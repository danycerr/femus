#include "Equations_conf.h"
// ======================================================================================
#ifdef DS_EQUATIONS // 3D-2D Displacement equation
// ======================================================================================
// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverDS.h"

// config file --------------------------------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"                     // Mesh class
#include "MGSystem.h"                 // System class
#include "MGEquationsSystem.h"        // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ======================================================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif


// ======================================================
/// This function constructs the 3d-2D MGSolDSX_Y_Z class
MGSolDS::MGSolDS ( MGEquationsSystem& mg_equations_map_in,
                   const int vars_in[],
                   std::string eqname_in,
                   std::string varname_in ) :
    MGSolDA ( mg_equations_map_in,vars_in,eqname_in,varname_in ),
    // mesh params ------------
    _offset ( _mgmesh._NoNodes[_NoLevels-1] ), // mesh nodes
    // phys params ------------
    _dt ( _mgutils.get_par ( "dt" ) ), // parameter  dt
    _rhof ( mg_equations_map_in.get_par ( "rho0" ) ), // parameter  density reference
    _uref ( mg_equations_map_in.get_par ( "Uref" ) ), // parameter  vel reference
    _lref ( mg_equations_map_in.get_par ( "Lref" ) ), // parameter  length reference
    _Tref ( mg_equations_map_in.get_par ( "Tref" ) ), // parameter  temperature reference
    _rhos ( mg_equations_map_in.get_par ( "rhos" ) ) { // parameter solid density
    //  =================================================

    _FSI_idx=-1; //  navier-stokes-equations
    _DS_idx=-1;   //  energy-equations
    // class variable names
    _nDSdim=DIMENSION;
    _dir=0;

    if ( !varname_in.compare ( "dx" ) ) {
        _dir=0;
    }
    if ( !varname_in.compare ( "dy" ) ) {
        _dir=1;
    }
    if ( !varname_in.compare ( "dz" ) ) {
        _dir=2;
    }
#if DIMENSION ==2
    if ( _dir==2 ) {
        std::cout<<"Too many Dimension!!\n";
    }
#endif
    _var_names[0]=varname_in;
    _refvalue[0]=_lref;

    for ( int l=0; l<NDOF_FEM; l++ ) {
        _bc_vol[l]=-1;
        _bc_bd[l]=-1;
    }

// class solver type (SOLVERT  in MGSTconf.h)
//   for (int l=0;l<_NoLevels;l++) _solver[l]->set_solver_type(CGM);
    for ( int l=0; l<_NoLevels; l++ ) {
        _solver[l]->set_solver_type ( GMRESM );
    }
    return;
}

//  =====================================================================================
/// This function move the mesh points according to the coarser grid
//  =====================================================================================
void  MGSolDS::MoveMesh (
    const int Level  // Level <-
) {  // =================================================================================

    const int n_nodes=_mgmesh._NoNodes[Level];                    //number nodes
    const int n_elem=_mgmesh._NoElements[0][Level];               //number of elements
    const int  offset = _mgmesh._NoNodes[Level];                  // mesh nodes
    int  el_conn[NDOF_FEM];                                       // element connectivity
    int el_mat_nrows =0;
    int el_ndof[3];
    el_ndof[0]=0;
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];             // element nodes
    int        sur_toply[NDOF_FEMB];                              // boundary topology

    // coupling  basic system fields ------------------------------------------------------
    int idx_fsi= _data_eq[2].tab_eqs[FS_F];// Fluid Structure Interaction System-----------
    _FSI_idx= ( idx_fsi>=0 ) ? _data_eq[2].indx_ub[idx_fsi]:-1;

    int idx_ds= _data_eq[2].tab_eqs[SDSX_F];// Displacemnte in the x direction index--------
    _DS_idx= ( idx_ds>=0 ) ? _data_eq[2].indx_ub[idx_ds]:-1;
// ------------------------ computing ----------------------------------------------------
    for ( int ideg=1; ideg<3; ideg++ ) {
        el_ndof[ideg]= ( ( _nvars[ideg]>0 ) ?    _fe[ideg]->_NoShape[DIMENSION-1]:0 );
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    int el_mat_ncols = el_mat_nrows;                     //square matrixbc_Neum[0]=0;
    std::vector<int> el_dof_indices ( el_mat_ncols );    // element dof vector

    int Level_c=0;//_NoLevels-1;
    DenseVectorM LoD;//Local Displacement
    LoD.resize ( el_ndof[2] );
    for ( int proc=0; proc<_mgmesh._n_subdom; proc++ ) {
        int ndof_lev=0;
        for ( int pr=0; pr <proc; pr++ ) {
            int delta =_mgmesh._off_el[0][pr*_NoLevels+Level_c+1]
                       -_mgmesh._off_el[0][pr*_NoLevels+Level_c];
            ndof_lev +=delta;
        }
// ----------------------------------------------------------------------------------------
//---------------------------------find the coarse nodes-----------------------------------
// ----------------------------------------------------------------------------------------
        const int coarse_nel_e =_mgmesh._off_el[0][Level_c+_NoLevels*proc+1]; // start element
        const int coarse_nel_b =_mgmesh._off_el[0][Level_c+_NoLevels*proc];   // stop element
        int        coarse_node[NDOF_FEM];                                     // element connectivity
        std::cout;
        //coarse elemnt loop
        for ( int iel=0; iel < ( coarse_nel_e - coarse_nel_b ); iel++ ) {
            LoD.zero();
            //getting connectivity and coordinates
            _mgmesh.get_el_nod_conn ( 0,Level_c,iel,el_conn,_xx_qnds,proc );
            //getting boundary contions
            get_el_dof_bc ( Level_c,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd );

            int count=0;
            int written[NDOF_FEM]; //needed to have fine displacement on the solid
            for ( int inode=0; inode<el_ndof[2]; inode++ )     {
                coarse_node[inode]=  _mgmesh._node_map[Level_c][el_dof_indices[inode]];
                int glob_indx=_mgmesh._node_map[Level_c][el_dof_indices[inode]];
                if ( _bc_vol[inode]>1.5 ) {
                    LoD ( inode ) = ( *x_old[Level] ) ( glob_indx )- ( *x_oold[Level] ) ( glob_indx );
                    written[count]=inode;
                    count++;
                }//end solid
                else { //liquid
                    if ( inode < NDOF_P ) {
                        LoD ( inode ) = ( *x_old[Level] ) ( glob_indx )- ( *x_oold[Level] ) ( glob_indx );
                    }
                }

            }//end loop node
#if DIMENSION==2
            for ( int iside=0; iside< el_sides; iside++ )  {
                int check=0;
                for ( int idof=0; idof<NDOF_FEMB; idof++ ) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];     // local nodes
                }
                check=0;
                for ( int i=0; i<count; i++ ) if ( written[i]==sur_toply[NDOF_FEMB-1] ) {
                        check =1;
                    }
                if ( check==0 ) {
                    LoD ( sur_toply[NDOF_FEMB-1] ) = 0.5* ( LoD ( sur_toply[NDOF_FEMB-3] ) )
                                                     + 0.5* ( LoD ( sur_toply[NDOF_FEMB-2] ) );
                }

            }//end iside
#endif
#if DIMENSION==3
            for ( int iside=0; iside< el_sides; iside++ )  {
//                 int check=0;
                for ( int idof=0; idof<NDOF_FEMB; idof++ ) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];     // local nodes
                }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
//                 if (check==0) {
                if ( _bc_vol[sur_toply[4]]<1.5 )  LoD ( sur_toply[4] ) = 0.5* ( LoD ( sur_toply[0] ) )
                            + 0.5* ( LoD ( sur_toply[1] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
//                 if (check==0) {
                if ( _bc_vol[sur_toply[5]]<1.5 ) LoD ( sur_toply[5] ) = 0.5* ( LoD ( sur_toply[1] ) )
                            + 0.5* ( LoD ( sur_toply[2] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
//                 if (check==0) {
                if ( _bc_vol[sur_toply[6]]<1.5 )  LoD ( sur_toply[6] ) = 0.5* ( LoD ( sur_toply[2] ) )
                            + 0.5* ( LoD ( sur_toply[3] ) );
                if ( _bc_vol[sur_toply[7]]<1.5 )  LoD ( sur_toply[7] ) = 0.5* ( LoD ( sur_toply[3] ) )
                            + 0.5* ( LoD ( sur_toply[0] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
//                 if (check==0) {
                if ( _bc_vol[sur_toply[8]]<1.5 )   LoD ( sur_toply[8] ) = 0.25* ( LoD ( sur_toply[0] ) )
                            +0.25* ( LoD ( sur_toply[1] ) )
                            + 0.25* ( LoD ( sur_toply[2] ) )
                            + 0.25* ( LoD ( sur_toply[3] ) );
//                 }

            }//end iside
#endif
            int check=0;
            for ( int i=0; i<count; i++ ) if ( written[i]==NDOF_FEM-1 ) {
                    check =1;     //last node
                }
            if ( check==0 ) {
                LoD ( NDOF_FEM-1 ) =0;
                for ( int jnode=0; jnode<NDOF_P; jnode++ ) {
                    LoD ( NDOF_FEM-1 ) +=LoD ( jnode ) /NDOF_P;
                }
            }

            for ( int inode=0; inode < el_ndof[2] ; inode++ ) {
                disp[Level]->set ( _mgmesh._node_map[Level_c][el_dof_indices[inode]],LoD ( inode ) );
            }  //setting the global displacement
        }//end loop iel
        el_dof_indices.clear();
    }//end iproc
//      end coarse elemnt loop
//              end find the coarse nodes
// ----------------------------------------------------------------------------------------
//--------------------------------- finer level -------------------------------------------
// ----------------------------------------------------------------------------------------
    for ( int ilev=Level_c; ilev<_NoLevels; ilev++ ) {
        for ( int proc=0; proc<_mgmesh._n_subdom; proc++ ) {

            int ndof_lev=0;
            for ( int pr=0; pr <proc; pr++ ) {
                int delta =_mgmesh._off_el[0][pr*_NoLevels+ilev+1]-_mgmesh._off_el[0][pr*_NoLevels+ilev];
                ndof_lev +=delta;
            }
            const int nel_e =_mgmesh._off_el[0][ilev+_NoLevels*proc+1]; // start element
            const int nel_b =_mgmesh._off_el[0][ilev+_NoLevels*proc];   // stop element
            for ( int iel=0; iel < ( nel_e - nel_b ); iel++ ) {
                LoD.zero();
                //getting connectivity and coordinates
                _mgmesh.get_el_nod_conn ( 0,ilev,iel,el_conn,_xx_qnds,proc );
                //getting boundary contions
                get_el_dof_bc ( ilev,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd );

                int count=0;
                int written[NDOF_FEM]; //needed to have fine displacement on the solid
                int countg=0;
                for ( int inode=0; inode<el_ndof[2]; inode++ )     {
                    int dof_idx=_mgmesh._node_map[ilev][el_dof_indices[inode]];
                    if ( _bc_vol[inode]>1.5 ) {
                        written[count]=inode;
                        count++;
                        {
                            LoD ( inode ) = ( *x_old[Level] ) ( dof_idx )- ( *x_oold[Level] ) ( dof_idx );
                        }
                    }//end solid
                    else { //liquid
                        if ( inode < NDOF_P ) {
                            LoD ( inode ) = ( *disp[Level] ) ( dof_idx );
                        }
                    }

                }//end loop node
#if DIMENSION==2
                for ( int iside=0; iside< el_sides; iside++ )  {
                    int check=0;
                    for ( int idof=0; idof<NDOF_FEMB; idof++ ) {
                        sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];     // local nodes
                    }
                    check=0;
                    for ( int i=0; i<count; i++ ) if ( written[i]==sur_toply[NDOF_FEMB-1] ) {
                            check =1;
                        }
                    if ( check==0 ) {
                        LoD ( sur_toply[NDOF_FEMB-1] ) = 0.5* ( LoD ( sur_toply[NDOF_FEMB-3] ) )
                                                         + 0.5* ( LoD ( sur_toply[NDOF_FEMB-2] ) );
                    }
                }
#endif

#if DIMENSION==3
                for ( int iside=0; iside< el_sides; iside++ )  {
//                 int check=0;
                    for ( int idof=0; idof<NDOF_FEMB; idof++ ) {
                        sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];     // local nodes
                    }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
//                 if (check==0) {
                    if ( _bc_vol[sur_toply[4]]<1.5 )  LoD ( sur_toply[4] ) = 0.5* ( LoD ( sur_toply[0] ) )
                                + 0.5* ( LoD ( sur_toply[1] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
//                 if (check==0) {
                    if ( _bc_vol[sur_toply[5]]<1.5 ) LoD ( sur_toply[5] ) = 0.5* ( LoD ( sur_toply[1] ) )
                                + 0.5* ( LoD ( sur_toply[2] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
//                 if (check==0) {
                    if ( _bc_vol[sur_toply[6]]<1.5 )  LoD ( sur_toply[6] ) = 0.5* ( LoD ( sur_toply[2] ) )
                                + 0.5* ( LoD ( sur_toply[3] ) );
                    if ( _bc_vol[sur_toply[7]]<1.5 )  LoD ( sur_toply[7] ) = 0.5* ( LoD ( sur_toply[3] ) )
                                + 0.5* ( LoD ( sur_toply[0] ) );
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
//                 if (check==0) {
                    if ( _bc_vol[sur_toply[8]]<1.5 )   LoD ( sur_toply[8] ) = 0.25* ( LoD ( sur_toply[0] ) )
                                +0.25* ( LoD ( sur_toply[1] ) )
                                + 0.25* ( LoD ( sur_toply[2] ) )
                                + 0.25* ( LoD ( sur_toply[3] ) );
//                 }

                }//end iside
#endif

                int check=0;
                for ( int i=0; i<count; i++ ) if ( written[i]==NDOF_FEM-1 ) {
                        check =1;     //last node
                    }
                if ( check==0 ) {
                    for ( int jnode=0; jnode<NDOF_P; jnode++ ) {
                        LoD ( NDOF_FEM-1 ) +=LoD ( jnode ) /NDOF_P;
                    }

                }
                int check2=0;
                for ( int inode=0; inode < el_ndof[2] ; inode++ ) {
                    disp[Level]->set ( _mgmesh._node_map[ilev][el_dof_indices[inode]],LoD ( inode ) );
                }
            }//end loop iel

            el_dof_indices.clear();
        }//end iproc
    }//end ilev
    disp[Level]->close();
    return;
}
//  =====================================================================================
/// This function assembles the matrix and the rhs:
//  =====================================================================================
void  MGSolDS::GenMatRhs (
    const double time, // time  <-
    const int Level,  // Level <-
    const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // =================================================================================

// ----------------------------------------------------------------------------------------
//--------------------------------- Set Up ------------------------------------------------
// ----------------------------------------------------------------------------------------
    const int  ndim = DIMENSION;                                           //dimension
    int        el_conn[NDOF_FEM];                                          // element connectivity
    int        el_neigh[NDOF_FEM];                                         // element connectivity
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
    int        sur_toply[NDOF_FEMB];                                       // boundary topology
    int        flag_group[NDOF_FEM];
    int        flag_sur_group[NDOF_FEMB];
    int fl_int;    int phase;
    double l_old[DIMENSION*NDOF_FEM];
    double  p_old;
    // gauss integration  ------------------
    const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];            //elem gauss points
    int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];            //elem gauss points
    double det[3],JxW_g[3];                                     // Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION];                             // global derivatives at g point
    double x_m[DIMENSION]; double normal[DIMENSION];
    //number of constant[0]-linear[1]-quadratic[2] element dof
    int el_ndof[3]; el_ndof[0]=0;
    int elb_ndof[3];elb_ndof[0]=0;
    int el_mat_nrows =0;
    
    for ( int ideg=1; ideg<3; ideg++ ) {                            //     ...
        el_ndof[ideg]= ( ( _nvars[ideg]>0 ) ?    _fe[ideg]->_NoShape[ndim-1]:0 );              //   computing
        elb_ndof[ideg]= ( ( _nvars[ideg]>0 ) ?_fe[ideg]->_NoShape[ndim-2]:0 );            //     ...
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };

    int el_mat_ncols = el_mat_nrows;                     //square matrixbc_Neum[0]=0;
    std::vector<int> el_dof_indices ( el_mat_ncols );   // element dof vector
    MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *> ( &_mgmesh );
    _FSI_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[FS_F]];                    // Fluid-Structure Interaction equation
//    const int fsi_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[TVX_F]];        // Energy equation
    _DS_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[SDSX_F]];                   // Displacement equation
    const int p_idx=_data_eq[1].indx_ub[_data_eq[1].tab_eqs[P_F]];              // Pressure equation

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) -----------
    A[Level]->zero();
    if ( mode ==1 ) {
        b[Level]->zero();
    }
    DenseMatrixM KeM;
    DenseVectorM FeM;
    KeM.resize ( el_mat_nrows,el_mat_ncols );
    FeM.resize ( el_mat_nrows );
    int ndof_lev=0;
    for ( int pr=0; pr <_mgmesh._iproc; pr++ ) {
        int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
        ndof_lev +=delta;
    }
#if FSI_EQUATION%2==0
    int  el_ndofp=_fe[1]->_NoShape[ndim-1];
#endif
    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for ( int iel=0; iel < ( nel_e - nel_b ); iel++ ) {

        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        // geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (xx_qnds)
        _mgmesh.get_el_nod_conn ( 0,Level,iel,el_conn,_xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides,0,Level,iel,el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd );
        for ( int idim=0; idim<DIMENSION; idim++ ) {
            x_m[idim] =0.;
            for ( int d=0; d< NDOF_FEM; d++ ) { // element nodes xxg (DIM)
                x_m[idim] +=_xx_qnds[idim*NDOF_FEM+d]/NDOF_FEM;
                flag_group[d]=ext_mesh->_bc_id[el_dof_indices[d]];
            }
        }
// ----------------------------------------------------------------------------------------
//--------------------------------- external field ----------------------------------------
// ----------------------------------------------------------------------------------------
        for ( int deg=2; deg<3; deg++ ) {
            for ( int eq=0; eq<_data_eq[deg].n_eqs; eq++ ) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub );
            }
        }
#if FSI_EQUATIONS%2==0 // pressure as external field (projection or splitting)
        _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol ( 0,1,el_ndofp,el_conn, offset,0,_data_eq[1].ub );
        _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol ( 0,1,el_ndofp,el_conn, offset,1,_data_eq[1].ub );
#endif
        //  external cell properties -------------------------------------
        for ( int d=0; d< _nDSdim*NDOF_FEM; d++ ) {
            _bc_el[d]=1;
            l_old[d]= _data_eq[2].ub[_DS_idx*NDOF_FEM+d];
        }
        phase= ( ext_mesh->_mat_id[iel+nel_b]==2 ) ?0:1;

// ----------------------------------------------------------------------------------------
//--------------------------------- Fluid and Solid boundary ------------------------------
// ----------------------------------------------------------------------------------------
        for ( int iside=0; iside< el_sides; iside++ )  { 
            if ( el_neigh[iside] == -1 ) {

                for ( int idof=0; idof<NDOF_FEMB; idof++ ) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    int idofb=sur_toply[idof];//flag_sur_group[idof]=ext_mesh->_bc_id[el_dof_indices[idofb]];
                    for ( int idim=0; idim<DIMENSION; idim++ ) {
                        _xxb_qnds[idim*NDOF_FEMB+idof]=_xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                    }
                }
                // det hyperplane
                double  det2= _fe[2]->JacSur ( elb_ngauss-1,_xxb_qnds, _InvJac2 ); // jacobian
                double Ipenalty2=det2/_dt;                            // Dirichlet bc flag
                // normal
                _fe[2]->normal_g ( _xxb_qnds,x_m,normal );
                int dir_maxnormal = ( fabs ( normal[0] ) >fabs ( normal[1] ) ) ?0:1 ;
                dir_maxnormal= ( fabs ( normal[dir_maxnormal] ) >fabs ( normal[DIMENSION-1] ) ) ? dir_maxnormal:DIMENSION-1;


                set_bc_matrix ( KeM,FeM,dir_maxnormal,sur_toply,
                                el_ndof,elb_ndof, elb_ngauss,
                                l_old, normal,Ipenalty2 ); //  boundary conditions
                //std::cout << "KeM" << endl; cout << KeM << endl;cout << "FeM" << endl; cout << FeM << endl;
            } //end if side
        } //   end for fluid boundary ---------------------------------------------------------------
//        cout << "KeM" << endl; cout << KeM << endl;cout << "FeM" << endl; cout << FeM << endl;

// --------------------------------------------------------------------------------------
//--------------------------------- Fluid and Solid Volume   ----------------------------
// --------------------------------------------------------------------------------------
        if ( phase==0 ) { // ------------------------ fluid -----------------------------
            matrixrhsvol_liq_ds ( KeM,FeM, el_ndof,flag_group );
        } // ----------------------------------------------------------------------------

        else { // --------------------  solid -------------------------------------------
            matrixrhsvol_sol_ds ( KeM,FeM, el_ndof );
        }// ------------------------------------------------------------------------------------------

        /// e) Global assemblying energy equation
        A[Level]->add_matrix ( KeM,el_dof_indices );               // global matrix
        if ( mode == 1 ) {
            b[Level]->add_vector ( FeM,el_dof_indices );           // global rhs
        }
//       std::cout << "KeM  vol" << endl;cout << KeM << endl;cout << "FeM vol " << endl;cout << FeM << endl;
    } // end of element loop
    // clean
    el_dof_indices.clear();
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled (DS)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

    return;
}



// =======================================================================================
/// This function controls the assembly and the solution of the DS_equation system:
// =======================================================================================
void MGSolDS::MGTimeStep ( const double time, const int /*iter*/ ) { // ------------------

    std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
    /// B) Assemblying of the rhs (top level in surface and volume with MGSolNS::GenRhs,GenRhsB),
#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif
    GenMatRhs ( time,_NoLevels-1,1 );
    A[_NoLevels-1]->close();
    /// C) Assemblying of the  MGmatrices (MGSolNS::GenMatrix),
    for ( int Level = 0 ; Level < _NoLevels-1; Level++ ) {
        GenMatRhs ( time,Level,0 ); // matrix
        A[Level]->close();
    }
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
    /// E) Solution of the linear MGsystem (MGSolNS::MGSolve).
    MGSolve ( 1.e-6,40 );
#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
    x_oold[_NoLevels-1]->zero();
    x_old[_NoLevels-1]->localize ( *x_oold[_NoLevels-1] );
    x[_NoLevels-1]->localize ( *x_old[_NoLevels-1] );

// ----------------------------------------------------------------------------------------
//--------------------------------- moving mesh -------------------------------------------
// ----------------------------------------------------------------------------------------
    const int flag_moving_mesh = _mgutils.get_par ( "moving_mesh" );
    const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
    const int n_elem=_mgmesh._NoElements[0][_NoLevels-1];
    const int offsetp=_dir*n_nodes;
    if ( flag_moving_mesh ) {
#ifdef FINE_MOVEMENT
        /// E) mesh update
        const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
        const int offsetp=_dir*n_nodes;
        for ( int inode=0; inode<n_nodes; inode++ ) {
            double disp= ( *x_old[_NoLevels-1] ) ( inode )- ( *x_oold[_NoLevels-1] ) ( inode );
            _mgmesh._xyz[inode+offsetp] += disp;
            _mgmesh._dxdydz[inode+offsetp] = disp;
        }
#endif
#ifdef COARSE_MOVEMENT
        MoveMesh ( _NoLevels-1 );
        for ( int inode=0; inode<n_nodes; inode++ ) {
            _mgmesh._xyz[inode+offsetp] += ( *disp[_NoLevels-1] ) ( inode );
            _mgmesh._dxdydz[inode+offsetp]= ( *disp[_NoLevels-1] ) ( inode );
        }
        MoveMesh ( _NoLevels-1 );
#endif
    }
    // ==================================================================================
    return;
}
// ======================================================================================
// ======================================================================================


// ======================================================================================
/// This function controls the assembly and the solution of the DS_equation system:
/// In separete routines
// ======================================================================================
void MGSolDS::MGTimeStep_nl_setup ( const double time, const int /*iter*/ ) { // --------

    return;
}
int MGSolDS::MGTimeStep_nl_iter ( const double time, const int /*iter*/ ) { // ----------

    return 0;
}
void MGSolDS::MGTimeStep_nl_sol_up ( const double time, const int /*iter*/ ) { // -------

    return;
}
// ======================================================================================
// ======================================================================================

#endif //DS_EQUATIONS


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
