// ===============================================================
// --------------   DS system [SDSX_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef DS_EQUATIONS
// ==============================================================
// DS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// DS_EQUATIONS==1 coupled    solver (u,v,w,p)
// DS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// local Femus class include -----------------------------------
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverDS.h"       // Navier-Stokes class header file
// #include "MGMesh.h"          // Mesh class
#include "MGFE.h"          // Mesh class



#define PRESS_0 (2.)



void MGSolDS::set_bc_matrix(
  DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
  int dir_maxnormal,                          ///<  normal dir
  int sur_toply[],                            ///< boundary topology map
  int el_ndof[],                              ///< number of volume dofs
  int elb_ndof[],                             ///< number of boundary dofs
  int elb_ngauss,                             ///<  number of surface gaussian points
  double normal[],                            ///< normal
  const int iaxis                              ///< axisymmetric
) {
//   const double delta=_FSI_parameter.DIST_FIX;
  double eff_stress=1.; double pressure;
  double vel_bound[1]; double xyz_g[DIMENSION];
//   const int elb_dof0=(NDOF_K>1)?elb_ndof[_pres_order]*NDOF_K:elb_ndof[_pres_order];
  double det2= _fe[2]->JacSur(elb_ngauss-1,_xxb_qnds, _InvJac2);// jacobian
  double Ipenalty=det2;                            // Dirichlet bc flagdouble xyz_bg[DIMENSION];


  //  ============================================================================================
  //  =================== Dirichlet boundary condition (Dirichlet=0; Neumann=1) ====================
  //  ============================================================================================

//     enum bound_cond {
//     interior=          11,     ///<
//
//     outflow=           22,     ///<  5(n) \f$  \tau_{nn}=p=0  \f$  +  5(t)  \f$ \tau_{nt}=0   \f$
//     pressure_outlet=   28,     ///<  5(n) \f$  \tau_{nn}=p=0  \f$  + 3(t) tg velocity inlet  \f$ {\bf v}.{\bf t}= 0  \f$
//
//     outflow_p=         32,     ///<  6(n)    \f$ \tau_{nn}=p=p_0  \f$ +  5(t)  \f$  \tau_{nt}=0 \f$
//     pressure_inlet=    38,     ///<  3(n)  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$ + 3(t) tg velocity inlet  \f$ {\bf v}.{\bf t}= 0  \f$
//
//     velocity_in0=      88,     ///<  3(n)  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$ + 3(t) tg velocity inlet  \f$ {\bf v}.{\bf t}= 0  \f$
//     velocity_tg0=      88,     ///<  3(n) normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$ +  3(t) tg velocity inlet  \f$ {\bf v}.{\bf t}= 0  \f$
//     wall=              66,     ///<  1(n) u=0  + 1(t) v0, w=0
//     velocity_in=       98,     ///<  4(n)  normal velocity inlet  \f$ {\bf v}.{\bf n}=  {\bf v}_0.{\bf n}  \f$ + 3(t) tg velocity inlet  \f$ {\bf v}.{\bf t}= 0  \f$
//     velocity_tg=       89,     ///<  3(n)  normal velocity inlet  \f$ {\bf v}.{\bf t}=0  \f$+ 4(t) tg velocity inlet  \f$ {\bf v}.{\bf t}=  {\bf v}_0.{\bf t}  \f$
//     velocity=          99,     ///<  4(n)  normal velocity inlet  \f$ {\bf v}.{\bf n}=  {\bf v}_0.{\bf n}  \f$ +4(t) tg velocity inlet  \f$ {\bf v}.{\bf t}=  {\bf v}_0.{\bf t}  \f$
//     wall_turb=         84,     ///< 3(n)  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$ + 7(t)  \f$  \tau_{nt}=\alpha u_t  \f$
//
//
//
//     simmx=              82    ///< 21  simmetry      \f$ (p{\bf n}++  5(t)  \f$  \tau_{nt}=0 \f$
// };




//      const int bc_el=_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM];
  int integral=0;
  // Dirichlet
  for (int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) { // Dirichlet ****************************************************
    for (int ivar=0; ivar< _nvars[2]; ivar++)    {

      const int  indx_row=sur_toply[lbnode]+ivar*el_ndof[2]; // solution indx
      const int  indx_sol= indx_row+_dir*el_ndof[2];// matrix index
      const int  bc_var=_bc_bd[ indx_row];
      const int  bc_bc=bc_var%10;
      const int  bc_rhs=((bc_var%2)== 1)? 1:0;


      if (bc_bc >5)  { // *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
        if (bc_var/10 ==1 )  { // normal  16-19 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (ivar+_dir==dir_maxnormal) {
            _bc_el[ indx_row]=0;   // set flag to Dirichlet
            const double valn=bc_rhs*_data_eq[2].ub[indx_sol+(_FF_idx[SDSX_F])*el_ndof[2]];
            // along axis  u_x u or v or w---------------------------------- 16 or 17 -------------------------------------------
            if (bc_bc<8) { KeM(indx_row,indx_row) = Ipenalty;     FeM(indx_row) = Ipenalty*valn; }
            else { // normal  u.n n(u*nx+v*ny+w*nz)-----------------------------------18  19 -----------------------------------------------
              const double nii=Ipenalty*normal[ivar+_dir];
              KeM(indx_row,indx_row) = nii;  FeM(indx_row) = nii*valn; //diagonal
//               for (int   jvar=ivar+_dir+1; jvar< _nDSdim+ivar+_dir; jvar++)    { // non-diagonal part  ->  u_n=n(u*nx+v*ny+w*nz)
//                 const int    jvar1=jvar%DIMENSION;    const  double nij= Ipenalty*normal[jvar1];
//                 const  double valj=bc_rhs*_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
// //                 KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += nij;     FeM(indx_row) += nij*valj; // unsplitted
// 	 	FeM(indx_row) +=0; //splitted
//               }// jvar loop
            } //   end u.n -----------------------------18 19----------------------------------------------------
          }// ivar==dir_maxnormal
        } // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else  {// Diriclet tangential  6-9 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (ivar+_dir != dir_maxnormal) {
            _bc_el[indx_row]=0;   // set flag to Dirichlet
            const double vart= bc_rhs*_data_eq[2].ub[indx_sol+_FF_idx[SDSX_F]*el_ndof[2]];
            if (bc_bc<8) {  // tangential axial  u_tx,uty ---------------------------------- 6 or 7 -------------------------------------------
              KeM(indx_row,indx_row) = Ipenalty;    FeM(indx_row) = Ipenalty*vart;//diagonal
            } // --------------------------------------------6  7--------------------------------------------------
            else { //  tangential u.t -----------------------------------8   9 -----------------------------------------------
              const double nii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir];
              KeM(indx_row,indx_row) +=  Ipenalty*(1-nii);        FeM(indx_row) +=  Ipenalty*(1-nii)*vart;
              // non-diagonal part Dirichlet tang  homogeneous only
//               for (int   jvar=ivar+_dir+1; jvar< _nDSdim+ivar+_dir; jvar++)    {
//                 const  int   jvar1=jvar%DIMENSION;       const double nij=Ipenalty*normal[jvar1]*normal[ivar+_dir];
//                 const double valtj=bc_rhs*_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
// //                 KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += -nij;       FeM(indx_row) += -nij*valtj;// unsplitted
// // 		FeM(indx_row) +=0; //splitted
//               }// jvar loop
            }//  end  tangential u.t -----------------------------------8   9 ----------------------------------------------
          }//ivar != dir_maxnormal
        } //   end u.t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      } // bc_bc >5  *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
      else {if (bc_bc==4) { integral=1; }}
    }// ebd ivar
  }// end loop lbnode Dirichlet

  // Numann pressure *********************************************************************************************
  // Pressure dirichlet boundary condition  ----------------  12  13 -------------------------
//   int iiop=0;
//   for (int i=0; i< elb_ndof[2]; i++)  if (_bc_bd[sur_toply[i]+ dir_maxnormal*NDOF_FEM]==12 || _bc_bd[sur_toply[i]+ dir_maxnormal*NDOF_FEM]==13) {
//       iiop=1;   if ( _bc_bd[sur_toply[i]+ dir_maxnormal*NDOF_FEM]==13)  integral=2;
//     }
// //   if(_bc_bd[surbc_bd[sur_toply[1]+ _nDSdim*NDOF_FEM]==13) iiop=1;
//   if (iiop==1) {
//     for (int i=0; i< elb_ndof[1]; i++) {
//       const int bc_rhs= ( integral==2) ? 1:0;   const int   indx=sur_toply[i]+ _nDSdim*NDOF_FEM;// volume dof index
//       _bc_el[indx]=0;     KeM(indx,indx)  +=Ipenalty;      FeM(indx)  += bc_rhs* Ipenalty* _data_eq[1].ub[sur_toply[i]];
//     } // end -------------------------------------------------------------------------------
//     if (integral==2 && _FSI_parameter.PINTP==1) {// ----------------------------------------------------
//       for (int ivar=0; ivar< _nvars[2]; ivar++)    {
//         for (int  qp=0; qp< elb_ngauss; qp++) { //gaussian integration loop (n_gauss)
//           double det     = _fe[2]->JacSur(qp,_xxb_qnds, _InvJac2);   // local coord _phi_g and jac
//           double  JxW_g=det*_fe[2]->_weight1[_nDSdim-2][qp];  // weight
//           _fe[2]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[2]);       // global coord _phi_g
//           _fe[1]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[1]);       // global coord _phi_g
// //             interp_el_bd_sol(_xx_qnds,sur_toply,el_ndof[2],0,_nDSdim,_phi_g[2],elb_ndof[2],_xyzg);
// //             if(iaxis==1) { JxW_g  *=_xyzg[0]; } // axisymmetric  (index ->0)
//           // Assemblying NS equation
//           for (int i=0; i< elb_ndof[2]; i++) if (ivar == dir_maxnormal) { // Assemblying NS equation
//               const int   indx_row = sur_toply[i]+ivar*el_ndof[2];      // volume dof index
//               const int   indx_var =  indx_row+_dir*el_ndof[2]; // volume dof index
//               double press=0.;   for (int j=0; j<elb_ndof[1]; j++)   press   += _phi_g[1][j]*_data_eq[1].ub[(_FF_idx[SDSX_F])*NDOF_FEM+sur_toply[j]];
//               FeM(indx_row)  += -1.*_bc_el[indx_var]* JxW_g* _phi_g[2][i]*normal[ivar+_dir]*press;
//             }// i-loop
//         }// end gaussian  integration
//       } // ivar
//     } // if integral--------------------------------------------------------------------------------------------
// 
//   }  //   iiop if
// end Neumann  pressure conditions *********************************************************************

  if (integral==1) {

    eff_stress = 1.;      //-- 4 or 5-
    for (int ivar=0; ivar< _nvars[2]; ivar++)    {

      for (int  qp=0; qp< elb_ngauss; qp++) { //gaussian integration loop (n_gauss)
        // quad/linear  [2]=quad [1]=linear------------------------------------
        double det     = _fe[2]->JacSur(qp,_xxb_qnds, _InvJac2);   // local coord _phi_g and jac
        double  JxW_g=det*_fe[2]->_weight1[_nDSdim-2][qp];  // weight
        _fe[2]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[2]);       // global coord _phi_g
        _fe[1]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[1]);       // global coord _phi_g
//             interp_el_bd_sol(_xx_qnds,sur_toply,el_ndof[2],0,_nDSdim,_phi_g[2],elb_ndof[2],_xyzg);
//             if(iaxis==1) { JxW_g  *=_xyzg[0]; } // axisymmetric  (index ->0)
        // Assemblying NS equation
        for (int i=0; i< elb_ndof[2]; i++) if (ivar+_dir  != dir_maxnormal) { // Assemblying NS equation

            const int  bc_var=_bc_bd[sur_toply[i]+ivar*NDOF_FEM];
            const int  bc_bc=bc_var%10;
            const int  bc_rhs=((bc_var%2)== 1)? 1:0;
            // set up row i
            const double phii_g  = _phi_g[2][i];
            const int   indx_var = sur_toply[i]+(ivar+_dir)*el_ndof[2]; // volume dof index
            const int   indx_row = sur_toply[i]+(ivar)*el_ndof[2];      // volume dof index
            const double dtJxW_g = (_bc_el[indx_row]== 0) ? 0.: JxW_g;
            double nii= normal[ivar+_dir]*normal[ivar+_dir];

            for (int j=0; j<elb_ndof[2]; j++) { // quad -quad ------------------------------------------------------------------------------------------
              const double phij_g= _phi_g[2][j];
              double vali=_bc_el[indx_row]*_data_eq[2].ub[indx_var+_FF_idx[SDSX_F]*NDOF_FEM];
              KeM(indx_row,indx_row) += dtJxW_g*(1.-nii)*eff_stress*phii_g*phij_g;
              FeM(indx_row) += bc_rhs*dtJxW_g*(1.-nii)*eff_stress*phii_g*vali;

//               for (int  jvar=ivar+_dir+1; jvar< _nDSdim+ivar+_dir; jvar++) {     // other components than ivar
//                 int jvar1 = jvar%DIMENSION;
//                 double nij =_bc_el[sur_toply[j]+(ivar*_dir)*NDOF_FEM]*normal[jvar1]*normal[ivar+_dir];
//                 double valj=_bc_el[sur_toply[j]+(ivar*_dir)*NDOF_FEM]*_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[i]];
// //                 FeM(indx_row) += bc_rhs*dtJxW_g*eff_stress*nij*valj;                KeM(indx_row,sur_toply[j]+jvar1*el_ndof[2]) += dtJxW_g*eff_stress*nij*phij_g*phii_g     ;// unsplitted
// //                   FeM(indx_row) += 0.   ;// splitted
//               }// jvar-loop
            }//  j-loop ---------------------------------------------------------------------------------------------------------
          }// i-loop
      }// end gaussian  integration
    } // ivar


  } // integral






  return;
}

// void MGSolDS::set_bc_matrix(
//   DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
//   int dir_maxnormal,                          ///<  normal dir
//   int sur_toply[],                            ///< boundary topology map
//   int el_ndof[],                              ///< number of volume dofs
//   int elb_ndof[],                             ///< number of boundary dofs
//   int elb_ngauss,                             ///<  number of surface gaussian points
//   double normal[],                            ///< normal
//   const int iaxis                              ///< axisymmetric
// ) {
// 
//   double eff_stress=1.; double pressure;
//   double vel_bound[1]; double xyz_g[DIMENSION];
// //   const int elb_dof0=(NDOF_K>1)?elb_ndof[_pres_order]*NDOF_K:elb_ndof[_pres_order];
//   double det2= _fe[2]->JacSur(elb_ngauss-1,_xxb_qnds, _InvJac2);// jacobian
//   double Ipenalty=det2;                            // Dirichlet bc flagdouble xyz_bg[DIMENSION];
// 
// 
//   //  ============================================================================================
//   //  =================== Dirichlet boundary condition (Dirichlet=0; Neumann=1) ====================
//   //  ============================================================================================
//   // enum bound_cond
//   // velocity_in0=      1,   velocity_tg0=      2,   wall=              3,
//   // velocity_in=       5,   velocity_tg=       6,   velocity=          7,
//   //  ============================================================================================
// 
//   for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
// 
//     int bc_dn=(_bc_vol[sur_toply[lbnode]]/10);
//     if(bc_dn==0) {  // dirichlet -------------------------------------------------------------------------------------------------------------------
// 
//       int bc_s=_bc_vol[sur_toply[lbnode]]%10;
//       //  (0;??1) normal              (0;?1?) tg  // (0;1??) nonhomoge neous
//       int bc_normal=(bc_s%2); int bc_tg=((bc_s&2)>>1); int bc_rhs=((bc_s&4)>>2);
//      
//         const int  indx_sol=sur_toply[lbnode]+_dir*el_ndof[2];//ivar=0;
//         const int  indx_row=sur_toply[lbnode];//ivar=0;
// 
//         if(_dir == dir_maxnormal)  {  // normal ---------------------------------------------------------------------------------------------------
//        
//           _bc_el[sur_toply[lbnode]]=0;
//           //          Dirichlet (homogenoues and non  homogenoues) Assemblying
//           const double nii=Ipenalty*normal[_dir]/**normal[ivar+_dir]*/;
//           KeM(indx_row,indx_row) = nii; //diagonal ->  u_n=n(u*nx+v*ny+w*nz)
//           FeM(indx_row) = bc_normal*bc_rhs*nii*_data_eq[2].ub[indx_sol+(_FF_idx[SDSX_F]+_dir)*el_ndof[2]]; // only normal non-homogeneous
//           for(int   jvar=_dir+1; jvar< _nDSdim+_dir; jvar++)    {  // non-diagonal part  ->  u_n=n(u*nx+v*ny+w*nz)
//             int      jvar1=jvar%DIMENSION;
//             const double nij= Ipenalty*normal[jvar1]/**normal[ivar+_dir]*/;
//             const double valj=_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
//             FeM(indx_row) += (bc_rhs*bc_normal-1)*nij*valj; //  KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += aij;
//           }// jvar loop
// 
//         } //  end  Dirichlet normal --------------------------------------------------------------------------------------------------------------------------
// 
//         else  {  // Non normal direction  ---------------------------------------------------------------------------------------------------------------
// 
//           _bc_el[sur_toply[lbnode]]=0; // storage bc for vol part
// 
//           //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
//           // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
//           //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
//           // diagonal part
//           const double nii= normal[_dir]*normal[_dir]*fabs(bc_tg-bc_normal);
//           KeM(indx_row,indx_row) =  Ipenalty*(1-nii);
//           FeM(indx_row) =  Ipenalty*bc_rhs*bc_tg*(1-nii)*_data_eq[2].ub[indx_sol+(_FF_idx[SDSX_F]+_dir)*el_ndof[2]];  
//           // non-diagonal part
//           for(int   jvar=_dir+1; jvar< _nDSdim+_dir; jvar++)    {
//             int      jvar1=jvar%_nDSdim;
//             const double nij=Ipenalty*normal[jvar1]*normal[_dir]*fabs(bc_tg-bc_normal);
//             const double valj=_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
//             FeM(indx_row) -= (bc_rhs*bc_tg-1)*nij*valj;
// //             FeM(indx_row) += aij*valj*fabs(bc_tg-bc_normal);
// 
// //             KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += -aij;
//           }// jvar loop
//         } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------
// 
//       }  // dc
//  
// 
//     // ==============================================================================================================
//     // simmetry for vector field
//     // simmx=            21,   simmy=            22,   simmz=            23,
//     // simmxy=           24,   simmxz=           25,   simmyz=           26 <-  not used
//     // ===============================================================================================================
//     if(bc_dn==2) {  // simmetry--------------------------------------------------------------------------------------------
//       int bc_s=_bc_vol[sur_toply[lbnode]]%10;
// 
//       if((bc_s%3) <4) { // simmx=21,24 (x-plane),simmy=22,25(y-plane);simmz=23,26(z-plane);
//         int sim_dir= (bc_s%3)-1 ;
//         if( sim_dir == dir_maxnormal && _dir==dir_maxnormal)  {  // normal ---------------------------------------------------------------------------------------------------
// //           const int  indx_sol=sur_toply[lbnode]+( sim_dir+_dir)*el_ndof[2];//ivar=0;
//           const int  indx_row=sur_toply[lbnode];//ivar=0;
//           _bc_el[sur_toply[lbnode]]=0;
//           //          Dirichlet (homogenoues and non  homogenoues) Assemblying
//           KeM(indx_row,indx_row) = 1.;FeM(indx_row) = 0.; //diagonal ->  u_n=n(u*nx+v*ny+w*nz)
//           
//         } //  end  ivar == dir_maxnormal------------------------------------------------------------------------------------------------------------
//       }
// 
//       
//       // -----------------------------------------------------------
//     }// dc  end dc==0
//   } // lnopde
// 
// 
// 
// 
// 
//   //  ============================================================================================
//   //  ====================================== Neumann  ============================================
//   // outflow=          10,   interior=         11,   pressure_outlet=  12,
//   // wall_turb=        13,   outflow_p=        14,   pressure_inlet=   16,
//   //  ============================================================================================
//   int ivar=0;
// //   for(int  ivar=0; ivar< _nvars[2]; ivar++)
//   {
//     if(_bc_vol[sur_toply[NDOF_FEMB-1]]/10 ==1) {  /*bc_dn=1; }*/
// //     if(bc_dn==1) {
//       // Neumann  -------------------------------------------------------------------------------------------------------------------
//       int bc_s=_bc_vol[sur_toply[NDOF_FEMB-1]]%10;
//       int bc_normal=(bc_s%2);  int bc_tg=((bc_s&2)>>1); int bc_rhs=((bc_s&4)>>2); // (0;?1?) tg (0;??1) normal (0;1??) nonhomogeneous
//    
// //       if(bc_s !=3)    for(int k=0; k<elb_dof0 ; k++) { _bc_el[_pres_order*sur_toply[k]+ _nDSdim*NDOF_FEM]= 0  ; }// pressure bc
//       int flag_normal=(ivar+_dir == dir_maxnormal)?1:0 ;
//       if(bc_s == 3) {
//         flag_normal= !flag_normal; // turb condition  bc_tg=1;bc_normal=1;
//       }
//       if(flag_normal)  {  // normal ---------------------------------------------------------------------------------------------------
//         double   det;
//         double sign=1.; if(normal[dir_maxnormal]<0) {sign=-1.;}
//         for(int  qp=0; qp< elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)
//           // quad/linear  [2]=quad [1]=linear------------------------------------
//           det     = _fe[2]->JacSur(qp,_xxb_qnds, _InvJac2);   // local coord _phi_g and jac
//           double  JxW_g=det*_fe[2]->_weight1[_nDSdim-2][qp]; // weight
//           _fe[2]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[2]);   // global coord _phi_g
//           _fe[1]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[1]);   // global coord _phi_g
//           interp_el_bd_sol(_xx_qnds,sur_toply,el_ndof[2],0,_nDSdim,_phi_g[2],elb_ndof[2],_xyz_g);
// 
//           if(iaxis==1) { JxW_g  *=_xyz_g[0]; } // axisymmetric  (index ->0)
// 
//           // Assemblying NS equation
//           for(int i=0; i< elb_ndof[2]; i++) {  // Assemblying NS equation
//             // set up row i
//             const double phii_g=_phi_g[2][i];
//             const int   indx_var=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
//             const int   indx_row=sur_toply[i]+(ivar)*el_ndof[2];// volume dof index
//             // boundary flag
//             int bc_v1=(int)(_bc_vol[sur_toply[i]]/10);
//             //                 int bc_s1=(int)_bc_bd[indx_row];
//             //                  bc_rhs   =((bc_s1&4)>>2); // (1??) -> nonhomogeneous
//             //                  bc_tg    =((bc_s1&2)>>1); // (?1?) -> tg
//             //                  bc_normal=(bc_s1%2);      // (??1) -> normal
//             double dtJxW_g=JxW_g*_bc_el[sur_toply[i]+(ivar+_dir)*NDOF_FEM];
// 
//             // wall_turb boundary conditions
//             double un=0;  for(int  kdim=0; kdim<_nDSdim; kdim++) { un +=_data_eq[2].ub[sur_toply[i]+(_FF_idx[SDSX_F]+kdim)*NDOF_FEM];}
// //               u_old[sur_toply[i]+kdim*NDOF_FEM]*normal[kdim]; }
//             double ut=0;  for(int  jvar=0; jvar<_nDSdim; jvar++) {
//               const double tmp=_data_eq[2].ub[sur_toply[i]+(_FF_idx[SDSX_F]+jvar)*NDOF_FEM]-normal[jvar]*un; ut +=tmp*tmp;
//             }
//             vel_bound[0] =sqrt(ut);    double u_tau= eval_var2(vel_bound);
// #ifdef TBK_EQUATIONS
//             eff_stress =u_tau*u_tau/fabs(vel_bound[0]+1.e-20);
// #endif
//             for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad ------------------------------------------------------------------------------------------
//               const double phij_g= _phi_g[2][j];
// //               if(j<elb_ndof[1]) {
// //                 double pressure_ctrl=1*_data_eq[_pres_order].ub[_pres_order*sur_toply[j]];
// //                 FeM(indx_row)  += -1.*sign*dtJxW_g*phii_g*_phi_g[1][j]*pressure_ctrl;//_refvalue[3]);
// //               }
//               KeM(indx_row,indx_row) +=   dtJxW_g*bc_tg*bc_normal*eff_stress*phij_g*phii_g;
//               for(int  jvar=ivar+_dir+1; jvar< _nDSdim+_dir; jvar++)    {  // u
// //                 KeM(indx_row,sur_toply[j]+(jvar%DIMENSION)*el_ndof[2]) +=    dtJxW_g*eff_stress*phij_g*phii_g*bc_normal*bc_tg;
//               }// jvar-loop
//             }//  j-loop ---------------------------------------------------------------------------------------------------------
//           }// i-loop
//         }// end gaussian  integration
// //         for(int i=0; i< elb_ndof[_pres_order]; i++) {  // Assemblying NS equation
// //           const int   indx=_pres_order*sur_toply[i]+ _nDSdim*NDOF_FEM;// volume dof index
// //           KeM(indx,indx)  =det*(1-_bc_el[indx]);
// //            if(bc_s >3) {
// // //             std::cout << xyz_g[0]  << "  " << xyz_g[1]    << "  "<<std::endl;
// //             double pressure_ctrl= _data_eq[_pres_order].ub[_pres_order*sur_toply[i]];
// //             FeM(indx)  =  det* pressure_ctrl*(1-_bc_el[indx]);
// //           }
// //         }
//       }//  end  normal direction
//       //  -------------------------------------------------------------------------------------------------------------------
//       else  {  // Non normal direction  -------------------------------------------------------------------------------
//         for(int lbnode=0; lbnode< elb_ndof[2]; lbnode++) {  // Assemblying NS equation
//           // set up row                int  iflag=0;
//           bc_s=(_bc_vol[sur_toply[lbnode]]%10)%4;
//           bc_tg=((bc_s&2)>>1);   // (0;?1?) tg
//           bc_normal=(bc_s%2);    // (0;??1) normal
//           bc_rhs=((bc_s&4)>>2);  // (0;1??) nonhomogeneous
//           if(_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM] != 0) {
//             const int  indx_sol=sur_toply[lbnode]+(ivar+_dir)*el_ndof[2];//ivar=0;
//             const int  indx_row=sur_toply[lbnode]+ivar*el_ndof[2];//ivar=0;
//             _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= ((bc_tg+bc_normal) == 0) ? 1:0 ;  // storage bc for vol part
//             //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
//             // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
//             //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
//             // diagonal part
//             KeM(indx_row,indx_row) += ((bc_tg+bc_normal)%2)*Ipenalty;
//             double aii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*((bc_tg+bc_normal));
//             KeM(indx_row,sur_toply[lbnode]+ivar*el_ndof[2]) += -Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(bc_tg+bc_normal);
//             FeM(indx_row) +=  bc_rhs*bc_tg*Ipenalty*_data_eq[2].ub[_FF_idx[SDSX_F]*NDOF_FEM+indx_sol]
//                               -bc_tg*bc_rhs*aii*_data_eq[2].ub[_FF_idx[SDSX_F]*NDOF_FEM+indx_sol];   // Dirichlet tang  homogeneous only
//             // non-diagonal part
//             for(int   jvar=ivar+_dir+1; jvar< _nDSdim+ivar+_dir; jvar++)    {
//               int      jvar1=jvar%DIMENSION;
//               double aij=_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]*Ipenalty*normal[jvar1]*normal[ivar+_dir]*(bc_tg+bc_normal);
//               double valj=_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]*_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
//               FeM(indx_row) += -bc_rhs*bc_tg*aij*valj;
// //               KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += -aij;
//               FeM(indx_row) +=aij*valj;
//             }// jvar loop
//             //                     bc_rhs*Ipenalty*u_tau[i]   + bc_yplus*Ipenalty*y1_plus[i];
//           } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------
//         }
//       }
//     }  //  ivar |||||
//   }  //  bc_dn |||||
// 
// 
// 
// 
// 
// 
// 
// //   } // l
// 
//   // ========================================= Dirichlet boundary conditions ==============================
//   return;
// }





// #endif


// #endif




#endif
