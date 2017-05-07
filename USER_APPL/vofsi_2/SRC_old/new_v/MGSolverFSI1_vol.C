// ===============================================================
// --------------   NAVIER-STOKES system [FS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef FSI_EQUATIONS
#if FSI_EQUATIONS==1
// ==============================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverFSI.h"       // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#ifdef   TWO_PHASE
 #include "MGSolverCC.h" 
#endif

// Thermodinamical Properties ==========================================================================================
// constant
double  _FSIdensity(double T,double T_r) {return 1.;} // water (t K) T_ref
double  _FSIkviscosity(double T,double T_r) {return 1.;} // lead T K
// water
//     double  _NSdensity(double T,double T_r){return (1. - (T-3.9863)*(T-3.9863)*(T+288.9414)/(508929.2*(T+68.12963)))/
//         (1. - (T_r-3.9863)*(T_r-3.9863)*(T_r+288.9414)/(508929.2*(T_r+68.12963)));} // water (t C) T_ref
//      double  _NSkviscosity(double T){return 1.;} //  T K
//  // Lead
//     double  _NSdensity(double T,double T_r){return (11441-1.2795*T)/(11441-1.2795*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return return 4.55e-4*exp(1069/T);} // lead T K
//       // LBE
//     double  _NSdensity(double T,double T_r){return (11065-1.293*T)/(11065-1.293*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return 4.94e-4*exp(754.1/T);} // lead T K
//
//  ======================================================================================================================

void MGSolFSI::get_el_field_data(
  int iel, int Level,
  int el_conn [], int offset,int el_ndof[],int ndof_lev,
  double u_old[],  double u_oold[],   double u_nl[],
  double p_proj[], double dp_proj[]
) {

  int el_ndofp=el_ndof[1];

  // element nodes coordinates ----------------------------------------------------------------------------------------
//   for(int idim=0; idim<DIMENSION; idim++) for(int d=0; d< NDOF_FEM; d++) {
//       _data_eq[2].ub[idim*NDOF_FEM+d]=_xx_qnds[idim*NDOF_FEM+d];
//     }
  // external fields (from constant 0 to quadratic 2) -------------------------------------------------------------------
  for(int deg=0; deg<3; deg++)
    for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
      _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                           el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
    }

  //  internal field data (NS) -------------------------------------------------------------------------------------------

#if FSI_EQUATIONS==0 // pressure as external field (projection)// --------------------------------------------------------------
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_nonl_sol(0,_nNSdim,el_ndof[2],el_conn, offset,0,u_nl);
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0,1,el_ndofp,el_conn, offset,0,p_proj);
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndofp,el_conn, offset,0,dp_proj);

  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_sol(0,_nNSdim,el_ndof[2],el_conn, offset,0,u_old);
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_oldsol(0,_nNSdim,el_ndof[2],el_conn, offset,0,u_oold);
#endif

#if FSI_EQUATIONS==1 //   coupled ----------------------------------------------------------------------------------------
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_sol(_nNSdim,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub);    // pressure
   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_sol_piece(_nNSdim,1,el_ndof[0],iel, offset,0,_data_eq[0].ub);    // pressure
 
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_sol(0,_nNSdim,el_ndof[2],el_conn, offset,0,u_old);     // old vel
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[FS_F]]->get_el_nonl_sol(0,_nNSdim,el_ndof[2],el_conn, offset,0,u_nl);   //  non linear vel
#endif

#if FSI_EQUATIONS==2 // pressure as external field (splitting) -----------------------------------------------------------
  for(int kdim=0; kdim< _nNSdim; kdim++) {
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[_FF_idx[FS_F]+kdim]]->get_el_nonl_sol(0,1,el_ndof[2],el_conn,offset,kdim,u_nl);
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[_FF_idx[FS_F]+kdim]]->get_el_sol(0,1,el_ndof[2],el_conn, offset,kdim,u_old);
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[_FF_idx[FS_F]+kdim]]->get_el_oldsol(0,1,el_ndof[2],el_conn, offset,kdim,u_oold);
  }
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0,1,el_ndofp,el_conn, offset,0,p_proj);       // old pressure
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndofp,el_conn, offset,0,dp_proj);   // dp pressure
#endif

  return;


}






// ==============================================================================================
void MGSolFSI::matrixrhs_liq_vol(
  DenseMatrixM &KeM, DenseVectorM &FeM,
  int el_ndof[],
  double u_old[],  double u_oold[],   double u_nl[], double p_proj[], double dp_proj[],
  const int unsteady_flag, const int axysim,
  const double les, const int mode, int flag_group[]
) {

  double dphijdx_g2[DIMENSION],dphijdx_g1[DIMENSION];
  double dphiidx_g2[DIMENSION],dphiidx_g1[DIMENSION];
  double val_tbg[1];
  double vel_gddx[DIMENSION*DIMENSION*DIMENSION];
  double vel_gdx[DIMENSION*DIMENSION];
  double vela_g[DIMENSION],u_nlg[DIMENSION]; 
 
  int el_ndofp=el_ndof[1];
  double rho=1.;
  int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim-1];                //elem gauss points
  const int el_ndof2=el_ndof[2];
  double ff[NDOF_FEM*DIMENSION];
   int ord[NDOF_FEM*DIMENSION];
   int block=CCLEV;
  
  
  for(int qp=0; qp<  el_ngauss; qp++) {
    // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
    // quadratic continuous (2)  (velocity)
    const double det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);     // quadratic Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[ _nNSdim-1][qp];           // quadratic weight
    _fe[2]->get_phi_gl_g(_nNSdim,qp,_phi_g[2]);                     // quadratic shape function
    _fe[2]->get_dphi_gl_g(_nNSdim,qp,_InvJac2,_dphi_g[2]);          // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nNSdim,qp,_InvJac2,_ddphi_g[2]);        // global second derivatives
    // continuous linear (1) (linear pressure)
    const double det1      = _fe[1]->Jac(qp,_xx_qnds,_InvJac1);     // linear Jacobian
    double JxW_g1 =det1*_fe[1]->_weight1[ _nNSdim-1][qp];           // linear weight
    _fe[1]->get_phi_gl_g(_nNSdim,qp,_phi_g[1]);                     // linear shape funct
    _fe[1]->get_dphi_gl_g(_nNSdim,qp,_InvJac1,_dphi_g[1]);          // global coord deriv
    // discontinuous (0) (disc pressure)
    if(_nvars[0]>0) { _fe[0]->get_phi_gl_g(_nNSdim,qp,_phi_g[0]);}  // piecewise shape function

    interp_el_sol(_xx_qnds,0,_nNSdim,_phi_g[2],el_ndof2,_xyzg);
        _msolcc->get_2surten(_xyzg,ff,ord);
    double cc=_msolcc->get_2phase(block, _xyzg); // new cc
    rho = 100.*cc*1+(1-cc)*1.;            // new cc 
    double mu_tp= 1000.*cc+(1-cc);       // new cc
    
    // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------

    // quadratic fields (velocity)
    interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2]);    // field _ub_g[2][DIM]
    interp_el_gdx(u_old,0, _nNSdim,_dphi_g[2],el_ndof2,vel_gdx);    // derivatives  vel_gdx[DIM][DIM]
    interp_el_gddx(u_old,0, _nNSdim,_ddphi_g[2],el_ndof2,vel_gddx);
    // linear fields (pressure) -> projection and segregated
// #if (FSI_EQUATIONS%2==0)
//     interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndofp,_ub_g[1]);
// #endif
    
    if(axysim==1) {JxW_g2  *=_xyzg[0]; JxW_g1  *=_xyzg[0];}  // std::cout<<_ub_g[2][0]<<"   "<<std::endl;

    // Velocity, Reynolds and upwind  --------------------------------------------------------------------------
    double mod2_vel=1.e-20;        _sP=1.e-20;    double mu=1.; // Velocity
    for(int idim=0; idim<  _nNSdim; idim++) {
      vela_g[idim]=0.;   u_nlg[idim]=0.; // old and  non linear Velocity at gaussina point qp
      for(int k=0; k< NDOF_FEM; k++) {
        u_nlg[idim] += u_nl[k+idim*NDOF_FEM]*_phi_g[2][k];
        vela_g[idim] += u_old[k+idim*NDOF_FEM]*_phi_g[2][k];
      }
      mod2_vel += u_nlg[idim]*u_nlg[idim]; // module velocity
      for(int jdim=0; jdim<  _nNSdim; jdim++) {
        _sP += 0.25*(vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim])*
               (vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim]);
      }
    }
    mod2_vel =sqrt(mod2_vel);
    // upwind term -> f_upwind
    // h_eff=2/sum|s.dN| and f_upwind ---------------------------------------------------------------------------------
    double h_eff=1.e-21;
    for(int i=0; i<el_ndof2; i++) {
      double hh=1.e-20; for(int idim=0; idim< _nNSdim; idim++) { hh += u_nlg[idim]*_dphi_g[2][i+idim*el_ndof2]/mod2_vel; }
      h_eff += fabs(hh);
    }
    h_eff=2./h_eff; if(h_eff<1.e-10) {h_eff=1. ; std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";   }

    double Pe_h=0.5*mod2_vel*h_eff/_IRe;
    double a_opt=(1./tanh(Pe_h)-1./Pe_h); if(a_opt >1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
    double f_upwind=rho*0.5*a_opt*h_eff/(mod2_vel);   // upwind

    // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
    if(_FF_idx[T_F]>=0) {
      double Temp_g=_ub_g[2][_FF_idx[T_F]];   rho *= _FSIdensity(Temp_g,0);  mu *=  _FSIkviscosity(Temp_g,0);
    }

    _mu_turb=0.;// eff visc turb at g pt
    // -------------------- Turbulence [K_F] -> (quad,_indx_eqs[K_F]) -------------------------------------
    if(_FF_idx[K_F]>=0) {// turbulence _mu_turb evaluation
      _kappa_g[0]= _ub_g[2][_FF_idx[K_F]];  _kappa_g[1]= _ub_g[2][_FF_idx[K_F]+1];
//       _mu_turb=eval_var1(val_tbg);
      _mu_turb= _kappa_g[0]/_kappa_g[1];
    }
    double IRe_eff  =_IRe*mu_tp/*_mu_turb*/ /*+_IRe*les*h_eff*h_eff*sqrt(2*_sP) */ ;   //  visc_eff  at g pt

    // ===========================================================================================================
    //                                       D) Assembling NS equation
    // ===========================================================================================================
    for(int i=0; i< el_ndof2; i++) {    // +++++++++++
      // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
      const double phii_g=_phi_g[2][i];
       int i_flag_i=1;
//        if(fabs(flag_group[i])>999){
//          i_flag_i=0;
//          }
      // Advection and Laplacian operator (supg)
      for(int  ivar=0; ivar< _nvars[2]; ivar++)    {
        int    indx=i+ivar*el_ndof2; // ivar=0;
        double dtxJxW_g=JxW_g2*fabs(_bc_el[indx])*i_flag_i;
        double Adv_expl=0.,Lap_expl=0.,supglapexpl=0.; double  Phi_supg=0.;
        for(int idim=0; idim< _nNSdim; idim++) {
          int ix=(ivar+_dir)*_nNSdim+idim;
          dphiidx_g2[idim]=_dphi_g[2][i+idim*el_ndof2];
         if(_bc_el[indx]>0) Phi_supg +=_FSI_parameter.SUPG*f_upwind*u_nlg[idim]*dphiidx_g2[idim]; //f_upwind
          Adv_expl    +=  vela_g[idim]*vel_gdx[ix];
          Lap_expl    +=  IRe_eff*vel_gdx[ix]*dphiidx_g2[idim];
          supglapexpl +=  IRe_eff*vel_gddx[ix*_nNSdim+idim];
        }
        // -------------------------------------- Assemblying rhs ------------------------------------------------------
//      KeM(indx,indx) +=JxW_g2*(1-i_flag_i)*fabs(_bc_el[indx]);  
//      FeM(indx) +=ivar*1.*JxW_g2*(1-i_flag_i)*fabs(_bc_el[indx]);
        if(mode == 1) {
          FeM(indx)  +=  dtxJxW_g* (
                           unsteady_flag* rho*vela_g[ivar+_dir]*(phii_g+Phi_supg) /_dt     // time
                           + rho*_IFr*_dirg[ivar+_dir]* (phii_g+Phi_supg)   // x-gravity
#ifdef Nat_Conv
                           + rho* (_ub_g[2][T_idx]-_T_nc) *_grav[ivar+_dir]*_beta_mat* (phii_g+Phi_supg)
#endif
                         );
        }
//
        //-------------------------------------- Assemblying matrix ---------------------------------------------------
//
        for(int j=0; j<el_ndof2; j++) {
          const double phij_g= _phi_g[2][j];
            int i_flag_j=1; //if(fabs(flag_group[j])>999) i_flag_j=0;
         
          // set up
          double Lap_g=0.,Adv_g=0., Div_g= 0., Supg_lap=0.;
          if(axysim==1) {  Lap_g =2.*(1-(ivar+_dir))*IRe_eff*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]); }      // axysimmetry
         
          for(int kdim=0; kdim< _nNSdim; kdim++) {
            dphijdx_g2[kdim] =_dphi_g[2][j+kdim*el_ndof2];
            Adv_g    += u_nlg[kdim]*dphijdx_g2[kdim]; // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
            Lap_g    += (IRe_eff+_FSI_parameter.UPWIND*f_upwind*u_nlg[kdim]*u_nlg[kdim]) *dphijdx_g2[kdim]*dphiidx_g2[kdim];
            Supg_lap += IRe_eff*_ddphi_g[2][j* _nNSdim* _nNSdim+kdim* _nNSdim+kdim];
          }
           FeM(indx)  +=  dtxJxW_g*0.072*(ff[j+ivar*9]*phij_g)* (phii_g+Phi_supg);
          //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
          KeM(indx,j+ivar*el_ndof2) +=dtxJxW_g*rho*i_flag_j* (
                                        (Adv_g+ unsteady_flag*phij_g/_dt)* (phii_g+Phi_supg) // time + advection
                                        + Lap_g                       // viscous Laplacian
                                        -Supg_lap*Phi_supg            // SUPG viscous Laplacian
                                         + IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir] // viscous tensor
                                      );
          // --------------------------------- Block +1 [2-6-7] --------------------------------------------------------------
          // out of diagonal viscous tensor
          const int idimp1= (ivar+1+_dir) % _nNSdim;
           KeM(indx,j+idimp1*el_ndof2) +=   i_flag_j* dtxJxW_g*rho*(IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp1]

                                                        );
#if DIMENSION==3  // --------------------------------------------------------------------------------------
          //--------------------------------- Block +2 [3-4-8] --------------------------------------------------------------
          const int idimp2= (ivar+2+_dir) % _nNSdim;
          KeM(indx,j+idimp2*el_ndof2) +=
           i_flag_j* dtxJxW_g*rho*(IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp2]

                         );
#endif    // ----------------------------------------------------------------------
        } // end A element matrix quad -quad (end loop on j)---------------------------------------------

        // ------------------------------------------------------------------
// #if FSI_EQUATIONS==1      // B^T element matrix ( p*div(v) )--------------------
        for(int  ikl=0; ikl<2; ikl++) {    // ikl=0 discontinuous ikl=1 continuous pressure
          for(int  ivarl=0; ivarl<_nvars[ikl]; ivarl++) {
            for(int  j=0; j<el_ndof[ikl]; j++) {
              const double psij_g= _phi_g[ikl][j];
              
              if(axysim==1)  KeM(indx,j+ _nNSdim*el_ndof2) -=dtxJxW_g*(1- (ivar + _dir)) *psij_g*phii_g/_ub_g[2][0];
              KeM(indx,j+ _nNSdim*el_ndof2) +=dtxJxW_g* (    // MPascal
                                             -_FSI_parameter.PINTP*psij_g*dphiidx_g2[ivar]
                                             +(1-_FSI_parameter.PINTP)*_dphi_g[ikl][j+ivar*el_ndof[1]]*phii_g
                                                +_dphi_g[ikl][j+ivar*el_ndof[1]]*Phi_supg
                                              );
            } // j
          } // ivarl
        }
// #endif  // end B^T element matrix ------------------------------------
      } // end loop ivar
    } // end loop i

//------------------    QL    -----------------------------------------------------------
// #if FSI_EQUATIONS==1   // only coupled Assemblying Matrix linear ------------------------
    for(int  ikl=0; ikl<2; ikl++) {    // ikl=0 discontinuous ikl=1 continuous pressure
      for(int   ivarl=0; ivarl< _nvars[ikl]; ivarl++) {
        for(int   i=0; i< el_ndof[ikl]; i++) {    // +++++++++++
          // set up row i
          int  indx=i+el_ndof2*_nvars[2];//ivar=0;
          const double psii_g=_phi_g[ikl][i];
//           int i_flag_i= (flag_group[i]>999)?0:1;
//              if (mode == 1)   FeM(indx)  +=  JxW_g2*psii_g*_ub_g[1][0]*KOMP_FSI; //    _u_old[DIMENSION]*KOMP_FSI;
          double dtxJxWp_g=JxW_g2*fabs(_bc_el[indx]);

          // linear-quadratic Assemblying Matrix -----------------------------
          for(int j=0; j<el_ndof2; j++)     {
            const double phij_g= _phi_g[2][j];
//               int i_flag_j=1;if(fabs(flag_group[j])>999) i_flag_j=0;
            for(int  jvar=0; jvar< _nvars[2]; jvar++) {    // linear -quad
                               // p-equation
              if(axysim==1) KeM(indx,j+jvar*el_ndof2) +=  dtxJxWp_g*rho* (1- (jvar + _dir)) *psii_g*phij_g/_ub_g[2][0];
              KeM(indx,j+jvar*el_ndof2) +=  /*i_flag_j**/dtxJxWp_g*rho* (
                                              psii_g*_dphi_g[2][j+jvar*el_ndof2]  // div=0
                                            );

            }// jvar
          }  // j end linear-quad --------------------------------------------

          // linear-linear Assemblying Matrix ------------------------------
//           for(int j=0; j<el_ndof[1]; j++)  {
//             const double psij_g=_phi_g[1][i];
//                 KeM(indx,j+ _nNSdim*el_ndof2)  += JxW_g2*psii_g*psij_g*0.0* KOMP_FSI;
//           } // end linear liner -----------------------------------------------------
        }  // i
      }// ivarl end linear +++++++++++
    }// ikl=0 discontinuous ikl=1 continuous pressure
// #endif  // -------------------- FSI_EQUATIONS==1 ----------------------------------------
  } // end of the quadrature point qp-loop

// ====================== end volume (element) =======================================
  return;
}






// ==============================================================================================
void MGSolFSI::matrixrhs_sol_vol(
  DenseMatrixM &KeM, DenseVectorM &FeM,
  int el_ndof[],
  double u_old[],  double u_oold[],   double u_nl[], double p_proj[], double dp_proj[],
  const int unsteady_flag, const int axysim,
  const double les, const int mode, int flag_group[]
) {

  double dphijdx_g2[DIMENSION],dphijdx_g1[DIMENSION];
  double dphiidx_g2[DIMENSION],dphiidx_g1[DIMENSION];
  double val_tbg[1];
  double vel_gddx[DIMENSION*DIMENSION*DIMENSION];
  double vel_gdx[DIMENSION*DIMENSION];
  double vela_g[DIMENSION],u_nlg[DIMENSION]; 
  double l_old[DIMENSION*NDOF_FEM];
 
  int el_ndofp=el_ndof[1];
  double rho=1.;
  int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim-1];                //elem gauss points
  const int el_ndof2=el_ndof[2];
  
  
  
  for(int qp=0; qp<  el_ngauss; qp++) {
    // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
    // quadratic continuous (2)  (velocity)
    const double det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);     // quadratic Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[ _nNSdim-1][qp];           // quadratic weight
    _fe[2]->get_phi_gl_g(_nNSdim,qp,_phi_g[2]);                     // quadratic shape function
    _fe[2]->get_dphi_gl_g(_nNSdim,qp,_InvJac2,_dphi_g[2]);          // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nNSdim,qp,_InvJac2,_ddphi_g[2]);        // global second derivatives
    // continuous linear (1) (linear pressure)
    const double det1      = _fe[1]->Jac(qp,_xx_qnds,_InvJac1);     // linear Jacobian
    double JxW_g1 =det1*_fe[1]->_weight1[ _nNSdim-1][qp];           // linear weight
    _fe[1]->get_phi_gl_g(_nNSdim,qp,_phi_g[1]);                     // linear shape funct
    _fe[1]->get_dphi_gl_g(_nNSdim,qp,_InvJac1,_dphi_g[1]);          // global coord deriv
    // discontinuous (0) (disc pressure)
    if(_nvars[0]>0) { _fe[0]->get_phi_gl_g(_nNSdim,qp,_phi_g[0]);}  // piecewise shape function

    interp_el_sol(_xx_qnds,0,_nNSdim,_phi_g[2],el_ndof2,_xyzg);
    // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------

    // quadratic fields (velocity)
    interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2]);    // field _ub_g[2][DIM]
    interp_el_gdx(u_old,0, _nNSdim,_dphi_g[2],el_ndof2,vel_gdx);                                           // derivatives vel_gdx[DIM][DIM]
    interp_el_gddx(u_old,0, _nNSdim,_ddphi_g[2],el_ndof2,vel_gddx);
    // linear fields (pressure) -> projection and segregated
// #if (FSI_EQUATIONS%2==0)
//     interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndofp,_ub_g[1]);
// #endif
    
    if(axysim==1) {JxW_g2  *=_xyzg[0]; JxW_g1  *=_xyzg[0];}  // std::cout<<_ub_g[2][0]<<"   "<<std::endl;

    // Velocity, Reynolds and upwind  --------------------------------------------------------------------------
    double mod2_vel=1.e-20;        _sP=1.e-20;    double mu=1.; // Velocity
    for(int idim=0; idim<  _nNSdim; idim++) {
      vela_g[idim]=0.;   u_nlg[idim]=0.; // old and  non linear Velocity at gaussina point qp
      for(int k=0; k< NDOF_FEM; k++) {
        u_nlg[idim] += u_nl[k+idim*NDOF_FEM]*_phi_g[2][k];
        vela_g[idim] += u_old[k+idim*NDOF_FEM]*_phi_g[2][k];
         
      }
      mod2_vel += u_nlg[idim]*u_nlg[idim]; // module velocity
      for(int jdim=0; jdim<  _nNSdim; jdim++) {
        _sP += 0.25*(vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim])*
               (vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim]);
      }
    }
    mod2_vel =sqrt(mod2_vel);
    // upwind term -> f_upwind
    // h_eff=2/sum|s.dN| and f_upwind ---------------------------------------------------------------------------------
    double h_eff=1.e-21;
    for(int i=0; i<el_ndof2; i++) {
      double hh=1.e-20; for(int idim=0; idim< _nNSdim; idim++) { hh += u_nlg[idim]*_dphi_g[2][i+idim*el_ndof2]/mod2_vel; }
      h_eff += fabs(hh);
    }
    h_eff=2./h_eff; if(h_eff<1.e-10) {h_eff=1. ; std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";   }

    double Pe_h=0.5*mod2_vel*h_eff/ _mus;
    double a_opt=(1./tanh(Pe_h)-1./Pe_h); if(a_opt >1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
    double f_upwind=rho*0.5*a_opt*h_eff/(mod2_vel);   // upwind

    // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
    if(_FF_idx[T_F]>=0) {
      double Temp_g=_ub_g[2][_FF_idx[T_F]];   rho *= _FSIdensity(Temp_g,0);  mu *=  _FSIkviscosity(Temp_g,0);
    }

    _mu_turb=0.;// eff visc turb at g pt
    // -------------------- Turbulence [K_F] -> (quad,_indx_eqs[K_F]) -------------------------------------
    if(_FF_idx[K_F]>=0) {// turbulence _mu_turb evaluation
      _kappa_g[0]= _ub_g[2][_FF_idx[K_F]];  _kappa_g[1]= _ub_g[2][_FF_idx[K_F]+1];
//       _mu_turb=eval_var1(val_tbg);
      _mu_turb= _kappa_g[0]/_kappa_g[1];
    }
    double IRe_eff  =_mus*mu+_mu_turb  ;   //  visc_eff  at g pt
    double lambda  =_lambda*_FSI_parameter.COMPRESSIBLE; 
//        if( fabs(_xyzg[1])> 0.2   )  IRe_eff *=3.;
    // ===========================================================================================================
    //                                       D) Assembling NS equation
    // ===========================================================================================================
    for(int i=0; i< el_ndof2; i++) {    // +++++++++++
      // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
         
      
      const double phii_g=_phi_g[2][i];
      // Advection and Laplacian operator (supg)
      for(int  ivar=0; ivar< _nvars[2]; ivar++)    {
        int    indx=i+ivar*el_ndof2; // ivar=0;
          int i_flag_i= (fabs(flag_group[i])>999)?0:1;
        double dtxJxW_g=JxW_g2*fabs(_bc_el[indx]);
        double Adv_expl=0.,Lap_expl=0.,supglapexpl=0.; double  Phi_supg=0.;
        for(int idim=0; idim< _nNSdim; idim++) {
          int ix=(ivar+_dir)*_nNSdim+idim;
          dphiidx_g2[idim]=_dphi_g[2][i+idim*el_ndof2];
          if(_bc_el[indx] >0 ) Phi_supg +=_FSI_parameter.SUPG*f_upwind*u_nlg[idim]*dphiidx_g2[idim]; //f_upwind
          Adv_expl    +=  vela_g[idim]*vel_gdx[ix];
          Lap_expl    +=  IRe_eff*vel_gdx[ix]*dphiidx_g2[idim];
          supglapexpl +=  IRe_eff*vel_gddx[ix*_nNSdim+idim];
        }
        // -------------------------------------- Assemblying rhs ------------------------------------------------------
//       if ( _bc_el[indx] != 0){ KeM(indx,indx)  =dtxJxW_g; FeM(indx) =1.;}
//        KeM(indx,indx)  +=1*fabs(_bc_el[indx])*i_flag_i; 
//       if((flag_group[i])<999 )   
//         KeM(indx,indx)  +=dtxJxW_g/_dt;
//      FeM(indx) +=dtxJxW_g/_dt*ivar;
//        if ( _bc_el[indx] == 0){ KeM(indx,indx)  = 1; FeM(indx) =2.;}
//       if(fabs(flag_group[i])>999 && _bc_el[indx] != 0)KeM(indx,indx)  +=  dtxJxW_g;
        if(mode == 1) {
        
          FeM(indx)  +=  dtxJxW_g* (
                           unsteady_flag* rho*vela_g[ivar+_dir]* (phii_g+Phi_supg) /_dt     // time
                           + rho*_IFr*_dirg[ivar+_dir]* (phii_g+Phi_supg)   // x-gravity
#ifdef Nat_Conv
                           + rho* (_ub_g[2][T_idx]-_T_nc) *_grav[ivar+_dir]*_beta_mat* (phii_g+Phi_supg)
#endif
                         );
        }

        //-------------------------------------- Assemblying matrix ---------------------------------------------------
//
        for(int j=0; j<el_ndof2; j++) {
          const double phij_g= _phi_g[2][j];
          // set up
          double Lap_g=0.,Adv_g=0., Div_g= 0., Supg_lap=0.;
          if(axysim==1) {  Lap_g =2.*(1-(ivar+_dir))*IRe_eff*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]); }      // axysimmetry
         
          for(int kdim=0; kdim< _nNSdim; kdim++) {
            dphijdx_g2[kdim] =_dphi_g[2][j+kdim*el_ndof2];
            Adv_g    += u_nlg[kdim]*dphijdx_g2[kdim]; // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
            Lap_g    += (IRe_eff+_FSI_parameter.UPWIND*f_upwind*u_nlg[kdim]*u_nlg[kdim]) *dphijdx_g2[kdim]*dphiidx_g2[kdim];
            Supg_lap += IRe_eff*_ddphi_g[2][j* _nNSdim* _nNSdim+kdim* _nNSdim+kdim];
          }
             
          
          //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
          KeM(indx,j+ivar*el_ndof2) +=dtxJxW_g*rho*(
                                        (Adv_g+ unsteady_flag*(1./_dt+200.)*phij_g)* (phii_g+Phi_supg) // time + advection
                                        + Lap_g*_dt                        // viscous Laplacian
                                        -Supg_lap*Phi_supg *_dt            // SUPG viscous Laplacian
                                         + IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir]*_dt  // viscous tensor
                                           +lambda*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir]*_dt
                                      );
           FeM(indx) -=dtxJxW_g*rho*  _data_eq[2].ub[_FF_idx[SDSX_F+ivar+_dir]*NDOF_FEM+j]*(
                                        + Lap_g                       // viscous Laplacian
                                        -Supg_lap*Phi_supg            // SUPG viscous Laplacian
                                         + IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir] // viscous tensor
                                         +lambda*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir]
                                      );
          // --------------------------------- Block +1 [2-6-7] --------------------------------------------------------------
          // out of diagonal viscous tensor
          const int idimp1= (ivar+1+_dir) % _nNSdim;
           KeM(indx,j+idimp1*el_ndof2) +=    dtxJxW_g*rho*_dt*(
             IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp1]+
           lambda*dphijdx_g2[idimp1+_dir]*dphiidx_g2[ivar+_dir]);
          FeM(indx)-=dtxJxW_g*rho* _data_eq[2].ub[_FF_idx[SDSX_F+idimp1]*NDOF_FEM+j]*(
            IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp1]+
            lambda*dphijdx_g2[idimp1+_dir]*dphiidx_g2[ivar+_dir]
          );
#if DIMENSION==3  // --------------------------------------------------------------------------------------
          //--------------------------------- Block +2 [3-4-8] --------------------------------------------------------------
          const int idimp2= (ivar+2+_dir) % _nNSdim;
          KeM(indx,j+idimp2*el_ndof2) += dtxJxW_g*rho*(
            IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp2]+
            lambda*dphijdx_g2[idimp2+_dir]*dphiidx_g2[ivar+_dir] 
          );
          
          FeM(indx) -=dtxJxW_g*rho* _data_eq[2].ub[_FF_idx[SDSX_F+idimp2*NDOF_FEM+j]]*(
            IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp2]+
            lambda*dphijdx_g2[idimp2+_dir]*dphiidx_g2[ivar+_dir]
            
          );
#endif    // ----------------------------------------------------------------------
        } // end A element matrix quad -quad (end loop on j)---------------------------------------------

        // ------------------------------------------------------------------
// #if FSI_EQUATIONS==1      // B^T element matrix ( p*div(v) )--------------------
        for(int  ikl=0; ikl<2; ikl++) {    // ikl=0 discontinuous ikl=1 continuous pressure
          for(int  ivarl=0; ivarl<_nvars[ikl]; ivarl++) {
            for(int  j=0; j<el_ndof[ikl]; j++) {
              const double psij_g=_phi_g[ikl][j];
              
              if(axysim==1)  KeM(indx,j+ _nNSdim*el_ndof2) -=dtxJxW_g*(1- (ivar + _dir)) *psij_g*phii_g/_ub_g[2][0];
              KeM(indx,j+ _nNSdim*el_ndof2) +=dtxJxW_g*(1-_FSI_parameter.COMPRESSIBLE)*(    // MPascal
                                            -_FSI_parameter.PINTP*psij_g*dphiidx_g2[ivar]
                                             +(1-_FSI_parameter.PINTP)*_dphi_g[ikl][j+ivar*el_ndof[1]]*phii_g
                                                +_dphi_g[ikl][j+ivar*el_ndof[1]]*Phi_supg
                                              );
            } // j
          } // ivarl
        }
// #endif  // end B^T element matrix ------------------------------------
      } // end loop ivar
    } // end loop i

//------------------    QL    -----------------------------------------------------------
// #if FSI_EQUATIONS==1   // only coupled Assemblying Matrix linear ------------------------
    for(int  ikl=0; ikl<2; ikl++) {    // ikl=0 discontinuous ikl=1 continuous pressure
      for(int   ivarl=0; ivarl< _nvars[ikl]; ivarl++) {
        for(int   i=0; i< el_ndof[ikl]; i++) {    // +++++++++++
          // set up row i
          int  indx=i+el_ndof2*_nvars[2];//ivar=0;
          const double psii_g=_phi_g[ikl][i];
//              if (mode == 1)   FeM(indx)  +=  JxW_g2*psii_g*_ub_g[1][0]*KOMP_FSI; //    _u_old[DIMENSION]*KOMP_FSI;
           int i_flag_i= (fabs(flag_group[i])>999 && ikl>0 )?0:1;
          double dtxJxWp_g=JxW_g2*fabs(_bc_el[indx])*(1-_FSI_parameter.COMPRESSIBLE);
            
           KeM(indx,indx) +=JxW_g2*_FSI_parameter.COMPRESSIBLE*i_flag_i;
          // linear-quadratic Assemblying Matrix -----------------------------
          for(int j=0; j<el_ndof2; j++)     {
            const double phij_g= _phi_g[2][j];
            for(int  jvar=0; jvar< _nvars[2]; jvar++) {    // linear -quad
                               // p-equation
              if(axysim==1) KeM(indx,j+jvar*el_ndof2) +=  dtxJxWp_g*rho* (1- (jvar + _dir)) *psii_g*phij_g/_ub_g[2][0];
              KeM(indx,j+jvar*el_ndof2) +=  dtxJxWp_g*rho* (
                                              psii_g*_dphi_g[2][j+jvar*el_ndof2]  // div=0
                                            );

            }// jvar
          }  // j end linear-quad --------------------------------------------

          // linear-linear Assemblying Matrix ------------------------------
//           for(int j=0; j<el_ndof[1]; j++)  {
//             const double psij_g=_phi_g[1][i];
//                 KeM(indx,j+ _nNSdim*el_ndof2)  += JxW_g2*psii_g*psij_g*0.0* KOMP_FSI;
//           } // end linear liner -----------------------------------------------------
        }  // i
      }// ivarl end linear +++++++++++
    }// ikl=0 discontinuous ikl=1 continuous pressure
// #endif  // -------------------- FSI_EQUATIONS==1 ----------------------------------------
  } // end of the quadrature point qp-loop

// ====================== end volume (element) =======================================
  return;
}



// #ifdef TBK_EQUATIONS
// =====================================================
double  MGSolFSI::eval_var1(double val[]) {
  // Turbulent viscosity

  double kappa = _kappa_g[0];    double omega = _kappa_g[1];

#ifdef LOG_W
  omega = exp(_kappa_g[1]);
#endif

#ifdef LOG_K
  kappa = exp(_kappa_g[0]);
#endif

//    if(omega < 1.e-10) { omega= 1.e-10; }  // epsilon/omega
//    if(kappa< 1.e-10) { kappa= 1.e-10; }  // kappa
  double tau_k=1.;  // turbulent time    double tau_k=1.;  // turbulent time

#if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
  tau_k=_y_bcout/sqrt(kappa);
#endif  // end kappa -------------------------------------------------------

#if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ----------
  tau_k= CMU*kappa/omega; // tau=1/omega=CMU*kappa/epsilon
#endif   // kappa-epsilon --------------------------------------------------

#if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
  tau_k=1/omega; // tau=1/omega

#endif  // -----------------------------------------------------------------
// turbulent viscosity
  double Re_t= kappa*tau_k/_IRe;
  double mu_turb= fabs(Re_t);
//
  /* ================================== Nagano Model ==================================== */

#ifdef NAGANO /// Nagano model with low-Re correction, valid for k-e and k-w

  const int _NS_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[FS_F]];
  double vel = _data_eq[2].ub[(_NS_idx+AX_DIR)*NDOF_FEM+(NDOF_FEM-1)];
  double tildey = sqrt(_y_bcout*vel/_IRe);


  const double R_t = Re_t/CMU;                                             // turbulent Reynolds number
  const double R_d = _y_bcout*sqrt(kappa/sqrt(R_t))/_IRe;                   // non-dimensional Kolmogorov-based distance from the wall
  const double f_mu = (1.-exp(-1.*R_d/14.))*(1.-exp(-1.*R_d/14.));
  const double f_corr = 1. + 5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.);

  mu_turb *= f_corr*f_mu;    // Nagano low Re correction

// #ifdef WALL_FUNC_APP
//   double yp = _y_bcout * sqrt(_dirg[AX_DIR]*0.03025*9.81) / _IRe;
//   if(yp < 5.) { mu_turb *= f_mu; }
// #else
//   mu_turb *= f_mu;
// #endif

#endif

  /* ================================== SST Model ======================================= */

#ifdef SST

  // F1 and F2 coeff calculation
  double F1,F2;
  double F1_first  = 500.*_IRe*tau_k/(_y_bcout*_y_bcout);
  double F2_first= 2.*sqrt(kappa)*tau_k/(BETASTAR*_y_bcout);
  if(F1_first > F2_first) { F2_first=F1_first; }
  F2=tanh(F2_first*F2_first);

  // nu_t calculation
//       double alpha_star        = 1.;//(0.024+ Re_t/6.)/(1.+ Re_t/6.);
  double alpha_starb = sqrt(_sP)*F2*tau_k/0.31;
  if(alpha_starb < 1.) {
    alpha_starb=1.; /*printf("alpha star is %f and and %f and %f\n",alpha_starb, 1./tau_k, F2);*/
  }
  mu_turb /= alpha_starb;

//  mu_turb = Re_t;
  if(mu_turb > MU_TOP) { mu_turb =MU_TOP; }
  if(mu_turb < MU_LOW) { mu_turb =MU_LOW; }
//       printf("mu_turb is %f \n",mu_turb);

#endif
  return mu_turb;
}
/******************************************************************************************************/
//
//


#endif
#endif


