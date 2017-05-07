// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
#if NS_EQUATIONS==2
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class


// Thermodinamical Properties ==========================================================================================
// constant
double  _NSdensity(double T,double T_r) {return 1.;} // water (t K) T_ref
double  _NSkviscosity(double T,double T_r) {return 1.;} // lead T K
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

void MGSolNS::get_el_field_data(
  int iel, int Level,
  int el_conn [], int offset,int el_ndof[],int ndof_lev,
  double u_old[],  double u_oold[],   double u_nl[],
  double p_proj[], double dp_proj[]
) {

  int el_ndofp=el_ndof[1];

  // element nodes coordinates ----------------------------------------------------------------------------------------
  for (int idim=0; idim<DIMENSION; idim++) for (int d=0; d< NDOF_FEM; d++) {
      _data_eq[2].ub[idim*NDOF_FEM+d]=_xx_qnds[idim*NDOF_FEM+d];
    }
  // external fields (from constant 0 to quadratic 2) -------------------------------------------------------------------
  for (int deg=0; deg<3; deg++)
    for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
      _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                           el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
    }

  //  internal field data (NS) -------------------------------------------------------------------------------------------
// pressure as external field (splitting) -----------------------------------------------------------
  for (int kdim=0; kdim< _nNSdim; kdim++) {
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F+kdim]]->get_el_nonl_sol(0,1,el_ndof[2],el_conn,offset,kdim,u_nl);
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F+kdim]]->get_el_sol(0,1,el_ndof[2],el_conn, offset,kdim,u_old);
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F+kdim]]->get_el_oldsol(0,1,el_ndof[2],el_conn, offset,kdim,u_oold);
  }
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0,1,el_ndofp,el_conn, offset,0,p_proj);       // old pressure
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndofp,el_conn, offset,0,dp_proj);   // dp pressure

  return;


}






// ==============================================================================================
void MGSolNS::matrixrhsvol(DenseMatrixM &KeM, DenseVectorM &FeM,
                           int el_ndof[],const int mode,const double euler_impl,const double h_eff,
                          
                           double u_old[],  double u_oold[],   double u_nl[], double p_proj[], double dp_proj[]
                          ) {


// ==================================  Volume ===============================================

  double det2,det1,JxW_g2,JxW_g1;           // Jac, Jac*w Jacobean
  double dphijdx_g2[DIMENSION],dphijdx_g1[DIMENSION];
  double dphiidx_g2[DIMENSION],dphiidx_g1[DIMENSION];
  double val_tbg[1];
  double vel_gddx[DIMENSION*DIMENSION*DIMENSION];
  double vel_gdx[DIMENSION*DIMENSION];
  double vel_g[DIMENSION],u_nlg[DIMENSION];
  int el_ndofp=el_ndof[1];
//     euler_impl=1;
  double rho=1.;
  int  el_ngauss = _fe[2]->_NoGauss1[ _nNSdim-1];                //elem gauss points

// ==================================  Volume ===============================================
  int el_ndof2=el_ndof[2];
  double sumk=0;
  for (int qp=0; qp<  el_ngauss; qp++) {
    // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
    // quadratic continuous (2)  (velocity)
    const double det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);         // quadratic Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[ _nNSdim-1][qp];             // quadratic weight
    _fe[2]->get_phi_gl_g(_nNSdim,qp,_phi_g[2]);                     // quadratic shape function
    _fe[2]->get_dphi_gl_g(_nNSdim,qp,_InvJac2,_dphi_g[2]);             // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nNSdim,qp,_InvJac2,_ddphi_g[2]);         // global second derivatives
    // continuous linear (1) (linear pressure)
    const double det1      = _fe[1]->Jac(qp,_xx_qnds,_InvJac1);       // linear Jacobian
    double JxW_g1 =det1*_fe[1]->_weight1[ _nNSdim-1][qp];           // linear weight
    _fe[1]->get_phi_gl_g(_nNSdim,qp,_phi_g[1]);                     // linear shape funct
    _fe[1]->get_dphi_gl_g(_nNSdim,qp,_InvJac1,_dphi_g[1]);           // global coord deriv
    // discontinuous (0) (disc pressure)
    if (_nvars[0]>0) { _fe[0]->get_phi_gl_g(_nNSdim,qp,_phi_g[0]);  }       // piecewise shape function


    // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------
    // quadratic fields (velocity)
    interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2]);    // field _ub_g[2][DIM]
    interp_el_gdx(u_old,0, _nNSdim,_dphi_g[2],el_ndof2,vel_gdx);    // derivatives  vel_gdx[DIM][DIM]
    interp_el_gddx(u_old,0, _nNSdim,_ddphi_g[2],el_ndof2,vel_gddx);
    // linear fields (pressure) -> projection and segregated
    interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndofp,_ub_g[1]);
    
#ifdef AXISYM
    JxW_g2  *=_ub_g[2][0]; JxW_g1  *=_ub_g[2][0];  // std::cout<<_ub_g[2][0]<<"   "<<std::endl;
#endif

    // Velocity, Reynolds and upwind  --------------------------------------------------------------------------
    double mod2_vel=1.e-20;        _sP=1.e-20;    double mu=1.; // Velocity
    for (int idim=0; idim<  _nNSdim; idim++) {
      vel_g[idim]=0.;   u_nlg[idim]=0.; // old and  non linear Velocity at gaussina point qp
      for (int k=0; k< NDOF_FEM; k++) {
        u_nlg[idim] += u_nl[k+idim*NDOF_FEM]*_phi_g[2][k];
        vel_g[idim] += u_old[k+idim*NDOF_FEM]*_phi_g[2][k];
      }
      mod2_vel += u_nlg[idim]*u_nlg[idim]; // module velocity
         for(int jdim=0; jdim<  _nNSdim; jdim++) {  
         _sP += 0.25*(vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim] )*
                                         (vel_gdx[ jdim*_nNSdim+idim] +vel_gdx[ idim*_nNSdim+jdim] );
	 }
    }
    mod2_vel =sqrt(mod2_vel);
    // upwind term -> f_upwind
    double Pe_h=0.5*mod2_vel*h_eff/_IRe;
    double f_upwind=NS_parameter.NS_UPWIND*rho*0.25* (1./tanh(Pe_h)-1./Pe_h)*h_eff/(mod2_vel);      // SUPG

    // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
    if (_T_idx>=0) { double Temp_g=_ub_g[2][_T_idx];   rho *= _NSdensity(Temp_g,0);  mu *=  _NSkviscosity(Temp_g,0);}
//     double IRe_eff=_IRe*mu;
    // -------------------- Turbulence [K_F] -> (quad,_indx_eqs[K_F]) -------------------------------------
    _mu_turb=0.;// turbulence _mu_turb evaluation/
    if (_K_idx>=0) { _kappa_g[0]= _ub_g[2][_K_idx];  _kappa_g[1]= _ub_g[2][_K_idx+1]; _mu_turb=eval_var1(val_tbg); }
//     IRe_eff  *= (1.+_mu_turb);   //  visc_eff  at g pt
double IRe_eff=_IRe*mu* (1.+_mu_turb)+  _IRe*  NS_parameter.NS_LES*0.01*h_eff*h_eff*sqrt(2*_sP);
    // ===========================================================================================================
    //                                       D) Assembling NS equation
    // ===========================================================================================================
    for (int i=0; i< el_ndof2; i++) {   // +++++++++++
      // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
      const double phii_g=_phi_g[2][i];
      double  Phi_supg=0.;
      for (int idim=0; idim< _nNSdim; idim++) {
        dphiidx_g2[idim]=_dphi_g[2][i+idim*el_ndof2];
        Phi_supg += NS_parameter.NS_SUPG*_dt*u_nlg[idim]*dphiidx_g2[idim];
      }
      // Advection and Laplacian operator (supg)
      for (int  ivar=0; ivar< _nvars[2]; ivar++)    {
        int    indx=i+ivar*el_ndof2; // ivar=0;
        double dtxJxW_g=JxW_g2*_bc_el[indx];
        double Adv_expl=0.,Lap_expl=0.,supglapexpl=0.;
        for (int idim=0; idim< _nNSdim; idim++) {
          int ix=(ivar+_dir)*_nNSdim+idim;
          Adv_expl    +=   vel_g[idim]*vel_gdx[ix];
          Lap_expl    +=   IRe_eff*vel_gdx[ix]*dphiidx_g2[idim];
          supglapexpl +=  IRe_eff*vel_gddx[ix*_nNSdim+idim];
        }

        // -------------------------------------- Assemblying rhs ------------------------------------------------------
        if (mode == 1) {
          FeM(indx)  +=  dtxJxW_g* (
                           rho*vel_g[ivar+_dir]* (phii_g+Phi_supg) /_dt     // time
                           + rho*_IFr*_dirg[ivar+_dir]* (phii_g+Phi_supg)   // x-gravity
//                                 -1.*(sum)*dphiidx_g2[ivar+_dir]*penalty_f2      // rhs penalty
// #ifdef Nat_Conv
//                                       + rho* ( _ub_g[2][T_idx]-_T_nc ) *_grav[ivar+_dir]*_beta_mat* ( phii_g+Phi_supg )
// #endif
                         );
        }
        //-------------------------------------- Assemblying matrix ---------------------------------------------------
        for (int j=0; j<el_ndof2; j++) {
          const double phij_g= _phi_g[2][j];
          // set up
          double Lap_g=0.,Adv_g=0., Div_g= 0., Supg_lap=0.;
#ifdef  AXISYM
          Lap_g =2.*(1-(ivar+_dir))*IRe_eff*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);       // axysimmetry
#endif
          for (int kdim=0; kdim< _nNSdim; kdim++) {
            dphijdx_g2[kdim] =_dphi_g[2][j+kdim*el_ndof2];
            Adv_g    += u_nlg[kdim]*dphijdx_g2[kdim]; // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
            Lap_g    += (IRe_eff+f_upwind*u_nlg[kdim]*u_nlg[kdim]) *dphijdx_g2[kdim]*dphiidx_g2[kdim];
            Supg_lap += IRe_eff*_ddphi_g[2][j* _nNSdim* _nNSdim+kdim* _nNSdim+kdim];
          }
          //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
          // ===============================> Termine di pressione per il metodo proiettivo <==============================
          if (j<el_ndofp) {
            FeM(indx) += dtxJxW_g* (
                 -(p_proj[j]+dp_proj[j])*_dphi_g[1][j+(ivar+_dir)*el_ndof[1]]*phii_g
                 -(p_proj[j]+dp_proj[j])*_dphi_g[1][j+(ivar+_dir)*el_ndof[1]]*Phi_supg
//                            (p_proj[j]+dp_proj[j])*_phi_g[1][j]*dphiidx_g2[ivar+_dir]
//                             -(p_proj[j]+dp_proj[j])*_dphi_g[1][j+(ivar+_dir)*el_ndof[1]]*Phi_supg // ?????????????
                         );
          }

          KeM(indx,j+ivar*el_ndof2) +=dtxJxW_g*rho* (
                                       (Adv_g+ phij_g /_dt)* (phii_g+Phi_supg)            // time + advection
//                                         + euler_impl*Adv_g*phii_g                 //
                                        + /* euler_impl* */Lap_g                        // viscous Laplacian
//                                         + euler_impl*Adv_g*Phi_supg               // SUPG advection
                                        - /* euler_impl* */Supg_lap*Phi_supg            // SUPG viscous Laplacian
                                        + /*euler_impl* */IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[ivar+_dir] // viscous tensor
//                                             + impl_pen_f*penalty_f2*dphijdx_g2[ivar]*dphiidx_g2[ivar+_dir]
//  *(IRe_eff+0.*f_upwind[ivar+_dir]*u_nlg[ivar+_dir]*u_nlg[ivar+_dir]
                                      );
//
          // --------------------------------- Block +1 [2-6-7] --------------------------------------------------------------
//                        // out of diagonal viscous tensor
          const int idimp1= (ivar+1+_dir) % _nNSdim;
          FeM(indx) += -1.*u_nl[j+idimp1*el_ndof2]* dtxJxW_g*rho/* *euler_impl */* (IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp1]
                       //      +impl_pen_f*penalty_f2*dphijdx_g2[idimp1]*dphiidx_g2[ivar+_dir]
                                                                             );
#if DIMENSION==3  // --------------------------------------------------------------------------------------
          //--------------------------------- Block +2 [3-4-8] --------------------------------------------------------------
          const int idimp2= (ivar+2+_dir) % _nNSdim;
          FeM(indx) +=  -1.*u_nl[j+idimp2*el_ndof2]* dtxJxW_g*rho*/* euler_impl*  */(IRe_eff*dphijdx_g2[ivar+_dir]*dphiidx_g2[idimp2]
                        //       +impl_pen_f*penalty_f2*dphijdx_g2[idimp2]*dphiidx_g2[ivar+_dir]
                                                                              );
#endif    // ----------------------------------------------------------------------
        } // end A element matrix quad -quad (end loop on j)---------------------------------------------
        // ------------------------------------------------------------------
      } // end loop ivar
    } // end loop i
//------------------    QL    -----------------------------------------------------------
  } // end of the quadrature point qp-loop

// ====================== end volume (element) =======================================
  return;
}




// #ifdef TBK_EQUATIONS
// =====================================================
double  MGSolNS::eval_var1(double val[]) {
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
  double mu_turb= Re_t;
//
  /* ================================== Nagano Model ==================================== */

#ifdef NAGANO /// Nagano model with low-Re correction, valid for k-e and k-w

  const int _NS_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];
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

