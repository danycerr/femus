#include "Equations_conf.h"

// ============================================
#ifdef COLOR_EQUATIONS // 3D-2D Energy equation
// ============================================



#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverCOL.h"       // Navier-Stokes class header file

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MGFE.h"          // Mesh class,double vel_g[]

#ifdef   TWO_PHASE
 #include "MGSolverCC.h" 
#endif
// ===============================================================================================
void MGSolCOL::vol_integral(
  DenseMatrixM &KeM, DenseVectorM &FeM,
  const int el_ndof2,const int el_ngauss,
  double xx_qnds[],
  const int unsteady,const int mode
) { // ==============================================================================================

  double vel_g[DIMENSION];
  double C_der[DIMENSION],T_secder[DIMENSION*DIMENSION];
  double xyz_g[DIMENSION];
  double rhocp=1.; int block=CCLEV;
// --------------------------------------------------------------------------------------------------------------------
  /// c) gaussian integration loop (n_gauss)
  // ------------------------------------------------------------------------------------------------------------------
  for(int qp=0; qp< el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2      = _fe[2]->Jac(qp,xx_qnds,_InvJac2);     // Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[_nTdim-1][qp];       // weight
    _fe[2]->get_phi_gl_g(_nTdim,qp,_phi_g[2]);               // shape funct
    _fe[2]->get_dphi_gl_g(_nTdim,qp,_InvJac2,_dphi_g[2]); // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nTdim,qp,_InvJac2,_ddphi_g[2]); // local second deriv

    //  fields --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds,0,_nTdim,_phi_g[2],el_ndof2,xyz_g);
    interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2]); // quadratic
    //  derivatives
    interp_el_gdx(_data_eq[2].ub,_FF_idx[CO_F],1,_dphi_g[2],el_ndof2,C_der);
    interp_el_gddx(_data_eq[2].ub,_FF_idx[CO_F],1,_ddphi_g[2],el_ndof2,T_secder);
    
    
// //     color function 
//      _msolcc->get_2surten(xyz_g,ff,ord);
    double cc=0; double cc_st=0;
    if(_dir==0)
    {
   cc=_msolcc->get_2phase(block,xyz_g); // new cc
      if (cc!=cc || cc < 1.e-15)   cc=0.;
    }
    else{
       cc_st=_msolcc->get_2phase(block,xyz_g); // new cc
      if (cc_st!=cc_st || cc_st < 1.e-3 || cc_st > .5-1.e-3  )   cc_st=0.;
     cc= _ub_g[2][_FF_idx[CO_F]];
    }
//     double cc=_msolcc->get_2phase(block,xyz_g); // new cc
#ifdef AXISYM   // axisymmetric (index -> 0)
    JxW_g[2]  *=xyz_g[0];
#endif

    // h_eff=2/sum|s.dN| and f_upwind ---------------------------------------------------------------------------------
//     double h_eff=1.e-21;
//     for(int i=0; i<el_ndof2; i++) {
//       double hh=1.e-20; for(int idim=0; idim< _nTdim; idim++) hh += vel_g[idim]*_dphi_g[2][i+idim*el_ndof2]/mod2_vel;
//       h_eff += fabs(hh);
//     }
//     h_eff=2./h_eff; if(h_eff<1.e-10) {h_eff=1. ; std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";   }
// 
//     double Pe_h=0.5*mod2_vel*h_eff/(_IRe*_IPrdl);
//     double a_opt=(1./tanh(Pe_h)-1./Pe_h); if(a_opt >1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
//     double f_upwind=rhocp*0.5*a_opt*h_eff/(mod2_vel);   // upwind

    /// d) Local (element) assemblying energy equation
    // =====================================================================================================================
    for(int i=0; i<el_ndof2; i++)     {
 
      double dtxJxW_g=JxW_g2*_bc_el[i];           // area with bc and weight
      double Phi_supg=0.,Lap_expl=0.,Lap_supg=0.,Adv_expl=0.;  // supg, Lap explicit , Lap explicit supg
//       for(int idim=0; idim< _nTdim; idim++) {
//         const double  dphiidxg=_dphi_g[2][i+idim*el_ndof2];
//         Adv_expl += vel_g[idim]*C_der[idim];
//         Lap_supg += alpha_eff*T_secder[idim*_nTdim+idim];       // explicit Laplacian supg
//         Lap_expl += alpha_eff*C_der[idim]* dphiidxg;            // explicit Laplacian
// //         Phi_supg += _C_parameter.SUPG*_bc_el[i]* f_upwind*vel_g[idim]* dphiidxg; // phii_g+
//       }

      const double phii_g=_phi_g[2][i];

      // Rhs Assemblying  ---------------------------------------------------------------------------------------------------
      if(mode == 1) { // rhs
        if(_dir==0){ FeM(i) += dtxJxW_g*(cc*(phii_g+Phi_supg)  );}
        else{
          double sum=1.e-5;
         for(int kdim=0; kdim<_nTdim; kdim++ ) sum +=C_der[kdim]*C_der[kdim];
         for(int kdim=0; kdim<_nTdim; kdim++ ) FeM(i) += dtxJxW_g*(cc_st*C_der[kdim]*_dphi_g[2][i+kdim*el_ndof2]/sqrt(sum));  //    ((cc<0.001 || cc>0.2 -0.001 )? 1.:1)*
        }
      }
      
      // Matrix Assemblying ------------------------------------------------------------------------------------------------
      for(int j=0; j<el_ndof2; j++) {
        double Adv=0.,Lap=0.,Lap_supgi=0.;
        for(int idim=0; idim<  _nTdim; idim++) {
          const double  dphiidxg=_dphi_g[2][i+idim*el_ndof2];
          const double  dphijdxg=_dphi_g[2][j+idim*el_ndof2];
          Lap += dphijdxg*dphiidxg;                                                        // diffusion
        }

        const double phij_g= _phi_g[2][j];
        KeM(i,j) +=dtxJxW_g*( // energy-equation matrix
                     phij_g*(phii_g)// time term
                     +Lap*(_dir*0.e-0+(1-_dir)*1.e-9 ) //diffusion term
                   );
//         std::cout << "Matrix \n"  <<KeM  << "FeM\n"  <<FeM;
      }
    } // ----------------------------------------
  } // end of the quadrature point qp-loop ***********************

  return;
}

#endif
