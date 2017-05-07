// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
// ==============================================================
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#ifdef NS_EQUATIONS
#if NS_EQUATIONS==2

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#define PRESS_0 (2.)
#include "MGSolverNS.h"       // Navier-Stokes class header file

// local Femus class include -----------------------------------
#include "MGFE.h"          // Mesh class



 //  ============================================================================================
 //  ======     set the boundary condition (bc row in KeM and FeM)
 //  ============================================================================================
 /// This function sets the boundary condition on the local element matrix KeM and rhs FeM
 //  ============================================================================================
 //  Boundary conditions into local element matrix (Diriclet normal=0; Neumann=1) 
 //  Diriclet<10 (enum bound_cond):
 //                 simm=0,velocity_in0=1,velocity_tg0=2,wall=3,
 //                 simm1=4,velocity_in=5,velocity_tg=6,velocity=7;
 //  Neumann >= 10 (enum bound_cond):
 //                 outflow=10,pressure_outlet=12,interior=11,wall_turb=13,
 //                 outflow_p=14,pressure_inlet=16;
 //  =============================================================================================
void MGSolNS::set_bc_matrix(
  DenseMatrixM &KeM, DenseVectorM &FeM,
  int dir_maxnormal,int sur_toply[], 
  int el_ndof[],  int elb_ndof[],  int elb_ngauss, 
  double u_old[], double normal[], double Ipenalty
  ) {



  double eff_stress=1.; double pressure;
  double vel_bound[1];

  //  ============================================================================================
  //  ====================================== Neumann  ============================================
  //  ============================================================================================
  for(int  ivar=0; ivar< _nvars[2]; ivar++)    {
    int bc_dn=0;  if(_bc_vol[sur_toply[NDOF_FEMB-1]]/10>0) { bc_dn=1; }
    if(bc_dn>0) {   // Neumann  -------------------------------------------------------------------------------------------------------------------
      int bc_s=_bc_vol[sur_toply[NDOF_FEMB-1]]%10;
      int bc_tg=((bc_s&2)>>1);   // (0;?1?) tg
      int bc_normal=(bc_s%2);    // (0;??1) normal
      int bc_rhs=((bc_s&4)>>2);  // (0;1??) nonhomogeneous
      //              int bc_yplus=((bc_s&8)>>3);  // (0;1??) nonhomogeneous
      int flag_normal=(ivar+_dir == dir_maxnormal)?1:0 ;
      if(bc_s == 3) { flag_normal= !flag_normal; } // turb condition

     
      if(flag_normal)  {  // normal ---------------------------------------------------------------------------------------------------
        for(int  qp=0; qp< elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)
          // quad/linear  [2]=quad [1]=linear------------------------------------
          double   det     = _fe[2]->JacSur(qp,_xxb_qnds, _InvJac2);   // local coord _phi_g and jac
          double  JxW_g=det*_fe[2]->_weight1[_nNSdim-2][qp]; // weight
          _fe[2]->get_phi_gl_g(_nNSdim-1,qp,_phi_g[2]);   // global coord _phi_g
          _fe[1]->get_phi_gl_g(_nNSdim-1,qp,_phi_g[1]);   // global coord _phi_g
//                interp_el_sol(pressure,0,1,_phi_g[1],elb_ndof,_ub_g[1]);
#ifdef AXISYM   // axisymmetric  (index ->0)
          interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
          JxW_g  *=_ub_g[2][0];
#endif
          // Assemblying NS equation
          for(int i=0; i< elb_ndof[2]; i++) {  // Assemblying NS equation
            // set up row i
            const double phii_g=_phi_g[2][i];
            const int   indx_var=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
            const int   indx_row=sur_toply[i]+(ivar)*el_ndof[2];// volume dof index
            // boundary flag
           
            //                 int bc_s1=(int)_bc_bd[indx_row];
            //                  bc_rhs   =((bc_s1&4)>>2); // (1??) -> nonhomogeneous
            //                  bc_tg    =((bc_s1&2)>>1); // (?1?) -> tg
            //                  bc_normal=(bc_s1%2);      // (??1) -> normal
            int bc_v1=(int)(_bc_vol[sur_toply[i]]/10);
            double dtJxW_g=JxW_g*bc_v1;
            //  if(mode == 1)   {
//             FeM(indx_row)  += -1.*bc_rhs*dtJxW_g*phii_g*(normal[ivar+_dir])*20;//(_ub_g[1][0]);
            //    }
            // wall_turb boundary conditions
            double un=0;  for(int  kdim=0; kdim<_nNSdim; kdim++) { un +=u_old[sur_toply[i]+kdim*NDOF_FEM]*normal[kdim]; }
            double ut=0;  for(int  jvar=0; jvar<_nNSdim; jvar++) {
              const double tmp=u_old[sur_toply[i]+jvar*NDOF_FEM]-normal[jvar]*un; ut +=tmp*tmp;
            }
            vel_bound[0] =sqrt(ut);
            double u_tau= eval_var2(vel_bound);
#ifdef TBK_EQUATIONS
            eff_stress =u_tau*u_tau/fabs(vel_bound[0]+1.e-20);
#endif
            //                 printf("\t\t\t\t\t\t\t log region vel_bound %f u_tau %f y1_plus %f  Re_tau %f \n", vel_bound[i], u_tau[i],  y1_plus[i], u_tau[i]*0.03025/_IRe);

            for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad ------------------------------------------------------------------------------------------
              const double phij_g= _phi_g[2][j];
              KeM(indx_row,indx_row) +=   dtJxW_g*bc_tg*bc_normal*eff_stress*phij_g*phii_g
                                          // (?1?) -> tg     /* applied when bc_flag = 2, bc_Neum = 1*/
                                          // + dtJxW_g*eff_stress*normal[ivar+_dir]*normal[ivar+_dir]*phij_g*phii_g*(
                                          //  bc_normal // (??1) -> normal
                                          //  -bc_tg )    // (?1?) -> tg
                                          ;
              for(int  jvar=ivar+_dir+1; jvar< _nNSdim+_dir; jvar++)    {  // u
                int jvar1=jvar%DIMENSION;
                FeM(indx_row) +=   -1.*_data_eq[2].ub[_NS_idx*NDOF_FEM+sur_toply[j]+jvar1*el_ndof[2]]*
                                        dtJxW_g*eff_stress*phij_g*phii_g*bc_normal*bc_tg;
              }// jvar-loop
            }//  j-loop ---------------------------------------------------------------------------------------------------------

          }// i-loop
        }// end gaussian  integration
      }//  end  normal direction
      //  -------------------------------------------------------------------------------------------------------------------
      else  {  // Non normal direction  -------------------------------------------------------------------------------
        for(int lbnode=0; lbnode< elb_ndof[2]; lbnode++) {  // Assemblying NS equation
          // set up row                int  iflag=0;
          bc_s=(_bc_vol[sur_toply[lbnode]]%10)%4;
          bc_tg=((bc_s&2)>>1);   // (0;?1?) tg
          bc_normal=(bc_s%2);    // (0;??1) normal
          bc_rhs=((bc_s&4)>>2);  // (0;1??) nonhomogeneous
//                 if((_bc_vol[sur_toply[lbnode]]%10)%4 == 3){iflag=1;_bc_in[sur_toply[lbnode]+ivar*NDOF_FEM]= 1;}

          const int  indx_sol=sur_toply[lbnode]+(ivar+_dir)*el_ndof[2];//ivar=0;
          const int  indx_row=sur_toply[lbnode]+ivar*el_ndof[2];//ivar=0;
//             _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= ((_bc_vol[sur_toply[lbnode]]%10)%4 == 0 ) ? 1:0 ; // storage bc for vol part
          _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= ((bc_tg+bc_normal) == 0) ? 1:0 ;  // storage bc for vol part
//                  _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= ((_bc_vol[sur_toply[lbnode]]%10) == 1 ) ? 1:0 ; // storage bc for vol part
          //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
          // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
          //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
          // diagonal part
//           Ipenalty=1.;
          KeM(indx_row,indx_row) += ((bc_tg+bc_normal)%2)*Ipenalty;
          double aii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*((bc_tg+bc_normal));
          KeM(indx_row,sur_toply[lbnode]+ivar*el_ndof[2]) += -Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(bc_tg+bc_normal);
          FeM(indx_row) +=  bc_rhs*bc_tg*Ipenalty*u_old[indx_sol]
                            -bc_tg*bc_rhs*aii*u_old[indx_sol];   // Dirichlet tang  homogeneous only
          // non-diagonal part
          for(int   jvar=ivar+_dir+1; jvar< _nNSdim+ivar+_dir; jvar++)    {
            int      jvar1=jvar%DIMENSION;
            double aij=Ipenalty*normal[jvar1]*normal[ivar+_dir]*(bc_tg+bc_normal);
            double valj=_data_eq[2].ub[(_NS_idx+jvar1)*NDOF_FEM+sur_toply[lbnode]];
            FeM(indx_row) += -bc_rhs*bc_tg*aij*valj;
            FeM(indx_row) += aij*valj;

          }// jvar loop
          //                     bc_rhs*Ipenalty*u_tau[i]   + bc_yplus*Ipenalty*y1_plus[i];
        } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------
      }

    }  //  ivar |||||
  }  //  bc_dn |||||


  //  ============================================================================================
  //  =================== Dirichlet boundary condition (Dirichlet=0; Neumann=1) ====================
  //  ============================================================================================
  // enum bound_cond{simm=0,velocity_in0=1,velocity_tg0=2,wall=3,
  // simm1=4,velocity_in=5,velocity_tg=6,velocity=7,
  // outflow=10,pressure_outlet=12,interior=11,turb=11,outflow_p=14,pressure_inlet=16
  //  ============================================================================================

  for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {

    int bc_dn=(_bc_vol[sur_toply[lbnode]]/10)%2;
    if(bc_dn==0) {  // dirichlet -------------------------------------------------------------------------------------------------------------------

      int bc_s=_bc_vol[sur_toply[lbnode]]%10;
      int bc_tg=((bc_s&2)>>1);   // (0;?1?) tg
      int bc_normal=(bc_s%2);    // (0;??1) normal
      int bc_rhs=((bc_s&4)>>2);  // (0;1??) nonhomogeneous

      for(int ivar=0; ivar< _nvars[2]; ivar++)    {
        const int  indx_sol=sur_toply[lbnode]+(ivar+_dir)*el_ndof[2];//ivar=0;
        const int  indx_row=sur_toply[lbnode]+(ivar)*el_ndof[2];//ivar=0;
//             int bc_yplus=((bc_s&8)>>3);  // (0;1??) nonhomogeneous
        
        if(ivar+_dir == dir_maxnormal)  {  // normal ---------------------------------------------------------------------------------------------------

          _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]=bc_dn;
          //          Dirichlet (homogenoues and non  homogenoues) Assemblying
          double aii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir];
          KeM(indx_row,indx_row) += aii; //diagonal ->  u_n=n(u*nx+v*ny+w*nz)
          FeM(indx_row) += bc_normal*bc_rhs*aii*u_old[indx_sol]; // only normal non-homogeneous
          for(int   jvar=ivar+_dir+1; jvar< _nNSdim+_dir+ivar; jvar++)    {  // non-diagonal part  ->  u_n=n(u*nx+v*ny+w*nz)
            int      jvar1=jvar%DIMENSION;
            double aij= Ipenalty*normal[jvar1]*normal[ivar+_dir];
            double valj=_data_eq[2].ub[(_NS_idx+jvar1)*NDOF_FEM+sur_toply[lbnode]];
            FeM(indx_row) += bc_rhs*bc_normal*aij*valj;
            FeM(indx_row) += -aij*valj;

          }// jvar loop
          //    bc_rhs*Ipenalty*u_tau[i]   + bc_yplus*Ipenalty*y1_plus[i];
        } //  end  Dirichlet normal --------------------------------------------------------------------------------------------------------------------------

        else  {  // Non normal direction  ---------------------------------------------------------------------------------------------------------------

          _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= (_bc_vol[sur_toply[lbnode]]%10 == 0) ? 1:0 ; // storage bc for vol part
          //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
          // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
          //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
          // diagonal part
          KeM(indx_row,indx_row) +=  Ipenalty*(1-(1-bc_tg)*(1-bc_normal));
          double aii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*fabs(bc_tg-bc_normal);
          KeM(indx_row,sur_toply[lbnode]+ivar*el_ndof[2]) += -aii;
          FeM(indx_row) +=  bc_rhs*bc_tg*Ipenalty*u_old[indx_sol]
                            -bc_tg*bc_rhs*aii*u_old[indx_sol];   // Dirichlet tang  homogeneous only
          // non-diagonal part
          for(int   jvar=ivar+_dir+1; jvar< _nNSdim+ivar+_dir; jvar++)    {
            int      jvar1=jvar%DIMENSION;
            double aij=Ipenalty*normal[jvar1]*normal[ivar+_dir]*fabs(bc_tg-bc_normal);
            double valj=_data_eq[2].ub[(_NS_idx+jvar1)*NDOF_FEM+sur_toply[lbnode]];
            FeM(indx_row) += -bc_rhs*bc_tg*aij*valj;
            FeM(indx_row) += aij*valj;

          }// jvar loop
          //       bc_rhs*Ipenalty*u_tau[i]   + bc_yplus*Ipenalty*y1_plus[i];
        } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------

      }  // dc
    }// lnode

  } // l

  // ========================================= Dirichlet boundary conditions ==============================
}


//********************************************************
// Function for the calculation of u_tau
double MGSolNS::eval_var2(double var[]) { // var2=utau

  // var[0]=


  const double  u_b=var[0]; // boundary velocity
  double utau2 = 0.;  // linear-log intersection
  double utau3;       // buffer region

  // Calculation of linear region utau   u+=u/utau=y+ (linear) ; y+=y*utau/nu ---------------------
  double utau1 = sqrt(u_b *_IRe /(_y_bcout+1e-20));//  ---------------> sqrt(u_b * _IRe / (_y_bcout+1e-20)); // utau (linear) <- var[0],_y_bcout
  double yp = utau1*_y_bcout/_IRe;// --------------->utau1*_y_bcout/_IRe;               // y+ (linear)

  // Calculation of logarithmic region utau --------------------------------------------------------
  if(yp > 5.) { // non linear
    double diff = 100.; // diff=fabs(utau(k)-utau(k-1))
    double utau2_old =11.6/(u_b+1.e-20);
    double E=exp(5.2*0.41);
    while(diff > 1.e-6) { // y+=y*utau/nu
      utau2 = u_b*0.41/(log(E*_y_bcout*utau2_old/_IRe));
      diff = fabs(utau2_old - utau2);
      utau2_old = utau2;
    }
    yp = utau2*_y_bcout/_IRe;//utau2*_y_bcout/_IRe;
  } else {utau2 = 0.;}

  // Calculation of utau through the musker relation ------------------------------------------------
  if(yp > 5. && yp < 40.) { // musker region
    double var_in[1];
    double diff = 100.; // diff=fabs(utau(k)-utau(k-1))
    int cont = 0;
    double utau3_old = utau2;
    while(diff > 1.e-6) {
      var_in[0]= utau3_old;
      utau3 = u_b/eval_var3(var_in);
      diff = fabs(utau3_old - utau3);
      utau3_old = utau3;
      cont ++;   if(cont > 3000) {utau3 = 0.; break;}
    }
  }

  // final choice utau1 (linear), utau2 (log)
  double utau;
  if(yp > 5.)  { // the max between the log and musker region
    utau = ((utau2>utau3)? utau2:utau3);
    if(utau3>2.*utau2) { utau = utau2; } // if twice
  } else { utau = utau1; }
// // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// utau = utau1;
// // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   if(_y_bcout * utau/_IRe < 45.) { // from tau/a=-dP/dx (half channel width)
//     printf("ulog %g umusk %g utau %1.9f up %1.9f yp %f \n",(utau2-sqrt(_dirg[1]*9.81*0.03025))/sqrt(_dirg[1]*9.81*0.03025),
//            (umusk - sqrt(_dirg[1]*9.81*0.03025))/sqrt(_dirg[1]*9.81*0.03025), utau, sqrt(_dirg[1]*9.81*0.03025), _y_bcout * utau/_IRe);
//   }
  return utau;
}
//********************************************************


double MGSolNS::eval_var3(double var[]) { // var3=vel

  // var[0]= dist, var[1]=  utau, var[2]=  nu
  double yplus = _y_bcout*var[0]/_IRe;
  double vel = 5.424*atan((2.*yplus-8.15)/16.7)
               + 4.1693*log(yplus + 10.6)
               - 0.8686*log(yplus*yplus - 8.15*yplus+86)
               - 3.52;
  return vel;
}



// #endif



#endif
#endif