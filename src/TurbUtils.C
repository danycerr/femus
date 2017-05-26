#include "TurbUtils.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

TurbUtils::TurbUtils(){
  FillModelMap();
  FillParameters(); // DEFAULT VALUES
}

TurbUtils::TurbUtils(double wall_distance, double TurbModel[], double nu, double alpha){
  FillModelMap();
  FillParameters(wall_distance, TurbModel, nu, alpha);
}

void TurbUtils::FillParameters(){
  _BoundWallDist = 1.e-3;
  _nu = 1.e-4;
  _alpha = 1.e-4;
  _klog  = 0.;
  _wlog  = 0.;
  _khlog = 0.;
  _whlog = 0.;
  return;
}
void TurbUtils::FillParameters(double wall_dist, std::string FirstDynEq, std::string SecDynEq, std::string FirstThermEq, std::string SecThermEq, double nu, double alpha)
{
  std::cout<<"\n ========================================================= \n"
  << "\033[38;5;118m \t \t SETTING THE TURBULENCE MODEL \033[0m \n";
  _BoundWallDist = wall_dist;
  _nu = nu;
  _alpha = alpha;
  
  int FDE = _DynamicModel[FirstDynEq];
  int SDE = _DynamicModel[SecDynEq];
  int FTE = _ThermalModel[FirstThermEq];
  int STE = _ThermalModel[SecThermEq];
  
  _Nagano = (FDE/10)%2;
  _Wilcox = ((FDE/10)&2)>>1;
  
  _klog  = (FDE%10)%2; 
  _wlog  = (SDE%10)%2;
  _emod  = ((SDE%10)&2)>>1;
  _numod = ((FDE%10)&2)>>1;

  _khlog = _ThermalModel[FirstThermEq]%10;
  _whlog = _ThermalModel[SecThermEq]%10;

  _nu    = nu;
  _alpha = alpha;  
  __IPr  = _alpha/_nu; 
  _IsFilled = true;
  
  return;
}

void TurbUtils::FillParameters(double wall_dist, 
			       std::vector< std::string > TurbModel, 
			       std::vector< std::string > SolveEqs,
			       std::vector<std::string> Constrain, 
			       double UnderRelaxation[2],
			       double nu, 
			       double alpha)
{
  std::cout<<"\n ========================================================= \n"
  << "\033[38;5;118m \t \t SETTING THE TURBULENCE MODEL \033[0m \n";
  _BoundWallDist = wall_dist;
  _nu = nu;
  _alpha = alpha;
  __IPr  = _alpha/_nu; 
  
  _DynUnderRel = (UnderRelaxation[0]==NULL)? 0.:UnderRelaxation[0];
  _ThermUnderRel = (UnderRelaxation[1]==NULL)? 0.:UnderRelaxation[1];
  
  int FDE = _DynamicModel[TurbModel[0]];
  int SDE = _DynamicModel[TurbModel[1]];
  int FTE = _ThermalModel[TurbModel[2]];
  int STE = _ThermalModel[TurbModel[3]];
  
  _Nagano = (FDE/10)%2;
  _Wilcox = ((FDE/10)&2)>>1;
  
  _klog  = (FDE%10)%2; 
  _wlog  = (SDE%10)%2;
  _emod  = ((SDE%10)&2)>>1;
  _numod = ((FDE%10)&2)>>1;

  _khlog = FTE%10;
  _whlog = STE%10;
  
  std::map<std::string, bool> YesNo;
  YesNo["yes"] = true;
  YesNo["no"]  = false;
  
  _SolveNS   = YesNo[SolveEqs[0]];
  _SolveT    = YesNo[SolveEqs[1]];
  _SolveTBK  = YesNo[SolveEqs[2]];
  _SolveTTBK = YesNo[SolveEqs[3]];  
  
  _Durbin    = YesNo[Constrain[0]]; 
  _YapCorr   = YesNo[Constrain[1]]; 
  
  _IsFilled = true;
  PrintStatus(TurbModel);
  return;
}


void TurbUtils::FillModelMap(){

  _DynamicModel["nagano_k"]      = nagano_k    ;
  _DynamicModel["nagano_logk"]   = nagano_logk ;
  _DynamicModel["nagano_w"]      = nagano_w    ;
  _DynamicModel["nagano_logw"]   = nagano_logw ;
  _DynamicModel["nagano_e"]      = nagano_e    ;
  _DynamicModel["wilcox_k"]      = wilcox_k    ;
  _DynamicModel["wilcox_logk"]   = wilcox_logk ;
  _DynamicModel["wilcox_nut"]    = wilcox_nut  ;
  _DynamicModel["wilcox_w"]      = wilcox_w    ;
  _DynamicModel["wilcox_logw"]   = wilcox_logw ;
  
  
  _ThermalModel["natural_kh"] = natural_kh;
  _ThermalModel["logarithmic_kh"] = logarithmic_kh;
  _ThermalModel["natural_omegah"] = natural_omegah;
  _ThermalModel["logarithmic_omegah"] = logarithmic_omegah; 
  return;
}

void TurbUtils::FillParameters(double wall_distance, double TurbModel[], double nu, double alpha){
  _BoundWallDist = wall_distance;
  _nu = nu;
  _alpha = alpha;
  _klog  = TurbModel[0];
  _wlog  = TurbModel[1];
  _khlog = TurbModel[2];
  _whlog = TurbModel[3];
  __IPr  = _alpha/_nu; 
  return;
}

double TurbUtils::CalcUtau(double vel_bound, double dist) {
    double umusk, ulog, ulin, utau, ulold, diff;
    double beta, vk;

    double yp;

    beta = 5.2;
    vk = 0.41;

    ulog = 0.;
    ulold = 11.6/vel_bound;

    // Calculation of utau through the linear relation
    ulin = sqrt(vel_bound * _nu / dist);
    yp = ulin*dist/_nu;

    // Calculation of utau through the logarithmic relation
    diff = 100.;
    if(yp > 5.) {
        while (diff > 1.e-6) {
            ulog = vel_bound*vk/(log(exp(beta*vk)*dist*ulold/_nu));
            diff = fabs(ulold - ulog);
            ulold = ulog;
        }
        yp = ulog*dist/_nu;
    }
    else {
        ulog = 0.;
    }

    // Calculation of utau through the musker relation

    if(yp > 5. && yp < 40.) {
        diff = 100.;
        int cont = 0;
        double umuskold = ulog;

        while (diff > 1.e-6) {
            umusk = vel_bound/Musker(dist, umuskold);
            diff = fabs(umuskold - umusk);
            umuskold = umusk;
            cont ++;
            if(cont > 3000) {
                umusk = 0.;
                break;
            }
        }
    }

    if(yp > 5.)  {
        utau = max (ulog,umusk);
        if(umusk>2.*ulog) utau = ulog;
    }
    else         utau = ulin;

    return utau;
}

double TurbUtils::Musker (double dist, double utau) {
    double yplus = dist*utau/_nu;
    double vel = 5.424*atan((2.*yplus-8.15)/16.7)
                 + 4.1693*log(yplus + 10.6)
                 - 0.8686*log(yplus*yplus - 8.15*yplus+86)
                 - 3.52;
    return vel;
}


double TurbUtils::CalcMuTurb(double KappaAndOmega[], double dist, double vel_sp){

  if(_numod) _MuTurb = KappaAndOmega[0]/_nu;
  else{// WILCOX AND NAGANO TURBULENCE MODELS - NOT FOR NUT EQUATION    
    double kappa, omega, epsilon;    
    {// REAL K, OMEGA AND EPSILON CALCULATION
      kappa = (1.-_klog)* KappaAndOmega[0] + _klog*exp(_klog*KappaAndOmega[0]);    
      if(_emod){ 
        epsilon = KappaAndOmega[1];
        omega   = epsilon/(kappa*__CMU);
      }
      else omega = (1.-_wlog)* KappaAndOmega[1] + _wlog*exp(_wlog*KappaAndOmega[1]);
      kappa = (kappa>1.e-20)? kappa:1.e-20;
      omega = (omega>1.e-20)? omega:1.e-20;
    }
    
    {// EDDY VISCOSITY CALCULATION
      __Ret   = kappa/(_nu*omega);               //  viscosity ratio
      __Rt    = __Ret/__CMU;                     // turbulent Reynolds number
      __Rd    = dist*sqrt(kappa/sqrt(__Rt))/_nu;          
      _MuTurb = __Ret;                           // turbulent viscosity ratio
      if(_Nagano){
        __fmu   = (1.-exp(-1.*__Rd/14.))*(1.-exp(-1.*__Rd/14.));
        __fcorr = 1. + 5./pow(__Rt,0.75)*exp(-1.*__Rt*__Rt/40000.);
        _MuTurb *= __fcorr*__fmu;
      }
    }
    if(_Durbin){
     if (_MuTurb > kappa/(sqrt(vel_sp + 1.e-10)*_nu)) {
        _MuTurb = kappa/(sqrt(vel_sp + 1.e-10)*_nu); 
     }
//     if (muturb > sqrt(2./3.)*kappa/(sqrt(0.5*vel_sp)*_nu)) muturb = sqrt(2./3.)*kappa/(sqrt(0.5*vel_sp)*_nu);  // 3D
    }    
  }
  return _MuTurb;
}

void TurbUtils::CalcDynTurSourceAndDiss(double KappaAndOmega[], double dist, double vel_sp ,double &muturb, double source[2], double diss[2]){
  
  // REAL KAPPA, OMEGA AND EPSILON
  double kappa,omega,epsilon;
  kappa =  omega = epsilon =1.;
  kappa = KappaAndOmega[0];
  omega = KappaAndOmega[1];
  epsilon = KappaAndOmega[1];
  
  if(!_numod) kappa = (1.-_klog)* KappaAndOmega[0] + _klog*exp(_klog*KappaAndOmega[0]);
  if(!_emod)  omega = (1.-_wlog)* KappaAndOmega[1] + _wlog*exp(_wlog*KappaAndOmega[1]);  

  kappa = (kappa>1.e-20)? kappa:1.e-20;
  omega = (omega>1.e-20)? omega:1.e-20;
  if(_emod) epsilon = (epsilon>1.e-20)? epsilon:1.e-20;

  const double kCorr = (_klog/kappa + (1.-_klog));
  const double wCorr = (_wlog/omega + (1.-_wlog));  
  
  // Values of __Rt and __Rd are assigned within TurbUtils::CalcMuTurb
  muturb =  CalcMuTurb(KappaAndOmega,dist);
  double prod_k = 0.5 * _nu * muturb * vel_sp;
  prod_k = min(prod_k, sqrt(8./3.)*kappa*sqrt(vel_sp));  // k production limitation -> Park
  
  if(_Nagano){// NAGANO TURBULENCE MODEL-----------------------------------------------------------------------
    const double f_exp  = (1.-exp(-1.*__Rd/(3.1)))*(1.-exp(-1.*__Rd/(3.1)))*(1.- 0.3*exp(-1.*__Rt*__Rt/42.25));
    if(_emod){// KAPPA - EPSILON TURBULENCE MODEL
      diss[0]   = __CMU*epsilon/kappa;      // correggere
      diss[1]   = epsilon*__C20*f_exp/kappa;
      source[0] = prod_k*kCorr;
      source[1] = __C10*prod_k*epsilon/kappa;         
    }
    else{// KAPPA - OMEGA AND LOG FORMULATION TURBULENCE MODEL
      
      diss[0]   = __CMU*omega;
      diss[1]   = __CMU*omega*(__C20*f_exp-1);
      source[0] = prod_k*kCorr;
      source[1] = (__C10-1.)*prod_k*omega/kappa*wCorr;     
      if(_YapCorr){
        double yap_term = 0.83*kappa*kappa*(sqrt(kappa)*wCorr/(__CMU*omega*2.44*dist)-1.)/__CMU;
        source[1] += max(yap_term,0.);
      }
    }
  }//----------------------------------------------------------------------------------------------------------
  if(_Wilcox){// WILCOX TURBULENCE MODEL-----------------------------------------------------------------------
    // FIRST EQUATION SOURCE AND DISSIPATION    
    if(_numod){// NUT EQUATION
      diss[0]   = __BETAN*omega;
      source[0] = 0.5*vel_sp*__AN/omega;      
    }
    else{// KAPPA AND LOG(KAPPA) EQUATION  
      diss[0]   = __BETAS*omega;
      source[0] = prod_k*kCorr;      
    }
    // SECOND EQUATION SOURCE AND DISSIPATION
    diss[1]   = __BETAW*omega;
    source[1] = __AW * 0.5 * vel_sp*wCorr;
  }//----------------------------------------------------------------------------------------------------------  
  return;
}

double TurbUtils::CalcAlphaTurb(double KappaAndOmega[], double TKappaAndOmega[], double dist)
{
  // REAL OMEGA AND OMEGAT
  double omega  = (1.-_wlog)* KappaAndOmega[1] + _wlog*exp(_wlog*KappaAndOmega[1]);
  double omegaT = (1.-_whlog)* TKappaAndOmega[1] + _whlog*exp(_whlog*TKappaAndOmega[1]);
  
  // MU_TURB calculation
  // __fmu __fcorr __Ret __Rt __Rd
  double nut = CalcMuTurb(KappaAndOmega, dist);
  nut /= (__fmu*__fcorr); // -> nut/nu
  __rT       = omega/omegaT;
  __F1t     = (1.-exp(-__Rd/(sqrt(__IPr)*14.)))*(1.-exp(-__Rd/14.));
  __F2at    = exp(-4.e-6*__Rt*__Rt);
  __F2bt    = exp(-2.e-5*__Rt*__Rt);  
    
  const double  IPrdlT = (0.1/__CMU)*__F1t*(
                           __Prdl_inf                                            /* Asymptotic contribution */
                           + 2.*__rT/(__rT+0.3)*__F2at                           /* Contribution far from wall */
                           + 1.3*__IPr*sqrt(2.*__rT)/(pow(__Rt,0.75))*__F2bt     /* Near Wall contribution */
                         );  
   const double alphaT = nut * _nu * IPrdlT;  
   return alphaT; 
}

void TurbUtils::CalcThermTurSourceAndDiss(double KappaAndOmega[],   // Dynamic Turbulence
					  double TKappaAndOmega[],  // Thermal Turbulence
					  double dist,              // Wall Distance
					  double sp,                // Squared value of velocity derivatives
					  double st,                // Modulus of temperature gradient
					  double &alphaturb,        // Eddy thermal diffusivity
					  double source[],          // Source terms
 					  double diss[],            // Dissipation terms
					  double meccterm[])        // Mechanical contribution
{
  // REAL K AND OMEGA, MECHANICAL AND THERMAL
  double kappa  = (1.-_klog)* KappaAndOmega[0] + _klog*exp(_klog*KappaAndOmega[0]);
  double omega  = (1.-_wlog)* KappaAndOmega[1] + _wlog*exp(_wlog*KappaAndOmega[1]);
  double kappaT = (1.-_khlog)* TKappaAndOmega[0] + _khlog*exp(_khlog*TKappaAndOmega[0]);
  double omegaT = (1.-_whlog)* TKappaAndOmega[1] + _whlog*exp(_whlog*TKappaAndOmega[1]);
  
  // ALPHA_TURB CALCULATION
  // __fmu __fcorr __Ret __Rt __Rd __rT __F1t __F2at __F2bt
  alphaturb = CalcAlphaTurb(KappaAndOmega, TKappaAndOmega, dist);

  const double TkCorr = ((1.-_khlog) + _khlog/kappaT);
  const double TwCorr = ((1.-_whlog) + _whlog/omegaT);
  const double f_exp  = (__CD2*(1.-0.3*exp(-__Rt*__Rt/42.25))-1.)*(1. - exp(-__Rd/5.7))*(1.-exp(-__Rd/5.7));

  double muturb = __fmu*__fcorr*kappa/omega;
  
    if(_Durbin){
      if (muturb > kappa/(sqrt(sp)*_nu))  muturb = kappa/(sqrt(sp)*_nu); 
    }
  
  const double prod_k = 0.5*sp*muturb;
  const double prod_kt = alphaturb*st;
   
  // THERMAL SOURCE
  source[0] = prod_kt*TkCorr;
  source[1] = (__CP1-1.)*prod_kt*TwCorr*omegaT/kappaT;
   
  // THERMAL DISSIPATION
  diss[0]   = __CMU*omegaT;
  diss[1]   = (__CD1-1.)*__CMU*omegaT;
  
  // MECHANICAL DISSIPATION AND SOURCE
  meccterm[0]  = f_exp*__CMU*omega;                      // DISSIPATION
  meccterm[1]  = __CP2*prod_k*omegaT*TwCorr/kappa;       // SOURCE

  if(_YapCorr){
    const double yap_term = 0.83*kappa*kappa*(sqrt(kappa)/(__CMU*omega*2.44*dist)-1.)/(__CMU * (_wlog + (1.-_wlog)*omega));
    meccterm[1]  += max(yap_term,0.);
  }
  return;
}


void TurbUtils::DynTurInitValues(double & kappa, double & omega, double WallDist, double AvVel, double Diameter, bool FlatProfile){
  
  if(FlatProfile){    
    const double Re_h     = AvVel*Diameter/_nu;
    const double In       = 0.16*pow(Re_h,-0.125);
    const double len      = 0.07*Diameter; 
    const double k_in     = 1.5*(In*AvVel)*(In*AvVel);
    const double e_in     = __CMU*k_in*sqrt(k_in)/len;
    const double w_in     = sqrt(k_in)/len;
    const double n_in    = k_in/w_in;
    
    if(_numod) kappa = n_in;
    else kappa = k_in;
    if(_emod) omega = e_in;
    else omega = w_in;    
  }
  else{
    // FIRST MESH POINT ASSUMED TO BE AT y+ = 10.
    const double utau = 10.*_nu/_BoundWallDist;
    
    double yp     = WallDist*utau/_nu;    
    double nut    = _nu/(1./(0.001093*yp*yp*yp) + 1./(0.41*yp));
    double w_lin  = 2.*_nu/(0.09*WallDist*WallDist);    
    double w_log  = utau/(sqrt(__CMU)*0.41*WallDist);    
    double w      = sqrt(w_lin*w_log);
    double k      = w * nut;        
    if(_numod) kappa = nut;
    else kappa = k;
    if(_emod) omega = utau*utau*utau/(0.41*WallDist);
    else omega = w;    
  }
  
  return;
}

void TurbUtils::TherTurInitValues(double& kappaT, double& omegaT, double WallDist, double NormFlux, double AvVel, double Diameter, bool FlatProfile)
{
  
  double kappa, omega, yp;
  DynTurInitValues(kappa, omega, WallDist, AvVel, Diameter, FlatProfile);
  
  if(FlatProfile){
    kappaT = kappa*(_alpha/_nu);
    omegaT = omega*(_alpha/_nu);
    kappaT = kappa;
    omegaT = omega;
    
  }
  else{    
    // FIRST MESH POINT ASSUMED TO BE AT y+ = 10.
    const double utau  = 10.*_nu/_BoundWallDist;
    const double T_tau = NormFlux/utau;            // NormFlux needs to be _qs/_rhof*_cp0    
    yp     = utau*WallDist/_nu;

    double wh1, wh2;
    wh2    = omega*sqrt(_alpha/_nu);
    wh1    = 2.*_alpha/(__CMU*WallDist*WallDist);    
    omegaT     = sqrt(wh1*wh2);
    omegaT     = omega*4;                                   // Assumed value of R = omega/omegaT = 0.25
    kappaT     = kappa*T_tau*T_tau/(utau*utau*_alpha/_nu);    
    const double yPr = yp*_nu/_alpha;
    const double khlin = 0.0001*pow(yPr,1.5);
    const double khlog = 1./(yPr*yPr*yPr);
    kappaT     = khlin*khlog/(khlin+khlog);
  }  
  
  return;
}




void TurbUtils::PrintStatus(std::vector<std::string>TurbModels)
{

  ofstream output;
  output.open("SimulationSpecs.txt");
  output << " MESSAGE PRINTED BY TURBUTILS CLASS "<<endl;
  output << endl << " TurbUtils filled: "<<_IsFilled<<endl;
  output << "---------------------------------------" <<endl;
  output << " Solved Equations: "<<endl;
  output << "\t Navier Stokes:        "<<_SolveNS<<endl;
  output << "\t Temperature:          "<<_SolveT<<endl;
  output << "\t Dynamical Turbulence: "<<_SolveTBK<<endl;
  output << "\t Thermal Turbulence:   "<<_SolveTTBK<<endl;
  output << "---------------------------------------" <<endl;
  output << " Turbulence Models: "<<endl;
  output << "\t Dynamic, first equation:  "<<TurbModels[0]<<endl;
  output << "\t Dynamic, second equation: "<<TurbModels[1]<<endl;
  output << "\t Thermal, first equation:  "<<TurbModels[2]<<endl;
  output << "\t Thermal, second equation: "<<TurbModels[3]<<endl; 
  output << "---------------------------------------" <<endl;
  output << " Realizability and Constrains:"<<endl;
  output << "\t Durbin Realizability Constrain: "<<_Durbin<<endl;
  output << "\t Yap Correction:                 "<<_YapCorr<<endl;
  output << "---------------------------------------" <<endl;
  output << " Geometry Parameters: "<<endl;
  output << "\t Fixed Wall Distance: "<<_BoundWallDist <<endl;
  output << "---------------------------------------" <<endl;
  output << " Physical Properties: "<<endl;
  output << "\t Kinematic Viscosity: "<<_nu<<endl;
  output << "\t Thermal diffusivity: "<<_alpha<<endl;
  
  return;
}
