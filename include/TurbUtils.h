#include <sstream>
#include <vector>
#include <map>
#ifndef __TurbUtils__
#define __TurbUtils__

class TurbUtils{
private:  
  // NAGANO-DERIVED TURBULENCE MODEL COEFFICIENTS
  const double __CMU = 0.09;
  const double __CP1 = 1.025;
  const double __CP2 = 0.9;
  const double __CD1 = 1.1;
  const double __CD2 = 1.9;
  const double __C20 = 1.9;
  const double __C10 = 1.5;
  const double __Prdl_inf = 0.75;
  
  // WILCOX-DERIVED TURBULENCE MODEL COEFFICIENTS
  const double __AW = 13./25.;
  const double __AN = 1. - __AW;
  const double __BETAS = 8./100.;
  const double __BETAW = 3./40.;
  const double __BETAN = __BETAS - __BETAW;
  
  // DYNAMICAL TURBULENCE FUNCTIONS
  double __Rt    ;
  double __Rd    ;
  double __fmu   ;
  double __fcorr ;
  double __Ret   ;
  
  // THERMAL TURBULENCE FUNCTIONS
  double __IPr      ;  
  double __rT       ;
  double __F1t      ;
  double __F2at     ;
  double __F2bt     ;    
   
public: 
  double _utau;
  double _kwall;
  double _wwall;
  double _BoundWallDist;
  double _MuTurb;
  double _AlphaTurb;
  double _nu;
  double _alpha;
  int _wlog, _klog, _khlog, _whlog; 
  int _Nagano, _Wilcox;
  int _numod, _emod;
  double _vmid, _diameter;
  bool   _IsFilled = false;
  bool   _SolveNS, _SolveT, _SolveTBK, _SolveTTBK, _YapCorr, _Durbin;
  double _DynUnderRel, _ThermUnderRel;
  enum DynTurbModel{
//     natural_k = 10,
//     logarithmic_k = 11,
//     natural_omega = 20,
//     logarithmic_omega = 21
    // NAGANO BASED MODELS
    nagano_k    =10,
    nagano_logk =11,
    nagano_w    =10,
    nagano_logw =11,
    nagano_e    =12,
    // WILCOX BASED MODELS
    wilcox_k    =20,
    wilcox_logk =21,
    wilcox_nut  =22,
    wilcox_w    =20,
    wilcox_logw =21
  };
  enum ThermTurbModel{
    natural_kh = 10,
    logarithmic_kh = 11,
    natural_omegah = 20,
    logarithmic_omegah = 21
  };
  
  std::map<std::string, DynTurbModel> _DynamicModel;
  
  std::map<std::string, ThermTurbModel> _ThermalModel; 
  
  
public:
  TurbUtils();  
  TurbUtils( double wall_dist, double TurbModel[], double nu, double alpha = 1.e-4);
  ~TurbUtils(){}
  
  double CalcUtau(double vel_bound, 
		  double dist);
  double CalcMuTurb(double KappaAndOmega[], 
		    double dist,
		   double vel_sp = 1.);
  void CalcDynTurSourceAndDiss(double KappaAndOmega[],
			       double dist, 
			       double vel_sp ,
			       double &muturb, 
			       double source[2], 
			       double diss[2]);
  double CalcAlphaTurb(double KappaAndOmega[],
		       double TKappaAndOmega[],
		       double dist);
  void CalcThermTurSourceAndDiss(double KappaAndOmega[], 
				 double TKappaAndOmega[], 
				 double dist, 
				 double sp, 
				 double st, 
				 double &alphaturb, 
				 double source[2], 
				 double diss[2], 
				 double meccterm[2]);
  
  void FillParameters(double wall_dist, 
		      double TurbModel[],
		      double nu, 
		      double alpha=1.e-4);
  void FillParameters(double wall_dist, 
		      std::string val1, 
		      std::string val2, 
		      std::string val3, 
		      std::string val4, 
		      double nu=1.e-4, 
		      double alpha=1.e-4);
  void FillParameters(double wall_dist, 
		      std::vector<std::string> TurbModel,
		      std::vector<std::string> SolveEqs,
		      std::vector<std::string> Constrain,
		      double UnderRelaxation[],
		      double nu=1.e-4, 
		      double alpha=1.e-4
		      );
  void FillParameters();
  void FillModelMap();
  
  void PrintStatus(std::vector<std::string> TurbModel);
  
  void DynTurInitValues(double & kappa, double & omega, double WallDist, double AvVel, double Diameter, bool FlatProfile);
  void TherTurInitValues(double & kappaT, double & omegaT, double WallDist, double NormFlux, double AvVel, double Diameter, bool FlatProfile);
  
  inline bool CheckIfFilled(){return _IsFilled;}
  inline void SetAvVel(double vel){ _vmid = vel;}
  inline void SetDiameter(double diameter){ _diameter = diameter;}
  double Musker (double dist, double utau) ;  
};


#endif