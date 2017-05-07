#ifndef _config_cfg
#define _config_cfg

#include "Domain_conf.h"


#define VOF_SUBSTEP   1
// #define RESTART 10
// #define SPLIT
// cc level
#define CCLEV 3


// ***************************************************
// 2D 	( DIM2= defined = two-dimensional simulations )
// ***************************************************
#if DIMENSION==2 // 2D dim=2 --------------------------------
// #define DIMENSION    2
#define NX (40)        // number of nx element
#define NY (40)        // number of ny element
#define NZ (0)         // number of nz element
#define AX (0.)        // x-interval  (AX,BX)
#define AY (0.)        // y-interval  (AY,BY)
#define AZ (0.)        // z-interval  (AZ,BZ)
#define BX (1.)        // x-interval  (AX,BX)
#define BY (1.)        // y-interval  (AY,BY)
#define BZ (0.)        // z-interval  (AZ,BZ)

#endif


// *******************************************************
// 3D 	( DIM2= undefined = three-dimensional simulations )
// *******************************************************
#if DIMENSION==3  // 3D dim=3 ------------------------------
// #define DIMENSION    3
#define NX (40)        // number of nx element
#define NY (40)        // number of ny element
#define NZ (40)         // number of nz element
#define AX (0.)        // x-interval  (AX,BX)
#define AY (0.)        // y-interval  (AY,BY)
#define AZ (0.)        // z-interval  (AZ,BZ)
#define BX (1.)        // x-interval  (AX,BX)
#define BY (1.)        // y-interval  (AY,BY)
#define BZ (1.)        // z-interval  (AZ,BZ)



#endif


#endif
