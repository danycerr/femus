#ifndef _domainconf
#define _domainconf

//-------------------------------------------------------------
//Dimension of the problem
// DIMENSION = 2 or 3 (D)
#define DIMENSION  (2)
//  #define AXISYM
#define NUM_MESH (1)
//-----------------------------------------------------------

#define  BDRY_TOLL  1.e-12 //tolerance for setting the BCs

// libemesh no this 
 #define MATBC_INTERFACE       // boundary conditions and material

// #define HAVE_GROUP (1) // see in the gencase file only for med

//Physical dimensions of the mesh:

const int AX_DIR = 1;

#if DIMENSION==2  // =============================

  #define LXB (0.)
  #define LXE (.3)
  #define LYB (0.)
  #define LYE (1.)
  #define LZB (0.)
  #define LZE (0.)


//   #define LXB (0.)
//   #define LXE (0.03025)
//   #define LYB (0.)
//   #define LYE (0.1)
//   #define LZB (0.)
//   #define LZE (0.)

#endif

#if DIMENSION==3 // =========================

  #define LXB (0.)
  #define LXE (0.03025)
  #define LYB (0.)
  #define LYE (0.01)
  #define LZB (0.)
  #define LZE (0.01)

#endif
#endif
