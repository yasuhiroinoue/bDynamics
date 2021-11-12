/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef _CTRL_H
#define _CTRL_H

#define BROWNIAN

//#define USE_MPI

//////////#ifdef USE_MPI
//////////const int NX=4;
//////////const int NY=1;
//////////const int NZ=1;
////////#define NX 4
////////#define NY 1
////////#define NZ 1
//////////ÇÌÇØÇ™ÇÌÇ©ÇËÇ‹ÇπÇÒÇ≈Ç∑ÇÀÅD
////////#else
//////////ïœÇ¶ÇƒÇÕë ñ⁄
////////#define NX 1
////////#define NY 1
////////#define NZ 1
extern int NX;// 4
extern int NY;// 1
extern int NZ;// 1
//////////#endif//USE MPI

#define POLYMERIZATION
#define DEPOLYMERIZATION
//#define SEVERING

#define WALL_BOUNDARY

#ifdef WALL_BOUNDARY
// #define SHAKE_BOTH_SIDE
// #define WALL_CYCLIC

#define WALL_BIND
#endif

#define DEBUG_PERIODIC


//#define SRD

#ifdef SRD
//#define COOLDOWN
#endif


#ifdef WALL_BOUNDARY
#define COUETTE_BOTH_SIDE
#endif

#define ARP2_3

//#define INITIAL_CONDITION_DEF

#endif
