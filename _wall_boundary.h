/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include"_ctrl.h"
#ifdef WALL_BOUNDARY
#ifndef __WALL_BOUNDARY_H
#define __WALL_BOUNDARY_H

#include "_class.h"

namespace WALL{

extern double ene_pot;
// extern double WALL_Z;
extern struct _WALL_Z WALL_Z[2];
void wall_bind_init();
void Wall_ptcl_frc(int);
void bind_cal(int myrank,int numprocs);
void bind_frc(int,int);
void Wall_Location_FirstStep(int myrank,int numprocs);
void Wall_Location_SecondStep(int myrank,int numprocs);

}//namespace WALL
#endif //__WALL_BOUNDARY_H
#endif //WALL_BOUNDARY
