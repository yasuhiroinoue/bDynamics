/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include"_ctrl.h"
#ifdef USE_MPI

#ifndef ___MYMPI_H
#define ___MYMPI_H

#include"_class2.h"

namespace mympi{
extern packer bit_sender;
extern int pair_pre;
extern int pair_bck;

void mympi_init(int&,char**&,int&,int&);
void mympi_finalize();
void mympi_ptcl_init(const int myrank,const int numprocs);

void mympi_send_pre(const int myrank,const int numprocs);
void mympi_send_bck(const int myrank,const int numprocs);
void mympi_send_bck_second(const int myrank,const int numprocs);

void mympi_conformation_change_nucleation(const int myrank,const int numprocs);
void mympi_conformation_change(const int myrank,const int numprocs);
void mympi_location_gather(const int myrank,const int numprocs);
void mympi_replenish_ptcl(const int myrank,const int numprocs);
void mympi_flow_wallinf_bcast(const int myrank,const int numprocs);

}//namespace mympi

#endif //___MYMPI_H

#endif //USE_MPI

