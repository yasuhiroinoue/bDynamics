/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef FUNCTION_H 
#define FUNCTION_H 

#include "_system.h"
#include "_ctrl.h"

#ifdef SRD
#include "_srd.h"
#endif

//void Periodic(_vec<double>& a);

void ptcl_initial_condition();
void PTCL_output(double t);
void ene_format();//初期化
void Frc();//力
void Frc(const int,const int);//力
void PTCL_energy();//エネルギ
void PTCL_Location();//位置
void PTCL_Velocity();//速さ

namespace _Actin{

#ifdef POLYMERIZATION
void polymerization(const int myrank,const int numprocs);
void nucleation(const int myrank,const int numprocs);
#endif //POLYMERIZATION
#ifdef DEPOLYMERIZATION
void depolymerization(const int myrank,const int numprocs);
#endif //DEPOLYMERIZATION
#ifdef SEVERING
void severing(const int myrank,const int numprocs);
#endif //SEVERING

}//namespace _Actin{


namespace PV_out{
	void PVWin_init_bin(const int,const int);//PVWin用初期化
	void PVWin_output_bin(const int,const int);//PVWin用アウトプット
	//void PVWin_init_restart_bin();
}
void file_init();
#ifdef WALL_BOUNDARY
void frc_wall_output(const int, const int);
void OutputPlusEnd(const int myrank,const int);
#endif
//void loc_format();


// void p_level_cal();
// void list_cal();

void systemmonitor_ptcl(int step);
void len_out();
void arg_out();

void struct_init(const int myrank,const int numprocs);
void struct_finalize();

// void bind_frc();
// void bind_cal();
//void ptcl_debug(double t);

void restart_binary_output(int restart_step);
void restart_binary_input(int restart_step);
void restart_binary_input(const int ,const int ,int);
void ptcl_fila_count();
void history(string c);
void dir_init();


void debug();
//void numbering_debug();



void Ptcl_Location_FirstStep(const int myrank,const int numprocs);
void Ptcl_Location_SecondStep(const int myrank,const int numprocs);


//------------------------------------------------------------------------------
void pdf_strain(int start_step , int end_step=STEP_END ,int del_n=50);

_list* list_calc(const _vec<double>&);

void replenish_ptcl(const int myrank,const int numprocs);


void output_cc(const int myrank, const int numprocs);
void output_polymerization_info(const int myrank, const int numprocs);
void output_frc_wall(const int myrank, const int numprocs);

void ptcl_wall_interaction(const int myrank, const int numprocs);
void ptcl_number_output(const int myrank, const int numprocs);

void OutputClassUpdate(const int myrank, const int numprocs);
void OutputForceCount(const int myrank,const int);
void Flow_Output(const int myrank,const int);

void InputInitefile(const int& argc,char**& argv,const int myrank,const int numprocs);

#endif

