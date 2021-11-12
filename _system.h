/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef SYSTEM_H 
#define SYSTEM_H 

//#define N 2


#include <cmath>
#include "mymath.h"
#include "_ctrl.h"

#define INTVL_COL 100
//MPC_dt/PTCL_dt
#define INTVL_step 1
//(100/INTVL_COL)




#define SYS_kT 1.0//0.5//1.0//2.0//0.5	//0.001//0.50//ガウス分布定数 温度 0.10 - 1.0


#define SYS_X 32			//48//60//32//64//32
#define SYS_Y 32			//48//60//32//64//32
#define SYS_Z 64//48			//48//60//32//64//32
const double WALL_Z_INIT=33.0;
//const double WALL_Z_INIT=50-8;//SYS_Z-8;
const int list_X(12);
const int list_Y(12);
const int list_Z(24);


#ifdef WALL_BOUNDARY
namespace WALL{
//const double jijuu=140.0;//1.0;
extern double Wall_Friction;//=100.0;
}

#define		SVX		0.0000		//Z平面 速度x
#define		SVY		0.00		//Z平面 速度y
#endif //#ifdef WALL_BOUNDARY


const int Number_Ptcl_Active=4096;
//const int Number_Ptcl_Active=4596;
const int Number_Ptcl_Spare=4096;
const int Number_Ptcl_Spare_Dens=512;
//const int Number_Ptcl_Active=(64*16);//1024;
//const int Number_Ptcl_Spare=0;
//const int Number_Ptcl_Spare_Dens=0;
const int Number_Ptcl = Number_Ptcl_Active+Number_Ptcl_Spare;//4096;//2048//2572//3072//2048//10000//2048//1//2048
//14//2048//0//4096//4092//4096//4000//0//4000//2000//6000//2000//16//4096//11264//4096//2048//16//2000//2000	//2000	//4000	//2000	//0//4096//11264//1331 //ptcl粒子数

const int offset_impound=3;

#define STEP_EVRY 100
extern int STEP_START;// 0
extern int STEP_END;// 104001//44001//34001//100001//34001//1000//34001
#define AVR_TIME 2000
#define MPC_OUTPUT_TIME 2000
#define restart_out 10000

#define P_monitor 10
#define STEP_MONITOR  100
#define PVWin_OUT 10
#define PTCL_OUT 10
#define WALL_FRC_OUT 1//100

#define fila_count_out 1000
#define fila_arg 1000
#define fila_arg1 200
#define fila_out 10
const int replenish_slice=1;

#define _GETA 1E-5 //ゲタ

//ptcl
//#define feature_max 6000
//#define number_max 10000
const int feature_max = Number_Ptcl/2+1;//2000;//Number_Ptcl/2+1;
const int number_max = Number_Ptcl+1;
const int segment_max = 1000;

//#define PTCL_MASS 5.0			//20//4.0//20//4.4	//10.0

const double PTCL_dt=0.000005;//0.000007											//0.000002			//0.00001//
//const double Ptcl_Friction=1.0;													//0.02//0.1//1000.0
extern double Ptcl_Friction;//=1.0;													//0.02//0.1//1000.0
const double sigma_LJ=1.0;													//1.2		//1.0	//ptcl 大きさ
//const double Ptcl_Rot_Friction=(Ptcl_Friction*4.0/3.0*sigma_LJ*sigma_LJ*0.25);		//0.01//100.0
extern double Ptcl_Rot_Friction;//=(Ptcl_Friction*4.0/3.0*sigma_LJ*sigma_LJ*0.25);		//0.01//100.0
const double epsilon_LJ=0.187804878;											//0.000187805	//0.25
const double r_cut = sigma_LJ*pow(2.0,1.0/6.0);//2.50 * sigma_LJ;					//cutoff
const double r_cut2 = r_cut*r_cut;					//cutoff

const double cutoffuper_sub=0.8;
const double cutoffuper=cutoffuper_sub*sigma_LJ;
const double cutoffuper2=cutoffuper*cutoffuper;


const double k_bane=2000.0;								//2100.00	//100.0		//バネ定数
const double rb=(sigma_LJ*0.5);								//0.0			//バネ釣り合い位置
const double k_angl=500.0;//1000.0;//10.0//1.0//15.0		//13.00									//500.0	//角バネ定数
const double eq_angl = PI;

//const double lence_eq = pow(4.0*epsilon_LJ*12.0*pow(sigma_LJ,12.0)/k_bane,1.0/14.0);

const double r_por_f_g = 1.02*sigma_LJ;			//重合範囲・cut_off以下・listを使って計算しているため．
const double r_por_g_g = 1.2*sigma_LJ;			//重合範囲・cut_off以下・listを使って計算しているため．
// const double r1 = 7.0;							//重合範囲.エネルギ
const double por_angl = PI*60.0/180.0;			//重合範囲

const double probability_deporimerize  =0.00007;			//0.00015//0.00045
//const double probability_porimerize_f_p=1e-3;//1e-4;
extern double probability_porimerize_f_p;//=1e-3;//1e-4;
//5e-5;//1e-4;
//5e-5;
//5.27E-08;//5.0E-09;
//0.1;//1.0;				//0.005//1.0//0.01//0.8
const double probability_porimerize_p_p=0.0;//0.00000005//0.0000001//1.0				//0.0000001//0.00001
//#define probability_severing 0.0000002
const double pro_sev=0.00000025;					//0.0000001;//0.000005;//0.0000001;//0.001;
const double l_eq=1.0440;							//1.0;	//1.0440;//流速０の時の<l>
const double e_c=0.05;								//10000;//0.01;//0.05;



#define p_time 150.0




#define K_Torsion k_angl		//50.0//1.0//20		//100.00
#define hanshu 13.5
#define k_Bend_angl 0.1*k_angl//0.0//k_angl	//13.00
//k_angl

//const double lattice_cut=_ARP2_3::SIGMA_LJ_ARP_ARP;//=3.20 * sigma_LJ;//r_cut;




//#ifdef SRD
//
//
//const double MPCMASS = 1.0;
//const double MPC_dt=(INTVL_COL*PTCL_dt);
////1.0
//
//const int cell_size = SYS_X*SYS_Y*SYS_Z;
//
//const double VX=0.0;		//ガウス分布定数 平均速度x
//const double VY=0.0;		//ガウス分布定数 平均速度y
//const double VZ=0.0;		//ガウス分布定数 平均速度z
//
//const double number_dens=5;
//		//mpc粒子数 密度5-10になるように.rho-5
//const int Number_mpc=(number_dens*SYS_X*SYS_Y*SYS_Z);
//	//327680		//163840		//327680		//262144		//163840		//1310720		//163840 //mpc粒子数 密度5-10になるように.rho-5
//
//const double rho_mpc = (double)number_dens*MPCMASS;
//
//const double G=0.00050;
////PI*0.50;//RAND*PI;//衝突角度を変えるときは、ここを変える。(例) tc = PI*0.50;// 90度
//const double col_angl = PI*130.0/180.0;//PI*45.0/180.0;//PI*130.0/180.0; //2.2689280275926284500007979990352;
//const double col_cell = 1.0;
//
//const int SYS_X_col=SYS_X/col_cell;
//const int SYS_Y_col=SYS_Y/col_cell;
//const int SYS_Z_col=SYS_Z/col_cell;
//const int cell_size_col = SYS_X_col*SYS_Y_col*SYS_Z_col;
//
//const double number_dens_cell = (double)number_dens * col_cell * col_cell * col_cell;
//
//
//
//#endif

#endif
