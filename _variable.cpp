/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include <sstream>
#include <fstream>
#include "_variable.h"
#include "_system.h"

int N(0);
int RESTART(0);// 0//0//1//0//1	//1=on
int RESTART_NUMBER(4000);// 4000//14000
int NX(0);
int NY(0);
int NZ(0);

double probability_porimerize_f_p(1.0);
int STEP_START(0);// 0
int STEP_END(104001);// 104001//44001//34001//100001//34001//1000//34001
double Ptcl_Friction(1.0);													//0.02//0.1//1000.0
double Ptcl_Rot_Friction(Ptcl_Friction*4.0/3.0*sigma_LJ*sigma_LJ*0.25);		//0.01//100.0
namespace WALL{
	double jijuu(100);
	double Wall_Friction(100);
}

double ene_pot(0.0);
double ene_kine(0.0);
double ene(0.0);
_vec<double> p(0.0,0.0,0.0);

//_ptcl* ptcl(NULL);
_PTCL ptcl;//(Number_Ptcl);
#ifdef SRD
_mpc* mpc;
_cell* cell;
_node* node;
#endif //SRD
_list* list(NULL);

_list_sub list_sub;
_list_sub_sev list_sub_sev;
// _numbering* numbering;

//namespace _Actin{
// 	_Actin::_Bond* bond=NULL;
// 	_Actin::_Bend* bend=NULL;
//	int counter_bond=0;
// 	int counter_bend=0;
	_Actin::_Plus_End* plus_end;
	_Actin::_Minus_End* minus_end;	

//}


int* counter_num(0);

//int* bind_wall;


//int list_X(0);
//int list_Y(0);
//int list_Z(0);

_calcArea calcArea;
	
int calc_length_X(0);
int calc_length_Y(0);
int calc_length_Z(0);

int list_size(0);

#ifdef WALL_BOUNDARY
namespace WALL{
_bind_wall* bind_wall;
int bind_X = 0;
int bind_Y = 0;
int bind_size =0;
}
#endif


double t = 0.0;
int step(STEP_START);

std::ofstream fout_err;

//ofstream fout_sev;
//stringstream ss_fout_sev;
//_check_sev* check_sev;
int his_sev_num[16];
int his_sev_the[180];

int no_fila = 0;
double T;
int counter_pol;
int counter_depol;
int counter_severing=0;
int counter_severing_max=0;


// _vec<double> average_randam_force(0.0,0.0,0.0);
// _vec<double> veriance_randam_force(0.0,0.0,0.0);
// const _vec<double> Inite_X=_vec<double>((double)SYS_X*0.5,(double)SYS_X*0.5,(double)SYS_X*0.5);
// _vec<double> veriance_X(0.0,0.0,0.0);

#ifdef USE_MPI
#include<vector>
//#if NX>1
std::vector<_ptcl*> send_list_ptcl_x_pre;
std::vector<_ptcl*> send_list_ptcl_x_pre_inner;
std::vector<_ptcl*> send_list_ptcl_x_pre_inner_new;
std::vector<_ptcl*> send_list_ptcl_x_bck;
std::vector<_ptcl*> send_list_ptcl_x_bck_inner;
std::vector<_ptcl*> send_list_ptcl_x_bck_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck_inner_new;
//#endif// ZX
//#if NY>1
std::vector<_ptcl*> send_list_ptcl_y_pre;
std::vector<_ptcl*> send_list_ptcl_y_pre_inner;
std::vector<_ptcl*> send_list_ptcl_y_pre_inner_new;
std::vector<_ptcl*> send_list_ptcl_y_bck;
std::vector<_ptcl*> send_list_ptcl_y_bck_inner;
std::vector<_ptcl*> send_list_ptcl_y_bck_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck_inner_new;
//#endif//NY>1
//#if NZ>1
std::vector<_ptcl*> send_list_ptcl_z_pre;
std::vector<_ptcl*> send_list_ptcl_z_pre_inner;
std::vector<_ptcl*> send_list_ptcl_z_pre_inner_new;
std::vector<_ptcl*> send_list_ptcl_z_bck;
std::vector<_ptcl*> send_list_ptcl_z_bck_inner;
std::vector<_ptcl*> send_list_ptcl_z_bck_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre_inner_new;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck_inner;
std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck_inner_new;
//#endif//NZ>1
#endif// USE_MPI

_vec<double> replenish_loc[Number_Ptcl_Spare_Dens];
#ifdef USE_MPI
int replenish_index[Number_Ptcl_Spare_Dens];
std::vector<int> replenish_changelist;
#endif //#ifdef USE_MPI


int Active_Ptcl(Number_Ptcl_Active+Number_Ptcl_Spare_Dens);


double* filament_force;
double* filament_force_sub;
double Wall_RandamF;


double force_wall_mono(0);
double force_wall_fila(0);

