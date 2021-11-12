/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef VARIABLE_H 
#define VARIABLE_H 

#include "_class.h"
#include "_system.h"
#include "_ctrl.h"

extern int N;		//アウトプットナンバー(ファイル読み込み)
extern int RESTART;	//リスタートするかどうか(ファイル読み込み，0=off)// 0//0//1//0//1	//1=on
extern int RESTART_NUMBER;//リスタートステップ(ファイル読み込み)// 4000//14000
#define RESTART_N N//0000


extern double ene_pot;	//LJポテンシャル
extern double ene_kine;	//運動エネルギ？慣性項が合ったときの名残
extern double ene;		//エネルギ総和
extern _vec<double> p;	//運動量？慣性項が合ったときの名残

#ifdef SRD
extern _mpc* mpc;
extern _cell* cell;
extern _node* node;
#endif
//extern _ptcl* ptcl;
extern _PTCL ptcl;		//粒子情報
extern _list* list;
extern _list_sub list_sub;			//重合・脱重合のとき使うリスト
extern _list_sub_sev list_sub_sev;	//切断のとき使うリスト
// extern _numbering* numbering;

extern _Actin::_Plus_End* plus_end;		//フィラメントの端情報保持
extern _Actin::_Minus_End* minus_end;	


extern int* counter_num;

//extern int* bind_wall;


//extern int list_X;// = 0;
//extern int list_Y;// = 0;
//extern int list_Z;// = 0;
	
extern _calcArea calcArea;
extern int calc_length_X;
extern int calc_length_Y;
extern int calc_length_Z;

extern int list_size;// = 0;

#ifdef WALL_BOUNDARY
namespace WALL{
	extern _bind_wall* bind_wall;
	extern int bind_X;// = 0;
	extern int bind_Y;// = 0;
	extern int bind_size;// =0;
const int del_bind=4;
}
#endif


extern double t;
extern int step;
extern ofstream fout_err;

//extern ofstream fout_sev;
//extern stringstream ss_fout_sev;
//extern _check_sev* check_sev;
extern int his_sev_num[16];
extern int his_sev_the[180];

extern int no_fila;
extern double T;
extern int counter_pol;
extern int counter_depol;
extern int counter_severing;
extern int counter_severing_max;


#ifdef USE_MPI
#include<vector>
//#if NX>1
extern std::vector<_ptcl*> send_list_ptcl_x_pre;
extern std::vector<_ptcl*> send_list_ptcl_x_pre_inner;
extern std::vector<_ptcl*> send_list_ptcl_x_pre_inner_new;
extern std::vector<_ptcl*> send_list_ptcl_x_bck;
extern std::vector<_ptcl*> send_list_ptcl_x_bck_inner;
extern std::vector<_ptcl*> send_list_ptcl_x_bck_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_pre_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_x_bck_inner_new;
//#endif// ZX
//#if NY>1
extern std::vector<_ptcl*> send_list_ptcl_y_pre;
extern std::vector<_ptcl*> send_list_ptcl_y_pre_inner;
extern std::vector<_ptcl*> send_list_ptcl_y_pre_inner_new;
extern std::vector<_ptcl*> send_list_ptcl_y_bck;
extern std::vector<_ptcl*> send_list_ptcl_y_bck_inner;
extern std::vector<_ptcl*> send_list_ptcl_y_bck_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_pre_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_y_bck_inner_new;
//#endif//NY>1
//#if NZ>1
extern std::vector<_ptcl*> send_list_ptcl_z_pre;
extern std::vector<_ptcl*> send_list_ptcl_z_pre_inner;
extern std::vector<_ptcl*> send_list_ptcl_z_pre_inner_new;
extern std::vector<_ptcl*> send_list_ptcl_z_bck;
extern std::vector<_ptcl*> send_list_ptcl_z_bck_inner;
extern std::vector<_ptcl*> send_list_ptcl_z_bck_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_pre_inner_new;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck_inner;
extern std::vector<_ARP2_3::_arp2_3*> send_list_arp_z_bck_inner_new;
//#endif//NZ>1
#endif// USE_MPI

extern int Active_Ptcl;//=Number_Ptcl_Active+Number_Ptcl_Spare_Dens;

extern _vec<double> replenish_loc[Number_Ptcl_Spare_Dens];
#ifdef USE_MPI
extern int replenish_index[Number_Ptcl_Spare_Dens];
extern std::vector<int> replenish_changelist;
#endif //#ifdef USE_MPI



extern double* filament_force;
extern double* filament_force_sub;
extern double Wall_RandamF;


extern double force_wall_mono;
extern double force_wall_fila;


namespace WALL{
extern double jijuu;
}



#endif

