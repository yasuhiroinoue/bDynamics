/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef _CLASS_H
#define _CLASS_H

#include "vec.h"
#include "quaternion.h"
#include "_system.h"
#include "_ctrl.h"
#include "_output_class.h"

struct _list;
class _config;
namespace _ARP2_3{struct _arp2_3;}
// namespace _Actin{struct _Bond;}
#ifdef BROWNIAN
struct _ptcl{
	int index;
	_vec<double> loc;		//posision
	_vec<double> loc_pre;	//previus posision
	_vec<double> force;		//external force

	_quat<double> quat;		//角運動方程式用
	_quat<double> quat_pre;
	_vec<double> adhesion_loc;
	_vec<double> adhesion_loc_inite;
	_vec<double> adhesion_frc;


	_vec<double> udloc;//フィラメントのバネの方向(tmp)
	double norm;		//フィラメントのバネの長さ

	_config* cnf;//状態管理(どの粒子とつながっているか，など.)
	
	_ptcl();
	void add(_list*);
	_list* list_calc();
	bool Copy(const _ptcl& p){//_PTCLのCopyから以外の読み出し禁止
//indexとconfigは変えない
		loc.IN(p.loc);
		loc_pre.IN(p.loc_pre);
		force.IN(p.force);
		quat.IN(p.quat);
		quat_pre.IN(p.quat_pre);
		adhesion_loc.IN(p.adhesion_loc);
		adhesion_loc_inite.IN(p.adhesion_loc_inite);
		adhesion_frc.IN(p.adhesion_frc);
		return true;
	};
};

class _config{
public:
	_list* li;
	int list_num;
	void add(_list*);
	
	int feature;			//主に単量体かどうかの判断 モノマー=0
// 	int num;
	_ptcl* ptcl;//相互参照
	_ptcl* plus;
	_ptcl* minus;
	
	_ARP2_3::_arp2_3* edge_arp2_3;
	_ARP2_3::_arp2_3* end_arp2_3;

	int bind;				//boolian でいいんだけどね//_bind_wall

	_config():li(NULL),list_num(0),
			feature(0),ptcl(NULL),
			plus(NULL),minus(NULL),
			bind(-1){}

	bool Copy(const _config& p);//基本的に相互参照
};
	
class _PTCL{
public:
	_ptcl* pt;
	_config* con;
	int max_size;
	int used_size;
	
	_PTCL():pt(NULL),con(NULL),max_size(0),used_size(0){}
	_PTCL(const int& size){
		mynew(size);
	}
	~_PTCL(){
		delete[] pt;
		delete[] con;
	}
	void mynew(const int& size){
		used_size=0;
		max_size=size;
		{
			try{
				pt = new _ptcl[size];
			}catch(const std::bad_alloc&){std::cerr<<"allocerr"<<std::endl;exit(1);}
			try{
				con= new _config[size];
			}catch(const std::bad_alloc&){std::cerr<<"allocerr"<<std::endl;exit(1);	}
		}
		for(int i(0);i<size;i++){
			pt[i].index=i;
			pt[i].cnf=&con[i];
			con[i].ptcl=&pt[i];
		}
	};
	_ptcl& operator[] (const int& index){
		//if(index < used_size)
		return pt[index];
//		return NULL;
	}
	bool Copy(_ptcl& p1,_ptcl& p2){
		if(p1.index==p2.index)return false;
//		p2.cnf->li->rm(p2,p2.cnf->list_num);
		con[p1.index].Copy(con[p2.index]);
		p1.Copy(p2);
		return true;
	}
	bool ReSize(const int& size){
		//だるぃ
		return true;
	}
};














#else
struct _ptcl{
	_vec<double> loc;
	_vec<double> vel;
	_vec<double> frc_old;
	_vec<double> frc_new;
	double mass;
	int feature;//_numbering* number
	int num;
	double time;
	int p_level;
	int bind;
	_list* li;
};
#endif

#ifdef SRD
struct _mpc{
	_vec<double> loc;
	_vec<double> vel;
	double mass;
};

struct _cell{
	double mass;
	_vec<double> vel;
	_quat<double> rot;
};

struct _node{
//	double col;
	int count;
	double mass;
	_vec<double> mom;
	double s[6];
};
#endif

struct _total{
	public:
		double   mas;
		_vec<double> mom;
		double   ene;
};
/*struct _list{
	int num;
	_ptcl* index[500];
};*/

namespace _ARP2_3{
	struct _arp2_3{
		public:
		_vec<double> loc;
		_vec<double> loc_pre;
		_vec<double> force;
		int edge_ptcl_index;
		int end_ptcl_index;
		_ptcl* edge_ptcl();
		_ptcl* end_ptcl();
		void edge_ptcl_in(int);
		void end_ptcl_in(int);
		
		int list_index;
		_list* li();
		int list_num;
		
		int index;
		
		_arp2_3():loc(),loc_pre(),force(),edge_ptcl_index(-1),end_ptcl_index(-1),list_index(-1){}
		~_arp2_3(){
			loc=loc_pre=force=_vec<double>(0.0,0.0,0.0);
			edge_ptcl_index=end_ptcl_index=NULL;
		}
		_list* list_calc();
		void add(_list*);
	};
}




struct _list{//一セルないの情報
	int num;	//保持しているptclの総数
	int num_arp;	//保持しているarp23の総数
	int myrank;		//計算領域の確認
	_ptcl* ptcl_index[484];		//セル内粒子のindex，４８４は一セルないに入る最大粒子数’(適当50くらいでもいい？)ここを変えると下も同様に変える
	_ARP2_3::_arp2_3* arp_index[100];
	_list* list_index[27];		//近接のセルのポインター//２７は隣接のセル数，値は固定
//	_vec<int> index;
	_list(){//初期化(NULL埋め)
		num=0;
		num_arp=0;
		for(int i=0;i<484;++i)ptcl_index[i]=NULL;
		for(int i=0;i<100;++i)arp_index[i]=NULL;
		for(int i=0;i<27; ++i)list_index[i]=0;
	}
	void init(){
		num=0;
		for(int i=0;ptcl_index[i];i++)ptcl_index[i]=0;
	}
	void init_all(){
		num=0;
		num_arp=0;
		for(int i=0;i<484;++i)ptcl_index[i]=NULL;
		for(int i=0;i<100;++i)arp_index[i]=NULL;
	}
	void add(_ptcl* p){
		ptcl_index[num]=p;
		++num;
	}
	void add(_ARP2_3::_arp2_3* p){
		arp_index[num_arp]=p;
		++num_arp;
	}
	void rm(_ptcl* p){
		--num;
		int i=0;
		if(num){
			for(i=0;p!=ptcl_index[i];++i);
			if(i!=num){
				ptcl_index[i]=ptcl_index[num];
				ptcl_index[i]->cnf->list_num=i;
			}
		}
		ptcl_index[num]=NULL;
	}
	void rm(_ptcl*,int i){
		--num;
		if((num!=0) && (i!=num)){
			ptcl_index[i]=ptcl_index[num];
			ptcl_index[i]->cnf->list_num=i;
		}
		ptcl_index[num]=NULL;
		
	}
	void rm(_ARP2_3::_arp2_3*,int i){
		--num_arp;
		if((num_arp) && (i!=num_arp)){
			arp_index[i]=arp_index[num_arp];
			arp_index[i]->list_num=i;
		}
		arp_index[num_arp]=NULL;
	}
};

struct _list_sub{
	int num;
	_ptcl* index[500];
	
};
struct _list_sub_sev{
	int num;
	int index[500];
	
};
struct _bind_wall{
	_vec<double> loc;
	_ptcl* index;
};
struct _WALL_Z{
	_vec<double> loc;
	_vec<double> loc_pre;
	_vec<double> LJ;
	_vec<double> adhe;
	_vec<double> gaus;
#ifdef SRD
	_vec<double> gaus_srd;
#endif
	_WALL_Z(){
		loc=loc_pre=LJ=adhe=gaus=_vec<double>(0.0,0.0,0.0);
#ifdef SRD
		gaus_srd.IN(0.0,0.0,0.0);
#endif
	}
};

struct _total_sub{
	int counter;
	double   mas;
	_vec<double> mom;
	double   ene;
};

//struct _check_sev{
//	int n;
//	double t;
//	double the;
//};

// struct _numbering{
// 	int num;
// 	_ptcl* ptcl_n[number_max];
// 	_ARP2_3::_arp2_3* end_arp2_3;
// 	_numbering(){
// 		num=0;
// 		for(int i=0;i<number_max;i++){ptcl_n[i]=0;}
// 		return;
// 	}
// 	void init(){
// 		num=0;
// 		for(int i=0;i<number_max;i++){ptcl_n[i]=0;}
// 		return;
// 	}
// 	void sort_p(_ptcl* p){
// 		bool flag=false;
// 		for(int i=0 ; ptcl_n[i] ; i++){
// 			if(flag){
// 				ptcl_n[i-1]=ptcl_n[i];
// 				ptcl_n[i-1]->num--;
// 				ptcl_n[i]=0;
// 			}
// 			else{
// 				if(p==ptcl_n[i]){flag=true;	ptcl_n[i]=0;}//flag=(p==ptcl_n[i]);if(flag)ptcl_n[i]=0;
// 			}
// 		}
// 		//if(!flag){cout<<"error_numbering"<<endl;}
// 		return;
// 	}
// 	void sort_n(int n){
// 		ptcl_n[n]=0;
// 		for(int i=n;ptcl_n[i+1];i++){
// 			ptcl_n[i]=ptcl_n[i+1];
// 			ptcl_n[i]->num=i;
// 			ptcl_n[i+1]=0;
// 		}return;
// 	}
// };
// 

namespace _Actin{
	struct _Plus_End{
		int num;
		int end_index[feature_max];
		_ptcl* end(const int);
		_ptcl* end();
//		void end_in(const int);
		
		_Plus_End(){
			num=0;
			for(int i=0;i<feature_max;++i){end_index[i]=-1;}
		}
		~_Plus_End(){
			num=0;
			for(int i=0;i<feature_max;++i){end_index[i]=-1;}
		}
		void add(_ptcl* b);
		void rm(_ptcl* b);
		void polymerization_fila(int i,_ptcl* pa);
		void polymerization_mono(_ptcl* pa,_ptcl* pb);
		void depolymerization(_ptcl* b);
	};
	struct _Minus_End{
		int num;
		int end_index[feature_max];
		_ptcl* end(const int);
		_ptcl* end();
		
		_Minus_End(){
			num=0;
			for(int i=0;i<feature_max;++i){end_index[i]=-1;}
		}
		~_Minus_End(){
			num=0;
			for(int i=0;i<feature_max;++i){end_index[i]=-1;}
		}
		void add(_ptcl* b);
		void rm(_ptcl* b);
		void depolymerization(_ptcl* b,int i);
	};
}

#ifdef USE_MPI
enum conf_ch {no=0,pol_mono,pol_fila,depol,sev,wall_bind,wall_debind,arp_edge,arp_end,arp_debind};
namespace mympi{
struct _conformation_change{								//packing class	を作ったからあんまり意味がないなぁ...
	conf_ch cc;
	int index[2];
	_conformation_change(const conf_ch& a,const int& b,const int& c){
		cc=a;
		index[0]=b;
		index[1]=c;
	}
};
}
#endif//#ifdef USE_MPI

class _calcArea{
public:
	int nx,ny,nz;
	int x_pre,x_bck,y_pre,y_bck,z_pre,z_bck;
	int pair_pre_x,pair_bck_x,pair_pre_y,pair_bck_y,pair_pre_z,pair_bck_z;
	
	
//	_calcArea():x_pre(0),x_bck(0),y_pre(0),y_bck(0),z_pre,z_bck(0){}
	_calcArea():x_pre(0),x_bck(list_X),y_pre(0),y_bck(list_Y),z_pre(0),z_bck(list_Z),
				pair_pre_x(0),pair_bck_x(0),pair_pre_y(0),pair_bck_y(0),pair_pre_z(0),pair_bck_z(0){}
	void add(const int myrank,const int);
};




#endif
