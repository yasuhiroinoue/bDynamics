/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<cstring>
#include "_system.h"
#include "_class.h"
#include "_arp2_3.h"
#include "_variable.h"
#include "_wall_boundary.h"

using namespace std;

namespace PV_out{
string pv_filename;//outputFileName
const int MOLTYPE=0;		//1
//#define NMOL (Number_Ptcl*2+Number_Arp2_3*2)
#ifdef WALL_BOUNDARY
const int NMOL=(Number_Ptcl*2+Number_Arp2_3+18+9);//18+9は壁の線を書く用の粒子
#else
const int NMOL=(Number_Ptcl*2+Number_Arp2_3);
#endif
const int NBOND=Number_Ptcl+36;
//#if RESTART==1
//const int NDATE=((STEP_END-RESTART_NUMBER)/PVWin_OUT);
//#else
//const int NDATE=STEP_END/PVWin_OUT;
//#endif

/*------------------PVWin-------------------*/

void PVWin_init_bin(const int myrank,const int numprocs){
int NDATE;
if(RESTART==1)NDATE=((STEP_END-RESTART_NUMBER)/PVWin_OUT);
else NDATE=STEP_END/PVWin_OUT;

class header{ 
public:
	char ch[32];
	unsigned short moltype,nmol,nbond,ndata;
	double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
	header(int NDATE):moltype(MOLTYPE),nmol(NMOL),nbond(NBOND),
		vlx1(0.0),vlx2((double)SYS_X),vly1(0.0),vly2((double)SYS_Y),vlz1(0.0),vlz2((double)SYS_Z),
		time0(0.0),dt((double)PVWin_OUT){strncpy(ch,"PV-01 /Shoji-Maruyama Laboratory",32);ndata=(NDATE);}
	~header(){}
}h(NDATE);
	{//output名の代入
		stringstream ss;
		ss.str("");
//		ss<<"actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ss<<"actinpv/_actin"<<setw(4)<<setfill('0')<<N<<'_'<<myrank<<".pv";
		pv_filename=ss.str();
	}
	if(NMOL){//粒子直径andヘッダ書き出し
		double *s;
		s=new double[NMOL];
		for(int i=0;i<NMOL;++i){
			if(i<Number_Ptcl)s[i]=sigma_LJ;
			else if(i<Number_Ptcl*2) s[i]=0.25;
			else if(i<Number_Ptcl*2+Number_Arp2_3) s[i]=_ARP2_3::SIGMA_LJ_ARP_ARP;
//			else s[i]=0.25;
#ifdef WALL_BOUNDARY
			else s[i]=0.0;
#endif
			
		}
		ofstream fout(pv_filename.c_str(),ios::out | ios::binary);
		fout.write((char*)&h,sizeof(header));
		fout.write((char*)s,sizeof(double)*NMOL);
		delete[] s;
		fout.close();
	}return;
}

#pragma pack(1)
struct atom{
	unsigned char attr;
	short x,y,z;
	atom():attr(0),x(0),y(0),z(0){}
	~atom(){attr=0;x=y=z=0;}
	template <class TEMPLATE> void IN(const TEMPLATE& loc_x,const TEMPLATE& loc_y,const TEMPLATE& loc_z){x=(short)(loc_x*100);y=(short)(loc_y*100);z=(short)(loc_z*100);}
	template <class TEMPLATE> void IN(const _vec<TEMPLATE>& loc){IN(loc.x ,loc.y ,loc.z);}
	template <class TEMPLATE> void IN(const TEMPLATE (&loc)[3]){IN(loc[0] ,loc[1] ,loc[2]);}
	void IN(const _ptcl& p){IN(p.loc);}//基本的に定義は少なく..inline展開されるはずだから
	void IN(const unsigned char& atr ,const _ptcl& p){attr=atr;IN(p.loc);}
	template <class TEMPLATE> void IN(const unsigned char& atr ,const _vec<TEMPLATE>& loc){attr=atr;IN(loc);}
	template <class TEMPLATE> void IN(const unsigned char& atr ,const TEMPLATE (&loc)[3]){attr=atr;IN(loc);}
};
#pragma pack()
#pragma pack(1)
struct Bond{
	unsigned char Attr;
	unsigned short p1,p2;
	Bond():Attr(0xff),p1(0),p2(0){}//非表示で確保
	~Bond(){Attr=0;p1=p2=0;}
};
#pragma pack()

void PVWin_output_bin(const int myrank,const int numprocs){
	if(NMOL){
		atom* at;
		Bond* bond;
		at = new atom[NMOL];
		bond = new Bond[NBOND];
//		for(int i=0 ; i<Number_Ptcl ; i++){
		for(int i=0 ; i<Active_Ptcl ; i++){
			_ptcl& p = ptcl[i];
			if(p.cnf->feature)	at[i].attr=(unsigned char)((p.cnf->feature-1)%5)+1;
			else at[i].attr=(unsigned char)0;
// 			at[i].attr=(unsigned char)p.li->myrank;
			at[i].IN(p);
		}
//		for(int i=Number_Ptcl ; i<Number_Ptcl*2 ; i++){
		for(int i=Number_Ptcl ; i<Number_Ptcl+Active_Ptcl ; i++){
			_ptcl& p = ptcl[i-Number_Ptcl];
			if(p.cnf->feature)	at[i].attr=(unsigned char)13;
			else at[i].attr=(unsigned char)14;
			at[i].IN(p.loc+p.adhesion_loc);
		}
		int in(0);
//		for(_ptcl** end=&plus_end->end[0];*end;++end){
		for(int i(0);i<plus_end->num;++i){
			in++;
			_ptcl* end=plus_end->end(i);
			at[end->index].attr=in;
			if(end->index>=64*32)at[end->index].attr=8;
			for(_ptcl* p=end;p=p->cnf->minus;){
				at[p->index].attr=in;
				if(p->index>=64*32)at[p->index].attr=8;
			}
			if(in>=5)in=0;
		}
//		for(atom* p=at+(Number_Ptcl_Active+Number_Ptcl_Spare_Dens);p!=at+Active_Ptcl;p++)p->attr=6;
		for(atom* p=at+(Active_Ptcl);p!=at+Number_Ptcl;p++)p->attr=7;
		
		
#ifdef ARP2_3
		for(int i=Number_Ptcl*2 ; i<Number_Ptcl*2+Number_Arp2_3 ; i++){
			_ARP2_3::_arp2_3* p = &_ARP2_3::arp2_3[i-2*Number_Ptcl];
			if(p->edge_ptcl())at[i].attr=(unsigned char)11;
			else at[i].attr=(unsigned char)12;
// 			at[i].attr=(unsigned char)p->li->myrank;
			at[i].IN(p->loc);
		}
//		for(int i=Number_Ptcl*2+Number_Arp2_3 ; i<Number_Ptcl*2+Number_Arp2_3*2 ; i++){
//			_ARP2_3::_arp2_3* p = &_ARP2_3::arp2_3[i-(2*Number_Ptcl+Number_Arp2_3)];
//			at[i].attr=(unsigned char)12;
//			at[i].x=(short)((p->loc.x+p->adhesion_loc.x)*100);
//			at[i].y=(short)((p->loc.y+p->adhesion_loc.y)*100);
//			at[i].z=(short)((p->loc.z+p->adhesion_loc.z)*100);
//		}
#endif
#ifdef WALL_BOUNDARY
		{
			int i=Number_Ptcl*2+Number_Arp2_3;
			for(int z=0;z<=1;z++){
			for(int y=0;y<=2;y++){
			for(int x=0;x<=2;x++){
	 			at[i].attr=(unsigned char)15;
				at[i].x=(short)(x*SYS_X*50-x);
				at[i].y=(short)(y*SYS_Y*50-y);
// #ifdef SHAKE_BOTH_SIDE
// 				if(!z)at[i].z=(short)(WALL::WALL_Z*100);
// #else
// 				if(!z)at[i].z=(short)0;
// #endif
// 				else at[i].z=(short)((SYS_Z-WALL::WALL_Z)*100);
				at[i].z=(short)(WALL::WALL_Z[z].loc.z*100);
				++i;
			}}}
			for(int y=0;y<=2;y++){
			for(int x=0;x<=2;x++){
	 			at[i].attr=(unsigned char)15;
				at[i].x=(short)(x*SYS_X*50-x);
				at[i].y=(short)(y*SYS_Y*50-y);
				at[i].z=(short)(offset_impound*SYS_Z/(double)list_Z*100);
				++i;
			}}
		}
#endif
		
		int counter_bond=0;
// 		for(int i=1 ; numbering[i].num ; i++){
// 			_ptcl* p1=numbering[i].ptcl_n[0];
// 			_ptcl* p2;
// 			for(int j=0 ; p2=numbering[i].ptcl_n[j+1] ; j++){
// 				bond[counter_bond].Attr=0x07;
// 				bond[counter_bond].p1=p1->index + 1;
// 				bond[counter_bond].p2=p2->index + 1;
// 				counter_bond++;
// 				if(counter_bond==NBOND)break;
// 				p1=p2;
// 			}if(counter_bond==NBOND)break;
// #ifdef ARP2_3
// 			_ARP2_3::_arp2_3 *arp=NULL;
// 			for(int j=0 ; p1=numbering[i].ptcl_n[j] ; j++){
// 				if(arp=p1->edge_arp2_3){
// 					bond[counter_bond].Attr=0x37;
// 					bond[counter_bond].p1=p1->index + 1;
// 					bond[counter_bond].p2=arp->index +Number_Ptcl*2 + 1;
// 					counter_bond++;
// 					if(counter_bond==NBOND)break;
// 					if(p2=arp->end_ptcl){
// 						bond[counter_bond].Attr=0x57;
// 						bond[counter_bond].p1=p2->index + 1;
// 						bond[counter_bond].p2=arp->index +Number_Ptcl*2 + 1;
// 						counter_bond++;
// 						if(counter_bond==NBOND)break;
// 					}
// 				}
// 			}if(counter_bond==NBOND)break;
// #endif
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl];++p){
			if(p->cnf->plus){
				bond[counter_bond].Attr=0x07;
				bond[counter_bond].p1=(p - &ptcl[0]) + 1;
				bond[counter_bond].p2=(p->cnf->plus - &ptcl[0]) + 1;
				counter_bond++;
				if(counter_bond==NBOND)break;
			}
			if(p->cnf->edge_arp2_3){
				bond[counter_bond].Attr=0x37;
				bond[counter_bond].p1=(p - &ptcl[0]) + 1;
				bond[counter_bond].p2=(p->cnf->edge_arp2_3 - _ARP2_3::arp2_3) + 1+Number_Ptcl*2;
				counter_bond++;
				if(counter_bond==NBOND)break;
				
			}
			if(p->cnf->end_arp2_3){
				bond[counter_bond].Attr=0x57;
				bond[counter_bond].p1=(p - &ptcl[0]) + 1;
				bond[counter_bond].p2=(p->cnf->end_arp2_3 - _ARP2_3::arp2_3) + 1+Number_Ptcl*2;
				counter_bond++;
				if(counter_bond==NBOND)break;
				
			}
		}
#ifdef WALL_BOUNDARY
		{
			int i=Number_Ptcl*2+Number_Arp2_3;
			for(int z=0;z<=2;z++){
			for(int y=0;y<=2;y++){
			for(int x=0;x<=2;x++){
				if(x<=1){
					bond[counter_bond].Attr=0xf7;
					bond[counter_bond].p1=i + 1;
					bond[counter_bond].p2=i + 2;
					counter_bond++;
					if(counter_bond==NBOND)break;
				}
				if(y<=1){
					bond[counter_bond].Attr=0xf7;
					bond[counter_bond].p1=i + 1;
					bond[counter_bond].p2=i + 4;
					counter_bond++;
					if(counter_bond==NBOND)break;
				}
				++i;
			}if(counter_bond==NBOND)break;
			}if(counter_bond==NBOND)break;
			}
		}
#endif
		
		
		ofstream fout(pv_filename.c_str(),ios::out | ios::app | ios::binary);
		fout.write((char*)at,NMOL*sizeof(atom));
		fout.write((char*)bond,NBOND*sizeof(Bond));
		fout.close();
		delete[] at;
		delete[] bond;
	}return;
}
}	//namespace PV_out
















/*

void PVWin_init_bin(){
	if(Number_Ptcl){
		stringstream ss;
		ss.str("");
		ss<<"./actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ofstream fout(ss.str().c_str(),ios::out | ios::binary);
		struct header{ 
			char ch[32];
			int moltype,nmol,nbond,ndata;
			double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
		}h;
		memset(&h,0,sizeof(header));
		strncpy(h.ch,"PV-32 /Shoji-Maruyama Laboratory",32);
		h.moltype = (int)0;
		h.nmol = (int)NMOL;
		h.nbond = (int)NBOND;
		h.ndata = (int)NDATE;
		h.vlx1=0.0;	h.vlx2=(double)SYS_X;
		h.vly1=0.0;	h.vly2=(double)SYS_Y;
		h.vlz1=0.0;	h.vlz2=(double)SYS_Z;
		h.time0=0.0;	h.dt=(double)PVWin_OUT;
		fout.write((char*)&h,sizeof(header));
		
		double *s;
		s=new double[NMOL];
		for(int i=0;i<Number_Ptcl;i++){
			s[2*i]=1.0;	s[2*i+1]=0.25;
		}
		fout.write((char*)s,sizeof(double)*NMOL);
		
		fout.close();
	}
}
#pragma pack(1)
struct atom{
	unsigned char attr;
	int x,y,z;
//コンストラクタ
	atom(){
		attr=0;
		x=y=z=0;
	}
//デストラクタ
	~atom(){
		attr=0;
		x=y=z=0;
	}
};
#pragma pack()
#pragma pack(1)
struct Bond{
	unsigned char Attr;
	int p1,p2;
//コンストラクタ
	Bond(){	Attr=0xff;	p1=p2=0;}
	~Bond(){Attr=0;	p1=p2=0;	}
//デストラクタ
};
#pragma pack()

void PVWin_output_bin(){
	if(Number_Ptcl){
		atom* at;
		Bond* bond;
		at = new atom[NMOL];
		bond = new Bond[NBOND];
		for(int i=0 ; i<Number_Ptcl ; i++){
			_ptcl* p = &ptcl[i];
			at[2*i].attr=(unsigned char)1;
			at[2*i].x=(int)(p->loc.x*100);
			at[2*i].y=(int)(p->loc.y*100);
			at[2*i].z=(int)(p->loc.z*100);
			at[2*i+1].attr=(unsigned char)14;
			at[2*i+1].x=(int)((p->loc.x+p->adhesion_loc.x)*100);
			at[2*i+1].y=(int)((p->loc.y+p->adhesion_loc.y)*100);
			at[2*i+1].z=(int)((p->loc.z+p->adhesion_loc.z)*100);
		}
		for(int i=0; i<Number_Ptcl-1 ; i++){
			bond[2*i].Attr=0x05;
			bond[2*i].p1=2*i+1;
			bond[2*i].p2=2*i+3;
			bond[2*i+1].Attr=0x53;
			bond[2*i+1].p1=2*i+1;
			bond[2*i+1].p2=2*i+2;
		}
			bond[2*(Number_Ptcl-1)].Attr=0x53;
			bond[2*(Number_Ptcl-1)].p1=2*(Number_Ptcl-1)+1;
			bond[2*(Number_Ptcl-1)].p2=2*Number_Ptcl;
		
		
		stringstream ss;
		ss.str("");
		ss<<"./actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ofstream fout(ss.str().c_str(),ios::out | ios::app | ios::binary);
		fout.write((char*)at,NMOL*sizeof(atom));
		fout.write((char*)bond,NBOND*sizeof(Bond));
		fout.close();
		delete[] at;
	}
}
*/
/*void PVWin_init_bin(){
	if(Number_Ptcl){
		stringstream ss;
		ss.str("");
		ss<<"./actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ofstream fout(ss.str().c_str(),ios::out | ios::binary);
		struct header{ 
			char ch[32];
			unsigned short moltype,nmol,nbond,ndata;
			double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
		};
		if(Number_Ptcl>=65535 || STEP_END/PVWin_OUT>=65535){	fout_err<<"pvwin_err\n";	}
			// >=ffffのチェック
		header h;
		memset(&h,0,sizeof(header));
		strncpy(h.ch,"PV-01 /Shoji-Maruyama Laboratory",32);
			//difference between "pv-01" & "pv-32" is difference of headertype...?
		h.moltype = (unsigned short)1;
		h.nmol = (unsigned short)Number_Ptcl;
		h.nbond = (unsigned short)0;
		h.ndata = (unsigned short)(STEP_END/PVWin_OUT);
		h.vlx1=0.0;	h.vlx2=(double)SYS_X;
		h.vly1=0.0;	h.vly2=(double)SYS_Y;
		h.vlz1=0.0;	h.vlz2=(double)SYS_Z;
		h.time0=0.0;	h.dt=(double)PVWin_OUT;
		fout.write((char*)&h,sizeof(header));
		
		double s[3];
		memset(s,0,3*sizeof(double));
		switch(h.moltype){
			case 1:
				s[0]=1.0;
				break;
			case 2:
				s[0]=1.2;
				s[1]=0.6;
				break;
			case 3:
				s[0]=1.2;
				s[1]=0.6;
				s[2]=0.6;
			break;
		}
		fout.write((char*)s,sizeof(double)*h.moltype);
		fout.close();
	}
}
void PVWin_init_restart_bin(){
	if(Number_Ptcl!=0){
		stringstream ss;
		ss.str("");
		ss<<"./actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ofstream fout(ss.str().c_str(),ios::out | ios::binary);
		struct header{ 
			char ch[32];
			unsigned short moltype,nmol,nbond,ndata;
			double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
		};
		if(Number_Ptcl>=65535 || STEP_END/PVWin_OUT>=65535){	fout_err<<"pvwin_err\n";	}
			// >=ffffのチェック
		header h;
		memset(&h,0,sizeof(header));
		strncpy(h.ch,"PV-01 /Shoji-Maruyama Laboratory",32);
		h.moltype = (unsigned short)1;
		h.nmol = (unsigned short)Number_Ptcl;
		h.nbond = (unsigned short)0;
		h.ndata = (unsigned short)((STEP_END-RESTART_NUMBER)/PVWin_OUT);
		h.vlx1=0.0;	h.vlx2=(double)SYS_X;
		h.vly1=0.0;	h.vly2=(double)SYS_Y;
		h.vlz1=0.0;	h.vlz2=(double)SYS_Z;
		h.time0=(double)RESTART_NUMBER;	h.dt=(double)PVWin_OUT;
		fout.write((char*)&h,sizeof(header));
		
		double s[3];
		memset(s,0,3*sizeof(double));
		switch(h.moltype){
			case 1:
				s[0]=1.0;
				break;
			case 2:
				s[0]=1.2;
				s[1]=0.6;
				break;
			case 3:
				s[0]=1.2;
				s[1]=0.6;
				s[2]=0.6;
			break;
		}
		fout.write((char*)s,sizeof(double)*h.moltype);
		fout.close();
	}
}
#pragma pack(1)
struct atom{
  unsigned char attr;
  short x,y,z;
};
#pragma pack()
void PVWin_output_bin(){
	if(Number_Ptcl!=0){
		atom* at;
		at = new atom[Number_Ptcl];
		memset(at,0,Number_Ptcl*sizeof(atom));
		for(int i=0 ; i<Number_Ptcl ; i++){
#if 0
			_ptcl* p = &ptcl[i];
			int col = p->feature;
			int p_l = p->p_level;
			int bind =p->bind;
			int pvwin = 14;//15 - col;
			
			if(col != 0 && p_l ==1) pvwin = 14;
			if(col != 0 && p_l ==0) pvwin = 0;
			if(bind ==1)pvwin = 15;
			if( col == 0 ) pvwin = 1;
			
			at[i].attr=(unsigned char)pvwin;
			at[i].x=(short)(ptcl[i].loc.x*100);
			at[i].y=(short)(ptcl[i].loc.y*100);
			at[i].z=(short)(ptcl[i].loc.z*100);
#endif
#if 1
			_ptcl* p = &ptcl[i];
			
			int pvwin = p->feature % 10 + 1;
			if(p->bind == 1 )pvwin = 15;
			if(p->feature == 0 ) pvwin = 14;
			
			at[i].attr=(unsigned char)pvwin;
			at[i].x=(short)(p->loc.x*100);
			at[i].y=(short)(p->loc.y*100);
			at[i].z=(short)(p->loc.z*100);
#endif
#if 0
			_ptcl* p = &ptcl[i];
			
			int pvwin = p->feature % 10 + 1;
			if(p->bind == 1 )pvwin = 15;
			if(p->feature == 0 ) pvwin = 14;
			if(i==3588)pvwin = 11;
			if(i==0)pvwin = 0;
			
			at[i].attr=(unsigned char)pvwin;
			at[i].x=(short)(p->loc.x*100);
			at[i].y=(short)(p->loc.y*100);
			at[i].z=(short)(p->loc.z*100);
#endif
		}
		stringstream ss;
		ss.str("");
		ss<<"./actinpv/_actin"<<setw(4)<<setfill('0')<<N<<".pv";
		ofstream fout(ss.str().c_str(),ios::out | ios::app | ios::binary);
		fout.write((char*)at,Number_Ptcl*sizeof(atom));
		fout.close();
		memset(at,0,Number_Ptcl*sizeof(atom));
		delete[] at;
	}
}*/
