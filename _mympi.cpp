/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include"_ctrl.h"
#ifdef USE_MPI



#include<iostream>
#include<vector>
#include <algorithm>
#include "mpi.h"
#include"_class.h"
#include"_class2.h"
#include "_variable.h"
#include "_arp2_3.h"
#include "_periodic.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include "_wall_boundary.h"
#include "_function.h"
void struct_finalize();

//template<class T ,int N_ARRAY> inline std::ostream& operator<<(std::ostream& stream,const T (&in)[N_ARRAY]){
//	stream<<in[0];
//	for(int i=1;i<N_ARRAY;++i)stream<<'\t'<<in[i];
//	return stream;
//}





using _ARP2_3::_arp2_3;
using _ARP2_3::arp2_3;
using WALL::WALL_Z;
namespace mympi{

int pair_pre;
int pair_bck;

//_ptcl* ptcl_sender;

struct _pointer_list{
public:
	_ptcl* p;
	_ARP2_3::_arp2_3* ar;
	_list* li;
};
//vector<_pointer_list> pointer_list;
_pointer_list* pointer_list;
packer bit_sender(10485760);//10485760=1024^2*10=10Mばいと
packer bit_recver(10485760);

_WALL_Z* wall_tmp;
int* displs;
int* displs_sub;
//void debug(const int myrank){
//// 	std::stringstream ss;
//// 	std::ofstream fout;
//// 	ss.str("");
//// 	ss<<"./output"<<std::setw(4)<<std::setfill('0')<<N<<"/"<<"ptcl"<<myrank<<".xls";
//// 	fout.open(ss.str().c_str(),std::ios::out );
//// 	fout<<"step=\t"<<step<<"\tt=\t"<<t<<std::endl;
//// 	fout<<"\tlocation\t\t\tforce\t\t\tfila_number\tnumber"<<endl;
//// 	for( int i = 0; i < Number_Ptcl; i++){
//// 		_ptcl* p = &ptcl[i];
//// 		fout<<i<<'\t'<<p->loc<<'\t'<<p->force<<'\t'
//// 			<<p->cnf->feature<<'\t'
//// //			<<p->num
//// 			<<std::endl;
//// 	}
//// 	fout.close();
//// 	fout.clear();
//
//	
//	
//	
//	
//	MPI_Finalize();
//	struct_finalize();
//	exit(0);
//}
//

void mympi_init(int& argc,char**& argv,int& myrank,int& numprocs){
	MPI_Init(&argc,&argv);						// 初期化
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);	// プロセス数の取得
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);		// マイランク取得
	
	if(!(myrank))pair_pre=numprocs-1;
	else pair_pre=myrank-1;
	if(myrank==numprocs-1) pair_bck=0;
	else pair_bck=myrank+1;
	
//	pointer_list.resize(numprocs);
	pointer_list=new _pointer_list[numprocs];
	
	send_list_ptcl_x_pre.reserve( Number_Ptcl/numprocs );
	send_list_ptcl_x_pre_inner.reserve( Number_Ptcl/numprocs );
	send_list_ptcl_x_bck.reserve( Number_Ptcl/numprocs );
	send_list_arp_x_pre.reserve( Number_Arp2_3/numprocs );
	send_list_arp_x_pre_inner.reserve( Number_Arp2_3/numprocs );
	send_list_arp_x_bck.reserve( Number_Arp2_3/numprocs );
	
	init_genrand(myrank + numprocs*3);
	wall_tmp = new _WALL_Z[2*numprocs];
	displs = new int[numprocs];
	displs_sub = new int[numprocs];
	return;
}

void mympi_finalize(){
	MPI_Finalize();
	delete[] wall_tmp;
	delete[] displs;
	delete[] displs_sub;
}


void mympi_ptcl_init(const int myrank,const int numprocs){
	{
		_pointer_list pointer_list_send;
		pointer_list_send.p=&ptcl[0];
		pointer_list_send.li=list;
	//	pointer_list_send.b=bond;
		pointer_list_send.ar=_ARP2_3::arp2_3;
		
		MPI_Allgather(&pointer_list_send,sizeof(_pointer_list),MPI_BYTE,pointer_list,sizeof(_pointer_list),MPI_BYTE,MPI_COMM_WORLD);
//		MPI_Allgather(&pointer_list_send,sizeof(_pointer_list),MPI_BYTE,&pointer_list,sizeof(_pointer_list),MPI_BYTE,MPI_COMM_WORLD);
	}
	MPI_Bcast(ptcl.pt,((Number_Ptcl_Active+Number_Ptcl_Spare_Dens)*sizeof(_ptcl)), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_ARP2_3::arp2_3,Number_Arp2_3*sizeof(_ARP2_3::_arp2_3), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ptcl.con,((Number_Ptcl_Active+Number_Ptcl_Spare_Dens)*sizeof(_config)), MPI_BYTE, 0, MPI_COMM_WORLD);
	if(myrank){
	for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl_Active+Number_Ptcl_Spare_Dens];++p){
		p->cnf=&ptcl.con[p->index];
		p->cnf->ptcl=p;
		if(p->cnf->plus) p->cnf->plus =&ptcl[ (p->cnf->plus  - pointer_list[0].p)];
		if(p->cnf->minus)p->cnf->minus=&ptcl[ (p->cnf->minus - pointer_list[0].p)];
		if(p->cnf->edge_arp2_3)p->cnf->edge_arp2_3=arp2_3 + (p->cnf->edge_arp2_3 - pointer_list[0].ar);
		if(p->cnf->end_arp2_3) p->cnf->end_arp2_3 =arp2_3 + (p->cnf->end_arp2_3  - pointer_list[0].ar);
	}}
#ifdef WALL_BOUNDARY
	MPI_Bcast(WALL::bind_wall,WALL::bind_size*sizeof(_bind_wall), MPI_BYTE, 0, MPI_COMM_WORLD);
	if(myrank){
	for(_bind_wall* p=WALL::bind_wall;p!=WALL::bind_wall+WALL::bind_size;++p){
		if(p->index)p->index =&ptcl[(p->index  - pointer_list[0].p)];
	}}
#endif
	if(myrank){
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl_Active+Number_Ptcl_Spare_Dens];++p){
			_list* li_sub=p->list_calc();
			p->add(li_sub);
			li_sub->add(p);
		}
		for(_ARP2_3::_arp2_3* p=_ARP2_3::arp2_3; p<_ARP2_3::arp2_3+Number_Arp2_3;++p){
			_list* li_sub=p->list_calc();
			p->add(li_sub);
			li_sub->add(p);
		}
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl_Active+Number_Ptcl_Spare_Dens];++p){
			if((p->cnf->plus) && !(p->cnf->minus) && !(p->cnf->end_arp2_3))minus_end->add(p);
			else if(!(p->cnf->plus) && (p->cnf->minus))plus_end->add(p);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return;
}
//#if NZ>1
void mympi_send_Z_pre(const int myrank,const int numprocs){
	send_list_ptcl_z_pre.clear();
	send_list_ptcl_z_pre_inner.clear();
	send_list_ptcl_z_pre_inner_new.clear();
	send_list_ptcl_z_bck.clear();
	send_list_ptcl_z_bck_inner.clear();
	send_list_ptcl_z_bck_inner_new.clear();
	send_list_arp_z_pre.clear();
	send_list_arp_z_pre_inner.clear();
	send_list_arp_z_pre_inner_new.clear();
	send_list_arp_z_bck.clear();
	send_list_arp_z_bck_inner.clear();
	send_list_arp_z_bck_inner_new.clear();
	
	{
		{
			int z[4]={(calcArea.nz!=0   )?(calcArea.z_pre-1):list_Z-1,
					  (calcArea.nz!=NY-1)?(calcArea.z_bck)  :0,
					  calcArea.z_pre,
					  calcArea.z_bck-1
				};
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z[0] + y*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z[1] + y*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z[2] + y*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z[3] + y*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck_inner.push_back(arp);}
				}
			}}
if(NX>1){
			int x[2]={(calcArea.nx!=0   )?(calcArea.x_pre-1):list_X-1,
					  (calcArea.nx!=NX-1)?(calcArea.x_bck)  :0,
				};
			for(int i(0);i<2;++i){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z[0] + y*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z[1] + y*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z[2] + y*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z[3] + y*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck_inner.push_back(arp);}
				}
			}}
}
if (NY>1){
			int y[2]={(calcArea.ny!=0   )?(calcArea.y_pre-1):list_Y-1,
					  (calcArea.ny!=NY-1)?(calcArea.y_bck)  :0,
				};
			for(int i(0);i<2;++i){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z[0] + y[i]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z[1] + y[i]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z[2] + y[i]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z[3] + y[i]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck_inner.push_back(arp);}
				}
			}}
if(NX>1){
			for(int i(0);i<2;++i){
			for(int k(0);k<2;++k){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				int x[2]={(calcArea.nx!=0   )?(calcArea.x_pre-1):list_X-1,
						  (calcArea.nx!=NX-1)?(calcArea.x_bck)  :0,
					};
				{//のりしろ前外側
					li=list+(z[0] + y[i]*list_Z + x[k]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z[1] + y[i]*list_Z + x[k]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z[2] + y[i]*list_Z + x[k]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z[3] + y[i]*list_Z + x[k]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_z_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_z_bck_inner.push_back(arp);}
				}
			}}
}
}
		}
	}
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_z_pre.size() + send_list_arp_z_pre.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_pre.begin();p!=send_list_ptcl_z_pre.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_pre.begin();p!=send_list_arp_z_pre.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&(index));
		bit_sender.pack(&((*p)->force));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_z,myrank+4*numprocs             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_z,calcArea.pair_bck_z+4*numprocs,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_z_bck.size() + send_list_arp_z_bck.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_bck.begin();p!=send_list_ptcl_z_bck.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_bck.begin();p!=send_list_arp_z_bck.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->force));
	}
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_z, myrank+5*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_z, calcArea.pair_pre_z+5*numprocs ,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NZ>1
//#if NY>1
void mympi_send_Y_pre(const int myrank,const int numprocs){
	send_list_ptcl_y_pre.clear();
	send_list_ptcl_y_pre_inner.clear();
	send_list_ptcl_y_pre_inner_new.clear();
	send_list_ptcl_y_bck.clear();
	send_list_ptcl_y_bck_inner.clear();
	send_list_ptcl_y_bck_inner_new.clear();
	send_list_arp_y_pre.clear();
	send_list_arp_y_pre_inner.clear();
	send_list_arp_y_pre_inner_new.clear();
	send_list_arp_y_bck.clear();
	send_list_arp_y_bck_inner.clear();
	send_list_arp_y_bck_inner_new.clear();
	{
		{
			int y[4]={(calcArea.ny!=0   )?(calcArea.y_pre-1):list_Y-1,
					  (calcArea.ny!=NY-1)?(calcArea.y_bck)  :0,
					  calcArea.y_pre,
					  calcArea.y_bck-1
				};
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z + y[0]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z + y[1]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z + y[2]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z + y[3]*list_Z + x*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_bck_inner.push_back(arp);}
				}
			}}
if(NX>1){
			int x[2]={(calcArea.nx!=0   )?(calcArea.x_pre-1):list_X-1,
					  (calcArea.nx!=NX-1)?(calcArea.x_bck)  :0,
				};
			for(int i(0);i<2;++i){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z + y[0]*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z + y[1]*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z + y[2]*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z + y[3]*list_Z + x[i]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_y_bck_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_y_bck_inner.push_back(arp);}
				}
			}}
}
		}
	}
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_y_pre.size() + send_list_arp_y_pre.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_pre.begin();p!=send_list_ptcl_y_pre.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_pre.begin();p!=send_list_arp_y_pre.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&(index));
		bit_sender.pack(&((*p)->force));
	}
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_y,myrank+2*numprocs             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_y,calcArea.pair_bck_y+2*numprocs,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_y_bck.size() + send_list_arp_y_bck.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_bck.begin();p!=send_list_ptcl_y_bck.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_bck.begin();p!=send_list_arp_y_bck.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->force));
	}
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_y, myrank+3*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_y, calcArea.pair_pre_y+3*numprocs ,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NY>1
//#if NX>1
void mympi_send_X_pre(const int myrank,const int numprocs){
	send_list_ptcl_x_pre.clear();
	send_list_ptcl_x_pre_inner.clear();
	send_list_ptcl_x_pre_inner_new.clear();
	send_list_ptcl_x_bck.clear();
	send_list_ptcl_x_bck_inner.clear();
	send_list_ptcl_x_bck_inner_new.clear();
	send_list_arp_x_pre.clear();
	send_list_arp_x_pre_inner.clear();
	send_list_arp_x_pre_inner_new.clear();
	send_list_arp_x_bck.clear();
	send_list_arp_x_bck_inner.clear();
	send_list_arp_x_bck_inner_new.clear();
	
	{
		{
//			int x[4]={(pair_pre+1)*calc_length_X-1,
//					  pair_bck*calc_length_X,
//					  myrank*calc_length_X,
//					  (myrank+1)*calc_length_X-1
//					  };

			int x[4]={(calcArea.nx!=0   )?(calcArea.x_pre-1):list_X-1,
					  (calcArea.nx!=NX-1)?(calcArea.x_bck)  :0,
					  calcArea.x_pre,
					  calcArea.x_bck-1
				};
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				_ptcl* pa;
				_ARP2_3::_arp2_3* arp;
				_list* li;
				{//のりしろ前外側
					li=list+(z + y*list_Z + x[0]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_x_pre.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_x_pre.push_back(arp);}
				}{//のりしろ後外側
					li=list+(z + y*list_Z + x[1]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_x_bck.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_x_bck.push_back(arp);}
				}{//のりしろ内
					li=list+(z + y*list_Z + x[2]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){send_list_ptcl_x_pre_inner.push_back(pa);}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_x_pre_inner.push_back(arp);}
				}{//のりしろ内側
					li=list+(z + y*list_Z + x[3]*(list_Y*list_Z));
					for(int j=0;pa=li->ptcl_index[j];j++){
						send_list_ptcl_x_bck_inner.push_back(pa);
						if(pa->index>=Active_Ptcl){
							cout<<myrank<<" "<<pa->index<<" list_err "<<endl;
//							cout<<li-list<<" "<<(li-list)/(list_Y*list_Z)<<" "<<((li-list)%(list_Y*list_Z))/list_Z<<" "<<(li-list)%list_Z<<endl;
//							cout<<li->num<<" "<<j<<endl;
							
						}
					}
					for(int j=0;arp=li->arp_index[j];j++){send_list_arp_x_bck_inner.push_back(arp);}
				}
			}}
		}
	}
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_x_pre.size() + send_list_arp_x_pre.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_pre.begin();p!=send_list_ptcl_x_pre.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_pre.begin();p!=send_list_arp_x_pre.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&(index));
		bit_sender.pack(&((*p)->force));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_x,myrank             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_x,calcArea.pair_bck_x,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int send_size_ptcl=send_list_ptcl_x_bck.size() + send_list_arp_x_bck.size();
		bit_sender.pack(&send_size_ptcl);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_bck.begin();p!=send_list_ptcl_x_bck.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->force));
		if((*p)->cnf->feature)bit_sender.pack(&((*p)->adhesion_frc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_bck.begin();p!=send_list_arp_x_bck.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->force));
	}
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_x, myrank+numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_x, calcArea.pair_pre_x+numprocs ,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				ptcl[index_p].force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				if(ptcl[index_p].cnf->feature){
					ptcl[index_p].adhesion_frc+=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				(arp2_3+(index_p-Number_Ptcl))->force+=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NX>1
void mympi_send_pre(const int myrank,const int numprocs){
if(NZ>1){
	mympi_send_Z_pre(myrank,numprocs);
}
if(NY>1){
	mympi_send_Y_pre(myrank,numprocs);
}
if(NX>1){
	mympi_send_X_pre(myrank,numprocs);
}
#ifdef WALL_BOUNDARY
	MPI_Gather(WALL::WALL_Z, 2*sizeof(_WALL_Z), MPI_BYTE, wall_tmp,2*sizeof(_WALL_Z), MPI_BYTE,0, MPI_COMM_WORLD);
	static double* fm_ff_gather(new double[numprocs*2]);
	double fm_ff[2]={force_wall_mono ,force_wall_fila};
	MPI_Gather(fm_ff, 2, MPI_DOUBLE, fm_ff_gather, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(!myrank){
		for(int j=1;j<numprocs;++j){
			for(int i=0;i<2;i++){
				WALL_Z[i].LJ  +=wall_tmp[j*2+i].LJ;
				WALL_Z[i].adhe+=wall_tmp[j*2+i].adhe;
				WALL_Z[i].gaus+=wall_tmp[j*2+i].gaus;
#ifdef SRD                  
				WALL_Z[i].gaus_srd+=wall_tmp[j*2+i].gaus_srd;
#endif                      
			}
			force_wall_mono+=fm_ff_gather[j*2];
			force_wall_fila+=fm_ff_gather[j*2+1];
//			cout<<fm_ff_gather[j*2]<<'\t'<<fm_ff_gather[j*2+1]<<endl;
		}
	}
// 	if(!myrank){
// 	if(step>10000){
// 		
// 		for(int j=0;j<numprocs;++j){
// 		cout<<step<<' '<<t<<' '<<myrank<<" "<<j;
// 		for(int i=0;i<2;i++){
// 			cout<<wall_tmp[j*2+i].loc.z<<' '<<wall_tmp[j*2+i].LJ<<" "<<wall_tmp[j*2+i].adhe<<" "<<wall_tmp[j*2+i].gaus<<endl;
// 		}}
// 		
// 	}
// 	}
#endif
}
//#if NZ>1
void mympi_send_Z_bck(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_z_bck_inner.size() + send_list_arp_z_bck_inner.size() + send_list_ptcl_z_bck_inner_new.size() + send_list_arp_z_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_bck_inner.begin();p!=send_list_ptcl_z_bck_inner.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_bck_inner.begin();p!=send_list_arp_z_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_bck_inner_new.begin();p!=send_list_ptcl_z_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_bck_inner_new.begin();p!=send_list_arp_z_bck_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_z, myrank+10*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_z, calcArea.pair_pre_z+10*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3) {
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	
	MPI_Wait(&srequest, &sstatus);
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_z_pre_inner.size() + send_list_arp_z_pre_inner.size() + send_list_ptcl_z_pre_inner_new.size() + send_list_arp_z_pre_inner_new.size();
		bit_sender.pack(&sendsize);
// 		std::cout<<myrank<<" "<<sendsize<<std::endl;
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_pre_inner.begin();p!=send_list_ptcl_z_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_pre_inner.begin();p!=send_list_arp_z_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_pre_inner_new.begin();p!=send_list_ptcl_z_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_pre_inner_new.begin();p!=send_list_arp_z_pre_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_z, myrank+11*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_z, calcArea.pair_bck_z+11*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3){
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NZ>1
//#if NY>1
void mympi_send_Y_bck(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_y_bck_inner.size() + send_list_arp_y_bck_inner.size() + send_list_ptcl_y_bck_inner_new.size() + send_list_arp_y_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_bck_inner.begin();p!=send_list_ptcl_y_bck_inner.end();p++){
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_bck_inner.begin();p!=send_list_arp_y_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_bck_inner_new.begin();p!=send_list_ptcl_y_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_bck_inner_new.begin();p!=send_list_arp_y_bck_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_y, myrank+8*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_y, calcArea.pair_pre_y+8*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3) {
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	
	MPI_Wait(&srequest, &sstatus);
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_y_pre_inner.size() + send_list_arp_y_pre_inner.size() + send_list_ptcl_y_pre_inner_new.size() + send_list_arp_y_pre_inner_new.size();
		bit_sender.pack(&sendsize);
// 		std::cout<<myrank<<" "<<sendsize<<std::endl;
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_pre_inner.begin();p!=send_list_ptcl_y_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_pre_inner.begin();p!=send_list_arp_y_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_pre_inner_new.begin();p!=send_list_ptcl_y_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_pre_inner_new.begin();p!=send_list_arp_y_pre_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_y, myrank+9*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_y, calcArea.pair_bck_y+9*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3){
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NY>1
//#if NX>1
void mympi_send_X_bck(const int myrank,const int numprocs){

//std::ofstream fout;
//static bool flag(true);
//stringstream ss;
//ss.str("");
//ss<<myrank<<".dat";
//if(flag){
//	fout.open(ss.str().c_str(),std::ios::out);
//	flag=false;
//}else{
//	fout.open(ss.str().c_str(),std::ios::out | std::ios::app);
//}
//std::ofstream fout1;
//static bool flag1(true);
//stringstream ss1;
//ss1.str("");
//ss<<myrank<<"1.dat";
//if(flag1){
//	fout1.open(ss.str().c_str(),std::ios::out);
//	flag1=false;
//}else{
//	fout1.open(ss.str().c_str(),std::ios::out | std::ios::app);
//}
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_x_bck_inner.size() + send_list_arp_x_bck_inner.size() + send_list_ptcl_x_bck_inner_new.size() + send_list_arp_x_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_bck_inner.begin();p!=send_list_ptcl_x_bck_inner.end();p++){
//				if(((*p)->index)>=Active_Ptcl){
//					cout<<myrank<<" "<<((*p)->index)<<" ptcl_index_Err "<<endl;
////					exit(0);
//				}
		bit_sender.pack(&((*p)->index));
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_bck_inner.begin();p!=send_list_arp_x_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_bck_inner_new.begin();p!=send_list_ptcl_x_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_bck_inner_new.begin();p!=send_list_arp_x_bck_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
//fout<<"mpi send 0  "<<myrank<<" "<<numprocs<<" "<<calcArea.pair_bck_x<<" "<<calcArea.pair_pre_x<<" "<<bit_sender.size()<<" "<<bit_recver.size()<<endl;
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_x, myrank+6*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_x, calcArea.pair_pre_x+6*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
//fout<<"mpi send 1  "<<myrank<<" "<<numprocs<<" "<<calcArea.pair_bck_x<<" "<<calcArea.pair_pre_x<<" "<<bit_sender.size()<<" "<<bit_recver.size()<<endl;
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
//fout1<<*send_size_ptcl<<" send_size_ptcl"<<endl;
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
//fout1<<index_p<<endl;
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
//				if(index_p>=Active_Ptcl){
//					cout<<index_p<<" ptcl_index_Err "<<endl;
//					exit(0);
//				}
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3) {
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	
	MPI_Wait(&srequest, &sstatus);
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_x_pre_inner.size() + send_list_arp_x_pre_inner.size() + send_list_ptcl_x_pre_inner_new.size() + send_list_arp_x_pre_inner_new.size();
		bit_sender.pack(&sendsize);
// 		std::cout<<myrank<<" "<<sendsize<<std::endl;
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_pre_inner.begin();p!=send_list_ptcl_x_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_pre_inner.begin();p!=send_list_arp_x_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_pre_inner_new.begin();p!=send_list_ptcl_x_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
			bit_sender.pack(&(*p)->quat_pre);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_pre_inner_new.begin();p!=send_list_arp_x_pre_inner_new.end();p++){
		int index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
		bit_sender.pack(&((*p)->loc_pre));
	}
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_x, myrank+7*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_x, calcArea.pair_bck_x+7*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < Number_Ptcl+Number_Arp2_3) {
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else if(index_p < 2*Number_Ptcl+Number_Arp2_3) {
				p=&ptcl[(index_p-Number_Ptcl-Number_Arp2_3)];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p < 2*Number_Ptcl+2*Number_Arp2_3){
				arp=arp2_3+(index_p-2*Number_Ptcl-Number_Arp2_3);
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				_list* li=arp->list_calc();
				arp->force.IN(0.0,0.0,0.0);
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else {
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NX>1
void mympi_send_bck(const int myrank,const int numprocs){
if(NX>1){
	mympi_send_X_bck(myrank, numprocs);
}
if(NY>1){
	mympi_send_Y_bck(myrank, numprocs);
}
if(NZ>1){
	mympi_send_Z_bck(myrank, numprocs);
}
#ifdef WALL_BOUNDARY
	MPI_Bcast(WALL::WALL_Z,2*sizeof(_WALL_Z), MPI_BYTE,0, MPI_COMM_WORLD);
#endif
}

//#if NZ>1
void mympi_send_Z_bck_second(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_z_bck_inner.size() + send_list_arp_z_bck_inner.size() + send_list_ptcl_z_bck_inner_new.size() + send_list_arp_z_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_bck_inner.begin();p!=send_list_ptcl_z_bck_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_bck_inner_new.begin();p!=send_list_ptcl_z_bck_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_bck_inner.begin();p!=send_list_arp_z_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_bck_inner_new.begin();p!=send_list_arp_z_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_z ,myrank+12*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_z ,calcArea.pair_pre_z+12*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_z_pre_inner.size() + send_list_arp_z_pre_inner.size() + send_list_ptcl_z_pre_inner_new.size() + send_list_arp_z_pre_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_pre_inner.begin();p!=send_list_ptcl_z_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_z_pre_inner_new.begin();p!=send_list_ptcl_z_pre_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_pre_inner.begin();p!=send_list_arp_z_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_z_pre_inner_new.begin();p!=send_list_arp_z_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_z,myrank+13*numprocs             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_z,calcArea.pair_bck_z+13*numprocs,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else { 
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NZ>1
//
//#if NY>1
void mympi_send_Y_bck_second(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_y_bck_inner.size() + send_list_arp_y_bck_inner.size() + send_list_ptcl_y_bck_inner_new.size() + send_list_arp_y_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_bck_inner.begin();p!=send_list_ptcl_y_bck_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_bck_inner_new.begin();p!=send_list_ptcl_y_bck_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_bck_inner.begin();p!=send_list_arp_y_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_bck_inner_new.begin();p!=send_list_arp_y_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_y ,myrank+12*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_y ,calcArea.pair_pre_y+12*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_y_pre_inner.size() + send_list_arp_y_pre_inner.size() + send_list_ptcl_y_pre_inner_new.size() + send_list_arp_y_pre_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_pre_inner.begin();p!=send_list_ptcl_y_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_y_pre_inner_new.begin();p!=send_list_ptcl_y_pre_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_pre_inner.begin();p!=send_list_arp_y_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_y_pre_inner_new.begin();p!=send_list_arp_y_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_y,myrank+13*numprocs             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_y,calcArea.pair_bck_y+13*numprocs,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else { 
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NY>1
//#if NX>1
void mympi_send_X_bck_second(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_x_bck_inner.size() + send_list_arp_x_bck_inner.size() + send_list_ptcl_x_bck_inner_new.size() + send_list_arp_x_bck_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_bck_inner.begin();p!=send_list_ptcl_x_bck_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_bck_inner_new.begin();p!=send_list_ptcl_x_bck_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_bck_inner.begin();p!=send_list_arp_x_bck_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_bck_inner_new.begin();p!=send_list_arp_x_bck_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	MPI_Status sstatus;
	MPI_Status rstatus;
	MPI_Request srequest;
	MPI_Request rrequest;
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_bck_x ,myrank+12*numprocs              ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_pre_x ,calcArea.pair_pre_x+12*numprocs ,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else{
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
	
	bit_sender.clear();
	bit_recver.clear();
	{
		int sendsize=send_list_ptcl_x_pre_inner.size() + send_list_arp_x_pre_inner.size() + send_list_ptcl_x_pre_inner_new.size() + send_list_arp_x_pre_inner_new.size();
		bit_sender.pack(&sendsize);
	}
	
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_pre_inner.begin();p!=send_list_ptcl_x_pre_inner.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_x_pre_inner_new.begin();p!=send_list_ptcl_x_pre_inner_new.end();p++){
		bit_sender.pack(&(*p)->index);
		bit_sender.pack(&((*p)->loc));
		if((*p)->cnf->feature){
			bit_sender.pack(&(*p)->adhesion_loc);
			bit_sender.pack(&(*p)->quat);
		}
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_pre_inner.begin();p!=send_list_arp_x_pre_inner.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	for(std::vector<_arp2_3*>::iterator p=send_list_arp_x_pre_inner_new.begin();p!=send_list_arp_x_pre_inner_new.end();p++){
		int index=(*p)->index+Number_Ptcl;
		bit_sender.pack(&index);
		bit_sender.pack(&((*p)->loc));
	}
	
	
	MPI_Isend(bit_sender[0],bit_sender.size()   ,MPI_BYTE,calcArea.pair_pre_x,myrank+13*numprocs             ,MPI_COMM_WORLD,&srequest);
	MPI_Irecv(bit_recver[0],bit_recver.maxsize(),MPI_BYTE,calcArea.pair_bck_x,calcArea.pair_bck_x+13*numprocs,MPI_COMM_WORLD,&rrequest);
	MPI_Wait(&rrequest, &rstatus);
	{
		int index=sizeof(int);
		int* send_size_ptcl=(int*)bit_recver[index];
		index+=sizeof(int);
		_ptcl* p;
		_arp2_3* arp;
		for(int i=0;i<*send_size_ptcl;i++){
			int index_p=*(int*)bit_recver[index];
			index+=sizeof(int);
			if(index_p<Number_Ptcl){
				p=&ptcl[index_p];
				p->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				p->loc_pre=p->loc;
				p->force.IN(0.0,0.0,0.0);
				if(p->cnf->feature){
					p->adhesion_loc=*(_vec<double>*)bit_recver[index];
					index+=sizeof(_vec<double>);
					p->adhesion_frc.IN(0.0,0.0,0.0);
					p->quat=*(_quat<double>*)bit_recver[index];
					index+=sizeof(_quat<double>);
					p->quat_pre=p->quat;
				}
				_list* li=p->list_calc();
				if(li!=p->cnf->li){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)!=2){		send_list_ptcl_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)!=0){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->cnf->li->myrank%3)!=2){		send_list_ptcl_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)!=0){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
					}
					p->cnf->li->rm(p,p->cnf->list_num);
					p->add(li);
					li->add(p);
				}
			}else if(index_p<Number_Ptcl+Number_Arp2_3){
				arp=(arp2_3+(index_p-Number_Ptcl));
				arp->loc=*(_vec<double>*)bit_recver[index];
				index+=sizeof(_vec<double>);
				arp->loc_pre=arp->loc;
				arp->force.IN(0.0,0.0,0.0);
				_list* li=arp->list_calc();
				if(li!=arp->li()){
					if(li->myrank!=-1){
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((arp->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(arp);
						}else if(((li->myrank/3)%3)==0 && ((arp->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(arp);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (arp->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(arp);
						}else if((li->myrank%3)==0 && (arp->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(arp);}
}
					}
					arp->li()->rm(arp,arp->list_num);
					arp->add(li);
					li->add(arp);
				}
			}else { 
				exit(0);
			}
		}
	}
	MPI_Wait(&srequest, &sstatus);
}
//#endif // NX>1
void mympi_send_bck_second(const int myrank,const int numprocs){
if(NX>1){
	mympi_send_X_bck_second(myrank, numprocs);
}
if(NY>1){
	mympi_send_Y_bck_second(myrank, numprocs);
}
if(NZ>1){
	mympi_send_Z_bck_second(myrank, numprocs);
}
#ifdef WALL_BOUNDARY
	MPI_Bcast(WALL::WALL_Z,2*sizeof(_WALL_Z), MPI_BYTE,0, MPI_COMM_WORLD);
#endif
	
	{
		double potential[3]={ene_pot,_ARP2_3::ene_pot,WALL::ene_pot};
		double potential_sub[3]={0.0, 0.0, 0.0};
		MPI_Reduce(potential,potential_sub,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if(!myrank){
			ene_pot			=potential_sub[0];
			_ARP2_3::ene_pot=potential_sub[1];
			WALL::ene_pot	=potential_sub[2];
		}
	}
	for(int i(0);i<Number_Ptcl;++i)filament_force_sub[i]=0.0;
	MPI_Reduce(filament_force,filament_force_sub,Number_Ptcl,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(!myrank){for(int i(0);i<Number_Ptcl;++i)filament_force[i]=filament_force_sub[i];}
	
	bit_sender.clear();
	bit_recver.clear();

}
inline bool nucleation(_ptcl* p_g1,_ptcl* p_g2,_ptcl* p_g1_sub,_ptcl* p_g2_sub ,const int myrank,const int numprocs,int i){//さてどうしたものか
	//nucleatがかぶったとき
	//方法1 myrankの若い方が優先
	//方法2 例外処理として0番機に選択させ，もう一度通信しなおす．
	//まぁ1のほうが簡単だよね．．
	//しかし，ねぇ？
	
	//polymerizationとかぶったとき→かぶらなくしました
	
	if(!(p_g1->cnf->feature) && !(p_g1->cnf->feature)){
		plus_end->polymerization_mono(p_g2,p_g1);
		p_g1->loc=p_g1_sub->loc;
		p_g2->loc=p_g2_sub->loc;
		p_g1->loc_pre=p_g1->loc;
		p_g2->loc_pre=p_g2->loc;
		p_g1->loc_pre=p_g1->loc;
		p_g2->loc_pre=p_g2->loc;
	p_g1->cnf->feature = 1;
	p_g2->cnf->feature = 1;
		p_g1->adhesion_loc=p_g1_sub->adhesion_loc;
		p_g2->adhesion_loc=p_g2_sub->adhesion_loc;
		p_g1->adhesion_loc_inite = p_g1->adhesion_loc;
		p_g2->adhesion_loc_inite = p_g2->adhesion_loc;
		p_g1->quat.IN(0.0 ,0.0 ,0.0);
		p_g1->quat_pre.IN(0.0 ,0.0 ,0.0);
		p_g2->quat.IN(0.0 ,0.0 ,0.0);
		p_g2->quat_pre.IN(0.0 ,0.0 ,0.0);
	}else{
		return false;
	}
	counter_pol++;
	return true;
}

inline bool polymerization(int i,_ptcl* p_f,_ptcl* p_g,_vec<double>& loc){
//取り合いになったらどうしようか
	plus_end->polymerization_fila(i,p_g);
	p_g->cnf->feature = 1;
	p_f->cnf->feature = 2;
	p_g->adhesion_loc = loc;
	p_g->adhesion_loc_inite = p_g->adhesion_loc;
	p_g->quat.IN(0.0 ,0.0 ,0.0);
	p_g->quat_pre.IN(0.0 ,0.0 ,0.0);

	counter_pol++;
	return true;
}
inline void depolymerization(_ptcl* p,int i){
	counter_depol++;
	minus_end->depolymerization(p,i);
	p->cnf->feature = 0;
	return;
}
inline void wall_boundary_bind(_ptcl* p,_bind_wall* b){
	p->cnf->bind= b - WALL::bind_wall;
	b->index=p;
	return;
}
inline void wall_boundary_debind(_ptcl* p,_bind_wall* b){
	b->index = NULL;
	p->cnf->bind = -1;
	return;
}
inline void arp2_3_edge_bind(_arp2_3* arp,_ptcl* p1){
	arp->edge_ptcl_in(p1->index);
	p1->cnf->edge_arp2_3=arp;//->index);
	return;
}
inline void arp2_3_end_bind(_arp2_3* arp,_ptcl* p,_vec<double>& loc){
	arp->end_ptcl_in(p->index);
	p->cnf->feature = 1;//no_fila;
	p->cnf->end_arp2_3=arp;//->index);
	plus_end->add(p);
		p->adhesion_loc = loc;
		p->adhesion_loc_inite = p->adhesion_loc;
		p->quat.IN(0.0 ,0.0 ,0.0);
		p->quat_pre.IN(0.0 ,0.0 ,0.0);
		
	return;
}
inline void arp2_3_debind(_arp2_3* arp){
	arp->edge_ptcl()->cnf->edge_arp2_3=NULL;
	if(arp->end_ptcl()){
		arp->end_ptcl()->cnf->end_arp2_3=NULL;
		if(arp->end_ptcl()->cnf->plus){
			minus_end->add(arp->end_ptcl());
		}else{
			plus_end->rm(arp->end_ptcl());
			arp->end_ptcl()->cnf->feature = 0;
		}
	}
	arp->edge_ptcl_in(-1);//=NULL;
	arp->end_ptcl_in(-1);//=NULL;
	return;
}
void mympi_conformation_change_nucleation(const int myrank,const int numprocs){
	{
		int size=bit_sender.size();
		MPI_Allgather(&size,1,MPI_INT,displs_sub,1,MPI_INT,MPI_COMM_WORLD);
		{
			displs[0]=0;
			for(int i=1;i<numprocs;++i){
				displs[i]=displs_sub[i-1]+displs[i-1];
			}
		}
		if((displs[numprocs-1]+displs_sub[numprocs-1])==sizeof(int)*numprocs)return;
		MPI_Allgatherv(bit_sender[0],bit_sender.size(),MPI_BYTE,bit_recver[0],displs_sub,displs,MPI_BYTE,MPI_COMM_WORLD);
	}
	int size_sub(0);
	int index(0);
// 	MPI_Allgather(bit_sender[0],bit_sender.size(),MPI_BYTE,bit_recver[0],bit_recver.maxsize(),MPI_BYTE,MPI_COMM_WORLD);
	for(int i=0;i<numprocs;++i){
		int size=*((int*)bit_recver[index]);
		size_sub+=size;
		index+=sizeof(int);
		if(i!=myrank){
		while(index<size_sub){
			_ptcl* pa=(_ptcl*)bit_recver[index];
			_ptcl* pb=(_ptcl*)bit_recver[index+sizeof(_ptcl)];
			index+=2*sizeof(_ptcl);
			if( !(nucleation(&ptcl[pa->index],&ptcl[pb->index],pa,pb,myrank,numprocs,i)) ){
				
			}
		}
		}
		index=size_sub;
	}
	bit_sender.clear();
	bit_recver.clear(index);
}
void mympi_conformation_change(const int myrank,const int numprocs){
	{
		int size=bit_sender.size();
		MPI_Allgather(&size,1,MPI_INT,displs_sub,1,MPI_INT,MPI_COMM_WORLD);
		{
			displs[0]=0;
			for(int i=1;i<numprocs;++i){
				displs[i]=displs_sub[i-1]+displs[i-1];
			}
		}
		if((displs[numprocs-1]+displs_sub[numprocs-1])==sizeof(int)*numprocs)return;
		MPI_Allgatherv(bit_sender[0],bit_sender.size(),MPI_BYTE,bit_recver[0],displs_sub,displs,MPI_BYTE,MPI_COMM_WORLD);
	}
	int size_sub(0);
	int index(0);
		for(int i=0;i<numprocs;++i){
			size_sub+=*((int*)bit_recver[index]);
			index+=sizeof(int);
			if(i!=myrank){
			while(index<size_sub){
				conf_ch* cc=(conf_ch*)bit_recver[index];
				index+=sizeof(conf_ch);
// cout<<myrank<<'\t'<<step<<'\t'<<i<<"\t"<<*cc<<" "<<*(int*)bit_recver[index]<<' '<<*(int*)bit_recver[index+sizeof(int)]<<'\n'<<endl;
				switch(*cc){
					case no:{
						break;
					}
// 					case pol_mono:{
// 						_ptcl* pa=(_ptcl*)bit_recver[index];
// 						_ptcl* pb=(_ptcl*)bit_recver[index+sizeof(_ptcl)];
// 						index+=2*sizeof(_ptcl);
// 						nucleation(ptcl[pa].index,ptcl+pb->index,pa,pb,myrank,numprocs,i);
// 						break;
// 					}
					case pol_fila:{//2
						int in=*(int*)bit_recver[index];
						int index_f=*(int*)bit_recver[index+sizeof(int)];
						int index_g=*(int*)bit_recver[index+2*sizeof(int)];
						_vec<double> adhesion_loc=*(_vec<double>*)bit_recver[index+3*sizeof(int)];
						index+=3*sizeof(int)+sizeof(_vec<double>);
						polymerization(in,&ptcl[index_f],&ptcl[index_g],adhesion_loc);
						break;
					}
					case depol:{//3
						int* index_p=(int*)bit_recver[index];
						int* in=(int*)bit_recver[index+sizeof(int)];
						index+=2*sizeof(int);
						depolymerization(&ptcl[*index_p],*in);
						break;
					}
					case wall_bind:{//5
						int* index_p=(int*)bit_recver[index];
						int* in=(int*)bit_recver[index+sizeof(int)];
						index+=2*sizeof(int);
						wall_boundary_bind(&ptcl[*index_p],WALL::bind_wall+*in);
						break;
					}
					case wall_debind:{//6
						int* index_p=(int*)bit_recver[index];
						int* in=(int*)bit_recver[index+sizeof(int)];
						index+=2*sizeof(int);
						wall_boundary_debind(&ptcl[*index_p],WALL::bind_wall+*in);
						break;
					}
					case arp_edge:{//7
						int* index_a=(int*)bit_recver[index];
						int* index_p=(int*)bit_recver[index+sizeof(int)];
						index+=2*sizeof(int);
						arp2_3_edge_bind(arp2_3+*index_a,&ptcl[*index_p]);
						break;
					}
					case arp_end:{//8
						int* index_a=(int*)bit_recver[index];
						int* index_p=(int*)bit_recver[index+sizeof(int)];
						_vec<double> adhesion_loc=*(_vec<double>*)bit_recver[index+2*sizeof(int)];
						index+=2*sizeof(int)+sizeof(_vec<double>);
						arp2_3_end_bind(arp2_3+*index_a,&ptcl[*index_p],adhesion_loc);
						break;
					}
					case arp_debind:{//8
						int* index_a=(int*)bit_recver[index];
						index+=2*sizeof(int);
						arp2_3_debind(arp2_3+*index_a);
						break;
					}
					default:{
						break;
					}
				}
			}
			}
			index=size_sub;
		}
		bit_sender.clear();
		bit_recver.clear(index);
}

void mympi_location_gather(const int myrank,const int numprocs){
	bit_sender.clear();
	bit_recver.clear();
	if(myrank){
		{
			int size(0);
//			for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++)size+=(list+i)->num+(list+i)->num_arp;
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				int i(z + y*list_Z + x*(list_Y*list_Z));
				size+=(list+i)->num+(list+i)->num_arp;
			}}}
			bit_sender.pack(&size);
		}
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li = list+i;
			_ptcl* pa;
			_arp2_3* arp;
			for(int j=0;pa=li->ptcl_index[j];j++){
				bit_sender.pack(&pa->index);
				bit_sender.pack(&pa->loc);
				if(pa->cnf->feature)	bit_sender.pack(&pa->quat);
			}
			for(int j=0;arp=li->arp_index[j];j++){
				int i=arp->index+Number_Ptcl;
				bit_sender.pack(&i);
				bit_sender.pack(&arp->loc);
			}
		}}}
	}
	{
		int size=bit_sender.size();
		MPI_Gather(&size,1,MPI_INT,displs_sub,1,MPI_INT,0,MPI_COMM_WORLD);
		if(!myrank){
			displs[0]=0;
			for(int i=1;i<numprocs;++i){
				displs[i]=displs_sub[i-1]+displs[i-1];
			}
		}
		MPI_Gatherv(bit_sender[0],bit_sender.size(),MPI_BYTE,bit_recver[0],displs_sub,displs,MPI_BYTE,0,MPI_COMM_WORLD);
	}
	
	
	if(!myrank){
		int index(sizeof(int));
		int size(0);
		int index_p;
		
		for(int i=1;i<numprocs;i++){
//cout<<index<<" "<<displs[i]<<" "<<*((int*)bit_recver[index])<<endl;
			index+=sizeof(int);
			size=*((int*)bit_recver[index]);
			index+=sizeof(int);
			for(int j=0;j<size;j++){
				index_p=*((int*)bit_recver[index]);
//cout<<i<<'\t'<<j<<'\t'<<index_p<<'\t'<<index<<endl;
				if(index_p<Number_Ptcl){
					_ptcl* p=&ptcl[index_p];
//cout<<p->cnf->feature<<endl;
					p->loc=*((_vec<double>*)bit_recver[index+sizeof(int)]);
					if(p->cnf->feature){
						p->quat=*((_quat<double>*)bit_recver[index+sizeof(int)+sizeof(_vec<double>)]);
						p->adhesion_loc= _quaternion::rot(p->quat,p->adhesion_loc_inite);
						index+=sizeof(_quat<double>);
					}
				}else{
					(arp2_3+index_p-Number_Ptcl)->loc=*((_vec<double>*)bit_recver[index+sizeof(int)]);
				}
				index+=sizeof(int)+sizeof(_vec<double>);
			}
		}
		bit_recver.clear(index);
	}
	bit_sender.clear();
	bit_recver.clear();

//	for(int i(0);i<Number_Ptcl;++i)filament_force_sub[i]=0.0;
//	MPI_Reduce(filament_force,filament_force_sub,Number_Ptcl,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//	if(!myrank){for(int i(0);i<Number_Ptcl;++i)filament_force[i]=filament_force_sub[i];}
	
	{
		static _output_PlusEndPosition* sub=new _output_PlusEndPosition[numprocs];
		MPI_Gather(&output_PlusEndPosition,sizeof(_output_PlusEndPosition),MPI_BYTE,sub,sizeof(_output_PlusEndPosition),MPI_BYTE,0,MPI_COMM_WORLD);
		output_PlusEndPosition.clear();
		if(!myrank){
			for(int i(0);i<numprocs;i++){
				output_PlusEndPosition.app(sub[i]);
			}
		}
	//	delete[] sub;
	}
//	{
//		static _PTCL_Ene_Output* sub=new _PTCL_Ene_Output[numprocs];
//		MPI_Gather(&PTCL_Ene_Output,sizeof(_PTCL_Ene_Output),MPI_BYTE,sub,sizeof(_PTCL_Ene_Output),MPI_BYTE,0,MPI_COMM_WORLD);
//		PTCL_Ene_Output.clear();
//		if(!myrank){
//			for(int i(0);i<numprocs;i++){
//				cout<<" "<<sub[i].ptcl_ene_pot2<<" "<<PTCL_Ene_Output.ptcl_ene_pot2<<endl;
//				PTCL_Ene_Output.app(sub[i],i);
//			}
//		}
//		//delete[] sub;
//	}
	
	
	if(step%1000==0){
		int size=PolymerizationInfo.size()*sizeof(_PolymerizationInfoSub);
		MPI_Allgather(&size,1,MPI_INT,displs_sub,1,MPI_INT,MPI_COMM_WORLD);
			displs[0]=0;
			for(int i=1;i<numprocs;++i){
				displs[i]=displs_sub[i-1]+displs[i-1];
			}
		
		if((displs[numprocs-1]+displs_sub[numprocs-1])!=0){
			_PolymerizationInfoSub *sub=new _PolymerizationInfoSub[(displs[numprocs-1]+displs_sub[numprocs-1])/sizeof(_PolymerizationInfoSub)+1];//new Type[1]は怖い
			
			//MPI_Gatherv(PolymerizationInfo.buf.begin(),size,MPI_BYTE,sub,displs_sub,displs,MPI_BYTE,0,MPI_COMM_WORLD);//PGIではコンパイルできる 2010.11.16
			MPI_Gatherv(&PolymerizationInfo.buf[0],size,MPI_BYTE,sub,displs_sub,displs,MPI_BYTE,0,MPI_COMM_WORLD);//Intel用に変更 2010.11.16
			PolymerizationInfo.clear();
			
			if(!myrank){
				for(int i(0);i<(displs[numprocs-1]+displs_sub[numprocs-1])/sizeof(_PolymerizationInfoSub);++i){
					PolymerizationInfo.add(sub[i]);
				}
			}
			delete[] sub;
		}
		
		{
			_flow_cell_element* sub_abs = new _flow_cell_element[numprocs * flow.size.x * flow.size.y * flow.size.z];
			_flow_cell_element* sub_rel = new _flow_cell_element[numprocs * flow.size.x * flow.size.y * flow.size.z];
			_flow_cell_element* sub_abs_mono = new _flow_cell_element[numprocs * flow.size.x * flow.size.y * flow.size.z];
			_flow_cell_element* sub_rel_mono = new _flow_cell_element[numprocs * flow.size.x * flow.size.y * flow.size.z];
			MPI_Gather(flow.abs,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,
						sub_abs,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,0,MPI_COMM_WORLD);
			MPI_Gather(flow.rel,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,
						sub_rel,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,0,MPI_COMM_WORLD);
			MPI_Gather(flow.abs_mono,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,
						sub_abs_mono,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,0,MPI_COMM_WORLD);
			MPI_Gather(flow.rel_mono,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,
						sub_rel_mono,sizeof(_flow_cell_element)*flow.size.x*flow.size.y*flow.size.z,MPI_BYTE,0,MPI_COMM_WORLD);
			double tmp_t=flow.start_step;
			flow.clear();
			
			if(!myrank){
				flow.start_step=tmp_t;
				for(int i(0);i<numprocs;i++){
					for(int j(0);j<flow.size.x*flow.size.y*flow.size.z;j++){
						flow.abs[j].add(sub_abs[i*flow.size.x*flow.size.y*flow.size.z + j]);
						flow.rel[j].add(sub_rel[i*flow.size.x*flow.size.y*flow.size.z + j]);
						flow.abs_mono[j].add(sub_abs_mono[i*flow.size.x*flow.size.y*flow.size.z + j]);
						flow.rel_mono[j].add(sub_rel_mono[i*flow.size.x*flow.size.y*flow.size.z + j]);
					}
				}
				
			}
			delete[] sub_abs;
			delete[] sub_rel;
			delete[] sub_abs_mono;
			delete[] sub_rel_mono;
		//	delete[] sub;
		}	
		
	}
	
	
	
}
void mympi_flow_wallinf_bcast(const int myrank,const int numprocs){
	MPI_Bcast(&flow.wall_loc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&flow.wall_vel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


//void mympi_replenish_ptcl(const int myrank,const int numprocs){
//	int size(bit_sender.size());
//	MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
//	
//	if(size!=sizeof(int)){
//		MPI_Bcast(bit_sender[0],size,MPI_BYTE,0,MPI_COMM_WORLD);
//		if(myrank){
//			int Active_Ptcl	=*((int*)bit_sender[sizeof(int)]);
//			int margin		=*((int*)bit_sender[2*sizeof(int)]);
//			if(margin>0){
//				_vec<double> *loc=(_vec<double>*)bit_sender[3*sizeof(int)];
//				for(int i=Active_Ptcl;(i<Active_Ptcl+margin)&&(i<Number_Ptcl);i++){
//					_ptcl* p=&ptcl[i];
//					p->loc=*loc;
//					p->loc_pre=p->loc;
////					p->p_level = 1;
//					p->cnf->feature = 0;
//					p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
//					p->adhesion_loc_inite = p->adhesion_loc;
//					p->quat.IN(0.0,0.0,0.0);
//					p->quat_pre=p->quat;
//					p->index=p - ptcl;
//					_list* li_sub=list_calc(*loc);
//					p->add(li_sub);
//					li_sub->add(p);
//					loc++;
//				}
//			}
//		}
//	}
//	bit_sender.clear();
//	return;
//}
void mympi_replenish_ptcl(const int myrank,const int numprocs){
struct _change_index{
	int first;
	int sec;
	_vec<double> loc;
}*change_index;

	bit_sender.clear();
	if(!myrank){
		bit_sender.pack(&Active_Ptcl);
		bit_sender.pack(replenish_loc,Number_Ptcl_Spare_Dens);
		bit_sender.pack(replenish_index,Number_Ptcl_Spare_Dens);
		{
			int size(replenish_changelist.size());
			bit_sender.pack(&size);
			for(int i(0);i<size;i+=2){
				bit_sender.pack(&replenish_changelist[i]);
				bit_sender.pack(&replenish_changelist[i+1]);
				bit_sender.pack(&(ptcl[replenish_changelist[i]].loc));
			}
		}
//		int size(bit_sender.size());
	}
	int size(bit_sender.size());
	MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
	
	MPI_Bcast(bit_sender[0],size,MPI_BYTE,0,MPI_COMM_WORLD);
	if(myrank){
		int index(sizeof(int));
		int Active_Ptcl_new	=*((int*)bit_sender[index]);
				index+=sizeof(int);
		_vec<double>* loc=(_vec<double>*)bit_sender[index];
				index+=Number_Ptcl_Spare_Dens*sizeof(_vec<double>);
		int* loc_index=(int*)bit_sender[index];
				index+=Number_Ptcl_Spare_Dens*sizeof(int);
		int change_size=*((int*)bit_sender[index]);
				index+=sizeof(int);
		change_index=(_change_index*)bit_sender[index];
//				index+=change_size*sizeof(int)+change_size*sizeof(_vec<double>)/2;
		
//		for(int x=0;x<list_X;x++){//粒子数の確認
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<offset_impound;z++){
//			(list+(z + y*list_Z + x*(list_Y*list_Z)))->init();
//		}}}
//		
		for(int i(0);i<Number_Ptcl_Spare_Dens;++i){
//cout<<myrank<<"\t"<<loc_index[i]<<"\t"<<loc[i]<<endl;
//{
//	stringstream ss;ss.str("");
//	ss<<myrank<<".xls";
//	static bool flag(true);
//	if(flag){
//		std::ofstream fout(ss.str().c_str(),std::ios::out);
//		fout.close();
//		flag=false;
//	}
//	std::ofstream fout(ss.str().c_str(),std::ios::out |std::ios::app);
//fout<<step<<"\t"<<myrank<<"\t"<<Active_Ptcl_new<<'\t'<<i<<'\t'<<loc_index[i]<<"\t"<<loc[i]<<endl;
			if(loc_index[i]==-1)continue;
			_ptcl& p =(ptcl[loc_index[i]]);
			p.loc.IN(*(loc+i));
			p.loc_pre.IN(p.loc);
//				p.feature = 0;
			p.adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
			p.adhesion_loc_inite.IN(p.adhesion_loc);
			p.quat.IN(0.0,0.0,0.0);
			p.quat_pre.IN(p.quat);
			//list更新
			_list* li_sub=p.list_calc();
//if(step>=362)fout<<p.list_index<<'\t'<<p.index<<'\t'<<Active_Ptcl<<'\t'<<Active_Ptcl_new<<endl;
//			if(li_sub!=p.li()){
				if(p.index<Active_Ptcl)p.cnf->li->rm(&p,p.cnf->list_num);
//if(step>=362)fout<<p.list_index<<endl;
				p.add(li_sub);
//if(step>=362)fout<<p.list_index<<endl;
				li_sub->add(&p);
//if(step>=362)fout<<p.list_index<<endl;
//			}
//	fout.close();
//}
		}
//if(!(myrank-1))std::cout<<myrank<<" "<<step<<" "<<change_size<<" "<<Active_Ptcl_new<<std::endl;
		if(change_size){
			for(int i(0);i<change_size;i+=2){
				_ptcl* pa=&ptcl[change_index->first];
				_ptcl* pb=&ptcl[change_index->sec];
//				pa->loc.IN(change_index->loc);
//if(!(myrank-1))std::cout<<myrank<<" "<<change_index->first<<" "<<change_index->sec<<endl;				
				_list* li=list_calc(change_index->loc);
				pa->cnf->li->rm(pa,pa->cnf->list_num);
//				li->rm(pb,pb->list_num);
//				pa->Copy(pb);
				ptcl.Copy(*pa,*pb);
//				pb->cnf->feature=0;
				pa->add(li);
				li->add(pa);
				pa->loc.IN(change_index->loc);
				
				change_index++;
			}
		}
		if(Active_Ptcl_new < Active_Ptcl){for(_ptcl* p(&ptcl[Active_Ptcl_new]);p!=&ptcl[Active_Ptcl];++p){p->cnf->li->rm(p,p->cnf->list_num);} }
		Active_Ptcl =  Active_Ptcl_new;
	}
	bit_sender.clear();
	return;
}




















// void mympi_send_pre(int myrank,int numprocs){
// 	send_list_ptcl_pre.clear();
// 	send_list_ptcl_pre_inner.clear();
// 	send_list_ptcl_pre_inner_new.clear();
// 	send_list_ptcl_bck.clear();
// 	send_list_ptcl_bck_inner.clear();
// 	send_list_ptcl_bck_inner_new.clear();
// 	send_list_arp_pre.clear();
// 	send_list_arp_pre_inner.clear();
// 	send_list_arp_pre_inner_new.clear();
// 	send_list_arp_bck.clear();
// 	send_list_arp_bck_inner.clear();
// 	send_list_arp_bck_inner_new.clear();
// 	
// // 	std::cout<<t<<'\t'
// // 		<<send_list_ptcl_pre.capacity()<<'\t'
// // 		<<send_list_ptcl_pre_inner.capacity()<<'\t'
// // 		<<send_list_ptcl_bck.capacity()<<'\t'
// // 		<<send_list_arp_pre.capacity()<<'\t'
// // 		<<send_list_arp_pre_inner.capacity()<<'\t'
// // 		<<send_list_arp_bck.capacity()<<'\t'
// // 		<<ptcl_sender.capacity()<<'\t'
// // 		<<ptcl_recver.capacity()<<std::endl;
// 
// 	int x;
// 	if(myrank)x=myrank*calc_length_X-1;
// 	else x=list_X - 1;
// 	for(int y=0;y<list_Y;y++){
// 	for(int z=0;z<list_Z;z++){
// 		_list* li=list+(z + y*list_Z + x*(list_Y*list_Z));
// 		_ptcl* pa;
// 		_ARP2_3::_arp2_3* arp;
// 		for(int j=0;pa=li->ptcl_index[j];j++){
// 			send_list_ptcl_pre.push_back(pa);
// 		}
// 		for(int j=0;arp=li->arp_index[j];j++){
// 			send_list_arp_pre.push_back(arp);
// 		}
// 	}}
// 	if(myrank==numprocs-1)x=0;
// 	else x=myrank*calc_length_X;
// 	for(int y=0;y<list_Y;y++){
// 	for(int z=0;z<list_Z;z++){
// 		_list* li=list+(z + y*list_Z + x*(list_Y*list_Z));
// 		_ptcl* pa;
// 		_ARP2_3::_arp2_3* arp;
// 		for(int j=0;pa=li->ptcl_index[j];j++){
// 			send_list_ptcl_bck.push_back(pa);
// 		}
// 		for(int j=0;arp=li->arp_index[j];j++){
// 			send_list_arp_bck.push_back(arp);
// 		}
// 	}}
// 	x=myrank*calc_length_X;
// 	for(int y=0;y<list_Y;y++){
// 	for(int z=0;z<list_Z;z++){
// 		_list* li=list+(z + y*list_Z + x*(list_Y*list_Z));
// 		_ptcl* pa;
// 		_ARP2_3::_arp2_3* arp;
// 		for(int j=0;pa=li->ptcl_index[j];j++){
// 			send_list_ptcl_pre_inner.push_back(pa);
// 		}
// 		for(int j=0;arp=li->arp_index[j];j++){
// 			send_list_arp_pre_inner.push_back(arp);
// 		}
// 	}}
// 	x=(myrank+1)*calc_length_X-1;
// 	for(int y=0;y<list_Y;y++){
// 	for(int z=0;z<list_Z;z++){
// 		_list* li=list+(z + y*list_Z + x*(list_Y*list_Z));
// 		_ptcl* pa;
// 		_ARP2_3::_arp2_3* arp;
// 		for(int j=0;pa=li->ptcl_index[j];j++){
// 			send_list_ptcl_bck_inner.push_back(pa);
// 		}
// 		for(int j=0;arp=li->arp_index[j];j++){
// 			send_list_arp_bck_inner.push_back(arp);
// 		}
// 	}}
// 	
// 	ptcl_sender.clear();
// 	ptcl_recver.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_pre.begin();p!=send_list_ptcl_pre.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->force;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_pre.begin();p!=send_list_arp_pre.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->force;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	MPI_Status sstatus;
// 	MPI_Status rstatus;
// 	MPI_Request srequest;
// 	MPI_Request rrequest;
// 
// 	
// // 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank,MPI_COMM_WORLD,&request);
// // 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck,MPI_COMM_WORLD,&request);//requestは上と同じだとダメですね
// // 	MPI_Send(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank,MPI_COMM_WORLD);
// // 	MPI_Recv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck,MPI_COMM_WORLD);
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	int recvCount;
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// //	std::cout<<myrank<<'\t'<<ptcl_recver.size()<<'\t'<<recvCount<<std::endl;
// //注　ptcl_recver.size()=0です．vectorの関数を使うとひどいめに会いそう
// 	
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			(ptcl+ptcl_recver[i].index)->force+=ptcl_recver[i].vec;
// 		}else {
// 			(arp2_3+(ptcl_recver[i].index-Number_Ptcl))->force+=ptcl_recver[i].vec;
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// 	
// 	ptcl_sender.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_bck.begin();p!=send_list_ptcl_bck.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->force;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_bck.begin();p!=send_list_arp_bck.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->force;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+30,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+30,MPI_COMM_WORLD,&rrequest);//requestは上と同じだとダメですね
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			(ptcl+ptcl_recver[i].index)->force+=ptcl_recver[i].vec;
// 		}else {
// 			(arp2_3+(ptcl_recver[i].index-Number_Ptcl))->force+=ptcl_recver[i].vec;
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// 	
// }
// 
// void mympi_send_bck(int myrank,int numprocs){
// 	ptcl_sender.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_bck_inner.begin();p!=send_list_ptcl_bck_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_bck_inner.begin();p!=send_list_arp_bck_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_bck_inner_new.begin();p!=send_list_ptcl_bck_inner_new.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl+Number_Arp2_3;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 		tmp.index=(*p)->index+Number_Ptcl+Number_Arp2_3;
// 		tmp.vec=(*p)->loc_pre;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_bck_inner_new.begin();p!=send_list_arp_bck_inner_new.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 		tmp.index=(*p)->index+2*Number_Ptcl+Number_Arp2_3;
// 		tmp.vec=(*p)->loc_pre;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	MPI_Status sstatus;
// 	MPI_Status rstatus;
// 	MPI_Request srequest;
// 	MPI_Request rrequest;
// 	
// // 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD,&request);
// // 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD,&request);
// // 	MPI_Send(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD);
// // 	MPI_Recv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD);
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD,&rrequest);
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	int recvCount;
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// 	
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			_ptcl* p=(ptcl+ptcl_recver[i].index);
// 			p->loc=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else if(ptcl_recver[i].index<Number_Ptcl+Number_Arp2_3){
// 			_arp2_3* p=(arp2_3+(ptcl_recver[i].index-Number_Ptcl));
// 			p->loc=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else if(ptcl_recver[i].index<2*Number_Ptcl+Number_Arp2_3){
// 			_ptcl* p=(ptcl+(ptcl_recver[i].index-Number_Ptcl-Number_Arp2_3));
// 			p->loc=ptcl_recver[i].vec;
// 			i++;
// 			p->loc_pre=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else if(ptcl_recver[i].index<2*Number_Ptcl+2*Number_Arp2_3){
// 			_arp2_3* p=(arp2_3+(ptcl_recver[i].index-2*Number_Ptcl-Number_Arp2_3));
// 			p->loc=ptcl_recver[i].vec;
// 			i++;
// 			p->loc_pre=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// 	
// 	ptcl_sender.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_pre_inner.begin();p!=send_list_ptcl_pre_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_pre_inner_new.begin();p!=send_list_ptcl_pre_inner_new.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_pre_inner.begin();p!=send_list_arp_pre_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_pre_inner_new.begin();p!=send_list_arp_pre_inner_new.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	
// // // 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD,&request);
// // 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD,&request);
// // 	MPI_Send(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD);
// // 	MPI_Recv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD);
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD,&rrequest);
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// 	
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			_ptcl* p=(ptcl+ptcl_recver[i].index);
// 			p->loc=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else {
// 			_arp2_3* p=(arp2_3+(ptcl_recver[i].index-Number_Ptcl));
// 			p->loc=ptcl_recver[i].vec;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// }
// void mympi_send_bck_second(int myrank,int numprocs){
// 	ptcl_sender.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_bck_inner.begin();p!=send_list_ptcl_bck_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_bck_inner.begin();p!=send_list_arp_bck_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	MPI_Status sstatus;
// 	MPI_Status rstatus;
// 	MPI_Request srequest;
// 	MPI_Request rrequest;
// 	
// // 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD,&request);
// // 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD,&request);
// // 	MPI_Send(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD);
// // 	MPI_Recv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD);
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,myrank+10,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,pair_pre+10,MPI_COMM_WORLD,&rrequest);
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	int recvCount;
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// 	
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			_ptcl* p=(ptcl+ptcl_recver[i].index);
// 			p->loc=ptcl_recver[i].vec;
// 			p->loc_pre=p->loc;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else {
// 			_arp2_3* p=(arp2_3+(ptcl_recver[i].index-Number_Ptcl));
// 			p->loc=ptcl_recver[i].vec;
// 			p->loc_pre=p->loc;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// 	
// 	ptcl_sender.clear();
// 	for(std::vector<_ptcl*>::iterator p=send_list_ptcl_pre_inner.begin();p!=send_list_ptcl_pre_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	for(std::vector<_arp2_3*>::iterator p=send_list_arp_pre_inner.begin();p!=send_list_arp_pre_inner.end();p++){
// 		_ptcl_sender tmp;
// 		tmp.index=(*p)->index+Number_Ptcl;
// 		tmp.vec=(*p)->loc;
// 		ptcl_sender.push_back(tmp);
// 	}
// 	
// 	
// // // 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD,&request);
// // 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD,&request);
// // 	MPI_Send(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD);
// // 	MPI_Recv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD);
// 	MPI_Isend(ptcl_sender.begin(),ptcl_sender.size()*sizeof(_ptcl_sender),MPI_BYTE,pair_pre,myrank+20,MPI_COMM_WORLD,&srequest);
// 	MPI_Irecv(ptcl_recver.begin(),ptcl_recver.capacity()*sizeof(_ptcl_sender),MPI_BYTE,pair_bck,pair_bck+20,MPI_COMM_WORLD,&rrequest);
// 	MPI_Wait(&rrequest, &rstatus);
// 	MPI_Wait(&srequest, &sstatus);
// 	MPI_Get_count(&rstatus, MPI_BYTE, &recvCount);
// 	recvCount/=sizeof(_ptcl_sender);
// 	
// 	for(int i=0;i<recvCount;i++){
// 		if(ptcl_recver[i].index<Number_Ptcl){
// 			_ptcl* p=(ptcl+ptcl_recver[i].index);
// 			p->loc=ptcl_recver[i].vec;
// 			p->loc_pre=p->loc;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}else {
// 			_arp2_3* p=(arp2_3+(ptcl_recver[i].index-Number_Ptcl));
// 			p->loc=ptcl_recver[i].vec;
// 			p->loc_pre=p->loc;
// 			p->force.IN(0.0,0.0,0.0);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->cnf->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 	}
// 	memset(ptcl_recver.begin(),0,recvCount*sizeof(_ptcl_sender));
// 	
// 
// }

}//namespace mympi

#endif //USE_MPI
