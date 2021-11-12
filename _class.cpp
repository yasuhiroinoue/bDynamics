/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include<iostream>
#include<fstream>
#include <cmath>
#include "_system.h"
#include "_class.h"
#include "_variable.h"
#include "_periodic.h"
#include"_arp2_3.h"

//const double lattice_cut=_ARP2_3::Cut_Off_Arp2_3_Arp2_3;//=3.20 * sigma_LJ;//r_cut;
//const int cell_number_X=12;
//const int cell_number_Y=12;
//const int cell_number_Z=24;
//const double lattice_cut_X=(double)SYS_X/(double)cell_number_X;
//const double lattice_cut_Y=(double)SYS_Y/(double)cell_number_Y;
//const double lattice_cut_Z=(double)SYS_Z/(double)cell_number_Z;
const double lattice_cut_X=(double)SYS_X/(double)list_X;
const double lattice_cut_Y=(double)SYS_Y/(double)list_Y;
const double lattice_cut_Z=(double)SYS_Z/(double)list_Z;

_ptcl::_ptcl(){}
bool _config::Copy(const _config& p){//基本的に相互参照
//		listは別管理でよろしく
//		li=p.li;
//		list_num=p.list_num;	
	feature=p.feature;
	
	plus=p.plus;
	minus=p.minus;
	
	edge_arp2_3=p.edge_arp2_3;
	end_arp2_3=p.end_arp2_3;
	
	bind=p.bind;
	
#ifdef WALL_BOUNDARY
	if(bind!=-1){
		(WALL::bind_wall+bind)->index=ptcl;
	}
#endif// WALL_BOUNDARY
	if(plus){
		plus->cnf->minus=ptcl;
	}
	if(minus){
		minus->cnf->plus=ptcl;
	}
//cout<<this->ptcl->index<<" "<<((plus)?plus->index:-1)<<" "<<((minus)?minus->index:-1)<<endl;
	if(plus==NULL &&minus!=NULL){
		for(int i(0) ;i<plus_end->num ;++i){
			_ptcl* end(plus_end->end(i));
			if(end==p.ptcl){
				plus_end->end_index[i]=this->ptcl->index;
				break;
			}
		}
	}else if(plus!=NULL &&minus==NULL){
		for(int i(0) ;i<minus_end->num ;++i){
			_ptcl* end(minus_end->end(i));
			if(end==p.ptcl){
				minus_end->end_index[i]=this->ptcl->index;
				break;
			}
		}
	}
//	if(plus!=NULL &&minus==NULL){
//		for(int i(0) ;i<plus_end->num ;++i){
//			_ptcl* end(plus_end->end(i));
//			if(end==p.ptcl){
//				end=ptcl;
//				break;
//			}
//		}
//	}else if(plus==NULL &&minus!=NULL){
//		for(int i(0) ;i<minus_end->num ;++i){
//			_ptcl* end(minus_end->end(i));
//			if(end==p.ptcl){
//				end=ptcl;
//				break;
//			}
//		}
//	}
	if(edge_arp2_3){
		edge_arp2_3->edge_ptcl_index=ptcl->index;
	}
	if(end_arp2_3){
		end_arp2_3->end_ptcl_index=ptcl->index;
	}
	return true;
}





//_list* _ptcl::li()   {return  (list_index==-1)?NULL:list+list_index;}
//_ptcl* _ptcl::plus() {return  (plus_index==-1)?NULL:ptcl+plus_index;}
//_ptcl* _ptcl::minus(){return (minus_index==-1)?NULL:ptcl+minus_index;}
////void _ptcl::plus_in(_ptcl* p)  {(p==NULL)?plus_index=-1 : plus_index=p-ptcl;}
////void _ptcl::minus_in(_ptcl* p) {(p==NULL)?minus_index=-1 : minus_index=p-ptcl;}
//void _ptcl::plus_in(int index) {plus_index=index;}
//void _ptcl::minus_in(int index){minus_index=index;}
//
//_ARP2_3::_arp2_3* _ptcl::edge_arp2_3(){return  (edge_arp2_3_index==-1)?NULL:_ARP2_3::arp2_3 + edge_arp2_3_index;}
//_ARP2_3::_arp2_3* _ptcl::end_arp2_3() {return  (end_arp2_3_index==-1)? NULL:_ARP2_3::arp2_3 + end_arp2_3_index;};
//void _ptcl::edge_arp2_3_in(int index){edge_arp2_3_index=index;}
//void _ptcl::end_arp2_3_in(int index){end_arp2_3_index=index;}
_list* _ARP2_3::_arp2_3::li()   {return  (list_index==-1)?NULL:list+list_index;}

_ptcl* _ARP2_3::_arp2_3::edge_ptcl(){return  (edge_ptcl_index==-1)?NULL:&ptcl[edge_ptcl_index];}
_ptcl* _ARP2_3::_arp2_3::end_ptcl() {return  (end_ptcl_index==-1) ?NULL:&ptcl[end_ptcl_index];}
void _ARP2_3::_arp2_3::edge_ptcl_in(int index){edge_ptcl_index=index;}
void _ARP2_3::_arp2_3::end_ptcl_in( int index){end_ptcl_index=index;}

inline int list_index_calc(_vec<double>& loc){
	int lx = (int)std::floor(loc.x/lattice_cut_X);
	int ly = (int)std::floor(loc.y/lattice_cut_Y);
	int lz = (int)std::floor(loc.z/lattice_cut_Z);
	if(lx==list_X){lx -= 1;}
	if(ly==list_Y){ly -= 1;}
	if(lz==list_Z){lz -= 1;}
	return (lz + ly*list_Z + lx*(list_Y*list_Z));
}
inline int list_index_calc( const _vec<double>& loc){
	int lx = (int)std::floor(loc.x/lattice_cut_X);
	int ly = (int)std::floor(loc.y/lattice_cut_Y);
	int lz = (int)std::floor(loc.z/lattice_cut_Z);
	if(lx==list_X){lx -= 1;}
	if(ly==list_Y){ly -= 1;}
	if(lz==list_Z){lz -= 1;}
	return (lz + ly*list_Z + lx*(list_Y*list_Z));
}
void list_calc(_vec<double>& loc,int& lx,int& ly,int& lz){
	lx = (int)std::floor(loc.x/lattice_cut_X);
	ly = (int)std::floor(loc.y/lattice_cut_Y);
	lz = (int)std::floor(loc.z/lattice_cut_Z);
	if(lx==list_X){lx -= 1;}
	if(ly==list_Y){ly -= 1;}
	if(lz==list_Z){lz -= 1;}
	return;
}
void list_calc( const _vec<double>& loc,int& lx,int& ly,int& lz){
	lx = (int)std::floor(loc.x/lattice_cut_X);
	ly = (int)std::floor(loc.y/lattice_cut_Y);
	lz = (int)std::floor(loc.z/lattice_cut_Z);
	if(lx==list_X){lx -= 1;}
	if(ly==list_Y){ly -= 1;}
	if(lz==list_Z){lz -= 1;}
	return;
}
_list* _ptcl::list_calc(){
	int li=list_index_calc(loc);
	return (list+li);
}
_list* _ARP2_3::_arp2_3::list_calc(){
	int li=list_index_calc(loc);
	return (list+li);
}
_list* list_calc(const _vec<double>& loc){
	int li=list_index_calc(loc);
	return (list+li);
}
void _ptcl::add(_list* list1){
	cnf->li=list1;
	cnf->list_num=list1->num;//list内の配列のindex
}
void _ARP2_3::_arp2_3::add(_list* list1){
	list_index=list1 - list;
	list_num=list1->num_arp;
}




void list_initialize(const int myrank,const int numprocs){
// 	list_X = (int) std::floor((double)SYS_X/(double)lattice_cut);
// 	list_Y = (int) std::floor((double)SYS_Y/(double)lattice_cut);
// 	list_Z = (int) std::floor((double)SYS_Z/(double)lattice_cut);
	if(lattice_cut_X>=_ARP2_3::Cut_Off_Arp2_3_Arp2_3 
	&& lattice_cut_Y>=_ARP2_3::Cut_Off_Arp2_3_Arp2_3 
	&& lattice_cut_Z>=_ARP2_3::Cut_Off_Arp2_3_Arp2_3
	&& !(list_X%NX) && !(list_Y%NY) && !(list_Z%NZ)){
//		list_X = cell_number_X;
//		list_Y = cell_number_Y;
//		list_Z = cell_number_Z;
//		calc_length_X=list_X/numprocs;
		calc_length_X=list_X/NX;
		calc_length_Y=list_Y/NY;
		calc_length_Z=list_Z/NZ;
		
	}else{fout_err<<"lattice_err"<<std::endl;exit(0);}
	list_size = list_X *list_Y *list_Z;
	list = new _list[list_size];
	for(int x=0;x<list_X;x++){
	for(int y=0;y<list_Y;y++){
	for(int z=0;z<list_Z;z++){
		int list_index =(z + y*list_Z + x*(list_Y*list_Z));
		int counter = 0;
		for(int lx1=-1;lx1<=1;lx1++){
		for(int ly1=-1;ly1<=1;ly1++){
		for(int lz1=-1;lz1<=1;lz1++){
		int lz=z+lz1;	int ly=y+ly1;	int lx=x+lx1;	
			if(		lx==-1		){lx = list_X-1;}
			else if(lx==list_X	){lx = 0;}
			if(		ly==-1		){ly = list_Y-1;}
			else if(ly==list_Y	){ly = 0;}
#ifndef WALL_BOUNDARY
			if(		lz==-1		){lz = list_Z-1;}
			else if(lz==list_Z	){lz = 0;}
#endif
#ifdef WALL_BOUNDARY
			if(lz==-1	||	lz==list_Z	){
				list[list_index].list_index[counter]=NULL;
				counter++;
				continue;
			}
#endif
			int list_index_sub =lz + ly*list_Z + lx*(list_Y*list_Z);
			list[list_index].list_index[counter] = &list[list_index_sub];
			counter++;
		}
		}
		}
	}
	}
	}
//	for(int i=0;i<numprocs;++i){
//		for(int x=i*calc_length_X+1;x<(i+1)*calc_length_X-1;x++){
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<list_Z;z++){
//			_list* li=list + (z + y*list_Z + x*(list_Y*list_Z));
//			li->myrank=i;
//		}}}
//		int x=i*calc_length_X;
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<list_Z;z++){
//			_list* li=list + (z + y*list_Z + x*(list_Y*list_Z));
//			li->myrank=i+numprocs;
//		}}
//		x=(i+1)*calc_length_X-1;
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<list_Z;z++){
//			_list* li=list + (z + y*list_Z + x*(list_Y*list_Z));
//			li->myrank=i+2*numprocs;
//		}}
//	}
	for(int i(0);i<list_size;++i)(list+i)->myrank=-1;
	{
		for(int x=0;x<list_X;x++){
		for(int y=0;y<list_Y;y++){
		for(int z=0;z<list_Z;z++){
			int dx,dy,dz;
			if(x%calc_length_X==0)dx=0;
			else if((x+1)%calc_length_X==0)dx=2;
			else dx=1;
			if(y%calc_length_Y==0)dy=0;
			else if((y+1)%calc_length_Y==0)dy=2;
			else dy=1;
			if(z%calc_length_Z==0)dz=0;
			else if((z+1)%calc_length_Z==0)dz=2;
			else dz=1;
			
			(list + (z + y*list_Z + x*(list_Y*list_Z)))->myrank=dz+dy*3+dx*9;

		}}}
	}

//	{
//		for(int x=calcArea.x_pre; x<calcArea.x_bck; x++){
//		for(int y=calcArea.y_pre; y<calcArea.y_bck; y++){
//		for(int z=calcArea.z_pre; z<calcArea.z_bck; z++){
//			int dx,dy,dz;
//			if(x==calcArea.x_pre)dx=0;
//			else if(x==calcArea.x_bck-1)dx=2;
//			else dx=1;
//			if(y==calcArea.y_pre)dy=0;
//			else if(y==calcArea.y_bck-1)dy=2;
//			else dy=1;
//			if(z==calcArea.z_pre)dz=0;
//			else if(z==calcArea.z_bck-1)dz=2;
//			else dz=1;
//			
//			(list + (z + y*list_Z + x*(list_Y*list_Z)))->myrank=dz+dy*3+dx*9;
//
//		}}}
//	}
	return;
}
//_ptcl* _Actin::_Plus_End::end(const int index){return (end_index[0]==-1)? NULL:(ptcl + end_index[index]);}
_ptcl* _Actin::_Plus_End::end(const int index){return (end_index[index]==-1)? NULL:(&ptcl[end_index[index]]);}
_ptcl* _Actin::_Plus_End::end(){return (end_index[0]==-1)? NULL:(&ptcl[end_index[0]]);}

void _Actin::_Plus_End::add(_ptcl* b){
//	end[num]=b;
	end_index[num]=b->index;
	++num;
}
void _Actin::_Plus_End::rm(_ptcl* b){
	--num;
	for(int i=0;i<num;++i){
		if(end_index[i]==b->index){
			end_index[i]=end_index[num];
			break;
		}
	}
	end_index[num]=-1;
}
void _Actin::_Plus_End::polymerization_fila(int i,_ptcl* pa){
	pa->cnf->minus=(end(i));//結合情報の代入
	end(i)->cnf->plus=(pa);	//結合情報の代入
	end_index[i] = pa->index;//配列の更新
//	pa->minus=end[i];
//	end[i]->plus=pa;
//	end[i] = pa;
}
void _Actin::_Plus_End::polymerization_mono(_ptcl* pa,_ptcl* pb){
	pa->cnf->minus=(pb);
	pb->cnf->plus=(pa);
	add(pa);
//	end[num] = pa;
//	++num;
	minus_end->add(pb);//良ろしくないなぁ;;
//	pa->minus=pb;
//	pb->plus=pa;
//	end[num] = pa;
//	++num;
//	minus_end->add(pb);//良ろしくないなぁ;;
}
void _Actin::_Plus_End::depolymerization(_ptcl* b){
	rm(b);
//	--num;
//	for(int i=0;i<num;++i){
//		if(end[i]==b){
//			end[i]=end[num];
//			break;
//		}
//	}
//	end[num]=NULL;
	return;
}


_ptcl* _Actin::_Minus_End::end(const int index){return (end_index[index]==-1)? NULL:(&ptcl[end_index[index]]);}
_ptcl* _Actin::_Minus_End::end(){return (end_index[0]==-1)? NULL:(&ptcl[end_index[0]]);}

void _Actin::_Minus_End::add(_ptcl* b){
// 	end[num]=b;
	end_index[num]=b->index;
	++num;
}
void _Actin::_Minus_End::rm(_ptcl* b){
	--num;
	for(int i=0;i<num;++i){
		if(end_index[i]==b->index){
			end_index[i]=end_index[num];
			break;
		}
	}
	end_index[num]=-1;
}
void _Actin::_Minus_End::depolymerization(_ptcl* b,int i){
	_ptcl* b_sub;
	if(b_sub=b->cnf->plus){
		if(!b_sub->cnf->plus){//2粒子のとき
			plus_end->depolymerization(b_sub);//良ろしくないなぁ;;
			--num;
			if(i!=num){
				end_index[i]=end_index[num];
			}
			end_index[num]=-1;
		}else{
			end_index[i]=b->cnf->plus->index;
		}
		{
			b->cnf->plus=NULL;
			b_sub->cnf->minus=NULL;//_in(-1);
//			b->plus_in(NULL);
//			b_sub->minus_in(NULL);
		}
	}
	return;
}


void _calcArea::add(const int myrank,const int numprocs)
{
	if(list_X%NX!=0 || list_Y%NY!=0 || list_Z%NZ!=0)exit(0);
	nx=((int)floor((double)myrank/(double)(NY*NZ)));
	ny=((int)floor((double)myrank/(double)NZ)%NY);
	nz=(myrank%NZ);
	
	x_pre=(nx  )*list_X/NX;
	x_bck=(nx+1)*list_X/NX;
	y_pre=(ny  )*list_Y/NY;
	y_bck=(ny+1)*list_Y/NY;
	z_pre=(nz  )*list_Z/NZ;
	z_bck=(nz+1)*list_Z/NZ;
	
	pair_pre_x=((int)floor((double)myrank/(double)(NY*NZ))==0   )? (myrank+(NX-1)*NY*NZ):(myrank-NY*NZ);
	pair_bck_x=((int)floor((double)myrank/(double)(NY*NZ))==NX-1)? (myrank-(NX-1)*NY*NZ):(myrank+NY*NZ);
	pair_pre_y=((int)floor((double)myrank/(double)NZ)%NY==0   )? (myrank+(NY-1)*NZ):(myrank-NZ); 
	pair_bck_y=((int)floor((double)myrank/(double)NZ)%NY==NY-1)? (myrank-(NY-1)*NZ):(myrank+NZ); 
	pair_pre_z=(myrank%NZ==0   )? (myrank+NZ-1):(myrank-1); 
	pair_bck_z=(myrank%NZ==NZ-1)? (myrank-NZ+1):(myrank+1);
//
//std::cout<<myrank<<"\t"<<pair_pre_x<<"\t"<<pair_bck_x<<"\t"<<pair_pre_y<<"\t"<<pair_bck_y<<"\t"<<pair_pre_z<<"\t"<<pair_bck_z<<std::endl;
//exit(0);
}

