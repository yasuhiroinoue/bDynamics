/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include <iostream>
#include "_class.h"
#include "_variable.h"
#include "_ctrl.h"
#include "_arp2_3.h"
#include "_wall_boundary.h"
using namespace std;

#ifdef WALL_BOUNDARY
_output_variable_wall output_variable_wall;
#endif// WALL_BOUNDARY
_output_PlusEndPosition output_PlusEndPosition;
_PTCL_Ene_Output PTCL_Ene_Output;
_PolymerizationInfo PolymerizationInfo;
_flow flow;
vector<_force_count> Force_Count(32);


void _output_PlusEndPosition::update(){
	_ptcl* pa;
	++counter_sub;
#ifdef WALL_BOUNDARY
	wall_loc+=WALL::WALL_Z[1].loc.z;
	wall_loc2+=WALL::WALL_Z[1].loc.z*WALL::WALL_Z[1].loc.z;
#endif// WALL_BOUNDARY
	for(int i=0;pa=plus_end->end(i);++i){
#ifdef USE_MPI
	if(pa->loc.x>=(double)calcArea.nx*SYS_X/NX && pa->loc.x<(calcArea.nx+1.0)*SYS_X/NX
	&& pa->loc.y>=(double)calcArea.ny*SYS_Y/NY && pa->loc.y<(calcArea.ny+1.0)*SYS_Y/NY
	&& pa->loc.z>=(double)calcArea.nz*SYS_Z/NZ && pa->loc.z<(calcArea.nz+1.0)*SYS_Z/NZ){
#endif
		IN(pa->loc.z);
#ifdef USE_MPI
	}
#endif
	}
}
void _output_PlusEndPosition::IN(const double& location){
	++counter;
#ifdef WALL_BOUNDARY
	loc+=WALL::WALL_Z[1].loc.z-location;
	loc2+=(WALL::WALL_Z[1].loc.z-location)*(WALL::WALL_Z[1].loc.z-location);
#endif//def WALL_BOUNDARY
	loc_sub+=location;
	loc_sub2+=location*location;
	
}
void _PTCL_Ene_Output::update(){
	++counter;
	ptcl_ene_pot +=ene_pot;
	ptcl_ene_pot2+=ene_pot*ene_pot;
	arp_ene_pot  +=_ARP2_3::ene_pot;
	arp_ene_pot2 +=_ARP2_3::ene_pot*_ARP2_3::ene_pot;
#ifdef WALL_BOUNDARY
	wall_ene_pot +=WALL::ene_pot;
	wall_ene_pot2+=WALL::ene_pot*WALL::ene_pot;
#endif//def WALL_BOUNDARY

}


void _PolymerizationInfo::add(const int& fila_index,const int& mono_index){
#ifdef WALL_BOUNDARY
	_PolymerizationInfoSub tmp(step ,t ,fila_index ,mono_index ,filament_force[fila_index] ,filament_force[mono_index] ,
			WALL::WALL_Z[1].loc.z ,WALL::WALL_Z[1].loc.z - ptcl[fila_index].loc.z ,WALL::WALL_Z[1].loc.z - ptcl[fila_index].loc.z,
			ptcl[fila_index].loc ,ptcl[fila_index].loc);
#else
	_PolymerizationInfoSub tmp(step ,t ,fila_index ,mono_index ,filament_force[fila_index] ,filament_force[mono_index] ,
			SYS_Z ,SYS_Z - ptcl[fila_index].loc.z ,SYS_Z - ptcl[fila_index].loc.z,
			ptcl[fila_index].loc ,ptcl[fila_index].loc);

#endif//def WALL_BOUNDARY
	buf.push_back(tmp);
}

void Force_Count_Update(){
	for(int i(0);i<plus_end->num;++i){
		_ptcl* end=plus_end->end(i);
		Force_Count[i].add(filament_force[end->index]);
	}
	return;
}
	void _force_count::add(const double& frc){
		if(myAbs(frc)>1.0e-30){
			if(checker){
				duration.back()++;
			}else{
				duration.push_back(1);
				time.push_back(t);
				checker=1;
			}
			force.push_back(frc);
		}else{
			if(checker){
				checker=0;
			}
		}
	}

	void _force_count::Output(const std::string& ss1, const std::string& ss2, const bool flag, const bool flag_force){
		if(flag){
			std::ofstream fout1(ss1.c_str(),std::ios::out);
			std::ofstream fout2(ss2.c_str(),std::ios::out);
			fout1.close();
			fout2.close();
		}
		if(duration.size()!=0 &&( (checker!=0) || flag_force)){
			int i(0);
			std::ofstream fout1(ss1.c_str(),std::ios::out | std::ios::app);
			int index(0);
			for(std::vector<int>::iterator p=duration.begin(); p!=duration.end(); p++,i++){
				for(int j(0);j<(*p);j++){
					fout1.precision(10);
					fout1<<time[i]+j*PTCL_dt<<'\t'<<j<<'\t'<<force[index]<<std::endl;
					index++;
				}
			}
			fout1.close();
			
			i=0;
			std::ofstream fout2(ss2.c_str(),std::ios::out | std::ios::app);
			for(std::vector<int>::iterator p=duration.begin(); p!=duration.end(); p++,i++){fout2<<time[i]<<'\t'<<(*p)*PTCL_dt<<std::endl;}
			fout2.close();
			
			
			force.clear();
			duration.clear();
			time.clear();
		}return;
	}





_flow_cell_element::_flow_cell_element():count(0),vel(0.0,0.0,0.0){}
void _flow_cell_element::clear(){
	count=0;
//	start_step=t;
	vel.IN(0.0,0.0,0.0);
}
void _flow::clear(){
	start_step=t;
	for(int i(0);i<size.x*size.y*size.z;++i){
		abs[i].clear();
		rel[i].clear();
		abs_mono[i].clear();
		rel_mono[i].clear();
	}
}
_flow::_flow(){
	start_step=t;
	size_length.x=1.0;
	size_length.y=1.0;
	size_length.z=1.0;
	if( SYS_X - ((int)(SYS_X/size_length.x))*size_length.x!=0 
	  ||SYS_Y - ((int)(SYS_Y/size_length.y))*size_length.y!=0 
	  ||SYS_Z - ((int)(SYS_Z/size_length.z))*size_length.y!=0 ){
		cout<<"class _flow init err"<<endl;
	}
	size.x=SYS_X/size_length.x;
	size.y=SYS_Y/size_length.y;
	size.z=SYS_Z/size_length.z;
	
	abs = new _flow_cell_element[size.x*size.y*size.z];
	rel = new _flow_cell_element[size.x*size.y*size.z];
	abs_mono = new _flow_cell_element[size.x*size.y*size.z];
	rel_mono = new _flow_cell_element[size.x*size.y*size.z];
	clear();
}

void _flow_cell_element::IN(const _vec<double>& v){
	count++;
	vel+=v;
}
void _flow::IN(const _vec<double>& l,const _vec<double>& v,const int& feature){
	int index = calc_index(l);
	abs[index].IN(v);
	if(!feature){abs_mono[index].IN(v);}

	_vec<double> locsub(l);
	_vec<double> velsub(v);
#ifdef WALL_BOUNDARY
	locsub.z += SYS_Z - wall_loc;
	velsub.z -= wall_vel;
#endif// WALL_BOUNDARY
	index = calc_index(locsub);
	rel[index].IN(velsub);
	if(!feature){rel_mono[index].IN(velsub);}
}
void _flow::Wall_IN(const double& l,const double& v){
	wall_loc=l;
	wall_vel=v;
}







void OutputClassUpdate(const int myrank, const int numprocs){
	output_PlusEndPosition.update();
	if(!myrank)PTCL_Ene_Output.update();
//	if(!myrank)Force_Count_Update();
}
