/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef _OUTPUT_CLASS_H
#define _OUTPUT_CLASS_H

#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>

#ifdef WALL_BOUNDARY
class _output_variable_wall{
public:
	int counter;
	double loc,loc2;
	double vel,vel2;
	double frc,frc2;
	double frc_ran,frc_ran2;
	double frc_pot,frc_pot2;
	_vec<double> frc_bane,frc_bane2;
	double frc_mono,frc_mono2;
	double frc_fila,frc_fila2;
	
	_output_variable_wall():loc(0.0) ,loc2(0.0) ,vel(0.0) ,vel2(0.0) ,
			frc(0.0) ,frc2(0.0) ,frc_ran(0.0) ,frc_ran2(0.0) ,frc_pot(0.0) ,frc_pot2(0.0) ,
			frc_bane(0.0,0.0,0.0),frc_bane2(0.0,0.0,0.0),
			frc_mono(0.0) ,frc_mono2(0.0) ,frc_fila(0.0) ,frc_fila2(0.0)
			{}
	void IN(const double& l ,const double& v ,const double& f ,const double& fr ,const double& fp ,const _vec<double>& fb ,const double& fm ,const double& ff){
		loc_in(l);
		vel_in(v);
		frc_in(f);
		frc_ran_in(fr);
		frc_pot_in(fp);
		frc_bane_in(fb);
		frc_mono_in(fm);
		frc_fila_in(ff);
		++counter;
		return;
	}
	void loc_in(const double& in){
		loc+=in;
		loc2+=in*in;
	}
	void vel_in(const double& in){
		vel+=in;
		vel2+=in*in;
	}
	void frc_in(const double& in){
		frc+=in;
		frc2+=in*in;
	}
	void frc_ran_in(const double& in){
		frc_ran+=in;
		frc_ran2+=in*in;
	}
	void frc_pot_in(const double& in){
		frc_pot+=in;
		frc_pot2+=in*in;
	}
	void frc_bane_in(const _vec<double>& in){
		frc_bane+=in;
		//frc_bane2+=_vec<double>(in.x*in.x ,in.y*in.y ,in.z*in.z);
		frc_bane2.x+=in.x*in.x;
		frc_bane2.y+=in.y*in.y;
		frc_bane2.z+=in.z*in.z;
	}
	void frc_mono_in(const double& in){
		frc_mono+=in;
		frc_mono2+=in*in;
	}
	void frc_fila_in(const double& in){
		frc_fila+=in;
		frc_fila2+=in*in;
	}
	void clear(){
		counter=0;
		loc=loc2=vel=vel2=0.0;
		frc=frc2=frc_ran=frc_ran2=0.0;
		frc_pot=frc_pot2=0.0;
		frc_bane.IN(0.0 ,0.0 ,0.0);
		frc_bane2.IN(0.0 ,0.0 ,0.0);
		frc_mono=frc_mono2=frc_fila=frc_fila2=0.0;
	}
	void clean(){clear();}
}extern output_variable_wall;
inline std::ostream& operator<<(std::ostream& stream, _output_variable_wall& tmp)
{//friendéŒ¾‚·‚é‚Ì‚ÍŒ™‚¢‚Å‚·D‚»‚à‚»‚àprivate‚ÅéŒ¾‚·‚é‚Ì‚ªk(ry
	return (tmp.counter)? stream<<tmp.loc/(double)tmp.counter<<	'\t'<<(tmp.loc2/(double)tmp.counter - tmp.loc*tmp.loc/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.vel/(double)tmp.counter<<		'\t'<<(tmp.vel2/(double)tmp.counter - tmp.vel*tmp.vel/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc/(double)tmp.counter<<		'\t'<<(tmp.frc2/(double)tmp.counter - tmp.frc*tmp.frc/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc_ran/(double)tmp.counter<<	'\t'<<(tmp.frc_ran2/(double)tmp.counter  - tmp.frc_ran*tmp.frc_ran/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc_pot/(double)tmp.counter<<	'\t'<<(tmp.frc_pot2/(double)tmp.counter  - tmp.frc_pot*tmp.frc_pot/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc_bane/(double)tmp.counter<<'\t'
		<<(tmp.frc_bane2/(double)tmp.counter - _vec<double>(tmp.frc_bane.x*tmp.frc_bane.x ,tmp.frc_bane.y*tmp.frc_bane.y ,tmp.frc_bane.z*tmp.frc_bane.z)/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc_mono/(double)tmp.counter<<'\t'<<(tmp.frc_mono2/(double)tmp.counter - tmp.frc_mono*tmp.frc_mono/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.frc_fila/(double)tmp.counter<<'\t'<<(tmp.frc_fila2/(double)tmp.counter - tmp.frc_fila*tmp.frc_fila/(double)(tmp.counter*tmp.counter))<<'\t'
		<<tmp.counter
		:stream<<"NO DATE"<<endl;
}
#endif// WALL_BOUNDARY

class _output_PlusEndPosition{
public:
	int counter;
	int counter_sub;
	double loc,loc2;
	double wall_loc,wall_loc2,loc_sub,loc_sub2;
	_output_PlusEndPosition():counter(0),wall_loc(0),loc_sub(0),loc_sub2(0),loc(0),loc2(0){}
	void update();
	void IN(const double&);
	friend std::ostream& operator<<(std::ostream& stream, _output_PlusEndPosition& tmp){
		return (tmp.counter)? stream<<tmp.loc/(double)tmp.counter<<	'\t'<<(tmp.loc2/(double)tmp.counter - tmp.loc*tmp.loc/(double)(tmp.counter*tmp.counter))<<'\t'
				<<tmp.wall_loc/(double)tmp.counter_sub<<	'\t'<<(tmp.wall_loc2/(double)tmp.counter_sub - tmp.wall_loc*tmp.wall_loc/(double)(tmp.counter_sub*tmp.counter_sub))<<'\t'
				<<tmp.loc_sub/(double)tmp.counter<<	'\t'<<(tmp.loc_sub2/(double)tmp.counter - tmp.loc_sub*tmp.loc_sub/(double)(tmp.counter*tmp.counter))<<'\t'
				<<tmp.counter<<'\t'<<tmp.counter_sub
			:stream<<"NO DATE"<<endl;
	}
	void clear(){
		counter=counter_sub=0;
		loc=loc2=wall_loc=wall_loc2=loc_sub=loc_sub2=0.0;
//		memset(this,0,sizeof(_output_PlusEndPosition));
	}
	void clean(){clear();}
	void app(const _output_PlusEndPosition& tmp){
		counter    +=tmp.counter;
		counter_sub+=tmp.counter_sub;
		loc        +=tmp.loc;
		loc2       +=tmp.loc2;
		wall_loc   +=tmp.wall_loc;
		wall_loc2  +=tmp.wall_loc2;
		loc_sub    +=tmp.loc_sub;
		loc_sub2   +=tmp.loc_sub2;
	}
}extern output_PlusEndPosition;



class _PTCL_Ene_Output{
public:
	int counter;
	clock_t timer;
	double speed;
	double ptcl_ene_pot,ptcl_ene_pot2;
//	double ene_kine;
//	_vec<double> p;
	double arp_ene_pot,arp_ene_pot2;
	double wall_ene_pot,wall_ene_pot2;
	
	_PTCL_Ene_Output():counter(0),ptcl_ene_pot(0.0),ptcl_ene_pot2(0.0),
			arp_ene_pot(0.0),arp_ene_pot2(0.0),
			wall_ene_pot(0.0),wall_ene_pot2(0.0),
			timer(clock())
			{}
	void update();
	friend std::ostream& operator<< (std::ostream& stream,_PTCL_Ene_Output& tmp){
		return (tmp.counter)? stream<<tmp.ptcl_ene_pot/(double)tmp.counter<<'\t'
				<<(tmp.ptcl_ene_pot2/(double)tmp.counter - tmp.ptcl_ene_pot*tmp.ptcl_ene_pot/(double)(tmp.counter*tmp.counter))<<'\t'
				<<tmp.arp_ene_pot/(double)tmp.counter<< '\t'<<(tmp.arp_ene_pot2/(double)tmp.counter  - tmp.arp_ene_pot * tmp.arp_ene_pot/(double)(tmp.counter*tmp.counter))<<'\t'
#ifdef WALL_BOUNDARY
				<<tmp.wall_ene_pot/(double)tmp.counter<<'\t'<<(tmp.wall_ene_pot2/(double)tmp.counter - tmp.wall_ene_pot*tmp.wall_ene_pot/(double)(tmp.counter*tmp.counter))<<'\t'
#endif
				<<tmp.speed<<'\t'<<tmp.counter
			:stream<<"NO DATE"<<endl;
	}
	void calc_speed(){
		clock_t timer_sub(timer);
		timer=clock();
		speed=(timer - timer_sub)/(double)CLOCKS_PER_SEC/(double)counter*1000.0;
		return;
	}
	void clear(){
		counter=0;
		ptcl_ene_pot=ptcl_ene_pot2=arp_ene_pot=arp_ene_pot2=wall_ene_pot=wall_ene_pot2=0.0;
	}
	void clean(){clear();}
	void app(const _PTCL_Ene_Output& tmp,const int myrank){
		if(!myrank)counter+=tmp.counter;
		ptcl_ene_pot +=tmp.ptcl_ene_pot;
		ptcl_ene_pot2+=tmp.ptcl_ene_pot2;
		arp_ene_pot  +=tmp.arp_ene_pot;
		arp_ene_pot2 +=tmp.arp_ene_pot2;
		wall_ene_pot +=tmp.wall_ene_pot;
		wall_ene_pot2+=tmp.wall_ene_pot2;
	}
}extern PTCL_Ene_Output;

//	fout_ploy_for<<step<<'\t'<<p_f->index<<'\t'<<p_g->index<<'\t'<<filament_force[p_f->index]<<'\t'<<filament_force[p_g->index]<<'\t'
//		<<WALL::WALL_Z[1].loc.z<<'\t'<<WALL::WALL_Z[1].loc.z - p_f->loc.z<<'\t'<<WALL::WALL_Z[1].loc.z - p_g->loc.z<<'\t'<<endl;
class _PolymerizationInfoSub{
public:
	int step;
	double t;
	int fila_index,mono_index;
	_vec<double> locf;
	_vec<double> locg;
	double fila_frc_form_wall,mono_frc_form_wall;
	double wall_loc;
	double fila_loc,mono_loc;
	_PolymerizationInfoSub(const int& a,const double& b,const int& c,const int& d,const double& e,const double& f,const double& g,const double& h,const double& i
							,const _vec<double>& j,const _vec<double>& k):
		step(a) ,t(b) ,fila_index(c) ,mono_index(d) ,fila_frc_form_wall(e) ,mono_frc_form_wall(f) ,wall_loc(g) ,fila_loc(h) ,mono_loc(i) ,locf(j) ,locg(k)
		{}
	_PolymerizationInfoSub(){}
	const std::string header(){
		std::string tmp="step\tt\tindex\t\tfrc_from_wall\t\tposition\t\t\tfila\t\t\tmono\n";
					tmp+="\t\tfila\tmono\tfila\tmono\twall\tfila(from wall)\tmono(from wall)\tx\ty\tz\tx\ty\tz";
		return tmp;
	}
	friend std::ostream& operator<< (std::ostream& stream,_PolymerizationInfoSub& tmp){
		return stream<<tmp.step<<'\t'<<tmp.t<<'\t'<<tmp.fila_index<<'\t'<<tmp.mono_index<<'\t'
			<<tmp.fila_frc_form_wall<<'\t'<<tmp.mono_frc_form_wall<<'\t'
			<<tmp.wall_loc<<'\t'<<tmp.fila_loc<<'\t'<<tmp.locf<<'\t'<<tmp.locg<<'\t'<<tmp.mono_loc;
	}
};
inline bool _PolymerizationInfoSub_name_desc( const _PolymerizationInfoSub& left, const _PolymerizationInfoSub& right){
	return left.t < right.t;
}


class _PolymerizationInfo{
public:
	vector<_PolymerizationInfoSub> buf;
	
	_PolymerizationInfo(){}
	void add(const int& fila_index,const int& mono_index);
	void add(const _PolymerizationInfoSub& tmp){buf.push_back(tmp);}
	void push_back(const int& fila_index,const int& mono_index){add(fila_index ,mono_index);}
	void push_back(const _PolymerizationInfoSub& tmp){add(tmp);}
	void clear(){buf.clear();}
	void clean(){clear();}
	void sort(){std::sort(buf.begin(), buf.end(), _PolymerizationInfoSub_name_desc);}
	const std::string header(){return _PolymerizationInfoSub().header();}
	int size() {return(int)buf.size();}
}extern PolymerizationInfo;

#include <fstream>
#include <string>

template <class TEMPLATE> inline TEMPLATE myAbs(const TEMPLATE& a){return (a<0)?(-1.0*a):a;}

class _force_count{
public:
	vector<double> force;
	vector<int> duration;
	vector<double> time;
	int checker;
	_force_count():checker(0){}
//	void update();
	void add(const double& frc);//{
//		if(myAbs(frc)>1.0e-10){
//			if(checker){
//				duration.back()++;
//			}else{
//				duration.push_back(1);
//				time.pust_back(t);
//				checker=1;
//			}
//			force.push_back(frc);
//		}else{
//			if(checker){
//				checker=0;
//			}
//		}
//	}
	
	void Output(const std::string& ss1, const std::string& ss2, const bool flag, const bool flag_force);//{
//		if(flag){
//			std::ofstream fout1(ss1.c_str(),std::ios::out);
//			std::ofstream fout2(ss2.c_str(),std::ios::out);
//			fout1.close();
//			fout2.close();
//		}
//		if(!checker){
//			int i(0);
//			std::ofstream fout1(ss1.c_str(),std::ios::out | std::ios::app);
//			int index(0);
//			for(std::vector<int>::iterator p=duration.begin(); p!=duration.end(); p++,i++){
//				for(int j(0);j<(*p);j++){
//					fout1<<time[i]+j*PTCL_dt<<"\t"<<force[index]<<std::endl;
//					index++;
//				}
//			}
//			fout1.close();
//			
//			i=0;
//			std::ofstream fout2(ss2.c_str(),std::ios::out | std::ios::app);
//			for(std::vector<int>::iterator p=duration.begin(); p!=duration.end(); p++,i++){fout2<<time[i]<<'\t'<<(*p)<<std::endl;}
//			fout2.close();
//			
//			
//			force.clear();
//			duration.clear();
//			time.clear();
//		}return;
//	}
};
extern vector<_force_count> Force_Count;


class _flow_cell_element{
public:
	int count;//–§“x‚æ‚¤
//	double start_step;
	_vec<double> vel;
	
	_flow_cell_element();//:count(0),vel(0.0,0.0,0.0){start_step=t;}
	void clear();//{
	
	void IN(const _vec<double>& v);
	void add(const _flow_cell_element& t){
		count+=t.count;
		vel+=t.vel;
	}
};

class _flow{
public:
	_vec<double> size_length;
	_vec<int> size;
	double start_step;
	double wall_loc;
	double wall_vel;
	_flow_cell_element* abs;
	_flow_cell_element* rel;
	_flow_cell_element* abs_mono;
	_flow_cell_element* rel_mono;
	
	_flow();
	int calc_index(const _vec<double>& l){
		_vec<int> index((int)(l.x/size_length.x),(int)(l.y/size_length.y) ,(int)(l.z/size_length.z) );
		if(index.x==size.x)index.x--;
		if(index.y==size.y)index.y--;
		if(index.z==size.z)index.z--;
		return (index.x + index.y*size.x + index.z*size.x*size.y);
	}
	void IN(const _vec<double>& l,const _vec<double>& v,const int&);
	void Wall_IN(const double& l,const double& v);
	void clear();//{for(int i(0);i<size.x*size.y*size.z;++i){abs[i].clear();rel[i].clear();}};
}extern flow;


#endif //_OUTPUT_CLASS_H
