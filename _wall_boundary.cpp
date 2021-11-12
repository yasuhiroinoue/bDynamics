/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include"_ctrl.h"
#ifdef WALL_BOUNDARY
#include<iostream>
#include<fstream>
#include<cmath>

#include"_class.h"
#include"_variable.h"
#include"_periodic.h"
#include"_arp2_3.h"
#include"_mympi.h"
#include "_function.h"
using std::cout;
using std::endl;
using _ARP2_3::_arp2_3;
using _ARP2_3::arp2_3;
using _ARP2_3::EPSILON_LJ_ACTIN_ARP;
using _ARP2_3::SIGMA_LJ_ACTIN_ARP;
using _ARP2_3::CutOffUperActinArp2_3;
using _ARP2_3::CutOffUperActinArp2_3_2;

namespace WALL{

double ene_pot=0.0;
//double WALL_Z=0.0;
struct _WALL_Z WALL_Z[2];
double omega=PI;
const double r_wall_ptcl = r_cut;//pow(2.0,1.0/6.0)*sigma_LJ;//2.50 * sigma_LJ;//sigma_LJ*pow(2.0,1.0/6.0)*2.0;
const double r_wall_ptcl2 = r_wall_ptcl*r_wall_ptcl;
const double r_wall_arp = _ARP2_3::Cut_Off_Actin_Arp2_3;//pow(2.0,1.0/6.0)*SIGMA_LJ_ACTIN_ARP;
const double r_wall_arp2 = r_wall_arp*r_wall_arp;
const double r_bind_wall = 1.2*sigma_LJ;
// const double r_bind_wall = 1.2*sigma_LJ;

#ifdef BROWNIAN


void wall_bind_init(){
//	WALL_Z[0].loc.z=offset_impound*SYS_Z/(double)list_Z;//SYS_Z;
	WALL_Z[0].loc.z=0.0;
	WALL_Z[0].loc_pre=WALL_Z[0].loc;
	WALL_Z[1].loc.z=WALL_Z_INIT+offset_impound*SYS_Z/(double)list_Z;//SYS_Z;
	WALL_Z[1].loc_pre=WALL_Z[1].loc;
	bind_X = (int) std::floor((double)SYS_X / (double)del_bind);
	bind_Y = (int) std::floor((double)SYS_Y / (double)del_bind);
	bind_size = bind_X * bind_Y * 2;
	
	//bind_wall = new int[bind_size];
	bind_wall = new _bind_wall[bind_size];
	memset(bind_wall,0,bind_size*sizeof(_bind_wall));
	for(int z=0 ; z<2 ; z++){
	for(int y=0 ; y<bind_Y ; y++){
	for(int x=0 ; x<bind_X ; x++){
		int bi = x + y*bind_X + z*(bind_X*bind_Y);
		bind_wall[bi].loc.IN(del_bind*(x+0.5),del_bind*(y+0.5),z*SYS_Z+offset_impound*SYS_Z/(double)list_Z);
//		if(!z)bind_wall[bi].loc.z+=1.0;
	}}}
}


inline _vec<double> LJ_ptcl(_vec<double>* dloc,double* potential){
	static const double s06=(double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ
						  * (double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ;
	_vec<double> F=_vec<double>(0.0,0.0,0.0);
	double r2=dloc->sqr();
	if(r2>=r_wall_ptcl2){
		*potential=0.0;
		return F;
	}
	if(r2<cutoffuper2){
		r2=cutoffuper2;
		double r06=s06/(cutoffuper2*cutoffuper2*cutoffuper2);
		double r12=r06*r06;
		*potential=4.0*epsilon_LJ*(r12 - r06) * (cutoffuper - dloc->norm() + (double)sigma_LJ) + epsilon_LJ;
		F=4.0*epsilon_LJ*(12.0*r12 - 6.0*r06)**dloc/r2;
		return F;
	}
	double r06=s06/(r2*r2*r2);
	double r12=r06*r06;
	{
		*potential=4.0*epsilon_LJ*(r12 - r06) + epsilon_LJ;
		F=4.0*epsilon_LJ*(12.0*r12 - 6.0*r06)**dloc/r2;
		return F;
	}
}
inline _vec<double> LJ_arp(_vec<double>* dloc,double* potential){
	static const double s06=SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP
						  * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP;
	_vec<double> F=_vec<double>(0.0,0.0,0.0);
	double r2=dloc->sqr();
	if(r2>=r_wall_arp2){
		*potential=0.0;
		return F;
	}
	if( r2<CutOffUperActinArp2_3_2){
		r2=CutOffUperActinArp2_3_2;
		double r06=s06/(CutOffUperActinArp2_3_2*CutOffUperActinArp2_3_2*CutOffUperActinArp2_3_2);
		double r12=r06*r06;
		*potential=4.0*EPSILON_LJ_ACTIN_ARP*(r12 - r06) * (CutOffUperActinArp2_3 - dloc->norm() + SIGMA_LJ_ACTIN_ARP) + EPSILON_LJ_ACTIN_ARP;
		F=4.0*EPSILON_LJ_ACTIN_ARP*(12.0*r12 - 6.0*r06)**dloc/r2;
		return F;
	}
	double r06=s06/(r2*r2*r2);
	double r12=r06*r06;
	{
		*potential=4.0*EPSILON_LJ_ACTIN_ARP*(r12 - r06) + EPSILON_LJ_ACTIN_ARP;
		F=4.0*EPSILON_LJ_ACTIN_ARP*(12.0*r12 - 6.0*r06)**dloc/r2;
		return F;
	}
}


void Wall_ptcl_frc(int myrank){
	for(int i(0);i<Number_Ptcl;++i)filament_force[i]=0.0;
	_ptcl* p;
	_ARP2_3::_arp2_3* arp;
//	for(int x=myrank*calc_length_X;x<(myrank+1)*calc_length_X;x++){
//	for(int y=0;y<list_Y;y++){
//	for(int z=0;z<list_Z;z++){
// 	for(int z=0;z<list_Z;z+=list_Z-1){
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		_list* li = list+(z + y*list_Z + x*(list_Y*list_Z));
		for(int i=0;p=li->ptcl_index[i];++i){
			if( p->loc.z<r_wall_ptcl+WALL_Z[0].loc.z){
				_vec<double> r = _vec<double>(0.0,0.0,p->loc.z-WALL_Z[0].loc.z);
				if(r.z>0){
					double potential;
					_vec<double> force_sub = LJ_ptcl(&r,&potential);
					p->force += force_sub;
					WALL_Z[0].LJ -= force_sub;
					ene_pot += potential;
					
					filament_force[p->index]=force_sub.z;
					
				}
				else{fout_err<<"wall_error_up"<<endl<<" ptcl_number = "<<p-&ptcl[0]<<" step = "<<step<<endl;}
			}
			if(p->cnf->feature){
			if( p->loc.z<r_wall_ptcl+(offset_impound*SYS_Z/(double)list_Z)){
				_vec<double> r = _vec<double>(0.0,0.0,p->loc.z-(offset_impound*SYS_Z/(double)list_Z));
				if(r.z>0){
					double potential;
					_vec<double> force_sub = LJ_ptcl(&r,&potential);
					p->force += force_sub;
					WALL_Z[0].LJ -= force_sub;
					ene_pot += potential;
					
					filament_force[p->index]=force_sub.z;
					
				}
				else{fout_err<<"wall_error_up"<<endl<<" ptcl_number = "<<p-&ptcl[0]<<" step = "<<step<<endl;}
			}
			}
			if( p->loc.z>WALL_Z[1].loc.z - r_wall_ptcl){
				_vec<double> r = _vec<double>(0.0,0.0,p->loc.z - WALL_Z[1].loc.z);
				if(r.z<0){
					double potential;
					_vec<double> force_sub = LJ_ptcl(&r,&potential);
					p->force += force_sub;
					WALL_Z[1].LJ -= force_sub;
					ene_pot += potential;
					
					
					(p->cnf->feature!=0)?force_wall_fila+=force_sub.z:force_wall_mono+=force_sub.z;
					filament_force[p->index]=force_sub.z;
					
				}
				else{fout_err<<"wall_error_down"<<endl<<"ptcl_number = "<<p-&ptcl[0]<<" step = "<<step<<endl;}
			}
		}
		for(int i=0;arp=li->arp_index[i];++i){
			if( arp->loc.z<r_wall_arp+WALL_Z[0].loc.z){
				_vec<double> r = _vec<double>(0.0,0.0,arp->loc.z-WALL_Z[0].loc.z);
				if(1)		r = _vec<double>(0.0,0.0,arp->loc.z-(offset_impound*SYS_Z/(double)list_Z));
				if(r.z>0){
					double potential;
					_vec<double> force_sub = LJ_arp(&r,&potential);
					arp->force += force_sub;
					WALL_Z[0].LJ -= force_sub;
					ene_pot += potential;
				}
				else{fout_err<<"wall_error_up"<<endl<<" ptcl_number = "<<arp-_ARP2_3::arp2_3<<" step = "<<step<<endl;}
			}
			if( arp->loc.z> WALL_Z[1].loc.z - r_wall_arp){
				_vec<double> r = _vec<double>(0.0,0.0,arp->loc.z - WALL_Z[1].loc.z);
				if(r.z<0){
					double potential;
					_vec<double> force_sub = LJ_arp(&r,&potential);
					arp->force += force_sub;
					WALL_Z[1].LJ -= force_sub;
					ene_pot += potential;
				}
				else{fout_err<<"wall_error_down"<<endl<<"ptcl_number = "<<arp-_ARP2_3::arp2_3<<" step = "<<step<<endl;}
			}
		}
	}}}
}
#else //BROWNIAN
/*
void Wall_ptcl(_ptcl& p, _vec<double>& loc){
	
	double mas = p.mass;
	if( loc.z < _GETA ){
		
		_vec<double> mv = p.mass * p.vel;
		
		double tau = (loc.z - _GETA)/ p.vel.z;
		loc.x -= p.vel.x * tau;
		loc.y -= p.vel.y * tau;

#ifndef COUETTE_BOTH_SIDE
		p.vel.x = (double)gauss((SYS_kT/mas),0.0);
		p.vel.y = (double)gauss((SYS_kT/mas),0.0);
#endif
#ifdef COUETTE_BOTH_SIDE
		p.vel.x = (double)gauss((SYS_kT/mas),(-SVX));
		p.vel.y = (double)gauss((SYS_kT/mas),(-SVY));
#endif
		p.vel.z = (double)gamma_noflow((SYS_kT/mas));
		
		
		tau = RAND*PTCL_dt;
		loc.x += p.vel.x * tau;
		loc.y += p.vel.y * tau;
		loc.z  = p.vel.z * tau + _GETA;
		
		mv -= p.mass * p.vel;
		//mv /= PTCL_dt;
		
		fout_err<<"Wall_ptcl_down_gauss ptcl_number= "<<&p-ptcl<<" step = "<<step<<endl;
		
	}else if( loc.z >= (SYS_Z-_GETA) ){
		
		_vec<double> mv = p.mass * p.vel;
		
		double tau = (loc.z - SYS_Z + _GETA) / p.vel.z;
		loc.x -= p.vel.x * tau;
		loc.y -= p.vel.y * tau;
		
		p.vel.x =  (double)gauss((SYS_kT/mas),SVX);
		p.vel.y =  (double)gauss((SYS_kT/mas),SVY);
		p.vel.z = -(double)gamma_noflow((SYS_kT/mas));
		
		tau = RAND*PTCL_dt;
		loc.x += p.vel.x * tau;
		loc.y += p.vel.y * tau;
		loc.z  = p.vel.z * tau + SYS_Z - _GETA;
		
		mv -= p.mass * p.vel;
		//mv /= PTCL_dt;
		
		fout_err<<"Wall_ptcl_up_gauss   ptcl_number= "<<&p-ptcl<<" step = "<<step<<endl;
		
	}
	
	p.loc = loc;
	(void)Periodic(p.loc);
}
void Wall_ptcl_frc(_ptcl& p){
	if( p.loc.z<r_wall){
		_vec<double> r = _vec<double>(0.0,0.0,p.loc.z);
//		if(p.loc.z<0){ r.IN(0.0,0.0,p.loc.z+0.5);}
		if(r.z>0){
			double norm = r.norm();// + _GETA;
			_vec<double> frc = 4.0*epsilon_LJ*12.0/sigma_LJ*pow(sigma_LJ/norm,13.0)*r/norm;
			p.frc_new += frc;
			
		}
		else{fout_err<<"wall_error_up"<<endl<<" ptcl_number = "<<&p-ptcl<<" step = "<<step<<endl;}
	}
	if( p.loc.z>SYS_Z - r_wall){
		_vec<double> r = _vec<double>(0.0,0.0,p.loc.z - (double)SYS_Z);
//		if(p.loc.z>SYS_Z){ r.IN(0.0,0.0,p.loc.z - (double)SYS_Z-0.5);}
		if(r.z<0){
		double norm = r.norm();// + _GETA;
		_vec<double> frc = 4.0*epsilon_LJ*12.0/sigma_LJ*pow(sigma_LJ/norm,13.0)*r/norm;
		p.frc_new += frc;
		
		}
		else{fout_err<<"wall_error_down"<<endl<<"ptcl_number = "<<&p-ptcl<<" step = "<<step<<endl;}
	}
	
}
*/
#endif //BROWNIAN
//
//const double jijuu=0.0;//1.0;
//const double Wall_Friction=100.0;
void Wall_Location_FirstStep(int myrank,int numprocs){
	if(!myrank){
		for(int i=0;i<2;++i){
			_WALL_Z* p=WALL_Z+i;
#ifdef WALL_CYCLIC
			if(!i)p->loc.z=4.0*(1.0-cos(omega*t));
			else p->loc.z=SYS_Z - 4.0*(1.0-cos(omega*t));
#ifndef SHAKE_BOTH_SIDE
			if(!i)p->loc.z=0.0;
			else p->loc.z=SYS_Z - 8.0*(1.0-cos(omega*t));
#endif

#else //WALL_CYCLIC

#ifdef SHAKE_BOTH_SIDE
// 			if(!i)p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z + jijuu)/Wall_Friction );
// 			else  p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z - jijuu)/Wall_Friction );
// #ifdef SRD
// 			if(!i)p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z + jijuu)/Wall_Friction );
// 			else  p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z - jijuu)/Wall_Friction );
// #endif
#else //SHAKE_BOTH_SIDE
			if(i)  p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z - jijuu)/Wall_Friction );
// if(i)cout<<p->LJ.z + p->adhe.z + p->gaus.z - jijuu<<std::endl;
// #ifdef SRD
// 			if(!i)p->loc.z = 0.0;
// 			else  p->loc.z+= (0.5*PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z - jijuu)/Wall_Friction );
// #endif
#endif //SHAKE_BOTH_SIDE

#endif //WALL_CYCLIC

			p->LJ.IN(0.0,0.0,0.0);
			p->adhe.IN(0.0,0.0,0.0);
			p->gaus.IN(0.0,0.0,0.0);
#ifdef SRD
			p->gaus_srd.IN(0.0,0.0,0.0);
#endif
			if(p->loc.z>SYS_Z)p->loc.z=SYS_Z;
//			if(i)p->loc.z=50;
			if(p->loc.z<0)p->loc.z=0.0;
		}
	}
	force_wall_fila=0;
	force_wall_mono=0;
	return;
}
void Wall_Location_SecondStep(int myrank,int numprocs){
	if(!myrank){
	static const double amplitude=sqrt(2.0*SYS_kT*Wall_Friction/PTCL_dt);
	static const double root3=sqrt(3.0)*2.0;
	_vec<double> RandamF(0.0,0.0,0.0);
		for(int i=0;i<2;++i){
			_WALL_Z* p=WALL_Z+i;
#ifdef WALL_CYCLIC
			if(!i)p->loc.z=4.0*(1.0-cos(omega*t));
			else p->loc.z=SYS_Z - 4.0*(1.0-cos(omega*t));
#ifndef SHAKE_BOTH_SIDE
			if(i)p->loc.z=0.0;
			else p->loc.z=SYS_Z - 8.0*(1.0-cos(omega*t));
#endif

#else //WALL_CYCLIC

			RandamF= (_vec<double>(0.0,0.0,RAND-0.5) * (root3 * amplitude));
		if(i){Wall_RandamF=RandamF.z;}
#ifdef SHAKE_BOTH_SIDE
// 			if(!i)p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z + RandamF.z + jijuu)/Wall_Friction );
// 			else  p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z + RandamF.z - jijuu)/Wall_Friction );
// #ifdef SRD
// 			if(!i)p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z + RandamF.z + jijuu)/Wall_Friction );
// 			else  p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z + RandamF.z - jijuu)/Wall_Friction );
// #endif
#else //SHAKE_BOTH_SIDE
			if(i) p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z  + RandamF.z - jijuu)/Wall_Friction );
//if(i)cout<<p->LJ.z + p->adhe.z + p->gaus.z  + RandamF.z - jijuu<<" "<<p->LJ.z<<" "<<p->adhe.z<<" "<<p->gaus.z<<" "<<RandamF.z<<" "<<jijuu<<endl;
// #ifdef SRD
// 			if(!i)p->loc_pre.z = 0.0;
// 			else  p->loc_pre.z+= (PTCL_dt*( p->LJ.z + p->adhe.z + p->gaus.z +p->gaus_srd.z + RandamF.z - jijuu)/Wall_Friction );
// #endif
#endif //SHAKE_BOTH_SIDE

#endif //WALL_CYCLIC


			if(!i){
				
			}else{
				
				output_variable_wall.IN(p->loc_pre.z ,( p->LJ.z + p->adhe.z + p->gaus.z + RandamF.z - jijuu)/Wall_Friction ,
										p->LJ.z + p->adhe.z + p->gaus.z ,RandamF.z ,
										p->LJ.z ,p->adhe ,
										force_wall_mono ,force_wall_fila );
				flow.Wall_IN(p->loc_pre.z,( p->LJ.z + p->adhe.z + p->gaus.z + RandamF.z - jijuu)/Wall_Friction );
			}
			p->LJ.IN(0.0,0.0,0.0);
			p->adhe.IN(0.0,0.0,0.0);
			p->gaus.IN(0.0,0.0,0.0);
#ifdef SRD
			p->gaus_srd.IN(0.0,0.0,0.0);
#endif
			if(p->loc_pre.z>SYS_Z)p->loc_pre.z=SYS_Z;
//			if(i)p->loc_pre.z=50;
			if(p->loc_pre.z<0)p->loc_pre.z=0.0;
			p->loc.z=p->loc_pre.z;
		}
//		std::cout<<RandamF<<std::endl;
//		if(step>0)exit(0);
	}
	force_wall_fila=0;
	force_wall_mono=0;
	return;
}



void bind_cal(int myrank,int numprocs){
	for(int z=0 ; z<2 ; z++){
	for(int y=0 ; y<bind_Y ; y++){
	for(int x=0 ; x<bind_X ; x++){
		int bi = x + y*bind_X + z*(bind_X*bind_Y);
		if(bind_wall[bi].index==NULL){
			memset(&list_sub , 0 , sizeof(_list_sub));
			_vec<double> r=bind_wall[bi].loc;
//		if(r.x>=(double)myrank*SYS_X/numprocs && r.x<(myrank+1.0)*SYS_X/numprocs){
		if(r.x>=(double)calcArea.nx*SYS_X/NX && r.x<(calcArea.nx+1.0)*SYS_X/NX){
		if(r.y>=(double)calcArea.ny*SYS_Y/NY && r.y<(calcArea.ny+1.0)*SYS_Y/NY){
		if(r.z>=(double)calcArea.nz*SYS_Z/NZ && r.z<=(calcArea.nz+1.0)*SYS_Z/NZ){
			
			_list* li=list_calc(r);
			_list* li_sub=0;
			_ptcl* p =0;
			for(int i=0;i<=27;i++){
				if(li_sub=li->list_index[i]){
					for(int k=0 ; p=li_sub->ptcl_index[k] ; k++){
						if(p->cnf->feature){
// 							if(numbering[p->cnf->feature].num!=1){
// 								if(p->num==0 || p->num==numbering[p->cnf->feature].num-1){
							if((p->cnf->plus && !p->cnf->minus) || (!p->cnf->plus && p->cnf->minus)){
								if(!p->cnf->edge_arp2_3 && p->cnf->end_arp2_3){
									_vec<double> delloc = r - p->loc;
									distance_boundary(delloc);
									if(delloc.norm()<=r_bind_wall){
										list_sub.num +=1;
										list_sub.index[list_sub.num]=p;
									}
								}
							}
						}
					}
				}
			}
			if(list_sub.num==1){
				_ptcl* p = list_sub.index[1];
				p->cnf->bind=bi;
				bind_wall[bi].index=p;
#ifdef USE_MPI
				{
					mympi::_conformation_change cc(wall_bind,p->index,bi);
					mympi::bit_sender.pack(&cc);
				}
#endif//USE_MPI
				fout_err<<"bind "<<p->cnf->feature<<' '<<x*del_bind<<' '<<y*del_bind<<' '<<z*SYS_Z<<endl;
			}
			else if(list_sub.num > 1){
				int num = (int)ceil(list_sub.num*RAND);
				while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
				_ptcl* p = list_sub.index[num];
				p->cnf->bind=bi;
				bind_wall[bi].index=p;
#ifdef USE_MPI
				{
					mympi::_conformation_change cc(wall_bind,p->index,bi);
					mympi::bit_sender.pack(&cc);
				}
#endif//USE_MPI
				fout_err<<"bind "<<p->cnf->feature<<' '<<x*del_bind<<' '<<y*del_bind<<' '<<z*SYS_Z<<endl;
			}
		}}}
		}
	}}}
	return;
}
void de_bind(_ptcl* p , _bind_wall* bind ){
	if(p->cnf->feature==0){
		bind->index = 0;
		p->cnf->bind = -1;
		//bind.feature = 0;
#ifdef USE_MPI
		{
			mympi::_conformation_change cc(wall_debind,p->index,bind - bind_wall);
			mympi::bit_sender.pack(&cc);
		}
#endif//USE_MPI
		fout_err<<"de_bind_mono. ptcl_number = "<<p-&ptcl[0]<<" step = "<<step<<endl;
	}
	else{
		_vec<double> delloc = p->loc - bind->loc;
		distance_boundary(delloc);
		if(delloc.norm()>3.0 * sigma_LJ){
			bind->index = 0;
			p->cnf->bind = -1;
#ifdef USE_MPI
			{
				mympi::_conformation_change cc(wall_debind,p->index,bind - bind_wall);
				mympi::bit_sender.pack(&cc);
			}
#endif//USE_MPI
			fout_err<<"de_bind_fila. ptcl_number = "<<p-&ptcl[0]<<" step = "<<step<<endl;
		}
	}
	return;
}
#ifdef BROWNIAN
void bind_loc(){
	for(int z=0 ; z<2 ; z++){
	for(int y=0 ; y<bind_Y ; y++){
	for(int x=0 ; x<bind_X ; x++){
		int bi = x + y*bind_X + z*(bind_X*bind_Y);
		_bind_wall* b=bind_wall+bi;
		_vec<double> loc =b->loc;
		_vec<double> vel(SVX,SVY,0.0);
#ifdef COUETTE_BOTH_SIDE
		if(z==0)	{loc -= 0.5*vel*PTCL_dt;}
#endif
		if(z!=0)	{loc += 0.5*vel*PTCL_dt;}
		while(loc.x<0)		{loc.x += (double)SYS_X;}
		while(loc.x>=SYS_X)	{loc.x -= (double)SYS_X;}
		while(loc.y<0)		{loc.y += (double)SYS_Y;}
		while(loc.y>=SYS_Y)	{loc.y -= (double)SYS_Y;}
		b->loc = loc;
		if(!z){
//			b->loc.z=WALL_Z[0].loc.z;
		}else{
			b->loc.z=SYS_Z-WALL_Z[1].loc.z;
		}
	}}}
	return;
}
void bind_frc(int myrank,int numprocs){
	bind_loc();
	for(int z=0 ; z<2 ; z++){
		for(int y=0 ; y<bind_Y ; y++){
			for(int x=0 ; x<bind_X ; x++){
				int bi = x + y*bind_X + z*(bind_X*bind_Y);
// 				bind_loc(bind_wall+bi);
//				if(bind_wall[bi].loc.x>=(double)myrank*(double)SYS_X/(double)numprocs && bind_wall[bi].loc.x<(myrank+1.0)*(double)SYS_X/(double)numprocs){
				if(bind_wall[bi].loc.x>=(double)calcArea.nx*SYS_X/NX && bind_wall[bi].loc.x<(calcArea.nx+1.0)*SYS_X/NX){
				if(bind_wall[bi].loc.y>=(double)calcArea.ny*SYS_Y/NY && bind_wall[bi].loc.y<(calcArea.ny+1.0)*SYS_Y/NY){
				if(bind_wall[bi].loc.z>=(double)calcArea.nz*SYS_Z/NZ && bind_wall[bi].loc.z<=(calcArea.nz+1.0)*SYS_Z/NZ){
					_ptcl* p;
					if(p=bind_wall[bi].index){
						de_bind(p , bind_wall+bi );
						if(bind_wall[bi].index){
							const _vec<double> &r=bind_wall[bi].loc;
							
							_ptcl* pc;
							if(pc = p->cnf->plus);
							else if(pc = p->cnf->minus);
							else continue;
							
							_vec<double> dloc1 = p->loc - r;
							_vec<double> dloc2 = pc->loc - p->loc;
							_vec<double> dloc3;
							if(z==0){dloc3.IN(0.0,0.0,1.0);}
							else{dloc3.IN(0.0,0.0,-1.0);}
							
							
							distance_boundary(dloc1);
							distance_boundary(dloc2);
							
							double norm1 =dloc1.norm()+1.0e-9;
							double norm2 =dloc2.norm()+1.0e-9;
							dloc1/=norm1;
							dloc2/=dloc2.norm();
							double cosphi = dloc1*dloc2;
							double cosphi2 = dloc1*dloc3;
		//						cout<<angl/PI*180.0<<'\t'<<angl2/PI*180.0<<endl;
							
							_vec<double> force_sub = -1.0*k_bane*(norm1-0.5)*dloc1;
							ene_pot += k_bane*(norm1-0.5)*(norm1-0.5)*0.5;
							
							p->force  += force_sub;
							if(z==0){
								WALL_Z[0].adhe -= force_sub;
							}
							else{
								WALL_Z[1].adhe -= force_sub;
							}
							
//							if(norm1<1.0e-3){
//								if(z==0){dloc1.IN(0.0,0.0,1.0);}
//								else{dloc1.IN(0.0,0.0,-1.0);}
//								norm1=1.0;
//							}
							double F=	k_angl;//=k_angl*(sinphi*cos(eq_angl)-sin(eq_angl)*cosphi)/sinphi
							_vec<double> force_sub1 = F*(dloc1-cosphi*dloc2)/norm2;
							_vec<double> force_sub2 = F*(dloc2-cosphi*dloc1)/norm1;
							_vec<double> force_sub3 = F*(dloc3-cosphi2*dloc1)/norm1;
							ene_pot += k_angl*(1-cosphi);//=k_angl*(1-(cosphi*cos(eq_angl)-sinphi*sin(eq_angl)));
							ene_pot += k_angl*(1-cosphi2);//=k_angl*(1-(cosphi*cos(eq_angl)-sinphi*sin(eq_angl)));
							
							p->force  += force_sub2 + force_sub3 - force_sub1;
							pc->force += force_sub1;
							
							if(z==0){
								WALL_Z[0].adhe -= (force_sub2 + force_sub3);
							}
							else{
								WALL_Z[1].adhe -= (force_sub2 + force_sub3);
							}
						}
					}
				}}}
			}
		}
	}
	return;
}
#else //BROWNIAN
/*void bind_loc(_bind_wall* b){
	_vec<double> loc =b->loc;
	_vec<double> vel(SVX,SVY,0.0);
#ifdef COUETTE_BOTH_SIDE
	if((int)loc.z==0)	{loc -= vel*PTCL_dt;}
#endif
	if((int)loc.z!=0)	{loc += vel*PTCL_dt;}
	while(loc.x<0)		{loc.x += (double)SYS_X;}
	while(loc.x>=SYS_X)	{loc.x -= (double)SYS_X;}
	while(loc.y<0)		{loc.y += (double)SYS_Y;}
	while(loc.y>=SYS_Y)	{loc.y -= (double)SYS_Y;}
	b->loc = loc;
	return;
}
void bind_frc(){
	for(int z=0 ; z<2 ; z++){
		for(int y=0 ; y<bind_Y ; y++){
			for(int x=0 ; x<bind_X ; x++){
				int bi = x + y*bind_X + z*(bind_X*bind_Y);
				bind_loc(bind_wall+bi);
				_ptcl* p = 0;
				if(p=bind_wall[bi].index){
//					_ptcl* p = bind_wall[bi].index;
					de_bind(p , bind_wall+bi );
					if(bind_wall[bi].index){
						_vec<double> r=bind_wall[bi].loc;
						
						_ptcl* pc;
						if(p->num==0){pc = numbering[p->cnf->feature].ptcl_n[1];}
						else{pc = numbering[p->cnf->feature].ptcl_n[p->num -1];}
						
						_vec<double> dloc1 = r - p->loc;
						_vec<double> dloc2 = pc->loc - p->loc;
						_vec<double> dloc3;
						if(z==0){dloc3.IN(0.0,0.0,1.0);}
						else{dloc3.IN(0.0,0.0,-1.0);}
						
						distance_boundary(dloc1);
						distance_boundary(dloc2);
						
						double angl = acos( (dloc1/dloc1.norm()) * (dloc2/dloc2.norm())) + 1.0E-9;
						double angl2 = acos( (dloc1/dloc1.norm()) * (dloc3/dloc3.norm())) + 1.0E-9;
//						cout<<angl/PI*180.0<<'\t'<<angl2/PI*180.0<<endl;
						
						_vec<double> force = frce_bane(dloc1);
						p->frc_new -= force;
						p->frc_new -= frce_angl( dloc1 , dloc2 , angl )+frce_angl( dloc2 , dloc1 , angl );
						p->frc_new += frce_angl( dloc1 , dloc3 , angl2 );
						pc->frc_new += frce_angl( dloc2 , dloc1 , angl );
						ene_pot += pot_bane(dloc1);//ƒ|ƒeƒ“ƒVƒƒƒ‹
						ene_pot += pot_angl( angl ) + pot_angl( angl2 );
					}
				}
				
			}
		}
	}
	return;
}*/
#endif //BROwNIAN

}//namespace wall

#endif //WALL_BOUNDARY
