/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstring>
#include "vec.h"
#include "mymath.h"
#include "_system.h"
#include "_class.h"
#include "_function.h"
#include "_variable.h"
#include "_ctrl.h"
#include "_periodic.h"
#include "_arp2_3.h"
#include "_wall_boundary.h"
#include "_mympi.h"

using namespace std;


/*				初期化・出力				*/
void list_initialize(int myrank,int numprocs);
void struct_init(const int myrank,const int numprocs){
	
#ifdef SRD
	mpc = new _mpc[Number_mpc];
	cell = new _cell[(cell_size_col)];
	node = new _node[(cell_size_col)];
	
	memset(mpc,0,Number_mpc*sizeof(_mpc));
	memset(cell,0,cell_size_col*sizeof(_cell));
	memset(node,0,cell_size_col*sizeof(_node));
#endif //SRD
//	ptcl = new _ptcl[Number_Ptcl];
	ptcl.mynew(Number_Ptcl);	//この構造体はカオスです
	filament_force = new double[Number_Ptcl];
	filament_force_sub = new double[Number_Ptcl];
	counter_num = new int[number_max];
	
	plus_end = new _Actin::_Plus_End;
	minus_end = new _Actin::_Minus_End;
	
//	memset(ptcl,0,Number_Ptcl*sizeof(_ptcl));
//	memset(list,0,list_size*sizeof(_list));
	
	list_initialize(myrank,numprocs);
#ifdef WALL_BOUNDARY
	WALL::wall_bind_init();
#endif// WALL_BOUNDARY
	_ARP2_3::arp_init();
	
	return;
}
void struct_finalize(){
#ifdef SRD
	delete[] mpc;
	delete[] cell;
	delete[] node;
#endif //SRD
	delete[] list;
#ifdef WALL_BOUNDARY
	delete[] WALL::bind_wall;
#endif// WALL_BOUNDARY
	delete[] counter_num;
// 	delete[] numbering;
//	delete[] bond;
	delete plus_end;
	delete minus_end;
	_ARP2_3::arp_finalize();
	return;
}
void adhetion_init(_ptcl* p,int num){
	if(num>1){
		_ptcl* pa=p-1;//フィラメントが粒子番号連続している設定
		_vec<double> delloc = p->loc - pa->loc;
		distance_boundary(delloc);	delloc/=delloc.norm();
		_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
		quat/=quat.norm();
		p->adhesion_loc = _quaternion::rot( quat , pa->adhesion_loc);
		p->adhesion_loc_inite = p->adhesion_loc;
		p->quat.IN(0.0 ,0.0 ,0.0);
		p->quat_pre.IN(0.0 ,0.0 ,0.0);
	}else if(num==1){
		_ptcl* pa=p-1;
		_vec<double> delloc = p->loc - pa->loc;	distance_boundary(delloc);	delloc/=delloc.norm();
		_vec<double> JIKU(RAND-0.5,RAND-0.5,RAND-0.5);
		while(1){
			while( JIKU.norm() <1.0e-6){
				JIKU.x = RAND-0.5;
				JIKU.y = RAND-0.5;
				JIKU.z = RAND-0.5;
			}
			JIKU/=JIKU.norm();
			JIKU=JIKU%delloc;
			if(JIKU.norm()>1.0e-6){
				JIKU/=JIKU.norm();
				break;
			}
		}
		JIKU*=0.5 * sigma_LJ;
		pa->adhesion_loc.IN(JIKU);
		pa->adhesion_loc_inite = pa->adhesion_loc;
		pa->quat.IN(0.0 ,0.0 ,0.0);
		pa->quat_pre.IN(0.0 ,0.0 ,0.0);
		_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
		quat/=quat.norm();
		p->adhesion_loc = _quaternion::rot( quat , pa->adhesion_loc);
		p->adhesion_loc_inite = p->adhesion_loc;
		p->quat.IN(0.0 ,0.0 ,0.0);
		p->quat_pre.IN(0.0 ,0.0 ,0.0);
	}//num==0はスルー
	return;
}

bool ptcl_ini_g_actin(_ptcl* p, const _vec<double>& loc, int feature,int num){
	bool flag=true;
	_list* li_sub=list_calc(loc);
	_list* li_sub2=li_sub;
	
	_ptcl* pa;
	for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
		if(li_sub2=li_sub->list_index[a]){
#else
		li_sub2=li_sub->list_index[a];
#endif //WALL_BOUNDARY
		for(int j=0;pa=li_sub2->ptcl_index[j];j++){
			_vec<double>delloc = pa->loc - loc;
			distance_boundary(delloc);
			if((p-pa!=1 && p-pa!=-1)&& delloc.norm()<sigma_LJ)flag=false;
		}
#ifdef WALL_BOUNDARY
		}
#endif //WALL_BOUNDARY
	}
	if(flag){
		p->loc.IN(loc);
		p->loc_pre.IN(p->loc);
// 		p->time=0.0;
//		p->p_level = 1;
		p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
		p->adhesion_loc_inite = p->adhesion_loc;
		p->quat.IN(0.0,0.0,0.0);
		p->quat_pre=p->quat;
		p->cnf->feature=feature;
		if(feature!=0){
			adhetion_init(p,num);
		}
		p->add(li_sub);
		li_sub->add(p);
	}
	return flag;
}

#ifdef BROWNIAN
void ptcl_initial_condition(){
//	for( int i = 0; i < Number_Ptcl; i++){
//		ptcl[i].mass = PTCL_MASS;
//	}
	int i=0;
	if(Number_Ptcl){	
		for(i=0;i<Number_Ptcl_Active;i++){
			_vec<double>loc;
			while(!i){
				int t(0);
				for(int j=0;j<WALL::bind_X*WALL::bind_Y;j+=2){//フィラメントの固定のためにかなりややこしくなっている
					if(t==4 && (j%2)){
						--j;
						t=0;
					}else if(t==4 && !(j%2)){
						++j;
						t=0;
					}
//					if(t==2 && (j%8)){
//						j+=6;
//						t=0;
//					}else if(t==2 && !(j%8)){
//						j+=10;
//						t=0;
//					}
//					if(t==2 && (j%8)){
//						j-=2;
//						t=0;
//					}else if(t==2 && !(j%8)){
//						j+=2;
//						t=0;
//					}
					++t;
					for(int k=0;k<64;k++){
						if(!k){
#ifndef WALL_BOUNDARY
							loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
							//loc.IN(18.0,18.0,1.0);
							loc.IN((double)(j%WALL::bind_X)*WALL::del_bind+2.0,(double)((j/WALL::bind_X))*WALL::del_bind+2.0,1.0 + offset_impound*SYS_Z/(double)list_Z);
							WALL::bind_wall[j].index=&ptcl[i];		//固定端
							ptcl[i].cnf->bind=j;//&WALL::bind_wall[j];
#endif
							while( !(ptcl_ini_g_actin((&ptcl[i]),loc,1,k)) ){
#ifndef WALL_BOUNDARY
								loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
								exit(0);
#endif
							}
							loc.IN(0.0,0.0,1.0);
							loc/=loc.norm();
						}else{
							_vec<double>loc1 = ptcl[i-k].loc+ loc*(double)k*0.500;
							ptcl[i].cnf->minus=&ptcl[i-1];
							ptcl[i-1].cnf->plus=&ptcl[i];
#ifdef WALL_BOUNDARY
							if(loc1.z<1.0 || loc1.z>SYS_Z-1.0 || loc1.z>WALL_Z_INIT + offset_impound*SYS_Z/(double)list_Z){
								exit(0);
								i=i-k;
								--j;
								break;
							}
#endif
							Periodic( loc1 );
							if( !(ptcl_ini_g_actin(&ptcl[i],loc1,1,k)) ){
								fout_err<<"ptcl_init_error\ni="<<i<<endl;
								exit(0);
								i=i-k;
								--j;
								break;
							}
						}
						i++;
						if(i==Number_Ptcl_Active)break;
					}if(i==Number_Ptcl_Active)break;
				}if(i==Number_Ptcl_Active)break;
			}if(i==Number_Ptcl_Active)break;
			{//モノマー(粒子浴以外)
#ifndef WALL_BOUNDARY
//				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
// 				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)SYS_Z-2.0)*RAND+1.0);
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)WALL_Z_INIT-1.0)*RAND + offset_impound*SYS_Z/(double)list_Z);
#endif
				while( !(ptcl_ini_g_actin(&(ptcl[i]),loc,0,0)) ){
#ifndef WALL_BOUNDARY
					loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
// 					loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)SYS_Z-2.0)*RAND+1.0);
					loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)WALL_Z_INIT-1.0)*RAND + offset_impound*SYS_Z/(double)list_Z);
#endif
				}
			}
		}
	}
	{//モノマー(粒子浴)	
		for(i=Number_Ptcl_Active;i<Number_Ptcl_Active+Number_Ptcl_Spare_Dens;i++){
			_vec<double>loc;
			{
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
				while( !(ptcl_ini_g_actin(&(ptcl[i]),loc,0,0)) ){
					loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
				}
			}
		}
	}
	for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl_Active];p++){
		if(p->cnf->plus && !p->cnf->minus && (p->cnf->end_arp2_3==NULL))minus_end->add(p);
		else if(!p->cnf->plus && p->cnf->minus)plus_end->add(p);
	}
	for(int k=0;k<Number_Ptcl_Active+Number_Ptcl_Spare_Dens;k++){Periodic( ptcl[k].loc );}
	for(int k=0;k<Number_Ptcl_Active+Number_Ptcl_Spare_Dens;k++){ ptcl[k].loc_pre.IN(ptcl[k].loc);}
	fout_err<<"no_fila = "<<no_fila<<endl;
	fout_err<<"ptcl_init_ok "<<endl;
//	debug();
#ifdef ARP2_3
	_ARP2_3::arp_posion_init();
#endif //ARP2_3
	return;
}
#else //BROWNIAN
/*
void ptcl_initial_condition(){
	for( int i = 0; i < Number_Ptcl; i++){
		ptcl[i].mass = PTCL_MASS;
	}
	int i=0;
	int j=1;
	int im=0;
	if(Number_Ptcl!=0){	
		for(int x=0;x<SYS_X;x++){
			for(int y=0;y<SYS_Y;y++){
				for(int z=0;z<SYS_Z;z++){
					if(i<32){
						if( (x%16== 0) && (y%2==0) && (z%2==0) ){
							for(int k=0 ; k<32 ; k++){
								ptcl[i+k].loc.IN((double)x+0.5+k,(double)y+0.5,(double)z+1.5);
								//ptcl[i+k].loc.IN((double)x+0.5+k*l_eq,(double)y+16.0,(double)z+16.0);
						//		ptcl[i+k].vel.IN(0.0 , 0.0 , 0.0);
						//		ptcl[i+k].vel.IN((double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0));
						//		ptcl[i+k].vel.IN(2.0*(double)SVX/(double)SYS_Z*ptcl[i].loc.z - (double)SVX , 0.0 , 0.0);
								ptcl[i+k].vel.IN((double)gauss(SYS_kT/PTCL_MASS,2.0*(double)SVX/(double)SYS_Z*ptcl[i].loc.z - (double)SVX ) ,
											(double)gauss(SYS_kT/PTCL_MASS,0.0),(double)gauss(SYS_kT/PTCL_MASS,0.0));
// 								ptcl[i+k].time = 0.0;
								ptcl[i+k].p_level = 1;
								
								ptcl[i+k].feature = j;
								ptcl[i+k].num = k;
								numbering[j].num++;
								numbering[j].ptcl_n[k]=&ptcl[i+k];
							}
							i += 32;
							j++;
						}
					}
					else{
						if( (x%2== 0) && (y%2==0) && (z%2==0) && (y+1)*(z+1)!=1 ){
							ptcl[i].loc.IN((double)x+0.5,(double)y+0.5,(double)z+1.0);
							//ptcl[i].vel.IN((double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0));
						//	ptcl[i].vel.IN(0.0 , 0.0 , 0.0);
						//	ptcl[i].vel.IN((double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0),(double)gauss(SYS_kT,0.0));
						//	ptcl[i].vel.IN(2.0*(double)SVX/(double)SYS_Z*ptcl[i].loc.z - (double)SVX , 0.0 , 0.0);
							ptcl[i].vel.IN((double)gauss(SYS_kT/PTCL_MASS,2.0*(double)SVX/(double)SYS_Z*ptcl[i].loc.z - (double)SVX ) ,
								(double)gauss(SYS_kT/PTCL_MASS,0.0),(double)gauss(SYS_kT/PTCL_MASS,0.0));
							
// 							ptcl[i].time = 0.0;
							ptcl[i].p_level = 1;
							
							ptcl[i].feature = 0;
							ptcl[i].num = im;
							numbering[0].num++;
							numbering[0].ptcl_n[im]=&ptcl[i];
							
							i++;
							im++;
						}
					}
					if(i==Number_Ptcl)break;
				}
				if(i==Number_Ptcl)break;
			}
			if(i==Number_Ptcl)break;
		}
		no_fila=j-1;
	}
	else{	no_fila=0;}
	
	if(i != Number_Ptcl){	fout_err<<"ptcl_init_error\n";	exit(0);	}
	fout_err<<"no_fila = "<<no_fila<<endl;
	fout_err<<"ptcl_init_ok "<<endl;
//	for(int i=1 ; i<no_fila+1 ; i++){
//		counter_severing_max+=countfea[i]-1;
//	}
//	cout<<"ok"<<endl;exit(0);
	for(int k=0 ; k<Number_Ptcl ; k++){Periodic( ptcl[k].loc );}
	list_cal();
	numbering_debug();
	return;
}*/
#endif//BROWNIAN


// inline double pot( double r){//ポテンシャル(スカラー)
// 	return (4.0*epsilon_LJ)*pow(sigma_LJ/r,12.0);
// }
// inline double dif_pot( double r){
// 	return -1.0*(4.0*epsilon_LJ)*12.0/sigma_LJ*pow(sigma_LJ/r,13.0);
// }
// inline double d_dif_pot( double r){
// 	return (4.0*epsilon_LJ)*12.0/sigma_LJ*13.0/sigma_LJ*pow(sigma_LJ/r,14.0);
// }
// 
// inline double repot( _vec<double> r){
// 	if(r.norm()<=r_cut){return pot(r.norm())-pot(r_cut)+dif_pot(r_cut)*(r.norm() - r_cut)+1/2.0*d_dif_pot(r_cut)*pow(r.norm() - r_cut,2);}
// 	else{return 0.0;}
// }
// inline double pot_bane( _vec<double> r){
// 	return k_bane*(r.norm()-rb)*(r.norm()-rb)/2.0;
// }
// inline double pot_angl( double a ){
// 	return k_angl*(a-eq_angl)*(a-eq_angl)/2.0;
// }
// inline double dif_pot_angl( double a ){
// 	return k_angl*(a-eq_angl);
// }
// inline _vec<double> refrce( _vec<double> r){//LJポテンシャル
// 	if(r.norm()<=r_cut){return -1.0*(dif_pot(r.norm())+dif_pot( r_cut )+d_dif_pot( r_cut )*(r.norm() - r_cut))*r/r.norm();}
// 	else{return _vec<double>(0.0,0.0,0.0);}
// }
// inline _vec<double> frce_bane( _vec<double> r){//バネ
// 	return -1.0*k_bane*(r.norm()-rb)*r/r.norm();
// }
// inline _vec<double> frce_angl( _vec<double> r1 , _vec<double> r2 , double a ){//角、バネ
// 	return -1.0*dif_pot_angl( a )*(-1.0*r2/r2.norm()+r1/r1.norm()*cos(a))/(r1.norm()*sin(a));
// }

_vec<double> LJ(_ptcl* p1,_ptcl* p2,double* potential){//26
	static const double s06=(double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ
						  * (double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ;
	_vec<double> F(0.0,0.0,0.0);
	_vec<double> dloc=p1->loc - p2->loc;					//---
	distance_boundary(dloc);
	double r2=dloc.sqr();									//***++
	if(r2>=r_cut2){
		*potential=0.0;
		return F;
	}
	if(r2<cutoffuper2){
		r2=cutoffuper2;
		double r06=s06/(cutoffuper2*cutoffuper2*cutoffuper2);
		double r12=r06*r06;
		*potential=4.0*epsilon_LJ*(r12 - r06) * (cutoffuper - dloc.norm() + (double)sigma_LJ) + epsilon_LJ;
		F=4.0*epsilon_LJ*(12.0*r12 - 6.0*r06)*dloc/r2;
		return F;
	}
	double r06=s06/(r2*r2*r2);								//**/
	double r12=r06*r06;										//*
	{
		*potential=4.0*epsilon_LJ*(r12 - r06) + epsilon_LJ;	//-**+
		F=(4.0*epsilon_LJ*(12.0*r12 - 6.0*r06)/r2)*dloc;		//**-**/***
		return F;
	}
}

#ifdef BROWNIAN

void Adhesion_Torsion_Force(_ptcl* p2,_ptcl* p3){
	_vec<double> r23 = p3->loc - p2->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
	_vec<double> n123 = p2->adhesion_loc % r23;
	_vec<double> n234 = r23 % p3->adhesion_loc;
	double norm1 = n123.norm();
	double norm2 = n234.norm();
	n123/=norm1;
	n234/=norm2;
	
	double cosph = n123*n234;
	double sinph = -1.0*((n123%n234)*ur23);
	//cout<<std::asin(sinph)*180/PI<<endl;
//	sinph *= -1.0;
//	cout<<std::acos(cosph)*180/PI<<' '<<std::asin(sinph)*180/PI<<endl;
	if(sinph>=0.0)sinph+=1.0e-9;
	else sinph-=1.0e-9;
	
	static const double cosph0 = n123*n234;
	static const double sinph0 = -1.0*((n123%n234)*ur23);
	//static double cosph0 = std::cos(hannshu_angle);
	//static double sinph0 = std::sin(hannshu_angle);
//	cout<<std::acos(cosph0)*180/PI<<' '<<std::asin(sinph0)*180/PI<<endl;
	
	double F = -1.0 * K_Torsion*(sinph*cosph0 - cosph*sinph0)/sinph;
//	cout<<F<<endl;
	
	_vec<double> N1 = (n234 - cosph*n123)/norm1;
	_vec<double> N2 = (n123 - cosph*n234)/norm2;
	
	_vec<double> force_sub1 = F * (N1 % r23);//K_Torsion*(sinph*cosph0 - cosph*sinph0)*r23.norm()/norm1*n123;//F * (N1 % r23);
	_vec<double> force_sub2 = F * (N1 % p2->adhesion_loc) -  F * (N2 % p3->adhesion_loc);
	_vec<double> force_sub3 = F * (N2 % r23);
//	cout << "f1 " << force_sub1 << endl;
//	cout << "f2 " << force_sub2 << endl;
//	cout << "f3 " << force_sub3 << endl;
//	if(t>0.1)exit(0);
	p2->adhesion_frc += force_sub1;
	p2->force += force_sub2 - force_sub1;
	p3->force += force_sub3 - force_sub2;
	p3->adhesion_frc -= force_sub3;
	
	
	
	return;
}

// void Adhesion_Torsion_Force(_Actin::_Bond* b){
// 	//_vec<double> r23 = p3->loc - p2->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
// 	_vec<double> ur23 = b->udloc;
// 	_vec<double> n123 = (b->minus->adhesion_loc % ur23) * b->norm;
// 	_vec<double> n234 = (ur23 % b->plus->adhesion_loc) * b->norm;
// 	double norm1 = n123.norm();
// 	double norm2 = n234.norm();
// 	n123/=norm1;
// 	n234/=norm2;
// 	
// 	double cosph = n123*n234;
// 	double sinph = -1.0*((n123%n234)*ur23);
// 	//cout<<std::asin(sinph)*180/PI<<endl;
// //	sinph *= -1.0;
// //	cout<<std::acos(cosph)*180/PI<<' '<<std::asin(sinph)*180/PI<<endl;
// 	if(sinph>=0.0)sinph+=1.0e-9;
// 	else sinph-=1.0e-9;
// 	
// 	static const double cosph0 = n123*n234;
// 	static const double sinph0 = -1.0*((n123%n234)*ur23);
// 	//static double cosph0 = std::cos(hannshu_angle);
// 	//static double sinph0 = std::sin(hannshu_angle);
// //	cout<<std::acos(cosph0)*180/PI<<' '<<std::asin(sinph0)*180/PI<<endl;
// 	
// 	double F = -1.0 * K_Torsion*(sinph*cosph0 - cosph*sinph0)/sinph;
// //	cout<<F<<endl;
// 	
// 	_vec<double> N1 = (n234 - cosph*n123)/norm1;
// 	_vec<double> N2 = (n123 - cosph*n234)/norm2;
// 	
// 	_vec<double> force_sub1 = F * (N1 % ur23)* b->norm;//K_Torsion*(sinph*cosph0 - cosph*sinph0)*r23.norm()/norm1*n123;//F * (N1 % r23);
// 	_vec<double> force_sub2 = F * (N1 % b->minus->adhesion_loc) -  F * (N2 % b->plus()->adhesion_loc);
// 	_vec<double> force_sub3 = F * (N2 % ur23)* b->norm;
// //	cout << "f1 " << force_sub1 << endl;
// //	cout << "f2 " << force_sub2 << endl;
// //	cout << "f3 " << force_sub3 << endl;
// //	if(t>0.1)exit(0);
// 	b->minus->adhesion_frc += force_sub1;
// 	b->minus->force += force_sub2 - force_sub1;
// 	b->plus()->force += force_sub3 - force_sub2;
// 	b->plus()->adhesion_frc -= force_sub3;
// 	
// 	return;
// }
void Adhesion_Torsion_Force(_ptcl* pa){//二面角ポテンシャルもっと簡単にできる．というか不安定点があるのが解消できる．
	//_vec<double> r23 = p3->loc - p2->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
	_vec<double> ur23 = pa->udloc;
	_vec<double> n123 = (pa->adhesion_loc % ur23) * pa->norm;
	_vec<double> n234 = (ur23 % pa->cnf->plus->adhesion_loc) * pa->norm;
	double norm1 = n123.norm();
	double norm2 = n234.norm();
	n123/=norm1;
	n234/=norm2;
	
	double cosph = n123*n234;
	double sinph = -1.0*((n123%n234)*ur23);
	//cout<<std::asin(sinph)*180/PI<<endl;
//	sinph *= -1.0;
//	cout<<std::acos(cosph)*180/PI<<' '<<std::asin(sinph)*180/PI<<endl;
	if(sinph>=0.0)sinph+=1.0e-9;
	else sinph-=1.0e-9;
	
	static const double cosph0 = n123*n234;
	static const double sinph0 = -1.0*((n123%n234)*ur23);
	//static double cosph0 = std::cos(hannshu_angle);
	//static double sinph0 = std::sin(hannshu_angle);
//	cout<<std::acos(cosph0)*180/PI<<' '<<std::asin(sinph0)*180/PI<<endl;
	
	double F = -1.0 * K_Torsion*(sinph*cosph0 - cosph*sinph0)/sinph;
//	cout<<F<<endl;
	
	_vec<double> N1 = (n234 - cosph*n123)/norm1;
	_vec<double> N2 = (n123 - cosph*n234)/norm2;
	
	_vec<double> force_sub1 = F * (N1 % ur23)* pa->norm;//K_Torsion*(sinph*cosph0 - cosph*sinph0)*r23.norm()/norm1*n123;//F * (N1 % r23);
	_vec<double> force_sub2 = F * (N1 % pa->adhesion_loc) -  F * (N2 % pa->cnf->plus->adhesion_loc);
	_vec<double> force_sub3 = F * (N2 % ur23)* pa->norm;
//	cout << "f1 " << force_sub1 << endl;
//	cout << "f2 " << force_sub2 << endl;
//	cout << "f3 " << force_sub3 << endl;
//	if(t>0.1)exit(0);
	pa->adhesion_frc += force_sub1;
	pa->force += force_sub2 - force_sub1;
	pa->cnf->plus->force += force_sub3 - force_sub2;
	pa->cnf->plus->adhesion_frc -= force_sub3;
	
//if(step>10){
//	cout<<force_sub1<<'\n'
//		<<force_sub2<<'\n'
//		<<force_sub3<<'\n'<<endl;
//	
//	cout<<K_Torsion*(sinph*cosph0 - cosph*sinph0)*pa->norm/norm1*n123<<'\n'
//		<<K_Torsion*(sinph*cosph0 - cosph*sinph0)*( (pa->adhesion_loc*pa->udloc)/norm1*n123 + (pa->udloc*pa->plus()->adhesion_loc)/norm2*n234)<<'\n'
//		<<-1.0*K_Torsion*(sinph*cosph0 - cosph*sinph0)*pa->norm/norm2*n234<<"\n\n\n"<<endl;
//	
//	exit(0);
//}		
	return;
}


void Adhesion_Force(_ptcl* pa){
	double norm_p,norm_m;
	_vec<double> udloc_p,udloc_m;
	_vec<double> force_sub_p1,force_sub_p2,force_sub_m1,force_sub_m2,force_sub1;
//	for(_Actin::_Bond* b=bond;b->plus;++b){
		norm_p=pa->cnf->plus->adhesion_loc.norm();
		udloc_p=pa->cnf->plus->adhesion_loc/norm_p;
		norm_m=pa->adhesion_loc.norm();
		udloc_m=pa->adhesion_loc/norm_m;
		
		double cosphi_p = pa->udloc * udloc_p;
		double sinphi_p = std::sqrt(1.0 - cosphi_p*cosphi_p);//0~pi
		double cosphi_m = pa->udloc * udloc_m;
		double sinphi_m = std::sqrt(1.0 - cosphi_m*cosphi_m);//0~pi
		
		double F_p = -1.0*k_Bend_angl*(cosphi_p/sinphi_p);
		double F_m = -1.0*k_Bend_angl*(cosphi_m/sinphi_m);
		force_sub_p1 = F_p * (pa->udloc - cosphi_p*udloc_p )/norm_p;
//		force_sub_p2 = F_p * (udloc_p - cosphi_p*pa->udloc )/pa->norm;
//		force_sub_m1 = F_m * (udloc_m - cosphi_m*pa->udloc )/pa->norm;
		force_sub_m2 = F_m * (pa->udloc - cosphi_m*udloc_m )/norm_m;
		force_sub1=F_p * (udloc_p - cosphi_p*pa->udloc )/pa->norm + F_m * (udloc_m - cosphi_m*pa->udloc )/pa->norm;
		
		pa->cnf->plus->adhesion_frc	+= force_sub_p1 ;
		//pa->plus()->force			+= (-force_sub_p1+force_sub_p2+force_sub_m1);
		//pa->force			-= (force_sub_p2+force_sub_m1+force_sub_m2);
		pa->cnf->plus->force			+= (force_sub1-force_sub_p1);
		pa->force			-= (force_sub1+force_sub_m2);
		pa->adhesion_frc	+= force_sub_m2 ;
		
		
	//	Adhesion_Torsion_Force(b->minus , b->plus());
		Adhesion_Torsion_Force(pa);
	
// 	}
	return;
}
// void Adhesion_Force(){
// 	double norm_p,norm_m;
// 	_vec<double> udloc_p,udloc_m;
// 	_vec<double> force_sub_p1,force_sub_p2,force_sub_m1,force_sub_m2,force_sub1;
// 	for(_Actin::_Bond* b=bond;b->plus;++b){
// 		norm_p=b->plus->adhesion_loc.norm();
// 		udloc_p=b->plus->adhesion_loc/norm_p;
// 		norm_m=b->minus->adhesion_loc.norm();
// 		udloc_m=b->minus->adhesion_loc/norm_m;
// 		
// 		double cosphi_p = b->udloc * udloc_p;
// 		double sinphi_p = std::sqrt(1.0 - cosphi_p*cosphi_p);//0~pi
// 		double cosphi_m = b->udloc * udloc_m;
// 		double sinphi_m = std::sqrt(1.0 - cosphi_m*cosphi_m);//0~pi
// 		
// 		double F_p = -1.0*k_Bend_angl*(cosphi_p/sinphi_p);
// 		double F_m = -1.0*k_Bend_angl*(cosphi_m/sinphi_m);
// 		force_sub_p1 = F_p * (b->udloc - cosphi_p*udloc_p )/norm_p;
// //		force_sub_p2 = F_p * (udloc_p - cosphi_p*b->udloc )/b->norm;
// //		force_sub_m1 = F_m * (udloc_m - cosphi_m*b->udloc )/b->norm;
// 		force_sub_m2 = F_m * (b->udloc - cosphi_m*udloc_m )/norm_m;
// 		force_sub1=F_p * (udloc_p - cosphi_p*b->udloc )/b->norm + F_m * (udloc_m - cosphi_m*b->udloc )/b->norm;
// 		
// 		b->plus->adhesion_frc	+= force_sub_p1 ;
// 		//b->plus->force			+= (-force_sub_p1+force_sub_p2+force_sub_m1);
// 		//b->minus->force			-= (force_sub_p2+force_sub_m1+force_sub_m2);
// 		b->plus->force			+= (force_sub1-force_sub_p1);
// 		b->minus->force			-= (force_sub1+force_sub_m2);
// 		b->minus->adhesion_frc	+= force_sub_m2 ;
// 		
// 		
// 	//	Adhesion_Torsion_Force(b->minus , b->plus);
// 		Adhesion_Torsion_Force(b);
// 	
// 	}
// 	return;
// }
// 	

// void Adhesion_Force(){
// 	for(int i = 1 ; numbering[i].num ; i++){
// 		_ptcl *pc_0, *pc_1, *pc_2;
// 		_vec<double> force_sub0,force_sub1,force_sub2,force_sub3;
// 		_vec<double> del_f_0,del_f_1;
// 		if(numbering[i].num!=1){
// 			for(int j=0 ; ; j++){
// 				if(j==0){
// 					pc_1 = numbering[i].ptcl_n[j];
// 					pc_2 = numbering[i].ptcl_n[j+1];
// 					
// 					del_f_1 = pc_2->loc - pc_1->loc;
// 					distance_boundary(del_f_1);
// 					double d1 = acos( (del_f_1/del_f_1.norm()) * (pc_1->adhesion_loc / sigma_LJ * 2.0) ) + 1.0E-9;
// 					force_sub2 = Adhesion_Bending_Force( pc_1->adhesion_loc , del_f_1 , d1 );
// 					force_sub3 = Adhesion_Bending_Force( del_f_1 , pc_1->adhesion_loc , d1 );
// 					
// 					pc_1->adhesion_frc += force_sub2;
// 					pc_1->force -= force_sub2+force_sub3;
// 					pc_2->force += force_sub3;
// 					
// 					del_f_0 = -1.0*del_f_1;
// 					pc_0=pc_1;
// 					pc_1=pc_2;
// 					
// 					continue;
// 				}
// 				if( !numbering[i].ptcl_n[j+1] ){
// 					
// 					double d0 = acos( (del_f_0/del_f_0.norm()) * (pc_1->adhesion_loc / sigma_LJ * 2.0) ) + 1.0E-9;
// 					
// 					force_sub0 = Adhesion_Bending_Force( del_f_0 , pc_1->adhesion_loc , d0 );
// 					force_sub1 = Adhesion_Bending_Force( pc_1->adhesion_loc , del_f_0 , d0 );
// 					
// 					pc_0->force += force_sub0;
// 					pc_1->force -= force_sub0+force_sub1;
// 					pc_1->adhesion_frc += force_sub1;
// 					
// 					break;
// 				}
// 				
// 				pc_2 = numbering[i].ptcl_n[j+1];
// 				del_f_1 = pc_2->loc - pc_1->loc;
// 				distance_boundary(del_f_1);
// 				double d0 = acos( (del_f_0/del_f_0.norm()) * (pc_1->adhesion_loc / sigma_LJ * 2.0) );
// 				double d1 = acos( (del_f_1/del_f_1.norm()) * (pc_1->adhesion_loc / sigma_LJ * 2.0) );
// 				
// 				force_sub0 = Adhesion_Bending_Force( del_f_0 , pc_1->adhesion_loc , d0 );
// 				force_sub1 = Adhesion_Bending_Force( pc_1->adhesion_loc , del_f_0 , d0 ) + Adhesion_Bending_Force( pc_1->adhesion_loc , del_f_1 , d1 );
// 			//	force_sub2 = Adhesion_Bending_Force( pc_1->adhesion_loc , del_f_1 , d1 );
// 				force_sub3 = Adhesion_Bending_Force( del_f_1 , pc_1->adhesion_loc , d1 );
// 				
// 				pc_0->force += force_sub0;
// 			//	pc_1->force -= force_sub0+force_sub1+force_sub2+force_sub3;;
// 				pc_1->force -= force_sub0+force_sub1+force_sub3;;
// 				pc_2->force += force_sub3;
// 			//	pc_1->adhesion_frc += force_sub1 + force_sub2;
// 				pc_1->adhesion_frc += force_sub1;
// 				
// 				del_f_0 = -1.0*del_f_1;
// 				pc_0=pc_1;
// 				pc_1=pc_2;
// 				
// 				
// 			}
// 			for(int j=0 ; numbering[i].ptcl_n[j+1] ; j++){
// 				Adhesion_Torsion_Force(numbering[i].ptcl_n[j] , numbering[i].ptcl_n[j+1]);
// 			}
// 		}
// 	}
// 	return;
// }
	


void calc_norm_udloc(_list* const li){
	_ptcl* pa(NULL);
	for(int j(0);pa=li->ptcl_index[j];j++){
		if(pa->cnf->plus){
			_vec<double> dloc=pa->cnf->plus->loc - pa->loc;		//---
			distance_boundary(dloc);						
			pa->norm = dloc.norm();							//***++sqrt
			pa->udloc = dloc/pa->norm;						// ///
		}
	}
}
//void calc_norm_udloc(_ptcl* pa){
//	_vec<double> dloc=pa->plus()->loc - pa->loc;		//---
//	distance_boundary(dloc);						
//	pa->norm = dloc.norm();							//***++sqrt
//	pa->udloc = dloc/pa->norm;						// ///
//}
//void linear_spring(_ptcl* pa){
//	_vec<double> force_sub = -1.0*k_bane*(pa->norm-rb)*pa->udloc;	//-*****
//	ene_pot += k_bane*(pa->norm-rb)*(pa->norm-rb)*0.5;				//--***
//	pa->plus()->force += force_sub;									//+++
//	pa->force -= force_sub;											//---
//
//}
//void bend_spring(_ptcl* pa, _ptcl* pb){
//	double cosphi = (pa->udloc) * (pb->udloc);					//***++
////			double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
////			double d=acos(cosphi);
//	//double F=	k_angl*(d -eq_angl)/sinphi;
//	double F=	k_angl;//=k_angl*(sinphi*cos(eq_angl)-sin(eq_angl)*cosphi)/sinphi
//	
//	_vec<double> force_sub1 = (F/(pa->norm))*(pb->udloc-cosphi*pa->udloc);		//***---/***
//	_vec<double> force_sub2 = (F/(pb->norm))*(pa->udloc-cosphi*pb->udloc);		//***---/***
//	ene_pot += k_angl*(1-cosphi);//=k_angl*(1-(cosphi*cos(eq_angl)-sinphi*sin(eq_angl)));	//+-*
//	
//	pa->plus()->force += force_sub1;												//+++
//	pa->force -= force_sub1-force_sub2;											//------
//	pb->force -= force_sub2;													//---
//	
//}

void Frc(int myrank,int numprocs){
//	for(int i=0;i<list_size;i++){
//	for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){//計算領域内でセルリストを回す
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li(list+i);
		_ptcl* pa(NULL);
		_ptcl* pb(NULL);
		_ptcl* pa_p(NULL);
		_ptcl* pa_m(NULL);
		for(int j=0;pa=li->ptcl_index[j];j++){				//Lennard-Jones potential actin to actin in mycell
//pa=NULLに成るとforを抜ける．頭から粒子のポインターが入っていて最後がNULLに成っている
// 			if(pa->b_plus){
// 				pa_p=pa->b_plus->plus();
// 			}if(pa->b_minus){
// 				pa_m=pa->b_minus->minus;
// 			}
			if(pa->cnf->plus){
				pa_p=pa->cnf->plus;
			}if(pa->cnf->minus){
				pa_m=pa->cnf->minus;
			}
			for(int k=j+1;pb=li->ptcl_index[k];k++){		//Lennard-Jones potential actin to actin
				if((pa->cnf->feature!=0) 
					&& (pa_p==pb || pa_m==pb)
					){continue;}
				double potential;
				_vec<double> force_sub = LJ(pa,pb,&potential);
				pa->force += force_sub;						//+++
				pb->force -= force_sub;						//---
				ene_pot += potential;  						//+
			}                          
			pa_p = NULL;
			pa_m = NULL;
		}
		
		_list* li_sub = 0;
		
		for(int a=0;a<13;a++){//同種の粒子同士は半分でいい
#ifdef WALL_BOUNDARY
			if(li_sub=li->list_index[a]){
#else
			li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
			for(int j=0;pa=li->ptcl_index[j];j++){
// 				if(pa->b_plus){
// 					pa_p=pa->b_plus->plus();
// 				}if(pa->b_minus){
// 					pa_m=pa->b_minus->minus;
// 				}
				if(pa->cnf->plus){
					pa_p=pa->cnf->plus;
				}if(pa->cnf->minus){
					pa_m=pa->cnf->minus;
				}
				for(int k=0;pb=li_sub->ptcl_index[k];k++){		//Lennard-Jones potential actin to actin in enighbercell
					if((pa->cnf->feature!=0) 
						&& (pa_p==pb || pa_m==pb)
						){continue;}
					double potential;
					_vec<double> force_sub = LJ(pa,pb,&potential);
					pa->force += force_sub;						//+++
					pb->force -= force_sub;						//---
					ene_pot += potential;  						//+
				}                          
				pa_p = NULL;
				pa_m = NULL;
			}
#ifdef WALL_BOUNDARY
			}
#endif //WALL_BOUNDARY
		}
	}}}
	{
//		for(_Actin::_Bond* b=bond;b->plus;++b){
//			b->vec_cal();
//		}
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
//			int i(z + y*list_Z + x*(list_Y*list_Z));
			calc_norm_udloc( (list+(z + y*list_Z + x*(list_Y*list_Z))) );
//			_list* li = list+i;
//			_ptcl* pa = 0;
//			for(int j=0;pa=li->ptcl_index[j];j++){
//				if(pa->plus()){
////					calc_norm_udloc(pa);
//					_vec<double> dloc=pa->plus()->loc - pa->loc;		//---
//					distance_boundary(dloc);						
//					pa->norm = dloc.norm();							//***++sqrt
//					pa->udloc = dloc/pa->norm;						// ///
//				}
//			}
		}}}
#ifdef USE_MPI
		{
			int x1=(calcArea.nx!=0   )?(calcArea.x_pre-1):list_X-1;
			int x2=(calcArea.nx!=NX-1)?(calcArea.x_bck)  :0;
			int y1=(calcArea.ny!=0   )?(calcArea.y_pre-1):list_Y-1;
			int y2=(calcArea.ny!=NY-1)?(calcArea.y_bck)  :0;
			int z1=(calcArea.nz!=0   )?(calcArea.z_pre-1):list_Z-1;
			int z2=(calcArea.nz!=NZ-1)?(calcArea.z_bck)  :0;
if(NX>1){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				calc_norm_udloc( (list+(z + y*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z + y*list_Z + x2*(list_Y*list_Z))) );
			}}
}
if(NY>1){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				calc_norm_udloc( (list+(z + y1*list_Z + x*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z + y2*list_Z + x*(list_Y*list_Z))) );
			}}
}
if(NZ>1){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
				calc_norm_udloc( (list+(z1 + y*list_Z + x*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y*list_Z + x*(list_Y*list_Z))) );
			}}
}
if((NX-1)*(NY-1)!=0){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				calc_norm_udloc( (list+(z + y1*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z + y2*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z + y1*list_Z + x2*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z + y2*list_Z + x2*(list_Y*list_Z))) );
			}
}
if((NY-1)*(NZ-1)!=0){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
				calc_norm_udloc( (list+(z1 + y1*list_Z + x*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y1*list_Z + x*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z1 + y2*list_Z + x*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y2*list_Z + x*(list_Y*list_Z))) );
			}
}
if((NZ-1)*(NX-1)!=0){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
				calc_norm_udloc( (list+(z1 + y*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z1 + y*list_Z + x2*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y*list_Z + x2*(list_Y*list_Z))) );
			}
}
if((NX-1)*(NY-1)*(NZ-1)!=0){
			{
				calc_norm_udloc( (list+(z1 + y1*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y1*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z1 + y2*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y2*list_Z + x1*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z1 + y1*list_Z + x2*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y1*list_Z + x2*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z1 + y2*list_Z + x2*(list_Y*list_Z))) );
				calc_norm_udloc( (list+(z2 + y2*list_Z + x2*(list_Y*list_Z))) );
			}
}
		}


#endif
	}
	{
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li = list+i;
			_ptcl* pa = 0;
			_ptcl* pb = 0;
			for(int j=0;pa=li->ptcl_index[j];j++){//引張・曲げばね
				if(pa->cnf->plus){
//					linear_spring(pa);
					_vec<double> force_sub = -1.0*k_bane*(pa->norm-rb)*pa->udloc;	//-*****
					ene_pot += k_bane*(pa->norm-rb)*(pa->norm-rb)*0.5;				//--***
					pa->cnf->plus->force += force_sub;									//+++
					pa->force -= force_sub;											//---
					if(pb=pa->cnf->minus){
//						bend_spring(pa,pb);
						double cosphi = (pa->udloc) * (pb->udloc);					//***++ 5
			//			double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
			//			double d=acos(cosphi);
						//double F=	k_angl*(d -eq_angl)/sinphi;
						double F=	k_angl;//=k_angl*(sinphi*cos(eq_angl)-sin(eq_angl)*cosphi)/sinphi
						
						_vec<double> force_sub1 = (F/(pa->norm))*(pb->udloc-cosphi*pa->udloc);		//***---/*** 10
						_vec<double> force_sub2 = (F/(pb->norm))*(pa->udloc-cosphi*pb->udloc);		//***---/*** 10
						ene_pot += k_angl*(1-cosphi);//=k_angl*(1-(cosphi*cos(eq_angl)-sinphi*sin(eq_angl)));	//+-* 3
						
						pa->cnf->plus->force += force_sub1;												//+++ 3
						pa->force -= force_sub1-force_sub2;											//------ 6
						pb->force -= force_sub2;													//--- 3
					}
					Adhesion_Force(pa);//ねじりばね
				}
			}
		}}}
	}
// 	{
// 		_Actin::_Bond* ba;
// 		for(_Actin::_Bond* b=bond;b->plus;++b){
// 			_vec<double> force_sub = -1.0*k_bane*(b->norm-rb)*b->udloc;
// 			ene_pot += k_bane*(b->norm-rb)*(b->norm-rb)*0.5;
// 			b->plus->force += force_sub;
// 			b->minus->force -= force_sub;
// 			if(ba=b->b_plus){
// 				double cosphi = (b->udloc) * (ba->udloc);
// 	//			double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
// 	//			double d=acos(cosphi);
// 				//double F=	k_angl*(d -eq_angl)/sinphi;
// 				double F=	k_angl;//=k_angl*(sinphi*cos(eq_angl)-sin(eq_angl)*cosphi)/sinphi
// 				
// 				_vec<double> force_sub1 = F*(b->udloc-cosphi*ba->udloc)/(ba->norm);
// 				_vec<double> force_sub2 = F*(ba->udloc-cosphi*b->udloc)/(b->norm);
// 				ene_pot += k_angl*(1-cosphi);//=k_angl*(1-(cosphi*cos(eq_angl)-sinphi*sin(eq_angl)));
// 				
// 				ba->plus->force += force_sub1;
// 				ba->minus->force -= force_sub1-force_sub2;
// 				b->minus->force -= force_sub2;
// 			
// 				
// 			}
// 		}
// 	}
#ifdef WALL_BOUNDARY
	WALL::Wall_ptcl_frc(myrank);
#endif
#ifdef WALL_BIND
	WALL::bind_frc(myrank,numprocs);
#endif
	
#ifdef ARP2_3
	_ARP2_3::_ARP_ARP::Frc(myrank);
// if(t>0.001915)cout<<"frc5 ac-ar"<<endl;
	_ARP2_3::_ACTIN_ARP::Frc(myrank);
#endif
//	Adhesion_Force();
	
	
	
	return;
}

/*
void Frc(){
	for(int i=0;i<list_size;i++){
		_list* li = list+i;
		_ptcl* pa = 0;
		_ptcl* pb = 0;
		for(int j=0;pa=li->ptcl_index[j];j++){
			for(int k=j+1;pb=li->ptcl_index[k];k++){
				if(pa->feature==pb->feature && (pa->num-pb->num==1 || pa->num-pb->num==-1)){continue;}
				_vec<double> dloc = pa->loc - pb->loc;
				//distance_boundary(dloc);
				_vec<double> force_sub = refrce(dloc);
				pa->force += force_sub;
				pb->force -= force_sub;
				ene_pot += repot(dloc);
			}
		}
		
		_list* li_sub = 0;
		
		for(int a=0;a<13;a++){
#ifdef WALL_BOUNDARY
			if(li_sub=li->list_index[a]){
#else
			li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
			for(int j=0;pa=li->ptcl_index[j];j++){
				for(int k=0;pb=li_sub->ptcl_index[k];k++){
					if(pa->feature==pb->feature && (pa->num-pb->num==1 || pa->num-pb->num==-1)){continue;}
					_vec<double> dloc = pa->loc - pb->loc;
					distance_boundary(dloc);
					_vec<double> force_sub = refrce(dloc);
					pa->force += force_sub;
					pb->force -= force_sub;
					ene_pot += repot(dloc);
				}
			}
#ifdef WALL_BOUNDARY
			}
#endif //WALL_BOUNDARY
		}
	}
	for(int i = 1 ; numbering[i].num ; i++){
		if(numbering[i].num!=1){
			if(numbering[i].num == 2){
				_vec<double> dloc = numbering[i].ptcl_n[0]->loc - numbering[i].ptcl_n[1]->loc;
				distance_boundary(dloc);
				_vec<double> force_sub = frce_bane(dloc);
				numbering[i].ptcl_n[0]->force += force_sub;
				numbering[i].ptcl_n[1]->force -= force_sub;
				ene_pot += pot_bane(dloc);//ポテンシャル
				continue;
			}
			for(int j = 1 ; numbering[i].ptcl_n[j+1] ; j++){
				_vec<double> dloc1 = numbering[i].ptcl_n[j-1]->loc - numbering[i].ptcl_n[j]->loc;
				_vec<double> dloc2 = numbering[i].ptcl_n[j+1]->loc - numbering[i].ptcl_n[j]->loc;
				distance_boundary(dloc1);
				distance_boundary(dloc2);
				if(j==1){
					_vec<double> force_sub = frce_bane(dloc1);
					numbering[i].ptcl_n[j-1]->force += force_sub;
					numbering[i].ptcl_n[j]->force -= force_sub;
					ene_pot += pot_bane(dloc1);
				}
				_vec<double> force_sub_B = frce_bane(dloc2);
				numbering[i].ptcl_n[j]->force -= force_sub_B;
				numbering[i].ptcl_n[j+1]->force += force_sub_B;
				ene_pot += pot_bane(dloc2);
				
				double d = acos( (dloc1/dloc1.norm()) * (dloc2/dloc2.norm())) + 1.0E-9;
				
				_vec<double> force_sub1 = frce_angl( dloc1 , dloc2 , d );
				_vec<double> force_sub2 = frce_angl( dloc2 , dloc1 , d );
				
				numbering[i].ptcl_n[j-1]->force += force_sub1;
				numbering[i].ptcl_n[j]->force -= force_sub1+force_sub2;
				numbering[i].ptcl_n[j+1]->force += force_sub2;
				ene_pot += pot_angl( d );
			}
		}
	}
#ifdef WALL_BOUNDARY
	for(int i =0; i<Number_Ptcl;i++){
		WALL::Wall_ptcl_frc(*(ptcl + i));
	}
#endif
#ifdef WALL_BIND
	bind_frc();
#endif
	
#ifdef ARP2_3
	_ARP2_3::_ARP_ARP::Frc();
	_ARP2_3::_ACTIN_ARP::Frc();
#endif
//	Adhesion_Force();
	
	ptcl[0].force-=_vec<double>(0.0,0.0,0.263414634);
	ptcl[13].force+=_vec<double>(0.0,0.0,0.263414634);
	
	
	
	return;
}*/
#else //BROWNIAN
/*void Frc(){
	for(int i=0;i<list_size;i++){
		_list* li = list+i;
		_ptcl* pa = 0;
		_ptcl* pb = 0;
		for(int j=0;pa=li->ptcl_index[j];j++){
			for(int k=j+1;pb=li->ptcl_index[k];k++){
				_vec<double> dloc = pa->loc - pb->loc;
				//distance_boundary(dloc);
				_vec<double> force_sub = refrce(dloc);
				pa->frc_new += force_sub;
				pb->frc_new -= force_sub;
				ene_pot += repot(dloc);
			}
		}
		
		_list* li_sub = 0;
		
		for(int a=0;a<13;a++){
#ifdef WALL_BOUNDARY
			if(li_sub=li->list_index[a]){
#else
			li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
			for(int j=0;pa=li->ptcl_index[j];j++){
				for(int k=0;pb=li_sub->ptcl_index[k];k++){
					_vec<double> dloc = pa->loc - pb->loc;
					distance_boundary(dloc);
					_vec<double> force_sub = refrce(dloc);
					pa->frc_new += force_sub;
					pb->frc_new -= force_sub;
					ene_pot += repot(dloc);
				}
			}
#ifdef WALL_BOUNDARY
			}
#endif //WALL_BOUNDARY
		}
	}
	for(int i = 1 ; numbering[i].num ; i++){
		if(numbering[i].num!=1){
		if(numbering[i].num == 2){
			_vec<double> dloc = numbering[i].ptcl_n[0]->loc - numbering[i].ptcl_n[1]->loc;
			distance_boundary(dloc);
			numbering[i].ptcl_n[0]->frc_new += frce_bane(dloc);
			numbering[i].ptcl_n[1]->frc_new -= frce_bane(dloc);
			ene_pot += pot_bane(dloc);//ポテンシャル
			continue;
		}
		for(int j = 1 ; numbering[i].ptcl_n[j+1] ; j++){
			_vec<double> dloc1 = numbering[i].ptcl_n[j-1]->loc - numbering[i].ptcl_n[j]->loc;
			_vec<double> dloc2 = numbering[i].ptcl_n[j+1]->loc - numbering[i].ptcl_n[j]->loc;
			distance_boundary(dloc1);
			distance_boundary(dloc2);
			if(j==1){
				numbering[i].ptcl_n[j-1]->frc_new += frce_bane(dloc1);
				numbering[i].ptcl_n[j]->frc_new -= frce_bane(dloc1);
				ene_pot += pot_bane(dloc1);
			}
			numbering[i].ptcl_n[j]->frc_new -= frce_bane(dloc2);
			numbering[i].ptcl_n[j+1]->frc_new += frce_bane(dloc2);
			ene_pot += pot_bane(dloc2);
			
			double d = acos( (dloc1/dloc1.norm()) * (dloc2/dloc2.norm())) + 1.0E-9;
			numbering[i].ptcl_n[j-1]->frc_new += frce_angl( dloc1 , dloc2 , d );
			numbering[i].ptcl_n[j]->frc_new -= frce_angl( dloc1 , dloc2 , d )+frce_angl( dloc2 , dloc1 , d );
			numbering[i].ptcl_n[j+1]->frc_new += frce_angl( dloc2 , dloc1 , d );
			ene_pot += pot_angl( d );
		}}
	}
#ifdef WALL_BOUNDARY
	for(int i =0; i<Number_Ptcl;i++){
		Wall_ptcl_frc(*(ptcl + i));
	}
#endif
	//bind_frc();
	
	return;
}*/
/*void Frc(){
	for(int i=0 ; i<Number_Ptcl ; i++){
		_list* li = ptcl[i].li;
		_ptcl* p = 0;
		for(int j=1 ; p=li->index[j] ; j++){
//			_ptcl* p=li->index[j];
			if(p - &ptcl[i] > 0){
				_vec<double> dloc = ptcl[i].loc - p->loc;
				distance_boundary(dloc);
				if(dloc.norm()<=r_cut){
					ptcl[i].frc_new += refrce(dloc);
					p->frc_new -= refrce(dloc);
					ene_pot += repot(dloc);
				}
			}
		}
	}
	for(int i = 1 ; numbering[i].num ; i++){
		if(numbering[i].num!=1){
		if(numbering[i].num == 2){
			_vec<double> dloc = numbering[i].ptcl_n[0]->loc - numbering[i].ptcl_n[1]->loc;
			distance_boundary(dloc);
			numbering[i].ptcl_n[0]->frc_new += frce_bane(dloc);
			numbering[i].ptcl_n[1]->frc_new -= frce_bane(dloc);
			ene_pot += pot_bane(dloc);//ポテンシャル
			continue;
		}
		for(int j = 1 ; numbering[i].ptcl_n[j+1] ; j++){
			_vec<double> dloc1 = numbering[i].ptcl_n[j-1]->loc - numbering[i].ptcl_n[j]->loc;
			_vec<double> dloc2 = numbering[i].ptcl_n[j+1]->loc - numbering[i].ptcl_n[j]->loc;
			distance_boundary(dloc1);
			distance_boundary(dloc2);
			if(j==1){
				numbering[i].ptcl_n[j-1]->frc_new += frce_bane(dloc1);
				numbering[i].ptcl_n[j]->frc_new -= frce_bane(dloc1);
				ene_pot += pot_bane(dloc1);
			}
			numbering[i].ptcl_n[j]->frc_new -= frce_bane(dloc2);
			numbering[i].ptcl_n[j+1]->frc_new += frce_bane(dloc2);
			ene_pot += pot_bane(dloc2);
			
			double d = acos( (dloc1/dloc1.norm()) * (dloc2/dloc2.norm())) + 1.0E-9;
			numbering[i].ptcl_n[j-1]->frc_new += frce_angl( dloc1 , dloc2 , d );
			numbering[i].ptcl_n[j]->frc_new -= frce_angl( dloc1 , dloc2 , d )+frce_angl( dloc2 , dloc1 , d );
			numbering[i].ptcl_n[j+1]->frc_new += frce_angl( dloc2 , dloc1 , d );
			ene_pot += pot_angl( d );
		}}
	}
#ifdef WALL_BOUNDARY
	for(int i =0; i<Number_Ptcl;i++){
		WALL::Wall_ptcl_frc(*(ptcl + i));
	}
#endif
	return;
}*/
#endif
#ifdef BROWNIAN
// 
// void Ptcl_Adhesion_Location_FirstStep(){
// //	for(int i=0 ;i<Number_Ptcl;i++ ){
// //		_ptcl* p = (ptcl+i);
// //		p->quat_pre = p->quat;
// //		_vec<double> torque = p->adhesion_loc % p->adhesion_frc;
// //		_vec<double> omega = torque/Ptcl_Rot_Friction;
// //		
// //		p->quat += (0.25 * (omega * p->quat) * PTCL_dt);//dq/dt=0.5*ω*q
// //		p->quat /= (p->quat.norm());
// //		
// //		p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc);
// //		p->adhesion_frc.IN(0.0,0.0,0.0);
// //	}
// 	_ptcl* p;
// 	for(int i=1; numbering[i].num ;i++){
// 		for(int j=0;p=numbering[i].ptcl_n[j];j++){
// 			p->quat_pre = p->quat;
// 			_vec<double> torque = (p->adhesion_loc % p->adhesion_frc);
// 			_vec<double> omega = (torque/Ptcl_Rot_Friction);
// 			
// 			p->quat += (0.25 * (omega * p->quat) * PTCL_dt);//dq/dt=0.5*ω*q
// 			p->quat /= (p->quat.norm());
// 			
// 			p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc_inite);
// 			p->adhesion_frc.IN(0.0,0.0,0.0);
// 		}
// 	}
// 	return;
// }
// void Ptcl_Adhesion_Location_SecondStep(){
// //	static const double amplitude=sqrt(2.0*SYS_kT*Ptcl_Rot_Friction/PTCL_dt);
// //	static const double root3=sqrt(3.0)*2.0;
// //	_vec<double> RandamF(0.0,0.0,0.0);
// //	for(int i=0 ;i<Number_Ptcl;i++ ){
// //		RandamF=amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3;
// //		_ptcl* p = (ptcl+i);
// //		_vec<double> torque = p->adhesion_loc % p->adhesion_frc + RandamF;
// //		_vec<double> omega = torque/Ptcl_Rot_Friction;
// //	//	if(torque.norm()*0.5>=RandamF.norm())cout<<"torque"<<endl;
// //	//	else cout<<"randam"<<endl;
// //		
// //		p->quat = p->quat_pre + (0.5 * (omega * p->quat) * PTCL_dt);//dq/dt=0.5*ω*q
// //		p->quat /= (p->quat.norm());
// //		
// //		p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc);
// //		p->adhesion_frc.IN(0.0,0.0,0.0);
// //	}
// 	static const double amplitude=sqrt(2.0*SYS_kT*Ptcl_Rot_Friction/PTCL_dt);
// 	static const double root3=sqrt(3.0)*2.0;
// 	_vec<double> RandamF(0.0,0.0,0.0);
// 	_ptcl* p;
// 	for(int i=1; numbering[i].num ;i++){
// 		for(int j=0;p=numbering[i].ptcl_n[j];j++){
// 			RandamF=amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3;
// 			_vec<double> torque = ((p->adhesion_loc % p->adhesion_frc) + RandamF);
// 			_vec<double> omega = torque/Ptcl_Rot_Friction;
// 
// 			p->quat = (p->quat_pre + (0.5 * (omega * p->quat) * PTCL_dt));//dq/dt=0.5*ω*q
// 			p->quat /= (p->quat.norm());
// 			
// 			p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc_inite);
// 			p->adhesion_frc.IN(0.0,0.0,0.0);
// 		}
// 	}
// 	return;
// }
// void Ptcl_Location_FirstStep(){
// 	for(int i=0 ;i<Number_Ptcl;i++ ){
// 		_ptcl* p = (ptcl+i);
// 		p->loc+= (0.5*PTCL_dt*( p->force + p->adhesion_frc )/Ptcl_Friction );
// 		p->force.IN(0.0,0.0,0.0);
// 		Periodic(p->loc);
// 		_list* li=p->list_calc();
// 		if(li!=p->li){
// 			p->li->rm(p,p->list_num);
// 			p->add(li);
// 			li->add(p);
// // 			p->li=li;
// // 			p->list_num=li->num;
// 		}
// 	}
// 	Ptcl_Adhesion_Location_FirstStep();
// #ifdef ARP2_3
// 	_ARP2_3::Arp_Location_FirstStep();
// #endif
// //	list_cal();
// 	return;
// }
// void Ptcl_Location_SecondStep(){
// 	static const double amplitude=sqrt(2.0*SYS_kT*Ptcl_Friction/PTCL_dt);
// 	static const double root3=sqrt(3.0)*2.0;
// 	_vec<double> RandamF(0.0,0.0,0.0);
// 	for(int i=0 ;i<Number_Ptcl;i++ ){
// 		RandamF= (_vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * (root3 * amplitude));
// 		_ptcl* p = (ptcl+i);
// // 		if(i<20)cout<<i<<'\t'<<p->loc<<"\t"<<p->loc_pre<<endl;
// // 		if(i<20)cout<<i<<'\t'<<p->force<<"\t"<<p->adhesion_frc<<"\t"<<RandamF<<endl;
// 		p->loc_pre+= ( PTCL_dt * (p->force + p->adhesion_frc + RandamF) / Ptcl_Friction );
// //		average_randam_force+=PTCL_dt * RandamF;
// //		veriance_randam_force+=(_vec<double>(RandamF.x*RandamF.x ,RandamF.y*RandamF.y ,RandamF.z*RandamF.z)*PTCL_dt*PTCL_dt);
// // 		_vec<double> tmp= Inite_X-p->loc;
// // 		veriance_X +=_vec<double>(tmp.x*tmp.x ,tmp.y*tmp.y ,tmp.z*tmp.z);
// 		p->force.IN(0.0,0.0,0.0);
// 		Periodic(p->loc_pre);
// 		p->loc=p->loc_pre;
// 		_list* li=p->list_calc();
// 		if(li!=p->li){
// 			p->li->rm(p,p->list_num);
// 			p->add(li);
// 			li->add(p);
// 		}
// 	}
// 	Ptcl_Adhesion_Location_SecondStep();
// #ifdef ARP2_3
// 	_ARP2_3::Arp_Location_SecondStep();
// #endif
// //	list_cal();
// 	
// // 	{
// // 		
// // 		veriance_X/=Number_Ptcl;
// // 		stringstream ss;
// // 		ss.str("");
// // 		ss<<OUTPUT_N<<"veriance_x.dat";
// // 		ofstream fout_X(ss.str().c_str(),ios::out | ios::app);
// // 		fout_X<<t<<"\t\t"<<veriance_X<<endl;
// // 		fout_X.close();
// // 		
// // 		veriance_X.IN(0.0,0.0,0.0);
// // 		
// // // 		ss.str("");
// // // 		ss<<OUTPUT_N<<"veriance_x_binary";
// // // 		ofstream fout_X_bin(ss.str().c_str(),ios::out | ios::app | ios::binary);
// // // 		fout_X_bin.write((char*)&t,sizeof(double));
// // // 		fout_X_bin.write((char*)&tmp,sizeof(_vec<double>));
// // // 		fout_X_bin.write((char*)&veriance_X,sizeof(_vec<double>));
// // // 		fout_X_bin.close();
// // 		
// // 	}
// 	
// 	return;
// }
void Ptcl_Adhesion_Location_FirstStep(_ptcl* p){
			p->quat_pre = p->quat;
			_vec<double> torque = (p->adhesion_loc % p->adhesion_frc);
			_vec<double> omega = (torque/Ptcl_Rot_Friction);
			
			p->quat += (0.25 * (omega * p->quat) * PTCL_dt);//dq/dt=0.5*ω*q
			p->quat /= (p->quat.norm());
			
			p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc_inite);
			p->adhesion_frc.IN(0.0,0.0,0.0);
	return;
}
void Ptcl_Adhesion_Location_SecondStep(_ptcl* p){
	static const double amplitude=sqrt(2.0*SYS_kT*Ptcl_Rot_Friction/PTCL_dt);
	static const double root3=sqrt(3.0)*2.0;
	_vec<double> RandamF(0.0,0.0,0.0);
			RandamF=amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3;
			_vec<double> torque = ((p->adhesion_loc % p->adhesion_frc) + RandamF);
			_vec<double> omega = torque/Ptcl_Rot_Friction;

			p->quat = (p->quat_pre + (0.5 * (omega * p->quat) * PTCL_dt));//dq/dt=0.5*ω*q
			p->quat /= (p->quat.norm());
			
			p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc_inite);
			p->adhesion_frc.IN(0.0,0.0,0.0);
	return;
}
void Ptcl_Location_FirstStep(const int myrank,const int numprocs){
//	for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li_sub = list+i;
		_ptcl* p =NULL;
		for(int j=0;p=li_sub->ptcl_index[j];j++){
			p->loc+= (0.5*PTCL_dt*( p->force + p->adhesion_frc )/Ptcl_Friction );
			p->force.IN(0.0,0.0,0.0);
			Periodic(p->loc);
			if(p->cnf->feature)Ptcl_Adhesion_Location_FirstStep(p);
		}
	}}}
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){//list更新
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li_sub = list+i;
		_ptcl* p =NULL;
		for(int j=0;p=li_sub->ptcl_index[j];j++){
			_list* li=p->list_calc();//list計算どこのlistか
			if(li!=p->cnf->li){//セルを出て行った場合
#ifdef USE_MPI
				if(li->myrank!=-1){		//糊代にいれば転送リストに
if(NX>1){
					if((li->myrank/9)==2 && (p->cnf->li->myrank/9)==1){		send_list_ptcl_x_bck_inner_new.push_back(p);
					}else if((li->myrank/9)==0 && (p->cnf->li->myrank/9)==1){	send_list_ptcl_x_pre_inner_new.push_back(p);}
}
if(NY>1){
					if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)==1){		send_list_ptcl_y_bck_inner_new.push_back(p);
					}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)==1){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
					if((li->myrank%3)==2 && (p->cnf->li->myrank%3)==1){		send_list_ptcl_z_bck_inner_new.push_back(p);
					}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)==1){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
				}
#endif//USE_MPI
				p->cnf->li->rm(p,p->cnf->list_num);
				p->add(li);
				li->add(p);
				j--;
			}
		}
	}}}
//	Ptcl_Adhesion_Location_FirstStep();
#ifdef ARP2_3
	_ARP2_3::Arp_Location_FirstStep(myrank,numprocs);
#endif
#ifdef WALL_BOUNDARY
	WALL::Wall_Location_FirstStep(myrank,numprocs);
#endif// WALL_BOUNDARY
//	list_cal();
	return;
}
void Ptcl_Location_SecondStep(const int myrank,const int numprocs){
#ifdef WALL_BOUNDARY
	WALL::Wall_Location_SecondStep(myrank,numprocs);
#endif// WALL_BOUNDARY
#ifdef USE_MPI
	mympi::mympi_flow_wallinf_bcast(myrank,numprocs);
#endif// USE_MPI
	static const double amplitude=sqrt(2.0*SYS_kT*Ptcl_Friction/PTCL_dt);
	static const double root3=sqrt(3.0)*2.0;
	_vec<double> RandamF(0.0,0.0,0.0);
//	for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li_sub = list+i;
		_ptcl* p =NULL;
		for(int j=0;p=li_sub->ptcl_index[j];j++){
			RandamF= (_vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * (root3 * amplitude));
			p->loc_pre+= ( PTCL_dt * (p->force + p->adhesion_frc + RandamF) / Ptcl_Friction );
			
			flow.IN(p->loc,((p->force + p->adhesion_frc + RandamF) / Ptcl_Friction),p->cnf->feature);
			
			p->force.IN(0.0,0.0,0.0);
			Periodic(p->loc_pre);
			p->loc=p->loc_pre;
			
			if(p->cnf->feature)Ptcl_Adhesion_Location_SecondStep(p);
		}
	}}}
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li_sub = list+i;
		_ptcl* p =NULL;
		for(int j=0;p=li_sub->ptcl_index[j];j++){
			
			_list* li=p->list_calc();
			if(li!=p->cnf->li){
#ifdef USE_MPI
				if(li->myrank!=-1){		//糊代にいれば転送リストに
if(NX>1){
					if((li->myrank/9)==2 && (p->cnf->li->myrank/9)==1){		send_list_ptcl_x_bck_inner_new.push_back(p);
					}else if((li->myrank/9)==0 && (p->cnf->li->myrank/9)==1){	send_list_ptcl_x_pre_inner_new.push_back(p);}
}
if(NY>1){
					if(((li->myrank/3)%3)==2 && ((p->cnf->li->myrank/3)%3)==1){		send_list_ptcl_y_bck_inner_new.push_back(p);
					}else if(((li->myrank/3)%3)==0 && ((p->cnf->li->myrank/3)%3)==1){	send_list_ptcl_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
					if((li->myrank%3)==2 && (p->cnf->li->myrank%3)==1){		send_list_ptcl_z_bck_inner_new.push_back(p);
					}else if((li->myrank%3)==0 && (p->cnf->li->myrank%3)==1){	send_list_ptcl_z_pre_inner_new.push_back(p);}
}
				}
#endif
				p->cnf->li->rm(p,p->cnf->list_num);
				p->add(li);
				li->add(p);
				j--;
			}
		}
	}}}
//////debug
////	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
////	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
////	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
////		int i(z + y*list_Z + x*(list_Y*list_Z));
////		_list* li_sub = list+i;
////		_ptcl* p =NULL;
////		for(int j=0;p=li_sub->ptcl_index[j];j++){
////			
////			_list* li=p->list_calc();
////			if(li!=p->cnf->li){
////				cout<<"list error "<<endl;exit(0);
////			}
////		}
////	}}}
#ifdef ARP2_3
	_ARP2_3::Arp_Location_SecondStep(myrank,numprocs);
#endif
	return;
}
void ene_format(){//初期化
	ene_kine=0.0;
	ene_pot=0.0;
	_ARP2_3::ene_pot=0.0;
#ifdef WALL_BOUNDARY
	WALL::ene_pot=0.0;
#endif// WALL_BOUNDARY
	p.IN(0.0,0.0,0.0);
	return;
}

// void p_level_cal(){
// 	for(int i=0 ;i<Number_Ptcl;i++ ){
// 		if(t-ptcl[i].time>p_time){
// 			if(ptcl[i].feature){ptcl[i].p_level=0;}
// 			else{ptcl[i].p_level=1;}
// 		}
// 	}
// 	return;
// }
#else //BROWNIAN
/*inline double kinetic( double m , _vec<double> v){
	return 1.0/2.0*m*v.sqr();
}
inline _vec<double> momentum( double m , _vec<double> v){
	return m*v;
}
void PTCL_energy(){
	for(int i=0 ; i<Number_Ptcl ; i++){
		ene_kine += kinetic(ptcl[i].mass , ptcl[i].vel);//運動エネルギ
		p += momentum(ptcl[i].mass , ptcl[i].vel);//運動量
	}
	ene = ene_kine + ene_pot;//エネルギ
	return;
}
inline _vec<double> Ptcl_location_sub(_ptcl* p){//calculus of finite differences with second-order
	return p->loc+PTCL_dt*(p->vel+PTCL_dt/(2.0*p->mass)*p->frc_new);
}
void PTCL_Location(){
	for(int i=0 ; i<Number_Ptcl ; i++){
		_vec<double> loc = Ptcl_location_sub(ptcl+i);
#ifdef WALL_BOUNDARY
		WALL::Wall_ptcl(ptcl[i],loc);
#endif
#ifndef WALL_BOUNDARY
		ptcl[i].loc = loc;
		Periodic( ptcl[i].loc );
#endif
		ptcl[i].frc_old = ptcl[i].frc_new;
	}
	list_cal();
	return;
}
inline _vec<double> Ptcl_velocity_sub(_ptcl* p){//calculus of finite differences with second-order
	return p->vel + PTCL_dt/(2.0*p->mass)*(p->frc_old + p->frc_new);
}
void PTCL_Velocity(){//速さ
	for(int i=0;i < Number_Ptcl; i++){
		ptcl[i].vel = Ptcl_velocity_sub(ptcl+i);
	}
	return;
}
void ene_format(){//初期化
	for( int i=0 ; i<Number_Ptcl ; i++){
		ptcl[i].frc_new = _vec<double>(0.0,0.0,0.0);
	}
	ene_kine=0.0;
	ene_pot=0.0;
	p=_vec<double>(0.0,0.0,0.0);
	return;
}
void p_level_cal(){
	for(int i=0 ;i<Number_Ptcl;i++ ){
		if(t-ptcl[i].time>p_time){
			if(ptcl[i].feature){ptcl[i].p_level=0;}
			else{ptcl[i].p_level=1;}
		}
	}
	return;
}*/
#endif //BROWNIAN

namespace _Actin{

#ifdef DEPOLYMERIZATION
void depolymerization(const int myrank,const int numprocs){
	double probability = probability_deporimerize;
	_ptcl* p=NULL;
	for(int i=0 ; p=minus_end->end(i) ; i++){
#ifdef USE_MPI
// 	if(p->cnf->li->myrank==myrank || p->cnf->li->myrank==myrank+numprocs || p->cnf->li->myrank==myrank+2*numprocs ){
	if(p->loc.x>=(double)calcArea.nx*SYS_X/NX && p->loc.x<(calcArea.nx+1.0)*SYS_X/NX
	&& p->loc.y>=(double)calcArea.ny*SYS_Y/NY && p->loc.y<(calcArea.ny+1.0)*SYS_Y/NY
	&& p->loc.z>=(double)calcArea.nz*SYS_Z/NZ && p->loc.z<(calcArea.nz+1.0)*SYS_Z/NZ){
#endif
		double a = RAND;
		if(a<probability 
//			&& p->p_level==0 
			&& p->cnf->bind==-1)
		{
			counter_depol++;
			minus_end->depolymerization(p,i);
			p->cnf->feature = 0;
// 			p->time = t;
#ifdef USE_MPI
			{
				mympi::_conformation_change cc(depol,p->index,i);
				mympi::bit_sender.pack(&cc);
			}
#endif//USE_MPI
		}
#ifdef USE_MPI
	}
#endif
	}
	return;
}
#endif //DEPOLYMERIZATION
#ifdef POLYMERIZATION
void polymerization_mono(_ptcl* p_g1, _ptcl* p_g2 ,const int myrank,const int numprocs){
	
	plus_end->polymerization_mono(p_g2,p_g1);
	
	p_g1->cnf->feature = 1;
// 	p_g1->time = t;
	p_g2->cnf->feature = 1;
// 	p_g2->time = t;
	
	_vec<double> delloc = p_g2->loc - p_g1->loc;	distance_boundary(delloc);	delloc/=delloc.norm();
	_vec<double> JIKU(RAND-0.5,RAND-0.5,RAND-0.5);
	while(1){
		while( JIKU.norm() <1.0e-6){
			JIKU.x = RAND-0.5;
			JIKU.y = RAND-0.5;
			JIKU.z = RAND-0.5;
		}
		JIKU/=JIKU.norm();
		JIKU=JIKU%delloc;
		if(JIKU.norm()>1.0e-6){
			JIKU/=JIKU.norm();
			break;
		}
	}
	JIKU*=0.5 * sigma_LJ;
//	p_g1->adhesion_loc.IN(delloc.y-delloc.z , delloc.z-delloc.x , delloc.x-delloc.y);
//	p_g1->adhesion_loc = (p_g1->adhesion_loc / p_g1->adhesion_loc.norm()) * 0.5 * sigma_LJ;
	p_g1->adhesion_loc.IN(JIKU);
	p_g1->adhesion_loc_inite = p_g1->adhesion_loc;
	p_g1->quat.IN(0.0 ,0.0 ,0.0);
	p_g1->quat_pre.IN(0.0 ,0.0 ,0.0);
	_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
	quat/=quat.norm();
	p_g2->adhesion_loc = _quaternion::rot( quat , p_g1->adhesion_loc);
	p_g2->adhesion_loc_inite = p_g2->adhesion_loc;
	p_g2->quat.IN(0.0 ,0.0 ,0.0);
	p_g2->quat_pre.IN(0.0 ,0.0 ,0.0);
#ifdef USE_MPI
	{
		mympi::bit_sender.pack(p_g1);
		mympi::bit_sender.pack(p_g2);
	}
#endif//USE_MPI
	
	return;
}
void nucleation(const int myrank,const int numprocs){
// 	if(numbering[0].num>1){
		int check_list=0;
		_ptcl* list_up[100][2];
		_ptcl* pa = 0;
		_ptcl* pb = 0;
		_list* lia=0;
		_list* lib=0;
//		for(int x=myrank*calc_length_X;x<(myrank+1)*calc_length_X;x++){
//// 		for(int x=0;x<list_X;x++){
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<list_Z;z++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			lia = list+(z + y*list_Z + x*(list_Y*list_Z));
// 		for(int i=0;i<list_size;i++){
// 			lia=list+i;
			for(int j=0;pa=lia->ptcl_index[j];j++){
//				if(pa->feature==0 && pa->p_level==1){
				if(pa->cnf->feature==0){
					for(int k=j+1;pb=lia->ptcl_index[k];k++){
//						if(pb->cnf->feature==0 && pb->p_level==1){
						if(pb->cnf->feature==0){
							_vec<double> delloc = pa->loc - pb->loc;
#ifndef BROWNIAN
//							_vec<double> delvel = pa->vel - pb->vel;
#endif
							distance_boundary(delloc);
							
							if(delloc.norm() < r_por_g_g 
#ifndef BROWNIAN
//							&& delvel.sqr() < r1 
#endif
							){
								double probability = probability_porimerize_p_p;
								double c = RAND;
								
								if(c<probability){
									list_up[check_list][0]=pa;
									list_up[check_list][1]=pb;
									check_list++;
								}
							}
						}
					}
				}
			}
			for(int in=0;in<13;in++){
				if(lib=lia->list_index[in]){
					for(int j=0;pa=lia->ptcl_index[j];j++){
//						if(pa->cnf->feature==0 && pa->p_level==1){
						if(pa->cnf->feature==0){
							for(int k=0;pb=lib->ptcl_index[k];k++){
//								if(pb->cnf->feature==0 && pb->p_level==1){
								if(pb->cnf->feature==0){
									_vec<double> delloc = pa->loc - pb->loc;
#ifndef BROWNIAN
//									_vec<double> delvel = pa->vel - pb->vel;
#endif
									distance_boundary(delloc);
									
									if(delloc.norm() < r_por_g_g 
#ifndef BROWNIAN
//									&& delvel.sqr() < r1 
#endif
									){
										double probability = probability_porimerize_p_p;
										double c = RAND;
										
										if(c<probability){
											list_up[check_list][0]=pa;
											list_up[check_list][1]=pb;
											check_list++;
										}
									}
								}
							}
						}
					}
				}
			}
		}}}
		if(check_list==1){
			polymerization_mono( list_up[0][0] , list_up[0][1] , myrank, numprocs);
			counter_pol++;
		}
		else if(check_list!=0){
			_ptcl* pc=0;
			_ptcl* pd=0;
			for(int ia=0;ia<check_list;ia++){
				int check_kaburi_num=0;
				int check_kaburi[50];
				for(int ib=0;ib<2;ib++){
					pc=list_up[ia][ib];
					for(int ja=ia+1;ja<check_list;ja++){
						for(int jb=0;jb<2;jb++){
							pd=list_up[ja][jb];
							if(pc==pd){
								check_kaburi_num++;
								check_kaburi[0]=ia;
								check_kaburi[check_kaburi_num]=ja;
							}
						}
					}
				}
				if(!check_kaburi_num){
					if(list_up[ia][0]->cnf->feature==0 &&list_up[ia][1]->cnf->feature==0 ){
						polymerization_mono( list_up[ia][0] , list_up[ia][1] , myrank, numprocs);
						counter_pol++;
					}
				}
				else{
					int index=(int)std::floor((check_kaburi_num+1)*RAND);
					while(index==(check_kaburi_num+1))index=(int)std::floor((check_kaburi_num+1)*RAND);
					if(list_up[index][0]->cnf->feature==0 &&list_up[index][1]->cnf->feature==0 ){
						int index_sub=check_kaburi[index];
						polymerization_mono( list_up[index_sub][0] , list_up[index_sub][1] , myrank, numprocs);
						counter_pol++;
					}
				}
			}
		}
	return;
}

void polymerization_fila(int i,_ptcl* p_f,_ptcl* p_g ,const int myrank,const int numprocs){
	plus_end->polymerization_fila(i,p_g);
	p_g->cnf->feature = 1;
// 	p_g->time = t;
	p_f->cnf->feature = 2;
	_vec<double> delloc = p_f->loc - p_g->loc;
	distance_boundary(delloc);	delloc/=delloc.norm();
	_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
	quat/=quat.norm();
	p_g->adhesion_loc = _quaternion::rot( quat , p_f->adhesion_loc);
	p_g->adhesion_loc_inite = p_g->adhesion_loc;
	p_g->quat.IN(0.0 ,0.0 ,0.0);
	p_g->quat_pre.IN(0.0 ,0.0 ,0.0);
#ifdef USE_MPI
	{
		conf_ch cc=pol_fila;
		mympi::bit_sender.pack(&cc);
		mympi::bit_sender.pack(&i);
		mympi::bit_sender.pack(&p_f->index);
		mympi::bit_sender.pack(&p_g->index);
		mympi::bit_sender.pack(&p_g->adhesion_loc);
	}
#endif//USE_MPI
	
	
	PolymerizationInfo.add(p_f->index , p_g->index);
	
	return;
}
void polymerization(const int myrank,const int numprocs){
	_ptcl* pa=NULL;
	for(int i=0;pa=plus_end->end(i);++i){//フィラメント端でループをまわす
#ifdef USE_MPI
//	if(pa->cnf->li->myrank==myrank || pa->cnf->li->myrank==myrank+numprocs || pa->cnf->li->myrank==myrank+2*numprocs ){
	if(pa->loc.x>=(double)calcArea.nx*SYS_X/NX && pa->loc.x<(calcArea.nx+1.0)*SYS_X/NX
	&& pa->loc.y>=(double)calcArea.ny*SYS_Y/NY && pa->loc.y<(calcArea.ny+1.0)*SYS_Y/NY
	&& pa->loc.z>=(double)calcArea.nz*SYS_Z/NZ && pa->loc.z<(calcArea.nz+1.0)*SYS_Z/NZ){//計算領域内かどうか//自分が計算するかどうか
#endif
		memset(&list_sub , 0 , sizeof(_list_sub));//十号の可能性のある粒子listの初期化
		if(pa->cnf->bind!=-1)continue;//bind==-1ならば自由端
		_list* li = pa->cnf->li;
		_list* li_sub = 0;
		_ptcl* pb=0;
		for(int in=0;in<27;in++){
			if(li_sub=li->list_index[in]){
				for(int k=0 ; pb=li_sub->ptcl_index[k] ; k++){
//					if(pb->cnf->feature == 0 && pb->p_level==1 ){
					if(pb->cnf->feature == 0){//0=モノマー，1=フィラメント
						_vec<double> delloc = pb->loc - pa->loc;
#ifndef BROWNIAN
//						_vec<double> delvel = pb->vel - pa->vel;
#endif
						distance_boundary(delloc);
						if(delloc.norm() < r_por_f_g //距離判定
#ifndef BROWNIAN
//						&& delvel.norm()*delvel.norm() < r1 
#endif
						){
							_vec<double> delloc2;
							{//フィラメントの軸の計算
								if(pa->cnf->minus){//一個後ろに粒子があれば
									delloc2 = pa->loc - pa->cnf->minus->loc;
								}
								else {//arp23がフィラメント端(マイナス)になっていて，かつ，arp23にアクチンが一つしか付いていないとき
									delloc2 = pa->loc - pa->cnf->end_arp2_3->loc;
								}
							}
							distance_boundary(delloc2);
// 							double delangl = acos( (delloc/delloc.norm()) * (delloc2/delloc2.norm())) + 1.0E-9;
// 							if(delangl < por_angl){
							if((delloc/delloc.norm()) * (delloc2/delloc2.norm()) > 0.5){//前方60度．ハードコーディングに気づかなかった
								double probability = probability_porimerize_f_p;
								double c = RAND;
								if(c<probability){	//結合確率の判定
									list_sub.num +=1;
									list_sub.index[list_sub.num]=pb;
								}
							}
						}
					}
				}
			}
		}
		if(list_sub.num==1){//重合するものが一個だった場合
			polymerization_fila(i, pa , list_sub.index[1] , myrank, numprocs);
			counter_pol++;
		}
		else if(list_sub.num>1){//重合するものが複数だった場合，乱数で選ぶ
			int num = (int)ceil(list_sub.num*RAND);
			while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
			polymerization_fila(i, pa , list_sub.index[num] , myrank, numprocs);
			counter_pol++;
		}
#ifdef USE_MPI
	}
#endif//USE_MPI
	}
}
#endif //POLYMERIZATION
#ifdef SEVERING
void severing_sub(_ptcl* pa,_ptcl* pb,const int myrank,const int numprocs){
	if(pb->plus()){
		if(pa->minus()){
			plus_end->add(pa);
			minus_end->add(pb);
			pa->cnf->feature=1;
			pb->cnf->feature=1;
		}else{
			minus_end->add(pb);
			pb->cnf->feature=1;
			pa->cnf->feature=0;
// 			pa->time = t;
		}
	}else{
		if(pa->minus()){
			plus_end->add(pa);
			pb->cnf->feature=0;
// 			pb->time = t;
			pa->cnf->feature=1;
		}else{
			pb->cnf->feature=0;
// 			pb->time = t;
			pa->cnf->feature=0;
// 			pa->time = t;
		}
	}
	pa->plus_in(-1);//	pa->plus_in(NULL);
	pb->minus_in(-1);//	pb->minus_in(NULL);
	
	return;
}
inline double sev_pro(_ptcl* pa , _ptcl* pb){
	//_vec<double> delloc = pa->loc - pb->loc;
	//distance_boundary(delloc);
	//return pro_sev*exp(-1.0*(delloc.norm()/l_eq-1)/e_c);
	return pro_sev;
}
inline double sev_pro(_Actin::_Bond*){
	//_vec<double> delloc = pa->loc - pb->loc;
	//distance_boundary(delloc);
	//return pro_sev*exp(-1.0*(delloc.norm()/l_eq-1)/e_c);
	return pro_sev;
}
void severing(const int myrank,const int numprocs){
//	for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
		int i(z + y*list_Z + x*(list_Y*list_Z));
		_list* li = list+i;
		_ptcl* pa = 0;
		for(int j=0;pa=li->ptcl_index[j];j++){
			if(pa->plus()){
				double probability=sev_pro(b);
				double c = RAND;
				if(c<probability){
					severing_sub(pa,pa->plus(),myrank,numprocs);
				}
			}
		}
	}}}
// 	for(_Actin::_Bond* b=bond;b->plus;++b){
// 		double probability=sev_pro(b);
// 		double c = RAND;
// 		if(c<probability){
// 			severing_sub(b,myrank,numprocs);
// 		}
// 	}
	return;
}
#endif //SEVERING

}// namespace _Actin{

// #endif //BROWNIAN
// #ifdef DEPOLYMERIZATION
// void depolymerization(){
// 	double probability = probability_deporimerize;
// 	for(int i = 1 ; numbering[i].num ; i++){
// 		if(numbering[i].num!=1){
// 			
// 			double a = RAND;
// 			if(a<probability && numbering[i].ptcl_n[0]->p_level==0 && numbering[i].ptcl_n[0]->bind==-1){
// 				counter_depol++;
// 				if(numbering[i].num==2){
// 					for(int j=0 ; numbering[i].ptcl_n[j] ; j++){
// 						numbering[i].ptcl_n[j]->cnf->feature = 0;
// 						numbering[i].ptcl_n[j]->time = t;
// 						numbering[i].ptcl_n[j]->num = numbering[0].num + j;
// 						numbering[0].ptcl_n[numbering[0].num + j]=numbering[i].ptcl_n[j];
// 						numbering[i].ptcl_n[j]=0;
// 					}
// 					numbering[i].num=0;
// 					numbering[0].num+=2;
// 					for(int j=i+1 ; numbering[j].num ; j++){
// 						numbering[j-1]=numbering[j];
// 						memset(&numbering[j],0,sizeof(_numbering));
// 						for(int k=0 ; numbering[j-1].ptcl_n[k] ; k++){
// 							numbering[j-1].ptcl_n[k]->cnf->feature = j-1;
// 						}
// 					}
// 				//	if(no_fila - i!=0){
// 				//		memmove(&numbering[i],&numbering[i+1],(no_fila - i)*sizeof(_numbering));
// 				//		memset(&numbering[no_fila],0,sizeof(_numbering));
// 				//	}
// 					no_fila-=1;
// 				}
// 				else{
// 					_ptcl* p = numbering[i].ptcl_n[0];
// 					p->cnf->feature = 0;
// 					p->time = t;
// 					p->num = numbering[0].num;
// 					numbering[0].ptcl_n[numbering[0].num] = p;
// 					numbering[0].num++;
// 					numbering[i].ptcl_n[0]=0;
// 					numbering[i].sort_n(0);
// 					numbering[i].num--;
// 				}
// 			}
// 		}
// 	}
// 	return;
// }
// #endif //DEPOLYMERIZATION
// #ifdef POLYMERIZATION
// void polymerization_fila(_ptcl* p_f, _ptcl* p_g ){
// 	numbering[0].sort_n(p_g->num);
// 	numbering[0].num--;
// 	numbering[p_f->cnf->feature].ptcl_n[p_f->num +1] = p_g;
// 	numbering[p_f->cnf->feature].num++;
// 	p_g->cnf->feature = p_f->cnf->feature;
// 	p_g->num = p_f->num +1;
// 	p_g->time = t;
// 	
// 	//p_g->quat = p_f->quat;
// 	//p_g->quat_pre = p_g->quat;
// 	_vec<double> delloc = p_f->loc - p_g->loc;	distance_boundary(delloc);	delloc/=delloc.norm();
// 	_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
// 	quat/=quat.norm();
// 	p_g->adhesion_loc = _quaternion::rot( quat , p_f->adhesion_loc);
// 	p_g->adhesion_loc_inite = p_g->adhesion_loc;
// 	p_g->quat.IN(0.0 ,0.0 ,0.0);
// 	p_g->quat_pre.IN(0.0 ,0.0 ,0.0);
// 	
// 	return;
// }
// void polymerization_mono(_ptcl* p_g1, _ptcl* p_g2){
// 	numbering[0].sort_n(p_g1->num);
// 	numbering[0].sort_n(p_g2->num);
// 	numbering[0].num-=2;
// 	no_fila++;
// 	p_g1->cnf->feature = no_fila;
// 	p_g1->num = 0;
// 	p_g1->time = t;
// 	p_g2->cnf->feature = no_fila;
// 	p_g2->num = 1;
// 	p_g2->time = t;
// 	numbering[no_fila].ptcl_n[0]=p_g1;
// 	numbering[no_fila].ptcl_n[1]=p_g2;
// 	numbering[no_fila].num=2;
// 	
// 	_vec<double> delloc = p_g2->loc - p_g1->loc;	distance_boundary(delloc);	delloc/=delloc.norm();
// 	p_g1->adhesion_loc.IN(delloc.y-delloc.z , delloc.z-delloc.x , delloc.x-delloc.y);
// 	p_g1->adhesion_loc = (p_g1->adhesion_loc / p_g1->adhesion_loc.norm()) * 0.5 * sigma_LJ;
// 	p_g1->adhesion_loc_inite = p_g1->adhesion_loc;
// 	p_g1->quat.IN(0.0 ,0.0 ,0.0);
// 	p_g1->quat_pre.IN(0.0 ,0.0 ,0.0);
// 	_quat<double> quat((hanshu+1.0)*PI/hanshu , delloc);
// 	quat/=quat.norm();
// 	p_g2->adhesion_loc = _quaternion::rot( quat , p_g1->adhesion_loc);
// 	p_g2->adhesion_loc_inite = p_g2->adhesion_loc;
// 	p_g2->quat.IN(0.0 ,0.0 ,0.0);
// 	p_g2->quat_pre.IN(0.0 ,0.0 ,0.0);
// 	
// 	return;
// }
// void polymerization(){
// 	for( int j=1 ; numbering[j].num ; j++){
// 		memset(&list_sub , 0 , sizeof(_list_sub));
// 		_ptcl* pa = numbering[j].ptcl_n[numbering[j].num-1];
// 		if(pa->bind!=-1)continue;
// 		_list* li = pa->li;
// 		_list* li_sub = 0;
// 		_ptcl* pb=0;
// 		for(int in=0;in<27;in++){
// 			if(li_sub=li->list_index[in]){
// 				for(int k=0 ; pb=li_sub->ptcl_index[k] ; k++){
// 					if(pb->cnf->feature == 0 && pb->p_level==1 ){
// 						_vec<double> delloc = pb->loc - pa->loc;
// #ifndef BROWNIAN
// //						_vec<double> delvel = pb->vel - pa->vel;
// #endif
// 						distance_boundary(delloc);
// 						if(delloc.norm() < r_por_f_g 
// #ifndef BROWNIAN
// //						&& delvel.norm()*delvel.norm() < r1 
// #endif
// 						){
// 							_vec<double> delloc2;
// 							if(numbering[j].num!=1){
// 								delloc2 = pa->loc - numbering[j].ptcl_n[numbering[j].num-2]->loc;
// 							}
// 							else {
// 								delloc2 = pa->loc - numbering[j].end_arp2_3->loc;
// 							}
// 							distance_boundary(delloc2);
// 							double delangl = acos( (delloc/delloc.norm()) * (delloc2/delloc2.norm())) + 1.0E-9;
// 							if(delangl < por_angl){
// 								double probability = probability_porimerize_f_p;
// 								double c = RAND;
// 								if(c<probability){
// 									list_sub.num +=1;
// 									list_sub.index[list_sub.num]=pb;
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		if(list_sub.num==1){
// 			polymerization_fila( pa , list_sub.index[1] );
// 			counter_pol++;
// 		}
// 		else if(list_sub.num>1){
// 			int num = (int)ceil(list_sub.num*RAND);
// 			while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
// 			polymerization_fila( pa , list_sub.index[num] );
// 			counter_pol++;
// 		}
// 	}
// 	if(numbering[0].num>1){
// 		int check_list=0;
// 		_ptcl* list_up[100][2]={0};
// 		_ptcl* pa = 0;
// 		_ptcl* pb = 0;
// 		_list* lia=0;
// 		_list* lib=0;
// 		for(int i=0;i<list_size;i++){
// 			lia=list+i;
// 			for(int j=0;pa=lia->ptcl_index[j];j++){
// 				if(pa->cnf->feature==0 && pa->p_level==1){
// 					for(int k=j+1;pb=lia->ptcl_index[k];k++){
// 						if(pb->cnf->feature==0 && pb->p_level==1){
// 							_vec<double> delloc = pa->loc - pb->loc;
// #ifndef BROWNIAN
// //							_vec<double> delvel = pa->vel - pb->vel;
// #endif
// 							distance_boundary(delloc);
// 							
// 							if(delloc.norm() < r_por_g_g 
// #ifndef BROWNIAN
// //							&& delvel.sqr() < r1 
// #endif
// 							){
// 								double probability = probability_porimerize_p_p;
// 								double c = RAND;
// 								
// 								if(c<probability){
// 									list_up[check_list][0]=pa;
// 									list_up[check_list][1]=pb;
// 									check_list++;
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 			for(int in=0;in<13;in++){
// 				if(lib=lia->list_index[in]){
// 					for(int j=0;pa=lia->ptcl_index[j];j++){
// 						if(pa->cnf->feature==0 && pa->p_level==1){
// 							for(int k=0;pb=lib->ptcl_index[k];k++){
// 								if(pb->cnf->feature==0 && pb->p_level==1){
// 									_vec<double> delloc = pa->loc - pb->loc;
// #ifndef BROWNIAN
// //									_vec<double> delvel = pa->vel - pb->vel;
// #endif
// 									distance_boundary(delloc);
// 									
// 									if(delloc.norm() < r_por_g_g 
// #ifndef BROWNIAN
// //									&& delvel.sqr() < r1 
// #endif
// 									){
// 										double probability = probability_porimerize_p_p;
// 										double c = RAND;
// 										
// 										if(c<probability){
// 											list_up[check_list][0]=pa;
// 											list_up[check_list][1]=pb;
// 											check_list++;
// 										}
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		if(check_list==1){
// 			polymerization_mono( list_up[0][0] , list_up[0][1] );
// 			counter_pol++;
// 		}
// 		else if(check_list!=0){
// 			_ptcl* pc=0;
// 			_ptcl* pd=0;
// 			for(int ia=0;ia<check_list;ia++){
// 				int check_kaburi_num=0;
// 				int check_kaburi[50]={0};
// 				for(int ib=0;ib<2;ib++){
// 					pc=list_up[ia][ib];
// 					for(int ja=ia+1;ja<check_list;ja++){
// 						for(int jb=0;jb<2;jb++){
// 							pd=list_up[ja][jb];
// 							if(pc==pd){
// 								check_kaburi_num++;
// 								check_kaburi[0]=ia;
// 								check_kaburi[check_kaburi_num]=ja;
// 							}
// 						}
// 					}
// 				}
// 				if(!check_kaburi_num){
// 					if(list_up[ia][0]->cnf->feature==0 &&list_up[ia][1]->cnf->feature==0 ){
// 						polymerization_mono( list_up[ia][0] , list_up[ia][1] );
// 						counter_pol++;
// 					}
// 				}
// 				else{
// 					int index=(int)std::floor((check_kaburi_num+1)*RAND);
// 					while(index==(check_kaburi_num+1))index=(int)std::floor((check_kaburi_num+1)*RAND);
// 					if(list_up[index][0]->cnf->feature==0 &&list_up[index][1]->cnf->feature==0 ){
// 						int index_sub=check_kaburi[index];
// 						polymerization_mono( list_up[index_sub][0] , list_up[index_sub][1] );
// 						counter_pol++;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	/*
// 	if(numbering[0].num>1){
// 		_ptcl* pa = 0;
// 		for( int j = 0; pa=numbering[0].ptcl_n[j] ; j++){
// 			memset(&list_sub , 0 , sizeof(_list_sub));
// 			if(pa->p_level==1){
// 				_list* li=pa->li;
// 				_list* li_sub=0;
// 				_ptcl* pb = 0;
// 				for(int in=0;in<14;in++){
// 					if(li_sub=li->list_index[in]){
// 						for(int k=0 ; pb=li_sub->ptcl_index[k] ; k++){
// 							if(pb->cnf->feature==0 && pb->p_level==1 && pb!=pa){
// 								_vec<double> delloc = pa->loc - pb->loc;
// 								_vec<double> delvel = pa->vel - pb->vel;
// 								distance_boundary(delloc);
// 								
// 								if(delloc.norm() < r_por_g_g && delvel.norm()*delvel.norm() < r1 ){
// 									double probability = probability_porimerize_p_p;
// 									double c = RAND;
// 									
// 									if(c<probability){
// 										list_sub.num +=1;
// 										list_sub.index[list_sub.num]= pb;
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 			if(list_sub.num==1){
// 				polymerization_mono( pa , list_sub.index[1] );
// 				counter_pol++;
// 			}
// 			else if(list_sub.num>1){
// 				int num = (int)ceil(list_sub.num*RAND);
// 				while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
// 				polymerization_mono( pa , list_sub.index[num] );
// 				counter_pol++;
// 			}
// 		}
// 	}*/
// 	return;
// }
// 
// /*void polymerization(){
// 	for( int j=1 ; numbering[j].num ; j++){
// 		memset(&list_sub , 0 , sizeof(_list_sub));
// 		_ptcl* pa = numbering[j].ptcl_n[numbering[j].num-1];
// 		if(pa->bind!=-1)continue;
// -;
// 		_ptcl* pb=0;
// 		for(int k=1 ; pb=li->index[k] ; k++){
// //			_ptcl* pb = li->index[k];
// 			if(pb->cnf->feature == 0 && pb->p_level==1 ){
// 				_vec<double> delloc = pb->loc - pa->loc;
// 				_vec<double> delvel = pb->vel - pa->vel;
// 				distance_boundary(delloc);
// 				if(delloc.norm() < r_por_f_g && delvel.norm()*delvel.norm() < r1 ){
// 					_vec<double> delloc2 = pa->loc - numbering[j].ptcl_n[numbering[j].num-2]->loc;
// 					distance_boundary(delloc2);
// 					double delangl = acos( (delloc/delloc.norm()) * (delloc2/delloc2.norm())) + 1.0E-9;
// 					if(delangl < por_angl){
// 						double probability = probability_porimerize_f_p;
// 						double c = RAND;
// 						if(c<probability){
// 							list_sub.num +=1;
// 							list_sub.index[list_sub.num]=pb;
// 						}
// 					}
// 				}
// 			}
// 		}
// 		if(list_sub.num==1){
// 			polymerization_fila( pa , list_sub.index[1] );
// 			counter_pol++;
// 		}
// 		else if(list_sub.num>1){
// 			int num = (int)ceil(list_sub.num*RAND);
// 			while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
// 			polymerization_fila( pa , list_sub.index[num] );
// 			counter_pol++;
// 		}
// 	}
// 	
// 	if(numbering[0].num>1){
// 		_ptcl* pa = 0;
// 		for( int j = 0; pa=numbering[0].ptcl_n[j] ; j++){
// 			memset(&list_sub , 0 , sizeof(_list_sub));
// //			_ptcl* pa = numbering[0].ptcl_n[j];
// 			if(pa->p_level==1){
// 				_list* li=pa->li;
// 				_ptcl* pb = 0;
// 				for(int k=1 ; pb=li->index[k] ; k++){
// //					_ptcl* pb = li->index[k];
// 					if(pb->cnf->feature==0 && pb->p_level==1 && pb - pa>0){
// 						_vec<double> delloc = pa->loc - pb->loc;
// 						_vec<double> delvel = pa->vel - pb->vel;
// 						distance_boundary(delloc);
// 						
// 						if(delloc.norm() < r_por_g_g && delvel.norm()*delvel.norm() < r1 ){
// 							double probability = probability_porimerize_p_p;
// 							double c = RAND;
// 							
// 							if(c<probability){
// 								list_sub.num +=1;
// 								list_sub.index[list_sub.num]= pb;
// 							}
// 						}
// 					}
// 				}
// 			}
// 			if(list_sub.num==1){
// 				polymerization_mono( pa , list_sub.index[1] );
// 				counter_pol++;
// 			}
// 			else if(list_sub.num>1){
// 				int num = (int)ceil(list_sub.num*RAND);
// 				while(num<=0 || num>list_sub.num){num = (int)ceil(list_sub.num*RAND);}
// 				polymerization_mono( pa , list_sub.index[num] );
// 				counter_pol++;
// 			}
// 		}
// 	}
// 	return;
// }*/
// #endif //POLYMERIZATION
// #ifdef SEVERING
// void severing_sub(int feature,int num){
// 	if(numbering[feature].num==2){
// 		if(num!=0){cout<<"severing_error"<<endl;}
// 		for(int i=0 ; i<2 ; i++){
// 			numbering[feature].ptcl_n[i]->cnf->feature = 0;
// 			numbering[feature].ptcl_n[i]->time = t;
// 			numbering[feature].ptcl_n[i]->num = numbering[0].num + i;
// 			numbering[0].ptcl_n[numbering[0].num + i]=numbering[feature].ptcl_n[i];
// 			numbering[feature].ptcl_n[i]=0;
// 		}
// 		numbering[feature].num=0;
// 		numbering[0].num+=2;
// 		for(int i=feature+1 ; numbering[i].num ; i++){
// 			numbering[i-1]=numbering[i];
// 			memset(&numbering[i],0,sizeof(_numbering));
// 			for(int j=0 ; numbering[i-1].ptcl_n[j] ; j++){
// 				numbering[i-1].ptcl_n[j]->cnf->feature = i-1;
// 			}
// 		}
// 		no_fila--;
// 	}
// 	else if(num == 0){
// 		_ptcl* p = numbering[feature].ptcl_n[0];
// 		p->cnf->feature=0;
// 		p->time=t;
// 		p->num = numbering[0].num;
// 		numbering[0].ptcl_n[numbering[0].num]=p;
// 		numbering[0].num++;
// 		numbering[feature].sort_n(0);
// 		numbering[feature].num--;
// 	}
// 	else if(num==numbering[feature].num - 2){
// 		_ptcl* p = numbering[feature].ptcl_n[num+1];
// 		p->cnf->feature=0;
// 		p->time=t;
// 		p->num = numbering[0].num;
// 		numbering[0].ptcl_n[numbering[0].num]=p;
// 		numbering[0].num++;
// 		numbering[feature].num--;
// 		numbering[feature].ptcl_n[num+1]=0;
// 	}
// 	else{
// 		no_fila++;
// 		for(int i=num+1 ; numbering[feature].ptcl_n[i] ; i++){
// 			_ptcl* p = numbering[feature].ptcl_n[i];
// 			p->cnf->feature=no_fila;
// 			p->num=i-num-1;
// 			numbering[feature].ptcl_n[i]=0;
// 			numbering[feature].num--;
// 			numbering[no_fila].ptcl_n[numbering[no_fila].num]=p;
// 			numbering[no_fila].num++;
// 		}
// 	}
// 	return;
// }
// inline double sev_pro(_ptcl* pa , _ptcl* pb){
// 	//_vec<double> delloc = pa->loc - pb->loc;
// 	//distance_boundary(delloc);
// 	//return pro_sev*exp(-1.0*(delloc.norm()/l_eq-1)/e_c);
// 	return pro_sev;
// }
// void severing(){
// 	for(int i=1 ; numbering[i].num ; i++){
// 		if(numbering[i].num!=1){
// 			memset(&list_sub_sev , 0 , sizeof(_list_sub_sev));
// 			for(int j=1 ; numbering[i].ptcl_n[j] ; j++){
// 				double probability=sev_pro(numbering[i].ptcl_n[j-1] , numbering[i].ptcl_n[j]);
// 				double c = RAND;
// 				if(c<probability){
// 					list_sub_sev.num++;
// 					if(list_sub_sev.num>=500){
// 	//					fout_err_sub("list_sub_max_error\n");
// 						fout_err<<"list_sub_max_error"<<endl;
// 						list_sub_sev.num -= 1;
// 						break;
// 					}
// 					list_sub_sev.index[list_sub_sev.num]=j-1;
// 				}
// 				
// 			}
// 			if(list_sub_sev.num==1){
// 				int j = list_sub_sev.index[1];
// 				severing_sub(i,j);
// 				fout_err<<"severing "
// 					<<i<<' '<<j<<' '
// 					<<endl;
// 				counter_severing++;
// 			}
// 			else if(list_sub_sev.num>1){
// 				int num = (int)ceil(list_sub_sev.num*RAND);
// 				while(num<=0 || num>list_sub_sev.num){num = (int)ceil(list_sub_sev.num*RAND);}
// 				int j = list_sub_sev.index[num];
// 				severing_sub(i,j);
// 				fout_err<<"severing "
// 					<<i<<' '<<j<<' '
// 					<<endl;
// 				counter_severing++;
// 			}
// 		}
// 	}
// 	return;
// }
// #endif //SEVERING

// void fout_err_sub(string ss1 , int end_fout){
// 	static string ss_fout_err_sub[10];
// 	static unsigned int counter_four_err=0;
// 	ss_fout_err_sub[counter_four_err] = ss1;
// 	counter_four_err++;
// 	if(counter_four_err==10 || end_fout!=0){
// 		ofstream fout_err;
// 		stringstream ss_fout;
// 		ss_fout.str("");
// 		ss_fout<<OUTPUT_N<<"fout_err.dat";
// 		fout_err.open(ss_fout.str().c_str(),ios::out);
// 		for(int i=0 ; i<counter_four_err ; i++){
// 			fout_err<<ss_fout_err_sub[i];
// 		}
// 		counter_four_err=0;
// 		memset(ss_fout_err_sub,0,10*sizeof(string));
// 	}
// 	return;
// }


//void replenish_ptcl(const int myrank,const int numprocs){
////	static int Active_Ptcl=Number_Ptcl_Active+Number_Ptcl_Spare_Dens;
//	
//	if(Number_Ptcl-Active_Ptcl<0)return;
//	
//		for(_ptcl* p=ptcl;p!=ptcl+Active_Ptcl;p++){
//			_list* li=p->list_calc();
//			if(li!=p->cnf->li){
//				p->cnf->li->rm(p,p->list_num);
//				p->add(li);
//				li->add(p);
//			}
//		}
//		int counter_ptcl(0);
//		for(int x=0;x<list_X;x++){
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<offset_impound;z++){
//			counter_ptcl+=(list+(z + y*list_Z + x*(list_Y*list_Z)))->num;
//		}}}
//		int margin=Number_Ptcl_Spare_Dens-counter_ptcl;
//		
//#ifdef USE_MPI
//		mympi::bit_sender.pack(&Active_Ptcl);
//		mympi::bit_sender.pack(&margin);
//#endif //#ifdef USE_MPI
//		
//		if(margin>0){
//			for(int i=Active_Ptcl;(i<Active_Ptcl+margin)&&(i<Number_Ptcl);i++){
//				_vec<double>loc;
//				{
//					loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
//					while( !(ptcl_ini_g_actin((ptcl+i),loc,0,0)) ){
//						loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
//					}
//				}
//#ifdef USE_MPI
//				mympi::bit_sender.pack(&loc);
//#endif //#ifdef USE_MPI
//			}
//			Active_Ptcl+=margin;
//			if(Active_Ptcl > Number_Ptcl)Active_Ptcl=Number_Ptcl;
//		}
//		fout_err<<"step = "<<step<<"\t Active_Ptcl = "<<Active_Ptcl<<endl;
//	
//	return;
//}

double LJ_Potential(const _vec<double>& vec,_ptcl* p2){
	static const double s06=(double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ
						  * (double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ;
	_vec<double> dloc=vec - p2->loc;
	distance_boundary(dloc);
	double r2=dloc.sqr();
	if(r2>=r_cut2){
		return 0.0;
	}
	double r06=s06/(r2*r2*r2);
	double r12=r06*r06;
	return 4.0*epsilon_LJ*(r12 - r06) + epsilon_LJ;
}

bool Chack_Potensial(const _vec<double>& vec){
	double potential(0.0);
	_list* li(list_calc(vec));
	_ptcl* pa(NULL);
	
	_list* li_sub(NULL);
	for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
		if(li_sub=li->list_index[a]){
#else
		li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
			for(int k=0;pa=li_sub->ptcl_index[k];k++){		//Lennard-Jones potential actin to actin in enighbercell
				potential += LJ_Potential(vec,pa);
			}                          
#ifdef WALL_BOUNDARY
		}
#endif //WALL_BOUNDARY
	}
//	cout<<potential<<endl;
	return (RAND<exp(-potential/SYS_kT));
}

void debug_()
{
	_ptcl* p;
	std::ofstream fout("debug0526.dat",ios::out);
		int counter_ptcl(0);
		for(int x=0;x<list_X;x++){//粒子数の確認
		for(int y=0;y<list_Y;y++){
		for(int z=0;z<offset_impound;z++){
			counter_ptcl+=(list+(z + y*list_Z + x*(list_Y*list_Z)))->num;
		}}}
	fout<<counter_ptcl<<std::endl;
	for(_ptcl* p(&ptcl[0]);p!=&ptcl[Number_Ptcl];++p){
		fout<<p->index<<'\t'<<p->loc<<std::endl;
	}
	fout<<std::endl;
	for(int x=0;x<list_X;x++){//粒子数の確認
	for(int y=0;y<list_Y;y++){
	for(int z=0;z<list_Z;z++){
		fout<<x<<'\t'<<y<<'\t'<<z<<'\t'<<(list+z+y*list_Z+x*list_Z*list_Y)->num<<'\t';
		for(int i(0);p=(list+z+y*list_Z+x*list_Z*list_Y)->ptcl_index[i];++i){fout<<p->index<<'\t';}fout<<std::endl;
	}}}
	
	
	
	
	
	fout.close();
}
void replenish_ptcl(const int myrank,const int numprocs){
	static const double Z_check(offset_impound*SYS_Z/(double)list_Z);
//	cout<<"replinish "<<step<<" "<<Active_Ptcl<<std::endl;
//	cout<<step<<'\t'<<(ptcl+4614)->loc<<endl;
#ifdef USE_MPI
	for(_ptcl* p=&ptcl[0];p!=&ptcl[Active_Ptcl];p++){//位置はギャザーされているがリスト更新はされていない．
		_list* li=p->list_calc();
		if(li!=p->cnf->li){
			p->cnf->li->rm(p,p->cnf->list_num);
			p->add(li);
			li->add(p);
		}
	}
#endif //#ifdef USE_MPI
	int counter_ptcl(0);
	for(int x=0;x<list_X;x++){//粒子数の確認
	for(int y=0;y<list_Y;y++){
	for(int z=0;z<offset_impound;z++){
		counter_ptcl+=(list+(z + y*list_Z + x*(list_Y*list_Z)))->num;
	}}}
	int margin=Number_Ptcl_Spare_Dens-counter_ptcl;
	if(Active_Ptcl+margin > Number_Ptcl){
		margin=Number_Ptcl-Active_Ptcl;
	}
#ifdef USE_MPI
	replenish_changelist.clear();
#endif//#ifdef USE_MPI
	
#ifdef USE_MPI
	for(int i(0);i<Number_Ptcl_Spare_Dens;++i){*(replenish_index+i)=-1;}
#endif//#ifdef USE_MPI
	
	for(int x=0;x<list_X;x++){
	for(int y=0;y<list_Y;y++){
	for(int z=0;z<offset_impound;z++){
		(list+(z + y*list_Z + x*(list_Y*list_Z)))->init();
	}}}
	
	int counter(0);
	for(_ptcl* p(&ptcl[0]);p!=&ptcl[Active_Ptcl];++p)
	{
		if(p->index < Active_Ptcl+margin){
			if(counter < Number_Ptcl_Spare_Dens){
				if(p->loc.z < Z_check){//基本的に入れ替え
//				if((p->list_index%list_Z) < offset_impound){//基本的に入れ替え
					if(!(p->cnf->feature)){
						_vec<double>& location(*(replenish_loc+counter));
						location.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
						while(!(Chack_Potensial(location))){
							location.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
						}
						p->loc.IN(location);
						p->loc_pre.IN(p->loc);
//								p->cnf->feature = 0;
						p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
						p->adhesion_loc_inite.IN(p->adhesion_loc);
						p->quat.IN(0.0,0.0,0.0);
						p->quat_pre.IN(p->quat);
						//list更新
						_list* li_sub=p->list_calc();
						p->add(li_sub);
						li_sub->add(p);
#ifdef USE_MPI
						*(replenish_index+counter)=p->index;
#endif//#ifdef USE_MPI
					}else{//フィラメントが魔酔いこんだとき
						_list* li_sub=p->list_calc();
						p->add(li_sub);
						li_sub->add(p);
						
					}
					++counter;
				}
			}else{
				if(p->loc.z < Z_check){
//				if((p->list_index%list_Z) < offset_impound){//基本的に入れ替え//後から探して入れ替え//配列の整列
					int tmp(-1);
					for(_ptcl* q(&ptcl[Active_Ptcl+margin]);q!=&ptcl[Active_Ptcl];q++){
						if(q->loc.z >= Z_check || q->cnf->feature!=0){
							tmp=q->index;
							break;
						}
					}
					if(!(tmp+1)){
					cerr<<"repli err margin "<<margin<<" index "<<p->index<<"AM "<<Active_Ptcl+margin<<"cou "<<counter<<" "<<endl;
						debug_();
						exit(1);//とりあえず強制終了
					}
					_list* li=ptcl[tmp].cnf->li;
					li->rm(&ptcl[tmp],ptcl[tmp].cnf->list_num);
// cout<<"copy "<<tmp<<std::endl;
//					p->Copy(&ptcl[tmp]);
					ptcl.Copy(*p,ptcl[tmp]);
// cout<<"copy "<<std::endl;
					ptcl[tmp].loc.z=0;
					ptcl[tmp].cnf->feature=0;
					
					p->add(li);
					li->add(p);
					
#ifdef USE_MPI
					replenish_changelist.push_back(p->index);
					replenish_changelist.push_back(tmp);
#endif//USE_MPI
					++counter;
				}//else{
				//	
				//}
			}
		}else{// if(counter >= Number_Ptcl_Spare_Dens){//粒子削除
//			cout<<"rm "<<p->index<<'\t'<<p->loc.z<<" "<<p->list_index%list_Z<<endl;
//			if(p->loc.z >= Z_check){
//				cout<<"miss_take "<<step<<endl;
//				exit(0);
//			}
		}
	}
	if(counter < Number_Ptcl_Spare_Dens){//粒子追加
		_ptcl* p(&ptcl[Active_Ptcl]);
		while(counter < Number_Ptcl_Spare_Dens){
			_vec<double>& location(*(replenish_loc+counter));
			location.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
			while(!(Chack_Potensial(location))){
				location.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
			}
			p->loc.IN(location);
			p->loc_pre.IN(p->loc);
			p->cnf->feature = 0;
			p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
			p->adhesion_loc_inite.IN(p->adhesion_loc);
			p->quat.IN(0.0,0.0,0.0);
			p->quat_pre.IN(p->quat);
			//list更新
			_list* li_sub=p->list_calc();
			p->add(li_sub);
			li_sub->add(p);
			
			p->cnf->bind=-1;
			p->cnf->plus=NULL;
			p->cnf->minus=NULL;
			p->cnf->edge_arp2_3=NULL;
			p->cnf->end_arp2_3=NULL;
			
//					cout<<"add  "<<p->index<<'\t'<<p->loc.z<<" "<<p->list_index%list_Z<<endl;
			
#ifdef USE_MPI
			*(replenish_index+counter)=p->index;
#endif//#ifdef USE_MPI
			++counter;
			++p;
		}//std::cout<<p->index<<std::endl;
	}
	
	Active_Ptcl+=margin;
//	std::cout<<Active_Ptcl<<'\t'<<margin<<std::endl;
//	cout<<step<<'\t'<<(ptcl+4614)->loc<<endl;
//	fout_err<<"replinish "<<step<<" "<<Active_Ptcl<<std::endl;
	
	return;
}


//double LJ_Potential(_ptcl* p1,_ptcl* p2){
//	static const double s06=(double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ
//						  * (double)sigma_LJ * (double)sigma_LJ * (double)sigma_LJ;
//	_vec<double> dloc=p1->loc - p2->loc;
//	distance_boundary(dloc);
//	double r2=dloc.sqr();
//	if(r2>=r_cut2){
//		return 0.0;
//	}
//	double r06=s06/(r2*r2*r2);
//	double r12=r06*r06;
//	return 4.0*epsilon_LJ*(r12 - r06) + epsilon_LJ;
//}
//
//
//bool Calc_Potensial(void){
//	double potential(0.0);
//	for(int x=0;x<list_X;x++){
//	for(int y=0;y<list_Y;y++){
//	for(int z=0;z<offset_impound;z++){
//		_list* li(list+ (z + y*list_Z + x*(list_Y*list_Z)) );
//		_ptcl* pa(NULL);
//		_ptcl* pb(NULL);
//		for(int j=0;pa=li->ptcl_index[j];j++){				//Lennard-Jones potential actin to actin in mycell
//			for(int k=j+1;pb=li->ptcl_index[k];k++){		//Lennard-Jones potential actin to actin
//				potential += LJ_Potential(pa,pb);
////if(LJ_Potential(pa,pb)>1.0)cout<<LJ_Potential(pa,pb)<<" "<<(pa->loc-pb->loc).norm()<<"1 "<<pa->index<<'\t'<<pa->loc<<'\n'<<"2 "<<pb->index<<'\t'<<pb->loc<<endl;
////				if(RAND>exp(-LJ_Potential(pa,pb)/SYS_kT))return false;
//			}
//		}
//		
//		_list* li_sub = 0;
//		for(int a=0;a<13;a++){
//#ifdef WALL_BOUNDARY
//			if(li_sub=li->list_index[a]){
//#else
//			li_sub=li->list_index[a];
//#endif //WALL_BOUNDARY
//			for(int j=0;pa=li->ptcl_index[j];j++){
//				for(int k=0;pb=li_sub->ptcl_index[k];k++){		//Lennard-Jones potential actin to actin in enighbercell
//					potential += LJ_Potential(pa,pb);
////if(LJ_Potential(pa,pb)>1.0)cout<<LJ_Potential(pa,pb)<<" "<<(pa->loc-pb->loc).norm()<<"1 "<<pa->index<<'\t'<<pa->loc<<'\n'<<"2 "<<pb->index<<'\t'<<pb->loc<<endl;
////					if(RAND>exp(-LJ_Potential(pa,pb)/SYS_kT))return false;
//				}                          
//			}
//#ifdef WALL_BOUNDARY
//			}
//#endif //WALL_BOUNDARY
//		}
//	}}}
//	cout<<potential<<endl;
//	return (RAND<exp(-potential/SYS_kT));
//	return true;
//}
//void replenish_ptcl(const int myrank,const int numprocs){
//	static const double Z_check(offset_impound*SYS_Z/(double)list_Z);
//// 	cout<<"replinish "<<step<<" "<<Active_Ptcl<<std::endl;
//	
//#ifdef USE_MPI
//		for(_ptcl* p=ptcl;p!=ptcl+Active_Ptcl;p++){//何しているの？//位置はギャザーされているがリスト更新はされていない．
//			_list* li=p->list_calc();
//			if(li!=p->cnf->li){
//				p->cnf->li->rm(p,p->list_num);
//				p->add(li);
//				li->add(p);
//			}
//		}
//#endif //#ifdef USE_MPI
//		int counter_ptcl(0);
//		for(int x=0;x<list_X;x++){//粒子数の確認
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<offset_impound;z++){
//			counter_ptcl+=(list+(z + y*list_Z + x*(list_Y*list_Z)))->num;
//		}}}
//		int margin=Number_Ptcl_Spare_Dens-counter_ptcl;
//		if(Active_Ptcl+margin > Number_Ptcl){
//			margin=Number_Ptcl-Active_Ptcl;
//		}
//		bool flag(false);
//#ifdef USE_MPI
//			replenish_changelist.clear();
//#endif//#ifdef USE_MPI
//		while(!(flag))
//		{
//			
//			for(int i(0);i<Number_Ptcl_Spare_Dens;++i)
//			{
//				(replenish_loc+i)->IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(offset_impound*SYS_Z/(double)list_Z-1.0)*RAND+1.0);
//#ifdef USE_MPI
//				*(replenish_index+i)=-1;
//#endif//#ifdef USE_MPI
//			}
//			for(int x=0;x<list_X;x++){//粒子数の確認
//			for(int y=0;y<list_Y;y++){
//			for(int z=0;z<offset_impound;z++){
//				(list+(z + y*list_Z + x*(list_Y*list_Z)))->init();
//			}}}
//			
//			int counter(0);
//			for(_ptcl* p(ptcl);p!=ptcl+Active_Ptcl;++p)
//			{
//				if(p->index<Active_Ptcl+margin){
//					if(counter < Number_Ptcl_Spare_Dens){
//						if(p->loc.z < Z_check){//基本的に入れ替え
//							if(!(p->cnf->feature)){
//								p->loc.IN(*(replenish_loc+counter));
//								p->loc_pre.IN(p->loc);
////								p->cnf->feature = 0;
//								p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
//								p->adhesion_loc_inite.IN(p->adhesion_loc);
//								p->quat.IN(0.0,0.0,0.0);
//								p->quat_pre.IN(p->quat);
//								//list更新
//								_list* li_sub=p->list_calc();
//								p->add(li_sub);
//								li_sub->add(p);
//								
//#ifdef USE_MPI
//								*(replenish_index+counter)=p->index;
//#endif//#ifdef USE_MPI
//								++counter;
//							}else{//フィラメントが魔酔いこんだとき
//								_list* li_sub=p->list_calc();
//								p->add(li_sub);
//								li_sub->add(p);
//								
//								++counter;
//							}
//						}
//					}else{//後から探して入れ替え
//						if(p->loc.z < Z_check){
//							int tmp(0);
//							for(_ptcl* q(ptcl+Active_Ptcl);q!=ptcl+Active_Ptcl;q++){
//								if(q->cnf->feature){
//									tmp=q->index;
//									break;
//								}
//							}
//							if(!tmp)exit(1);//とりあえず強制終了
//							_list* li=(ptcl+tmp)->cnf->li;
//							li->rm((ptcl+tmp),(ptcl+tmp)->list_num);
//cout<<"copy "<<tmp<<std::endl;
//							p->Copy(ptcl+tmp);
//cout<<"copy "<<std::endl;
//							(ptcl+tmp)->cnf->feature=0;
//							p->add(li);
//							li->add(p);
//							
//#ifdef USE_MPI
//							replenish_changelist.pash_back(p->index);
//							replenish_changelist.pash_back(tmp);
//#endif//USE_MPI
//						}else{//フィラメントが魔酔いこんだとき
//							_list* li_sub=p->list_calc();
//							p->add(li_sub);
//							li_sub->add(p);
//						}
//					}
//				}else{// if(counter >= Number_Ptcl_Spare_Dens){//粒子削除
//					
//				}
//			}
//			if(counter < Number_Ptcl_Spare_Dens){//粒子追加
//				_ptcl* p(ptcl + Active_Ptcl);
//				while(counter<Number_Ptcl_Spare_Dens){
//					p->loc.IN(*(replenish_loc+counter));
//					p->loc_pre.IN(p->loc);
//					p->cnf->feature = 0;
//					p->adhesion_loc.IN(0.5 * sigma_LJ , 0.0 , 0.0);
//					p->adhesion_loc_inite.IN(p->adhesion_loc);
//					p->quat.IN(0.0,0.0,0.0);
//					p->quat_pre.IN(p->quat);
//					//list更新
//					_list* li_sub=p->list_calc();
//					p->add(li_sub);
//					li_sub->add(p);
//					
//					p->bind=-1;
//					p->plus_in(-1);
//					p->minus_in(-1);
//					p->edge_arp2_3_in(-1);
//					p->end_arp2_3_in(-1);
//					
//#ifdef USE_MPI
//					*(replenish_index+counter)=p->index;
//#endif//#ifdef USE_MPI
//					++counter;
//					++p;
//				}
//			}
//			flag=Calc_Potensial();
//		}
//		
//		Active_Ptcl+=margin;
//	cout<<"replinish "<<step<<" "<<Active_Ptcl<<std::endl;
//	
//	return;
//}
//
//
