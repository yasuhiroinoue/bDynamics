/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/

#include"_ctrl.h"
#ifdef ARP2_3

#include<iostream>
#include<fstream>
#include<cmath>
#include"_arp2_3.h"
#include"_periodic.h"
#include"_variable.h"
#include "quaternion.h"
#include"_mympi.h"
#include "_function.h"

using namespace std;

void list_calc(_vec<double>&,int&,int&,int&);
//_list* list_calc(const _vec<double>& loc);

namespace _ARP2_3{
	_arp2_3* arp2_3=NULL;
// 	_Bind* bind_arp2_3=NULL;
	 double ene_pot=0.0;
	bool Arp_ini_sub(_arp2_3* p, const _vec<double>& loc){
		bool flag=true;
		_list* li_sub=list_calc(loc);
		_list* li_sub2=li_sub;
		_ptcl* pa;
		_arp2_3* par;
		for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
			if(li_sub2=li_sub->list_index[a]){
#else
			li_sub2=li_sub->list_index[a];
#endif //WALL_BOUNDARY
			for(int j=0;pa=li_sub2->ptcl_index[j];j++){
				_vec<double>delloc = pa->loc - loc;
				distance_boundary(delloc);
				if(delloc.norm()<SIGMA_LJ_ACTIN_ARP)flag=false;
			}
			for(int j=0;par=li_sub2->arp_index[j];j++){
				_vec<double>delloc = par->loc - loc;
				distance_boundary(delloc);
				if(delloc.norm()<SIGMA_LJ_ARP_ARP)flag=false;
			}
#ifdef WALL_BOUNDARY
			}
#endif //WALL_BOUNDARY
		}
		if(flag){
			p->loc=loc;
			p->loc_pre=p->loc;
			p->index=p - arp2_3;
			p->add(li_sub);
			li_sub->add(p);
		}
		return flag;
	}
	void arp_init(void){
		_ARP2_3::arp2_3=new _ARP2_3::_arp2_3[Number_Arp2_3];
		
		for(int i=0;i<Number_Arp2_3;i++){(_ARP2_3::arp2_3 + i)->index=i;}
	}
	
	void arp_posion_init(void){
//		_ARP2_3::arp2_3=new _ARP2_3::_arp2_3[Number_Arp2_3];
		_vec<double> loc;
		for(int i=0;i<Number_Arp2_3;i++){
#ifndef WALL_BOUNDARY
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
// 				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)SYS_Z-2.0*SIGMA_LJ_ACTIN_ARP)*RAND+SIGMA_LJ_ACTIN_ARP);
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)WALL_Z_INIT-1.0*SIGMA_LJ_ACTIN_ARP)*RAND + offset_impound*SYS_Z/(double)list_Z);
#endif
			while(!(Arp_ini_sub((arp2_3+i),loc)) ){
#ifndef WALL_BOUNDARY
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,(double)SYS_Z*RAND);
#else
// 				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)SYS_Z-2.0*SIGMA_LJ_ACTIN_ARP)*RAND+SIGMA_LJ_ACTIN_ARP);
				loc.IN((double)SYS_X*RAND,(double)SYS_Y*RAND,((double)WALL_Z_INIT-1.0*SIGMA_LJ_ACTIN_ARP)*RAND + offset_impound*SYS_Z/(double)list_Z);
#endif
			}
		}
// 		arp2_3[0].loc.IN(16.0,16.0-SIGMA_LJ_ACTIN_ARP,16.0);arp2_3[0].loc_pre=arp2_3[0].loc;
		return;
	}
	void arp_finalize(){
		delete[] _ARP2_3::arp2_3;
	}
	
	namespace _ARP_ARP{
		_vec<double> LJ(_arp2_3* p1,_arp2_3* p2,double* potential){
			static const double s06=SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP
								  * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP;
			_vec<double> F=_vec<double>(0.0,0.0,0.0);
			_vec<double> dloc=p1->loc - p2->loc;
			distance_boundary(dloc);
			double r2=dloc.sqr();
			if(r2>=Cut_Off_Arp2_3_Arp2_3_2){
				*potential=0.0;
				return F;
			}
				if( r2<CutOffUperArp2_3Arp2_3_2){
					r2=CutOffUperArp2_3Arp2_3_2;
					double r06=s06/(CutOffUperArp2_3Arp2_3_2*CutOffUperArp2_3Arp2_3_2*CutOffUperArp2_3Arp2_3_2);
					double r12=r06*r06;
					*potential=4.0*EPSILON_LJ_ARP_ARP*(r12 - r06) * (CutOffUperArp2_3Arp2_3 - dloc.norm() +SIGMA_LJ_ARP_ARP) + EPSILON_LJ_ARP_ARP;
					F=4.0*EPSILON_LJ_ARP_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
					return F;
				}
			double r06=s06/(r2*r2*r2);
			double r12=r06*r06;
			{
				*potential=4.0*EPSILON_LJ_ARP_ARP*(r12 - r06) + EPSILON_LJ_ARP_ARP;
				F=4.0*EPSILON_LJ_ARP_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
				return F;
			}
		}

		void Frc(const int myrank){
//			for(int x=myrank*calc_length_X;x<(myrank+1)*calc_length_X;x++){
////			for(int x=0;x<list_X;x++){
//			for(int y=0;y<list_Y;y++){
//			for(int z=0;z<list_Z;z++){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				int list_index =z + y*list_Z + x*(list_Y*list_Z);
				
				_list* li = list+list_index;
				_arp2_3* pa = 0;
				_arp2_3* pb = 0;
				if(li->num_arp){
					for(int j=0;pa=li->arp_index[j];j++){
						for(int k=j+1;pb=li->arp_index[k];k++){
							double potential;
							_vec<double> force_sub = LJ(pa,pb,&potential);
							pa->force += force_sub;
							pb->force -= force_sub;
							ene_pot += potential;
						}
					}
					
					_list* li_sub = NULL;
// 					bool flag=false;
// 					for(int li_z_sub=-1;li_z_sub<=1;li_z_sub++){
// 					for(int li_y_sub=-1;li_y_sub<=1;li_y_sub++){
// 					for(int li_x_sub=-1;li_x_sub<=1;li_x_sub++){
// 						if(li_z_sub==0&&li_y_sub==0&&li_x_sub==0){flag=true;break;}
// 						int lx = x + li_x_sub;
// 						int ly = y + li_y_sub;
// 						int lz = z + li_z_sub;
// 						while(lx<0)			lx+=list_X;
// 						while(lx>=list_X)	lx-=list_X;
// 						while(ly<0)			ly+=list_Y;
// 						while(ly>=list_Y)	ly-=list_Y;
// #ifdef WALL_BOUNDARY
// 						if(lz<0 || lz>=list_Z)continue;
// #else //WALL_BOUNDARY
// 						while(lz<0)			lz+=list_Z;
// 						while(lz>=list_Z)	lz-=list_Z;
// #endif //WALL_BOUNDARY
// 						li_sub=(list+ (lx + ly*list_X + lz*(list_X*list_Y) ));
// 						for(int j=0;pa=li->arp_index[j];j++){
// 							for(int k=0;pb=li_sub->arp_index[k];k++){
// 								double potential;
// 								_vec<double> force_sub = LJ(pa,pb,&potential);
// 								pa->force += force_sub;
// 								pb->force -= force_sub;
// //								ene_pot += repot(dloc);
// 								ene_pot += potential;
// 							}
// 						}
// 					}if(flag)break;
// 					}if(flag)break;
// 					}
// 					
					for(int a=0;a<13;a++){
#ifdef WALL_BOUNDARY
						if(li_sub=li->list_index[a]){
#else
						li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
						for(int j=0;pa=li->arp_index[j];j++){
							for(int k=0;pb=li_sub->arp_index[k];k++){
								double potential;
								_vec<double> force_sub = LJ(pa,pb,&potential);
								pa->force += force_sub;
								pb->force -= force_sub;
								ene_pot += potential;
							}
						}
#ifdef WALL_BOUNDARY
						}
#endif //WALL_BOUNDARY
					}
				}
			}}}
			return;
			
		}
	}	//namespace _ARP_ARP
	
	
	namespace _ACTIN_ARP{
		_vec<double> LJ(_ptcl* p1,_arp2_3* p2,double* potential){
			static const double s06=SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP
								  * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP;
			_vec<double> F=_vec<double>(0.0,0.0,0.0);
			_vec<double> dloc=p1->loc - p2->loc;
			distance_boundary(dloc);
			double r2=dloc.sqr();
			if(r2>=Cut_Off_Actin_Arp2_3_2){
				*potential=0.0;
				return F;
			}
				if( r2<CutOffUperActinArp2_3_2){
					r2=CutOffUperActinArp2_3_2;
					double r06=s06/(CutOffUperActinArp2_3_2*CutOffUperActinArp2_3_2*CutOffUperActinArp2_3_2);
					double r12=r06*r06;
					*potential=4.0*EPSILON_LJ_ACTIN_ARP*(r12 - r06) * (CutOffUperActinArp2_3 - dloc.norm() + SIGMA_LJ_ACTIN_ARP) + EPSILON_LJ_ACTIN_ARP;
					F=4.0*EPSILON_LJ_ACTIN_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
					return F;
				}
			double r06=s06/(r2*r2*r2);
			double r12=r06*r06;
			if(r2>=Cut_Off_Actin_Arp2_3_2){
				*potential=0.0;
				return F;
			}else{
				*potential=4.0*EPSILON_LJ_ACTIN_ARP*(r12 - r06) + EPSILON_LJ_ACTIN_ARP;
				F=4.0*EPSILON_LJ_ACTIN_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
				return F;
			}
		}
		void Bind_frc(_arp2_3* arp){
			_ptcl* pa=NULL;
			_ptcl* pb=NULL;
			if(pa=arp->edge_ptcl()){
				_vec<double> dloc1 = (pa->loc + pa->adhesion_loc) - arp->loc;
				distance_boundary(dloc1);
				double norm1 = dloc1.norm();
				_vec<double> udloc1 = dloc1/norm1;
				{//親フィラメントと結合・張力
					_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(norm1-EQ_LEN_ACTIN_ARP)*udloc1;
					pa->adhesion_frc += force_sub1;
					arp->force -= force_sub1;
					ene_pot+=K_TENSILE_ACTIN_ARP*(norm1-EQ_LEN_ACTIN_ARP)*(norm1-EQ_LEN_ACTIN_ARP)/2.0;
				}{//曲げ
					double norm2 = pa->adhesion_loc.norm();
					_vec<double> udloc2 = pa->adhesion_loc/norm2;
					double cosphi = udloc1 * udloc2;
					//double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
					
					_vec<double> force_sub1 = K_ANGLE_ACTIN_ARP*(udloc1 - cosphi * udloc2)/(norm2);
					_vec<double> force_sub2 = K_ANGLE_ACTIN_ARP*(udloc2 - cosphi * udloc1)/(norm1);
					//dif_pot=K*sin(phi-phi0)=k*sin(phi-pi)=-k*sinphi
					pa->force += force_sub1;
					pa->adhesion_frc -= force_sub1 + force_sub2;
					arp->force += force_sub2;
					ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
				}
			}
			if(pb=arp->end_ptcl()){
				_vec<double> dloc1 = arp->loc - pb->loc;
				distance_boundary(dloc1);
				double norm1 = dloc1.norm();
				_vec<double> udloc1 = dloc1/norm1;
				{//子フィラメントと結合
					_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(norm1-SIGMA_LJ_ACTIN_ARP)*udloc1;
					arp->force += force_sub1;
					pb->force -= force_sub1;
					ene_pot+=K_TENSILE_ACTIN_ARP*(udloc1-SIGMA_LJ_ACTIN_ARP)*(udloc1-SIGMA_LJ_ACTIN_ARP)/2.0;
				}{//子フィラメントの角度保持・70度
					if(pa=arp->edge_ptcl()){
						_vec<double> dloc2 = arp->loc - (pa->loc + pa->adhesion_loc);
						distance_boundary(dloc2);
						double norm2 = dloc2.norm();
						_vec<double> udloc2 = dloc2/norm2;
						double cosphi = (udloc1) * (udloc2);
						double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
						
//							cout<<std::acos(cosphi)*180.0/PI<<'\t'<<std::asin(sinphi)*180.0/PI<<endl;
						
						static const double cosphi0 = std::cos(16.0/18.0*PI);
						static const double sinphi0 = std::sin(16.0/18.0*PI);
						
						double F=K_ANGLE_ACTIN_ARP*(cosphi0 - cosphi*sinphi0/sinphi);
						
						_vec<double> force_sub1 = F*(udloc1 - cosphi * udloc2)/(norm2);
						_vec<double> force_sub2 = F*(udloc2 - cosphi * udloc1)/(norm1);
						
						pa->adhesion_frc -= force_sub2;
						arp->force += force_sub1 + force_sub2;
						pb->force -= force_sub1;
						ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
					}
				}{//同一平面に
					if(pa=arp->edge_ptcl()){
						_vec<double> dloc_sub=-1.0*dloc1;
						_vec<double> r23 =  arp->loc - pa->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
						_vec<double> r21;
						if(pa->cnf->plus){
							r21 = pa->cnf->plus->loc - pa->loc;
							distance_boundary(r21);
						}else if(pa->cnf->minus){
							r21 = pa->loc - pa->cnf->minus->loc;
							distance_boundary(r21);
						}else if(pa->cnf->end_arp2_3){
							r21 = pa->loc - pa->cnf->end_arp2_3->loc;
							distance_boundary(r21);
						}else {
							fout_err<<"arp_err"<<endl;
//							continue;
						}
						_vec<double> n123 = r21 % r23;
						_vec<double> n234 = r23 % dloc_sub;
						double norm_n1 = n123.norm();
						double norm_n2 = n234.norm();
						n123/=norm_n1;
						n234/=norm_n2;
						
						double cosph = n123*n234;
//							double sinph = -1.0*(n123%n234)*ur23;
//							if(sinph>=0.0)sinph+=1.0e-9;
//							else sinph-=1.0e-9;
						
	//					double cosph0 = std::cos(PI);
	//					double sinph0 = std::sin(PI);
						
					//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
						double F = K_Torsion_ARP;
						
						_vec<double> N1 = (n234 - cosph*n123)/norm_n1;
						_vec<double> N2 = (n123 - cosph*n234)/norm_n2;
						
						_vec<double> force_sub1 = F * (N1 % r23);
						_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % dloc_sub);
						_vec<double> force_sub3 = F * (N2 % r23);
						
						if(pa->cnf->plus){
							pa->cnf->plus->force += force_sub1;
							pa->force += force_sub2 - force_sub1;
						}else if(pa->cnf->minus){
							pa->cnf->minus->force -= force_sub1;
							pa->force += force_sub2 + force_sub1;
						}else if(pa->cnf->end_arp2_3){
							pa->cnf->end_arp2_3->force -= force_sub1;
							pa->force += force_sub2 + force_sub1;
						}else {
							fout_err<<"arp_err"<<endl;
//							continue;
						}
						arp->force += force_sub3 - force_sub2;
						pb->force -= force_sub3;
						
					}
				}{//子フィラメントの結合点を外に向ける
					double norm2 = pb->adhesion_loc.norm();
					_vec<double> udloc2 = pb->adhesion_loc/norm2;
					double cosphi = udloc1 * udloc2;
					double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
					
					double F = -1.0*K_TENSILE_ACTIN_ARP*(cosphi/sinphi);
					
					_vec<double> force_sub1 = F*(udloc1 - cosphi * udloc2)/norm2;
					_vec<double> force_sub2 = F*(udloc2 - cosphi * udloc1)/norm1;
					
					pb->adhesion_frc += force_sub1;
					pb->force -= force_sub1 + force_sub2;
					arp->force += force_sub2;
					ene_pot+=1.0-sinphi;
				}{//子フィラメントの回転を止める
					if(pa=arp->edge_ptcl()){
// 						_vec<double> r23 = pb->loc - arp->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
						_vec<double> r23 = dloc1; _vec<double> ur23 = udloc1;
						_vec<double> r21 = pa->loc - arp->loc;	distance_boundary(r21);
						_vec<double> n123 = r21 % r23;
						_vec<double> n234 = r23 % pb->adhesion_loc;
						double norm_n1 = n123.norm();
						double norm_n2 = n234.norm();
						n123/=norm_n1;
						n234/=norm_n2;
						
						double cosph = n123*n234;
//						double sinph = -1.0*(n123%n234)*ur23;
//						if(sinph>=0.0)sinph+=1.0e-9;
//						else sinph-=1.0e-9;
						
//						double cosph0 = std::cos(PI);
//						double sinph0 = std::sin(PI);
						
					//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
						double F = K_Torsion_ARP;
						
						_vec<double> N1 = (n234 - cosph*n123)/norm_n1;
						_vec<double> N2 = (n123 - cosph*n234)/norm_n2;
						
						_vec<double> force_sub1 = F * (N1 % r23);
						_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % pb->adhesion_loc);
						_vec<double> force_sub3 = F * (N2 % r23);
						pa->force += force_sub1;
						arp->force += force_sub2 - force_sub1;
						pb->force += force_sub3 - force_sub2;
						pb->adhesion_frc -= force_sub3;
					}
				}
				{//子フィラメントの曲げ剛性
					_ptcl* pc=NULL;
					if(pb->cnf->plus){
						pc=pb->cnf->plus;
//					if(pc=numbering[pb->feature].ptcl_n[1]){
						_vec<double> dloc2 = arp->loc - pb->loc;
						_vec<double> dloc3 = pc->loc - pb->loc;
						distance_boundary(dloc2);
						distance_boundary(dloc3);
						
						double cosphi = (dloc2/dloc2.norm()) * (dloc3/dloc3.norm());
					//	double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
						
						_vec<double> force_sub1 = -k_angl*(dloc2/dloc2.norm() - cosphi * dloc3/dloc3.norm())/(dloc3.norm());
						_vec<double> force_sub2 = -k_angl*(dloc3/dloc3.norm() - cosphi * dloc2/dloc2.norm())/(dloc2.norm());
						
						pc->force += force_sub1;
						pb->force -= force_sub1 + force_sub2;
						arp->force += force_sub2;
						ene_pot+=k_angl*(1.0+cosphi);
					}
				}
			}
		}
		void Frc(const int myrank){
			_ptcl* pa=NULL;
			_arp2_3* arp=NULL;
//			for(int x=myrank*calc_length_X;x<(myrank+1)*calc_length_X;x++){
////			for(int x=0;x<list_X;x++){
//			for(int y=0;y<list_Y;y++){
//			for(int z=0;z<list_Z;z++){
			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
				_list* li = list+(z + y*list_Z + x*(list_Y*list_Z));
				_list* li_sub=NULL;
				if(li->num_arp){
					for(int j=0;arp=li->arp_index[j];j++){
						Bind_frc(arp);
					}
					for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
						if(li_sub=li->list_index[a]){
#else
						li_sub=li->list_index[a];
#endif //WALL_BOUNDARY
						for(int j=0;arp=li->arp_index[j];j++){
							for(int k=0;pa=li_sub->ptcl_index[k];k++){
								double potential;
								_vec<double> force_sub = LJ(pa,arp,&potential);
								pa->force += force_sub;
								arp->force -= force_sub;
								ene_pot += potential;
							}
						}
#ifdef WALL_BOUNDARY
						}
#endif //WALL_BOUNDARY
					}
				}
			}}}
			return;
		}
		
	}	//namespace _ACTIN_ARP
	
	//曲率接着しやすさ
//	void Calc_Strain(const int myrank,const int numprocs){
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
//			for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
//			for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
//			for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
//			_list* li(list+i);
//			_ptcl* pa(NULL);
//			_ptcl* pplus(NULL);
//			_ptcl* pminus(NULL);
//			for(int j=0;pa=li->ptcl_index[j];j++){
//				pplus=pa->plus;
//				pminus=pa->minus;
//				if(pplus!=NULL && pminus!=NULL){//フィラメントの間の粒子
//					_vec<double>dloc1= pa->loc - pminus->loc;
//					_vec<double>dloc2= pplus->loc - pa->loc;
//					distance_boundary(dloc1);
//					distance_boundary(dloc2);
//					_vec<double>kappa=(dloc2/dloc2.norm() - dloc1/dloc1.norm())/((dloc2 + dloc1).norm()*0.5);
//					pa->strain=-(kappa*pa->adhesion_loc);
//				}else if(pplus!=NULL || pminus!=NULL){//フィラメントの端
//					pa->strain=0.0;
//				}
//			}
//		}
//		return;
//	}
//	inline double arp_bind_pro(_ptcl *pa){
//		return pro_bind*exp(-sigma_LJ*sigma_LJ*0.125*(pa->strain*pa->strain+2.0*pa->strain));
//	}
	
	
	void Bind(const int myrank,const int numprocs){
//		Calc_Strain(myrank,numprocs);
		_ptcl* p=NULL;
		_arp2_3* arp=NULL;
		_list_sub list_sub;
		//memset(&list_sub,0,sizeof(_list_sub));
		//for(int m=0;m<Number_Arp2_3;m++){
// 			arp=(arp2_3+m);
//		for(int x=myrank*calc_length_X;x<(myrank+1)*calc_length_X;x++){
////		for(int x=0;x<list_X;x++){
//		for(int y=0;y<list_Y;y++){
//		for(int z=0;z<list_Z;z++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			_list* li_arp = list+(z + y*list_Z + x*(list_Y*list_Z));
			_list* li_sub1=NULL;
			for(int j=0;arp=li_arp->arp_index[j];j++){
				if(!arp->edge_ptcl()){
					memset(&list_sub,0,sizeof(_list_sub));
//					_list* li=NULL;
					for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
						if(li_sub1=li_arp->list_index[a]){
#else
						li_sub1=li_arp->list_index[a];
#endif //WALL_BOUNDARY
						for(int j=0;p=li_sub1->ptcl_index[j];j++){
							if((p->cnf->feature) && (!p->cnf->edge_arp2_3)){
								_vec<double> dloc = arp->loc - (p->loc + p->adhesion_loc);
								distance_boundary(dloc);
								if(dloc.norm()<=SIGMA_LJ_ARP_ARP*0.5
									&& (p->adhesion_loc/p->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
									
//									double c = RAND;
//									if(c<arp_bind_pro(p)){
										list_sub.index[list_sub.num]=p;
										list_sub.num++;
		//								cout<<list_sub.num<<endl;
									}
//								}
							}
						}
#ifdef WALL_BOUNDARY
						}
#endif //WALL_BOUNDARY
					}
					if(list_sub.num==1){
						_ptcl* p1=list_sub.index[0];
						arp->edge_ptcl_in(p1->index);
						p1->cnf->edge_arp2_3=arp;
#ifdef USE_MPI
						{
							mympi::_conformation_change cc(arp_edge,arp->index,p1->index);
							mympi::bit_sender.pack(&cc);
						}
#endif//USE_MPI
					}else if(list_sub.num>1){
						int num = (int)std::floor((double)list_sub.num * RAND);
						while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
						_ptcl* p1=list_sub.index[num];
						arp->edge_ptcl_in(p1->index);
						p1->cnf->edge_arp2_3 = arp;
#ifdef USE_MPI
						{
							mympi::_conformation_change cc(arp_edge,arp->index,p1->index);
							mympi::bit_sender.pack(&cc);
						}
#endif//USE_MPI
					}
				}
				else if(!arp->end_ptcl()){
					if(!arp->edge_ptcl()->cnf->feature){
						arp->edge_ptcl()->cnf->edge_arp2_3=NULL;//NULL;
						arp->edge_ptcl_in(-1);//=NULL;
						continue;
					}
					memset(&list_sub,0,sizeof(_list_sub));
					_vec<double> axis1;
					{//軸の計算
						_ptcl* pa=NULL;
						//_ptcl* pb=NULL;
						if(pa = arp->edge_ptcl()){
							
							_vec<double> dloc1 = arp->loc - pa->loc;
							_vec<double> dloc2;
							distance_boundary(dloc1);
							dloc1/=dloc1.norm();
// 							if(pa->num==0){
// 								if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 									dloc2 = pb->loc - pa->loc;
// 								}else if(numbering[pa->feature].end_arp2_3){
// 									dloc2 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// 								}else{fout_err<<"arp_err_bind"<<endl;	continue;}
// 							}else if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 								dloc2 = pb->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 							}else{
// 								dloc2 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 							}
// 							if(pa->b_plus){
// 								dloc2 = pa->b_plus->plus->loc - pa->loc;
// 							}else if(pa->b_minus){
// 								dloc2 =  pa->loc - pa->b_minus->minus->loc;
// 							}else if(pa->end_arp2_3){
// 								dloc2 = pa->loc - pa->end_arp2_3->loc;
// 							}
							if(pa->cnf->plus){
								if(pa->cnf->minus){
									dloc2 = pa->cnf->plus->loc - pa->cnf->minus->loc;
								}else{
									dloc2 = pa->cnf->plus->loc - pa->loc;
								}
							}else{
								if(pa->cnf->minus){
									dloc2 = pa->loc - pa->cnf->minus->loc;
								}else if(pa->cnf->end_arp2_3){
									dloc2 = pa->loc - pa->cnf->end_arp2_3->loc;
								}
							}
							
							distance_boundary(dloc2);
							dloc2/=dloc2.norm();
							
							double costhe = dloc2*dloc1;
							double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
							
							static const double costhe0 = std::cos(PI*70.0/180.0);
							static const double sinthe0 = std::sin(PI*70.0/180.0);
							
							axis1 = (costhe0 - sinthe0*costhe/sinthe)*dloc2 + (sinthe0/sinthe)*dloc1;
						}else{fout_err<<"arp_err_bind"<<endl;	continue;}
					}
					_list* li=NULL;
					for(int a=0;a<27;a++){
#ifdef WALL_BOUNDARY
						if(li=arp->li()->list_index[a]){
#else
						li=arp->li()->list_index[a];
#endif //WALL_BOUNDARY
						for(int j=0;p=li->ptcl_index[j];j++){
							if(!(p->cnf->feature)){
								_vec<double> dloc = p->loc - arp->loc;
								distance_boundary(dloc);
								if(dloc.norm()<=SIGMA_LJ_ACTIN_ARP && (axis1)*(dloc/dloc.norm())>0.5){
									list_sub.index[list_sub.num]=p;
									list_sub.num++;
		//							cout<<list_sub.num<<endl;
									
								}
							}
						}
#ifdef WALL_BOUNDARY
						}
#endif //WALL_BOUNDARY
					}
					if(list_sub.num){
						p=list_sub.index[0];
						if(list_sub.num>1){
							int num = (int)std::floor((double)list_sub.num * RAND);
							while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
							p=list_sub.index[num];
						}
						arp->end_ptcl_in(p->index);
// 						numbering[0].sort_n(p->num);
// 						numbering[0].num--;
// 						no_fila++;
						p->cnf->feature = 1;//no_fila;
// 						p->num = 0;
// 						p->time = t;
						p->cnf->end_arp2_3 = (arp);
						plus_end->add(p);
						
						{
							_vec<double> dloc1 = (arp->edge_ptcl()->loc + arp->edge_ptcl()->adhesion_loc) - arp->loc;
							_vec<double> dloc2 = p->loc - arp->loc;
							distance_boundary(dloc1);
							distance_boundary(dloc2);
							dloc1/=dloc1.norm();
							dloc2/=dloc2.norm();
							
							double costhe = dloc2*dloc1;
							double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
							
							_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
							p->adhesion_loc = axis * sigma_LJ * 0.5;
							p->adhesion_loc_inite = p->adhesion_loc;
							p->quat.IN(0.0 ,0.0 ,0.0);
							p->quat_pre.IN(0.0 ,0.0 ,0.0);
							
						}
						
#ifdef USE_MPI
						{
							conf_ch cc=arp_end;
							mympi::bit_sender.pack(&cc);
							mympi::bit_sender.pack(&arp->index);
							mympi::bit_sender.pack(&p->index);
							mympi::bit_sender.pack(&p->adhesion_loc);
						}
#endif//USE_MPI
						
	// 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
	// 					distance_boundary(axis);
	// 					axis = axis % arp->adhesion_loc;
	// 					_quat<double> quat(0.5*PI,axis);
	// 					quat/=quat.norm();
	// 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
	// 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
	// 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
						
// 					}else if(list_sub.num>1){
// 						int num = (int)std::floor((double)list_sub.num * RAND);
// 						while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 						p=list_sub.index[num];
// 						arp->end_ptcl=p;
// 						
// 						numbering[0].sort_n(p->num);
// 						numbering[0].num--;
// 						no_fila++;
// 						p->feature = no_fila;
// 						p->num = 0;
// 						p->time = t;
// 						numbering[no_fila].ptcl_n[0]=p;
// 						numbering[no_fila].num=1;
// 						numbering[no_fila].end_arp2_3=arp;
// 						
// 						{
// 							_vec<double> dloc1 = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 							_vec<double> dloc2 = p->loc - arp->loc;
// 							distance_boundary(dloc1);
// 							distance_boundary(dloc2);
// 							dloc1/=dloc1.norm();
// 							dloc2/=dloc2.norm();
// 							
// 							double costhe = dloc2*dloc1;
// 							double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 							
// 							_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
// 							p->adhesion_loc = axis * sigma_LJ * 0.5;
// 							p->adhesion_loc_inite = p->adhesion_loc;
// 							p->quat.IN(0.0 ,0.0 ,0.0);
// 							p->quat_pre.IN(0.0 ,0.0 ,0.0);
// 						}
// 						
// 	// 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 	// 					distance_boundary(axis);
// 	// 					axis = axis % arp->adhesion_loc;
// 	// 					_quat<double> quat(0.5*PI,axis);
// 	// 					quat/=quat.norm();
// 	// 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// 	// 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// 	// 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 	// 					
					}
				}
				if((arp->edge_ptcl()) && (!arp->edge_ptcl()->cnf->feature)){
					arp->edge_ptcl()->cnf->edge_arp2_3=NULL;
//					numbering[arp->end_ptcl->feature].end_arp2_3=NULL;
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
#ifdef USE_MPI
					{
						mympi::_conformation_change cc(arp_debind,arp->index,0);
						mympi::bit_sender.pack(&cc);
					}
#endif//USE_MPI
				}
			}
		}}}
		return;
	}
	
// 	void Bind(){
// 		_ptcl* p=NULL;
// 		_arp2_3* arp=NULL;
// 		_list_sub list_sub;
// 		//memset(&list_sub,0,sizeof(_list_sub));
// 		for(int m=0;m<Number_Arp2_3;m++){
// 			arp=(arp2_3+m);
// 			if(!arp->edge_ptcl){
// 				memset(&list_sub,0,sizeof(_list_sub));
// 				_list* li=NULL;
// 				for(int a=0;a<27;a++){
// #ifdef WALL_BOUNDARY
// 					if(li=arp->li->list_index[a]){
// #else
// 					li=arp->li->list_index[a];
// #endif //WALL_BOUNDARY
// 					for(int j=0;p=li->ptcl_index[j];j++){
// 						if((p->feature) && (!p->edge_arp2_3)){
// 							_vec<double> dloc = arp->loc - (p->loc + p->adhesion_loc);
// 							distance_boundary(dloc);
// 							if(dloc.norm()<=SIGMA_LJ_ARP_ARP*0.5
// 								&& (p->adhesion_loc/p->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
// 								list_sub.index[list_sub.num]=p;
// 								list_sub.num++;
// //								cout<<list_sub.num<<endl;
// 							}
// 						}
// 					}
// #ifdef WALL_BOUNDARY
// 					}
// #endif //WALL_BOUNDARY
// 				}
// 				if(list_sub.num==1){
// 					_ptcl* p1=list_sub.index[0];
// 					arp->edge_ptcl = p1;
// 					p1->edge_arp2_3 = arp;
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					_ptcl* p1=list_sub.index[num];
// 					arp->edge_ptcl=p1;
// 					p1->edge_arp2_3 = arp;
// 				}
// 			}
// 			else if(!arp->end_ptcl){
// 				if(!arp->edge_ptcl->feature){
// 					arp->edge_ptcl->edge_arp2_3=NULL;
// 					arp->edge_ptcl=NULL;
// 					continue;
// 				}
// 				memset(&list_sub,0,sizeof(_list_sub));
// 				_vec<double> axis1;
// 				{
// 					_ptcl* pa=NULL;
// 					_ptcl* pb=NULL;
// 					if(pa = arp->edge_ptcl){
// 						
// 						_vec<double> dloc1 = arp->loc - pa->loc;
// 						_vec<double> dloc2;
// 						distance_boundary(dloc1);
// 						dloc1/=dloc1.norm();
// 						if(pa->num==0){
// 							if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 								dloc2 = pb->loc - pa->loc;
// 							}else if(numbering[pa->feature].end_arp2_3){
// 								dloc2 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// 							}else{fout_err<<"arp_err_bind"<<endl;	continue;}
// 						}else if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 							dloc2 = pb->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 						}else{
// 							dloc2 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 						}
// 						distance_boundary(dloc2);
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						static const double costhe0 = std::cos(PI*70.0/180.0);
// 						static const double sinthe0 = std::sin(PI*70.0/180.0);
// 						
// 						axis1 = (costhe0 - sinthe0*costhe/sinthe)*dloc2 + (sinthe0/sinthe)*dloc1;
// 					}else{fout_err<<"arp_err_bind"<<endl;	continue;}
// 				}
// 				_list* li=NULL;
// 				for(int a=0;a<27;a++){
// #ifdef WALL_BOUNDARY
// 					if(li=arp->li->list_index[a]){
// #else
// 					li=arp->li->list_index[a];
// #endif //WALL_BOUNDARY
// 					for(int j=0;p=li->ptcl_index[j];j++){
// 						if(!(p->feature)){
// 							_vec<double> dloc = p->loc - arp->loc;
// 							distance_boundary(dloc);
// 							if(dloc.norm()<=SIGMA_LJ_ACTIN_ARP && (axis1)*(dloc/dloc.norm())>0.5){
// 								list_sub.index[list_sub.num]=p;
// 								list_sub.num++;
// 	//							cout<<list_sub.num<<endl;
// 								
// 							}
// 						}
// 					}
// #ifdef WALL_BOUNDARY
// 					}
// #endif //WALL_BOUNDARY
// 				}
// // 				for(int i=0; i<Number_Ptcl; ++i){
// // 					p=(ptcl+i);
// // 					if(!p->feature){
// // 						_vec<double> dloc = p->loc - arp->loc;
// // 						distance_boundary(dloc);
// // 						if(dloc.norm()<=SIGMA_LJ_ACTIN_ARP && (axis1)*(dloc/dloc.norm())>0.5){
// // 							list_sub.index[list_sub.num]=p;
// // 							list_sub.num++;
// // //							cout<<list_sub.num<<endl;
// // 							
// // 						}
// // 					}
// // 				}
// 				if(list_sub.num==1){
// 					p=list_sub.index[0];
// 					arp->end_ptcl=p;
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					
// 					{
// 						_vec<double> dloc1 = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 						_vec<double> dloc2 = p->loc - arp->loc;
// 						distance_boundary(dloc1);
// 						distance_boundary(dloc2);
// 						dloc1/=dloc1.norm();
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
// 						p->adhesion_loc = axis * sigma_LJ * 0.5;
// 						p->adhesion_loc_inite = p->adhesion_loc;
// 						p->quat.IN(0.0 ,0.0 ,0.0);
// 						p->quat_pre.IN(0.0 ,0.0 ,0.0);
// 						
// 					}
// 					
// 					
// // 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// // 					distance_boundary(axis);
// // 					axis = axis % arp->adhesion_loc;
// // 					_quat<double> quat(0.5*PI,axis);
// // 					quat/=quat.norm();
// // 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// // 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// // 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					p=list_sub.index[num];
// 					arp->end_ptcl=p;
// 					
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					{
// 						_vec<double> dloc1 = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 						_vec<double> dloc2 = p->loc - arp->loc;
// 						distance_boundary(dloc1);
// 						distance_boundary(dloc2);
// 						dloc1/=dloc1.norm();
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
// 						p->adhesion_loc = axis * sigma_LJ * 0.5;
// 						p->adhesion_loc_inite = p->adhesion_loc;
// 						p->quat.IN(0.0 ,0.0 ,0.0);
// 						p->quat_pre.IN(0.0 ,0.0 ,0.0);
// 					}
// 					
// // 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// // 					distance_boundary(axis);
// // 					axis = axis % arp->adhesion_loc;
// // 					_quat<double> quat(0.5*PI,axis);
// // 					quat/=quat.norm();
// // 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// // 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// // 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// // 					
// 				}
// 				
// 			}
// 			else if(!arp->edge_ptcl->feature){
// 				arp->edge_ptcl->edge_arp2_3=NULL;
// 				numbering[arp->end_ptcl->feature].end_arp2_3=NULL;
// 				arp->edge_ptcl=NULL;
// 				arp->end_ptcl=NULL;
// 			}
// 		}
// 		return;
// 	}
	
	void Arp_Location_FirstStep(const int myrank,const int numprocs){
//	cout<<"first"<<endl;
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li_sub = list+i;
			_arp2_3* p =NULL;
			for(int j=0;p=li_sub->arp_index[j];j++){
				p->loc+= (0.5*PTCL_dt* p->force /Arp2_3_Friction);
				p->force.IN(0.0,0.0,0.0);
				Periodic(p->loc);
			}
		}}}
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li_sub = list+i;
			_arp2_3* p =NULL;
			for(int j=0;p=li_sub->arp_index[j];j++){
				_list* li=p->list_calc();
				if(li!=p->li()){
#ifdef USE_MPI
					if(li->myrank!=-1){
if(NX>1){
						if((li->myrank/9)==2 && (p->li()->myrank/9)!=2){		send_list_arp_x_bck_inner_new.push_back(p);
						}else if((li->myrank/9)==0 && (p->li()->myrank/9)!=0){	send_list_arp_x_pre_inner_new.push_back(p);}
}
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(p);}
}
					}
#endif
					p->li()->rm(p,p->list_num);
					p->add(li);
					li->add(p);
					j--;
				}
			}
		}}}
		return;
	}
	void Arp_Location_SecondStep(const int myrank,const int numprocs){
//	cout<<"second"<<endl;
		static const double amplitude=sqrt(2.0*SYS_kT*Arp2_3_Friction/PTCL_dt);
		static const double root3=sqrt(3.0)*2.0;
		_vec<double> RandamF(0.0,0.0,0.0);
//		for(int i=myrank*calc_length_X*list_Y*list_Z;i<(myrank+1)*calc_length_X*list_Y*list_Z;i++){
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li_sub = list+i;
			_arp2_3* p =NULL;
			for(int j=0;p=li_sub->arp_index[j];j++){
				RandamF= (amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3);
				p->loc_pre+=PTCL_dt * (p->force + RandamF) / Arp2_3_Friction;
				p->force.IN(0.0,0.0,0.0);
				Periodic(p->loc_pre);
				p->loc=p->loc_pre;
			}
		}}}
		for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
		for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
		for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
			int i(z + y*list_Z + x*(list_Y*list_Z));
			_list* li_sub = list+i;
			_arp2_3* p =NULL;
			for(int j=0;p=li_sub->arp_index[j];j++){
				_list* li=p->list_calc();
				if(li!=p->li()){
#ifdef USE_MPI
					if(li->myrank!=-1){
if(NX>1){
						if((li->myrank/9)==2 && (p->li()->myrank/9)!=2){		send_list_arp_x_bck_inner_new.push_back(p);
						}else if((li->myrank/9)==0 && (p->li()->myrank/9)!=0){	send_list_arp_x_pre_inner_new.push_back(p);}
}
if(NY>1){
						if(((li->myrank/3)%3)==2 && ((p->li()->myrank/3)%3)!=2){		send_list_arp_y_bck_inner_new.push_back(p);
						}else if(((li->myrank/3)%3)==0 && ((p->li()->myrank/3)%3)!=0){	send_list_arp_y_pre_inner_new.push_back(p);}
}
if(NZ>1){
						if((li->myrank%3)==2 && (p->li()->myrank%3)!=2){		send_list_arp_z_bck_inner_new.push_back(p);
						}else if((li->myrank%3)==0 && (p->li()->myrank%3)!=0){	send_list_arp_z_pre_inner_new.push_back(p);}
}
					}
#endif
					p->li()->rm(p,p->list_num);
					p->add(li);
					li->add(p);
					j--;
				}
			}
		}}}
		return;
	}

}

//#endif //ARP2_3

#if 0
// template <class T> inline ostream& operator<< (ostream& stream,_quat<T> quat){
// 	return stream<<quat.n<<'\t'<<quat.x<<'\t'<<quat.y<<'\t'<<quat.z;
// }
// 
// namespace _ARP2_3{
// 	_arp2_3* arp2_3=NULL;
// 	 double ene_pot=0.0;
// 	bool Arp_ini_sub(_arp2_3* p, const _vec<double>& loc){
// 		bool flag=true;
// 		_list* li_sub=list_calc(loc);
// 		_list* li_sub2=li_sub;
// // 		int lx = (int)std::floor(loc.x/r_cut);
// // 		int ly = (int)std::floor(loc.y/r_cut);
// // 		int lz = (int)std::floor(loc.z/r_cut);
// // 		if(lx==list_X){lx -= 1;}
// // 		if(ly==list_Y){ly -= 1;}
// // 		if(lz==list_Z){lz -= 1;}
// // 		int li = lx + ly*list_X + lz*(list_X*list_Y);
// // 		_list* li_sub=list+li;
// // 		_list* li_sub2=list+li;
// 		
// 		_ptcl* pa;
// 		_arp2_3* par;
// 		for(int a=0;a<27;a++){
// #ifdef WALL_BOUNDARY
// 			if(li_sub2=li_sub->list_index[a]){
// #else
// 			li_sub2=li_sub->list_index[a];
// #endif //WALL_BOUNDARY
// 			for(int j=0;pa=li_sub2->ptcl_index[j];j++){
// 				_vec<double>delloc = pa->loc - loc;
// 				distance_boundary(delloc);
// 				if(delloc.norm()<SIGMA_LJ_ACTIN_ARP)flag=false;
// 			}
// 			for(int j=0;par=li_sub2->arp_index[j];j++){
// 				_vec<double>delloc = par->loc - loc;
// 				distance_boundary(delloc);
// 				if(delloc.norm()<SIGMA_LJ_ARP_ARP)flag=false;
// 			}
// #ifdef WALL_BOUNDARY
// 			}
// #endif //WALL_BOUNDARY
// 		}
// // 		for(int i=0;i<Number_Arp2_3;i++){
// // 			if(p==(arp2_3+i))break;
// // 			_vec<double>delloc = arp2_3[i].loc - loc;
// // 			distance_boundary(delloc);
// // 			if(delloc.norm()<SIGMA_LJ_ARP_ARP)flag=false;
// // 		}
// 		if(flag){
// 			p->loc=loc;
// 			p->loc_pre=p->loc;
// 			p->index=p - arp2_3;
// 			p->add(li_sub);
// 			li_sub->add(p);
// 		}
// 		return flag;
// 	}
// 	void arp_init(void){
// 		_ARP2_3::arp2_3=new _ARP2_3::_arp2_3[Number_Arp2_3];
// 		for(int i=0;i<Number_Arp2_3;i++){
// 			_vec<double> loc((double)SYS_X*RAND, (double)SYS_Y*RAND, (double)SYS_Z*RAND);
// 			while(!(Arp_ini_sub((arp2_3+i),loc)) ){loc.IN((double)SYS_X*RAND, (double)SYS_Y*RAND, (double)SYS_Z*RAND);}
// 		}
// // 		arp2_3[0].loc.IN(16.0,16.0-SIGMA_LJ_ACTIN_ARP,16.0);arp2_3[0].loc_pre=arp2_3[0].loc;
// 		return;
// 	}
// 	namespace _ARP_ARP{
// // 		inline double pot(const double r){
// // 			return (4.0*EPSILON_LJ_ARP_ARP)*pow(SIGMA_LJ_ARP_ARP/r,12.0);
// // 		}
// // 		inline double dif_pot( double r){
// // 			return -1.0*(4.0*EPSILON_LJ_ARP_ARP)*12.0/SIGMA_LJ_ARP_ARP*pow(SIGMA_LJ_ARP_ARP/r,13.0);
// // 		}
// // 		inline double d_dif_pot( double r){
// // 			return (4.0*EPSILON_LJ_ARP_ARP)*156.0/(SIGMA_LJ_ARP_ARP*SIGMA_LJ_ARP_ARP)*pow(SIGMA_LJ_ARP_ARP/r,14.0);
// // 		}
// // 		
// // 		inline double repot(_vec<double> r){
// // 			if(r.norm()<=Cut_Off_Arp2_3_Arp2_3){
// // 				return pot(r.norm())-pot((double)Cut_Off_Arp2_3_Arp2_3)
// // 					+dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)
// // 					+0.5*d_dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3);
// // 			}
// // 			else{return 0.0;}
// // 		}
// // 		inline _vec<double> refrce(_vec<double> r){//LJポテンシャル
// // 			if(r.norm()<=Cut_Off_Arp2_3_Arp2_3){
// // 				return -( dif_pot(r.norm())
// // 					+dif_pot((double)Cut_Off_Arp2_3_Arp2_3)
// // 					+d_dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)
// // 					)*r/r.norm();
// // 			}
// // 			else{return _vec<double>(0.0,0.0,0.0);}
// // 		}
// 		inline _vec<double> LJ(_arp2_3* p1,_arp2_3* p2,double* potential){
// 			static const double s06=SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP
// 								  * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP * SIGMA_LJ_ARP_ARP;
// 			_vec<double> dloc=p1->loc - p2->loc;
// 			distance_boundary(dloc);
// 			double r2=dloc.sqr();
// 			double r06=s06/(r2*r2*r2);
// 			double r12=r06*r06;
// 			_vec<double> F=_vec<double>(0.0,0.0,0.0);
// 			if(r2>=Cut_Off_Arp2_3_Arp2_3_2){
// 				*potential=0.0;
// 				return F;
// 			}else{
// 				*potential=4.0*EPSILON_LJ_ARP_ARP*(r12 - r06) + EPSILON_LJ_ARP_ARP;
// 				F=4.0*EPSILON_LJ_ARP_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
// 				return F;
// 			}
// 		}
// 
// 		void Frc(){
// // 			for(int i=0;i<Number_Arp2_3;i++){
// // 				for(int j=i+1;j<Number_Arp2_3;j++){
// // 					_vec<double> dloc = _ARP2_3::arp2_3[i].loc - _ARP2_3::arp2_3[j].loc;
// // 					distance_boundary(dloc);
// // 					_vec<double> force_sub = refrce(dloc);
// // 					_ARP2_3::arp2_3[i].force += force_sub;
// // 					_ARP2_3::arp2_3[j].force -= force_sub;
// // 					ene_pot += repot(dloc);
// // 				}
// // 			}
// //			for(int i=0;i<list_size;i++){
// 			for(int z=0;z<list_Z;z++){
// 			for(int y=0;y<list_Y;y++){
// 			for(int x=0;x<list_X;x++){
// 				int list_index =x + y*list_X + z*(list_X*list_Y);
// 				
// 				_list* li = list+list_index;
// 				_arp2_3* pa = 0;
// 				_arp2_3* pb = 0;
// 				if(li->num_arp){
// 					for(int j=0;pa=li->arp_index[j];j++){
// 						for(int k=j+1;pb=li->arp_index[k];k++){
// //							_vec<double> dloc = pa->loc - pb->loc;
// 							//distance_boundary(dloc);
// //							_vec<double> force_sub = refrce(dloc);
// 							double potential;
// 							_vec<double> force_sub = LJ(pa,pb,&potential);
// 							pa->force += force_sub;
// 							pb->force -= force_sub;
// //							ene_pot += repot(dloc);
// 							ene_pot += potential;
// 						}
// 					}
// 					
// 					_list* li_sub = NULL;
// 					bool flag=false;
// 					for(int li_z_sub=-2;li_z_sub<=2;li_z_sub++){
// 					for(int li_y_sub=-2;li_y_sub<=2;li_y_sub++){
// 					for(int li_x_sub=-2;li_x_sub<=2;li_x_sub++){
// 						if(li_z_sub==0&&li_y_sub==0&&li_x_sub==0){flag=true;break;}
// 						int lx = x + li_x_sub;
// 						int ly = y + li_y_sub;
// 						int lz = z + li_z_sub;
// 						while(lx<0)			lx+=list_X;
// 						while(lx>=list_X)	lx-=list_X;
// 						while(ly<0)			ly+=list_Y;
// 						while(ly>=list_Y)	ly-=list_Y;
// #ifdef WALL_BOUNDARY
// 						if(lz<0 || lz>=list_Z)continue;
// #else //WALL_BOUNDARY
// 						while(lz<0)			lz+=list_Z;
// 						while(lz>=list_Z)	lz-=list_Z;
// #endif //WALL_BOUNDARY
// 						li_sub=(list+ (lx + ly*list_X + lz*(list_X*list_Y) ));
// 						for(int j=0;pa=li->arp_index[j];j++){
// 							for(int k=0;pb=li_sub->arp_index[k];k++){
// // 								_vec<double> dloc = pa->loc - pb->loc;
// // 								distance_boundary(dloc);
// // 								_vec<double> force_sub = refrce(dloc);
// 								double potential;
// 								_vec<double> force_sub = LJ(pa,pb,&potential);
// 								pa->force += force_sub;
// 								pb->force -= force_sub;
// //								ene_pot += repot(dloc);
// 								ene_pot += potential;
// 							}
// 						}
// 					}if(flag)break;
// 					}if(flag)break;
// 					}
// 					
// 	// 				for(int a=0;a<13;a++){
// 	// #ifdef WALL_BOUNDARY
// 	// 					if(li_sub=li->list_index[a]){
// 	// #else
// 	// 					li_sub=li->list_index[a];
// 	// #endif //WALL_BOUNDARY
// 	// 					for(int j=0;pa=li->arp_index[j];j++){
// 	// 						for(int k=0;pb=li_sub->arp_index[k];k++){
// 	// // 							if((pa->feature!=0) 
// 	// // 								&& (pa->feature==pb->feature) 
// 	// // 		//						&& (pa->num-pb->num==1 || pa->num-pb->num==-1)
// 	// // 								){continue;}
// 	// 							_vec<double> dloc = pa->loc - pb->loc;
// 	// 							distance_boundary(dloc);
// 	// 							_vec<double> force_sub = refrce(dloc);
// 	// 							pa->force += force_sub;
// 	// 							pb->force -= force_sub;
// 	// 							ene_pot += repot(dloc);
// 	// 						}
// 	// 					}
// 	// #ifdef WALL_BOUNDARY
// 	// 					}
// 	// #endif //WALL_BOUNDARY
// 	// 				}
// 	//			}
// 				}}}
// 			}
// 			return;
// 			
// 		}
// 	}	//namespace _ARP_ARP
// 	
// 	
// 	namespace _ACTIN_ARP{
// // 		inline double pot(double r){//ポテンシャル(スカラー)
// // 			return (4.0*EPSILON_LJ_ACTIN_ARP)*pow(SIGMA_LJ_ACTIN_ARP/r,12.0);
// // 		}
// // 		inline double dif_pot(double r){
// // 			return -1.0*(4.0*EPSILON_LJ_ACTIN_ARP)*12.0/SIGMA_LJ_ACTIN_ARP*pow(SIGMA_LJ_ACTIN_ARP/r,13.0);
// // 		}
// // 		inline double d_dif_pot(double r){
// // 			return (4.0*EPSILON_LJ_ACTIN_ARP)*156.0/(SIGMA_LJ_ACTIN_ARP*SIGMA_LJ_ACTIN_ARP)*pow(SIGMA_LJ_ACTIN_ARP/r,14.0);
// // 		}
// // 
// // 		inline double repot(_vec<double> r){
// // 			if(r.norm()<=Cut_Off_Actin_Arp2_3){
// // 				return pot(r.norm())
// // 					-pot(Cut_Off_Actin_Arp2_3)
// // 					+dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)
// // 					+0.5*d_dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3);
// // 			}
// // 			else{return 0.0;}
// // 		}
// // 		inline _vec<double> refrce( _vec<double> r){//LJポテンシャル
// // 			if(r.norm()<=Cut_Off_Actin_Arp2_3){
// // 				return -( dif_pot(r.norm())
// // 					+dif_pot(Cut_Off_Actin_Arp2_3)
// // 					+d_dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)
// // 					)*r/r.norm();
// // 			}
// // 			else{return _vec<double>(0.0,0.0,0.0);}
// // 		}
// 		inline _vec<double> LJ(_ptcl* p1,_arp2_3* p2,double* potential){
// 			static const double s06=SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP
// 								  * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP * SIGMA_LJ_ACTIN_ARP;
// 			_vec<double> dloc=p1->loc - p2->loc;
// 			distance_boundary(dloc);
// 			double r2=dloc.sqr();
// 			double r06=s06/(r2*r2*r2);
// 			double r12=r06*r06;
// 			_vec<double> F=_vec<double>(0.0,0.0,0.0);
// 			if(r2>=Cut_Off_Actin_Arp2_3_2){
// 				*potential=0.0;
// 				return F;
// 			}else{
// 				*potential=4.0*EPSILON_LJ_ACTIN_ARP*(r12 - r06) + EPSILON_LJ_ACTIN_ARP;
// 				F=4.0*EPSILON_LJ_ACTIN_ARP*(12.0*r12 - 6.0*r06)*dloc/r2;
// 				return F;
// 			}
// 		}
// 		void Frc(){
// 			_ptcl* pa=NULL;
// 			_ptcl* pb=NULL;
// 			_arp2_3* arp=NULL;
// 			for(int i=0;i<Number_Arp2_3;i++){
// 				arp=(arp2_3+i);
// 				
// 				if(pa=arp->edge_ptcl){
// 					_vec<double> dloc1 = (pa->loc + pa->adhesion_loc) - arp->loc;
// 					distance_boundary(dloc1);
// 					double norm1 = dloc1.norm();
// 					_vec<double> udloc1 = dloc1/norm1;
// 					{//親フィラメントと結合・張力
// 						_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(norm1-EQ_LEN_ACTIN_ARP)*udloc1;
// 						pa->adhesion_frc += force_sub1;
// 						arp->force -= force_sub1;
// 						ene_pot+=K_TENSILE_ACTIN_ARP*(norm1-EQ_LEN_ACTIN_ARP)*(norm1-EQ_LEN_ACTIN_ARP)/2.0;
// 					}{//曲げ
// 						double norm2 = pa->adhesion_loc.norm();
// 						_vec<double> udloc2 = pa->adhesion_loc/norm2;
// 						double cosphi = udloc1 * udloc2;
// 						//double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						_vec<double> force_sub1 = K_ANGLE_ACTIN_ARP*(udloc1 - cosphi * udloc2)/(norm2);
// 						_vec<double> force_sub2 = K_ANGLE_ACTIN_ARP*(udloc2 - cosphi * udloc1)/(norm1);
// 						//dif_pot=K*sin(phi-phi0)=k*sin(phi-pi)=-k*sinphi
// 						pa->force += force_sub1;
// 						pa->adhesion_frc -= force_sub1 + force_sub2;
// 						arp->force += force_sub2;
// 						ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
// 					}
//  				}
// 				if(pb=arp->end_ptcl){
// 					_vec<double> dloc1 = arp->loc - pb->loc;
// 					distance_boundary(dloc1);
// 					double norm1 = dloc1.norm();
// 					_vec<double> udloc1 = dloc1/norm1;
// 					{//子フィラメントと結合
// 						_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(norm1-SIGMA_LJ_ACTIN_ARP)*udloc1;
// 						arp->force += force_sub1;
// 						pb->force -= force_sub1;
// 						ene_pot+=K_TENSILE_ACTIN_ARP*(udloc1-SIGMA_LJ_ACTIN_ARP)*(udloc1-SIGMA_LJ_ACTIN_ARP)/2.0;
// 					}{//子フィラメントの角度保持・70度
// 						if(pa=arp->edge_ptcl){
// 							_vec<double> dloc2 = arp->loc - (pa->loc + pa->adhesion_loc);
// 							distance_boundary(dloc2);
// 							double norm2 = dloc2.norm();
// 							_vec<double> udloc2 = dloc2/norm2;
// 							double cosphi = (udloc1) * (udloc2);
// 							double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
// 							
// //							cout<<std::acos(cosphi)*180.0/PI<<'\t'<<std::asin(sinphi)*180.0/PI<<endl;
// 							
// 							static const double cosphi0 = std::cos(16.0/18.0*PI);
// 							static const double sinphi0 = std::sin(16.0/18.0*PI);
// 							
// 							double F=K_ANGLE_ACTIN_ARP*(cosphi0 - cosphi*sinphi0/sinphi);
// 							
// 							_vec<double> force_sub1 = F*(udloc1 - cosphi * udloc2)/(norm2);
// 							_vec<double> force_sub2 = F*(udloc2 - cosphi * udloc1)/(norm1);
// 							
// 							pa->adhesion_frc -= force_sub2;
// 							arp->force += force_sub1 + force_sub2;
// 							pb->force -= force_sub1;
// 							ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
// 						}
// 					}{//同一平面に
// 						if(pa=arp->edge_ptcl){
// 							_vec<double> dloc_sub=-1.0*dloc1;
// 							_vec<double> r23 =  arp->loc - pa->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
// 							_vec<double> r21;
// 							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// 								r21 = numbering[pa->feature].ptcl_n[pa->num+1]->loc - pa->loc;
// 								distance_boundary(r21);
// 							}else if(pa->num!=0 && numbering[pa->feature].ptcl_n[pa->num-1]){
// 								r21 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 								distance_boundary(r21);
// 							}else if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// 								r21 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// 								distance_boundary(r21);
// 							}else {
// 								fout_err<<"arp_err"<<endl;
// 								continue;
// 							}
// 							_vec<double> n123 = r21 % r23;
// 							_vec<double> n234 = r23 % dloc_sub;
// 							double norm_n1 = n123.norm();
// 							double norm_n2 = n234.norm();
// 							n123/=norm_n1;
// 							n234/=norm_n2;
// 							
// 							double cosph = n123*n234;
// //							double sinph = -1.0*(n123%n234)*ur23;
// //							if(sinph>=0.0)sinph+=1.0e-9;
// //							else sinph-=1.0e-9;
// 							
// 		//					double cosph0 = std::cos(PI);
// 		//					double sinph0 = std::sin(PI);
// 							
// 						//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
// 							double F = K_Torsion_ARP;
// 							
// 							_vec<double> N1 = (n234 - cosph*n123)/norm_n1;
// 							_vec<double> N2 = (n123 - cosph*n234)/norm_n2;
// 							
// 							_vec<double> force_sub1 = F * (N1 % r23);
// 							_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % dloc_sub);
// 							_vec<double> force_sub3 = F * (N2 % r23);
// 							
// 							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// 								numbering[pa->feature].ptcl_n[pa->num+1]->force += force_sub1;
// 								pa->force += force_sub2 - force_sub1;
// 							}else if(pa->num!=0 && numbering[pa->feature].ptcl_n[pa->num-1]){
// 								numbering[pa->feature].ptcl_n[pa->num-1]->force -= force_sub1;
// 								pa->force += force_sub2 + force_sub1;
// 							}else if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// 								numbering[pa->feature].end_arp2_3->force -= force_sub1;
// 								pa->force += force_sub2 + force_sub1;
// 							}else {
// 								fout_err<<"arp_err"<<endl;
// 								continue;
// 							}
// 							arp->force += force_sub3 - force_sub2;
// 							pb->force -= force_sub3;
// 							
// 						}
// 					}{//子フィラメントの結合点を外に向ける
// 						double norm2 = pb->adhesion_loc.norm();
// 						_vec<double> udloc2 = pb->adhesion_loc/norm2;
// 						double cosphi = udloc1 * udloc2;
// 						double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						double F = -1.0*K_TENSILE_ACTIN_ARP*(cosphi/sinphi);
// 						
// 						_vec<double> force_sub1 = F*(udloc1 - cosphi * udloc2)/norm2;
// 						_vec<double> force_sub2 = F*(udloc2 - cosphi * udloc1)/norm1;
// 						
// 						pb->adhesion_frc += force_sub1;
// 						pb->force -= force_sub1 + force_sub2;
// 						arp->force += force_sub2;
// 						ene_pot+=1.0-sinphi;
// 					}{//子フィラメントの回転を止める
// 						if(pa=arp->edge_ptcl){
// 							_vec<double> r23 = pb->loc - arp->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
// 							_vec<double> r21 = pa->loc - arp->loc;	distance_boundary(r21);
// 							_vec<double> n123 = r21 % r23;
// 							_vec<double> n234 = r23 % pb->adhesion_loc;
// 							double norm_n1 = n123.norm();
// 							double norm_n2 = n234.norm();
// 							n123/=norm_n1;
// 							n234/=norm_n2;
// 							
// 							double cosph = n123*n234;
// 	//						double sinph = -1.0*(n123%n234)*ur23;
// 	//						if(sinph>=0.0)sinph+=1.0e-9;
// 	//						else sinph-=1.0e-9;
// 							
// 	//						double cosph0 = std::cos(PI);
// 	//						double sinph0 = std::sin(PI);
// 							
// 						//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
// 							double F = K_Torsion_ARP;
// 							
// 							_vec<double> N1 = (n234 - cosph*n123)/norm_n1;
// 							_vec<double> N2 = (n123 - cosph*n234)/norm_n2;
// 							
// 							_vec<double> force_sub1 = F * (N1 % r23);
// 							_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % pb->adhesion_loc);
// 							_vec<double> force_sub3 = F * (N2 % r23);
// 							pa->force += force_sub1;
// 							arp->force += force_sub2 - force_sub1;
// 							pb->force += force_sub3 - force_sub2;
// 							pb->adhesion_frc -= force_sub3;
// 						}
// 					}
// 					{//子フィラメントの曲げ剛性
// 						_ptcl* pc=NULL;
// 						if(pc=numbering[pb->feature].ptcl_n[1]){
// 							_vec<double> dloc2 = arp->loc - pb->loc;
// 							_vec<double> dloc3 = pc->loc - pb->loc;
// 							distance_boundary(dloc2);
// 							distance_boundary(dloc3);
// 							
// 							double cosphi = (dloc2/dloc2.norm()) * (dloc3/dloc3.norm());
// 						//	double sinphi = std::sqrt(1.0-cosphi*cosphi)+1.0E-9;//0～PI
// 							
// 							_vec<double> force_sub1 = -k_angl*(dloc2/dloc2.norm() - cosphi * dloc3/dloc3.norm())/(dloc3.norm());
// 							_vec<double> force_sub2 = -k_angl*(dloc3/dloc3.norm() - cosphi * dloc2/dloc2.norm())/(dloc2.norm());
// 							
// 							pc->force += force_sub1;
// 							pb->force -= force_sub1 + force_sub2;
// 							arp->force += force_sub2;
// 							ene_pot+=k_angl*(1.0+cosphi);
// 						}
// 					}
// 				}
// // // 				_vec<int> list_index=arp->li->index;
// // 				int li_x=(int)std::floor(arp->loc.x/r_cut);
// // 				int li_y=(int)std::floor(arp->loc.y/r_cut);
// // 				int li_z=(int)std::floor(arp->loc.z/r_cut);
// 				int li_x=(int)std::floor(arp->loc.x/(double)lattice_cut);
// 				int li_y=(int)std::floor(arp->loc.y/(double)lattice_cut);
// 				int li_z=(int)std::floor(arp->loc.z/(double)lattice_cut);
// // if(t>0.001915)cout<<arp->loc<<endl;
// // if(t>0.001915)cout<<li_x<<'\t'<<li_y<<'\t'<<li_z<<endl;
// 				if(li_x==list_X)li_x--;
// 				if(li_y==list_Y)li_y--;
// 				if(li_z==list_Z)li_z--;
// 				for(int li_z_sub=-2;li_z_sub<=2;li_z_sub++){
// 				for(int li_y_sub=-2;li_y_sub<=2;li_y_sub++){
// 				for(int li_x_sub=-2;li_x_sub<=2;li_x_sub++){
// 					int lx = li_x + li_x_sub;
// 					int ly = li_y + li_y_sub;
// 					int lz = li_z + li_z_sub;
// // 					int lx = list_index.x + li_x_sub;
// // 					int ly = list_index.y + li_y_sub;
// // 					int lz = list_index.z + li_z_sub;
// 					while(lx<0)			lx+=list_X;
// 					while(lx>=list_X)	lx-=list_X;
// 					while(ly<0)			ly+=list_Y;
// 					while(ly>=list_Y)	ly-=list_Y;
// #ifdef WALL_BOUNDARY
// 					if(lz<0 || lz>=list_Z)continue;
// #else //WALL_BOUNDARY
// 					while(lz<0)			lz+=list_Z;
// 					while(lz>=list_Z)	lz-=list_Z;
// #endif //WALL_BOUNDARY
// 					_list* lp=(list+ (lx + ly*list_X + lz*(list_X*list_Y) ));
// 					for(int j=0;pa=lp->ptcl_index[j];j++){
// // 						_vec<double> dloc = pa->loc - arp->loc;
// // 						distance_boundary(dloc);
// // 						_vec<double> force_sub = refrce(dloc);
// 						double potential;
// 						_vec<double> force_sub = LJ(pa,arp,&potential);
// 						pa->force += force_sub;
// 						arp->force -= force_sub;
// // 						ene_pot += repot(dloc);
// 						ene_pot += potential;
// 					}
// 				}
// 				}
// 				}
// 			}
// 			return;
// 		}
// 		
// 	}	//namespace _ACTIN_ARP
// 	
// 	void Bind(){
// 		_ptcl* p=NULL;
// 		_arp2_3* arp=NULL;
// 		_list_sub list_sub;
// 		//memset(&list_sub,0,sizeof(_list_sub));
// 		for(int m=0;m<Number_Arp2_3;m++){
// 			arp=(arp2_3+m);
// 			if(!arp->edge_ptcl){
// 				memset(&list_sub,0,sizeof(_list_sub));
// // 				for(int i=1; numbering[i].num; i++){
// // 					for(int j=0; p=numbering[i].ptcl_n[j]; j++){
// // 						if(!p->edge_arp2_3){
// // 							_vec<double> dloc = arp->loc - (p->loc + p->adhesion_loc);
// // 							distance_boundary(dloc);
// // 							if(dloc.norm()<=SIGMA_LJ_ARP_ARP*0.5
// // 								&& (p->adhesion_loc/p->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
// // 								list_sub.index[list_sub.num]=p;
// // 								list_sub.num++;
// // //								cout<<list_sub.num<<endl;
// // 							}
// // 						}
// // 						
// // 					}
// // 				}
// 				_list* li=NULL;
// 				for(int a=0;a<27;a++){
// #ifdef WALL_BOUNDARY
// 					if(li=arp->li->list_index[a]){
// #else
// 					li=arp->li->list_index[a];
// #endif //WALL_BOUNDARY
// 					for(int j=0;p=li->ptcl_index[j];j++){
// 						if((p->feature) && (!p->edge_arp2_3)){
// 							_vec<double> dloc = arp->loc - (p->loc + p->adhesion_loc);
// 							distance_boundary(dloc);
// 							if(dloc.norm()<=SIGMA_LJ_ARP_ARP*0.5
// 								&& (p->adhesion_loc/p->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
// 								list_sub.index[list_sub.num]=p;
// 								list_sub.num++;
// //								cout<<list_sub.num<<endl;
// 							}
// 						}
// 					}
// #ifdef WALL_BOUNDARY
// 					}
// #endif //WALL_BOUNDARY
// 				}
// 				if(list_sub.num==1){
// 					_ptcl* p1=list_sub.index[0];
// 					arp->edge_ptcl = p1;
// 					p1->edge_arp2_3 = arp;
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					_ptcl* p1=list_sub.index[num];
// 					arp->edge_ptcl=p1;
// 					p1->edge_arp2_3 = arp;
// 				}
// 			}
// 			else if(!arp->end_ptcl){
// 				if(!arp->edge_ptcl->feature){
// 					arp->edge_ptcl->edge_arp2_3=NULL;
// 					arp->edge_ptcl=NULL;
// 					continue;
// 				}
// 				memset(&list_sub,0,sizeof(_list_sub));
// 				_vec<double> axis1;
// 				{
// 					_ptcl* pa=NULL;
// 					_ptcl* pb=NULL;
// 					if(pa = arp->edge_ptcl){
// 						
// 						_vec<double> dloc1 = arp->loc - pa->loc;
// 						_vec<double> dloc2;
// 						distance_boundary(dloc1);
// 						dloc1/=dloc1.norm();
// 						if(pa->num==0){
// 							if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 								dloc2 = pb->loc - pa->loc;
// 							}else if(numbering[pa->feature].end_arp2_3){
// 								dloc2 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// 							}else{fout_err<<"arp_err_bind"<<endl;	continue;}
// 						}else if(pb=numbering[pa->feature].ptcl_n[pa->num+1]){
// 							dloc2 = pb->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 						}else{
// 							dloc2 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 						}
// 						distance_boundary(dloc2);
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						static const double costhe0 = std::cos(PI*70.0/180.0);
// 						static const double sinthe0 = std::sin(PI*70.0/180.0);
// 						
// 						axis1 = (costhe0 - sinthe0*costhe/sinthe)*dloc2 + (sinthe0/sinthe)*dloc1;
// 					}else{fout_err<<"arp_err_bind"<<endl;	continue;}
// 				}
// 				_list* li=NULL;
// 				for(int a=0;a<27;a++){
// #ifdef WALL_BOUNDARY
// 					if(li=arp->li->list_index[a]){
// #else
// 					li=arp->li->list_index[a];
// #endif //WALL_BOUNDARY
// 					for(int j=0;p=li->ptcl_index[j];j++){
// 						if(!(p->feature)){
// 							_vec<double> dloc = p->loc - arp->loc;
// 							distance_boundary(dloc);
// 							if(dloc.norm()<=SIGMA_LJ_ACTIN_ARP && (axis1)*(dloc/dloc.norm())>0.5){
// 								list_sub.index[list_sub.num]=p;
// 								list_sub.num++;
// 	//							cout<<list_sub.num<<endl;
// 								
// 							}
// 						}
// 					}
// #ifdef WALL_BOUNDARY
// 					}
// #endif //WALL_BOUNDARY
// 				}
// // 				for(int i=0; i<Number_Ptcl; ++i){
// // 					p=(ptcl+i);
// // 					if(!p->feature){
// // 						_vec<double> dloc = p->loc - arp->loc;
// // 						distance_boundary(dloc);
// // 						if(dloc.norm()<=SIGMA_LJ_ACTIN_ARP && (axis1)*(dloc/dloc.norm())>0.5){
// // 							list_sub.index[list_sub.num]=p;
// // 							list_sub.num++;
// // //							cout<<list_sub.num<<endl;
// // 							
// // 						}
// // 					}
// // 				}
// 				if(list_sub.num==1){
// 					p=list_sub.index[0];
// 					arp->end_ptcl=p;
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					
// 					{
// 						_vec<double> dloc1 = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 						_vec<double> dloc2 = p->loc - arp->loc;
// 						distance_boundary(dloc1);
// 						distance_boundary(dloc2);
// 						dloc1/=dloc1.norm();
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
// 						p->adhesion_loc = axis * sigma_LJ * 0.5;
// 						p->adhesion_loc_inite = p->adhesion_loc;
// 						p->quat.IN(0.0 ,0.0 ,0.0);
// 						p->quat_pre.IN(0.0 ,0.0 ,0.0);
// 						
// 					}
// 					
// 					
// // 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// // 					distance_boundary(axis);
// // 					axis = axis % arp->adhesion_loc;
// // 					_quat<double> quat(0.5*PI,axis);
// // 					quat/=quat.norm();
// // 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// // 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// // 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					p=list_sub.index[num];
// 					arp->end_ptcl=p;
// 					
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					{
// 						_vec<double> dloc1 = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 						_vec<double> dloc2 = p->loc - arp->loc;
// 						distance_boundary(dloc1);
// 						distance_boundary(dloc2);
// 						dloc1/=dloc1.norm();
// 						dloc2/=dloc2.norm();
// 						
// 						double costhe = dloc2*dloc1;
// 						double sinthe = std::sqrt(1-costhe*costhe)+1.0E-9;//0～PI
// 						
// 						_vec<double> axis= (dloc1 - costhe*dloc2) / sinthe;
// 						p->adhesion_loc = axis * sigma_LJ * 0.5;
// 						p->adhesion_loc_inite = p->adhesion_loc;
// 						p->quat.IN(0.0 ,0.0 ,0.0);
// 						p->quat_pre.IN(0.0 ,0.0 ,0.0);
// 					}
// 					
// // 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// // 					distance_boundary(axis);
// // 					axis = axis % arp->adhesion_loc;
// // 					_quat<double> quat(0.5*PI,axis);
// // 					quat/=quat.norm();
// // 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// // 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// // 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// // 					
// 				}
// 				
// 			}
// 			else if(!arp->edge_ptcl->feature){
// 				arp->edge_ptcl->edge_arp2_3=NULL;
// 				numbering[arp->end_ptcl->feature].end_arp2_3=NULL;
// 				arp->edge_ptcl=NULL;
// 				arp->end_ptcl=NULL;
// 			}
// 		}
// 		return;
// 	}
// 	
// 	void Arp_Location_FirstStep(){
// //	cout<<"first"<<endl;
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			_ARP2_3::_arp2_3* p = (arp2_3+i);
// 			p->loc+= (0.5*PTCL_dt* p->force /Arp2_3_Friction);
// 			p->force.IN(0.0,0.0,0.0);
// 			Periodic(p->loc);
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 		return;
// 	}
// 	void Arp_Location_SecondStep(){
// //	cout<<"second"<<endl;
// 		static const double amplitude=sqrt(2.0*SYS_kT*Arp2_3_Friction/PTCL_dt);
// 		static const double root3=sqrt(3.0)*2.0;
// 		_vec<double> RandamF(0.0,0.0,0.0);
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			_ARP2_3::_arp2_3* p = (arp2_3+i);
// 			RandamF= (amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3);
// 			p->loc_pre+=PTCL_dt * (p->force + RandamF) / Arp2_3_Friction;
// 			p->force.IN(0.0,0.0,0.0);
// 			Periodic(p->loc_pre);
// 			p->loc=p->loc_pre;
// 			_list* li=p->list_calc();
// 			if(li!=p->li){
// 				p->li->rm(p,p->list_num);
// 				p->add(li);
// 				li->add(p);
// 			}
// 		}
// 		return;
// 	}
// 
// }
// 
// //#endif //ARP2_3

#endif

#if 0
// 
// 
// #include"_ctrl.h"
// #ifdef ARP2_3
// 
// #include<cstdio>
// #include<cmath>
// #include"_arp2_3.h"
// #include"_periodic.h"
// #include"_variable.h"
// #include "quaternion.h"
// 
// using namespace std;
// template <class T> inline ostream& operator<< (ostream& stream,_vec<T> vec){
// 	return stream<<vec.x<<'\t'<<vec.y<<'\t'<<vec.z;
// }
// template <class T> inline ostream& operator<< (ostream& stream,_quat<T> quat){
// 	return stream<<quat.n<<'\t'<<quat.x<<'\t'<<quat.y<<'\t'<<quat.z;
// }
// 
// namespace _ARP2_3{
// 	_arp2_3* arp2_3=NULL;
// 	 double ene_pot=0.0;
// 	void arp_init(void){
// 		_ARP2_3::arp2_3=new _ARP2_3::_arp2_3[Number_Arp2_3];
// 		int i=0;
// 		for(int x=0;x<SYS_X;x++){
// 			for(int y=0;y<SYS_Y;y++){
// 				for(int z=0;z<SYS_Z;z++){
// 					if(x%16==0 && y%4==2 && z%3==0){
// 						arp2_3[i].loc.IN((double)x+0.5,(double)y+0.5,(double)z+1.0);
// 						arp2_3[i].loc_pre=arp2_3[i].loc;
// 						arp2_3[i].index=i;
// //						cout<<arp2_3[i].loc<<'\t'<<arp2_3[i].loc_pre<<endl;
// 						i++;
// 						if(i==Number_Arp2_3)break;
// 					}
// 				}if(i==Number_Arp2_3)break;
// 			}if(i==Number_Arp2_3)break;
// 		}
// 		
// 		return;
// 	}
// 	namespace _ARP_ARP{
// 		inline double pot(const double r){
// 			return (4.0*EPSILON_LJ_ARP_ARP)*pow(SIGMA_LJ_ARP_ARP/r,12.0);
// 		}
// 		inline double dif_pot( double r){
// 			return -1.0*(4.0*EPSILON_LJ_ARP_ARP)*12.0/SIGMA_LJ_ARP_ARP*pow(SIGMA_LJ_ARP_ARP/r,13.0);
// 		}
// 		inline double d_dif_pot( double r){
// 			return (4.0*EPSILON_LJ_ARP_ARP)*156.0/(SIGMA_LJ_ARP_ARP*SIGMA_LJ_ARP_ARP)*pow(SIGMA_LJ_ARP_ARP/r,14.0);
// 		}
// 		
// 		inline double repot(_vec<double> r){
// 			if(r.norm()<=Cut_Off_Arp2_3_Arp2_3){
// 				return pot(r.norm())-pot((double)Cut_Off_Arp2_3_Arp2_3)
// 					+dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)
// 					+0.5*d_dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3);
// 			}
// 			else{return 0.0;}
// 		}
// 		inline _vec<double> refrce(_vec<double> r){//LJポテンシャル
// 			if(r.norm()<=Cut_Off_Arp2_3_Arp2_3){
// 				return -( dif_pot(r.norm())
// 					+dif_pot((double)Cut_Off_Arp2_3_Arp2_3)
// 					+d_dif_pot((double)Cut_Off_Arp2_3_Arp2_3)*(r.norm()-(double)Cut_Off_Arp2_3_Arp2_3)
// 					)*r/r.norm();
// 			}
// 			else{return _vec<double>(0.0,0.0,0.0);}
// 		}
// 		void Frc(){
// 			for(int i=0;i<Number_Arp2_3;i++){
// 				for(int j=i+1;j<Number_Arp2_3;j++){
// 					_vec<double> dloc = _ARP2_3::arp2_3[i].loc - _ARP2_3::arp2_3[j].loc;
// 					distance_boundary(dloc);
// 					_vec<double> force_sub = refrce(dloc);
// 					_ARP2_3::arp2_3[i].force += force_sub;
// 					_ARP2_3::arp2_3[j].force -= force_sub;
// 					ene_pot += repot(dloc);
// 				}
// 			}
// 			return;
// 		}
// 	}	//namespace _ARP_ARP
// 	
// 	
// 	namespace _ACTIN_ARP{
// 		inline double pot(double r){//ポテンシャル(スカラー)
// 			return (4.0*EPSILON_LJ_ACTIN_ARP)*pow(SIGMA_LJ_ACTIN_ARP/r,12.0);
// 		}
// 		inline double dif_pot(double r){
// 			return -1.0*(4.0*EPSILON_LJ_ACTIN_ARP)*12.0/SIGMA_LJ_ACTIN_ARP*pow(SIGMA_LJ_ACTIN_ARP/r,13.0);
// 		}
// 		inline double d_dif_pot(double r){
// 			return (4.0*EPSILON_LJ_ACTIN_ARP)*156.0/(SIGMA_LJ_ACTIN_ARP*SIGMA_LJ_ACTIN_ARP)*pow(SIGMA_LJ_ACTIN_ARP/r,14.0);
// 		}
// 
// 		inline double repot(_vec<double> r){
// 			if(r.norm()<=Cut_Off_Actin_Arp2_3){
// 				return pot(r.norm())
// 					-pot(Cut_Off_Actin_Arp2_3)
// 					+dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)
// 					+0.5*d_dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3);
// 			}
// 			else{return 0.0;}
// 		}
// 		inline _vec<double> refrce( _vec<double> r){//LJポテンシャル
// 			if(r.norm()<=Cut_Off_Actin_Arp2_3){
// 				return -( dif_pot(r.norm())
// 					+dif_pot(Cut_Off_Actin_Arp2_3)
// 					+d_dif_pot(Cut_Off_Actin_Arp2_3)*(r.norm()-Cut_Off_Actin_Arp2_3)
// 					)*r/r.norm();
// 			}
// 			else{return _vec<double>(0.0,0.0,0.0);}
// 		}
// 		void Frc(){
// 			_ptcl* pa=NULL;
// 			_ptcl* pb=NULL;
// 			_arp2_3* arp=NULL;
// 			for(int i=0;i<Number_Arp2_3;i++){
// 				arp=(arp2_3+i);
// 				
// 				if(pa=arp->edge_ptcl){
// 					_vec<double> dloc1 = (pa->loc + pa->adhesion_loc) - arp->loc;
// 					distance_boundary(dloc1);
// 					{
// 						_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(dloc1.norm()-EQ_LEN_ACTIN_ARP)*dloc1/dloc1.norm();
// 						pa->adhesion_frc += force_sub1;
// 						arp->force -= force_sub1;
// 						ene_pot+=K_TENSILE_ACTIN_ARP*(dloc1.norm()-EQ_LEN_ACTIN_ARP)*(dloc1.norm()-EQ_LEN_ACTIN_ARP)/2.0;
// 					}
// 					{
// 						double cosphi = (dloc1/dloc1.norm()) * (pa->adhesion_loc/pa->adhesion_loc.norm());
// 						//double sinphi = std::sqrt(1-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						_vec<double> force_sub1 = K_ANGLE_ACTIN_ARP*(dloc1/dloc1.norm() - cosphi * pa->adhesion_loc/pa->adhesion_loc.norm())/(pa->adhesion_loc.norm());
// 						_vec<double> force_sub2 = K_ANGLE_ACTIN_ARP*(pa->adhesion_loc/pa->adhesion_loc.norm() - cosphi * dloc1/dloc1.norm())/(dloc1.norm());
// 						//dif_pot=K*sin(phi-phi0)=k*sin(phi-pi)=-k*sinphi
// 						pa->force += force_sub1;
// 						pa->adhesion_frc -= force_sub1 + force_sub2;
// 						arp->force += force_sub2;
// 						ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
// 					}
// 					{
// 						double cosphi = (dloc1/dloc1.norm()) * (arp->adhesion_loc/arp->adhesion_loc.norm());
// 						double sinphi = std::sqrt(1-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						static const double cosphi0 = std::cos(16.0/18.0*PI);
// 						static const double sinphi0 = std::sin(16.0/18.0*PI);
// 						
// 						_vec<double> force_sub1 = K_ANGLE_ACTIN_ARP*(sinphi*cosphi0 - cosphi*sinphi0)
// 								*(dloc1/dloc1.norm() - cosphi * arp->adhesion_loc/arp->adhesion_loc.norm())/(sinphi*arp->adhesion_loc.norm());
// 						_vec<double> force_sub2 = K_ANGLE_ACTIN_ARP*(sinphi*cosphi0 - cosphi*sinphi0)
// 								*(arp->adhesion_loc/arp->adhesion_loc.norm() - cosphi * dloc1/dloc1.norm())/(sinphi*dloc1.norm());
// 						
// 						arp->adhesion_frc += force_sub1;
// 						arp->force -= force_sub1 + force_sub2;
// 						pa->adhesion_frc += force_sub2;
// 						ene_pot+=K_ANGLE_ACTIN_ARP*(1.0-cosphi*cosphi0-sinphi*sinphi0);
// 					}
// 					{
// //						if(numbering[pa->feature].num!=1){
// 							_vec<double> r23 =  arp->loc - pa->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
// 							_vec<double> r21;
// //							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// //								r21 = numbering[pa->feature].ptcl_n[pa->num+1]->loc - pa->loc;
// //								distance_boundary(r21);
// //							}else {
// //								if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// //									r21 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// //								}
// //								else {r21 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;}
// //								distance_boundary(r21);
// //							}
// 							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// 								r21 = numbering[pa->feature].ptcl_n[pa->num+1]->loc - pa->loc;
// 								distance_boundary(r21);
// 							}else if(pa->num!=0 && numbering[pa->feature].ptcl_n[pa->num-1]){
// 								r21 = pa->loc - numbering[pa->feature].ptcl_n[pa->num-1]->loc;
// 								distance_boundary(r21);
// 							}else if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// 								r21 = pa->loc - numbering[pa->feature].end_arp2_3->loc;
// 								distance_boundary(r21);
// 							}else {
// 								fout_err<<"arp_err"<<endl;
// 								continue;
// 							}
// 							_vec<double> n123 = r21 % r23;
// 							_vec<double> n234 = r23 % arp->adhesion_loc;
// 							double norm1 = n123.norm();
// 							double norm2 = n234.norm();
// 							n123/=norm1;
// 							n234/=norm2;
// 							
// 							double cosph = n123*n234;
// //							double sinph = -1.0*(n123%n234)*ur23;
// //							if(sinph>=0.0)sinph+=1.0e-9;
// //							else sinph-=1.0e-9;
// 							
// 		//					double cosph0 = std::cos(PI);
// 		//					double sinph0 = std::sin(PI);
// 							
// 						//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
// 							double F = K_Torsion_ARP;
// 							
// 							_vec<double> N1 = (n234 - cosph*n123)/norm1;
// 							_vec<double> N2 = (n123 - cosph*n234)/norm2;
// 							
// 							_vec<double> force_sub1 = F * (N1 % r23);
// 							_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % arp->adhesion_loc);
// 							_vec<double> force_sub3 = F * (N2 % r23);
// 							
// //							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// //								numbering[pa->feature].ptcl_n[pa->num+1]->force += force_sub1;
// //								pa->force += force_sub2 - force_sub1;
// //							}else {
// //								if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// //									numbering[pa->feature].end_arp2_3->force -= force_sub1;
// //									pa->force += force_sub2 + force_sub1;
// //								}
// //								else {
// //									numbering[pa->feature].ptcl_n[pa->num-1]->force -= force_sub1;
// //									pa->force += force_sub2 + force_sub1;
// //								}
// //							}
// 							if(numbering[pa->feature].ptcl_n[pa->num+1]){
// 								numbering[pa->feature].ptcl_n[pa->num+1]->force += force_sub1;
// 								pa->force += force_sub2 - force_sub1;
// 							}else if(pa->num!=0 && numbering[pa->feature].ptcl_n[pa->num-1]){
// 								numbering[pa->feature].ptcl_n[pa->num-1]->force -= force_sub1;
// 								pa->force += force_sub2 + force_sub1;
// 							}else if(numbering[pa->feature].num==1 && numbering[pa->feature].end_arp2_3!=NULL){
// 								numbering[pa->feature].end_arp2_3->force -= force_sub1;
// 								pa->force += force_sub2 + force_sub1;
// 							}else {
// 								fout_err<<"arp_err"<<endl;
// 								continue;
// 							}
// 							arp->force += force_sub3 - force_sub2;
// 							arp->adhesion_frc -= force_sub3;
// 							
// 						}
// //					}
// 				}
// 				if(pb=arp->end_ptcl){
// 					_vec<double> dloc1 = (arp->loc + arp->adhesion_loc) - pb->loc;
// 					distance_boundary(dloc1);
// 					{
// 						_vec<double> force_sub1 = -1.0*K_TENSILE_ACTIN_ARP*(dloc1.norm()-0.5*sigma_LJ)*dloc1/dloc1.norm();
// 						arp->adhesion_frc += force_sub1;
// 						pb->force -= force_sub1;
// 						ene_pot+=K_TENSILE_ACTIN_ARP*(dloc1.norm()-0.5*sigma_LJ)*(dloc1.norm()-0.5*sigma_LJ)/2.0;
// 					}
// 					{
// 						double cosphi = (dloc1/dloc1.norm()) * (arp->adhesion_loc/arp->adhesion_loc.norm());
// 					//	double sinphi = std::sqrt(1-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						_vec<double> force_sub1 = K_ANGLE_ACTIN_ARP*(dloc1/dloc1.norm() - cosphi * arp->adhesion_loc/arp->adhesion_loc.norm())/(arp->adhesion_loc.norm());
// 						_vec<double> force_sub2 = K_ANGLE_ACTIN_ARP*(arp->adhesion_loc/arp->adhesion_loc.norm() - cosphi * dloc1/dloc1.norm())/(dloc1.norm());
// 						
// 						arp->force += force_sub1;
// 						arp->adhesion_frc -= force_sub1 + force_sub2;
// 						pb->force += force_sub2;
// 						ene_pot+=K_ANGLE_ACTIN_ARP*(1.0+cosphi);
// 					}
// 					{
// 						double cosphi = (dloc1/dloc1.norm()) * (pb->adhesion_loc/pb->adhesion_loc.norm());
// 						double sinphi = std::sqrt(1-cosphi*cosphi)+1.0E-9;//0～PI
// 						
// 						_vec<double> force_sub1 = K_TENSILE_ACTIN_ARP*( - cosphi)
// 								*(dloc1/dloc1.norm() - cosphi * pb->adhesion_loc/pb->adhesion_loc.norm())/(sinphi*pb->adhesion_loc.norm());
// 						_vec<double> force_sub2 = K_TENSILE_ACTIN_ARP*( - cosphi)
// 								*(pb->adhesion_loc/pb->adhesion_loc.norm() - cosphi * dloc1/dloc1.norm())/(sinphi*dloc1.norm());
// 						
// 						pb->adhesion_frc += force_sub1;
// 						pb->force -= force_sub1 + force_sub2;
// 						arp->adhesion_frc += force_sub2;
// 						ene_pot+=1.0-sinphi;
// 					}
// 					{
// 						_vec<double> r23 = pb->loc - arp->loc;	distance_boundary(r23);	_vec<double> ur23 = r23 / r23.norm();
// 						_vec<double> r21 = arp->edge_ptcl->loc - arp->loc;
// 						distance_boundary(r21);
// 						_vec<double> n123 = r21 % r23;
// 						_vec<double> n234 = r23 % pb->adhesion_loc;
// 						double norm1 = n123.norm();
// 						double norm2 = n234.norm();
// 						n123/=norm1;
// 						n234/=norm2;
// 						
// 						double cosph = n123*n234;
// //						double sinph = -1.0*(n123%n234)*ur23;
// //						if(sinph>=0.0)sinph+=1.0e-9;
// //						else sinph-=1.0e-9;
// 						
// //						double cosph0 = std::cos(PI);
// //						double sinph0 = std::sin(PI);
// 						
// 					//	double F = -1.0 * K_Torsion_ARP*(sinph*cosph0 - cosph*sinph0)/sinph;
// 						double F = K_Torsion_ARP;
// 						
// 						_vec<double> N1 = (n234 - cosph*n123)/norm1;
// 						_vec<double> N2 = (n123 - cosph*n234)/norm2;
// 						
// 						_vec<double> force_sub1 = F * (N1 % r23);
// 						_vec<double> force_sub2 = F * (N1 % r21) -  F * (N2 % pb->adhesion_loc);
// 						_vec<double> force_sub3 = F * (N2 % r23);
// 					//	cout << "f1 " << force_sub1 << endl;
// 					//	cout << "f2 " << force_sub2 << endl;
// 					//	cout << "f3 " << force_sub3 << endl;
// 					//	if(t>0.1)exit(0);
// 						arp->edge_ptcl->force += force_sub1;
// 						arp->force += force_sub2 - force_sub1;
// 						pb->force += force_sub3 - force_sub2;
// 						pb->adhesion_frc -= force_sub3;
// 						
// 					}
// 					{
// 						_ptcl* pc=NULL;
// 						if(pc=numbering[pb->feature].ptcl_n[1]){
// 							_vec<double> dloc2 = arp->loc - pb->loc;
// 							_vec<double> dloc3 = pc->loc - pb->loc;
// 							distance_boundary(dloc2);
// 							distance_boundary(dloc3);
// 							
// 							double cosphi = (dloc2/dloc2.norm()) * (dloc3/dloc3.norm());
// 						//	double sinphi = std::sqrt(1-cosphi*cosphi)+1.0E-9;//0～PI
// 							
// 							_vec<double> force_sub1 = -k_angl*(dloc2/dloc2.norm() - cosphi * dloc3/dloc3.norm())/(dloc3.norm());
// 							_vec<double> force_sub2 = -k_angl*(dloc3/dloc3.norm() - cosphi * dloc2/dloc2.norm())/(dloc2.norm());
// 							
// 							pc->force += force_sub1;
// 							pb->force -= force_sub1 + force_sub2;
// 							arp->force += force_sub2;
// 							ene_pot+=k_angl*(1.0+cosphi);
// 						}
// 					}
// 				}
// 				
// 				int li_x=(int)std::floor(arp->loc.x/r_cut);
// 				int li_y=(int)std::floor(arp->loc.y/r_cut);
// 				int li_z=(int)std::floor(arp->loc.z/r_cut);
// 				if(li_x==list_X)li_x--;
// 				if(li_y==list_Y)li_y--;
// 				if(li_z==list_Z)li_z--;
// 				for(int li_z_sub=-2;li_z_sub<=2;li_z_sub++){
// 				for(int li_y_sub=-2;li_y_sub<=2;li_y_sub++){
// 				for(int li_x_sub=-2;li_x_sub<=2;li_x_sub++){
// 					int lx = li_x + li_x_sub;
// 					int ly = li_y + li_y_sub;
// 					int lz = li_z + li_z_sub;
// 					while(lx<0)			lx+=list_X;
// 					while(lx>=list_X)	lx-=list_X;
// 					while(ly<0)			ly+=list_Y;
// 					while(ly>=list_Y)	ly-=list_Y;
// #ifdef WALL_BOUNDARY
// 					if(lz<0 || lz>=list_Z)continue;
// #else //WALL_BOUNDARY
// 					while(lz<0)			lz+=list_Z;
// 					while(lz>=list_Z)	lz-=list_Z;
// #endif //WALL_BOUNDARY
// 					_list* lp=(list+ (lx + ly*list_X + lz*(list_X*list_Y) ));
// 					for(int j=0;pa=lp->ptcl_index[j];j++){
// 						_vec<double> dloc = pa->loc - arp->loc;
// 						distance_boundary(dloc);
// 						_vec<double> force_sub = refrce(dloc);
// 						pa->force += force_sub;
// 						arp->force -= force_sub;
// 						ene_pot += repot(dloc);
// 					}
// 				}
// 				}
// 				}
// 			}
// 			return;
// 		}
// 		
// 	}	//namespace _ACTIN_ARP
// 	
// 	void Bind(){
// 		_ptcl* p=NULL;
// 		_arp2_3* arp=NULL;
// 		_list_sub list_sub;
// 		//memset(&list_sub,0,sizeof(_list_sub));
// 		for(int m=0;m<Number_Arp2_3;m++){
// 			arp=(arp2_3+m);
// 			if(!arp->edge_ptcl){
// 				memset(&list_sub,0,sizeof(_list_sub));
// 				for(int i=1; numbering[i].num; i++){
// 					for(int j=0; p=numbering[i].ptcl_n[j]; j++){
// 						if(!p->edge_arp2_3){
// 							_vec<double> dloc = arp->loc - (p->loc + p->adhesion_loc);
// 							distance_boundary(dloc);
// 							if(dloc.norm()<=SIGMA_LJ_ARP_ARP*0.5
// 								&& (p->adhesion_loc/p->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
// 								list_sub.index[list_sub.num]=p;
// 								list_sub.num++;
// //								cout<<list_sub.num<<endl;
// 							}
// 						}
// 						
// 					}
// 				}
// 				if(list_sub.num==1){
// 					_ptcl* p=list_sub.index[0];
// 					arp->edge_ptcl = p;
// 					p->edge_arp2_3 = arp;
// 					_vec<double> axis;
// 					_quat<double> quat;
// 					if(p->num==0){
// 						if(numbering[p->feature].num==1 && numbering[p->feature].end_arp2_3!=NULL){
// 							axis = (p->loc - numbering[p->feature].end_arp2_3->loc);
// 						}
// 						else {axis = (numbering[p->feature].ptcl_n[1]->loc - p->loc);}
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					}else if(!numbering[p->feature].ptcl_n[p->num+1]){
// 						axis = (p->loc - numbering[p->feature].ptcl_n[p->num-1]->loc);
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 						
// 					}else{
// 						axis = (numbering[p->feature].ptcl_n[p->num+1]->loc - numbering[p->feature].ptcl_n[p->num-1]->loc);
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 						
// 					}
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					_ptcl* p=list_sub.index[num];
// 					arp->edge_ptcl=p;
// 					p->edge_arp2_3 = arp;
// 					_vec<double> axis;
// 					_quat<double> quat;
// 					if(p->num==0){
// 						if(numbering[p->feature].num==1 && numbering[p->feature].end_arp2_3!=NULL){
// 							axis = (p->loc - numbering[p->feature].end_arp2_3->loc);
// 						}
// 						else {axis = (numbering[p->feature].ptcl_n[1]->loc - p->loc);}
// //						axis = (numbering[p->feature].ptcl_n[1]->loc - p->loc);
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					}else if(!numbering[p->feature].ptcl_n[p->num+1]){
// 						axis = (p->loc - numbering[p->feature].ptcl_n[p->num-1]->loc);
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 						
// 					}else{
// 						axis = (numbering[p->feature].ptcl_n[p->num+1]->loc - numbering[p->feature].ptcl_n[p->num-1]->loc);
// 						distance_boundary(axis);
// 						axis = axis % p->adhesion_loc;
// 						quat.IN(16/17*PI,axis);
// 						quat/=quat.norm();
// 						_vec<double> TMP_LOC=p->adhesion_loc - arp->loc;
// 						TMP_LOC*=0.5 * SIGMA_LJ_ARP_ARP / TMP_LOC.norm();
// 						arp->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 						
// 					}
// 				}
// 			}
// 			else if(!arp->end_ptcl){
// 				if(!arp->edge_ptcl->feature){
// 					arp->edge_ptcl->edge_arp2_3=NULL;
// 					arp->edge_ptcl=NULL;
// 					continue;
// 				}
// 				memset(&list_sub,0,sizeof(_list_sub));
// 				for(int i=0; i<Number_Ptcl; i++){
// 					p=(ptcl+i);
// 					if(!p->feature){
// 						_vec<double> dloc = p->loc - (arp->loc + arp->adhesion_loc);
// 						distance_boundary(dloc);
// 						if(dloc.norm()<=sigma_LJ*0.5&& (arp->adhesion_loc/arp->adhesion_loc.norm())*(dloc/dloc.norm())>0.5){
// 							list_sub.index[list_sub.num]=p;
// 							list_sub.num++;
// //							cout<<list_sub.num<<endl;
// 							
// 						}
// 					}
// 				}
// 				if(list_sub.num==1){
// 					_ptcl* p=list_sub.index[0];
// 					arp->end_ptcl=p;
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					
// 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 					distance_boundary(axis);
// 					axis = axis % arp->adhesion_loc;
// 					_quat<double> quat(0.5*PI,axis);
// 					quat/=quat.norm();
// 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					
// 				}else if(list_sub.num>1){
// 					int num = (int)std::floor((double)list_sub.num * RAND);
// 					while(num<0 || num>=list_sub.num){num = (int)std::floor((double)list_sub.num * RAND);}
// 					_ptcl* p=list_sub.index[num];
// 					arp->end_ptcl=p;
// 					
// 					numbering[0].sort_n(p->num);
// 					numbering[0].num--;
// 					no_fila++;
// 					p->feature = no_fila;
// 					p->num = 0;
// 					p->time = t;
// 					numbering[no_fila].ptcl_n[0]=p;
// 					numbering[no_fila].num=1;
// 					numbering[no_fila].end_arp2_3=arp;
// 					
// 					
// 					_vec<double> axis = (arp->edge_ptcl->loc + arp->edge_ptcl->adhesion_loc) - arp->loc;
// 					distance_boundary(axis);
// 					axis = axis % arp->adhesion_loc;
// 					_quat<double> quat(0.5*PI,axis);
// 					quat/=quat.norm();
// 					_vec<double> TMP_LOC = (arp->loc + arp->adhesion_loc) - p->loc;
// 					TMP_LOC*=0.5*sigma_LJ/TMP_LOC.norm();
// 					p->adhesion_loc = _quaternion::rot( quat , TMP_LOC);
// 					
// 				}
// 			}
// 			else if(!arp->edge_ptcl->feature){
// 				arp->edge_ptcl->edge_arp2_3=NULL;
// 				numbering[arp->end_ptcl->feature].end_arp2_3=NULL;
// 				arp->edge_ptcl=NULL;
// 				arp->end_ptcl=NULL;
// 			}
// 		}
// 		return;
// 	}
// 	
// 	void Arp_Adhesion_Location_FirstStep(){
// 		_ARP2_3::_arp2_3* arp=NULL;
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			arp=(arp2_3+i);
// 			if(arp->edge_ptcl!=NULL ||arp->end_ptcl!=NULL){
// 				arp->quat_pre = arp->quat;
// 				_vec<double> torque = arp->adhesion_loc % arp->adhesion_frc;
// 				_vec<double> omega = torque/Arp2_3_Ror_Friction;
// 				
// 				arp->quat += (0.25 * (omega * arp->quat) * PTCL_dt);//dq/dt=0.5*ω*q
// 				arp->quat /= (arp->quat.norm());
// 				
// 				arp->adhesion_loc = _quaternion::rot(arp->quat,arp->adhesion_loc);
// 				arp->adhesion_frc.IN(0.0,0.0,0.0);
// 			}
// 		}
// 		return;
// 	}
// 	void Arp_Adhesion_Location_SecondStep(){
// 		static const double amplitude=sqrt(2.0*SYS_kT*Arp2_3_Ror_Friction/PTCL_dt);
// 		static const double root3=sqrt(3.0)*2.0;
// 		_vec<double> RandamF(0.0,0.0,0.0);
// 		_ARP2_3::_arp2_3* p;
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			p = (arp2_3+i);
// 			if(p->edge_ptcl!=NULL ||p->end_ptcl!=NULL){
// 				RandamF=amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3;
// 				_vec<double> torque = p->adhesion_loc % p->adhesion_frc + RandamF;
// 				_vec<double> omega = torque/Arp2_3_Ror_Friction;
// 				
// 				p->quat = p->quat_pre + (0.5 * (omega * p->quat) * PTCL_dt);//dq/dt=0.5*ω*q
// 				p->quat /= (p->quat.norm());
// 				
// 				p->adhesion_loc = _quaternion::rot(p->quat,p->adhesion_loc);
// 				p->adhesion_frc.IN(0.0,0.0,0.0);
// 			}
// 		}
// 		return;
// 	}
// 	void Arp_Location_FirstStep(){
// //	cout<<"first"<<endl;
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			_ARP2_3::_arp2_3* p = (arp2_3+i);
// 			p->loc+=0.5*PTCL_dt*( p->force + p->adhesion_frc )/Arp2_3_Friction;
// 			p->force.IN(0.0,0.0,0.0);
// 			Periodic(p->loc);
// 		}
// 		Arp_Adhesion_Location_FirstStep();
// 		return;
// 	}
// 	void Arp_Location_SecondStep(){
// //	cout<<"second"<<endl;
// 		static const double amplitude=sqrt(2.0*SYS_kT*Arp2_3_Friction/PTCL_dt);
// 		static const double root3=sqrt(3.0)*2.0;
// 		_vec<double> RandamF(0.0,0.0,0.0);
// 		for(int i=0; i<Number_Arp2_3 ;i++){
// 			_ARP2_3::_arp2_3* p = (arp2_3+i);
// 			RandamF=amplitude * _vec<double>(RAND-0.5,RAND-0.5,RAND-0.5) * root3;
// 			p->loc_pre+=PTCL_dt * (p->force + p->adhesion_frc + RandamF) / Arp2_3_Friction;
// 			p->force.IN(0.0,0.0,0.0);
// 			Periodic(p->loc_pre);
// 			p->loc=p->loc_pre;
// 		}
// 		Arp_Adhesion_Location_SecondStep();
// 		return;
// 	}
// 
// }
// 
// #endif //ARP2_3

#endif

#endif //ARP2_3
