/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include<iostream>
#include<fstream>
#include<cmath>
#include "_system.h"
#include "_variable.h"
#include "_ctrl.h"
void debug();

void Periodic(_vec<double>& a){
	int nx = (int)std::floor(a.x/(double)SYS_X);
	int ny = (int)std::floor(a.y/(double)SYS_Y);
#ifndef WALL_BOUNDARY
	int nz = (int)std::floor(a.z/(double)SYS_Z);
#endif
	
	if(nx){a.x  -= (double)SYS_X*(double)nx;}
	if(ny){a.y  -= (double)SYS_Y*(double)ny;}
#ifndef WALL_BOUNDARY
	if(nz){a.z  -= (double)SYS_Z*(double)nz;}
#endif
	//a-=_vec<double>((double)SYS_X*floor(a.x/(double)SYS_X),
	//	(double)SYS_Y*floor(a.y/(double)SYS_Y),(double)SYS_Z*floor(a.z/(double)SYS_Z))
#ifdef DEBUG_PERIODIC	
	if(std::abs(nx)>2){
		int index = (_ptcl*)&a-(_ptcl*)&ptcl[0].loc;//((char*)&a-(char*)&ptcl[0].loc)/sizeof(_ptcl)
		if(index<0 ||index>=Number_Ptcl)index=-1;
		cout<<"\nMDptcl_fly_over "<<nx<<" systems_at_x.\n"
			<<"ptcl_number = "<<index<<" , step = "<<step<<"\nexit."<<endl;
		
		atexit(debug);
		exit(0);
	}
	if(std::abs(ny)>2){
		int index = (_ptcl*)&a-(_ptcl*)&ptcl[0].loc;
		if(index<0 ||index>=Number_Ptcl)index=-1;
		cout<<"\nMDptcl_fly_over "<<ny<<" systems_at_y.\n"
			<<"ptcl_number = "<<index<<" , step = "<<step<<"\nexit."<<endl;
		atexit(debug);
		exit(0);
	}
#ifndef WALL_BOUNDARY
	if(std::abs(nz)>2){
		int index = (_ptcl*)&a-(_ptcl*)&ptcl[0].loc;
		if(index<0 ||index>=Number_Ptcl)index=-1;
		cout<<"\nMDptcl_fly_over "<<nz<<" systems_at_z.\n"
			<<"ptcl_number = "<<index<<" , step = "<<step<<"\nexit."<<endl;
		atexit(debug);
		exit(0);
	}
#endif
#endif
	return;
}
void distance_boundary(_vec<double>& r){
	if(r.x < -SYS_X/2.0	){r.x += (double)SYS_X;}
	if(r.x > SYS_X/2.0	){r.x -= (double)SYS_X;}
	if(r.y < -SYS_Y/2.0	){r.y += (double)SYS_Y;}
	if(r.y > SYS_Y/2.0	){r.y -= (double)SYS_Y;}
#ifndef WALL_BOUNDARY
	if(r.z < -SYS_Z/2.0	){r.z += (double)SYS_Z;}
	if(r.z > SYS_Z/2.0	){r.z -= (double)SYS_Z;}
#endif
	return;
}
