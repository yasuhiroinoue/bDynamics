/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
//#include <stdio.h>
//#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "mymath.h"
#include "_system.h"
#include "_class.h"
#include "_function.h"
#include "_variable.h"
#include "_ctrl.h"
#include "_arp2_3.h"
#include "_mympi.h"

using namespace std;
std::ostream& OUTPUT_N(std::ostream& stream);

#ifdef USE_MPI
#include<mpi.h>
#endif// USE_MPI

#ifndef SRD
#ifdef BROWNIAN
int main(int argc,char** argv){//const を付けるとMPI_Initで怒られますね..
	int myrank(0);
	int numprocs(1);
#ifdef USE_MPI
	mympi::mympi_init(argc,argv,myrank,numprocs);
#endif //USE_MPI
	
	
	InputInitefile(argc,argv,myrank,numprocs);//設定ファイルの読み込み
	
	
	if((NX * NY * NZ)!=numprocs)return 1;
//InputInitefileで補正しているのでここで引っかかることはないはず？
//	init_genrand(myrank + numprocs*3 + N);
	
	calcArea.add(myrank,numprocs);//各計算機の計算領域がどこか計算
	if(!myrank){//0番期なら
		dir_init();//ディレクトリを作る
		std::stringstream ss_fout;
		ss_fout<<"./output"<<setw(4)<<setfill('0')<<N<<"/fout_err.dat";
		fout_err.open(ss_fout.str().c_str(),ios::out);
	}	
	struct_init(myrank,numprocs);			//構造体の初期化
	if(!myrank){
		file_init();		//ヘッダの書き出し
		std::stringstream ss;	ss.str("");
//		ss<<"Wall_Friction "<<WALL::Wall_Friction<<" - load "<<WALL::jijuu;
		history(ss.str());	//実行した時間の書き出し
		ptcl_initial_condition();			//0番機で配置後各機に送る．乱数が違うので
		PV_out::PVWin_init_bin(myrank,numprocs);	//pvwinの初期化
	}
//		PV_out::PVWin_init_bin(myrank,numprocs);
#ifdef USE_MPI
	mympi::mympi_ptcl_init(myrank,numprocs);
#endif //USE_MPI
if(RESTART){//リスタート
//cout<<myrank<<" "<<t<<endl;
		restart_binary_input(myrank,numprocs,RESTART_NUMBER);
//cout<<myrank<<" "<<t<<endl;
//		if(!myrank){
//				file_init();
//	//			PVWin_init_restart_bin();
//				PV_out::PVWin_init_bin(myrank,numprocs);
//				PV_out::PVWin_output_bin(myrank,numprocs);
//				systemmonitor_ptcl(step);
//				arg_out();	len_out();
//		}
//			}
}
#ifdef USE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif //USE_MPI
	
	for( ; step<STEP_END ; step++){
		for(int step_mpc=0 ; step_mpc<INTVL_step ; step_mpc++){//INTVL_step=1なので意味ない
			for(int i=0 ; i<(int)(INTVL_COL) ; i++){//INTVL_COL=100
				ene_format();			//
				Frc(myrank,numprocs);	//(力の計算)
#ifdef USE_MPI
				mympi::mympi_send_pre(myrank,numprocs);			//糊代をもらい
#endif //USE_MPI
				Ptcl_Location_FirstStep(myrank,numprocs);
#ifdef USE_MPI
				mympi::mympi_send_bck(myrank,numprocs);			//糊代を返す
#endif //USE_MPI
				ene_format();
				Frc(myrank,numprocs);
#ifdef USE_MPI
				mympi::mympi_send_pre(myrank,numprocs);
#endif //USE_MPI
				Ptcl_Location_SecondStep(myrank,numprocs);
#ifdef USE_MPI
				mympi::mympi_send_bck_second(myrank,numprocs);
#endif //USE_MPI
				
				if(step>14000){//初期状態を忘れるまで適当にまわす
#ifdef POLYMERIZATION
// if(t>=7.079)cout<<"7"<<endl;
//					_Actin::nucleation(myrank,numprocs);					//核生成
//#ifdef USE_MPI
//					mympi::mympi_conformation_change_nucleation(myrank,numprocs);
//#endif// USE_MPI
					_Actin::polymerization(myrank,numprocs);				//重合

#endif //POLYMERIZATION
#ifdef DEPOLYMERIZATION
// if(t>=7.079)cout<<"8"<<endl;
//					_Actin::depolymerization(myrank,numprocs);
#endif //DEPOLYMERIZATION
#ifdef SEVERING
//					_Actin::severing(myrank,numprocs);
#endif //SEVERING
					//bind_cal(myrank,numprocs);
#ifdef ARP2_3
//if(t>=7.0789)cout<<myrank<<" 9"<<endl;
//					_ARP2_3::Bind(myrank,numprocs);
#endif //ARP2_3
#ifdef USE_MPI
//if(t>=7.0789)cout<<myrank<<" 10"<<endl;
					mympi::mympi_conformation_change(myrank,numprocs);			//構造変化の共有
//if(t>=7.0789)cout<<myrank<<" 11"<<endl;
#endif// USE_MPI
				}
// if(t>=7.079)cout<<"10"<<endl;
// 				p_level_cal();
// if(t>=7.079)cout<<"11"<<endl;
				t += PTCL_dt;
				OutputClassUpdate(myrank,numprocs);
			}
			
		}
#ifdef USE_MPI
		mympi::mympi_location_gather(myrank,numprocs);
#endif//#ifdef USE_MPI
//			if( step % PVWin_OUT == 0 )		PV_out::PVWin_output_bin(myrank,numprocs);
		if(!myrank){
			if(step%100 == 0){fout_err<<"step"<<step<<endl;	}
			
			/*if( (step % AVR_TIME ) == 0 ){
				fout_err <<endl;
				fout_err << "avr_output" << endl;
				(void)AVR_output(AVR_TIME,step);
			}*/
			
			if(step%P_monitor == 0 )systemmonitor_ptcl(step);
			
			if( step % PVWin_OUT == 0 )		PV_out::PVWin_output_bin(myrank,numprocs);
			if( step % PTCL_OUT == 0 )		PTCL_output(t);
#ifdef WALL_BOUNDARY
			if( step % WALL_FRC_OUT == 0 )	frc_wall_output(myrank,numprocs);
//			if( step % WALL_FRC_OUT == 0 )	output_frc_wall(myrank,numprocs);
			OutputPlusEnd(myrank,numprocs);
#endif
			if( step % restart_out == 4000 )restart_binary_output(step);
			ptcl_fila_count();
			if( step % fila_out == 0 ){	arg_out();	}
			if( step % fila_arg == 0 ){	len_out();	}
			
			
			ptcl_number_output(myrank, numprocs);
			if(step%1000==0)output_polymerization_info(myrank, numprocs);
//cout<<step<<endl;
//			if((step%1000)==0)OutputForceCount(myrank, numprocs);	//なぜかうまくいかない．
//cout<<step<<endl;
			if((step % 10)==0)output_cc(myrank,numprocs);
			
			if((step % 200)==0)ptcl_wall_interaction(myrank,numprocs);
			if(step%5000==0)Flow_Output(myrank, numprocs);
			
			if((step % replenish_slice)==0)replenish_ptcl(myrank,numprocs);		//粒子浴
		}
#ifdef USE_MPI
		if((step % replenish_slice)==0)mympi::mympi_replenish_ptcl(myrank,numprocs);
#endif//#ifdef USE_MPI
		
	}
	if(!myrank){history("end");}
#ifdef USE_MPI
	mympi::mympi_finalize();
#endif //USE_MPI
	struct_finalize();

	//pdf_strain(8000);
	return 0;
}

#else //BROWNIAN


// int main(){
// 	dir_init();
// 	ss_fout<<"./output"<<setw(4)<<setfill('0')<<N<<"/fout_err.dat";
// 	fout_err.open(ss_fout.str().c_str(),ios::out);
// 	
// 	struct_init();
// 	file_init();
// 	history("start");
// 	ptcl_initial_condition();
// 	//PVWin_init();
// 	PVWin_init_bin();
// 	
// 	
// 	ene_format();
// 	Frc();
// 	
// 	PTCL_energy();
// 	PTCL_Location();
// 	
// 	//t += PTCL_dt;
// 	
// 	for(step =STEP_START ; step<STEP_END ; step++){
// 		for(int step_mpc=0 ; step_mpc<INTVL_step ; step_mpc++){
// 			for(int i=0 ; i<(int)(INTVL_COL) ; i++){
// 				ene_format();
// 				Frc();
// 				PTCL_energy();
// 				PTCL_Velocity();
// 				PTCL_Location();
// 				
// 				if(step>4000){			//+N*200){
// #ifdef POLYMERIZATION
// 					polymerization();
// #endif //POLYMERIZATION
// #ifdef DEPOLYMERIZATION
// 					depolymerization();
// #endif //DEPOLYMERIZATION
// #ifdef SEVERING
// 					severing();
// #endif //SEVERING
// 					//bind_cal();
// 				}
// 				p_level_cal();
// 				t += PTCL_dt;
// 				
// 			}
// 			
// 		}
// 		
// 		/*if( (step % AVR_TIME ) == 0 ){
// 			fout_err <<endl;
// 			fout_err << "avr_output" << endl;
// 			(void)AVR_output(AVR_TIME,step);
// 		}*/
// 		
// 		if(step%P_monitor == 0 )systemmonitor_ptcl(step);
// 		
// 		if(step%100 == 0){
// 			fout_err<<"step"<<step<<endl;
// 	//		//fout_err.close();
// 	//		//fout_err.open(ss_fout.str().c_str(),ios::out | ios::app);
// 	//		
// 		}
// 		//if( step % PVWin_OUT == 0 )		PVWin_output();
// 		if( step % PVWin_OUT == 0 )		PVWin_output_bin();
// 		if( step % PTCL_OUT == 0 )		PTCL_output(t);
// #ifdef WALL_BOUNDARY
// 		if( step % WALL_FRC_OUT == 0 )	frc_wall_output(step);
// #endif
// 		if( step % restart_out == 4000 )restart_binary_output(step);
// 		ptcl_fila_count();
// 		if( step % fila_out == 0 ){	arg_out();	}
// 		if( step % fila_arg == 0 ){	len_out();	}
// 		if( step % 10000 == 0 ){	numbering_debug();	}
// 		
// #if RESTART==1
// 		if(step ==STEP_START){
// 			restart_binary_input(RESTART_NUMBER);
// 			file_init();
// //			PVWin_init_restart_bin();
// 			PVWin_init_bin();
// 			PVWin_output_bin();
// 			systemmonitor_ptcl(step);
// 			arg_out();	len_out();
// 		}
// #endif
// #if 0
// 		{
// 			if(no_fila==0){
// 				static int check_step=0;
// 				check_step++;
// 				if(check_step>4000)exit(0);
// 			}
// 		}
// #endif
// 		
// 	}
// 	history("end");
// 	
// 	
// 	//pdf_strain(8000);
// 	return 0;
// }
#endif //BROWNIAN
#endif //NOT_SRD


//void list_check(){
//	for(int x=calcArea.x_pre;x<calcArea.x_bck;x++){
//	for(int y=calcArea.y_pre;y<calcArea.y_bck;y++){
//	for(int z=calcArea.z_pre;z<calcArea.z_bck;z++){
//		_list* li=list+(z + y*list_Z + x*(list_Y*list_Z));
//		_ptcl* pa;
//		for(int j=0;pa=li->ptcl_index[j];j++){
//			if(li!=pa->cnf->li){
//				cout<<x<<'\t'<<y<<'\t'<<z<<endl;
//				cout<<t<<" "<<j<<" "<<pa->index<<" "<<pa->cnf->feature<<" "<<pa->cnf->li->myrank<<" "<<pa->cnf->list_num<<endl;
//				exit(0);
//			}
//			if(j!=pa->cnf->list_num){
//				cout<<x<<'\t'<<y<<'\t'<<z<<endl;
//				cout<<t<<" "<<j<<" "<<pa->index<<" "<<pa->cnf->feature<<" "<<pa->cnf->li->myrank<<" "<<pa->cnf->list_num<<endl;
//				exit(0);
//			}
//			if(j>=li->num){
//				cout<<x<<'\t'<<y<<'\t'<<z<<endl;
//				cout<<t<<" "<<j<<" "<<pa->index<<" "<<pa->cnf->feature<<" "<<pa->cnf->li->myrank<<" "<<pa->cnf->list_num<<endl;
//				exit(0);
//			}
//		}
//	}}}
//}






#ifdef SRD
// int main(){
// 	dir_init();
// 	ss_fout<<"./output"<<setw(4)<<setfill('0')<<N<<"/fout_err.dat";
// 	fout_err.open(ss_fout.str().c_str(),ios::out);
// 	
// 	struct_init();
// 	file_init();
// 	history("start");
// 	
// 	MPC_Initialize();
// 	ptcl_initial_condition();
// 	//PVWin_init();
// 	PVWin_init_bin();
// 	
// 	
// 	ene_format();
// 	Frc();
// 	
// 	PTCL_energy();
// 	PTCL_Location();
// 	
// 	list_cal();
// 	//t += PTCL_dt;
// 	
// 	for(step =STEP_START ; step<STEP_END ; step++){
// 		for(int step_mpc=0 ; step_mpc<INTVL_step ; step_mpc++){
// 			for(int i=0 ; i<(int)(INTVL_COL) ; i++){
// 				ene_format();
// 				Frc();
// 				//bind_frc();
// 				PTCL_energy();
// 				PTCL_Velocity();
// 				PTCL_Location();
// 				list_cal();
// 				
// 				if(step>4000){			//+N*200){
// 					polymerization();
// 					depolymerization();
// 					severing();
// 					//bind_cal();
// 				}
// 				p_level_cal();
// 				t += PTCL_dt;
// 			}
// 			
// 			Propagation();
// 			CellCalculation();
// 			Collision();
// 		}
// 		
// 		//Mean();
// 		//if( (step % MPC_OUTPUT_TIME ) == 0 ){//mpc_output(step);}
// 		/*if( (step % AVR_TIME ) == 0 ){
// 			fout_err <<endl;
// 			fout_err << "avr_output" << endl;
// 			(void)AVR_output(AVR_TIME,step);
// 		}*/
// 		
// 		if(step%P_monitor == 0 )systemmonitor_ptcl(step);
// 		
// 		if(step%100 == 0){
// 			fout_err<<"step"<<step<<endl;
// 	//		//fout_err.close();
// 	//		//fout_err.open(ss_fout.str().c_str(),ios::out | ios::app);
// 	//		
// 		}
// 		if( step % STEP_MONITOR == 0 ) (void)SystemMonitor(step);
// 		//if( step % PVWin_OUT == 0 )		PVWin_output();
// 		if( step % PVWin_OUT == 0 )		PVWin_output_bin();
// 		if( step % PTCL_OUT == 0 )		PTCL_output(t);
// #ifdef WALL_BOUNDARY
// 		if( step % WALL_FRC_OUT == 0 )	frc_wall_output(step);
// #endif
// 		if( step % restart_out == 4000 )restart_binary_output(step);
// 		ptcl_fila_count();
// 		if( step % fila_out == 0 ){	arg_out();	}
// 		if( step % fila_arg == 0 ){	len_out();	}
// 		if( step % 10000 == 0 ){	numbering_debug();	}
// 		
// #if RESTART==1
// 		if(step ==STEP_START){
// 			restart_binary_input(RESTART_NUMBER);
// 			file_init();
// //			PVWin_init_restart_bin();
// 			PVWin_init_bin();
// 			PVWin_output_bin();
// 			systemmonitor_ptcl(step);
// 			arg_out();	len_out();
// 		}
// #endif
// #if 0
// 		{
// 			if(no_fila==0){
// 				static int check_step=0;
// 				check_step++;
// 				if(check_step>4000)exit(0);
// 			}
// 		}
// #endif
// 		
// 	}
// 	//{
// 	//	stringstream ss;
// 	//	ss<<"./output"<<N<<"/check_sev.dat";
// 	//	ofstream fout(ss.str().c_str());
// 	//	for(int i=0 ; i<Number_Ptcl ; i++){
// 	//		if(i%16==0){	fout<<endl;	}
// 	//		fout<<' '<<setw(8)<<check_sev[i].t;
// 	//	}
// 	//	fout<<endl;
// 	//	for(int i=0 ; i<Number_Ptcl ; i++){
// 	//		if(i%16==0){	fout<<endl;	}
// 	//		fout<<' '<<setw(4)<<check_sev[i].n;
// 	//	}
// 	//	fout<<endl;
// 	//	for(int i=0 ; i<Number_Ptcl ; i++){
// 	//		if(i%16==0){	fout<<endl;	}
// 	//		fout<<' '<<setw(4)<<check_sev[i].the/PI*180.0;
// 	//	}
// 	//}
// //	{
// //		stringstream ss;
// //		ss<<"./output"<<N<<"/pdf_sev_num.dat";
// //		ofstream fout(ss.str().c_str());
// //		unsigned int counter_his=0;
// //		for(int i=0 ; i<16 ; i++){
// //			counter_his += his_sev_num[i];
// //		}
// //		for(int i=0 ; i<16 ; i++){
// //			fout<<setw(3)<<i;
// //			fout<<'\t'<<setw(9)<<(double)his_sev_num[i]/(double)counter_his
// //				<<'\t'<<his_sev_num[i]<<'\t'<<counter_his<<endl;
// //		}
// //		fout.close();
// //		fout.clear();
// //		ss.str("");
// //		ss<<"./output"<<N<<"/pdf_sev_the.dat";
// //		fout.open(ss.str().c_str());
// //		for(int i=0 ; i<180 ; i++){
// //			fout<<setw(4)<<i;
// //			fout<<'\t'<<setw(9)<<(double)his_sev_the[i]/(double)counter_his<<endl;
// //		}
// //		fout.close();
// //		
// //	}
// 	//fout_err_sub("end.\n",1);
// 	history("end");
// 	
// 	
// 	//pdf_strain(8000);
// 	return 0;
// }
#endif //SRD
