/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<sys/stat.h>
#include <ctime>
#include<cstring>
#include "mymath.h"
#include "_system.h"
#include "_variable.h"
#include "_periodic.h"
#include "_arp2_3.h"
#include "_wall_boundary.h"


using namespace std;


double time1;
double time2;

std::ostream& OUTPUT_N(std::ostream& stream){//この関数はいけてない
	return stream<<"output"<<std::setw(4)<<std::setfill('0')<<N<<"/";
}
string OUTPUT_N(){
	stringstream ss;
	ss.str("");
	ss<<"output"<<std::setw(4)<<std::setfill('0')<<N<<"/";
	return ss.str();
}

inline bool MyMkdir(const std::string& ss){//mkdir
	ifstream fin(ss.c_str());
	if(!fin.is_open()){
		mkdir(ss.c_str(),0777);//0777はパーミション．8進数
	}else return false;
	fin.close();
	return true;
}
inline bool MyCopy(const std::string& ss_fromFile ,const std::string& ss_toFile){
	std::ifstream fromFile(ss_fromFile.c_str());
	std::ifstream toFileCheck(ss_toFile.c_str());
	std::ofstream toFile;
	if(fromFile.is_open() &&!(toFileCheck.is_open())){
		toFile.open(ss_toFile.c_str());
		if(fromFile.is_open() && toFile.is_open()){	toFile << fromFile.rdbuf();}
		else return false;
		toFile.close();
	}
	fromFile.close();
	toFileCheck.close();
	return true;
}
void rename_restart(std::string& ss){//fail name の頭に"_"を入れる
	string::size_type index = ss.find_last_of('\\');//カレントディレクトリの取得
	if(index == ss.npos)index=0;//発見されなければ
	else index++;
	string ss_tmp("_");
	ss.insert(index,ss_tmp);
}



//Read setting file
void InputInitefile(const int& argc,char**& argv,const int myrank,const int numprocs){
//	if(argc<2){for(int i(0);i<argc;i++)cout<<i<<" "<<argv[i]<<endl;exit(-1);}
	std::ifstream fin;
	if(argc>=2)fin.open(argv[1]);//settingファイル指定するならCommandLineから．
	if(!fin.is_open()){
		if(!myrank && (argc>=2))cerr<< "can't open CommandLine \""<<argv[1]<<'\"'<<endl;
		if(!myrank)cerr<< "try to Open \"system_init.dat\"."<<endl;//system_init.datを読みに．
		fin.clear();
		fin.open("system_init.dat");
		if(!fin.is_open()){
			if(!myrank)cerr<< "Calc under the Default Settings"<<endl;
			init_genrand(myrank);//乱数種
			NX=numprocs;
			NY=1;
			NZ=1;
			return ;
		}else {if(!myrank)cerr <<"Load Settings from \"system_init.dat\"."<<endl;}
	}
	string buf;
	string buf_sub;
	int Rand_Seed(0);
	int Rand_Seed_Switch(1);
	int line_counter(0);
	while(getline(fin,buf)){//1行読み込み
		std::stringstream ss;	ss.str("");//初期化
		ss<<buf;
//cout<<buf<<endl;
		ss>>buf_sub;//タブまたはスペース区切りの一つ目の文字
		if (buf_sub.empty() || buf.empty()){}
		else if (!buf_sub.compare(0,1,"#")){//commentouts
//			cout<<buf_sub<<endl;
		}else if(!buf_sub.compare(0,2,"//")){}//commentouts
		else if(buf_sub=="Output_Number" ){ss>>N;}				//出力ディレクトリナンバー
		else if(buf_sub=="Load"          ){ss>>WALL::jijuu;}	//荷重
		else if(buf_sub=="Pro_Poly"      ){ss>>probability_porimerize_f_p;}	//十号の結合確率
		else if(buf_sub=="Actin_Friction"){ss>>Ptcl_Friction;}				//抵抗係数
#ifdef WALL_BOUNDARY
		else if(buf_sub=="Wall_Friction" ){ss>>WALL::Wall_Friction;}
#endif// WALL_BOUNDARY
		else if(buf_sub=="Rand_Seed"     ){ss>>Rand_Seed;}
		else if(buf_sub=="Rand_Seed_Switch"){ss>>Rand_Seed_Switch;}			//乱数をNで変えるのかどうか
		else if(buf_sub=="Restart_Switch"){ss>>RESTART;}					//リスタートするかどうか
		else if(buf_sub=="Restart_Step"  ){ss>>RESTART_NUMBER;}
		else if(buf_sub=="Start_Step"    ){ss>>STEP_START;}
		else if(buf_sub=="End_Step"      ){ss>>STEP_END;}
		
		else if(buf_sub=="NX"){				ss>> NX;}						//X方向分割数
		else if(buf_sub=="NY"){				ss>> NY;}
		else if(buf_sub=="NZ"){				ss>> NZ;}
//		else if (buf_sub== ""){	ss>> ;}
		else{//defaults
			cout<<"line "<<line_counter<<" \""<<buf<<"\" is undef "<<endl;	//エラー
		}
		line_counter++;
		{//初期化
			buf.clear();
			buf_sub.clear();
		}
	}
	Ptcl_Rot_Friction=(Ptcl_Friction*4.0/3.0*sigma_LJ*sigma_LJ*0.25);		//回転の抵抗係数
	init_genrand(myrank + Rand_Seed + Rand_Seed_Switch*N);					//乱数種
	if((NX * NY * NZ)!=numprocs){//分割数とノード数が合わないとき
		NX=numprocs;//x軸のみで切る
		NY=1;
		NZ=1;
	}
//	cout<<"Output_Number "<<N<<endl
//	<<"Load "<<WALL::jijuu<<endl
//	<<"Pro_Poly "<<probability_porimerize_f_p<<endl
//	<<"Actin_Friction"<<Ptcl_Friction<<endl
//	<<"Wall_Friction "<<WALL::Wall_Friction<<endl
//	<<"Rand_Seed "<<Rand_Seed<<endl
//	<<"Rand_Seed_Switch "<<Rand_Seed_Switch<<endl
//	<<"Restart_Switch "<<RESTART<<endl
//	<<"Restart_Step "<<RESTART_NUMBER<<endl
//	<<"Start_Step "<<STEP_START<<endl
//	<<"End_Step "<<STEP_END<<endl
//	<<"NX "<<NX<<endl
//	<<"NY "<<NY<<endl
//	<<"NZ "<<NZ<<endl;
	
}



//const string Output_N(OUTPUT_N());
void dir_init(){
	MyMkdir(OUTPUT_N());
	MyMkdir("actinpv/");
	MyCopy("/home/deji/pv/ball.bmp" ,"actinpv/ball.bmp");
	MyCopy("/home/deji/pv/bondc.dat","actinpv/bondc.dat");
	
#ifdef INITIAL_CONDITION_DEF
	MyMkdir("initial_condition/");
#endif //INITIAL_CONDITION_DEF
	return;
}

void file_init(){
	stringstream ss;
	ofstream fout;
//	ss.str("");
//	ss<<OUTPUT_N<<"system_ptcl.dat";
//	fout.open(ss.str().c_str(),ios::out);
//	fout<<"step number_of_filament number_of_monomer atp adp"<<endl;
//	fout.close();
	
//	ss.str("");
//	ss<<OUTPUT_N<<"filament_len_ave_var.dat";
//	fout.open(ss.str().c_str(),ios::out);
//	fout<<"step	ave_len	var_len"<<endl;
//	fout.close();
//	
//	
//	ss.str("");
//	ss<<OUTPUT_N<<"filament_len_ave_var_ver1.dat";
//	fout.open(ss.str().c_str(),ios::out);
//	fout<<"step\tave_len\tvar_len\tave_len2\tvar_len2"<<endl;
//	fout.close();
	
#ifdef SRD
	double mu_kin;
	double mu_col;
	double gam = number_dens_cell;
	mu_kin = SYS_kT * MPC_dt * rho_mpc / (2.0*MPCMASS);
	mu_kin *= (( 5.0*gam - (gam - 1.0 + exp(-1.0*gam))*(2.0-cos(col_angl)-cos(2.0*col_angl)) )/( (gam-1.0+exp(-1.0*gam))*(2.0-cos(col_angl)-cos(2.0*col_angl)) ));
	mu_col = MPCMASS / (18.0*col_cell*MPC_dt)*( (gam-1.0+exp(-1.0*gam))*(1.0-cos(col_angl)) );
	
	double xi = 0;
	double xi_sub = 0;
	double lambda = number_dens_cell;
	for(int i=1 ; i<Number_mpc+1 ; i++){
		xi_sub=(exp(-1.0*lambda))*(MPCMASS*i/(PTCL_MASS+MPCMASS*i)*(1.0-cos(col_angl)));
		for(int j=i ; j>0 ; j--){
			xi_sub = xi_sub*lambda/j;
			if((int)(xi_sub*lambda/j*1.0E300)==0)break;
		}
		if(isinf(xi_sub)!=0 || isnan(xi_sub)!=0){exit(0);}
		xi+=xi_sub;
		if((int)(xi_sub*1.0E300)==0)break;
	}
	xi=xi*2.0*PTCL_MASS/(3.0*MPC_dt);
//	cout<<xi<<endl;
//	cout<<SYS_kT/xi<<endl;
	
	ss.str("");
	ss<<OUTPUT_N<<"viscosity.dat";
	fout.open(ss.str().c_str(),ios::out);
	fout<<"xi = "<<xi<<" T = "<<SYS_kT<<" D = "<<SYS_kT/xi<<endl;
	fout<<"no ptcl viscosity :"<<endl<<"  mu_kin ="<<mu_kin
		<<" ;"<<endl<<"  mu_col = "<<mu_col
		<<" ;"<<endl<<"  shear viscosity = "<<mu_kin+mu_col<<endl;
	fout<<"step  viscosity  T  xi?  D?  dt?"<<endl;
	fout.close();
#endif //SRD
	
	
	ss.str("");
	ss<<OUTPUT_N<<"sys.dat";
	fout.open(ss.str().c_str(),ios::out);
	fout
#ifdef SRD
		<<" SVX = "<<2.0*SVX
#endif
		<<" Number_Ptcl = "<<Number_Ptcl
#ifdef SRD
		<<" Number_mpc = "<<Number_mpc<<endl
#endif
		<<" L*L*L = "<<SYS_X<<'*'<<SYS_Y<<'*'<<SYS_Z<<endl
		<<" actin "<<endl/*<<" m = "<<PTCL_MASS*/<<" ; sigma = "<<sigma_LJ<<" ; dt = "<<PTCL_dt<<endl
#ifdef SRD
		<<" fluid "<<endl<<" m = "<<MPCMASS<<" ; rho = "<<number_dens<<" ; dt = "<<MPC_dt<<endl
		<<" col_angl = "<<col_angl*180.0/PI<<endl
		<<" col_cell = "<<col_cell<<endl<<endl
#endif //SRD
		<<" por_angl = "<<180.0-por_angl*180.0/PI<<endl
		<<" probability_porimerize_f_p = "<<probability_porimerize_f_p<<endl
		<<" probability_porimerize_p_p = "<<probability_porimerize_p_p<<endl
		<<" probability_deporimerize = "<<probability_deporimerize<<endl
		<<" pro_sev = "<<pro_sev<<endl
		<<" len_eq = "<<l_eq<<endl
#ifdef SRD
		<<" no (actin)ptcl viscos :\n  mu_kin ="<<mu_kin
		<<" ;\n  mu_col = "<<mu_col
		<<" ;\n  shear viscosity = "<<mu_kin+mu_col
#endif //SRD
		<<endl;
#ifdef SRD
	if((!isinf(xi)) && (!isnan(xi)) && (((int)(xi*1.0E7))!=0) ){
		fout<<" xi = "<<xi<<" T = "<<SYS_kT<<" D = "<<SYS_kT/xi<<endl;
	}
	else{	fout<<" xi = 0?or inf or nan. xi = "<<xi<<" T = "<<SYS_kT<<endl;}
#endif //SRD
	fout.close();
	
}

void history(string c){//ベンチ用　始まりと終わりの時間を出力
	std::ofstream fout("./#history.dat",ios::out | ios::app);
	time_t t=time(NULL);
	struct tm *gtm = localtime(&t);
	
	fout<<setw(4)<<gtm->tm_year + 1900<<'/'
		<<setw(2)<<setfill('0')<<gtm->tm_mon + 1<<'/'
		<<setw(2)<<setfill('0')<<gtm->tm_mday<<' '
		<<setw(2)<<setfill('0')<<gtm->tm_hour<<':'
		<<setw(2)<<setfill('0')<<gtm->tm_min<<':'
		<<setw(2)<<setfill('0')<<gtm->tm_sec
		<<" : N = "<<N
//		<<" : svx = "<<SVX*2.0<<"\n	: γ = "<<(double)SVX*2.0/SYS_Z
//		<<" : α = "<<col_angl*180.0/PI
		<<" n_p "<<Number_Ptcl<<" n_a "<<Number_Arp2_3<<"    "<<c<<endl;
		
	fout.close();
}
#ifdef WALL_BOUNDARY
void frc_wall_output(const int myrank ,const int numprocs){
	static string ss(OUTPUT_N()+"frc_wall.dat");
	static bool flag=true;
	if(flag){//header 書き出し
		std::ofstream fout;
if(RESTART)fout.open(ss.c_str(),ios::out | ios::app);
else fout.open(ss.c_str(),ios::out);
			if(!fout.is_open()){//頭に"_"を挿入します
				static int count_err(0);
				flag=true;
				if(count_err<10)count_err++;
				else{fout_err<<"frc_wall_output"<<endl;	exit(0);}
				rename_restart(ss);
				frc_wall_output(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
		fout<<"wall_force is "<<WALL::jijuu<<" friction is "<<WALL::Wall_Friction<<endl;
		fout<<"step\tt\twallposion\t\tvel\t\t"
			<<"NET_force(without Random)\t\t"
			<<"RANDOM_force\t\t"
			<<"force_from_LJpotential\t\t"
			<<"force_from_spring_potential\t\t\t\t\t\t"
			<<"force_from_monomers\t\t"
			<<"force_from_filaments\t\t"
			<<"counter"
			<<endl;
		fout<<"\t\t\t\tmean\tvariance\tmean\tvariance\t"
			<<"mean\tvariance\t"
			<<"mean\tvariance\t"
			<<"mean\tvariance\t"
			<<"mean\t\t\tvariance\t\t\t"
			<<"mean\tvariance\t"
			<<"mean\tvariance\t"
			<<endl;
		flag=false;
		fout.close();
	}
	static const int NVW(1000);
	static int counter_VW(0);
	static _output_variable_wall tmp_VW[NVW];
	static int tmp_VW_step[NVW];
	static double tmp_VW_t[NVW];
	{
		tmp_VW[counter_VW]=output_variable_wall;
		tmp_VW_step[counter_VW]=step;
		tmp_VW_t[counter_VW]=t;
		counter_VW++;
		output_variable_wall.clear();
	}
	if (counter_VW>=NVW || step>STEP_END-10){
		ofstream fout;
		fout.open(ss.c_str(),ios::out | ios::app);
			if(!fout.is_open()){//頭に"_"を挿入します
				static int count_err(0);
				flag=true;
				if(count_err<10)count_err++;else{fout_err<<"frc_wall_output"<<endl;	exit(0);}
				rename_restart(ss);
				frc_wall_output(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
		
		for(int i(0);i<counter_VW && i<NVW;++i){
			fout<<tmp_VW_step[i]<<"\t"<<tmp_VW_t[i]<<"\t"<<tmp_VW[i]<<endl;
		}
		fout.close();
		fout.clear();
		counter_VW=0;
	}
	
//#ifdef SRD
//	
//	double viscosity;
//	if(((int)(SVX*1.0E5)) != 0){viscosity=f_mp_g.x/(-1.0*SVX*2.0/SYS_Z);}
//	else{viscosity=f_mp_g.x/(-1.0*1E-3);}
//	ss.str("");
//	ss<<OUTPUT_N<<"viscosity.dat";
//	fout.open(ss.str().c_str(),ios::out | ios::app);
//	fout<<step<<"  "<<viscosity<<"  "<<T<<"  "<<6.0*PI*sigma_LJ*viscosity<<"  "
//		<<T/(6.0*PI*sigma_LJ*viscosity)<<"  "<<sigma_LJ*sigma_LJ/(T/(6.0*PI*sigma_LJ*viscosity))<<endl;
//	fout.close();
//	fout.clear();
//#endif
	
}
void OutputPlusEnd(const int myrank,const int numprocs){
	static string ss(OUTPUT_N()+"plus_end_posion.dat");
	static bool flag=true;
	if(flag){
//		std::ofstream fout(ss.c_str(),ios::out);
		ofstream fout;//(ss.c_str(),ios::out | ios::app);
if(RESTART)fout.open(ss.c_str(),ios::out | ios::app);
else fout.open(ss.c_str(),ios::out);
			if(!fout.is_open()){
				static int count_err(0);flag=true;
				if(count_err<10)count_err++;else{fout_err<<"OutputPlusEnd"<<endl;	exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				OutputPlusEnd(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
		fout<<"wall_force is "<<WALL::jijuu<<" friction is "<<WALL::Wall_Friction<<endl;
		fout<<"step\tt\tPlusEndPosion\t\tWallPosition\t\t"
			<<"PlusEndPosion_Sub\t\t"
			<<"counter\tcounterforwall"
			<<endl;
		fout<<"\t\tmean\tvariance\tmean\tvariance\t"
			<<"mean\tvariance\t"
			<<endl;
		flag=false;
		fout.close();
	}
	static const int NPEP(1000);
	static int counter_PEP(0);
	static _output_PlusEndPosition tmp_PEP[NPEP];
	static int tmp_PEP_step[NPEP];
	static double tmp_PEP_t[NPEP];
	{
		tmp_PEP[counter_PEP]=output_PlusEndPosition;
		tmp_PEP_step[counter_PEP]=step;
		tmp_PEP_t[counter_PEP]=t;
		counter_PEP++;
		output_PlusEndPosition.clear();
	}
	if (counter_PEP>=NPEP || step>STEP_END-10){
		ofstream fout;
		fout.open(ss.c_str(),ios::out | ios::app);
			if(!fout.is_open()){
				static int count_err(0);flag=true;
				if(count_err<10)count_err++;else{fout_err<<"OutputPlusEnd"<<endl;	exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				OutputPlusEnd(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
		
		for(int i(0);i<counter_PEP && i<NPEP;++i){
			fout<<tmp_PEP_step[i]<<"\t"<<tmp_PEP_t[i]<<"\t"<<tmp_PEP[i]<<endl;
		}
		fout.close();
		fout.clear();
		counter_PEP=0;
	}
}
void OutputForceCount(const int myrank,const int){
	{
		static bool flag=true;
		bool flag_force=false;
		if(step==STEP_END-1)flag_force=true;
//cout<<flag<<'\t'<<flag_force<<endl;
		for(int i(0);i<plus_end->num;++i){
			stringstream ss1;
			ss1.str("");
			ss1<<OUTPUT_N<<i<<"ForceCount_Force.dat";
			stringstream ss2;
			ss2.str("");
			ss2<<OUTPUT_N<<i<<"ForceCount_Time.dat";
			Force_Count[i].Output(ss1.str(),ss2.str(),flag,flag_force);
		}
		if(flag){flag=false;}
	}
}
#endif //WALL_BOUNDARY
/*void ptcl_debug(double t){
	{
		static bool flag(true);
		if(flag){
			std::stringstream ss;
			ss.str("");
			{
				ss<<OUTPUT_N<<"ptcl/";
				std::ifstream fin(ss.str().c_str());
				if(!fin.is_open()){	mkdir(ss.str().c_str(),0777);	}
				fin.close();
			}
			flag=false;
		}
	}
	
	stringstream ss;
	ss.str("");
	ss<<OUTPUT_N<<"ptcl/output"<<setw(8)<<setfill('0')<<(int)(t*100)<<".dat";
	ofstream fout(ss.str().c_str(),ios::out );
	for( int i = 0; i < Number_Ptcl; i++){
		_ptcl* p = &ptcl[i];
		fout<<t<<'\t'<<p->loc.x<<'\t'<<p->loc.y<<'\t'<<p->loc.z<<'\t'
			<<p->vel.x<<'\t'<<p->vel.y<<'\t'<<p->vel.z<<'\t'
			<<p->frc_new.x<<'\t'<<p->frc_new.y<<'\t'<<p->frc_new.z<<'\t'
			<<p->feature<<'\t'<<p->num<<endl;
	}
	
	fout.close();
}*/
void systemmonitor_ptcl(int step){
// 	string ss(Output_N+"system_ptcl.dat");
//	static bool flag(true);
//	if(flag){
//		std::ofstream fout.open(ss.c_str(),ios::out);
//		fout<<"step number_of_filament number_of_monomer atp adp"<<endl;
//		fout.close();
//		flag=false;
//	}
//	ofstream fout(ss.c_str(),ios::out | ios::app);
// 	int counter_p=0;
// // 	for(int i=0 ; numbering[0].ptcl_n[i] ; i++){
// // 		if(numbering[0].ptcl_n[i]->p_level==1){counter_p++;}
// // 	}
// 	for(_ptcl* p=ptcl;p<(ptcl+Number_Ptcl);++p){
// 		if(!(p->feature) && (p->p_level==1))counter_p++;
// 	}
// 	int list_bind_fila[no_fila+1];
// 	memset(list_bind_fila,0,(no_fila+1)*sizeof(int));
// 	for(int i=0 ; i<bind_size ; i++){
// 		if(bind_wall[i].index){
// 			list_bind_fila[bind_wall[i].index->feature] += 1;
// 		}
// 	}
// 	int count_p_bind=0;
// 	for(int i=1 ; i<=no_fila ; i++){
// 		if(list_bind_fila[i]){
// 			count_p_bind += numbering[i].num;
// 		}
// 	}
// 	
// 	
// 	fout<<step<<'\t'<<no_fila<<'\t'<<numbering[0].num<<'\t'<<counter_p
// 		<<'\t'<<numbering[0].num - counter_p<<'\t'
// 		<<counter_pol<<'\t'<<counter_depol<<'\t'<<counter_severing<<'\t'
// 		//<<(double)counter_severing/counter_severing_max<<'\t'
// 		//<<Number_Ptcl-numbering[0].num<<'\t'
// 		<<count_p_bind<<endl;
// 	fout.close();
}
void len_out(){
//	{
//		static bool flag(true);
//		if(flag){
//			std::stringstream ss;
//			ss.str("");
//			{
//				ss<<OUTPUT_N<<"counter/";
//				std::ifstream fin(ss.str().c_str());
//				if(!fin.is_open()){	mkdir(ss.str().c_str(),0777);	}
//				fin.close();
//			}
//			flag=false;
//		}
//	}
//	
// 	if(step % fila_arg == 0){
// 		double len_ave=((double)(Number_Ptcl - numbering[0].num))/((double)no_fila);
// 		double len_var=0;
// 		//if(no_fila==0){len_ave=0.0;}
// 		stringstream ss;
// 		ss.str("");
// 		ss<<OUTPUT_N<<"filament_len_ave_var.dat";
// 		ofstream fout_c;
// 		if(no_fila!=0){
// 			fout_c.open(ss.str().c_str(),ios::out | ios::app);
// 			for(int i=1 ; numbering[i].num ; i++){
// 				len_var += (numbering[i].num - len_ave) * (numbering[i].num - len_ave);
// 			}
// 			len_var /= no_fila;
// 			fout_c<<step<<'\t'<<len_ave<<'\t'<<len_var<<endl;
// 			
// 			fout_c.close();
// 		}
// 	
// 		if(step % fila_arg1 == 0){
// 			//stringstream ss;
// 			ss.str("");
// 			ss<<OUTPUT_N<<"counter/output_filament"<<setw(8)<<setfill('0')<<step<<".dat";
// 			if(no_fila!=0){
// 				ofstream fout_c;
// 				fout_c.open(ss.str().c_str(),ios::out);
// 				fout_c<<"number_of_firament = "<<no_fila<<endl
// 					  <<"average_length = "<<len_ave<<endl//((double)(Number_Ptcl-numbering[0].num))/((double)no_fila)<<endl
// 					  <<"variance_length = "<<len_var<<endl
// 					  <<"filament_number number_of_firament_particle"<<endl;
// 				for(int i=1 ; i<no_fila+1 ; i++){fout_c<<i<<'\t'<<numbering[i].num<<endl;}
// 				fout_c.close();
// 			}
// 		}
// 	}
// 	if(step % fila_arg == 0){
// 		double len_ave=0.0;
// 		double len_var=0.0;
// 		double len_ave2=0.0;
// 		double len_var2=0.0;
// 		stringstream ss;
// 		ss.str("");
// 		ss<<OUTPUT_N<<"filament_len_ave_var_ver1.dat";
// 		ofstream fout_c(ss.str().c_str(),ios::out | ios::app);;
// 		for(int i=1 ; numbering[i].num ; i++){
// 			double len_sub=0.0;
// 			if(numbering[i].num!=1){
// 				for(int j=0 ; numbering[i].ptcl_n[j+1] ; j++){
// 					_vec<double> del = numbering[i].ptcl_n[j]->loc - numbering[i].ptcl_n[j+1]->loc;
// 					distance_boundary(del);
// 					len_sub+=del.norm();
// 				}
// 			}
// 			len_sub+=sigma_LJ;
// 			len_ave+=len_sub;
// 			len_var+=len_sub*len_sub;
// 		}
// 		len_ave2=len_ave;	len_var2=len_var;
// 		len_ave+=(double)sigma_LJ*(double)numbering[0].num;
// 		len_var+=(double)sigma_LJ*(double)sigma_LJ*(double)numbering[0].num;
// 		if(no_fila!=0){
// 			len_ave2/=(double)no_fila;
// 			len_var2=len_var2/(double)no_fila - len_ave2*len_ave2;
// 		}
// 		else{len_ave2=0;	len_var2=0;}
// 		len_ave=len_ave/((double)(no_fila + numbering[0].num));
// 		len_var=len_var/((double)(no_fila + numbering[0].num)) - len_ave*len_ave;
// 		fout_c<<step<<'\t'<<len_ave<<'\t'<<len_var<<'\t'<<len_ave2<<'\t'<<len_var2<<endl;
// 		
// 		fout_c.close();
// 		
// 	}
	return;
}
void arg_out(){
	static string ss(OUTPUT_N()+"filament_arg_ave_var.dat");
	static string ss1(OUTPUT_N()+"counter/");
	stringstream ss2;
	ss2.str("");
	
	{
		static bool flag(true);
		if(flag){
			MyMkdir(ss1);
			ofstream fout_arg;//(ss.c_str(),ios::out | ios::app);
if(RESTART){			fout_arg.open(ss.c_str(),ios::out | ios::app);
}else{					fout_arg.open(ss.c_str(),ios::out);
}			

			fout_arg.close();
			flag=false;
		}
	}
//	if(step % fila_arg == 0){
		ofstream fout_arg;
		fout_arg.open(ss.c_str(),ios::out | ios::app);
			if(!fout_arg.is_open()){
				static int count_err(0);
				if(count_err<10)count_err++;else{fout_err<<"arg_out"<<endl;	exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				arg_out();//再帰で呼べばいいんじゃね
				return;
			}
		ofstream fout_polar;
		if(step % fila_arg == 0){
			ss2.str("");
			ss2<<OUTPUT_N<<"counter/seg_polar"<<setw(8)<<setfill('0')<<step<<".dat";
			fout_polar.open(ss2.str().c_str(),ios::out);
		}
		
		double ave_the(0.0);// = 0.0;
		double var_the(0.0);// = 0.0;
		double ave_phi(0.0);// = 0.0;
		double var_phi(0.0);// = 0.0;
		//double ave_len(0.0);// = 0.0;
		//double var_len(0.0);// = 0.0;
		double ave_strain(0.0);// = 0.0;
		double var_strain(0.0);// = 0.0;
		
		double ave_lxlz = 0.0;
		double var_lxlz = 0.0;
		int count_seg = 0;
		double fila_the = 0.0;
		double fila_phi = 0.0;
		
		int i(0);
//		for(_ptcl** end=&plus_end->end[0];*end;++end){
		for(int i_end(0);i_end<plus_end->num;++i_end){
			_ptcl* end = plus_end->end(i_end);
			++i;
			_ptcl* p=end;
			_ptcl* q=end;
			int j(0);
			while(p=p->cnf->minus){
				
				_vec<double> vec = q->loc - p->loc;
				distance_boundary(vec);
				if(vec.x<0){vec = -1.0*vec;}
				double the = acos(vec.z/vec.norm());
				double phi = asin(vec.y/(sqrt(vec.x*vec.x + vec.y*vec.y)));
				
				ave_the+=the;
				var_the+=the*the;
				ave_phi+=phi;
				var_phi+=phi*phi;
				//ave_len+=vec.norm();
				//var_len+=vec.norm()*vec.norm();
				ave_strain+=vec.norm()/l_eq-1.0;
				var_strain+=(vec.norm()/l_eq-1.0)*(vec.norm()/l_eq-1.0);
				
				ave_lxlz += vec.x*vec.z;
				var_lxlz += (vec.x*vec.z) * (vec.x*vec.z);
				
				count_seg++;
				
				if(step % fila_arg == 0){
					fout_polar<<i<<'\t'<<j<<'\t'
						<<vec.x<<'\t'<<vec.y<<'\t'<<vec.z<<"\t\t"
						<<vec.norm()<<'\t'<<the<<'\t'<<phi<<'\t'
						<<vec.norm()/l_eq-1.0<<'\t'<<vec.x*vec.z<<'\t'
						<<(p->loc.z + q->loc.z)*0.5<<endl;
				}
				q=p;
				++j;
			}
			_vec<double> vec = end->loc - q->loc;
			distance_boundary(vec);
			if(vec.x<0){vec *= -1.0;}
			double the = acos(vec.z/vec.norm());
			double phi = asin(vec.y/(sqrt(vec.x*vec.x + vec.y*vec.y)));
			fila_the += the;
			fila_phi += phi;
		}
		no_fila=i;
		if(count_seg!=0){
			ave_the /= (double)count_seg;	ave_phi /= (double)count_seg;	/*ave_len /= (double)count_seg;*/	ave_lxlz /= (double)count_seg;
			ave_strain /= (double)count_seg;
			var_the /= (double)count_seg;	var_phi /= (double)count_seg;	/*var_len /= (double)count_seg;*/	var_lxlz /= (double)count_seg;
			var_strain /= (double)count_seg;
			var_the -= ave_the*ave_the;		var_phi -= ave_phi*ave_phi;		/*var_len -= ave_len*ave_len;*/	var_lxlz -= ave_lxlz*ave_lxlz;
			var_strain -= ave_strain*ave_strain;
			fila_the /= no_fila;	fila_phi /= no_fila;
			fout_arg<<step<<'\t'<<ave_the<<'\t'<<var_the
				<<"\t\t"<<ave_phi<<'\t'<<var_phi
				//<<'\t'<<ave_len<<'\t'<<var_len
				<<"\t\t"<<ave_strain<<'\t'<<var_strain
				<<"\t\t"<<ave_lxlz<<'\t'<<var_lxlz
				<<"\t\t"<<fila_the<<'\t'<<fila_phi<<endl;
				
		}
		fout_arg.close();
		fout_polar.close();
//	}
}
void ptcl_fila_count(){
//	{
//		static bool flag(true);
//		if(flag){
//			std::stringstream ss;
//			ss.str("");
//			{
//				ss<<OUTPUT_N<<"counter/";
//				std::ifstream fin(ss.str().c_str());
//				if(!fin.is_open()){	mkdir(ss.str().c_str(),0777);	}
//				fin.close();
//			}
//			flag=false;
//		}
//	}
// 	for(int i=1 ; numbering[i].num ; i++ ){
// 		counter_num[numbering[i].num] += 1;
// 	}
// 	if(step % fila_count_out == 0){
// 		stringstream ss;
// 		ss.str("");
// 		ss<<OUTPUT_N<<"counter/output_counter_filament"<<setw(8)<<setfill('0')<<step<<".dat";
// 		ofstream fout(ss.str().c_str(),ios::out);
// // 		char fna[100];
// // 		FILE *fptr;
// // 		sprintf(fna,"./output%.4d/counter/output_counter_filament%07d.dat",N,step);
// // 		fptr = fopen(fna,"w");
// 		//int counter_num[number_max];
// 		for(int i=0;i<500/*number_max*/;i++){
// // 			fprintf(fptr,"%d %f\n",i,(double)counter_num[i]/(double)fila_count_out);
// 			fout<<i<<"\t"<<(double)counter_num[i]/(double)fila_count_out<<endl;;
// 		}
// 		memset(counter_num,0,number_max * sizeof(int));
// 		
// // 		fclose(fptr);
// 		fout.close();
// 	}
	
}
void PTCL_output(double t){
	static bool flag(true);
	static string ss(OUTPUT_N()+"output_ptcl.dat");
	if(flag){
		ofstream fout;//(ss.c_str(),ios::out | ios::app);
if(RESTART)fout.open(ss.c_str(),ios::out | ios::app);
else fout.open(ss.c_str(),ios::out);
//		ofstream fout(ss.c_str(),ios::out);
			if(!fout.is_open()){
				static int count_err(0);flag=true;
				if(count_err<10)count_err++;else{fout_err<<"PTCL_output"<<endl;	exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				PTCL_output(t);//再帰で呼べばいいんじゃね
				return;
			}
		fout<<"step\tt\tene\t\tarp\t\twall\t\tspeed[ms/dt]"<<endl;
		fout<<"\t\tmean\tvariance\tmean\tvariance\t"
			<<"mean\tvariance\t"
			<<endl;
		fout.close();
		flag=false;
	}
	
	static const int NPO(10);
	static int counter_PO(0);
	static _PTCL_Ene_Output tmp_PO[NPO];
	static int tmp_PO_step[NPO];
	static double tmp_PO_t[NPO];
	{
		PTCL_Ene_Output.calc_speed();
		tmp_PO[counter_PO]=PTCL_Ene_Output;
		tmp_PO_step[counter_PO]=step;
		tmp_PO_t[counter_PO]=t;
		counter_PO++;
		PTCL_Ene_Output.clear();
	}
	if (counter_PO>=NPO || step>STEP_END-10){
		ofstream fout;
		fout.open(ss.c_str(),ios::out | ios::app);
			if(!fout.is_open()){
				static int count_err(0);flag=true;
				if(count_err<10)count_err++;else{fout_err<<"PTCL_output"<<endl;	exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				PTCL_output(t);//再帰で呼べばいいんじゃね
				return;
			}

		for(int i(0);i<counter_PO && i<NPO;++i){
			fout<<tmp_PO_step[i]<<"\t"<<tmp_PO_t[i]<<"\t"<<tmp_PO[i]<<endl;
		}
		fout.close();
		fout.clear();
		counter_PO=0;
	}
//	ofstream fout(ss.c_str(),ios::out|ios::app);
//	fout<<step<<"\t"<<t<<'\t'<<PTCL_Ene_Output<<endl;
//	PTCL_Ene_Output.clear();
//	fout.close();
	
}

void Flow_Output(const int myrank,const int){
	static bool flag(true);
	std::stringstream ss;
	if(flag){
		{
			ss.str("");
			ss<<OUTPUT_N()<<"flow/";
			MyMkdir(ss.str());
		}
		flag=false;
	}
	ss.str("");
	ss<<OUTPUT_N<<"flow/absolute"<<setw(8)<<setfill('0')<<step<<".plt";
	ofstream fout_abs(ss.str().c_str(),ios::out);
	ss.str("");
	ss<<OUTPUT_N<<"flow/relative"<<setw(8)<<setfill('0')<<step<<".plt";
	ofstream fout_rel(ss.str().c_str(),ios::out);
	fout_abs<<"variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"u_mono\",\"v_mono\",\"w_mono\",\"rho\",\"rho_mono\""<<endl
			<<"zone t=\"3D zone\" i="<<flow.size.x <<" j="<<flow.size.y<<" k="<<flow.size.z <<"  f=point"<<endl;
	fout_rel<<"variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"u_mono\",\"v_mono\",\"w_mono\",\"rho\",\"rho_mono\""<<endl
			<<"zone t=\"3D zone\" i="<<flow.size.x <<" j="<<flow.size.y<<" k="<<flow.size.z <<"  f=point"<<endl;
	for(int k(0);k<flow.size.z;k++){
	for(int j(0);j<flow.size.y;j++){
	for(int i(0);i<flow.size.x;i++){
		double dt=(t - flow.start_step)/PTCL_dt;
//		if(!dt)continue;
		int index(i + j*flow.size.x + k*flow.size.x*flow.size.y);
		flow.abs[index].vel/=dt;
		flow.abs_mono[index].vel/=dt;
		fout_abs<<i<<' '<<j<<' '<<k<<' '<<flow.abs[index].vel.x<<' '<<flow.abs[index].vel.y<<' '<<flow.abs[index].vel.z<<' '
			<<flow.abs_mono[index].vel.x<<' '<<flow.abs_mono[index].vel.y<<' '<<flow.abs_mono[index].vel.z<<' '<<flow.abs[index].count/dt<<' '<<flow.abs_mono[index].count/dt<<endl;
		
		flow.rel[index].vel/=dt;
		flow.rel_mono[index].vel/=dt;
		fout_rel<<i<<' '<<j<<' '<<k<<' '<<flow.rel[index].vel.x<<' '<<flow.rel[index].vel.y<<' '<<flow.rel[index].vel.z<<' '
			<<flow.rel_mono[index].vel.x<<' '<<flow.rel_mono[index].vel.y<<' '<<flow.rel_mono[index].vel.z<<' '<<flow.rel[index].count/dt<<' '<<flow.rel_mono[index].count/dt<<endl;
		
	}}}
	flow.clear();
	fout_rel.close();
	fout_abs.close();
	
}

struct _pointer_set{
public:
	_ptcl* p;
#ifdef ARP2_3
	_ARP2_3::_arp2_3* ar;
#endif
	_list* li;
};

void restart_binary_output(int restart_step){
	{
		static bool flag(true);
		if(flag){
			string ss(OUTPUT_N()+"restart/");
			MyMkdir(ss);
		}
	}
	
	_pointer_set pointset;
	pointset.p=&ptcl[0];
#ifdef ARP2_3
	pointset.ar=_ARP2_3::arp2_3;
#endif
	pointset.li=list;
	
	stringstream ss;
	ofstream fout;
	
#ifdef SRD
	ss.str("");
	ss<<OUTPUT_N<<"restart/"<<setw(10)<<setfill('0')<<restart_step<<"_mpc.txt";
	fout.open(ss.str().c_str(),ios::out | ios::binary);
	fout.write(( char * ) mpc, Number_mpc*sizeof( _mpc ) );
	fout.close();
#endif
	int np=Number_Ptcl;
	ss.str("");
	ss<<OUTPUT_N<<"restart/"<<setw(10)<<setfill('0')<<restart_step<<"_ptcl.txt";
	fout.open(ss.str().c_str(),ios::out | ios::binary);
	fout.write(( char * )&np ,sizeof( int ) );
	fout.write(( char * )&Active_Ptcl ,sizeof( int ) );
	fout.write(( char * )&pointset ,sizeof( _pointer_set ) );
	fout.write(( char * )ptcl.pt ,Number_Ptcl*sizeof( _ptcl ) );
	fout.write(( char * )ptcl.con,Number_Ptcl*sizeof( _config ) );
	fout.write(( char * )plus_end  ,sizeof( _Actin::_Plus_End ) );
	fout.write(( char * )minus_end ,sizeof( _Actin::_Minus_End ) );
	fout.close();
	
#ifdef WALL_BOUNDARY
	ss.str("");
	ss<<OUTPUT_N<<"restart/"<<setw(10)<<setfill('0')<<restart_step<<"_bind_wall.txt";
	fout.open(ss.str().c_str(),ios::out | ios::binary);
	fout.write(( char * ) WALL::bind_wall,WALL::bind_size*sizeof( _bind_wall ) );
	fout.write(( char * ) WALL::WALL_Z ,2*sizeof( _WALL_Z ) );
	fout.close();
#endif //#ifdef WALL_BOUNDARY
#ifdef ARP2_3
	ss.str("");
	ss<<OUTPUT_N<<"restart/"<<setw(10)<<setfill('0')<<restart_step<<"_arp2_3.txt";
	fout.open(ss.str().c_str(),ios::out | ios::binary);
	fout.write(( char * ) _ARP2_3::arp2_3, Number_Arp2_3*sizeof( _ARP2_3::_arp2_3 ) );
	fout.close();
#endif

	
}
void restart_binary_input(const int myrank,const int numprocs,int restart_step){
	if(!myrank)fout_err<<"restart input "<<restart_step<<endl;
	_pointer_set pointset;
	stringstream ss;
	ss.str("");
	ifstream fin;
	
#ifdef SRD
	ss<<"./output"<<setw(4)<<setfill('0')<<RESTART_N<<"/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_mpc.txt";
//	ss<<"./output0/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_mpc.txt";
	fin.open(ss.str().c_str(),ios::in | ios::binary);
	fin.read(( char * ) mpc,Number_mpc*sizeof( _mpc ) );
	fin.close();
	fin.clear();
#endif
	
//cout<<myrank<<" 01"<<endl;
	int np;
	ss.str("");
	ss<<"./output"<<setw(4)<<setfill('0')<<RESTART_N<<"/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_ptcl.txt";
//	ss<<"./output0/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_ptcl.txt";
	fin.open(ss.str().c_str(),ios::in | ios::binary);
	fin.read(( char * )&np ,sizeof( int ) );
//cout<<myrank<<" 02 "<<np<<" "<<Number_Ptcl<<endl;
	if(np!=Number_Ptcl)exit(0);
//cout<<myrank<<" 03"<<endl;
	fin.read(( char * )&Active_Ptcl ,sizeof( int ) );
	fin.read(( char * )&pointset ,sizeof( _pointer_set ) );
	fin.read(( char * )ptcl.pt,Number_Ptcl*sizeof( _ptcl ) );
	fin.read(( char * )ptcl.con,Number_Ptcl*sizeof( _config ) );
	fin.read(( char * )plus_end  ,sizeof( _Actin::_Plus_End ) );
	fin.read(( char * )minus_end ,sizeof( _Actin::_Minus_End ) );
	fin.close();
	fin.clear();
//cout<<myrank<<" 04"<<endl;

#ifdef WALL_BOUNDARY
	ss.str("");
	ss<<"./output"<<setw(4)<<setfill('0')<<RESTART_N<<"/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_bind_wall.txt";
//	ss<<"./output0/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_bind_wall.txt";
	fin.open(ss.str().c_str(),ios::in | ios::binary);
	fin.read(( char * ) WALL::bind_wall, WALL::bind_size*sizeof( _bind_wall ) );
	fin.read(( char * ) WALL::WALL_Z ,2*sizeof( _WALL_Z ) );
	fin.close();
	fin.clear();
	for(_bind_wall* p=WALL::bind_wall;p!=WALL::bind_wall+WALL::bind_size;++p){
		if(p->index)p->index =&ptcl[0] + (p->index  - pointset.p);
	}
#endif //#ifdef WALL_BOUNDARY
#ifdef ARP2_3
	ss.str("");
	ss<<"./output"<<setw(4)<<setfill('0')<<RESTART_N<<"/restart/"<<setw(10)<<setfill('0')<<restart_step<<"_arp2_3.txt";
	fin.open(ss.str().c_str(),ios::out | ios::binary);
	fin.read(( char * ) _ARP2_3::arp2_3, Number_Arp2_3*sizeof( _ARP2_3::_arp2_3 ) );
	fin.close();
#endif
	{
//cout<<myrank<<" "<<1<<endl;
		for(_list* p=list;p!=list+list_size;++p){
			p->init_all();
		}
//cout<<myrank<<" "<<2<<endl;
		for(int i(0);i<Number_Ptcl;i++){
			ptcl.pt[i].index=i;
			ptcl.pt[i].cnf=&ptcl.con[i];
			ptcl.con[i].ptcl=&ptcl.pt[i];
		}
//cout<<myrank<<" "<<3<<endl;
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl];++p){//並列化組んでると違うことがありえる
			p->loc_pre=p->loc;
			p->quat_pre=p->quat;
		}
//cout<<myrank<<" "<<4<<endl;
//		for(_ptcl* p=ptcl;p!=ptcl+Number_Ptcl;++p){
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Active_Ptcl];++p){
			_list* li_sub=p->list_calc();
			p->add(li_sub);
			li_sub->add(p);
		}
//cout<<myrank<<" "<<5<<endl;

#ifdef ARP2_3
//		for(_ARP2_3::_arp2_3* p=_ARP2_3::arp2_3;p!=_ARP2_3::arp2_3+Number_Arp2_3;++p){
//			if(p->edge_ptcl)p->edge_ptcl=ptcl + (p->edge_ptcl - pointset.p);
//			if(p->end_ptcl) p->end_ptcl =ptcl + (p->end_ptcl  - pointset.p);
//		}
//cout<<myrank<<" "<<6<<endl;
		for(_ARP2_3::_arp2_3* p=_ARP2_3::arp2_3;p!=_ARP2_3::arp2_3+Number_Arp2_3;++p){
			_list* li_sub=p->list_calc();
			p->add(li_sub);
			li_sub->add(p);
		}
#endif
//cout<<myrank<<" "<<7<<endl;
		for(_ptcl* p=&ptcl[0];p!=&ptcl[Number_Ptcl];++p){
			if(p->cnf->plus) p->cnf->plus =&ptcl[0] + (p->cnf->plus  - pointset.p);
			if(p->cnf->minus)p->cnf->minus=&ptcl[0] + (p->cnf->minus - pointset.p);
			if(p->cnf->edge_arp2_3)p->cnf->edge_arp2_3=_ARP2_3::arp2_3 + (p->cnf->edge_arp2_3 - pointset.ar);
			if(p->cnf->end_arp2_3) p->cnf->end_arp2_3 =_ARP2_3::arp2_3 + (p->cnf->end_arp2_3  - pointset.ar);
		}
//		for(int i=0;plus_end->end[i];++i){
//			plus_end->end[i]=ptcl + (plus_end->end[i]  - pointset.p);
//		}
//		for(int i=0;minus_end->end[i];++i){
//			minus_end->end[i]=ptcl + (minus_end->end[i]  - pointset.p);
//		}
	}
//cout<<myrank<<" "<<8<<endl;
	
	step=restart_step;
	t += restart_step * PTCL_dt * INTVL_COL * INTVL_step;
	if(!myrank)fout_err<<"restart done"<<endl;
	return;
}

#ifdef BROWNIAN
void debug(){
	stringstream ss;
	ofstream fout;
	ss.str("");
	ss<<OUTPUT_N<<"ptcl.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"step=\t"<<step<<"\tt=\t"<<t<<endl;
	fout<<"\tlocation\t\t\tforce\t\t\tfila_number\tnumber"<<endl;
	for( int i = 0; i < Number_Ptcl; i++){
		_ptcl* p = &ptcl[i];
		fout<<i<<'\t'<<p->loc<<'\t'<<p->force<<'\t'
			<<p->cnf->feature<<'\t'
//			<<p->num
			<<endl;
	}
	fout.close();
	fout.clear();
#ifdef SRD
	ss.str("");
	ss<<OUTPUT_N<<"mpc.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"\tlocation\t\t\tvelocity\t\t\tforce\t\t\tmass"<<endl;
	for( int i = 0; i < Number_mpc; i++){
		_mpc* p = &mpc[i];
		fout<<i<<'\t'<<p->loc<<'\t'<<p->vel<<'\t'
			<<p->mass<<'\t'<<endl;
	}
	fout.close();
	fout.clear();
#endif
	ss.str("");
	ss<<OUTPUT_N<<"list.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"total\tptcl_index"<<endl;
	for( int i = 0; i < list_size ; i++){
		_list* p = &list[i];
		fout<<p->num;
		for(int j=0 ; j<p->num ; j++){
			fout<<'\t'<<setw(8)<<p->ptcl_index[j] - &ptcl[0];
		}
		fout<<endl;
	}
	fout.close();
	fout.clear();
	
// 	ss.str("");
// 	ss<<OUTPUT_N<<"bond.xls";
// 	fout.open(ss.str().c_str(),ios::out );
// 	fout<<"index\tplus\tminus\tb_plus\tb_minus"<<endl;
// 	for(_Actin::_Bond* b=bond;b->plus;++b){
// 		fout<<b-bond;
// 		fout<<'\t'<<setw(8)<<b->plus-ptcl<<'\t'<<b->minus-ptcl<<'\t'<<b->b_plus-bond<<'\t'<<b->b_minus-bond;
// 		fout<<endl;
// 	}
// 	fout.close();
// 	fout.clear();
}
#else //BROWNIAN
/*void debug(){
	stringstream ss;
	ofstream fout;
	ss.str("");
	ss<<OUTPUT_N<<"ptcl.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"\tlocation\t\t\tvelocity\t\t\tforce\t\t\tmass\tfila_number\tnumber"<<endl;
	for( int i = 0; i < Number_Ptcl; i++){
		_ptcl* p = &ptcl[i];
		fout<<i<<'\t'<<p->loc.x<<'\t'<<p->loc.y<<'\t'<<p->loc.z<<'\t'
			<<p->vel.x<<'\t'<<p->vel.y<<'\t'<<p->vel.z<<'\t'
			<<p->frc_new.x<<'\t'<<p->frc_new.y<<'\t'<<p->frc_new.z<<'\t'
			<<p->mass<<'\t'<<p->feature<<'\t'<<p->num<<endl;
	}
	fout.close();
	fout.clear();
#ifdef SRD
	ss.str("");
	ss<<OUTPUT_N<<"mpc.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"\tlocation\t\t\tvelocity\t\t\tforce\t\t\tmass"<<endl;
	for( int i = 0; i < Number_mpc; i++){
		_mpc* p = &mpc[i];
		fout<<i<<'\t'<<p->loc.x<<'\t'<<p->loc.y<<'\t'<<p->loc.z<<'\t'
			<<p->vel.x<<'\t'<<p->vel.y<<'\t'<<p->vel.z<<'\t'
			<<p->mass<<'\t'<<endl;
	}
	fout.close();
	fout.clear();
#endif
	ss.str("");
	ss<<OUTPUT_N<<"list.xls";
	fout.open(ss.str().c_str(),ios::out );
	fout<<"total\tptcl_index"<<endl;
	for( int i = 0; i < list_size ; i++){
		_list* p = &list[i];
		fout<<p->num;
		for(int j=0 ; j<p->num ; j++){
			fout<<'\t'<<setw(8)<<p->ptcl_index[j] - ptcl;
		}
		fout<<endl;
	}
	fout.close();
	fout.clear();
}*/

#endif //BROWNIAN
void output_cc(const int myrank, const int numprocs){
	static int counter_pol_sub(counter_pol);
	static int counter_depol_sub(counter_depol);
	static int counter_severing_sub(counter_severing);
	{
		static bool flag(true);
		if(flag){
			std::stringstream ss;
			ss<<OUTPUT_N<<"conf_chang.dat";
			std::ofstream fout;
if(RESTART){			fout.open(ss.str().c_str(),std::ios::out | std::ios::app);
}else{			fout.open(ss.str().c_str(),std::ios::out);
}			fout<<"step\tpolymerizasion\tdepolymerizasion\tsevering"<<endl;
			fout<<step<<'\t'<<counter_pol<<'\t'<<counter_depol<<'\t'<<counter_severing<<endl;
			fout.close();
			flag=false;
		}
	}
	
	if((counter_pol!=counter_pol_sub)||(counter_depol!=counter_depol_sub)||(counter_severing!=counter_severing_sub)){
		std::stringstream ss;
		ss<<OUTPUT_N<<"conf_chang.dat";
LOOP:
		ofstream fout(ss.str().c_str(),std::ios::out | std::ios::app);
			if(!fout.is_open()){
				static int count_err(0);
				if(count_err<10)count_err++;else{fout_err<<"PTCL_output"<<endl;	exit(0);}
				ss<<count_err;
				goto LOOP;
			}

		fout<<step<<'\t'<<counter_pol<<'\t'<<counter_depol<<'\t'<<counter_severing<<endl;
		
		fout.close();
		counter_pol_sub=counter_pol;
		counter_depol_sub=counter_depol;
		counter_severing_sub=counter_severing;
	}
	
	return;
}
void output_polymerization_info(const int myrank, const int numprocs){
	static string ss(OUTPUT_N()+"polymerization_info.dat");
	static bool flag(true);
	if(flag){
//		std::ofstream fout(ss.c_str(),std::ios::out);
		ofstream fout;//(ss.c_str(),ios::out | ios::app);
if(RESTART)fout.open(ss.c_str(),ios::out | ios::app);
else fout.open(ss.c_str(),ios::out);
		if(!fout.is_open()){
			static int count_err(0);
			if(count_err<10)count_err++;else{fout_err<<"output_polymerization_info_err1"<<endl;exit(0);}
			rename_restart(ss);//頭に"_"を挿入します
			output_polymerization_info(myrank,numprocs);
			return;
		}
		fout<<PolymerizationInfo.header()<<endl;
		fout.close();
		flag=false;
	}
	std::ofstream fout(ss.c_str(),std::ios::out |std::ios::app);
	if(!fout.is_open()){
		static int count_err(0);flag=true;
		if(count_err<10)count_err++;else{fout_err<<"output_polymerization_info_err"<<endl;exit(0);}
		rename_restart(ss);//頭に"_"を挿入します
		output_polymerization_info(myrank,numprocs);//再帰で呼べばいいんじゃね
		return;
	}

	
	PolymerizationInfo.sort();
	for(int i(0);i<PolymerizationInfo.size();++i){
		fout<<PolymerizationInfo.buf[i]<<endl;;
	}
	PolymerizationInfo.clear();
	return;
}


using _ARP2_3::_arp2_3;
using _ARP2_3::arp2_3;
using _ARP2_3::EPSILON_LJ_ACTIN_ARP;
using _ARP2_3::SIGMA_LJ_ACTIN_ARP;
using _ARP2_3::CutOffUperActinArp2_3;
using _ARP2_3::CutOffUperActinArp2_3_2;


const double r_wall_ptcl = r_cut;//pow(2.0,1.0/6.0)*sigma_LJ;//2.50 * sigma_LJ;//sigma_LJ*pow(2.0,1.0/6.0)*2.0;
const double r_wall_ptcl2 = r_wall_ptcl*r_wall_ptcl;
const double r_wall_arp = _ARP2_3::Cut_Off_Actin_Arp2_3;//pow(2.0,1.0/6.0)*SIGMA_LJ_ACTIN_ARP;
const double r_wall_arp2 = r_wall_arp*r_wall_arp;


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


void output_frc_wall(const int myrank, const int numprocs){
#ifdef WALL_BOUNDARY
	static string fna_frc;
	static string fna_arg;
	static string fna_arg_ave;
	{
		static bool flag(true);
		if(flag){
			std::stringstream ss;
			ss.str("");
			{
				ss<<OUTPUT_N<<"wall/";
				MyMkdir(ss.str());
			}
			ss.str("");
			ss<<OUTPUT_N<<"wall/frc.dat";
			std::ofstream fout(ss.str().c_str(),std::ios::out);
			fna_frc=ss.str();
			fout.close();
			
			ss.str("");
			ss<<OUTPUT_N<<"wall/arg.dat";
			fout.open(ss.str().c_str(),std::ios::out);
			fna_arg=ss.str();
			fout.close();
			flag=false;
			
			ss.str("");
			ss<<OUTPUT_N<<"wall/arg_ave.dat";
			fout.open(ss.str().c_str(),std::ios::out);
			fna_arg_ave=ss.str();
			fout.close();
			flag=false;
		}
	}
	ofstream fout(fna_frc.c_str(),std::ios::out | std::ios::app);
			if(!fout.is_open()){
				static int count_err(0);
				if(count_err<10)count_err++;else{fout_err<<"output_frc_wall"<<endl;exit(0);}
				rename_restart(fna_frc);//頭に"_"を挿入します
				output_frc_wall(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}

	fout<<step;
//	for(_ptcl** end=&plus_end->end[0];*end;++end){
	for(int i_end(0);i_end<plus_end->num;++i_end){
		_ptcl* end = plus_end->end(i_end);
		_ptcl* p = end;
		if( p->loc.z<r_wall_ptcl+(offset_impound*SYS_Z/(double)list_Z)){
			_vec<double> r = _vec<double>(0.0,0.0,p->loc.z-(offset_impound*SYS_Z/(double)list_Z));
			if(r.z>0){
				double potential;
				_vec<double> force_sub = LJ_ptcl(&r,&potential);
				
				fout<<'\t'<<force_sub.z;
			}
		}
		if( p->loc.z>WALL::WALL_Z[1].loc.z - r_wall_ptcl){
			_vec<double> r = _vec<double>(0.0,0.0,p->loc.z - WALL::WALL_Z[1].loc.z);
			if(r.z<0){
				double potential;
				_vec<double> force_sub = LJ_ptcl(&r,&potential);
				
				fout<<'\t'<<force_sub.z;
			}
		}
	}
	fout<<endl;
	fout.close();
	int count(0);
	double sum(0.0);
	fout.open(fna_arg.c_str(),std::ios::out | std::ios::app);
			if(!fout.is_open()){
				static int count_err(0);
				if(count_err<10)count_err++;else{fout_err<<"output_frc_wall"<<endl;exit(0);}
				rename_restart(fna_arg);//頭に"_"を挿入します
				output_frc_wall(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
	fout<<step;
//	for(_ptcl** end=&plus_end->end[0];*end;++end){
	for(int i_end(0);i_end<plus_end->num;++i_end){
		_ptcl* end = plus_end->end(i_end);
		_vec<double> vec;
		if((end)->cnf->minus)vec=(end)->cnf->minus->loc - (end)->loc;
		else if((end)->cnf->end_arp2_3)vec=(end)->cnf->end_arp2_3->loc - (end)->loc;
		else continue;
		distance_boundary(vec);
		double sub(vec.z/vec.norm());
		if(sub>=1.0)sub=1.0;
		if(sub<=-1.0)sub=-1.0;
//		if(vec.x<0){vec = -1.0*vec;}
		double the = acos(sub);
//		double phi = asin(vec.y/(sqrt(vec.x*vec.x + vec.y*vec.y)));
		fout<<'\t'<<the;
		count++;
		sum+=the;
	}
	fout<<endl;
	fout.close();
	
	fout.open(fna_arg_ave.c_str(),std::ios::out | std::ios::app);
			if(!fout.is_open()){
				static int count_err(0);
				if(count_err<10)count_err++;else{fout_err<<"output_frc_wall"<<endl;exit(0);}
				rename_restart(fna_arg_ave);//頭に"_"を挿入します
				output_frc_wall(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}
	fout<<step<<'\t'<<count<<'\t'<<sum/(double)count<<std::endl;
	fout.close();
	
	return;
#endif// WALL_BOUNDARY
}

void ptcl_wall_interaction(const int myrank, const int numprocs){
#ifdef WALL_BOUNDARY
	{
		static bool flag(true);
		if(flag){
			std::stringstream ss;
			ss.str("");
			{
				ss<<OUTPUT_N<<"ptcl/";
				MyMkdir(ss.str());
			}
			flag=false;
		}
	}
	
	stringstream ss;
	ss.str("");
	ss<<OUTPUT_N<<"ptcl/output"<<setw(8)<<setfill('0')<<step<<".dat";
	ofstream fout(ss.str().c_str(),ios::out );
	
	fout<<step<<'\t'<<Number_Ptcl<<'\t'<<Active_Ptcl<<'\t'<<WALL::WALL_Z[1].loc.z<<'\t'<<r_cut<<'\t'<<Wall_RandamF<<std::endl;
	
	for(int i(0); i<Active_Ptcl; i++){
		_ptcl& p = ptcl[i];
		fout<<i<<'\t'<<p.loc.x<<'\t'<<p.loc.y<<'\t'<<p.loc.z<<'\t'<<p.cnf->feature<<'\t';
		if(p.cnf->minus)fout<<p.cnf->minus->index<<"\t";
		else if(p.cnf->end_arp2_3)fout<<-2<<"\t";
		else fout<<-1<<"\t";
		if(p.cnf->plus)fout<<p.cnf->plus->index<<"\t";
		else fout<<-1<<"\t";
		fout<<filament_force[i];
		fout<<endl;
	}
	
	fout.close();
#endif// WALL_BOUNDARY
}
void ptcl_number_output(const int myrank, const int numprocs){
	if(myrank)return;
	static bool flag(true);
	static string ss(OUTPUT_N()+"ptcl_number.dat");
	{
		if(flag){
//			std::ofstream fout(ss.c_str());
			ofstream fout;//(ss.c_str(),ios::out | ios::app);
if(RESTART)fout.open(ss.c_str(),ios::out | ios::app);
else fout.open(ss.c_str(),ios::out);
				if(!fout.is_open()){
					static int count_err(0);flag=true;
					if(count_err<10)count_err++;else{fout_err<<"ptcl_number_output_err"<<endl;exit(0);}
					rename_restart(ss);//頭に"_"を挿入します
					ptcl_number_output(myrank,numprocs);//再帰で呼べばいいんじゃね
					return;
				}
			fout<<"load_is_"<<WALL::jijuu<<" NumberPtcl"<<Number_Ptcl<<endl;
			fout<<"step\tActivePtclNumber\n";
			fout.close();
			flag=false;
		}
	}
	static int tmp_s[100];
	static int tmp_AP[100];
	static int count_ap(0);
	{
		tmp_s[count_ap]=step;
		tmp_AP[count_ap]=Active_Ptcl;
		count_ap++;
	}
	if (count_ap>=100 || step>STEP_END-10){
		ofstream fout(ss.c_str(),ios::out |std::ios::app);
			if(!fout.is_open()){
				static int count_err(0);flag=true;
				if(count_err<10)count_err++;else{fout_err<<"ptcl_number_output_err"<<endl;exit(0);}
				rename_restart(ss);//頭に"_"を挿入します
				ptcl_number_output(myrank,numprocs);//再帰で呼べばいいんじゃね
				return;
			}

		for(int i(0);i<count_ap && i<100;++i){
			fout<<tmp_s[i]<<'\t'<<tmp_AP[i]<<std::endl;
		}
		fout.close();
		count_ap=0;
	}
}
