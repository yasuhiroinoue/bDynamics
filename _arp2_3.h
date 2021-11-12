/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef ___ARP2_3_H
#define ___ARP2_3_H

#include "_system.h"
#include "_class.h"

#define Number_Arp2_3 0//64//128//0//128//0//128//0//128//0//128//128


namespace _ARP2_3{
//---------------------------------ポテンシャル係数-------------------------------//

const double SIGMA_LJ_ARP_ARP		=(12.8/5.4);
const double SIGMA_LJ_ACTIN_ARP		=((SIGMA_LJ_ARP_ARP+sigma_LJ)*0.5);
const double EPSILON_LJ_ARP_ARP		=0.268292683;			//0.000268293
const double EPSILON_LJ_ACTIN_ARP = std::sqrt(EPSILON_LJ_ARP_ARP * epsilon_LJ);

const double K_TENSILE_ACTIN_ARP	=1000.000;
		//バネ定数
const double EQ_LEN_ACTIN_ARP		=(SIGMA_LJ_ARP_ARP*0.5);

const double K_ANGLE_ACTIN_ARP		=k_angl;		//10.0//20.0//1.000
const double EQ_ANGLE_ACTIN_ARP		=(PI*70.0/180.0);

const double K_Torsion_ARP			=0.1*k_angl;	//10.0
//0.50000
//1.0000

const double Cut_Off_Arp2_3_Arp2_3	=SIGMA_LJ_ARP_ARP*pow(2.0,1.0/6.0);//2.5*SIGMA_LJ_ARP_ARP;
const double Cut_Off_Actin_Arp2_3	=SIGMA_LJ_ACTIN_ARP*pow(2.0,1.0/6.0);//2.5*SIGMA_LJ_ACTIN_ARP;
const double Cut_Off_Arp2_3_Arp2_3_2=Cut_Off_Arp2_3_Arp2_3*Cut_Off_Arp2_3_Arp2_3;
const double Cut_Off_Actin_Arp2_3_2	=Cut_Off_Actin_Arp2_3 *Cut_Off_Actin_Arp2_3;


const double CutOffUperArp2_3Arp2_3	=cutoffuper_sub*SIGMA_LJ_ARP_ARP;
const double CutOffUperArp2_3Arp2_3_2=CutOffUperArp2_3Arp2_3*CutOffUperArp2_3Arp2_3;
const double CutOffUperActinArp2_3	=cutoffuper_sub*SIGMA_LJ_ACTIN_ARP;
const double CutOffUperActinArp2_3_2=CutOffUperActinArp2_3*CutOffUperActinArp2_3;


const double pro_bind=0.0025;					//0.0000001;//0.000005;//0.0000001;//0.001;

//--------------------------------------------------------------------------------//

const double  Arp2_3_Friction		=(Ptcl_Friction*3.0);
//#define Arp2_3_Ror_Friction Arp2_3_Friction*4.0/3.0*(SIGMA_LJ_ARP_ARP*0.5)*(SIGMA_LJ_ARP_ARP*0.5)

	extern _arp2_3* arp2_3;
// 	extern _Bind* bind_arp2_3;
	extern double ene_pot;
	namespace _ARP_ARP{
		void Frc(const int);
	}
	namespace _ACTIN_ARP{
		void Frc(const int);
	}
	void Bind(const int, const int);
	void arp_init(void);
	void arp_posion_init(void);
	void arp_finalize();
	void Arp_Location_FirstStep(const int myrank,const int numprocs);
	void Arp_Location_SecondStep(const int myrank,const int numprocs);
}






#endif //___ARP2_3_H
