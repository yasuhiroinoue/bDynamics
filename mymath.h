/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
//mylib.h



#ifndef MYLIB_H 

#define MYLIB_H 



//２次元行列に関する定義

//typedef double** MATRIX;

//MATRIX matrix_new(int, int);

//void matrix_delete(MATRIX);



//数学関数

static double HPI = 1.57079632679490;
static double PI =  3.14159265358979;
static double PI2 = 6.28318530717959;

#define TANE 10



#define SQR(x) (x * x)

#define ABS(x) sqrt(SQR(x))


double genrand_real3(void);

void init_genrand(unsigned long s);

double gauss( double t, double u);

double gamma_noflow( double t);

#define RAND (double)genrand_real3()





#endif

