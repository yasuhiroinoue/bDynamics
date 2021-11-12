/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
//mylib.cpp

#include <math.h>

#include "mymath.h"


/* Period parameters */  

#define N_MT 624

#define M_MT 397

#define MATRIX_A 0x9908b0dfUL   /* constant vector a */

#define UMASK 0x80000000UL /* most significant w-r bits */

#define LMASK 0x7fffffffUL /* least significant r bits */

#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )

#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))



static unsigned long state[N_MT]; /* the array for the state vector  */

static int left = 1;

static int initf = 0;

static unsigned long *next;


//ÇQéüå≥çsóÒÇ…ä÷Ç∑ÇÈíËã`

//MATRIX matrix_new(int col, int row){

//	MATRIX a = new double *[col];

//	for( int i = 0; i < row; i++)

//			a[i] = new double[row];

//

//	return a;

//}

//

//void matrix_delete(MATRIX a){

//	MATRIX b = a;

//	

//	while( *b != 0 ){

//		delete [] *b++;

//	}

//	delete [] a;

//}



//êîäwä÷êî

/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */

double gauss(double t, double u){
	double ran, gatai;
	
	do ran = RAND; while(ran == 0.0);
	gatai = (double)(sqrt(-2.0 * t * log(ran)) * sin(PI2 * RAND)) + u;
	return(gatai);
}

double gamma_noflow(double t){
	double ran, gamma;
	do ran = RAND; while(ran == 0.0);
	gamma = sqrt(-2.0 * t * log(ran));
	return(gamma);
}
/* initializes state[N] with a seed */

void init_genrand(unsigned long s){
	int j;
	state[0]= s & 0xffffffffUL;
	for (j=1; j<N_MT; j++) {
		state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */

        /* In the previous versions, MSBs of the seed affect   */

        /* only MSBs of the array state[].                        */

        /* 2002/01/09 modified by Makoto Matsumoto             */

		state[j] &= 0xffffffffUL;  /* for >32 bit machines */
	}
	left = 1; initf = 1;
}

static void next_state(void){
	unsigned long *p=state;
	int j;

    /* if init_genrand() has not been called, */

    /* a default initial seed is used         */

	if (initf==0) init_genrand(5489UL);
	left = N_MT;
	next = state;
	for (j=N_MT-M_MT+1; --j; p++) *p = p[M_MT] ^ TWIST(p[0], p[1]);
	for (j=M_MT; --j; p++) *p = p[M_MT-N_MT] ^ TWIST(p[0], p[1]);
	*p = p[M_MT-N_MT] ^ TWIST(p[0], state[0]);
}

/* generates a random number on (0,1)-real-interval */

double genrand_real3(void){
	unsigned long y;
	if (--left == 0) next_state();
	y = *next++;
    /* Tempering */

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);
	return ((double)y + 0.5) * (1.0/4294967296.0); 
    /* divided by 2^32 */

}
