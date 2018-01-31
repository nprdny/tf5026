#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "tf5026.h"

#define Nmx 20000

/* Place derivs_int() here ... */
void derivs_int(double x,double y[],double dydx[],int k, double f[], double r[],double dr[], int N)
{
	dydx[1]=f[k]*dr[k];
}

int main(){
	/* vektor-vektor grid */
	double *r, *dr, *f;
	
	/* variabel-variabel kerja */
	int i, N;
	double u;
	double *y, *dydx; //NR vector for dif. eq's
	double x, h;
	int k, dk;
	
	/* define pi */
	double pi;
	pi = 4.*atan(1);
	
	/* alokasi vektor hingga maksimum */
	r = dvector(0, Nmx);
	dr = dvector(0, Nmx);
	f = dvector(0, Nmx);
	y = dvector(1, 1);
	dydx = dvector(1, 1);
	
	/* loop over different integration mesh sizes */
	for(N=2; N<=Nmx; N*=10){
		/* konstruksi perubahan variabel dan fungsi integran untuk menghitung potensial akibat interaksi nukleus-elektron */
		for(i=0; i<=N; i++){
			u = ((double)i)/((double)N);
			r[i] = 1./(1.-u)-1;
			dr[i] = 1./(1.-u)/(1.-u);
			f[i] = (-1/r[i])*(exp(-2*r[i])/pi)*4.*pi*r[i]*r[i];
		}
		/* syarat-syarat khusus untuk menangani integrasi di titik-titik ujung */
		dr[0] = f[0] = 0;
		dr[N] = f[N] = 0;
		
		/* solusi runge-kutta menggunakan rk4tf5026()*/
		/* ukuran interval (h,dk)=2 krn algoritma RK memberi hasil hanya di titik-titik ujung */
		h = 2./N;
		dk = 2;
		y[1] = 0.;
		for(k=0; k<=N-2; k+=dk){
			x = k*h; //kalo remi k/2*h (?????)
			derivs_int(x,y,dydx,k,f,r,dr,N);
			rk4tf5026(y,dydx,1,x,h,y,derivs_int,k,dk,f,r,dr,N);
		}
		printf("N,T,S,Rk: %i %20.16f %20.16f %20.16f\n", N, trapint(f,r,dr,N), simpint(f,r,dr,N), y[1]);
		//printf("N,T,S,Rk: %6 %20.16f %20.16f\n", N, trapint(f,r,dr,N), simpint(f,r,dr,N));
	}
	/* deallocate NR vectors: should always clean up space at the end */
	free_dvector(r,0,Nmx);
	free_dvector(dr,0,Nmx);
	free_dvector(f,0,Nmx);
	free_dvector(y,1,1);
	free_dvector(dydx,1,1);
}
