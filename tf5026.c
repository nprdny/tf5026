#include <stdio.h>
#include <math.h>
#include "nrutil.h"

#define Nmx 20000

/* Fungsi derivs_int */
void derivs_int(double x,double y[],double dydx[],int k, double f[], double r[],double dr[], int N)
{
	dydx[1]=f[k]*dr[k];
}

/* Metode Runge Kutta */
void rk4tf5026(double y[], double dydx[], int n, double x, double h, double yout[], void(*derivs)(double x, double y[], double dydx[], int k, double f[],double r[], double dr[], int N), int k, int dk, double f[], double r[], double dr[], int N)

{
int i;
double xh,hh,h6,*dym,*dyt,*yt;

dym = dvector(1,n); // at this point: k_3
dyt = dvector(1,n); // at this point: k_2
yt = dvector(1,n); // sepertinya nilai awalnya = 0?
hh = h*0.5; // h = step size
h6 = h/6.0; 
xh = x+hh;

for (i=1;i<=n;i++)
	yt[i]=y[i]+hh*dydx[i]; // Array berisi y untuk bikin k_2
	(*derivs)(xh,yt,dyt,k,f,r,dr,N); // dyt=k_2

for (i=1;i<=n;i++)
	yt[i]=y[i]+hh*dyt[i]; // Array berisi y untuk bikin k_3
	(*derivs)(xh,yt,dym,k,f,r,dr,N); // dym=k_3

for (i=1;i<=n;i++)
{
    yt[i]=y[i]+h*dym[i]; // Array berisi y untuk bikin k_4
	dym[i] += dyt[i]; // dym[i]=dym[i]+dyt[i]
}

(*derivs)(x+h,yt,dyt,k,f,r,dr,N); //Fourth step.

for (i=1;i<=n;i++) // Accumulate increments with proper weights
yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

free_dvector(yt,1,n);
free_dvector(dyt,1,n);
free_dvector(dym,1,n);
}

/* Metode integrasi trapesium */
double trapint(double f[], double r[], double dr[], int N){
    double trapz=0;
    int i;
    for(i=0; i<N; i++){
        trapz = trapz+((r[i+1]-r[i])*((f[i]+f[i+1])/2))*dr[i];
    }
    return trapz;
}

/* Metode integrasi simpson */
double simpint(double f[], double r[], double dr[], int N){
    double sum2=0, sum4=0, simp;
    int i;
    for(i=1; i<N; i++){
        if(i%2==0){
            sum2=sum2+2*f[i];
        }
        else{
            sum4=sum4+4*f[i];
        }
        simp = ((((r[N]-r[0])/(3*N))*(f[0]+f[N]))+(((r[i+1]-r[i])/3)*(sum2+sum4)))*dr[i];
    }
    return simp;
}

/* Loop untuk trap & simp */
//int main(){
    /* Deklarasi untuk arrays, ingat bahwa ada N+1 titik*/
    /* Untuk ketelitian perhitungan, digunakan double precision */

//    double r[Nmx+1], dr[Nmx+1], f[Nmx+1];

    /* Deklarasi variabel kerja */
 /*    int i, N;
    double u; */

    /* Trik untuk mendapatkan seluruh digit pi*/
/*     double pi;
    pi = 4.*atan(1.); */

    /* Loop untuk ukuran N = 2, 20, 200, ..., Nmx */
//    for (N=2; N<=Nmx; N*=10){
            /* Untuk setiap ukuran, isi array*/
 //           for (i=0; i<=N; i++){
                /* Anggap u mempunyai selang yang sama antara 0 dan 1*/
 //               u = ((double)i)/((double)N);
                /*Definisi pengubahan variabel*/
/*                 r[i] = u;
                dr[i] = 1.; */
                /* Pernyataan fungsi*/
/*                 f[i] =  1/(1+r[i]*r[i])/(1+r[i]*r[i])/(0.25+0.125*pi); */

                //printf("Nilai N: %d \n", N);
                //printf("Nilai i: %d \n", i);
                //printf("Nilai r_%i: %f \n", i, r[i]);
                //printf("Fungsi f_%i: %f \n \n", i, f[i]);
           // }

            /* Kasus khusus di titik akhir */
            /* f[0] = f[0];
            r[0] = r[0];
            dr[0] = dr[0];

            f[N] = f[N];
            r[N] = r[N];
            dr[N] = dr[N]; */

            /* Memanggil routine integrasi dan hasil perhitungan */
           /*  printf("Trap, Simp: %d %20.14f %20.14f \n", N, trapint(f,r,dr,N), simpint(f,r,dr,N));
	}
} */


main(){
	/* Grid vectors */
	double *r,*dr,*f;
	
	/* Working variables */
	int i,N;
	double u;
	double *y,*dydx; /* NR vector for diff eqs */
	double x,h;
	int k,dk;
	
	/* Value of pi */
	double pi;
	pi=4.*atan(1.);
	
	/* Allocate NR vectors for maximum size used */
	r=dvector(0,Nmx);
	dr=dvector(0,Nmx);
	f=dvector(0,Nmx);
	y=dvector(1,1);
	dydx=dvector(1,1);

	/* Loop over different integration mesh sizes */
	for (N=2; N<=Nmx; N*=10) {
		
		/* Construct change of variables information and integrand function for computing V_el-nuc */
		for (i=0; i<=N; i++) {
			u=((double) i)/((double) N);
			r[i]=1./(1.-u)-1.;
			dr[i]=1./(1.-u)/(1.-u);
			f[i]=(-1/r[i])*(exp(-2*r[i])/pi)*4.*pi*r[i]*r[i];
		}
		
		/* Special conditions for handling end-points of integration */
		dr[0]=f[0]=0.;
		dr[N]=f[N]=0.;
			
		/* Runga-Kutta solution using rk4p480(). */
			
		/* Note step size (h,dk) is 2 because RK algorithm produces 7 results only at the end points. */
		h=2./N;
		dk=2;
		y[1]=0.;
			
		for (k=0; k<=N-2; k+=dk) {
			x=k/2*h;
			derivs_int(x,y,dydx,k,f,r,dr,N);
			rk4tf5026(y,dydx,1,x,h,y,derivs_int,k,dk,f,r,dr,N);
		}
		printf("T,S,RK: %6d %18.14f,%18.14f,%18.14f\n",N,trapint(f,r,dr,N),simpint(f,r,dr,N),y[1]);
	}
/* Deallocate NR vectors: Should always clean up space at the end. */
free_dvector(r,0,Nmx);
free_dvector(dr,0,Nmx);
free_dvector(f,0,Nmx);
free_dvector(y,1,2);
free_dvector(dydx,1,2);
}