#include <stdio.h>
#include <math.h>
#include "nrutil.h"

#define Nmx 20000

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
    double sumeven=0, sumodd=0, simp;
    int i;
    for(i=1; i<N; i++){
        if(i%2==0){
            sumeven=sumeven+2*f[i];
        }
        else{
            sumodd=sumodd+4*f[i];
        }
        //simp = ((((r[N]-r[0])/(3*N))*(f[0]+f[N]))+(((r[i+1]-r[i])/3)*(sumeven+sumodd)));
		simp = ((r[N]-r[0])/(3*N))*((f[0]+f[N])+(sumeven+sumodd));
    }
    return simp;
}

/* Metode Runge Kutta */
void rk4tf5026(double y[], double dydx[], int n, double x, double h, double yout[], void(*derivs_int)(double x, double y[], double dydx[], int k, double f[],double r[], double dr[], int N), int k, int dk, double f[], double r[], double dr[], int N){
	
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
	(*derivs_int)(xh,yt,dyt,k,f,r,dr,N); // dyt=k_2
	
	for (i=1;i<=n;i++)
		yt[i]=y[i]+hh*dyt[i]; // Array berisi y untuk bikin k_3
	(*derivs_int)(xh,yt,dym,k,f,r,dr,N); // dym=k_3
	
	for (i=1;i<=n;i++){
		yt[i]=y[i]+h*dym[i]; // Array berisi y untuk bikin k_4
		dym[i] += dyt[i]; // dym[i]=dym[i]+dyt[i]
		}
		
	(*derivs_int)(x+h,yt,dyt,k,f,r,dr,N); //Fourth step.
	
	for (i=1;i<=n;i++) // Accumulate increments with proper weights
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
		
	free_dvector(yt,1,n);
	free_dvector(dyt,1,n);
	free_dvector(dym,1,n);
	
	//return 0;
}
