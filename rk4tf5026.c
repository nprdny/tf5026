#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void rk4tf5026(double y[], double dydx[], int n, double x, double h, double yout[], void(*derivs)(double x, double y[], double dydx[], int k, double f[],double r[], double dr[], int N), int k, int dk, double f[], double r[], double dr[], int N)

/*Given values for the variables y[1..n] and their derivatives dydx[1..n] known at x, use the fourth-order Runge-Kutta method to advance the solution over an interval h and return the incremented variables as yout[1..n], which need not be a distinct array from y. The user supplies the routine derivs(x,y,dydx), which returns derivatives dydx at x.*/

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
	
	//derivs_int(double x, double y[], double dydx[], int k, double *f[], double r[], double dr[], int N)

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



