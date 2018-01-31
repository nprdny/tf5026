#ifndef _tf5026_H_
#define _tf5026_H_

double trapint(double f[], double r[], double dr[], int N);

double simpint(double f[], double r[], double dr[], int N);

void rk4tf5026(double y[], double dydx[], int n, double x, double h, double yout[], void(*derivs)(double x, double y[], double dydx[], int k, double f[], double r[], double dr[], int N), int k, int dk, double f[], double r[], double dr[], int N);
			   
#endif