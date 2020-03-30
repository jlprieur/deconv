/***********************************************************
* Test of minimisation routines from "Numerical Recipees"
* dcv_test0
*
* JLP
* Version 12/03/2001
************************************************************/
#include <jlp_dcv.h>

double func_test(double x[]);
void dfunc_test(double x[], double dx[]);

int main(int argc, char *argv[])
{
#define NN 2
int n=NN, iter; 
double p[NN+1], ftol, fret; 
double dx[NN+1];

printf("Test with Banana function \n");

/* Starting point: */
      p[1]= -1.9; 
      p[2]= 2.;
/* Tolerance */
      ftol=1.e-3;
dfunc_test(p,dx);
printf("value=%.2f gradient=%.2f %.2f\n",func_test(p),dx[1],dx[2]);
/* Calling conjudate gradient routine in multi-dimensions */
fret=-1;
 frprmn(p, n, ftol, &iter, &fret, func_test, dfunc_test);

/* Results: */
  printf(" Location of the minimum: %.3f %.3f \n", p[1],p[2]);
  printf(" Value of the minimum: %.5e\n", fret);
  printf(" Number of iterations: %d\n", iter);

return(0);
}
/*-------------------------------------------------------------------
* Function: 
* Banana function : start at X=[-1.9;2]. minimum at X=[1;1] : f(X)=0;
*-------------------------------------------------------------------*/
double func_test(double x[])
{
double ww;
  ww = (1. - x[1])*(1. - x[1]) + 100.*(x[2] - x[1])*(x[2] - x[1]);
  return(ww);
}
/*-------------------------------------------------------------------
* Gradient:
* BANANA FUNCTION
* z = 100*(x(2)-x(1))^2 + (1-x(1))^2
*  dx = [  -200*(x(2) - x(1))-2*(1-x(1)); 200*(x(2)-x(1)) ];
*-------------------------------------------------------------------*/
void dfunc_test(double x[], double dx[])
{
  dx[1] = -200. * (x[2] - x[1]) - 2. * (1 - x[1]); 
  dx[2] = 200. * (x[2] - x[1]);
return;
}
