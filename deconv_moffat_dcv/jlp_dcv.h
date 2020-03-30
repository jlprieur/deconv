/**************************************************************************
* Definitions needed by the deconvolution programs
*
* JLP
* Version 21/12/2006
**************************************************************************/
#include <stdio.h>
#include <stdlib.h> /* RAND_MAX is defined here! */
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <jlp_ftoc.h>
#include <jlp_fftw.h>

/******************* Prototypes: ****************************/

/* in "jlp_minimize.c": */ 
int jlp_minimize(double *x, int n, int nmax, double ftol, int *iter,
                 int iter_max, double *fret, double (*func)(double *),
                 void (*dfunc)(double *, double *), int positive, int lbfgs);

/* in "dcv_1D_sub.c": */ 
int dcv_1D_gmark(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_tikhonov(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_ggauss(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_sqrtrf(double alpha_0, double ss, double ftol, int *iter, int positive);
int dcv_1D_mem(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_spdiv();
int dcv_1D_wiener(double alpha_0);
int dcv_1D_banana(double ftol);
int dcv_1D_init(int nn1_0, int nn_0, int hh_nx_0);
void dcv_1D_free();
int dcv_1D_plot1(double *yy1, int nn, char *title);
int dcv_1D_plot1_log(double *yy1, int nn1, char *title);
int dcv_1D_plot2(double *yy1, double *yy2, int nn, char *title);
int dcv_1D_display_input(int nn);
int dcv_1D_transfert(double *hh);

/* in "dcv_2D_sub.c": */ 
int dcv_2D_gmark(double alpha_0, int positive, double ftol, int *iter,
                 double *L_phi, double *L_y);
int dcv_2D_tikhonov(double alpha_0, int positive, double ftol, int *iter,
                    double *L_phi, double *L_y);
int dcv_2D_ggauss(double alpha_0, int positive, double ftol, int *iter,
                  double *L_phi, double *L_y);
int dcv_2D_sqrtrf(double alpha_0, double ss, int positive, double ftol,
                  int *iter, double *L_phi, double *L_y);
int dcv_2D_mem(double alpha_0, double ftol, int *iter,
               double *L_phi, double *L_y);
int dcv_2D_spdiv();
int dcv_2D_wiener(double alpha_0);
int dcv_2D_init(int nx_0, int ny_0);
void dcv_2D_free();
int dcv_2D_noisy_signal(double snr);
int dcv_2D_transfert(double *hh);

/* splot/lib/jlp_display1.cpp */
int jlp_display1(double *xx, double *yy, int istart, int iend,
                 char *xlabel, char *ylabel, char *title, char *plotdev);
int jlp_display2(double *xx1, double *yy1, int istart1, int iend1,
                 double *xx2, double *yy2, int istart2, int iend2,
                 char *xlabel, char *ylabel, char *title, char *plotdev,
                 char *nchar1, char *nchar2);                                                

#define MAXI(x,y) ((x) > (y) ? (x) : (y))
/* Warning: error with "abs", so I define my own macro: */
#define ABS(x) ((x) > 0. ? (x) : (-(x)))
#define SQUARE(x) ((x) * (x))
#define MY_LOG10(x) (((x) > 0) ? log10(x) : 0.)
