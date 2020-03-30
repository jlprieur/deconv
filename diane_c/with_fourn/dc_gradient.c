/***************************************************************/
/*                                                             */
/* PROGRAM: diane_gradient                                     */
/*                                                             */
/* PURPOSE: Deconvolution with multiresolution constraint      */
/*          using a mean square minimization involving         */
/*          gradient method                                    */
/*                                                             */
/* INPUT:  argv[1] = generic_name                              */
/*         argv[2] = threshold for residual                    */
/*                                                             */
/* NB: *.PHIT  first approximation of the object       (input) */
/*     *.SNR   signal to noise ration in fourier domain (input)*/
/*     *.HR    synthetic aperture                      (input) */
/*     *.V     multiresolution support                (output) */
/*     *.G     interpolation function                 (output) */
/*     *.KR    stabilizing function                   (output) */
/*                                                             */
/* VERSION: December 1992                                      */
/*                                                             */
/* AUTHOR:  E. ANTERRIEU  modified by K. BOUYOUCEF             */ 
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

extern   void      fourn();
extern   float **  alloc_square_float();
extern   void      free_square_float();
extern   int       decomposition_mallat_ncs();
extern   int       recomposition_mallat_ncs();

/***************************************************************/
main(argc,argv)
int argc;
char *argv[];
/***************************************************************/
{

/*      DECLARATIONS
*/
auto     int      dim, size[4];
auto     char     file[100], type[16], mode[8], nature[16], comments[1024];
auto     float    **phit, **G, **V, **phin;
auto     float    **phin1;
auto     int      nbr_lin, nbr_col, err, d_order;
auto     float    eps;
register int      i, j, k;
auto     int      isign, nn[2], ndim;
auto     float    *phit_FT;

auto     float    **rn, **zn;
auto     float    epsilon;
auto     int      iteration;
auto     float    norm_rn, norm_rn_sav, norm_phi_n, norm_dphi_n, phi1, phi2, ps_dn_zn;
auto     float    omega_n, gamma_n;
auto     float    *rn_FT, *zn_FT, *Gr;

/**********************************************************************************/


/*      VERIFICATION DE LA LIGNE DE COMMANDE
*/
if (argc != 3)
  {
  (void) printf("\nUnexpected number of arguments\n");
  (void) printf("\nUSAGE:\n");
  (void) printf("diane_gradient  generic_name  eps\n\n");
  (void) printf("\t eps = threshold for residual\n");
  (void) printf("\t \n\n");
  exit(-1);
  }

/*       LECTURE DES DIFFERENTS FICHIERS ET DES DATAs
*/
dim = 2;
(void) strcpy(type,"float");
(void) strcpy(nature,"real");

sprintf(file,"%s.PHIT",argv[1]);
phit = (float **) read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];

sprintf(file,"%s.G",argv[1]);
G = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s.G & %s.PHIT] have different sizes\n",argv[1],argv[1]);
  exit(-1);
  }

sprintf(file,"%s.V",argv[1]);
V = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s.V & %s.PHIT] have different sizes\n",argv[1],argv[1]);
  exit(-1);
  }

err = sscanf(comments,"J=%d",&d_order);
if ((err != 1) || (d_order < 0) || (d_order > 10))
   d_order = 0;

err = sscanf(argv[2],"%f",&eps);
if ((err != 1) || (eps < 0.000000001) || (eps > 100.0))
  {
  fprintf(stderr,"\nWARNING : incorrect threshold for residual [%.10f]\n",eps);
  fprintf(stderr,"\ndecomposition order will be imposed = 0\n");
  }

fprintf(stderr,"\n*****************************");
fprintf(stderr,"\nPROGRAM :      diane_gradient");
fprintf(stderr,"\n*****************************");
fprintf(stderr,"\ndecomposition order of the support constraint = %d",d_order);
fprintf(stderr,"\nthreshold for residual error = %f\n\n",eps);

/*       FOURIER TRANSFORM OF phit USING NUMERICAL RECIPES SUBROUTINE
*/
phit_FT = (float *) malloc (2*nbr_lin*nbr_col*sizeof(float));
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phit_FT[2*(j+i*nbr_col)] = phit[i][j];
  phit_FT[2*(j+i*nbr_col)+1] = 0.0;
  }
ndim = 2;
isign = -1;
nn[0] = nbr_lin;
nn[1] = nbr_col;
(void) fourn(phit_FT-1,nn-1,ndim,isign);

/*        TRANSFORM G IN VECTOR FORM
*/
Gr = (float *) malloc (nbr_lin*nbr_col*sizeof(float));
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  Gr[j+i*nbr_col] = G[i][j];
(void) free_square_float(G,nbr_lin);
(void) free_square_float(phit,nbr_lin);


/*      METHODE DU GRADIENT
*/

rn = (float **) alloc_square_float(nbr_lin,nbr_col);
zn = (float **) alloc_square_float(nbr_lin,nbr_col);
phin = (float **) alloc_square_float(nbr_lin,nbr_col);
phin1 = (float **) alloc_square_float(nbr_lin,nbr_col);

rn_FT = (float*) malloc (2*nbr_lin*nbr_col*sizeof(float));
zn_FT = (float*) malloc (2*nbr_lin*nbr_col*sizeof(float));

iteration = 0;
epsilon = 1.0;

while (epsilon > eps)
  {
  iteration++;
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
      rn[i][j] = phin[i][j];


/*      PERFORMS OPERATION  rn = g * (phin - rn)
        IN FOURIER DOMAIN
*/
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    {
    rn_FT[2*(j+i*nbr_col)] = rn[i][j];
    rn_FT[2*(j+i*nbr_col)+1] = 0.0;
    }
  ndim = 2;
  isign = -1;
  nn[0] = nbr_lin;
  nn[1] = nbr_col;
  (void) fourn(rn_FT-1,nn-1,ndim,isign);

  for (k=0; k<2*nbr_lin*nbr_col; k++)
      rn_FT[k] = phit_FT[k] - rn_FT[k];
  for (k=0; k<nbr_lin*nbr_col; k++)
      {
      rn_FT[2*k] *= Gr[k];
      rn_FT[2*k+1] *= Gr[k];
      }
  isign = 1;
  (void) fourn(rn_FT-1,nn-1,ndim,isign);

  for (k=0; k<nbr_lin*nbr_col; k++)
    {
    i = k/nbr_col;
    j = k-i*nbr_col;
    rn[i][j] = rn_FT[2*k]/(nbr_lin*nbr_col);
    }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
     phin1[i][j] = rn[i][j];

  /********  projection sur des espaces multiresolution  ***********/
  err = (int) decomposition_mallat_ncs(rn,nbr_lin,nbr_col,d_order);
  if (err != 0)
     {
     fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
     exit(-1);
     }
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
      rn[i][j] *= V[i][j];
  err = (int) recomposition_mallat_ncs(rn,nbr_lin,nbr_col,d_order);
  if (err != 0)
     {
     fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
     exit(-1);
     }

  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
      zn[i][j] = rn[i][j];


/*      PERFORMS OPERATION  zn = g * zn
        IN FOURIER DOMAIN
*/
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    {
    zn_FT[2*(j+i*nbr_col)] = zn[i][j];
    zn_FT[2*(j+i*nbr_col)+1] = 0.0;
    }
  ndim = 2;
  isign = -1;
  nn[0] = nbr_lin;
  nn[1] = nbr_col;
  (void) fourn(zn_FT-1,nn-1,ndim,isign);

  for (k=0; k<nbr_lin*nbr_col; k++)
      {
      zn_FT[2*k] *= Gr[k];
      zn_FT[2*k+1] *= Gr[k];
      }
  isign = 1;
  (void) fourn(zn_FT-1,nn-1,ndim,isign);

  for (k=0; k<nbr_lin*nbr_col; k++)
    {
    i = k/nbr_col;
    j = k-i*nbr_col;
    zn[i][j] = zn_FT[2*k]/(nbr_lin*nbr_col);
    }

  /********  projection sur des espaces multiresolution  ***********/
  err = (int) decomposition_mallat_ncs(rn,nbr_lin,nbr_col,d_order);
  if (err != 0)
     {
     fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
     exit(-1);
     }
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
      rn[i][j] *= V[i][j];
  err = (int) recomposition_mallat_ncs(rn,nbr_lin,nbr_col,d_order);
  if (err != 0)
     {
     fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
     exit(-1);
     }

 
  norm_rn = 0.;
  ps_dn_zn = 0.;
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    {
    ps_dn_zn += rn[i][j]*zn[i][j];
    norm_rn += rn[i][j]*rn[i][j];
    }
  omega_n = norm_rn/ps_dn_zn;

  norm_phi_n = norm_dphi_n = 0.;
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    {
    phi1 = phin[i][j];
    phin1[i][j] = phin[i][j]+omega_n*phin1[i][j];
    phin[i][j] += omega_n*rn[i][j];
    phi2 = phin[i][j];
    norm_phi_n += phi2*phi2;
    norm_dphi_n += (phi2-phi1)*(phi2-phi1);
    }
  epsilon = norm_dphi_n/norm_phi_n;
  (void) fprintf(stderr,"it = %d  eps = %.10f\n",iteration,epsilon);
  }

(void) fprintf(stderr,"  Number of iterations: %d\n\n",iteration);
(void) free(phit_FT);
(void) free(rn_FT);
(void) free(zn_FT);
(void) free(Gr);
(void) free_square_float(rn,nbr_lin);
(void) free_square_float(zn,nbr_lin);
(void) free_square_float(V,nbr_lin);


/*      STORAGE OF THE FINAL SOLUTION
*/
dim = 2;
size[0] = nbr_lin;
size[1] = nbr_col;
strcpy(mode,"bin");
sprintf(file,"%s.PHIN",argv[1]);
write_data(file,phin,dim,size,type,mode,nature,comments);
 
sprintf(file,"%s.PHIN1",argv[1]);
write_data(file,phin1,dim,size,type,mode,nature,comments);

(void) free_square_float(phin,nbr_lin);
fprintf(stderr,"\n\n");
}
