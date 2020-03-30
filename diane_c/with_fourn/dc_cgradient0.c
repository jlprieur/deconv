/***************************************************************/
/*                                                             */
/* PROGRAM: gradient                                           */
/*                                                             */
/* PURPOSE: Deconvolution with multiresolution constraint      */
/*                                                             */
/* INPUT:  argv[1] = PHIT (first approximation of the solution)*/
/*         argv[2] = G (interpolation function)                */
/*         argv[3] = V (multiresolution support constraint)    */
/*         argv[4] = decomposition order for support constraint*/
/*         argv[5] = threshold for residual                    */
/*                                                             */
/* OUTPUT: argv[6] = clean map file                            */
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
auto     char     type[16], mode[8], nature[16], comments[1024];
auto     float    **phit, **G, **V, **phin;
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
if (argc != 7)
  {
  (void) printf("\a\nUnexpected number of arguments\n");
  (void) printf("\nUSAGE:\n");
  (void) printf("gradient  PHIT  G  V  d_order  eps  PHIN\n\n");
  (void) printf("\t PHIT = first approximation of the solution\n");
  (void) printf("\t G = interpolation function\n");
  (void) printf("\t V = multiresolution support constraint\n");
  (void) printf("\t d_order = decomposition order of the constraint\n");
  (void) printf("\t eps = threshold for residual\n");
  (void) printf("\t PHIN = final solution\n");
  (void) printf("\t \n\n");
  exit(-1);
  }

/*       LECTURE DES DIFFERENTS FICHIERS ET DES DATAs
*/
dim = 2;
(void) strcpy(type,"float");
(void) strcpy(nature,"real");

phit = (float **) read_data(argv[1],dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];

G = (float **) read_data(argv[2],dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s & %s] have different sizes\n",argv[1],argv[2]);
  exit(-1);
  }

V = (float **) read_data(argv[3],dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s & %s] have different sizes\n",argv[1],argv[3]);
  exit(-1);
  }

err = sscanf(argv[4],"%d",&d_order);
if ((err != 1) || (d_order < 0) || (d_order > 10))
  {
  fprintf(stderr,"\a\nFATAL ERROR : incorrect decomposition order [%d]\n",d_order);
  exit(-1);
  }

err = sscanf(argv[5],"%f",&eps);
if ((err != 1) || (eps < 0.000000001) || (eps > 100.0))
  {
  fprintf(stderr,"\a\nFATAL ERROR : incorrect threshold for residual [%f]\n",eps);
  exit(-1);
  }

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
    phin[i][j] += omega_n*rn[i][j];
/*
    if (phin[i][j] < 0.0) phin[i][j] = 0.0;
*/
    phi2 = phin[i][j];
    norm_phi_n += phi2*phi2;
    norm_dphi_n += (phi2-phi1)*(phi2-phi1);
    }
  epsilon = norm_dphi_n/norm_phi_n;
  (void) fprintf(stderr,"it = %d  eps = %f\n",iteration,epsilon);
  }

(void) printf("  Number of iterations: %d\n\n",iteration);
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
write_data(argv[6],phin,dim,size,type,mode,nature,comments);
 
(void) free_square_float(phin,nbr_lin);
}
