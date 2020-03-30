/*****************************************************************************/
/*                                                                           */
/* FUNCTION: recomposition_mallat_cs                                         */
/*                                                                           */
/* PURPOSE: rebuilts an image with Mallat algorithm using cubic splines.     */
/*          the decomposed image is supposed to be not normalized            */
/*                                                                           */
/* INPUT:  mat[0..nbr_lin-1][0..nbr_col-1] = image to be recomposed          */
/*         nbr_lin,nbr_col = dimensions of the image                         */
/*         order = order of decomposition                                    */
/*                                                                           */
/* OUTPUT: mat[0..nbr_lin-1][0..nbr_col-1] = image recomposed                */
/*                                                                           */
/* RETURN: 0 if every thing OK                                               */
/*         -1 if memory allocation failure                                   */
/*                                                                           */
/* VERSION: November 1992                                                    */
/*                                                                           */
/* AUTHORS: Karim BOUYOUCEF and Eric ANTERRIEU                               */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <malloc.h>
#define  NPT    41
#define  NPTS2  20   /*  (NPT-1)/2  */
/*****************************************************************************/
int     recomposition_mallat_ncs(mat,nbr_lin,nbr_col,order)
int     nbr_lin, nbr_col, order;
float   **mat;
/*****************************************************************************/
{
register  int    i, j, ii, jj, f, k;
auto      int    dim, iter, fmin;
auto      int    lin, lins2, imin, imax, i2;
auto      int    col, cols2, jmin, jmax, j2;
auto      float  *xvect, *yvect, H[NPT], G[NPT], yy, zz;

dim = nbr_lin;  if (nbr_col > nbr_lin) dim = nbr_col;
xvect = (float *) malloc(dim*sizeof(float));
if (xvect == NULL) return(-1);
yvect = (float *) malloc(dim*sizeof(float));
if (yvect == NULL) return(-1);

/* h filter
*/
H[0]  =  0.000000; H[1]  = -0.000164; H[2]  = -0.000202; H[3]  =  0.000327;
H[4]  =  0.000396; H[5]  = -0.000656; H[6]  = -0.000780; H[7]  =  0.001331;
H[8]  =  0.001546; H[9]  = -0.002745; H[10] = -0.003079; H[11] =  0.005799;
H[12] =  0.006141; H[13] = -0.012715; H[14] = -0.012145; H[15] =  0.029747;
H[16] =  0.022685; H[17] = -0.077808; H[18] = -0.035498; H[19] =  0.306830;
H[20] =  0.541736; H[21] =  0.306830; H[22] = -0.035498; H[23] = -0.077808;
H[24] =  0.022685; H[25] =  0.029747; H[26] = -0.012145; H[27] = -0.012715;
H[28] =  0.006141; H[29] =  0.005799; H[30] = -0.003079; H[31] = -0.002745;
H[32] =  0.001546; H[33] =  0.001331; H[34] = -0.000780; H[35] = -0.000656;
H[36] =  0.000396; H[37] =  0.000327; H[38] = -0.000202; H[39] = -0.000164;
H[40] =  0.000000;

/* g filter
*/
G[0]  =  0.000000; G[1]  =  0.000000; G[2]  =  0.000164; G[3]  = -0.000202;
G[4]  = -0.000327; G[5]  =  0.000396; G[6]  =  0.000656; G[7]  = -0.000780;
G[8]  = -0.001331; G[9]  =  0.001546; G[10] =  0.002745; G[11] = -0.003079;
G[12] = -0.005799; G[13] =  0.006141; G[14] =  0.012715; G[15] = -0.012145;
G[16] = -0.029747; G[17] =  0.022685; G[18] =  0.077808; G[19] = -0.035498;
G[20] = -0.306830; G[21] =  0.541736; G[22] = -0.306830; G[23] = -0.035498;
G[24] =  0.077808; G[25] =  0.022685; G[26] = -0.029747; G[27] = -0.012145;
G[28] =  0.012715; G[29] =  0.006141; G[30] = -0.005799; G[31] = -0.003079;
G[32] =  0.002745; G[33] =  0.001546; G[34] = -0.001331; G[35] = -0.000780;
G[36] =  0.000656; G[37] =  0.000396; G[38] = -0.000327; G[39] = -0.000202;
G[40] =  0.000164;

if (order == 0)   return(0);

for (iter=1; iter<=order; iter++)
  {
  lin = nbr_lin >> (order-iter);
  col = nbr_col >> (order-iter);
  lins2 = lin/2;
  cols2 = col/2;

  for (j=0; j<col; j++) /* transformation of columns */
    {
    for (i=0; i<lin; i++) {xvect[i] = 0.0; yvect[i] = mat[i][j];}
    for (i=0; i<lins2; i++)
      {
      i2 = i+i;
      imin = i2 - NPTS2; if (imin < 0)  imin = 0;
      imax = i2 + NPTS2; if (imax > (lin-1)) imax = lin-1;
      f = fmin = NPTS2 + (imin - i2);
      yy = yvect[i];
      zz = yvect[lins2+i];
      for (ii=imin; ii<=imax; ii++)
        {
        xvect[ii] += (yy*H[f] + zz*G[f]);
        f++;
        }
      }
    for (i=0; i<lin; i++)   mat[i][j] = xvect[i];
    }

  for (i=0; i<lin; i++) /* transformation of lines */
    {
    for (j=0; j<col; j++)  {xvect[j]=0.; yvect[j] = mat[i][j];}
    for (j=0; j<cols2; j++)
      {
      j2 = j+j;
      jmin = j2 - NPTS2; if (jmin < 0)  jmin = 0;
      jmax = j2 + NPTS2; if (jmax > (col-1)) jmax = col-1;
      f = fmin = NPTS2 + (jmin - j2);
      yy = yvect[j];
      zz = yvect[cols2+j];
      for (jj=jmin; jj<=jmax; jj++)
        {
        xvect[jj] += (yy*H[f] + zz*G[f]);
        f++;
        }
      }
    for (j=0; j<col; j++)   mat[i][j] = 4.0*xvect[j];
    }	
  }

free((char *) xvect);
free((char *) yvect);
return(0);
}
