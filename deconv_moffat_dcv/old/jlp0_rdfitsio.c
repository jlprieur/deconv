/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp0_rdfitsio.c
* To read FITS formatted 1-D, 2-D, and 3-D image files
* Using "FITSIO" C package.
* Formats supported are : FITS 8,16,32,-32,-64 
* (i.e. 1-byte, 2-byte or 4-byte integer, and 4-byte or 8-byte float values)
* JLP: comments and jlp_descriptors
*
* Contains:
* int JLP_VM_RDFITS(pntr_array,nx1,ny1,infile,comments,jlp_descr,dflag,istatus)
* int JLP_VM_RDFITS_3D(pntr_array,nx1,ny1,nz1,infile,comments,jlp_descr,
*                      dflag,istatus)
* int JLP_RDFITS(array,nx1,ny1,idim,infile,comments,jlp_descr,dflag,istatus)
* int jlp0_rdfits(pntr_array,array,nx1,ny1,nz1,idim,
*                infile,comments,jlp_descr,dflag,istatus,vm_flag)
* (and main program to test...)
*
* JLP
* Version 25-02-99
---------------------------------------------------------------------*/
#include   <jlp_ftoc.h>
#include   <../fitsio/fitsio.h>

/*
#define DEBUG 1
*/

static int jlp0_rdfits(INT_PNTR *pntr_array, float *array, INT4 *nx1, INT4 *ny1, 
                INT4 *nz1, INT4 *idim, char *infile, char *comments,
                char *jlp_descr, INT4 *dflag, INT4 *istatus, INT4 *vm_flag);

/* Main program to test JLP_RDFITS */
#ifdef TEST_PROGRAM
main()
{
  float        array[1000000];
  INT4         nx1, ny1, idim;
  char         infile[60], outfile[60];
  char         comments[81], jlp_descr[1024];
  INT4         istatus, dflag;
  register int i;

JLP_INQUIFMT();

  printf(" Test of jlp0_rdfits to read FITS files on disk\n");
  printf(" Version 26-01-99\n"); 

      printf(" Input FITS file   : ");
      scanf("%s",infile);
      printf(" Output MIDAS/FITS file : ");
      scanf("%s",outfile);

#if DEBUG
    printf(" Input FITS file   : >%s< \n",infile);
    printf(" Output FITSfile : >%s< \n",outfile);
#endif

idim = 1000; dflag = -1;
JLP_RDFITS(array,&nx1,&ny1,&idim,infile,comments,jlp_descr,&dflag,&istatus);
#if DEBUG
 printf(" JLP_RDFITS/istatus = %d \n",istatus);
#endif

#if DEBUG
 printf(" nx = %d, ny = %d \n",nx1,ny1);
 printf(" comments: %s \n",comments);
 printf(" image[0...20]: \n");
 for(i = 0; i < 20; i++) printf(" %f ",array[i]);
#endif

 JLP_WRITEIMAG(array,&nx1,&ny1,&idim,outfile,comments);
#if DEBUG
 printf(" JLP_WRFITS/istatus = %d \n",istatus);
#endif

JLP_END();
}
#endif

#define SIMPLE_VERSION
#ifdef SIMPLE_VERSION
int JLP_VM_READIMAG1(INT_PNTR *pntr_array, INT4 *nx1, INT4 *ny1, char *infile,
                     char *comments)
{
char jlp_descr[1024];
INT4 dflag, istatus;
dflag = 0;
JLP_VM_RDFITS(pntr_array, nx1, ny1, infile, comments, 
              jlp_descr, &dflag, &istatus);
return(0);
}
#endif

/**********************************************************************
* JLP_VM_RDFITS
*
* dflag = -1 (no error (warning) messages and no descriptors) 
*       = 0 (no descriptors) 
*       = 1 (descriptors) 
**********************************************************************/
int JLP_VM_RDFITS(INT_PNTR *pntr_array, INT4 *nx1, INT4 *ny1, char *infile,
                  char *comments, char *jlp_descr, INT4 *dflag, INT4 *istatus)
{
INT4 vm_flag, idim, nz;
float *array;
int istat;

vm_flag = 1; nz=0;
istat = jlp0_rdfits(pntr_array,array,nx1,ny1,&nz,&idim,
    infile,comments,jlp_descr,dflag,istatus,&vm_flag);

return(istat);
}
/**********************************************************************
* JLP_RDFITS
*
* dflag = -1 (no error (warning) messages and no descriptors) 
*       = 0 (no descriptors) 
*       = 1 (descriptors) 
**********************************************************************/
int JLP_RDFITS(float *array, INT4 *nx1, INT4 *ny1, INT4 *idim1,
               char *infile, char *comments, char *jlp_descr,
               INT4 *dflag, INT4 *istatus)
{
int istat;
INT4 vm_flag, nz;
INT_PNTR pntr_array;

vm_flag = 0; nz=0;
istat = jlp0_rdfits(&pntr_array,array,nx1,ny1,&nz,idim1,
    infile,comments,jlp_descr,dflag,istatus,&vm_flag);

return(istat);
}
/**********************************************************************
* JLP_VM_RDFITS_3D
*
* dflag = -1 (no error (warning) messages and no descriptors)
*       = 0 (no descriptors)
*       = 1 (descriptors)
**********************************************************************/
int JLP_VM_RDFITS_3D(INT_PNTR *pntr_array, INT4 *nx1, INT4 *ny1, INT4 *nz1, 
                  char *infile,
                  char *comments, char *jlp_descr, INT4 *dflag, INT4 *istatus)
{
INT4 vm_flag, idim, nz;
float *array;
int istat;

vm_flag = 1; 
istat = jlp0_rdfits(pntr_array,array,nx1,ny1,nz1,&idim,
    infile,comments,jlp_descr,dflag,istatus,&vm_flag);

return(istat);
}
/**********************************************************************
* jlp0_rdfits
*
* dflag = -1 (no error (warning) messages and no descriptors) 
*       = 0 (no descriptors) 
*       = 1 (descriptors output to screen) 
**********************************************************************/
static int jlp0_rdfits(INT_PNTR *pntr_array, float *array, INT4 *nx1, INT4 *ny1, 
                INT4 *nz1, INT4 *idim, char *infile, char *comments,
                char *jlp_descr, INT4 *dflag, INT4 *istatus, INT4 *vm_flag)
{
 char         filename[61], *pcc, err_message[81], buffer[81];
 int          fmt,num,type,nval,dsize,nz;
 char         devt;
 int          maxdim = 3, simple, bitpix, naxis, istat, extend, any_null_values;
 long         naxes[3], pcount, gcount, nelements;
 float        *array1;
 register int i, j;
 fitsfile *fptr;

*istatus = 0;

   strncpy(filename,infile,40);
   filename[40]='\0';
/* Check input characters strings (pb if fortran...) */
   pcc = filename;
   while(*pcc && *pcc != ' ') pcc++;
   if(*pcc == ' ') *pcc = '\0';

/* Check filename syntax and add ".fits" if no extension is present: */
   pcc = filename;
   while(*pcc && *pcc != '.') pcc++;
   if(*pcc != '.') strcpy(pcc,".fits");
   
   istat = 0;
   fits_open_file(&fptr,filename,READONLY,&istat);
   if (istat) {
     fits_read_errmsg(err_message);
     printf("jlp0_rdfits/ Cannot open input file : >%s< istat=%d\n %s \n",
            filename,istat,err_message);
     *istatus = -1; return(-1);
   }                     /* open file  */

/* decode header    */
    fits_read_imghdr(fptr, maxdim, &simple, &bitpix, &naxis, naxes, &pcount,
                     &gcount, &extend, &istat); 
    if (istat) {   
          printf("NOT supported FITS format! istat = %d\n",istat);
          *istatus = -2;  fits_close_file(fptr,&istat); return(-1);
        }

/* Try to fill "comments" with many keywords */

/* fits_read_key needs istat=0 in input ! */
    istat = 0;
    fits_read_key(fptr,TSTRING,"COMMENT",comments,buffer,&istat);
    istat = 0;
/* try from the beginning each time: */
    fits_read_record(fptr,0,buffer,&istat);
    istat = 0;
    fits_read_key(fptr,TSTRING,"DESCRIP",comments,buffer,&istat);
#ifdef DEBUG
printf(" DEBUG/ istat=%d DESCRIP = %s \n",istat,comments);
#endif
    if(istat == KEY_NO_EXIST) 
      {
      istat = 0;
      fits_read_record(fptr,0,buffer,&istat);
      istat = 0;
      fits_read_key(fptr,TSTRING,"OBJECT",comments,buffer,&istat);
#ifdef DEBUG
      printf(" DEBUG/ istat=%d OBJECT = %s \n",istat,comments);
#endif
      if(istat == KEY_NO_EXIST) 
         {
         istat = 0;
         fits_read_record(fptr,0,buffer,&istat);
         istat = 0;
         fits_read_key(fptr,TSTRING,"DESCRIP",comments,buffer,&istat);
#ifdef DEBUG
         printf(" DEBUG/ istat=%d DESCRIP = %s \n",istat,comments);
#endif
         if(istat == KEY_NO_EXIST) 
             {
             istat = 0;
             fits_read_record(fptr,0,buffer,&istat);
             istat = 0;
             fits_read_key(fptr,TSTRING,"HISTORY",comments,buffer,&istat);
             }
         }
      }

#ifdef DEBUG
printf("comments: >%s<\n",comments);
#endif

/* Descriptors if present */
    istat = 0;
    fits_read_record(fptr,0,buffer,&istat);
    istat = 0;
    fits_read_key(fptr,TSTRING,"HISTORY",jlp_descr,buffer,&istat);
/* Old programs:  jlp_descr + idescr * 62  */

/* Axes values: */
   *nx1 = naxes[0];
   if(naxis > 1 ) *ny1 = naxes[1];
/* In the case of 3-D images, keep the old value of nz, otherwise set it to 0: */
   if(naxis > 2 ) 
        *nz1 = naxes[2]; 
      else 
        *nz1 = 0;

/* Correction for one-dimensional images, if ny1 = 0 */
   if(*nx1 <= 0)
       {
          printf(" JLP_RDFITS/Fatal error: nx = %d \n",*nx1);
          *istatus = -2;  fits_close_file(fptr,&istat); return(-2);
       }

   nelements = naxes[0];
   if(naxis > 1) nelements *= naxes[1];
/* 3-D: */
   if(naxis > 2) nelements *= naxes[2];

/* Allocate memory space if needed: */
    if(*vm_flag)
    {
    dsize = nelements * sizeof(float);
    if((array1 = (float *) malloc(dsize)) == NULL)
     {
     printf("jlp0_rdfits/fatal error allocating memory, nel=%ld\n",nelements);
     *istatus = -2;  fits_close_file(fptr,&istat); return(-2);
     }
    }
    else
     array1 = array;

#ifdef DEBUG
   printf(" nx = %d, ny = %d, nz = %d nelements=%d \n",
         (int)*nx1, (int)*ny1, (int)*nz1, (int)nelements);
#endif

   istat = 0;

/* Check if size of image is consistent with buffer: */
   if(!(*vm_flag))
   {
     if(*nx1 > *idim)
      {
       printf(" JLP_RDFITS/Fatal error: input idim=%d is smaller than nx=%d \n",
               *idim, *nx1);
       *istatus = -2;
      }
   }

/* Read data:    */
    fits_read_img(fptr, TFLOAT, 1, nelements, 0, array1, 
                  &any_null_values, &istat); 
    if (istat) {
      fits_read_errmsg(err_message);
      printf("jlp0_rdfits/fits_read_img; error reading file : >%s<\
 istat=%d \n %s \n",
            filename,istat, err_message);
    *istatus = -3;
    }                     

    fits_close_file(fptr,&istat);

/* Check if idim is equal to nx1:
* if not transfer to good location */
   if(!(*vm_flag))
   {
   if(*idim != *nx1)
     {
       for(j = *ny1 -1; j >= 0; j--)
        for(i = 0; i < *nx1; i++)
            {
             array1[i + j * (*idim)] = array1[i + j * (*nx1)];
            } 
     }
   }
*pntr_array = (INT_PNTR)array1;
return(*istatus);
}
