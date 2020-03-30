/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp0_wrfits.c
* To write FITS formatted 2-D image files
* Formats supported are : FITS 32,-32 
* i.e., 4-byte integer (JLP_FORMAT=7) or 4-byte float values (JLP_FORMAT=8)
* JLP comments and jlp_descriptors
*
* June 1996: I remove as many calls as possible to MIDAS include files
* since it can change with MIDAS versions...
*
* Contains:
* int JLP_WRFITS(array,nx1,ny1,idim,filename,comments,jlp_descr,dflag,
*               istatus,out_type)
* (and main program to test...)
*
* JLP
* Version 27-06-96
---------------------------------------------------------------------*/
#include   <stdio.h>
#include   <jlp_ftoc.h>
#include   <../fitsio/fitsio.h>

/*
#define DEBUG 1
*/

static int jlp_wdescr_fits(char *jlp_descr, fitsfile *fptr);

/* Main program to test JLP_WRFITS */
#ifdef TEST_PROGRAM
main()
{
  float        array[1000000];
  INT4         nx1, ny1, idim1;
  char         infile[60], outfile[60];
  char         comments[81], jlp_descr[1024];
  INT4         istatus, dflag, out_type;

JLP_BEGIN();
JLP_INQUIFMT();

  printf(" jlp0_wrfits to convert MIDAS bdf files to FITS format\n");
  printf(" Version 11-11-91\n");

      printf(" Input fits file : ");
      scanf("%s",infile);
      printf(" Output file : ");
      scanf("%s",outfile);

#if DEBUG
    printf(" Input fits file : >%s< \n",infile);
    printf(" Output file     : >%s< \n",outfile);
#endif

idim1 = 1000; dflag = 0;
JLP_READIMAG(array,&nx1,&ny1,&idim1,infile,comments);
strcpy(jlp_descr,
      " It was really a long story isn't it. Although I thought it could be much shorter I couldn't help talking over and over...");

/*JLP_RDFITS(array,&nx1,&ny1,&idim1,infile,comments,jlp_descr,&dflag,&istatus);
*/
#if DEBUG
 printf(" JLP_RDFITS/istatus = %ld \n",istatus);
#endif

#if DEBUG
 printf(" nx = %ld, ny = %ld \n",nx1,ny1);
 printf(" comments: %s \n",comments);
#endif

out_type = 0;
JLP_WRFITS(array,&nx1,&ny1,&idim1,outfile,comments,jlp_descr,
           &dflag,&istatus,&out_type);
#if DEBUG
 printf(" JLP_WRFITS/istatus = %ld \n",istatus);
#endif

JLP_END();
}
#endif
#define SIMPLE_VERSION
#ifdef SIMPLE_VERSION
int JLP_WRITEIMAG(float *array, INT4 *nx1, INT4 *ny1, INT4 *idim1,
                  char *filename, char *comments)
{
char jlp_descr[1]; 
INT4 dflag, istatus, out_type;
out_type = 1;
dflag = 0;
jlp_descr[0] = '\0';
JLP_WRFITS(array, nx1, ny1, idim1, filename, comments, jlp_descr, &dflag,
           &istatus, &out_type);
return(0);
}
int JLP_D_WRITEIMAG(double *d_array, INT4 *nx1, INT4 *ny1, INT4 *idim1,
                  char *filename, char *comments)
{
char jlp_descr[1]; 
INT4 dflag, istatus, out_type;
float *f_array;
register int i, j;

if((f_array = (float *) malloc((*nx1) * (*ny1) * sizeof(float))) == NULL)
  {
   printf("JLP_D_WRITEIMAG/Fatal error allocating memory space: nx=%d ny1=%dn",
           *nx1, *ny1);
   exit(-1);
  }
for(j = 0; j < *ny1; j++) {
  for(i = 0; i < *nx1; i++) {
    f_array[i + j * (*nx1)] = d_array[i + j * (*idim1)];
  }
 }
out_type = 1;
dflag = 0;
jlp_descr[0] = '\0';
JLP_WRFITS(f_array, nx1, ny1, nx1, filename, comments, jlp_descr, &dflag,
           &istatus, &out_type);
free(f_array);
return(0);
}
#endif

/**********************************************************************
* JLP_WRFITS
*
* Output format:
* out_type = 0 : integer (32)
* out_type = 1 : real (-32) 
**********************************************************************/
int JLP_WRFITS(float *array, INT4 *nx1, INT4 *ny1, INT4 *idim1,
               char *filename, char *comments, char *jlp_descr, INT4 *dflag,
               INT4 *istatus, INT4 *out_type)
{
  char         *pcc, lhcuts[32], command[80], err_message[81];
  long int     naxes[2], nelements, nx, ny, idim; 
  float        work, lcut, hcut, *tmp_array;
  int          naxis, istat;
/* *fptr = pointer to FITS file, defined in fitsio.h */
  fitsfile *fptr;      
  register int i, j;

#if DEBUG
   printf("Beginning of JLP_WRFITS out_type=%d\n",*out_type);
   printf("Output file >%s< nx=%ld ny=%ld\n",filename,*nx1,*ny1);
#endif

*istatus = 0;

/* Check input characters strings (pb if fortran...) */
   pcc = filename;
   while(*pcc && *pcc != ' ') pcc++;
   if(*pcc == ' ') *pcc = '\0';

   comments[79]='\0';

/* Transfer to long int variables (necessary for Fortran interface with OSF1)*/
nx = *nx1; ny = *ny1; idim = *idim1;

/* Check filename syntax and add ".fits" if no extension is present: */
   pcc = filename;
   while(*pcc && *pcc != '.') pcc++;
   if(*pcc != '.') strcpy(pcc,".fits");

#if DEBUG
   printf("Opening output file >%s< nx=%ld ny=%ld idim=%ld\n",
          filename,nx,ny,idim);
#endif
   istat = 0;
   fits_create_file(&fptr, filename, &istat);
/* JLP99: overwrite file if already there: */
   if(istat)
    {
    printf("File already here: I delete it \n");
    sprintf(command,"rm %s",filename);
    jlp_system(command);
    istat = 0;
    fits_create_file(&fptr, filename, &istat);
    if(istat)
     {
     fits_read_errmsg(err_message);
     printf("JLP_WRFITS/Error when opening FITS file %s, istat=%d\n %s\n",
            filename, istat, err_message);
     *istatus = -1;
     return(-1);
     }
    }

/* compute max/min cuts of data */
  lcut = *array; hcut = lcut;
  for(j=0; j < ny; j++)
   {for(i=0; i < nx; i++)
      { 
        work = array[i + j * idim];
          if ( work < lcut) lcut = work;
           else if (hcut < work) hcut = work;
      }
   }

/* write FITS header      */
  naxis=2;
  naxes[0] = nx;
  naxes[1] = ny;
/* out_type = 0 : integer (32)
 out_type = 1 : real (-32) */ 
 if(out_type)
   fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &istat);
 else
   fits_create_img(fptr, LONG_IMG, naxis, naxes, &istat);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when writing header, istat=%d\n %s\n",
           istat, err_message);
   }

/* "LHCUTS  " */
  sprintf(lhcuts,"%12.5e %12.5e",lcut, hcut);
  fits_update_key(fptr, TSTRING, "LHCUTS  ", lhcuts, "Low and high cuts", &istat);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when writing keyword LHCUTS, istat=%d\n %s\n %s\n",
           istat, lhcuts, err_message);
   }

/* Copy comments to OBJECT for compatibility with other software*/
/* "OBJECT " */
  fits_update_key(fptr, TSTRING, "OBJECT ", comments, "", &istat);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when writing keyword OBJECT, istat=%d\n %s\n %s\n",
           istat, comments, err_message);
   }
  fits_update_key(fptr, TSTRING, "OBJECT ", comments, "", &istat);

/* Date */
  fits_write_date(fptr, &istat);

/* Copy jlp_descr to JLPDESCR */
/* "JLPDESCR " */
  if(*jlp_descr) jlp_wdescr_fits(jlp_descr, fptr);

/* "AUTHOR  " */
  fits_update_key(fptr, TSTRING, "AUTHOR  ", "Jean-Louis Prieur", 
                  "With cfitsio 1.4", &istat);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when writing keyword AUTHOR, istat=%d\n %s\n",
           istat, err_message);
   }


/* write prime data       */
  nelements = nx * ny;
  if((tmp_array = (float*) malloc(nelements*sizeof(float))) == NULL)
  { printf("JLP_WRFITS/Error when allocating memory for temporary array nel=%ld\n", 
            nelements); 
    return(-2);
  }
  for(j=0; j < ny; j++)
   for(i=0; i < nx; i++)
      tmp_array[i + j * nx] = array[i + j * idim]; 

/* Write the array to the image: (automatic conversion to the correct type,
 but should give the type of "tmp_array")*/
  fits_write_img(fptr, TFLOAT, 1, nelements, tmp_array, &istat);
  free(tmp_array);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when writing data to file, istat=%d\n %s\n",
           istat, err_message);
   }
  
/* Close file on disk: */
  fits_close_file(fptr, &istat);
  if(istat)
   {
   fits_read_errmsg(err_message);
   printf("JLP_WRFITS/Error when closing file %s, istat=%d\n %s\n",
           filename, istat, err_message);
   }

return(0); 
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Routine to write jlp_descriptors in FITS format
*
* JLP
---------------------------------------------------------------------*/
static int jlp_wdescr_fits(char *jlp_descr, fitsfile *fptr)
{
 register int i, k;
 int istat;
 char buffer[81], mydescr[1024], *pc, err_message[81];

  strcpy(mydescr,jlp_descr);

/********* First fills the end of descriptors with zeroes ********/
/* Look for the first zero: */
  mydescr[1023] = '\0';
  pc = mydescr;
  while(*pc) pc++;

/* Then replaces the blanks at the end with zeroes */
  pc--;
  while(*pc == ' ') {*pc = '\0'; pc--;}

/***** Then copies the descriptors to FITS header: (1024/62=17 lines maxi) ****/
  for (k=0; k<17; k++)
  {
  strncpy(buffer,&mydescr[k*62],62);
  buffer[62] = '\0';
#ifdef DEBUG
  printf(" Writing keyword JLPDESCR \n %s\n",buffer);
#endif
     istat = 0;
     fits_write_key_str(fptr, "JLPDESCR", buffer, "",&istat);
     if(istat)
      {
      fits_read_errmsg(err_message);
      printf("JLP_WRFITS/Error when writing keyword JLPDESCR, istat=%d \n %s\n %s\n",
              istat, buffer, err_message);
      }
/* Stop when first zero has been found */
    for (i=0; i<62; i++) if(buffer[i] == '\0') break;
    if(i < 62) break;
  }

/* end: */
 return(istat);
}
