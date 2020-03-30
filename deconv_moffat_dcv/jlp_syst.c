#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <jlp_ftoc.h>

/***************************************************************
* Use standard system routine:
***************************************************************/
int JLP_SYSTEM(char *string)
{
int status;
status = system(string);
return(status);
}
/*****************************************************************
* JLP_GETENV reads a symbol predefined with setenv :
* Problems since "runs" opens up a new session each time and therefore takes
* the values defined in ".cshrc" ...
*
* New version
******************************************************************/
 int JLP_GETENV(char *symbol, INT4 *length, char *string, INT4 *istat)
 {
   long int nochar, status;
   int len, jlp_system();
   char *getenv(), *pc;
   char mysymbol[41];
   char command[90], filename[31], pbuf[60];
   FILE *fp;

   len = *length;
   strncpy(mysymbol,symbol,len);
   mysymbol[len]='\0';

/* Removes blanks at the end: */
   pc = mysymbol; while(*pc && *pc !=' ')pc++; *pc = '\0';

   strcpy(filename,"jlp_symbol.tmp");
   sprintf(command,"printenv %s > %s",mysymbol,filename);
#if DEBUG
   printf("jlp_getenv: length = %d, mysymbol: >%s< \n",*length,mysymbol);
   printf("command: >%s< \n",command);
#endif

/* Actual command to the system  (be carefull with IBM to link
fortran programs with /lib/libc.a*/
   status = jlp_system(command);

/* Doesn't output an error message if error, since this simply
means the symbol is not there. The non-zero status should be handled
in the calling routine. */
/*
   if(status) {printf("jlp_getenv/ error in 'jlp_system' :\n >%s< \n",
                 command);
               *istat = -2; return(-2);}
*/
   if(status) { *istat = -2; return(-2);}

/* Open file in read mode (0): */
 if((fp = fopen(filename, "r")) == NULL)
 {
  *istat = -1;
  printf("JLP_GETENV/error opening file %s \n",filename);
  return(-1);
  }

/* Read it: */
 nochar = 60;
  status = fread(pbuf, sizeof(char), nochar, fp);
  if(status == 0)
  {*istat = 1; return(1);}
  else
  nochar = status;

#ifdef DEBUG
  printf("nochar %d\n",nochar);
  printf(" Buffer: \n >%s<\n",pbuf);
#endif

/* Note that for Unix, the last character is always '\n' (end of line) */
  nochar--;
  strncpy(string,pbuf,nochar);

/* Close it: */
 fclose(fp);

/* Removes the temporary file: */
 sprintf(command,"rm -f %s",filename);
#if DEBUG
   printf("command: >%s< \n",command);
#endif
 status = jlp_system(command);

   *istat = 0;
#if DEBUG
   printf("string: >%s< \n",string);
#endif
   return(0);
}
/*****************************************************************
* JLP_CTIME  to get the date and time:
******************************************************************/
int JLP_CTIME(char *string, INT4 *istat)
{
   long clock;
   char *p;
   *istat=1;
     time(&clock);
     if( (p=ctime(&clock)) != NULL)
          {*istat=0; strncpy(string,p,20);}
   return(0);
}
