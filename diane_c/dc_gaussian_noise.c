/*************************************************************/
/*                                                           */
/*           simulation d'un bruit gaussien                  */
/*                                                           */
/*************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>


extern  void      write_data();
extern  float  **alloc_square_float();


int nb_ligne = 512;
int nb_colonne = 512;
int Niteration = 10;

void main()

{

/*   declarations            
*/

int i,j,k;
float som;
auto float **image;
auto int size[4];
auto char comments[1024];

srand(9);
image = alloc_square_float( nb_ligne,nb_colonne);

for ( i=0; i<nb_ligne; i++ )
for ( j=0; j<nb_colonne; j++ )
    {
  for ( k=0; k<Niteration; k++ ) som += random();
  som = ((som/Niteration)-0.5)*sqrt(Niteration);
  image[i][j]=som;
    }
size[0]=nb_ligne;
size[1]=nb_colonne;
write_data( "bruit_gaus",image,2,size,"float","bin","real"," bruit gaussien ");
}


