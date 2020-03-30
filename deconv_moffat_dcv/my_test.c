#include <stdio.h>
#include <math.h>

float calculer(float (*pfonc)(float, float), float x,  float y);
float somme(float x, float y);

void main()
{
float a,b;
a=1.02; b=5.5;
printf(" somme=%f\n",calculer(somme, a, b));
}

float calculer(float (*pfonc)(float, float), float x, float y)
{
 return (*pfonc)(x,y);
}

float somme(float x, float y)
{
 return(x+y);
}
