
#ifndef _MAP
#define _MAP

#include <stdlib.h>

void mapd ( double , int, double *, int, int );    /* map x to y         */ //����������� ����� x �� ����������� ������� [-0.5,0.5] � ��������� 2^(-m), key = 1
void invmad ( int, double *, int , int *, double *, int , int );  /* map y to x         */
void xyd ( double *, int, double *, int );        /* get preimage       */

#endif
