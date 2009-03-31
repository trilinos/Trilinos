/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <ctype.h>

#include "zz_util_const.h"
#include "zoltan_mem.h"
#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Remove leading & trailing white space and convert to upper case. */

int Zoltan_Clean_String(
const char *string1,			/* original string */
char **pstring2) 		/* cleaned string to return */
{

    char     *string2;		/* cleaned up string */
    int       start, end;	/* indices bounding true string */
    int       length1;		/* length of string 1 */
    int       i;		/* loop counter */

    length1 = strlen(string1);
    start = 0;
    end = length1;
    while (start < length1 && isspace((int)(string1[start])))
	start++;
    while (end > start && isspace((int)(string1[end])))
	end--;

    string2 = (char *) ZOLTAN_MALLOC((end - start + 1) * sizeof(char));
    *pstring2 = string2;

    if (string2 == NULL)
	return (ZOLTAN_MEMERR);

    for (i = start; i < end; i++) {
	*string2++ = toupper(string1[i]);
    }
    *string2 = '\0';

    return (ZOLTAN_OK);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Version of strdup that uses Zoltan_Malloc */

char *Zoltan_Strdup(const char *str)
{
  /* char *yo = "Zoltan_Strdup"; */
  char *c = NULL;
  if (!str){
    return c;
  }

  c = (char *)ZOLTAN_MALLOC(strlen(str) + 1);

  if (c){
    strcpy(c, str);
  }

  return c;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
 /*
  * These transformations are used by the box assign functions when 
  * coordinates have been transformed due to degenerate geometries.
  */
void Zoltan_Transform_Point(
  double *p,         /* point to transform */
  double (*m)[3],    /* linear transformation */
  int *permute,      /* simplified transformation, permute coordinates */
  int d,             /* dimension of input (2 or 3) */
  int ndims,         /* dimension of output (1, 2 or 3) */
  double *v)         /* output */
{
  double tmp[3];

  tmp[0] = p[0];
  tmp[1] = p[1];
  tmp[2] = p[2];

  if (permute[0] >= 0){
    v[0] = p[permute[0]];
    v[1] = (ndims < 2) ? 0.0 : p[permute[1]];

    if (d == 3){
      v[2] = (ndims < 3) ? 0.0 : p[permute[2]];
    }
  }
  else{
    if (d == 2){  
      v[0] = m[0][0] * tmp[0]  +  m[0][1] * tmp[1];
   
      v[1] = (ndims < 2) ? 
              0.0 :
              m[1][0] * tmp[0]  +  m[1][1] * tmp[1];
    }
    else if (d == 3) {
      v[0] = m[0][0]*tmp[0] + m[0][1]*tmp[1] + m[0][2]*tmp[2];
   
      v[1] = (ndims < 2) ? 
              0.0 :
              m[1][0]*tmp[0] + m[1][1]*tmp[1] + m[1][2]*tmp[2];
   
      v[2] = (ndims < 3) ?
              0.0 :
              m[2][0]*tmp[0] + m[2][1]*tmp[1] + m[2][2]*tmp[2];
     
    }
  }
}

/* 
 * Input: the bounds of an axis-aligned box, a linear transformation, the
 * dimension (2 or 3) of the box, and the dimension of the transformed box
 * (1, 2 or 3).
 *
 * Output: the 4 or 8 vertices of the box obtained by applying the transformation
 * followed by:
 *  ndims is 1: projecting to the X-axis
 *  ndims is 2: projecting to the XY plane
 *  ndims is 3: no projection
 */

void Zoltan_Transform_Box_Points(
  double *lo, double *hi,    /* input: bounds of 2D or 3D axis-aligned box */
  double (*m)[3],            /* input: linear transformation     */
  int *perm,                 /* input: coordinate permutation */
  int d,                     /* dimension of box (2 or 3) */
  int ndims,                 /* dimension of transformed box (1, 2 or 3) */
  double (*v)[3])            /* output: 4 or 8 vertices of resulting box    */
{
     int i;

     if (d == 2){  
       v[0][0] = lo[0]; v[0][1] = lo[1];
       v[1][0] = lo[0]; v[1][1] = hi[1];
       v[2][0] = hi[0]; v[2][1] = hi[1]; 
       v[3][0] = hi[0]; v[3][1] = lo[1]; 

       for (i=0; i<4; i++){
         Zoltan_Transform_Point(v[i], m, perm, 2, ndims, v[i]);
       }
     }
     else if (d == 3) {
       v[0][0] = lo[0]; v[0][1] = lo[1]; v[0][2] = lo[2];
       v[1][0] = lo[0]; v[1][1] = hi[1]; v[1][2] = lo[2];
       v[2][0] = lo[0]; v[2][1] = hi[1]; v[2][2] = hi[2];
       v[3][0] = lo[0]; v[3][1] = lo[1]; v[3][2] = hi[2];
       v[4][0] = hi[0]; v[4][1] = lo[1]; v[4][2] = lo[2];
       v[5][0] = hi[0]; v[5][1] = hi[1]; v[5][2] = lo[2];
       v[6][0] = hi[0]; v[6][1] = hi[1]; v[6][2] = hi[2];
       v[7][0] = hi[0]; v[7][1] = lo[1]; v[7][2] = hi[2];

       for (i=0; i<8; i++){
         Zoltan_Transform_Point(v[i], m, perm, 3, ndims, v[i]);
       }
     }
}

/* 
 * Given the bounds of a 2D or 3D axis-aligned box and a linear
 * transformation, return the bounds of the transformed box,  Set
 * ndims as described above if the box is to be projected to a
 * lower dimension.
 *
 * This is used by the box assign functions when coordinates have been
 * transformed due to degenerate geometries.
 */

void Zoltan_Transform_Box( 
  double *lo, double *hi,  /* input: box bounds, output: bounds of transformed box */
  double (*m)[3],          /* 3x3 transformation */
  int *perm,               /* if transformation is simple coordinate permutation */
  int d,                   /* dimension of box */
  int ndims)               /* dimension of transformed box */
{          
     double v[8][3];
     int i, npoints;

     npoints = ((d == 2) ? 4 : 8);

     Zoltan_Transform_Box_Points(lo, hi, m, perm, d, ndims, v);

     lo[0] = hi[0] = v[0][0];
     lo[1] = hi[1] = 0.0;
     lo[2] = hi[2] = 0.0;

     if (ndims > 1){
       lo[1] = hi[1] = v[0][1];
       if (ndims > 2){
         lo[2] = hi[2] = v[0][2];
       }
     }

     for (i=1; i<npoints; i++){
       if (v[i][0] < lo[0] )     lo[0] = v[i][0];
       else if (v[i][0] > hi[0]) hi[0] = v[i][0];
       if (ndims > 1){
         if (v[i][1] < lo[1])      lo[1] = v[i][1];
         else if (v[i][1] > hi[1]) hi[1] = v[i][1];
         if (ndims > 2){
           if (v[i][2] < lo[2])      lo[2] = v[i][2];
           else if (v[i][2] > hi[2]) hi[2] = v[i][2];
         }
       }
     }
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
