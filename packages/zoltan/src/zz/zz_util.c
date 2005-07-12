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
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "zz_util_const.h"
#include "zoltan_mem.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Remove leading & trailing white space and convert to upper case. */

int Zoltan_Clean_String(
char *string1,			/* original string */
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

/* 
 * Input: the bounds of an axis-aligned box, a linear transformation, and a
 * count of dimensions (1, 2 or 3),
 *
 * Output: the 8 vertices of the box obtained by applying the transformation
 * followed by:
 *  ndims is 1: projecting to the X-axis
 *  ndims is 2: projecting to the XY plane
 *  ndims is 3: no projection
 *
 * This is used by the box assign functions when coordinates have been
 * transformed due to degenerate geometries.
 */

void Zoltan_Transform_Box_Points(
  double *lo, double *hi,    /* input: bounds of 3D axis-aligned box */
  double (*m)[3],            /* input: 3x3 linear transformation     */
  int ndims,                 /* input: 1, 2 or 3 dimensions in result */
  double (*v)[3])            /* output: 8 vertices of resulting box    */
{
     double temp[3];
     int i, j;

     v[0][0] = lo[0]; v[0][1] = lo[1]; v[0][2] = lo[2];
     v[1][0] = lo[0]; v[1][1] = hi[1]; v[1][2] = lo[2];
     v[2][0] = lo[0]; v[2][1] = hi[1]; v[2][2] = hi[2];
     v[3][0] = lo[0]; v[3][1] = lo[1]; v[3][2] = hi[2];
     v[4][0] = hi[0]; v[4][1] = lo[1]; v[4][2] = lo[2];
     v[5][0] = hi[0]; v[5][1] = hi[1]; v[5][2] = lo[2];
     v[6][0] = hi[0]; v[6][1] = hi[1]; v[6][2] = hi[2];
     v[7][0] = hi[0]; v[7][1] = lo[1]; v[7][2] = hi[2];

     for (i=0; i<8; i++){
       for (j=0; j<3; j++){
         temp[j] = v[i][j];
       }
       v[i][0] = m[0][0]*temp[0] + m[0][1]*temp[1] + m[0][2]*temp[2];
       if (ndims > 1 ){
         v[i][1] = m[1][0]*temp[0] + m[1][1]*temp[1] + m[1][2]*temp[2];
         if (ndims > 2){
           v[i][2] = m[2][0]*temp[0] + m[2][1]*temp[1] + m[2][2]*temp[2];
         }
         else{
           v[i][2] = 0.0;
         }
       }
       else{
         v[i][1] = 0.0;
         v[i][2] = 0.0;
       }
     }
}

/* 
 * Given the bounds of a 3D axis-aligned box and a 3x3 linear
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
  int ndims)               /* 1, 2 or 3 dimensional result */
{          
     double v[8][3];
     int i;

     Zoltan_Transform_Box_Points(lo, hi, m, ndims, v);

     lo[0] = hi[0] = v[0][0];

     if (ndims > 1){
       lo[1] = hi[1] = v[0][1];
       if (ndims > 2){
         lo[2] = hi[2] = v[0][2];
       }
     }

     for (i=1; i<8; i++){
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
