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


int
Zoltan_AllReduceInPlace(void *sndrcvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int ierr;

#ifndef MPI_IN_PLACE
  void * dummy;
  int size;

  MPI_Type_size(datatype, &size);

  dummy = ZOLTAN_MALLOC(size*count);
  if (dummy == NULL)
    return ZOLTAN_MEMERR;
  memcpy (dummy, sndrcvbuf, size*count);
  ierr = MPI_Allreduce(dummy, sndrcvbuf, count, datatype, op, comm);
  ZOLTAN_FREE(&dummy);
#else /* MPI_IN_PLACE */
  ierr = MPI_Allreduce(MPI_IN_PLACE, sndrcvbuf, count, datatype, op, comm);
#endif /* MPI_IN_PLACE */
  return (ierr);
}

MPI_Datatype _mpi_int_32_type;
MPI_Datatype _mpi_int_max_type;
MPI_Datatype _mpi_uint_32_type;
MPI_Datatype _mpi_uint_max_type;

char _z_int_specifier[16];
char _z_uint_specifier[16];
char _z_int_long_specifier[16];
char _z_uint_long_specifier[16];

int Zoltan_set_mpi_types()
{
  int rank;
  MPI_Datatype signed_types[16];
  MPI_Datatype unsigned_types[16];
  int size_short, size_unsigned_short, i;
  int size_int, size_unsigned_int;
  int size_long, size_unsigned_long;
  int size_long_long, size_unsigned_long_long;

  /* specifiers for scanf, printf, etc. */

  if (sizeof(Z_INT_L) == sizeof(int)){
    strcpy(_z_int_long_specifier,"%d");
    strcpy(_z_uint_long_specifier,"%u");
  }
  else if (sizeof(Z_INT_L) == sizeof(long)){
    strcpy(_z_int_long_specifier,"%ld");
    strcpy(_z_uint_long_specifier,"%lu");
  }
  else{   /* long long */
    strcpy(_z_int_long_specifier,"%Ld");
    strcpy(_z_uint_long_specifier,"%Lu");
  }

  if (sizeof(Z_INT) == sizeof(long)){
    strcpy(_z_int_specifier,"%ld");
    strcpy(_z_uint_specifier,"%lu");
  }
  else{   /* int or _int32 */
    strcpy(_z_int_specifier,"%d");
    strcpy(_z_uint_specifier,"%u");
  }

  MPI_Type_size(MPI_SHORT, &size_short);
  MPI_Type_size(MPI_UNSIGNED_SHORT, &size_unsigned_short);
  MPI_Type_size(MPI_INT, &size_int);
  MPI_Type_size(MPI_UNSIGNED, &size_unsigned_int);
  MPI_Type_size(MPI_LONG, &size_long);
  MPI_Type_size(MPI_UNSIGNED_LONG, &size_unsigned_long);
  MPI_Type_size(MPI_LONG_LONG, &size_long_long);
  MPI_Type_size(MPI_UNSIGNED_LONG_LONG, &size_unsigned_long_long);

  for (i=0; i < 16; i++){
    signed_types[i] = MPI_UNDEFINED;
    unsigned_types[i] = MPI_UNDEFINED;
  }

  signed_types[size_long_long] = MPI_LONG_LONG;
  signed_types[size_long] = MPI_LONG;
  signed_types[size_int] = MPI_INT;
  signed_types[size_short] = MPI_SHORT;

  unsigned_types[size_unsigned_long_long] = MPI_UNSIGNED_LONG_LONG;
  unsigned_types[size_unsigned_long] = MPI_UNSIGNED_LONG;
  unsigned_types[size_unsigned_int] = MPI_UNSIGNED;
  unsigned_types[size_unsigned_short] = MPI_UNSIGNED_SHORT;

  _mpi_int_32_type = signed_types[4];
  _mpi_uint_32_type = unsigned_types[4];
  
  if (signed_types[sizeof(Z_INT_L)] != MPI_UNDEFINED){
    _mpi_int_max_type = signed_types[sizeof(Z_INT_L)];
    _mpi_uint_max_type = unsigned_types[sizeof(Z_INT_L)];
  }
  else{
    _mpi_int_max_type = signed_types[4];
    _mpi_uint_max_type = unsigned_types[4];
  }

  if ((_mpi_int_32_type == MPI_UNDEFINED) ||
      (_mpi_uint_32_type == MPI_UNDEFINED) ||
      (_mpi_int_max_type == MPI_UNDEFINED) ||
      (_mpi_uint_max_type == MPI_UNDEFINED) ){

    return ZOLTAN_FATAL;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0){
    printf("Size of a Z_INT_L is %d bytes, Z_INT is %d bytes\n", sizeof(Z_INT_L), sizeof(Z_INT));
  } 

  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
