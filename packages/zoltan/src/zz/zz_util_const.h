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


#ifndef __ZOLTAN_UTIL_CONST_H
#define __ZOLTAN_UTIL_CONST_H

#include "zz_const.h"
#include "zoltan_types.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

extern unsigned int Zoltan_Hash(ZOLTAN_ID_PTR, int, unsigned int);
extern int Zoltan_Clean_String(const char *, char **);
extern char *Zoltan_Strdup(const char *);
void Zoltan_Transform_Point( double *p, double (*m)[3], int *a, int d,
  int ndims, double *v);
void Zoltan_Transform_Box(double *lo, double *hi, double (*m)[3], int *a, 
  int d, int ndims);
void Zoltan_Transform_Box_Points(double *lo, double *hi, double (*m)[3], 
  int *a, int d, int ndims, double (*v)[3]);
int Zoltan_AllReduceInPlace(void *, int , MPI_Datatype , MPI_Op , MPI_Comm );

/* A Zoltan_Map is like a C++ STL map.  It uses Zoltan_Hash.
 */

#define ZOLTAN_MAX_MAP 10 /* The max number of simultaneous sets created */

int Zoltan_Map_Create(ZZ *zz, int hash_range, int num_id_entries, int store_keys, int num_entries);
int Zoltan_Map_Destroy(ZZ *zz, int set_num);
int Zoltan_Map_Add(ZZ *zz, int set_num, int *key, void *data);
int Zoltan_Map_Find(ZZ *zz, int set_num, int *key, void **data);
int Zoltan_Map_Size(ZZ *zz, int set_num);
int Zoltan_Map_First(ZZ *zz, int set_num, int **key, void **data);
int Zoltan_Map_Next(ZZ *zz, int set_num, int **key, void **data);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
