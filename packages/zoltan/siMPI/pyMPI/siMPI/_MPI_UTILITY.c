/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************       _MPI_Utility.c     ***********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 26 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

void _MPI_safeFree(void* buffer,char* message) {
  _MPI_COVERAGE();
  if ( buffer ) free(buffer);
}

/*=============================================================================================*/
void *_MPI_safeMalloc(int size, char* message)
{
  void *temp;
  _MPI_COVERAGE();
  temp = malloc(size+128); /* HACK: Padding */
  _MPI_checkMemAlloc(temp,message);
  memset(temp,0,size);
  return temp;
}
/*=============================================================================================*/
void _MPI_checkMemAlloc (void* array, char* message)
{
  _MPI_COVERAGE();
  if (array == (void*) NULL) {
  _MPI_COVERAGE();
    fprintf(stderr, "Cannot allocate: %s\n", message);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}
/*=============================================================================================*/
void *_MPI_safeRealloc(void* oldBuffer, int size, char* message)
{
  _MPI_COVERAGE();
  oldBuffer = realloc (oldBuffer, size);
  _MPI_checkMemAlloc(oldBuffer, message);
  /* memset(oldBuffer+oldsize????,0,size); */
  return oldBuffer;

}
/*=============================================================================================*/
int _MPI_checkIntP (int *pt)
{
  _MPI_COVERAGE();
  if (pt == NULL)
    return _MPI_NOT_OK;
  return MPI_SUCCESS;
}
/*=============================================================================================*/
