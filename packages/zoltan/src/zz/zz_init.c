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


#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for initializing Zoltan.
 *  These functions are all callable by the application. 
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Initialize(int argc, char **argv, float *ver)
{
/*
 *  Function to initialize values needed in load balancing tools.
 *  The function should be called after MPI_Init if the application
 *  uses MPI.
 */

int mpi_flag;

  /* 
   *  Test whether MPI is already initialized.  If not, call MPI_Init.
   */

  MPI_Initialized(&mpi_flag);

  if (!mpi_flag) {
    MPI_Init(&argc, &argv);
  }

  /*
   * Now return the version so that the user knows which version of
   * the libarary is being used without having to get the source
   * code.
   */
  *ver = ZOLTAN_VERSION_NUMBER;

  return (ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
