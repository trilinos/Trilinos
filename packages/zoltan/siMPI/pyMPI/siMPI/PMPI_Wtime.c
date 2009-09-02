/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/**************************************************************************/
/* FILE   **************        PMPI_Wtime.c       ************************/
/**************************************************************************/
/* Author: Patrick Miller March  5 2004                                   */
/* Copyright (C) 2004 University of California Regents                    */
/**************************************************************************/
/* Simple "high resolution timer                                          */
/**************************************************************************/

#include "mpi.h"

#ifdef _MSC_VER
#include <gettimeofday.c>
#else
#include <sys/time.h>
#endif

/* STUB */
double PMPI_Wtime(void)
{
  struct timeval T;

  T.tv_sec = 0;
  T.tv_usec = 0;

  gettimeofday(&T,0);

  return T.tv_sec + 1e-6*T.tv_usec;
}
