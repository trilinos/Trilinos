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


#ifndef __HA_CONST_H
#define __HA_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Function prototypes */

extern int Zoltan_Get_Processor_Name(ZZ *, char *);
extern int Zoltan_Build_Machine_Desc(ZZ *);
extern int Zoltan_Divide_Machine(ZZ *, float *, int, MPI_Comm, int *, 
                                 int *, int *, int *, double *);

/* Misc. constants */
#define MAX_PROC_NAME_LEN (MPI_MAX_PROCESSOR_NAME+1)
#define MACHINE_DESC_FILE_DEFAULT "/etc/local/Zoltan_Machine_Desc"

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
