#ifndef _cfei_petsc_h_
#define _cfei_petsc_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/*
   This header defines the prototype for the PETSc-specific function that
   creates the LinSysCore struct pointer, which is used by FEI_create.
*/

#include "fei_LinSysCore_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

int PETSc_LinSysCore_create(LinSysCore** lsc, MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif /*_cfei_petsc_h_*/

