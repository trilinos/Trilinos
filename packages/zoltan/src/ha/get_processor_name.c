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
#include "ha_const.h"

int Zoltan_Get_Processor_Name(
   ZZ *zz,             /* The Zoltan structure.                */
   char *name          /* A string uniquely identifying the processor. 
                          We assume that at least MAX_PROC_NAME_LEN 
                          characters have been allocated.              */
)
{
/* This routine gets the name of the physical processor
   that an MPI process is running on.
 */

  int ierr = ZOLTAN_OK;
  int length;

  if (zz->Get_Processor_Name != NULL) {
    /* Use application-registered function */
    zz->Get_Processor_Name(zz->Get_Processor_Name_Data,
            name, &length, &ierr);
  }
  else {
    /* Use MPI_Get_processor_name by default */
    ierr = MPI_Get_processor_name(name, &length);
  }

  /* Add a trailing \0 to mark end of string */
  name[length] = '\0';

  return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
