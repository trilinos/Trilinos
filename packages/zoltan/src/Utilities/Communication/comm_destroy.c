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



#include <stdio.h>
#include "comm.h"
#include "zoltan_mem.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int       Zoltan_Comm_Destroy(
ZOLTAN_COMM_OBJ **plan)		/* communication data structure pointer */
{
    if (*plan == NULL) return ZOLTAN_OK;

    /* Free fields of the communication object. */
    ZOLTAN_FREE(&((*plan)->status));
    ZOLTAN_FREE(&((*plan)->request));
    ZOLTAN_FREE(&((*plan)->sizes));
    ZOLTAN_FREE(&((*plan)->sizes_to));
    ZOLTAN_FREE(&((*plan)->sizes_from));
    ZOLTAN_FREE(&((*plan)->starts_to_ptr));
    ZOLTAN_FREE(&((*plan)->starts_from_ptr));
    ZOLTAN_FREE(&((*plan)->indices_to_ptr));
    ZOLTAN_FREE(&((*plan)->indices_from_ptr));
    ZOLTAN_FREE(&((*plan)->indices_from));
    ZOLTAN_FREE(&((*plan)->indices_to));
    ZOLTAN_FREE(&((*plan)->lengths_from));
    ZOLTAN_FREE(&((*plan)->starts_to));
    ZOLTAN_FREE(&((*plan)->starts_from));
    ZOLTAN_FREE(&((*plan)->lengths_to));
    ZOLTAN_FREE(&((*plan)->procs_from));
    ZOLTAN_FREE(&((*plan)->procs_to));

    /* Free the communication object itself */
    ZOLTAN_FREE(plan);

    return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
