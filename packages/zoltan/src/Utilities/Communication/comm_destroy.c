/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
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


int       Zoltan_Comm_Destroy(
ZOLTAN_COMM_OBJ **plan)		/* communication data structure pointer */
{

    /* Free fields of the communication object. */
    ZOLTAN_FREE((void **) &((*plan)->status));
    ZOLTAN_FREE((void **) &((*plan)->request));
    ZOLTAN_FREE((void **) &((*plan)->sizes));
    ZOLTAN_FREE((void **) &((*plan)->sizes_to));
    ZOLTAN_FREE((void **) &((*plan)->sizes_from));
    ZOLTAN_FREE((void **) &((*plan)->starts_to_ptr));
    ZOLTAN_FREE((void **) &((*plan)->starts_from_ptr));
    ZOLTAN_FREE((void **) &((*plan)->indices_to_ptr));
    ZOLTAN_FREE((void **) &((*plan)->indices_from_ptr));
    ZOLTAN_FREE((void **) &((*plan)->indices_from));
    ZOLTAN_FREE((void **) &((*plan)->indices_to));
    ZOLTAN_FREE((void **) &((*plan)->lengths_from));
    ZOLTAN_FREE((void **) &((*plan)->starts_to));
    ZOLTAN_FREE((void **) &((*plan)->starts_from));
    ZOLTAN_FREE((void **) &((*plan)->lengths_to));
    ZOLTAN_FREE((void **) &((*plan)->procs_from));
    ZOLTAN_FREE((void **) &((*plan)->procs_to));

    /* Free the communication object itself */
    ZOLTAN_FREE((void **) plan);

    return(ZOLTAN_OK);
}
