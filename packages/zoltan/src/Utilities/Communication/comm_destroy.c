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
#include "comm_const.h"
#include "mem_const.h"


int       LB_Comm_Destroy(
struct Comm_Obj **plan)		/* communication data structure pointer */
{

    /* Free fields of the communication object. */
    LB_FREE((void **) &((*plan)->status));
    LB_FREE((void **) &((*plan)->request));
    LB_FREE((void **) &((*plan)->procs_from));
    LB_FREE((void **) &((*plan)->lengths_from));
    LB_FREE((void **) &((*plan)->indices_from));
    LB_FREE((void **) &((*plan)->procs_to));
    LB_FREE((void **) &((*plan)->lengths_to));
    LB_FREE((void **) &((*plan)->indices_to));

    /* Free the communication object itself */
    LB_FREE((void **) plan);

    return(LB_COMM_OK);
}
