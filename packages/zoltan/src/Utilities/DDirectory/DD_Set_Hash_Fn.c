/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * This software is distributed under the GNU Lesser General Public License. *
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
#include <stdlib.h>

#include "DD.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */






/*************  Zoltan_DD_Set_Hash_Fn()  ***********************/


int Zoltan_DD_Set_Hash_Fn (
 Zoltan_DD_Directory *dd,              /* directory state information */
 unsigned int (*hash) (ZOLTAN_ID_PTR, int, unsigned int))
     {
     char *yo = "Zoltan_DD_Set_Hash_Fn" ;

     /* input sanity checking */
     if (dd == NULL || hash == NULL)
        {
        ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
        return ZOLTAN_DD_INPUT_ERROR ;
        }

     dd->hash = hash ;

     if (dd->debug_level > 0)
        ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Successful") ;

     return ZOLTAN_DD_NORMAL_RETURN ;
     }

