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

#include "DD_Const.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */






/*************  Zoltan_DD_Set_Hash_Fn()  ***********************/


int Zoltan_DD_Set_Hash_Fn (Zoltan_DD_Directory *dd,
 unsigned int (*hash) (LB_ID_PTR, int, unsigned int))
     {

     /* input sanity checking */
     if (dd == NULL || hash == NULL)
        return ZOLTAN_DD_INPUT_ERROR ;

     dd->hash = hash ;

     if (dd->debug_level > 0)
        printf ("ZOLTAN_DD_SET_HASH_FN(%d): Successful completion\n",
         dd->my_proc) ;

     return ZOLTAN_DD_NORMAL_RETURN ;
     }

