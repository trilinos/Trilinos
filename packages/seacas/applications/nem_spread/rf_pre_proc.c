/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "rf_salsa.h"
#include "rf_fem_const.h"
#include "el_geom_const.h"
#include "rf_allo.h"
#include "el_elm.h"
#include "ps_pario_const.h"

/********* R O U T I N E S   D E F I N E D   I N   T H I S   F I L E **********
*
*
*      create_elem_types ()	 static void		pre_process ()
*
******************************************************************************/

void     create_elem_types(void);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void create_elem_types(void)

/*
 *      Function which creates a vector of element types for each element.
 *
 *       Author:          Scott Hutchinson (1421)
 */

{

  /* local variables */

  int i, j, iproc, ielem_type, ielem_count;

  /**************************** execution begins *****************************/

  Elem_Type = (int **)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2],
                                  sizeof(int *));

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {

    ielem_count = 0;

    Elem_Type[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1,
                                           Num_Internal_Elems[iproc] +
                                           Num_Border_Elems[iproc],
                                           sizeof(int));

  /*
   *     Loop through all the element blocks on the processor,
   *     setting the element type for each element in each block
   */

    for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++) {
      ielem_type = Proc_Elem_Blk_Types[iproc][i];
      for (j = 0; j < Proc_Num_Elem_In_Blk[iproc][i]; j++)
        Elem_Type[iproc][ielem_count++] = ielem_type;
    }

  }

} /* END of create_elem_types -----------------------------------------------*/

/*****************************************************************************/
/*                         END OF rf_pre_proc.c                              */
/*****************************************************************************/
