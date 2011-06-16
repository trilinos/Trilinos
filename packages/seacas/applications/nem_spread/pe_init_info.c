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
#include <stdio.h>
#include <string.h>

#include "exodusII.h"

#include "rf_salsa.h"
#include "rf_message.h"

#include "rf_io_const.h"
#include "rf_mp_const.h"
#include "el_geom_const.h"

#include "pe_init_struct.h"


/*************** R O U T I N E S   I N   T H I S   F I L E ********************
*
*           NAME                     TYPE           CALLED BY
* ------------------------------------------------------------------------
*	exch_init_info ()              void           main
*	  load_init_info ()     static void           exch_init_info
*	  extract_init_info ()  static void           exch_init_info
*
******************************************************************************/
static void load_init_info    (struct Brdcst_Struct *),
	    extract_init_info (struct Brdcst_Struct *);

/*****************************************************************************/
/*****************************************************************************/


void exch_init_info (void)

/*
*      exch_init_info -
*
*  Routine that broadcast global variables from Proc = 0 to all processors.
*  This routine is only called for the multiprocessor case.
*
*/
{

/* Local Variables */
  struct Brdcst_Struct             init_info;

/* External variable declarations */

 /*
  *    Load the communication structure in order to do a single broadcast
  *    of a lot of miscellaneous information to all processors.
  *      The loading is only done by processor zero.
  */

  if (Proc == 0) {
    load_init_info (&init_info);
    if (Debug_Flag) printf ("Communicate init info struct\n\n");
  }

 /*
  *    Broadcast the communications structure and then extract it
  */

  brdcst(Proc, Num_Proc, (char *)&init_info, sizeof(struct Brdcst_Struct), 0);

  if (Proc == 0 && Debug_Flag) (void) printf("Extract init info data\n\n");

  extract_init_info (&init_info);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void load_init_info (struct Brdcst_Struct *init_info)

/*
*	John N. Shadid Div. 1421 SNL
*	date:           12/18/92
*
*	Routine to load initializing info in init_info structure
*
*	NOTE: Any variable added here must appear in (1) an included
*	file above (or an external statement)  (2) in the init_info
*	struct (see rf_init_info.h) and (3) in extract_init_info() below.
*/

{


  /* rf_io.h */
  init_info->Debug_Flag      = Debug_Flag;
  init_info->Num_Nod_Var     = Num_Nod_Var;
  init_info->Num_Elem_Var    = Num_Elem_Var;
  init_info->Num_Glob_Var    = Num_Glob_Var;
  init_info->Num_Nset_Var    = Num_Nset_Var;
  init_info->Num_Sset_Var    = Num_Sset_Var;

  /* el_geom.h */
  init_info->Num_Dim     = Num_Dim;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void extract_init_info (struct Brdcst_Struct *init_info)

/*
*
*	Routine to extract initializing info from init_info structure
*/

{

  /* rf_io.h */
  Debug_Flag      = init_info->Debug_Flag;
  Num_Nod_Var     = init_info->Num_Nod_Var;
  Num_Elem_Var    = init_info->Num_Elem_Var;
  Num_Glob_Var    = init_info->Num_Glob_Var;
  Num_Nset_Var    = init_info->Num_Nset_Var;
  Num_Sset_Var    = init_info->Num_Sset_Var;

  /* el_geom.h */
  Num_Dim             = init_info->Num_Dim;

}
/*****************************************************************************/
