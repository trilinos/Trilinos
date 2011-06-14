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

#include "netcdf.h"

#include "rf_salsa.h"
#include "el_geom_const.h"
#include "rf_message.h"
#include "rf_io_const.h"
#include "rf_mp_const.h"

#include "exodusII.h"


/* need to hang on to this to write it out to the proc 0 file */
char    GeomTitle[MAX_LINE_LENGTH+1];

/************ R O U T I N E S   I N   T H I S   F I L E ***********************
*
*  Name_of_Routine		type		     Called by
*  ---------------------------------------------------------------
*
*  read_mesh_param ()			   	main:rf_salsa.c
*
******************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

extern void read_mesh_param (int *io_ws)

/*
*       Function which reads the sizing parameters from the EXODUS II database,
*       which contains the mesh.  This is used to cross-check against the
*       load balancer file.
*
*      ------------------------------------------------------------------------
*
*       Functions called:
*
*       check_exodus_error -- function which handles the error code returned by
*                             calls to EXODUS II API routines.
*
*      ------------------------------------------------------------------------
*/
{

/* Local variables */

  char   *yo = "read_mesh_param";
  char   *exofile;
  int     exoid, error, info[8];
  int     cpu_ws;
  float   version;

  extern void check_exodus_error (int error, char *);

/**************************** execution begins *******************************/
  cpu_ws = sizeof(float);

/*
 *         Generate the name of the parallel mesh file for this processor.
 *         Note: Default case is the scalar mesh file.
 */

  exofile = ExoFile;

  /* initialize the io word size on all of the processors */
  *io_ws = 0;

  /* Open the EXODUS II mesh file */

  if (Proc == 0) {

    exoid = ex_open (exofile, EX_READ, &cpu_ws, io_ws, &version);
    if (exoid == -1) {
      fprintf (stderr,
               "%s: ERROR openning up the mesh exoII file, %s\n", yo,
               exofile);
      exit (-1);
    }

    /* Read the initialization parameters */

    memset(GeomTitle, '\0', MAX_LINE_LENGTH*sizeof(char));
    error = ex_get_init (exoid, GeomTitle, &Num_Dim, &Num_Node, &Num_Elem,
                         &Num_Elem_Blk, &Num_Node_Set, &Num_Side_Set);
    check_exodus_error (error, "ex_get_init");

    if (Debug_Flag > 4) {
      printf ("\nKey definitions from regular exodus file (%s)\n", exofile);
      printf ("\tTitle of file: %s\n", GeomTitle);
      printf ("\tExodusII version number   = %f\n",   version);
      printf ("\tDimensionality of problem = %d\n",   Num_Dim);
      printf ("\tNumber of nodes           = %d\n",   Num_Node);
      printf ("\tNumber of elements        = %d\n",   Num_Elem);
      printf ("\tNumber of element blocks  = %d\n",   Num_Elem_Blk);
      printf ("\tNumber of node sets       = %d\n",   Num_Node_Set);
      printf ("\tNumber of side sets       = %d\n\n", Num_Side_Set);
    }

    info[0] = *io_ws;
    info[1] = Num_Dim;
    info[2] = Num_Node;
    info[3] = Num_Elem;
    info[4] = Num_Elem_Blk;
    info[5] = Num_Node_Set;
    info[6] = Num_Side_Set;
    /*
     * the processors need to know this, and this is the only place
     * that I can see to easily send it
     */
    info[7] = Restart_Info.Flag;

    /* Close the file */

    error = ex_close (exoid);  check_exodus_error (error, "ex_close");
  }

  brdcst (Proc, Num_Proc, (char *) info, (8 * sizeof(int)), 0);

  *io_ws            = info[0];
  Num_Dim           = info[1];
  Num_Node          = info[2];
  Num_Elem          = info[3];
  Num_Elem_Blk      = info[4];
  Num_Node_Set      = info[5];
  Num_Side_Set      = info[6];
  Restart_Info.Flag = info[7];

} /* END of routine read_mesh_params *****************************************/

