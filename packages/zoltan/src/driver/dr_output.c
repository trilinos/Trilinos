/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_output = "$Id$";
#endif

/*--------------------------------------------------------------------------*/
/* Purpose: Output the new element assignments.                             */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*    10 May 1999:       Date of creation.                                  */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"

int output_results(int Proc,
                   int Num_Proc,
                   PARIO_INFO_PTR pio_info,
                   ELEM_INFO elements[])
/*
 * For the first swipe at this, don't try to create a new
 * exodus/nemesis file or anything. Just get the global ids,
 * sort them, and print them to a new ascii file.
 */
{
  /* Local declarations. */
  char   par_out_fname[FILENAME_MAX+1], ctemp[FILENAME_MAX+1];

  int   *global_ids;
  int    i, j;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

  global_ids = (int *) malloc(Mesh.num_elems * sizeof(int));
  if (!global_ids) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  for (i = j = 0; i < Mesh.elem_array_len; i++) {
    if (elements[i].globalID >= 0) {
      global_ids[j] = elements[i].globalID;
      j++;
    }
  }

  sort_int(Mesh.num_elems, global_ids);

  /* generate the parallel filename for this processor */
  strcpy(ctemp, pio_info->pexo_fname);
  strcat(ctemp, ".out");
  gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc);

  fp = fopen(par_out_fname, "w");

  fprintf(fp, "Global element ids assigned to processor %d\n", Proc);
  for (i = 0; i < Mesh.num_elems; i++)
    fprintf(fp, "%d\n", global_ids[i]);

  fclose(fp);

  return 1;
}
