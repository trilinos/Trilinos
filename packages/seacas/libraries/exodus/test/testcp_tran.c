/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testcp - copy file test.exo created by testwt
 *
 * author - Sandia National Laboratories
 *          Larry A. Schoof - Original
 *
 *
 * environment - UNIX
 *
 * entry conditions -
 *   input parameters:
 *
 * exit conditions -
 *
 * revision history -
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int exoid, exoid1, error, idum;
  int CPU_word_size, IO_word_size;
  int num_nod_vars;
  int num_ele_vars;
  int i, j;

  float version;

  char *cdum = 0;

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use size in file */

  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  printf("\nafter ex_open\n");
  if (exoid < 0) {
    exit(1);
  }

  printf("test.exo is an EXODUSII file; version %4.2f\n", version);
  printf("         CPU word size %1d\n", CPU_word_size);
  printf("         I/O word size %1d\n", IO_word_size);
  ex_inquire(exoid, EX_INQ_API_VERS, &idum, &version, cdum);
  printf("EXODUSII API; version %4.2f\n", version);

  CPU_word_size = 4;
  IO_word_size  = 4;

  exoid1 = ex_create("testcp.exo",   /* filename */
                     EX_CLOBBER,     /* OK to overwrite */
                     &CPU_word_size, /* CPU float word size in bytes */
                     &IO_word_size); /* I/O float word size in bytes */

  printf("\nafter ex_create, exoid = %3d\n", exoid1);
  if (exoid1 < 0) {
    exit(1);
  }

  printf("         CPU word size %1d\n", CPU_word_size);
  printf("         I/O word size %1d\n", IO_word_size);

  error = ex_copy(exoid, exoid1);
  printf("\nafter ex_copy, error = %3d\n", error);

  /* See if any nodal transient variables on the input database */
  error = ex_get_variable_param(exoid, EX_NODE_BLOCK, &num_nod_vars);
  printf("\nafter ex_get_variable_param, error = %3d\n", error);
  error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_ele_vars);
  printf("\nafter ex_get_variable_param, error = %3d\n", error);

  if (num_nod_vars > 0) {
    num_nod_vars++;
    error = ex_put_variable_param(exoid1, EX_NODE_BLOCK, num_nod_vars);
    printf("\nafter ex_put_variable_param, error = %3d\n", error);
  }

  if (num_ele_vars > 0) {
    num_ele_vars++;
    error = ex_put_variable_param(exoid1, EX_ELEM_BLOCK, num_ele_vars);
    printf("\nafter ex_put_variable_param, error = %3d\n", error);
  }

  error = ex_copy_transient(exoid, exoid1);
  printf("\nafter ex_copy_transient, error = %3d\n", error);

  if (num_nod_vars > 0) {
    float node_vars[33]; /* Know that there are 33 nodes */
    error = ex_put_variable_name(exoid1, EX_NODE_BLOCK, num_nod_vars, "GregNode");
    printf("\nafter ex_put_variable_name, error = %3d\n", error);
    for (i = 0; i < 10; i++) {
      for (j = 0; j < 33; j++) {
        node_vars[j] = 3 + (float)(j + 1) * ((float)(i + 1) / 100.0);
      }
      error = ex_put_var(exoid1, i + 1, EX_NODE_BLOCK, num_nod_vars, 1, 33, node_vars);
      printf("\nafter ex_put_var, step %d, error = %3d\n", i + 1, error);
    }
  }

  if (num_ele_vars > 0) {
    error = ex_put_variable_name(exoid1, EX_ELEM_BLOCK, num_ele_vars, "GregElem");
    printf("\nafter ex_put_variable_name, error = %3d\n", error);

    /* Add the element variable to a single block.  1 element per block */
    float ele_vars[1];
    for (i = 0; i < 10; i++) { /* timesteps */
      for (j = 0; j < 1; j++) {
        ele_vars[j] = num_ele_vars + (float)(j + 1) * ((float)(i + 1) / 100.0);
      }
      error = ex_put_var(exoid1, i + 1, EX_ELEM_BLOCK, num_ele_vars, 10, 1, ele_vars);
      printf("\nafter ex_put_var, step %d, error = %3d\n", i + 1, error);
    }
  }
  error = ex_close(exoid);
  printf("\nafter ex_close, error = %3d\n", error);

  error = ex_close(exoid1);
  printf("\nafter ex_close, error = %3d\n", error);
  return 0;
}
