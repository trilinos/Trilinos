/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 EXPERIMENTAL -- DO NOT USE THIS CAPABILITY/CONCEPT FOR PRODUCTION WORK

 This is a test of an experimental capability to output only the transient
 data to an exodus file.  Eventually there will be a link between the
 transient file and the file containing the mesh data.  For now, we are
 only writing the transient data and verifying that the file has sufficient
 data to be understood by some of the seacas utilities.

 One of the ultimate goals is to make per-timestep output of results
 data faster since do not need the mesh definition bulk data
 (coordinates, connectivity, ...) for each step

 EXPERIMENTAL -- DO NOT USE THIS CAPABILITY/CONCEPT FOR PRODUCTION WORK
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exodusII.h"

int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_node_sets, num_side_sets, error;
  int  i, j, k, m;
  int  num_glo_vars, num_nod_vars, num_ele_vars;
  int *truth_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float  time_value;
  char  *var_names[3];
  char  *title = "This is a test";

  struct ex_block blocks[10];

  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */
  exoid = ex_create("test-results.exo", /* filename path */
                    EX_CLOBBER,         /* create mode */
                    &CPU_word_size,     /* CPU float word size in bytes */
                    &IO_word_size);     /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  num_dim       = 0;
  num_nodes     = 33;
  num_elem      = 7;
  num_elem_blk  = 7;
  num_node_sets = 0;
  num_side_sets = 0;

  ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
              num_side_sets);

  /* write element block parameters */
  for (i = 0; i < 10; i++) {
    blocks[i].type                = EX_ELEM_BLOCK;
    blocks[i].id                  = i + 10;
    blocks[i].num_entry           = 1;
    blocks[i].num_nodes_per_entry = 0;
    blocks[i].num_edges_per_entry = 0;
    blocks[i].num_faces_per_entry = 0;
    blocks[i].num_attribute       = 0;
  }

  ex_put_block_params(exoid, num_elem_blk, blocks);

  /* write results variables parameters and names */
  num_glo_vars = 1;

  var_names[0] = "glo_vars";

  error = ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  num_nod_vars = 2;
  /*              12345678901234567890123456789012 */
  var_names[0] = "node_variable_a_very_long_name_0";
  var_names[1] = "nod_var1";

  error = ex_put_variable_param(exoid, EX_NODAL, num_nod_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  num_ele_vars = 3;
  /*              0        1         2         3   */
  /*              12345678901234567890123456789012 */
  var_names[0] = "this_variable_name_is_short";
  var_names[1] = "this_variable_name_is_just_right";
  var_names[2] = "this_variable_name_is_tooooo_long";

  error = ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element variable truth table */
  truth_tab = (int *)calloc((num_elem_blk * num_ele_vars), sizeof(int));

  k = 0;
  for (i = 0; i < num_elem_blk; i++) {
    for (j = 0; j < num_ele_vars; j++) {
      truth_tab[k++] = 1;
    }
  }

  error = ex_put_truth_table(exoid, EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab);
  printf("after ex_put_elem_var_tab, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(truth_tab);

  /* for each time step, write the analysis results;
   * the code below fills the arrays glob_var_vals,
   * nodal_var_vals, and elem_var_vals with values for debugging purposes;
   * obviously the analysis code will populate these arrays
   */

  whole_time_step = 1;
  num_time_steps  = 10;

  glob_var_vals  = (float *)calloc(num_glo_vars, CPU_word_size);
  nodal_var_vals = (float *)calloc(num_nodes, CPU_word_size);
  elem_var_vals  = (float *)calloc(4, CPU_word_size);

  for (i = 0; i < num_time_steps; i++) {
    time_value = (float)(i + 1) / 100.0f;

    /* write time value */

    error = ex_put_time(exoid, whole_time_step, &time_value);
    printf("after ex_put_time, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write global variables */

    for (j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (float)(j + 2) * time_value;
    }

    error = ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals);
    printf("after ex_put_glob_vars, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write nodal variables */

    for (k = 1; k <= num_nod_vars; k++) {
      for (j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (float)k + ((float)(j + 1) * time_value);
      }

      error = ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals);
      printf("after ex_put_nodal_var, error = %d\n", error);
      if (error) {
        ex_close(exoid);
        exit(-1);
      }
    }

    /* write element variables */

    for (k = 1; k <= num_ele_vars; k++) {
      for (j = 0; j < num_elem_blk; j++) {
        for (m = 0; m < blocks[j].num_entry; m++) {
          elem_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          /* printf("elem_var_vals[%d]: %f\n",m,elem_var_vals[m]); */
        }
        error = ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, blocks[j].id,
                           blocks[j].num_entry, elem_var_vals);
        printf("after ex_put_elem_var, error = %d\n", error);
        if (error) {
          ex_close(exoid);
          exit(-1);
        }
      }
    }

    whole_time_step++;

    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    error = ex_update(exoid);
    printf("after ex_update, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
  }
  free(glob_var_vals);
  free(nodal_var_vals);
  free(elem_var_vals);

  /* close the EXODUS files
   */
  error = ex_close(exoid);
  printf("after ex_close, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }
  return 0;
}
