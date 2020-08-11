/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwt - test write an ExodusII database file
 *
 * author - Sandia National Laboratories
 *          Larry A. Schoof - Original
 *          Vic Yarberry    - Added headers and error logging
 *               7/7/93          Modified for use with Exodus 2.00
 *
 *
 * environment - UNIX
 *
 * entry conditions -
 *
 * exit conditions -
 *
 * revision history -
 *
 *  This is a test program for the C binding of the EXODUS II
 *  database write routines.
 *
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int num_elem_in_block[10];
  int num_node_sets, num_side_sets, error;
  int i, j, k, m;
  int node_list[100];
  int ebids[10];
  int num_qa_rec, num_info;
  int num_glo_vars, num_nod_vars, num_ele_vars;
  int whole_time_step, num_time_steps;
  int CPU_word_size, IO_word_size;
  int prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float  time_value;
  float  x[100], y[100], z[100];
  float  dist_fact[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */

  exoid = ex_create("test.exo",     /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */

  num_dim       = 3;
  num_nodes     = 33;
  num_elem      = 0;
  num_elem_blk  = 0;
  num_node_sets = 2;
  num_side_sets = 0;

  error = ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  x[0] = 0.0;
  y[0] = 0.0;
  z[0] = 0.0;
  x[1] = 1.0;
  y[1] = 0.0;
  z[1] = 0.0;
  x[2] = 1.0;
  y[2] = 1.0;
  z[2] = 0.0;
  x[3] = 0.0;
  y[3] = 1.0;
  z[3] = 0.0;

  /* Quad #2 */
  x[4] = 1.0;
  y[4] = 0.0;
  z[4] = 0.0;
  x[5] = 2.0;
  y[5] = 0.0;
  z[5] = 0.0;
  x[6] = 2.0;
  y[6] = 1.0;
  z[6] = 0.0;
  x[7] = 1.0;
  y[7] = 1.0;
  z[7] = 0.0;

  /* Hex #1 */
  x[8]  = 0.0;
  y[8]  = 0.0;
  z[8]  = 0.0;
  x[9]  = 10.0;
  y[9]  = 0.0;
  z[9]  = 0.0;
  x[10] = 10.0;
  y[10] = 0.0;
  z[10] = -10.0;
  x[11] = 1.0;
  y[11] = 0.0;
  z[11] = -10.0;
  x[12] = 1.0;
  y[12] = 10.0;
  z[12] = 0.0;
  x[13] = 10.0;
  y[13] = 10.0;
  z[13] = 0.0;
  x[14] = 10.0;
  y[14] = 10.0;
  z[14] = -10.0;
  x[15] = 1.0;
  y[15] = 10.0;
  z[15] = -10.0;

  /* Tetra #1 */
  x[16] = 0.0;
  y[16] = 0.0;
  z[16] = 0.0;
  x[17] = 1.0;
  y[17] = 0.0;
  z[17] = 5.0;
  x[18] = 10.0;
  y[18] = 0.0;
  z[18] = 2.0;
  x[19] = 7.0;
  y[19] = 5.0;
  z[19] = 3.0;

  /* Wedge #1 */
  x[20] = 3.0;
  y[20] = 0.0;
  z[20] = 6.0;
  x[21] = 6.0;
  y[21] = 0.0;
  z[21] = 0.0;
  x[22] = 0.0;
  y[22] = 0.0;
  z[22] = 0.0;
  x[23] = 3.0;
  y[23] = 2.0;
  z[23] = 6.0;
  x[24] = 6.0;
  y[24] = 2.0;
  z[24] = 2.0;
  x[25] = 0.0;
  y[25] = 2.0;
  z[25] = 0.0;

  /* Tetra #2 */
  x[26] = 2.7;
  y[26] = 1.7;
  z[26] = 2.7;
  x[27] = 6.0;
  y[27] = 1.7;
  z[27] = 3.3;
  x[28] = 5.7;
  y[28] = 1.7;
  z[28] = 1.7;
  x[29] = 3.7;
  y[29] = 0.0;
  z[29] = 2.3;

  /* 3d Tri */
  x[30] = 0.0;
  y[30] = 0.0;
  z[30] = 0.0;
  x[31] = 10.0;
  y[31] = 0.0;
  z[31] = 0.0;
  x[32] = 10.0;
  y[32] = 10.0;
  z[32] = 10.0;

  error = ex_put_coord(exoid, x, y, z);
  printf("after ex_put_coord, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  coord_names[0] = "xcoor";
  coord_names[1] = "ycoor";
  coord_names[2] = "zcoor";

  error = ex_put_coord_names(exoid, coord_names);
  printf("after ex_put_coord_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write individual node sets */

  error = ex_put_set_param(exoid, EX_NODE_SET, 20, 5, 5);
  printf("after ex_put_node_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  node_list[0] = 10;
  node_list[1] = 11;
  node_list[2] = 12;
  node_list[3] = 13;
  node_list[4] = 14;

  dist_fact[0] = 1.0;
  dist_fact[1] = 2.0;
  dist_fact[2] = 3.0;
  dist_fact[3] = 4.0;
  dist_fact[4] = 5.0;

  error = ex_put_set(exoid, EX_NODE_SET, 20, node_list, NULL);
  printf("after ex_put_node_set, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, 20, dist_fact);
  printf("after ex_put_node_set_dist_fact, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_param(exoid, EX_NODE_SET, 21, 3, 3);
  printf("after ex_put_node_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  node_list[0] = 20;
  node_list[1] = 21;
  node_list[2] = 22;

  dist_fact[0] = 1.1;
  dist_fact[1] = 2.1;
  dist_fact[2] = 3.1;

  error = ex_put_set(exoid, EX_NODE_SET, 21, node_list, NULL);
  printf("after ex_put_node_set, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, 21, dist_fact);
  printf("after ex_put_node_set_dist_fact, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5);
  printf("after ex_put_prop, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  error = ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array);
  printf("after ex_put_prop_array, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write QA records; test empty and just blank-filled records */

  num_qa_rec = 2;

  qa_record[0][0] = "TESTWT";
  qa_record[0][1] = "testwt";
  qa_record[0][2] = "07/07/93";
  qa_record[0][3] = "15:41:33";
  qa_record[1][0] = "";
  qa_record[1][1] = "                            ";
  qa_record[1][2] = "";
  qa_record[1][3] = "                        ";

  error = ex_put_qa(exoid, num_qa_rec, qa_record);
  printf("after ex_put_qa, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write information records; test empty and just blank-filled records */

  num_info = 3;

  info[0] = "This is the first information record.";
  info[1] = "";
  info[2] = "                                     ";

  error = ex_put_info(exoid, num_info, info);
  printf("after ex_put_info, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  num_ele_vars = 0;

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
    time_value = (float)(i + 1) / 100.;

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
        for (m = 0; m < num_elem_in_block[j]; m++) {
          elem_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          /* printf("elem_var_vals[%d]: %f\n",m,elem_var_vals[m]); */
        }
        error = ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, ebids[j], num_elem_in_block[j],
                           elem_var_vals);
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
