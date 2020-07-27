/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwt - test write an ExodusII database file (testwt-nfaced.exo)
 *
 *****************************************************************************/

#include "exodusII.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_elem_in_block[10], num_total_nodes_per_blk[10];
  int  num_face_in_block[10], num_total_faces_per_blk[10];
  int  num_node_sets, error;
  int  i, j, k, m, *connect;
  int  bids[10], nnpe[10];
  int  num_qa_rec, num_info;
  int  num_glo_vars, num_nod_vars, num_ele_vars;
  int *truth_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float  time_value;
  float  x[100], y[100], z[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char * block_names[10], *nset_names[10], *sset_names[10];
  char * prop_names[2], *attrib_names[2];
  char * title = "This is a test";
  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */

  exoid = ex_create("test-nfaced.exo", /* filename path */
                    EX_CLOBBER,        /* create mode */
                    &CPU_word_size,    /* CPU float word size in bytes */
                    &IO_word_size);    /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */
  {
    ex_init_params par;

    num_dim       = 3;
    num_nodes     = 14;
    num_elem      = 1;
    num_elem_blk  = 1;
    num_node_sets = 0;

    ex_copy_string(par.title, title, MAX_LINE_LENGTH + 1);
    par.num_dim       = num_dim;
    par.num_nodes     = num_nodes;
    par.num_edge      = 0;
    par.num_edge_blk  = 0;
    par.num_face      = 5;
    par.num_face_blk  = 1;
    par.num_elem      = num_elem;
    par.num_elem_blk  = num_elem_blk;
    par.num_node_sets = num_node_sets;
    par.num_edge_sets = 0;
    par.num_face_sets = 0;
    par.num_side_sets = 0;
    par.num_elem_sets = 0;
    par.num_node_maps = 0;
    par.num_edge_maps = 0;
    par.num_face_maps = 0;
    par.num_elem_maps = 0;

    error = ex_put_init_ext(exoid, &par);

    printf("after ex_put_init_ext, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }
  }

  /* write nodal coordinates values and names to database */
  x[0]  = 0.00000e+00;
  y[0]  = 0.00000e+00;
  z[0]  = 0.00000e+00;
  x[1]  = 2.00000e+00;
  y[1]  = 0.00000e+00;
  z[1]  = 0.00000e+00;
  x[2]  = 0.00000e+00;
  y[2]  = 2.00000e+00;
  z[2]  = 0.00000e+00;
  x[3]  = 2.00000e+00;
  y[3]  = 2.00000e+00;
  z[3]  = 0.00000e+00;
  x[4]  = 0.00000e+00;
  y[4]  = 0.00000e+00;
  z[4]  = 2.00000e+00;
  x[5]  = 2.00000e+00;
  y[5]  = 0.00000e+00;
  z[5]  = 2.00000e+00;
  x[6]  = 0.00000e+00;
  y[6]  = 2.00000e+00;
  z[6]  = 2.00000e+00;
  x[7]  = 2.00000e+00;
  y[7]  = 2.00000e+00;
  z[7]  = 2.00000e+00;
  x[8]  = 0.00000e+00;
  y[8]  = 3.50000e+00;
  z[8]  = 1.00000e+00;
  x[9]  = 2.00000e+00;
  y[9]  = 3.50000e+00;
  z[9]  = 1.00000e+00;
  x[10] = 0.00000e+00;
  y[10] = 3.00000e+00;
  z[10] = 1.50000e+00;
  x[11] = 2.00000e+00;
  y[11] = 3.00000e+00;
  z[11] = 1.50000e+00;
  x[12] = 0.00000e+00;
  y[12] = 3.00000e+00;
  z[12] = 0.50000e+00;
  x[13] = 2.00000e+00;
  y[13] = 3.00000e+00;
  z[13] = 0.50000e+00;

  error = ex_put_coord(exoid, x, y, z);
  printf("after ex_put_coord, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  coord_names[0] = "x";
  coord_names[1] = "y";
  coord_names[2] = "z";

  error = ex_put_coord_names(exoid, coord_names);
  printf("after ex_put_coord_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* Write the face block parameters */
  block_names[0]             = "face_block_1";
  num_face_in_block[0]       = 15;
  num_total_nodes_per_blk[0] = 54;
  bids[0]                    = 10;

  error = ex_put_block(exoid, EX_FACE_BLOCK, bids[0], "nsided", num_face_in_block[0],
                       num_total_nodes_per_blk[0], 0, 0, 0);
  printf("after ex_put_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write face connectivity */

  connect = (int *)calloc(num_total_nodes_per_blk[0], sizeof(int));

  i = 0;
  j = 0;

  connect[i++] = 5;
  connect[i++] = 6;
  connect[i++] = 8; /* connectivity of face 1 of element 1 */

  connect[i++] = 2;
  connect[i++] = 1;
  connect[i++] = 4; /* face 2 of element 1 */

  connect[i++] = 6;
  connect[i++] = 2;
  connect[i++] = 4;
  connect[i++] = 8; /* face 3 of element 1 */

  connect[i++] = 8;
  connect[i++] = 4;
  connect[i++] = 1;
  connect[i++] = 5; /* face 4 of element 1 */

  connect[i++] = 1;
  connect[i++] = 2;
  connect[i++] = 6;
  connect[i++] = 5; /*  face 5 of element 1 */

  connect[i++] = 5;
  connect[i++] = 8;
  connect[i++] = 7; /* connectivity of face 1 of element 2 */

  connect[i++] = 1;
  connect[i++] = 2;
  connect[i++] = 3;
  connect[i++] = 4;
  nnpe[j++]    = 4;

  connect[i++] = 5;
  connect[i++] = 3;
  connect[i++] = 4;
  connect[i++] = 6;
  nnpe[j++]    = 4;

  connect[i++] = 5;
  connect[i++] = 1;
  connect[i++] = 2;
  connect[i++] = 6;
  nnpe[j++]    = 4;

  connect[i++] = 6;
  connect[i++] = 2;
  connect[i++] = 4;
  nnpe[j++]    = 3;

  connect[i++] = 5;
  connect[i++] = 3;
  connect[i++] = 1;
  nnpe[j++]    = 3;

  assert(i == num_total_nodes_per_blk[0]);
  assert(j == num_face_in_block[0]);

  error = ex_put_conn(exoid, EX_FACE_BLOCK, bids[0], connect, NULL, NULL);
  printf("after ex_put_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(connect);
  connect = NULL;

  error = ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, bids[0], nnpe);
  printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element block parameters */
  block_names[0] = "nfaced_1";

  num_elem_in_block[0]       = 1;
  num_total_faces_per_blk[0] = 5;

  bids[0] = 10;

  error = ex_put_block(exoid, EX_ELEM_BLOCK, bids[0], "nfaced", num_elem_in_block[0], 0, 0,
                       num_total_faces_per_blk[0], 0);
  printf("after ex_put_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* Write face block names */
  error = ex_put_names(exoid, EX_FACE_BLOCK, block_names);
  printf("after ex_put_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* Write element block names */
  error = ex_put_names(exoid, EX_ELEM_BLOCK, block_names);
  printf("after ex_put_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element-face connectivity */
  connect = (int *)calloc(num_total_faces_per_blk[0], sizeof(int));

  i            = 0;
  j            = 0;
  connect[i++] = 1;
  connect[i++] = 2;
  connect[i++] = 3;
  connect[i++] = 4;
  connect[i++] = 5;
  nnpe[j++]    = 5; /* Number of faces per element */

  assert(i == num_total_faces_per_blk[0]);
  assert(j == num_elem_in_block[0]);

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, bids[0], NULL, NULL, connect);
  printf("after ex_put_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(connect);

  error = ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, bids[0], nnpe);
  printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write QA records; test empty and just blank-filled records */
  num_qa_rec = 2;

  qa_record[0][0] = "TESTWT-NFACED";
  qa_record[0][1] = "testwt-nfaced";
  qa_record[0][2] = "2010/02/15";
  qa_record[0][3] = "06:35:15";
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
  var_names[1] = EX_NODE_SET;

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

  var_names[0] = "ele_var0";
  var_names[1] = "ele_var1";
  var_names[2] = "ele_var2";

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
  elem_var_vals  = (float *)calloc(8, CPU_word_size);

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

      error = ex_put_nodal_var(exoid, whole_time_step, k, num_nodes, nodal_var_vals);
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
        error = ex_put_elem_var(exoid, whole_time_step, k, bids[j], num_elem_in_block[j],
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
