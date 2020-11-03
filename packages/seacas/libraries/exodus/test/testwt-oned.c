/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    if ((error = (funcall)) != NC_NOERR) {                                                         \
      fprintf(stderr, "ERROR Calling " #funcall ", error = %d\n", error);                          \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  int  error;
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_elem_in_block[10], num_nodes_per_elem[10];
  int  num_nodes_in_nset[10];
  int  num_node_sets, num_side_sets;
  int  i, j, k, m, *elem_map, *connect;
  int  node_list[100];
  int  ebids[10], nsids[10];
  int  num_qa_rec, num_info;
  int  num_nod_vars, num_ele_vars, num_nset_vars;
  int *truth_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;
  int  prop_array[2];

  int    elem_list[1], side_list[1];
  float *nodal_var_vals, *elem_var_vals;
  float *nset_var_vals;
  float  time_value;
  float  x[100];
  float  attrib[10], dist_fact[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char * block_names[10], *set_names[10];
  char * prop_names[2], *attrib_names[2];
  char * title = "This is a test";
  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS file */
  exoid = ex_create("test.exo",     /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for oned.e, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  EXCHECK(ex_set_max_name_length(exoid, 40));

  /* initialize file with parameters */

  num_dim       = 1;
  num_nodes     = 10;
  num_elem      = 10; /* 9 lines plus a point */
  num_elem_blk  = 3;
  num_node_sets = 2;
  num_side_sets = 2;

  EXCHECK(ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
                      num_side_sets));

  for (i = 0; i < num_nodes; i++) {
    x[i] = exp((float)i / 10.0f);
  }

  EXCHECK(ex_put_coord(exoid, x, NULL, NULL));

  coord_names[0] = "xcoor";
  EXCHECK(ex_put_coord_names(exoid, coord_names));

  /* Add nodal attributes */
  EXCHECK(ex_put_attr_param(exoid, EX_NODAL, 0, 1));

  EXCHECK(ex_put_one_attr(exoid, EX_NODAL, 0, 1, x));

  attrib_names[0] = "Node_attr_1";
  EXCHECK(ex_put_attr_names(exoid, EX_NODAL, 0, attrib_names));

  /* write element order map */
  elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = 10 * i;
  }

  EXCHECK(ex_put_map(exoid, elem_map));
  free(elem_map);

  /* write element block parameters */
  block_names[0] = "left_side";
  block_names[1] = "right_side";
  block_names[2] = "center";

  num_elem_in_block[0] = 4;
  num_elem_in_block[1] = 5;
  num_elem_in_block[2] = 1;

  num_nodes_per_elem[0] = 2;
  num_nodes_per_elem[1] = 2;
  num_nodes_per_elem[2] = 1;

  ebids[0] = 10;
  ebids[1] = 20;
  ebids[2] = 30;

  EXCHECK(ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "bar", num_elem_in_block[0],
                       num_nodes_per_elem[0], 0, 0, 1));
  EXCHECK(ex_put_block(exoid, EX_ELEM_BLOCK, ebids[1], "bar", num_elem_in_block[1],
                       num_nodes_per_elem[1], 0, 0, 1));
  EXCHECK(ex_put_block(exoid, EX_ELEM_BLOCK, ebids[2], "point", num_elem_in_block[2],
                       num_nodes_per_elem[2], 0, 0, 0));

  /* Write element block names */
  EXCHECK(ex_put_names(exoid, EX_ELEM_BLOCK, block_names));

  /* write element block properties */
  prop_names[0] = "DENSITY";
  EXCHECK(ex_put_prop_names(exoid, EX_ELEM_BLOCK, 1, prop_names));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], prop_names[0], 1));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], prop_names[0], 10));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], prop_names[0], 34));

  /* write element connectivity */
  connect = (int *)calloc(20, sizeof(int));
  for (i = 0; i < num_elem * 2; i += 2) {
    connect[i]     = i / 2 + 1;
    connect[i + 1] = i / 2 + 2;
  }

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL));
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[1], connect + 8, NULL, NULL));

  /* Circle */
  connect[0] = 5;
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[2], connect, NULL, NULL));
  free(connect);

  /* write element block attributes */
  for (i = 0; i < num_elem; i++) {
    attrib[i] = 3.14159 * i;
  }
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[1], attrib + num_elem_in_block[0]));

  attrib_names[0] = "THICKNESS";
  EXCHECK(ex_put_attr_names(exoid, EX_ELEM_BLOCK, ebids[0], attrib_names));
  attrib_names[0] = "WIDTH";
  EXCHECK(ex_put_attr_names(exoid, EX_ELEM_BLOCK, ebids[1], attrib_names));

  /* write individual node sets */
  num_nodes_in_nset[0] = 5;
  num_nodes_in_nset[1] = 3;

  nsids[0] = 20;
  nsids[1] = 21;

  EXCHECK(ex_put_set_param(exoid, EX_NODE_SET, nsids[0], 5, 5));

  node_list[0] = 1;
  node_list[1] = 3;
  node_list[2] = 5;
  node_list[3] = 7;
  node_list[4] = 9;

  dist_fact[0] = 1.0;
  dist_fact[1] = 2.0;
  dist_fact[2] = 3.0;
  dist_fact[3] = 4.0;
  dist_fact[4] = 5.0;

  EXCHECK(ex_put_set(exoid, EX_NODE_SET, nsids[0], node_list, NULL));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[0], dist_fact));

  EXCHECK(ex_put_set_param(exoid, EX_NODE_SET, nsids[1], 3, 3));

  node_list[0] = 2;
  node_list[1] = 4;
  node_list[2] = 6;

  dist_fact[0] = 1.0;
  dist_fact[1] = 2.0;
  dist_fact[2] = 3.0;

  EXCHECK(ex_put_set(exoid, EX_NODE_SET, nsids[1], node_list, NULL));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[1], dist_fact));

  /* Write node set names */
  set_names[0] = "all_odd_nodes";
  set_names[1] = "some_even_nodes";

  EXCHECK(ex_put_names(exoid, EX_NODE_SET, set_names));
  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, nsids[0], "FACE", 4));

  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, nsids[1], "FACE", 5));

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  EXCHECK(ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array));
  /* Add nodeset attributes */
  EXCHECK(ex_put_attr_param(exoid, EX_NODE_SET, nsids[0], 1));

  EXCHECK(ex_put_attr(exoid, EX_NODE_SET, nsids[0], x));

  attrib_names[0] = "Nodeset_attribute";
  EXCHECK(ex_put_attr_names(exoid, EX_NODE_SET, nsids[0], attrib_names));

  EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 1, 1, 1));
  EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 2, 1, 1));

  elem_list[0] = 1;
  side_list[0] = 1;
  dist_fact[0] = 2.0;
  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 1, elem_list, side_list));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 1, dist_fact));

  elem_list[0] = 9;
  side_list[0] = 2;
  dist_fact[0] = 3.0;
  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 2, elem_list, side_list));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 2, dist_fact));

  /* Write side set names */
  set_names[0] = "left_boundary";
  set_names[1] = "right_boundary";

  EXCHECK(ex_put_names(exoid, EX_SIDE_SET, set_names));
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

  EXCHECK(ex_put_qa(exoid, num_qa_rec, qa_record));

  /* write information records; test empty and just blank-filled records */
  num_info = 3;

  info[0] = "This is the first information record.";
  info[1] = "";
  info[2] = "                                     ";

  EXCHECK(ex_put_info(exoid, num_info, info));

  /* write results variables parameters and names */
  num_nod_vars = 2;
  /*              12345678901234567890123456789012 */
  var_names[0] = "node_variable_a_very_long_name_0";
  var_names[1] = "nod_var1";

  EXCHECK(ex_put_variable_param(exoid, EX_NODAL, num_nod_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names));

  num_ele_vars = 3;
  /*              0        1         2         3   */
  /*              12345678901234567890123456789012 */
  var_names[0] = "this_variable_name_is_short";
  var_names[1] = "this_variable_name_is_just_right";
  var_names[2] = "this_variable_name_is_tooooo_long";

  EXCHECK(ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names));

  num_nset_vars = 3;

  var_names[0] = "ns_var0";
  var_names[1] = "ns_var1";
  var_names[2] = "ns_var2";

  EXCHECK(ex_put_variable_param(exoid, EX_NODE_SET, num_nset_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names));

  /* write element variable truth table */
  truth_tab = (int *)calloc((num_elem_blk * num_ele_vars), sizeof(int));

  k = 0;
  for (i = 0; i < num_elem_blk; i++) {
    for (j = 0; j < num_ele_vars; j++) {
      truth_tab[k++] = 1;
    }
  }

  EXCHECK(ex_put_truth_table(exoid, EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab));

  free(truth_tab);

  /* for each time step, write the analysis results;
   * the code below fills the arrays glob_var_vals,
   * nodal_var_vals, and elem_var_vals with values for debugging purposes;
   */

  whole_time_step = 1;
  num_time_steps  = 10;

  nodal_var_vals = (float *)calloc(num_nodes, CPU_word_size);
  elem_var_vals  = (float *)calloc(num_elem, CPU_word_size);
  nset_var_vals  = (float *)calloc(10, CPU_word_size);

  for (i = 0; i < num_time_steps; i++) {
    time_value = (float)(i + 1) / 100.0f;

    /* write time value */
    EXCHECK(ex_put_time(exoid, whole_time_step, &time_value));

    /* write nodal variables */
    for (k = 1; k <= num_nod_vars; k++) {
      for (j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (float)k + ((float)(j + 1) * time_value);
      }
      EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals));
    }

    /* write element variables */
    for (k = 1; k <= num_ele_vars; k++) {
      for (j = 0; j < num_elem_blk; j++) {
        for (m = 0; m < num_elem_in_block[j]; m++) {
          elem_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, ebids[j], num_elem_in_block[j],
                           elem_var_vals));
      }
    }

    /* write nodeset variables */
    for (k = 1; k <= num_nset_vars; k++) {
      for (j = 0; j < num_node_sets; j++) {
        for (m = 0; m < num_nodes_in_nset[j]; m++) {
          nset_var_vals[m] = (float)(k + 3) + (float)(j + 4) + ((float)(m + 1) * time_value);
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODE_SET, k, nsids[j], num_nodes_in_nset[j],
                           nset_var_vals));
      }
    }

    whole_time_step++;

    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    EXCHECK(ex_update(exoid));
  }

  free(nodal_var_vals);
  free(elem_var_vals);
  free(nset_var_vals);

  /* close the EXODUS files */
  EXCHECK(ex_close(exoid));
  return 0;
}
