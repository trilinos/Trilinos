/*
 * Copyright(C) 1999-2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwt - test write an ExodusII database file
 *
 *  This is a test program for the C binding of the EXODUS II
 *  database write routines.
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    int error = (funcall);                                                                         \
    printf("after %s, error = %d\n", TOSTRING(funcall), error);                                    \
    if (error != EX_NOERR) {                                                                       \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  int  i, j, k, m;
  int  node_list[100], elem_list[100], side_list[100];
  int  num_glo_vars, num_nod_vars, num_ele_vars, num_nset_vars, num_sset_vars;
  int *truth_tab, *nset_tab, *sset_tab;
  int  whole_time_step, num_time_steps;
  int  prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals, *nset_var_vals, *sset_var_vals;
  float  time_value;
  float  dist_fact[100];
  char  *var_names[7];

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  int CPU_word_size = 0; /* sizeof(float) */
  int IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */

  int exoid = ex_create("test.exo",     /* filename path */
                        EX_CLOBBER,     /* create mode */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */

  int num_dim       = 3;
  int num_nodes     = 33;
  int num_elem      = 7;
  int num_elem_blk  = 7;
  int num_node_sets = 2;
  int num_side_sets = 5;

  EXCHECK(ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets));

  /* write QA records; test empty and just blank-filled records */
  int num_qa_rec = 2;

  char *qa_record[2][4];
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

  int num_info = 3;

  char *info[3];
  info[0] = "This is the first information record.";
  info[1] = "";
  info[2] = "                                     ";

  EXCHECK(ex_put_info(exoid, num_info, info));

  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  float x[100], y[100], z[100];
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

  EXCHECK(ex_put_coord(exoid, x, y, z));

  char *coord_names[] = {"xcoor", "ycoor", "zcoor"};
  EXCHECK(ex_put_coord_names(exoid, coord_names));

  /* write element order map */

  int *elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i;
  }

  EXCHECK(ex_put_id_map(exoid, EX_ELEM_MAP, elem_map));

  free(elem_map);

  /* write element block parameters */

  int num_elem_in_block[]  = {1, 1, 1, 1, 1, 1, 1};
  int num_nodes_per_elem[] = {4,  /* elements in block #1 are 4-node quads  */
                              4,  /* elements in block #2 are 4-node quads  */
                              8,  /* elements in block #3 are 8-node hexes  */
                              4,  /* elements in block #4 are 4-node tetras */
                              6,  /* elements in block #5 are 6-node wedges */
                              8,  /* elements in block #6 are 8-node tetras */
                              3}; /* elements in block #7 are 3-node tris   */

  int ebids[] = {10, 11, 12, 13, 14, 15, 16};
  int nattr[] = {1, 1, 1, 1, 1, 1, 1};

  char *eb_type[] = {"quad", "quad", "hex", "tetra", "wedge", "tetra", "tri"};
  EXCHECK(ex_put_concat_elem_block(exoid, ebids, eb_type, num_elem_in_block, num_nodes_per_elem,
                                   nattr, 0));

  /* write element block properties */
  char *prop_names[] = {"MATL", "DENSITY"};
  EXCHECK(ex_put_prop_names(exoid, EX_ELEM_BLOCK, 2, prop_names));

  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], "MATL", 10));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], "MATL", 20));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], "MATL", 30));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[3], "MATL", 40));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[4], "MATL", 50));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[5], "MATL", 60));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[6], "MATL", 70));

  /* write element connectivity */
  int *connect = (int *)calloc(8, sizeof(int));
  connect[0]   = 1;
  connect[1]   = 2;
  connect[2]   = 3;
  connect[3]   = 4;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL));

  connect[0] = 5;
  connect[1] = 6;
  connect[2] = 7;
  connect[3] = 8;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[1], connect, NULL, NULL));

  connect[0] = 9;
  connect[1] = 10;
  connect[2] = 11;
  connect[3] = 12;
  connect[4] = 13;
  connect[5] = 14;
  connect[6] = 15;
  connect[7] = 16;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[2], connect, NULL, NULL));

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[3], connect, NULL, NULL));

  connect[0] = 21;
  connect[1] = 22;
  connect[2] = 23;
  connect[3] = 24;
  connect[4] = 25;
  connect[5] = 26;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[4], connect, NULL, NULL));

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;
  connect[4] = 27;
  connect[5] = 28;
  connect[6] = 30;
  connect[7] = 29;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[5], connect, NULL, NULL));

  connect[0] = 31;
  connect[1] = 32;
  connect[2] = 33;

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[6], connect, NULL, NULL));

  free(connect);

  /* write element block attributes */

  float attrib[] = {3.14159};
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib));

  attrib[0] = 6.14159;
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[1], attrib));

  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[2], attrib));

  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[3], attrib));

  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[4], attrib));

  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[5], attrib));

  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[6], attrib));

  /* write individual node sets */
  int nsids[]             = {20, 21};
  int num_nodes_per_set[] = {5, 3};
  {
    int                 num_df_per_set[] = {5, 3};
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = nsids;
    set_specs.num_entries_per_set = num_nodes_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = NULL;
    set_specs.sets_dist_index     = NULL;
    set_specs.sets_entry_list     = NULL;
    set_specs.sets_extra_list     = NULL;
    set_specs.sets_dist_fact      = NULL;

    EXCHECK(ex_put_concat_sets(exoid, EX_NODE_SET, &set_specs));
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

  EXCHECK(ex_put_set(exoid, EX_NODE_SET, 20, node_list, NULL));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, 20, dist_fact));

  node_list[0] = 20;
  node_list[1] = 21;
  node_list[2] = 22;

  dist_fact[0] = 1.1;
  dist_fact[1] = 2.1;
  dist_fact[2] = 3.1;

  EXCHECK(ex_put_set(exoid, EX_NODE_SET, 21, node_list, NULL));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, 21, dist_fact));
  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4));
  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5));

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  EXCHECK(ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array));

  /* Define the sideset params at one time, then write individually */
  int ssids[]            = {30, 31, 32, 33, 34};
  int num_elem_per_set[] = {2, 2, 7, 8, 10};
  {
    int                 num_df_per_set[] = {4, 4, 0, 0, 0};
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = ssids;
    set_specs.num_entries_per_set = num_elem_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = NULL;
    set_specs.sets_dist_index     = NULL;
    set_specs.sets_entry_list     = NULL;
    set_specs.sets_extra_list     = NULL;
    set_specs.sets_dist_fact      = NULL;
    EXCHECK(ex_put_concat_sets(exoid, EX_SIDE_SET, &set_specs));
  }

  /* write individual side sets */

  /* side set #1  - quad */

  elem_list[0] = 2;
  elem_list[1] = 2;

  side_list[0] = 4;
  side_list[1] = 2;

  dist_fact[0] = 30.0;
  dist_fact[1] = 30.1;
  dist_fact[2] = 30.2;
  dist_fact[3] = 30.3;

  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 30, elem_list, side_list));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 30, dist_fact));

  /* side set #2  - quad, spanning 2 elements  */

  elem_list[0] = 1;
  elem_list[1] = 2;

  side_list[0] = 2;
  side_list[1] = 3;

  dist_fact[0] = 31.0;
  dist_fact[1] = 31.1;
  dist_fact[2] = 31.2;
  dist_fact[3] = 31.3;

  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 31, elem_list, side_list));
  EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 31, dist_fact));

  /* side set #3  - hex */

  elem_list[0] = 3;
  elem_list[1] = 3;
  elem_list[2] = 3;
  elem_list[3] = 3;
  elem_list[4] = 3;
  elem_list[5] = 3;
  elem_list[6] = 3;

  side_list[0] = 5;
  side_list[1] = 3;
  side_list[2] = 3;
  side_list[3] = 2;
  side_list[4] = 4;
  side_list[5] = 1;
  side_list[6] = 6;

  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 32, elem_list, side_list));

  /* side set #4  - tetras */

  elem_list[0] = 4;
  elem_list[1] = 4;
  elem_list[2] = 4;
  elem_list[3] = 4;
  elem_list[4] = 6;
  elem_list[5] = 6;
  elem_list[6] = 6;
  elem_list[7] = 6;

  side_list[0] = 1;
  side_list[1] = 2;
  side_list[2] = 3;
  side_list[3] = 4;
  side_list[4] = 1;
  side_list[5] = 2;
  side_list[6] = 3;
  side_list[7] = 4;

  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 33, elem_list, side_list));

  /* side set #5  - wedges and tris */

  elem_list[0] = 5;
  elem_list[1] = 5;
  elem_list[2] = 5;
  elem_list[3] = 5;
  elem_list[4] = 5;
  elem_list[5] = 7;
  elem_list[6] = 7;
  elem_list[7] = 7;
  elem_list[8] = 7;
  elem_list[9] = 7;

  side_list[0] = 1;
  side_list[1] = 2;
  side_list[2] = 3;
  side_list[3] = 4;
  side_list[4] = 5;
  side_list[5] = 1;
  side_list[6] = 2;
  side_list[7] = 3;
  side_list[8] = 4;
  side_list[9] = 5;

  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 34, elem_list, side_list));

  EXCHECK(ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100));

  EXCHECK(ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101));

  /* write results variables parameters and names */
  num_glo_vars  = 1;
  num_nod_vars  = 2;
  num_ele_vars  = 3;
  num_nset_vars = 4;
  num_sset_vars = 7;

  truth_tab = (int *)calloc((num_elem_blk * num_ele_vars), sizeof(int));
  nset_tab  = (int *)calloc((num_node_sets * num_nset_vars), sizeof(int));
  sset_tab  = (int *)calloc((num_side_sets * num_sset_vars), sizeof(int));

  k = 0;
  for (i = 0; i < num_elem_blk; i++) {
    for (j = 0; j < num_ele_vars; j++) {
      truth_tab[k++] = 1;
    }
  }

  k = 0;
  for (i = 0; i < num_node_sets; i++) {
    for (j = 0; j < num_nset_vars; j++) {
      if (k % 2 == 0) {
        nset_tab[k++] = 1;
      }
      else {
        nset_tab[k++] = 0;
      }
    }
  }

  k = 0;
  for (i = 0; i < num_side_sets; i++) {
    for (j = 0; j < num_sset_vars; j++) {
      if (k % 2 == 0) {
        sset_tab[k++] = 0;
      }
      else {
        sset_tab[k++] = 1;
      }
    }
  }

  EXCHECK(ex_put_all_var_param(exoid, num_glo_vars, num_nod_vars, num_ele_vars, truth_tab,
                               num_nset_vars, nset_tab, num_sset_vars, sset_tab));

  free(truth_tab);
  free(nset_tab);
  free(sset_tab);

  var_names[0] = "glo_vars";
  EXCHECK(ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names));

  /*              12345678901234567890123456789012 */
  var_names[0] = "node_variable_a_very_long_name_0";
  var_names[1] = "nod_var1";
  EXCHECK(ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names));

  var_names[0] = "ele_var0";
  var_names[1] = "ele_var1";
  var_names[2] = "ele_var2";
  EXCHECK(ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names));

  var_names[0] = "nset_var0";
  var_names[1] = "nset_var1";
  var_names[2] = "nset_var2";
  var_names[3] = "nset_var3";
  EXCHECK(ex_put_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names));

  var_names[0] = "sset_var0";
  var_names[1] = "sset_var1";
  var_names[2] = "sset_var2";
  var_names[3] = "sset_var3";
  var_names[4] = "sset_var4";
  var_names[5] = "sset_var5";
  var_names[6] = "sset_var6";
  EXCHECK(ex_put_variable_names(exoid, EX_SIDE_SET, num_sset_vars, var_names));

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
  nset_var_vals  = (float *)calloc(5, CPU_word_size);
  sset_var_vals  = (float *)calloc(10, CPU_word_size);

  for (i = 0; i < num_time_steps; i++) {
    time_value = (float)(i + 1) / 100.;

    /* write time value */

    EXCHECK(ex_put_time(exoid, whole_time_step, &time_value));

    /* write global variables */

    for (j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (float)(j + 2) * time_value;
    }

    EXCHECK(ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 0, num_glo_vars, glob_var_vals));

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

    int kk = 0;
    for (j = 0; j < num_node_sets; j++) {
      for (k = 0; k < num_nset_vars; k++) {
        if (kk++ % 2 == 0) {
          for (m = 0; m < num_nodes_per_set[j]; m++) {
            nset_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          }
          EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODE_SET, k + 1, nsids[j],
                             num_nodes_per_set[j], nset_var_vals));
        }
      }
    }

    /* write sideset variables */

    kk = 0;
    for (j = 0; j < num_side_sets; j++) {
      for (k = 0; k < num_sset_vars; k++) {
        if (kk++ % 2 != 0) {
          for (m = 0; m < num_elem_per_set[j]; m++) {
            sset_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          }
          EXCHECK(ex_put_var(exoid, whole_time_step, EX_SIDE_SET, k + 1, ssids[j],
                             num_elem_per_set[j], sset_var_vals));
        }
      }
    }

    whole_time_step++;
    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    EXCHECK(ex_update(exoid));
  }
  free(glob_var_vals);
  free(nodal_var_vals);
  free(elem_var_vals);
  free(nset_var_vals);
  free(sset_var_vals);

  /* close the EXODUS files
   */
  EXCHECK(ex_close(exoid));
  return 0;
}
