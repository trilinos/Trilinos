/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exodusII.h"
#include "exodusII_int.h"

#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    int f_error = (funcall);                                                                       \
    printf("after %s, error = %d\n", TOSTRING(funcall), f_error);                                  \
    if (f_error != EX_NOERR && f_error != EX_WARN) {                                               \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */
  int CPU_word_size = 8; /* sizeof(double) */
  int IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */
  int exoid = ex_create("test.exo",     /* filename path */
                        EX_CLOBBER,     /* create mode */
                        &CPU_word_size, /* CPU double word size in bytes */
                        &IO_word_size); /* I/O double word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */
  int num_dim       = 3;
  int num_nodes     = 33;
  int num_elem      = 7;
  int num_elem_blk  = 7;
  int num_node_sets = 2;
  int num_side_sets = 5;

  char *title = "This is a test";
  EXCHECK(ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
                      num_side_sets));

  /* clang-format off */
  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  double x[100], y[100], z[100];
  x[0] = 0.0;  y[0] = 0.0;  z[0] = 0.0;
  x[1] = 1.0;  y[1] = 0.0;  z[1] = 0.0;
  x[2] = 1.0;  y[2] = 1.0;  z[2] = 0.0;
  x[3] = 0.0;  y[3] = 1.0;  z[3] = 0.0;

  /* Quad #2 */
  x[4] = 1.0;  y[4] = 0.0;  z[4] = 0.0;
  x[5] = 2.0;  y[5] = 0.0;  z[5] = 0.0;
  x[6] = 2.0;  y[6] = 1.0;  z[6] = 0.0;
  x[7] = 1.0;  y[7] = 1.0;  z[7] = 0.0;

  /* Hex #1 */
  x[8]  = 0.0;   y[8]  = 0.0;  z[8]  = 0.0;
  x[9]  = 10.0;  y[9]  = 0.0;  z[9]  = 0.0;
  x[10] = 10.0;  y[10] = 0.0;  z[10] = -10.0;
  x[11] = 1.0;   y[11] = 0.0;  z[11] = -10.0;
  x[12] = 1.0;   y[12] = 10.0; z[12] = 0.0;
  x[13] = 10.0;  y[13] = 10.0; z[13] = 0.0;
  x[14] = 10.0;  y[14] = 10.0; z[14] = -10.0;
  x[15] = 1.0;   y[15] = 10.0; z[15] = -10.0;

  /* Tetra #1 */
  x[16] = 0.0;  y[16] = 0.0;  z[16] = 0.0;
  x[17] = 1.0;  y[17] = 0.0;  z[17] = 5.0;
  x[18] = 10.0; y[18] = 0.0;  z[18] = 2.0;
  x[19] = 7.0;  y[19] = 5.0;  z[19] = 3.0;

  /* Wedge #1 */
  x[20] = 3.0;  y[20] = 0.0;  z[20] = 6.0;
  x[21] = 6.0;  y[21] = 0.0;  z[21] = 0.0;
  x[22] = 0.0;  y[22] = 0.0;  z[22] = 0.0;
  x[23] = 3.0;  y[23] = 2.0;  z[23] = 6.0;
  x[24] = 6.0;  y[24] = 2.0;  z[24] = 2.0;
  x[25] = 0.0;  y[25] = 2.0;  z[25] = 0.0;

  /* Tetra #2 */
  x[26] = 2.7;  y[26] = 1.7;  z[26] = 2.7;
  x[27] = 6.0;  y[27] = 1.7;  z[27] = 3.3;
  x[28] = 5.7;  y[28] = 1.7;  z[28] = 1.7;
  x[29] = 3.7;  y[29] = 0.0;  z[29] = 2.3;

  /* 3d Tri */
  x[30] = 0.0;  y[30] = 0.0;  z[30] = 0.0;
  x[31] = 10.0; y[31] = 0.0;  z[31] = 0.0;
  x[32] = 10.0; y[32] = 10.0; z[32] = 10.0;
  /* clang-format on */

  EXCHECK(ex_put_coord(exoid, x, y, z));

  char *coord_names[] = {"xcoor", "ycoor", "zcoor"};
  EXCHECK(ex_put_coord_names(exoid, coord_names));

  /* Add nodal attributes */
  EXCHECK(ex_put_attr_param(exoid, EX_NODAL, 0, 2));
  EXCHECK(ex_put_one_attr(exoid, EX_NODAL, 0, 1, x));
  EXCHECK(ex_put_one_attr(exoid, EX_NODAL, 0, 2, y));

  {
    char *attrib_names[] = {"Node_attr_1", "Node_attr_2"};
    EXCHECK(ex_put_attr_names(exoid, EX_NODAL, 0, attrib_names));
  }

  /* write element id map */
  int *elem_map = (int *)calloc(num_elem, sizeof(int));
  for (int i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i * 10;
  }

  EXCHECK(ex_put_id_map(exoid, EX_ELEM_MAP, elem_map));

  free(elem_map);

  /* write element block parameters */
  struct ex_block blocks[10];
  for (int i = 0; i < 10; i++) {
    blocks[i].type                = EX_ELEM_BLOCK;
    blocks[i].id                  = 0;
    blocks[i].num_entry           = 0;
    blocks[i].num_nodes_per_entry = 0;
    blocks[i].num_edges_per_entry = 0;
    blocks[i].num_faces_per_entry = 0;
    blocks[i].num_attribute       = 0;
  }

  char *block_names[10];
  block_names[0] = "block_1";
  block_names[1] = "block_2";
  block_names[2] = "block_3";
  block_names[3] = "block_4";
  block_names[4] = "block_5";
  block_names[5] = "block_6";
  block_names[6] = "block_7";

  ex_copy_string(blocks[0].topology, "quad", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[1].topology, "quad", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[2].topology, "hex", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[3].topology, "tetra", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[4].topology, "wedge", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[5].topology, "tetra", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[6].topology, "tri", MAX_STR_LENGTH + 1);

  blocks[0].num_entry = 1;
  blocks[1].num_entry = 1;
  blocks[2].num_entry = 1;
  blocks[3].num_entry = 1;
  blocks[4].num_entry = 1;
  blocks[5].num_entry = 1;
  blocks[6].num_entry = 1;

  blocks[0].num_attribute = 1;
  blocks[1].num_attribute = 1;
  blocks[2].num_attribute = 1;
  blocks[3].num_attribute = 1;
  blocks[4].num_attribute = 1;
  blocks[5].num_attribute = 1;
  blocks[6].num_attribute = 1;

  blocks[0].num_nodes_per_entry = 4; /* elements in block #1 are 4-node quads  */
  blocks[1].num_nodes_per_entry = 4; /* elements in block #2 are 4-node quads  */
  blocks[2].num_nodes_per_entry = 8; /* elements in block #3 are 8-node hexes  */
  blocks[3].num_nodes_per_entry = 4; /* elements in block #4 are 4-node tetras */
  blocks[4].num_nodes_per_entry = 6; /* elements in block #5 are 6-node wedges */
  blocks[5].num_nodes_per_entry = 8; /* elements in block #6 are 8-node tetras */
  blocks[6].num_nodes_per_entry = 3; /* elements in block #7 are 3-node tris   */

  blocks[0].id = 10;
  blocks[1].id = 11;
  blocks[2].id = 12;
  blocks[3].id = 13;
  blocks[4].id = 14;
  blocks[5].id = 15;
  blocks[6].id = 16;

  /* Generate an error that name is not found since blocks have not
     yet been defined
  */
  int error = ex_put_name(exoid, EX_ELEM_BLOCK, blocks[0].id, block_names[0]);
  printf("after ex_put_name, error = %d\n", error);

  EXCHECK(ex_put_block_params(exoid, num_elem_blk, blocks));

  /* Write element block names */
  for (int i = 0; i < num_elem_blk; i++) {
    EXCHECK(ex_put_name(exoid, EX_ELEM_BLOCK, blocks[i].id, block_names[i]));
  }

  /* write element block properties */
  /*               12345678901234567890123456789012 */
  char *prop_names[2];
  prop_names[0] = "MATERIAL_PROPERTY_LONG_NAME_32CH";
  prop_names[1] = "DENSITY";

  EXCHECK(ex_put_prop_names(exoid, EX_ELEM_BLOCK, 2, prop_names));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[0].id, prop_names[0], 10));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[1].id, prop_names[0], 20));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[2].id, prop_names[0], 30));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[3].id, prop_names[0], 40));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[4].id, prop_names[0], 50));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[5].id, prop_names[0], 60));
  EXCHECK(ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[6].id, prop_names[0], 70));

  /* write element connectivity */
  {
    int connect[] = {1, 2, 3, 4};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[0].id, connect, NULL, NULL));
  }

  {
    int connect[] = {5, 6, 7, 8};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[1].id, connect, NULL, NULL));
  }

  {
    int connect[] = {9, 10, 11, 12, 13, 14, 15, 16};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[2].id, connect, NULL, NULL));
  }

  {
    int connect[] = {17, 18, 19, 20};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[3].id, connect, NULL, NULL));
  }

  {
    int connect[] = {21, 22, 23, 24, 25, 26};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[4].id, connect, NULL, NULL));
  }

  {
    int connect[] = {17, 18, 19, 20, 27, 28, 30, 29};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[5].id, connect, NULL, NULL));
  }

  {
    int connect[] = {31, 32, 33};
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[6].id, connect, NULL, NULL));
  }

  /* write element block attributes */
  double attrib[1];
  attrib[0] = 3.14159;
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[0].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[0].id, attrib));

  attrib[0] = 6.14159;
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[1].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[2].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[3].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[4].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[5].id, attrib));
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[6].id, attrib));

  {
    char *attrib_names[] = {"THICKNESS"};
    for (int i = 0; i < num_elem_blk; i++) {
      EXCHECK(ex_put_attr_names(exoid, EX_ELEM_BLOCK, blocks[i].id, attrib_names));
    }
  }

  /* write individual node sets */
  int num_nodes_in_nset[] = {5, 3};
  int nsids[]             = {20, 21};

  {
    EXCHECK(ex_put_set_param(exoid, EX_NODE_SET, nsids[0], 5, 5));

    int    node_list[] = {10, 11, 12, 13, 14};
    double dist_fact[] = {1.0, 2.0, 3.0, 4.0, 5.0};

    EXCHECK(ex_put_set(exoid, EX_NODE_SET, nsids[0], node_list, NULL));
    EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[0], dist_fact));
  }

  {
    EXCHECK(ex_put_set_param(exoid, EX_NODE_SET, nsids[1], 3, 3));

    int    node_list[] = {20, 21, 22};
    double dist_fact[] = {1.1, 2.1, 3.1};

    EXCHECK(ex_put_set(exoid, EX_NODE_SET, nsids[1], node_list, NULL));
    EXCHECK(ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[1], dist_fact));
  }

  /* Write node set names */
  char *nset_names[] = {"nset_1", "nset_2"};
  EXCHECK(ex_put_names(exoid, EX_NODE_SET, nset_names));
  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, nsids[0], "FACE", 4));
  EXCHECK(ex_put_prop(exoid, EX_NODE_SET, nsids[1], "FACE", 5));

  int prop_array[] = {1000, 2000};
  EXCHECK(ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array));

  /* Add nodeset attributes */
  EXCHECK(ex_put_attr_param(exoid, EX_NODE_SET, nsids[0], 1));
  EXCHECK(ex_put_attr(exoid, EX_NODE_SET, nsids[0], x));

  {
    char *attrib_names[] = {"Nodeset_attribute"};
    EXCHECK(ex_put_attr_names(exoid, EX_NODE_SET, nsids[0], attrib_names));
  }

  /* write individual side sets */
  int num_face_in_sset[] = {2, 2, 7, 8, 10};
  int ssids[]            = {30, 31, 32, 33, 34};

  {
    /* side set #1  - quad */
    EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, ssids[0], 2, 4));

    int    elem_list[] = {2, 2};
    int    side_list[] = {4, 2};
    double dist_fact[] = {30.0, 30.1, 30.2, 30.3};

    EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 30, elem_list, side_list));
    EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 30, dist_fact));
  }

  {
    /* side set #2  - quad, spanning 2 elements  */
    EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 31, 2, 4));

    int    elem_list[] = {1, 2};
    int    side_list[] = {2, 3};
    double dist_fact[] = {31.0, 31.1, 31.2, 31.3};

    EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 31, elem_list, side_list));
    EXCHECK(ex_put_set_dist_fact(exoid, EX_SIDE_SET, 31, dist_fact));
  }

  {
    /* side set #3  - hex */
    EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 32, 7, 0));

    int elem_list[] = {3, 3, 3, 3, 3, 3, 3};
    int side_list[] = {5, 3, 3, 2, 4, 1, 6};

    EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 32, elem_list, side_list));
  }

  {
    /* side set #4  - tetras */
    EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 33, 8, 0));

    int elem_list[] = {4, 4, 4, 4, 6, 6, 6, 6};
    int side_list[] = {1, 2, 3, 4, 1, 2, 3, 4};

    EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 33, elem_list, side_list));
  }

  {
    /* side set #5  - wedges and tris */
    EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 34, 10, 0));

    int elem_list[] = {5, 5, 5, 5, 5, 7, 7, 7, 7, 7};
    int side_list[] = {1, 2, 3, 4, 5, 1, 2, 3, 4, 5};

    EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 34, elem_list, side_list));
  }

  /* Write side set names */
  char *sset_names[] = {"sset_1", "sset_2", "sset_3", "sset_4", "sset_5"};
  EXCHECK(ex_put_names(exoid, EX_SIDE_SET, sset_names));
  EXCHECK(ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100));
  EXCHECK(ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101));

  /* write QA records; test empty and just blank-filled records */
  int   num_qa_rec = 2;
  char *qa_record[2][4];
  qa_record[0][0] = "TESTWT";
  qa_record[0][1] = "testwt";
  qa_record[0][2] = "07/07/93";
  qa_record[0][3] = "15:41:33";
  qa_record[1][0] = "Thirty-Two character QA Record|";
  qa_record[1][1] = "                            ";
  qa_record[1][2] = "";
  qa_record[1][3] = "                        ";

  EXCHECK(ex_put_qa(exoid, num_qa_rec, qa_record));

  /* write information records; test empty and just blank-filled records */
  char *info[3];
  info[0] = "This is the first information record.";
  info[1] = "";
  info[2] = "This info record is exactly 80 characters long.  last character should be pipe |";

  int num_info = 3;
  EXCHECK(ex_put_info(exoid, num_info, info));

  /* write results variables parameters and names */
  int num_glo_vars = 1;
  {
    char *var_names[] = {"glo_vars"};

    EXCHECK(ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars));
    EXCHECK(ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names));
  }

  int num_nod_vars = 2;
  {
    /*                    12345678901234567890123456789012 */
    char *var_names[] = {"node_variable_a_very_long_name_0", "nod_var1"};

    EXCHECK(ex_put_variable_param(exoid, EX_NODAL, num_nod_vars));
    EXCHECK(ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names));
  }

  int num_ele_vars = 3;
  {
    /* 0        1         2         3   */
    /* 12345678901234567890123456789012 */
    char *var_names[] = {"this_variable_name_is_short", "this_variable_name_is_just_right",
                         "this_variable_name_is_tooooo_long"};

    EXCHECK(ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars));
    EXCHECK(ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names));
  }

  int num_nset_vars = 3;
  {
    char *var_names[] = {"ns_var0", "ns_var1", "ns_var2"};

    EXCHECK(ex_put_variable_param(exoid, EX_NODE_SET, num_nset_vars));
    EXCHECK(ex_put_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names));
  }

  int num_sset_vars = 3;
  {
    char *var_names[] = {"ss_var0", "ss_var1", "ss_var2"};

    EXCHECK(ex_put_variable_param(exoid, EX_SIDE_SET, num_sset_vars));
    EXCHECK(ex_put_variable_names(exoid, EX_SIDE_SET, num_sset_vars, var_names));
  }

  /* write element variable truth table */
  int *truth_tab = (int *)calloc((num_elem_blk * num_ele_vars), sizeof(int));

  {
    int k = 0;
    for (int i = 0; i < num_elem_blk; i++) {
      for (int j = 0; j < num_ele_vars; j++) {
        truth_tab[k++] = 1;
      }
    }
  }

  EXCHECK(ex_put_truth_table(exoid, EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab));

  free(truth_tab);

  /* for each time step, write the analysis results;
   * the code below fills the arrays glob_var_vals,
   * nodal_var_vals, and elem_var_vals with values for debugging purposes;
   * obviously the analysis code will populate these arrays
   */

  double *glob_var_vals  = (double *)calloc(num_glo_vars, CPU_word_size);
  double *nodal_var_vals = (double *)calloc(num_nodes, CPU_word_size);
  double *elem_var_vals  = (double *)calloc(num_ele_vars, CPU_word_size);
  double *sset_var_vals  = (double *)calloc(10, CPU_word_size); /* max sides_in_sset */
  double *nset_var_vals  = (double *)calloc(5, CPU_word_size);  /* max nodes_in_nset */

  int num_time_steps = 10;
  for (int i = 0; i < num_time_steps; i++) {
    int    whole_time_step = i + 1;
    double time_value      = (double)(i + 1) / 100.;

    /* write time value */
    EXCHECK(ex_put_time(exoid, whole_time_step, &time_value));

    /* write global variables */
    for (int j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (double)(j + 2) * time_value;
    }
    EXCHECK(ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals));

    /* write nodal variables */
    for (int k = 1; k <= num_nod_vars; k++) {
      for (int j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (double)k + ((double)(j + 1) * time_value);
      }
      EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals));
    }

    /* write element variables */
    for (int k = 1; k <= num_ele_vars; k++) {
      for (int j = 0; j < num_elem_blk; j++) {
        for (int m = 0; m < blocks[j].num_entry; m++) {
          elem_var_vals[m] = (double)(k + 1) + (double)(j + 2) + ((double)(m + 1) * time_value);
          /* printf("elem_var_vals[%d]: %f\n",m,elem_var_vals[m]); */
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, blocks[j].id,
                           blocks[j].num_entry, elem_var_vals));
      }
    }

    /* write sideset variables */
    for (int k = 1; k <= num_sset_vars; k++) {
      for (int j = 0; j < num_side_sets; j++) {
        for (int m = 0; m < num_face_in_sset[j]; m++) {
          sset_var_vals[m] = (double)(k + 2) + (double)(j + 3) + ((double)(m + 1) * time_value);
          /* printf("sset_var_vals[%d]: %f\n",m,sset_var_vals[m]); */
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_SIDE_SET, k, ssids[j], num_face_in_sset[j],
                           sset_var_vals));
      }
    }

    /* write nodeset variables */
    for (int k = 1; k <= num_nset_vars; k++) {
      for (int j = 0; j < num_node_sets; j++) {
        for (int m = 0; m < num_nodes_in_nset[j]; m++) {
          nset_var_vals[m] = (double)(k + 3) + (double)(j + 4) + ((double)(m + 1) * time_value);
          /* printf("nset_var_vals[%d]: %f\n",m,nset_var_vals[m]); */
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODE_SET, k, nsids[j], num_nodes_in_nset[j],
                           nset_var_vals));
      }
    }

    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    EXCHECK(ex_update(exoid));
  }
  free(glob_var_vals);
  free(nodal_var_vals);
  free(elem_var_vals);
  free(sset_var_vals);
  free(nset_var_vals);

  /* close the EXODUS files
   */
  EXCHECK(ex_close(exoid));
  return 0;
}
