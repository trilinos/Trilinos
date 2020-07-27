/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwtm - test write many ExodusII database files
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
 *  database write routines. It tests multiple simultaneous output files.
 *
 *
 *****************************************************************************/

#define EX_TEST_INDIV_NODESET
#define EX_TEST_INDIV_SIDESET

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  exoidm[10], num_dim2, num_nodes2, num_elem2, num_elem_blk2;
  int  num_elem_in_block[10], num_node_sets, num_nodes_per_elem[10];
  int  num_elem_in_block2[10], num_node_sets2, num_nodes_per_elem2[10];
  int  num_side_sets, error;
  int  num_side_sets2, nexofiles = 5;
  int  i, j, k, m, n;
  int *elem_map, *connect, node_list[100], elem_list[100], side_list[100];
  int *elem_map2, *connect2, node_list2[100], elem_list2[100], side_list2[100];
  int  ebids[10], ids[10];
  int  ebids2[10], ids2[10];
  int  num_nodes_per_set[10], num_elem_per_set[10];
  int  num_nodes_per_set2[10], num_elem_per_set2[10];
  int  num_df_per_set[10], num_df_per_set2[10];
  int  df_ind[10], node_ind[10], elem_ind[10];
  int  df_ind2[10], node_ind2[10], elem_ind2[10];
  int  num_qa_rec, num_info;
  int  num_qa_rec2, num_info2;
  int  num_glo_vars, num_nod_vars, num_ele_vars;
  int  num_glo_vars2, num_nod_vars2, num_ele_vars2;
  int *truth_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;
  int  prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float  time_value;
  float  time_value2;
  float  x[100], y[100], z[100];
  float  attrib[1], dist_fact[1008];
  float  attrib2[1], dist_fact2[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char * coord_names2[3], *qa_record2[2][4], *info2[3], *var_names2[3];
  char   tmpstr[80];
  char * prop_names[2];
  char   exofname[256];

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II files */

  exoid = ex_create("test.exo",     /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  for (n = 0; n < nexofiles; n++) {
    sprintf(exofname, "test%d.exo", n);
    printf("test file name: %s\n", exofname);
    exoidm[n] = ex_create(exofname,       /* filename path */
                          EX_CLOBBER,     /* create mode */
                          &CPU_word_size, /* CPU float word size in bytes */
                          &IO_word_size); /* I/O float word size in bytes */
    printf("after ex_create for %s, exoid = %d\n", exofname, exoidm[n]);
  }

  /* initialize file with parameters */

  num_dim       = 3;
  num_nodes     = 26;
  num_elem      = 5;
  num_elem_blk  = 5;
  num_node_sets = 2;
  num_side_sets = 5;

  error = ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

  /* initialize file n with parameters */

  num_dim2       = 3;
  num_nodes2     = 26;
  num_elem2      = 5;
  num_elem_blk2  = 5;
  num_node_sets2 = 2;
  num_side_sets2 = 5;

  for (n = 0; n < nexofiles; n++) {
    sprintf(tmpstr, "This is test %d", n);

    error = ex_put_init(exoidm[n], tmpstr, num_dim2, num_nodes2, num_elem2, num_elem_blk2,
                        num_node_sets2, num_side_sets2);

    printf("after ex_put_init (%d), error = %d\n", n, error);
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

  error = ex_put_coord(exoid, x, y, z);
  printf("after ex_put_coord, error = %d\n", error);

  /* write nodal coordinates values and names to database */

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_coord(exoidm[n], x, y, z);
    printf("after ex_put_coord (%d), error = %d\n", n, error);
  }

  coord_names[0] = "xcoor";
  coord_names[1] = "ycoor";
  coord_names[2] = "zcoor";

  error = ex_put_coord_names(exoid, coord_names);
  printf("after ex_put_coord_names, error = %d\n", error);

  coord_names2[0] = "xcoor";
  coord_names2[1] = "ycoor";
  coord_names2[2] = "zcoor";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_coord_names(exoidm[n], coord_names2);
    printf("after ex_put_coord_names (%d), error = %d\n", n, error);
  }

  /* write element order map */

  elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i;
  }

  error = ex_put_map(exoid, elem_map);
  printf("after ex_put_map, error = %d\n", error);

  free(elem_map);

  elem_map2 = (int *)calloc(num_elem2, sizeof(int));

  for (i = 1; i <= num_elem2; i++) {
    elem_map2[i - 1] = i;
  }

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_map(exoidm[n], elem_map2);
    printf("after ex_put_map (%d), error = %d\n", n, error);
  }

  free(elem_map2);

  /* write element block parameters */

  num_elem_in_block[0] = 1;
  num_elem_in_block[1] = 1;
  num_elem_in_block[2] = 1;
  num_elem_in_block[3] = 1;
  num_elem_in_block[4] = 1;

  num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
  num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
  num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
  num_nodes_per_elem[3] = 4; /* elements in block #3 are 4-node tetras */
  num_nodes_per_elem[4] = 6; /* elements in block #3 are 6-node wedges */

  ebids[0] = 10;
  ebids[1] = 11;
  ebids[2] = 12;
  ebids[3] = 13;
  ebids[4] = 14;

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "quad", num_elem_in_block[0],
                       num_nodes_per_elem[0], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[1], "quad", num_elem_in_block[1],
                       num_nodes_per_elem[1], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[2], "hex", num_elem_in_block[2],
                       num_nodes_per_elem[2], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[3], "tetra", num_elem_in_block[3],
                       num_nodes_per_elem[3], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[4], "wedge", num_elem_in_block[4],
                       num_nodes_per_elem[4], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  /* write element block properties */

  prop_names[0] = "MATL";
  prop_names[1] = "DENSITY";
  error         = ex_put_prop_names(exoid, EX_ELEM_BLOCK, 2, prop_names);
  printf("after ex_put_prop_names, error = %d\n", error);

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], "MATL", 10);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], "MATL", 20);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], "MATL", 30);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[3], "MATL", 40);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[4], "MATL", 50);
  printf("after ex_put_prop, error = %d\n", error);

  /* files n */

  num_elem_in_block2[0] = 1;
  num_elem_in_block2[1] = 1;
  num_elem_in_block2[2] = 1;
  num_elem_in_block2[3] = 1;
  num_elem_in_block2[4] = 1;

  num_nodes_per_elem2[0] = 4; /* elements in block #1 are 4-node quads  */
  num_nodes_per_elem2[1] = 4; /* elements in block #2 are 4-node quads  */
  num_nodes_per_elem2[2] = 8; /* elements in block #3 are 8-node hexes  */
  num_nodes_per_elem2[3] = 4; /* elements in block #3 are 4-node tetras */
  num_nodes_per_elem2[4] = 6; /* elements in block #3 are 6-node wedges */

  ebids2[0] = 10;
  ebids2[1] = 11;
  ebids2[2] = 12;
  ebids2[3] = 13;
  ebids2[4] = 14;

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_block(exoidm[n], EX_ELEM_BLOCK, ebids2[0], "quad", num_elem_in_block2[0],
                         num_nodes_per_elem2[0], 0, 0, 1);
    printf("after ex_put_elem_block (%d), error = %d\n", n, error);

    error = ex_put_block(exoidm[n], EX_ELEM_BLOCK, ebids2[1], "quad", num_elem_in_block2[1],
                         num_nodes_per_elem2[1], 0, 0, 1);
    printf("after ex_put_elem_block (%d), error = %d\n", n, error);

    error = ex_put_block(exoidm[n], EX_ELEM_BLOCK, ebids2[2], "hex", num_elem_in_block2[2],
                         num_nodes_per_elem2[2], 0, 0, 1);
    printf("after ex_put_elem_block (%d), error = %d\n", n, error);

    error = ex_put_block(exoidm[n], EX_ELEM_BLOCK, ebids2[3], "tetra", num_elem_in_block2[3],
                         num_nodes_per_elem2[3], 0, 0, 1);
    printf("after ex_put_elem_block (%d), error = %d\n", n, error);

    error = ex_put_block(exoidm[n], EX_ELEM_BLOCK, ebids2[4], "wedge", num_elem_in_block2[4],
                         num_nodes_per_elem2[4], 0, 0, 1);
    printf("after ex_put_elem_block (%d), error = %d\n", n, error);

    /* write element block properties */

    prop_names[0] = "MATL";
    prop_names[1] = "DENSITY";
    error         = ex_put_prop_names(exoidm[n], EX_ELEM_BLOCK, 2, prop_names);
    printf("after ex_put_prop_names (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_ELEM_BLOCK, ebids2[0], "MATL", 100);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
    error = ex_put_prop(exoidm[n], EX_ELEM_BLOCK, ebids2[1], "MATL", 200);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
    error = ex_put_prop(exoidm[n], EX_ELEM_BLOCK, ebids2[2], "MATL", 300);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
    error = ex_put_prop(exoidm[n], EX_ELEM_BLOCK, ebids2[3], "MATL", 400);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
    error = ex_put_prop(exoidm[n], EX_ELEM_BLOCK, ebids2[4], "MATL", 500);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
  }

  /* write element connectivity */

  connect    = (int *)calloc(8, sizeof(int));
  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 5;
  connect[1] = 6;
  connect[2] = 7;
  connect[3] = 8;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[1], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 9;
  connect[1] = 10;
  connect[2] = 11;
  connect[3] = 12;
  connect[4] = 13;
  connect[5] = 14;
  connect[6] = 15;
  connect[7] = 16;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[2], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[3], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 21;
  connect[1] = 22;
  connect[2] = 23;
  connect[3] = 24;
  connect[4] = 25;
  connect[5] = 26;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[4], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  free(connect);

  for (n = 0; n < nexofiles; n++) {
    connect2    = (int *)calloc(8, sizeof(int));
    connect2[0] = 1;
    connect2[1] = 2;
    connect2[2] = 3;
    connect2[3] = 4;

    error = ex_put_conn(exoidm[n], EX_ELEM_BLOCK, ebids[0], connect2, NULL, NULL);
    printf("after ex_put_elem_conn (%d), error = %d\n", n, error);

    connect2[0] = 5;
    connect2[1] = 6;
    connect2[2] = 7;
    connect2[3] = 8;

    error = ex_put_conn(exoidm[n], EX_ELEM_BLOCK, ebids[1], connect2, NULL, NULL);
    printf("after ex_put_elem_conn (%d), error = %d\n", n, error);

    connect2[0] = 9;
    connect2[1] = 10;
    connect2[2] = 11;
    connect2[3] = 12;
    connect2[4] = 13;
    connect2[5] = 14;
    connect2[6] = 15;
    connect2[7] = 16;

    error = ex_put_conn(exoidm[n], EX_ELEM_BLOCK, ebids2[2], connect2, NULL, NULL);
    printf("after ex_put_elem_conn (%d), error = %d\n", n, error);

    connect2[0] = 17;
    connect2[1] = 18;
    connect2[2] = 19;
    connect2[3] = 20;

    error = ex_put_conn(exoidm[n], EX_ELEM_BLOCK, ebids2[3], connect2, NULL, NULL);
    printf("after ex_put_elem_conn (%d), error = %d\n", n, error);

    connect2[0] = 21;
    connect2[1] = 22;
    connect2[2] = 23;
    connect2[3] = 24;
    connect2[4] = 25;
    connect2[5] = 26;

    error = ex_put_conn(exoidm[n], EX_ELEM_BLOCK, ebids2[4], connect2, NULL, NULL);
    printf("after ex_put_elem_conn (%d), error = %d\n", n, error);

    free(connect2);
  }

  /* write element block attributes */

  attrib[0] = 3.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  attrib[0] = 6.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[1], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[2], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[3], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[4], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  for (n = 0; n < nexofiles; n++) {
    attrib2[0] = 3.;
    error      = ex_put_attr(exoidm[n], EX_ELEM_BLOCK, ebids[0], attrib2);
    printf("after ex_put_elem_attr (%d), error = %d\n", n, error);

    attrib2[0] = 6.;
    error      = ex_put_attr(exoidm[n], EX_ELEM_BLOCK, ebids[1], attrib2);
    printf("after ex_put_elem_attr (%d), error = %d\n", n, error);

    error = ex_put_attr(exoidm[n], EX_ELEM_BLOCK, ebids[2], attrib2);
    printf("after ex_put_elem_attr (%d), error = %d\n", n, error);

    error = ex_put_attr(exoidm[n], EX_ELEM_BLOCK, ebids[3], attrib2);
    printf("after ex_put_elem_attr (%d), error = %d\n", n, error);

    error = ex_put_attr(exoidm[n], EX_ELEM_BLOCK, ebids[4], attrib2);
    printf("after ex_put_elem_attr (%d), error = %d\n", n, error);
  }

#ifdef EX_TEST_INDIV_NODESET
  /* write individual node sets */

  error = ex_put_set_param(exoid, EX_NODE_SET, 20, 5, 5);
  printf("after ex_put_node_set_param, error = %d\n", error);

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
  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, 20, dist_fact);
  printf("after ex_put_node_set, error = %d\n", error);

  error = ex_put_set_param(exoid, EX_NODE_SET, 21, 3, 3);
  printf("after ex_put_node_set_param, error = %d\n", error);

  node_list[0] = 20;
  node_list[1] = 21;
  node_list[2] = 22;

  dist_fact[0] = 1.1;
  dist_fact[1] = 2.1;
  dist_fact[2] = 3.1;

  error = ex_put_set(exoid, EX_NODE_SET, 21, node_list, NULL);
  printf("after ex_put_node_set, error = %d\n", error);
  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, 21, dist_fact);
  printf("after ex_put_node_set, error = %d\n", error);

  error = ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4);
  printf("after ex_put_prop, error = %d\n", error);

  error = ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5);
  printf("after ex_put_prop, error = %d\n", error);

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  error = ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array);
  printf("after ex_put_prop_array, error = %d\n", error);

  /* file 2 */
  for (n = 0; n < nexofiles; n++) {
    error = ex_put_set_param(exoidm[n], EX_NODE_SET, 20, 5, 5);
    printf("after ex_put_node_set_param (%d), error = %d\n", n, error);

    node_list2[0] = 10;
    node_list2[1] = 11;
    node_list2[2] = 12;
    node_list2[3] = 13;
    node_list2[4] = 14;

    dist_fact2[0] = 1.0;
    dist_fact2[1] = 2.0;
    dist_fact2[2] = 3.0;
    dist_fact2[3] = 4.0;
    dist_fact2[4] = 5.0;

    error = ex_put_set(exoidm[n], EX_NODE_SET, 20, node_list2, NULL);
    printf("after ex_put_node_set (%d), error = %d\n", n, error);
    error = ex_put_set_dist_fact(exoidm[n], EX_NODE_SET, 20, dist_fact2);
    printf("after ex_put_node_set (%d), error = %d\n", n, error);

    error = ex_put_set_param(exoidm[n], EX_NODE_SET, 21, 3, 3);
    printf("after ex_put_node_set_param (%d), error = %d\n", n, error);

    node_list2[0] = 20;
    node_list2[1] = 21;
    node_list2[2] = 22;

    dist_fact2[0] = 1.1;
    dist_fact2[1] = 2.1;
    dist_fact2[2] = 3.1;

    error = ex_put_set(exoidm[n], EX_NODE_SET, 21, node_list2, NULL);
    printf("after ex_put_node_set (%d), error = %d\n", n, error);
    error = ex_put_set_dist_fact(exoidm[n], EX_NODE_SET, 21, dist_fact2);
    printf("after ex_put_node_set (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_NODE_SET, 20, "FACE", 4);
    printf("after ex_put_prop (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_NODE_SET, 21, "FACE", 5);
    printf("after ex_put_prop (%d), error = %d\n", n, error);

    prop_array[0] = 1000;
    prop_array[1] = 2000;

    error = ex_put_prop_array(exoidm[n], EX_NODE_SET, "VELOCITY", prop_array);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
  }

#else /* EX_TEST_INDIV_NODESET */
  /* write concatenated node sets; this produces the same information as
   * the above code which writes individual node sets
   */

  ids[0] = 20;
  ids[1] = 21;

  num_nodes_per_set[0] = 5;
  num_nodes_per_set[1] = 3;

  node_ind[0] = 0;
  node_ind[1] = 5;

  node_list[0] = 10;
  node_list[1] = 11;
  node_list[2] = 12;
  node_list[3] = 13;
  node_list[4] = 14;
  node_list[5] = 20;
  node_list[6] = 21;
  node_list[7] = 22;

  num_df_per_set[0] = 5;
  num_df_per_set[1] = 3;

  df_ind[0] = 0;
  df_ind[1] = 5;

  dist_fact[0] = 1.0;
  dist_fact[1] = 2.0;
  dist_fact[2] = 3.0;
  dist_fact[3] = 4.0;
  dist_fact[4] = 5.0;
  dist_fact[5] = 1.1;
  dist_fact[6] = 2.1;
  dist_fact[7] = 3.1;

  error = ex_put_concat_node_sets(exoid, ids, num_nodes_per_set, node_ind, node_list, dist_fact);
  printf("after ex_put_concat_node_sets, error = %d\n", error);

  error = ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5);
  printf("after ex_put_prop, error = %d\n", error);

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  error = ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array);
  printf("after ex_put_prop_array, error = %d\n", error);

  ids2[0] = 20;
  ids2[1] = 21;

  num_nodes_per_set2[0] = 5;
  num_nodes_per_set2[1] = 3;

  node_ind2[0] = 0;
  node_ind2[1] = 5;

  node_list2[0] = 10;
  node_list2[1] = 11;
  node_list2[2] = 12;
  node_list2[3] = 13;
  node_list2[4] = 14;
  node_list2[5] = 20;
  node_list2[6] = 21;
  node_list2[7] = 22;

  num_df_per_set2[0] = 5;
  num_df_per_set2[1] = 3;

  df_ind2[0] = 0;
  df_ind2[1] = 5;

  dist_fact2[0] = 1.0;
  dist_fact2[1] = 2.0;
  dist_fact2[2] = 3.0;
  dist_fact2[3] = 4.0;
  dist_fact2[4] = 5.0;
  dist_fact2[5] = 1.1;
  dist_fact2[6] = 2.1;
  dist_fact2[7] = 3.1;

  prop_array2[0] = 1000;
  prop_array2[1] = 2000;

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_concat_node_sets(exoidm[n], ids2, num_nodes_per_set2, num_df_per_set2, node_ind2,
                                    df_ind2, node_list2, dist_fact2);
    printf("after ex_put_concat_node_sets, error = %d\n", error);

    error = ex_put_prop(exoidm[n], EX_NODE_SET, 20, "FACE", 4);
    printf("after ex_put_prop, error = %d\n", error);
    error = ex_put_prop(exoidm[n], EX_NODE_SET, 21, "FACE", 5);
    printf("after ex_put_prop, error = %d\n", error);

    error = ex_put_prop_array(exoidm[n], EX_NODE_SET, "VELOCITY", prop_array2);
    printf("after ex_put_prop_array, error = %d\n", error);
  }

#endif /* EX_TEST_INDIV_NODESET */

#ifdef TEST_INDIV_SIDESET
  /* write individual side sets */

  /* side set #1  - quad */

  error = ex_put_set_param(exoid, EX_SIDE_SET, 30, 2, 4);
  printf("after ex_put_side_set_param, error = %d\n", error);

  elem_list[0] = 2;
  elem_list[1] = 2;

  side_list[0] = 4;
  side_list[1] = 2;

  dist_fact[0] = 30.0;
  dist_fact[1] = 30.1;
  dist_fact[2] = 30.2;
  dist_fact[3] = 30.3;

  error = ex_put_set(exoid, EX_SIDE_SET, 30, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  error = ex_put_set_dist_fact(exoid, EX_SIDE_SET, 30, dist_fact);
  printf("after ex_put_side_set_dist_fact, error = %d\n", error);

  /* side set #2  - quad  spanning elements */

  error = ex_put_set_param(exoid, EX_SIDE_SET, 31, 2, 4);
  printf("after ex_put_side_set_param, error = %d\n", error);

  elem_list[0] = 1;
  elem_list[1] = 2;

  side_list[0] = 2;
  side_list[1] = 3;

  dist_fact[0] = 31.0;
  dist_fact[1] = 31.1;
  dist_fact[2] = 31.2;
  dist_fact[3] = 31.3;

  error = ex_put_set(exoid, EX_SIDE_SET, 31, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  error = ex_put_set_dist_fact(exoid, EX_SIDE_SET, 31, dist_fact);
  printf("after ex_put_side_set_dist_fact, error = %d\n", error);

  /* side set #3  - hex */

  error = ex_put_set_param(exoid, EX_SIDE_SET, 32, 7, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

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

  error = ex_put_set(exoid, EX_SIDE_SET, 32, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  /* side set #4  - tetras */

  error = ex_put_set_param(exoid, EX_SIDE_SET, 33, 4, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

  elem_list[0] = 4;
  elem_list[1] = 4;
  elem_list[2] = 4;
  elem_list[3] = 4;

  side_list[0] = 1;
  side_list[1] = 2;
  side_list[2] = 3;
  side_list[3] = 4;

  error = ex_put_set(exoid, EX_SIDE_SET, 33, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  /* side set #5  - wedges */

  error = ex_put_set_param(exoid, EX_SIDE_SET, 34, 5, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

  elem_list[0] = 5;
  elem_list[1] = 5;
  elem_list[2] = 5;
  elem_list[3] = 5;
  elem_list[4] = 5;

  side_list[0] = 1;
  side_list[1] = 2;
  side_list[2] = 3;
  side_list[3] = 4;
  side_list[4] = 5;

  error = ex_put_set(exoid, EX_SIDE_SET, 34, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  error = ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100);
  printf("after ex_put_prop, error = %d\n", error);

  error = ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101);
  printf("after ex_put_prop, error = %d\n", error);

  /* file 2 */

  for (n = 0; n < nexofiles; n++) {
    /* side set 1 */

    error = ex_put_set_param(exoidm[n], EX_SIDE_SET, 30, 2, 4);
    printf("after ex_put_side_set_param (%d), error = %d\n", n, error);

    elem_list2[0] = 2;
    elem_list2[1] = 2;

    side_list2[0] = 4;
    side_list2[1] = 2;

    dist_fact2[0] = 30.0;
    dist_fact2[1] = 30.1;
    dist_fact2[2] = 30.2;
    dist_fact2[3] = 30.3;

    error = ex_put_set(exoidm[n], EX_SIDE_SET, 30, elem_list2, side_list2);
    printf("after ex_put_side_set (%d), error = %d\n", n, error);

    error = ex_put_set_dist_fact(exoidm[n], EX_SIDE_SET, 30, dist_fact2);
    printf("after ex_put_side_set_dist_fact (%d), error = %d\n", n, error);

    /* side set 2 */

    error = ex_put_set_param(exoidm[n], EX_SIDE_SET, 31, 2, 4);
    printf("after ex_put_side_set_param (%d), error = %d\n", n, error);

    elem_list2[0] = 1;
    elem_list2[1] = 2;

    side_list2[0] = 2;
    side_list2[1] = 3;

    dist_fact2[0] = 31.0;
    dist_fact2[1] = 31.1;
    dist_fact2[2] = 31.2;
    dist_fact2[3] = 31.3;

    error = ex_put_set(exoidm[n], EX_SIDE_SET, 31, elem_list2, side_list2);
    printf("after ex_put_side_set (%d), error = %d\n", n, error);

    error = ex_put_set_dist_fact(exoidm[n], EX_SIDE_SET, 31, dist_fact2);
    printf("after ex_put_side_set_dist_fact (%d), error = %d\n", n, error);

    /* side set #3  - hex */

    error = ex_put_set_param(exoidm[n], EX_SIDE_SET, 32, 7, 0);
    printf("after ex_put_side_set_param (%d), error = %d\n", n, error);

    elem_list2[0] = 3;
    elem_list2[1] = 3;
    elem_list2[2] = 3;
    elem_list2[3] = 3;
    elem_list2[4] = 3;
    elem_list2[5] = 3;
    elem_list2[6] = 3;

    side_list2[0] = 5;
    side_list2[1] = 3;
    side_list2[2] = 3;
    side_list2[3] = 2;
    side_list2[4] = 4;
    side_list2[5] = 1;
    side_list2[6] = 6;

    error = ex_put_set(exoidm[n], EX_SIDE_SET, 32, elem_list2, side_list2);
    printf("after ex_put_side_set (%d), error = %d\n", n, error);

    /* side set #4  - tetras */

    error = ex_put_set_param(exoidm[n], EX_SIDE_SET, 33, 4, 0);
    printf("after ex_put_side_set_param (%d), error = %d\n", n, error);

    elem_list2[0] = 4;
    elem_list2[1] = 4;
    elem_list2[2] = 4;
    elem_list2[3] = 4;

    side_list2[0] = 1;
    side_list2[1] = 2;
    side_list2[2] = 3;
    side_list2[3] = 4;

    error = ex_put_set(exoidm[n], EX_SIDE_SET, 33, elem_list2, side_list2);
    printf("after ex_put_side_set (%d), error = %d\n", n, error);

    /* side set #5  - wedges */

    error = ex_put_set_param(exoidm[n], EX_SIDE_SET, 34, 5, 0);
    printf("after ex_put_side_set_param (%d), error = %d\n", n, error);

    elem_list2[0] = 5;
    elem_list2[1] = 5;
    elem_list2[2] = 5;
    elem_list2[3] = 5;
    elem_list2[4] = 5;

    side_list2[0] = 1;
    side_list2[1] = 2;
    side_list2[2] = 3;
    side_list2[3] = 4;
    side_list2[4] = 5;

    error = ex_put_set(exoidm[n], EX_SIDE_SET, 34, elem_list2, side_list2);
    printf("after ex_put_side_set (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_SIDE_SET, 30, "COLOR", 100);
    printf("after ex_put_prop (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_SIDE_SET, 31, "COLOR", 101);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
  }

#else /* TEST_INDIV_SIDESET */
  /* write concatenated side sets; this produces the same information as
   * the above code which writes individual side sets
   */

  ids[0] = 30;
  ids[1] = 31;
  ids[2] = 32;
  ids[3] = 33;
  ids[4] = 34;

  node_list[0] = 8;
  node_list[1] = 5;
  node_list[2] = 6;
  node_list[3] = 7;

  node_list[4] = 2;
  node_list[5] = 3;
  node_list[6] = 7;
  node_list[7] = 8;

  node_list[8]  = 9;
  node_list[9]  = 12;
  node_list[10] = 11;
  node_list[11] = 10;

  node_list[12] = 11;
  node_list[13] = 12;
  node_list[14] = 16;
  node_list[15] = 15;

  node_list[16] = 16;
  node_list[17] = 15;
  node_list[18] = 11;
  node_list[19] = 12;

  node_list[20] = 10;
  node_list[21] = 11;
  node_list[22] = 15;
  node_list[23] = 14;

  node_list[24] = 13;
  node_list[25] = 16;
  node_list[26] = 12;
  node_list[27] = 9;

  node_list[28] = 14;
  node_list[29] = 13;
  node_list[30] = 9;
  node_list[31] = 10;

  node_list[32] = 16;
  node_list[33] = 13;
  node_list[34] = 14;
  node_list[35] = 15;

  node_list[36] = 17;
  node_list[37] = 18;
  node_list[38] = 20;

  node_list[39] = 18;
  node_list[40] = 19;
  node_list[41] = 20;

  node_list[42] = 20;
  node_list[43] = 19;
  node_list[44] = 17;

  node_list[45] = 19;
  node_list[46] = 18;
  node_list[47] = 17;

  node_list[48] = 25;
  node_list[49] = 24;
  node_list[50] = 21;
  node_list[51] = 22;

  node_list[52] = 26;
  node_list[53] = 25;
  node_list[54] = 22;
  node_list[55] = 23;

  node_list[56] = 26;
  node_list[57] = 23;
  node_list[58] = 21;
  node_list[59] = 24;

  node_list[60] = 23;
  node_list[61] = 22;
  node_list[62] = 21;

  node_list[63] = 24;
  node_list[64] = 25;
  node_list[65] = 26;

  node_ind[0] = 0;
  node_ind[1] = 4;
  node_ind[2] = 8;
  node_ind[3] = 36;
  node_ind[4] = 47;

  num_elem_per_set[0] = 2;
  num_elem_per_set[1] = 2;
  num_elem_per_set[2] = 7;
  num_elem_per_set[3] = 4;
  num_elem_per_set[4] = 5;

  num_nodes_per_set[0] = 4;
  num_nodes_per_set[1] = 4;
  num_nodes_per_set[2] = 28;
  num_nodes_per_set[3] = 12;
  num_nodes_per_set[4] = 18;

  elem_ind[0] = 0;
  elem_ind[1] = 2;
  elem_ind[2] = 4;
  elem_ind[3] = 11;
  elem_ind[4] = 15;

  elem_list[0]  = 2;
  elem_list[1]  = 2;
  elem_list[2]  = 1;
  elem_list[3]  = 2;
  elem_list[4]  = 3;
  elem_list[5]  = 3;
  elem_list[6]  = 3;
  elem_list[7]  = 3;
  elem_list[8]  = 3;
  elem_list[9]  = 3;
  elem_list[10] = 3;
  elem_list[11] = 4;
  elem_list[12] = 4;
  elem_list[13] = 4;
  elem_list[14] = 4;
  elem_list[15] = 5;
  elem_list[16] = 5;
  elem_list[17] = 5;
  elem_list[18] = 5;
  elem_list[19] = 5;

  error = ex_cvt_nodes_to_sides(exoid, num_elem_per_set, num_nodes_per_set, elem_ind, node_ind,
                                elem_list, node_list, side_list);
  printf("after ex_cvt_nodes_to_sides, error = %d\n", error);

  num_df_per_set[0] = 4;
  num_df_per_set[1] = 4;
  num_df_per_set[2] = 0;
  num_df_per_set[3] = 0;
  num_df_per_set[4] = 0;

  df_ind[0] = 0;
  df_ind[1] = 4;

  dist_fact[0] = 30.0;
  dist_fact[1] = 30.1;
  dist_fact[2] = 30.2;
  dist_fact[3] = 30.3;

  dist_fact[4] = 31.0;
  dist_fact[5] = 31.1;
  dist_fact[6] = 31.2;
  dist_fact[7] = 31.3;

  {
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = ids;
    set_specs.num_entries_per_set = num_elem_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = elem_ind;
    set_specs.sets_dist_index     = df_ind;
    set_specs.sets_entry_list     = elem_list;
    set_specs.sets_extra_list     = side_list;
    set_specs.sets_dist_fact      = dist_fact;
    error                         = ex_put_concat_sets(exoid, EX_SIDE_SET, &set_specs);
  }

  printf("after ex_put_concat_side_sets, error = %d\n", error);

  error = ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100);
  printf("after ex_put_prop, error = %d\n", error);

  error = ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101);
  printf("after ex_put_prop, error = %d\n", error);

  /* file 2 */

  ids2[0] = 30;
  ids2[1] = 31;
  ids2[2] = 32;
  ids2[3] = 33;
  ids2[4] = 34;

  node_list2[0] = 8;
  node_list2[1] = 5;
  node_list2[2] = 6;
  node_list2[3] = 7;

  node_list2[4] = 2;
  node_list2[5] = 3;
  node_list2[6] = 7;
  node_list2[7] = 8;

  node_list2[8]  = 9;
  node_list2[9]  = 12;
  node_list2[10] = 11;
  node_list2[11] = 10;

  node_list2[12] = 11;
  node_list2[13] = 12;
  node_list2[14] = 16;
  node_list2[15] = 15;

  node_list2[16] = 16;
  node_list2[17] = 15;
  node_list2[18] = 11;
  node_list2[19] = 12;

  node_list2[20] = 10;
  node_list2[21] = 11;
  node_list2[22] = 15;
  node_list2[23] = 14;

  node_list2[24] = 13;
  node_list2[25] = 16;
  node_list2[26] = 12;
  node_list2[27] = 9;

  node_list2[28] = 14;
  node_list2[29] = 13;
  node_list2[30] = 9;
  node_list2[31] = 10;

  node_list2[32] = 16;
  node_list2[33] = 13;
  node_list2[34] = 14;
  node_list2[35] = 15;

  node_list2[36] = 17;
  node_list2[37] = 18;
  node_list2[38] = 20;

  node_list2[39] = 18;
  node_list2[40] = 19;
  node_list2[41] = 20;

  node_list2[42] = 20;
  node_list2[43] = 19;
  node_list2[44] = 17;

  node_list2[45] = 19;
  node_list2[46] = 18;
  node_list2[47] = 17;

  node_list2[48] = 25;
  node_list2[49] = 24;
  node_list2[50] = 21;
  node_list2[51] = 22;

  node_list2[52] = 26;
  node_list2[53] = 25;
  node_list2[54] = 22;
  node_list2[55] = 23;

  node_list2[56] = 26;
  node_list2[57] = 23;
  node_list2[58] = 21;
  node_list2[59] = 24;

  node_list2[60] = 23;
  node_list2[61] = 22;
  node_list2[62] = 21;

  node_list2[63] = 24;
  node_list2[64] = 25;
  node_list2[65] = 26;

  node_ind2[0] = 0;
  node_ind2[1] = 4;
  node_ind2[2] = 8;
  node_ind2[3] = 36;
  node_ind2[4] = 47;

  num_elem_per_set2[0] = 2;
  num_elem_per_set2[1] = 2;
  num_elem_per_set2[2] = 7;
  num_elem_per_set2[3] = 4;
  num_elem_per_set2[4] = 5;

  num_nodes_per_set2[0] = 4;
  num_nodes_per_set2[1] = 4;
  num_nodes_per_set2[2] = 28;
  num_nodes_per_set2[3] = 12;
  num_nodes_per_set2[4] = 18;

  elem_ind2[0] = 0;
  elem_ind2[1] = 2;
  elem_ind2[2] = 4;
  elem_ind2[3] = 11;
  elem_ind2[4] = 15;

  elem_list2[0]  = 2;
  elem_list2[1]  = 2;
  elem_list2[2]  = 1;
  elem_list2[3]  = 2;
  elem_list2[4]  = 3;
  elem_list2[5]  = 3;
  elem_list2[6]  = 3;
  elem_list2[7]  = 3;
  elem_list2[8]  = 3;
  elem_list2[9]  = 3;
  elem_list2[10] = 3;
  elem_list2[11] = 4;
  elem_list2[12] = 4;
  elem_list2[13] = 4;
  elem_list2[14] = 4;
  elem_list2[15] = 5;
  elem_list2[16] = 5;
  elem_list2[17] = 5;
  elem_list2[18] = 5;
  elem_list2[19] = 5;

  num_df_per_set2[0] = 4;
  num_df_per_set2[1] = 4;
  num_df_per_set2[2] = 0;
  num_df_per_set2[3] = 0;
  num_df_per_set2[4] = 0;

  df_ind2[0] = 0;
  df_ind2[1] = 4;

  dist_fact2[0] = 30.0;
  dist_fact2[1] = 30.1;
  dist_fact2[2] = 30.2;
  dist_fact2[3] = 30.3;

  dist_fact2[4] = 31.0;
  dist_fact2[5] = 31.1;
  dist_fact2[6] = 31.2;
  dist_fact2[7] = 31.3;

  for (n = 0; n < nexofiles; n++) {
    error = ex_cvt_nodes_to_sides(exoidm[n], num_elem_per_set2, num_nodes_per_set2, elem_ind2,
                                  node_ind2, elem_list2, node_list2, side_list2);
    printf("after ex_cvt_nodes_to_sides (%d), error = %d\n", n, error);

    {
      struct ex_set_specs set_specs;

      set_specs.sets_ids            = ids2;
      set_specs.num_entries_per_set = num_elem_per_set2;
      set_specs.num_dist_per_set    = num_df_per_set2;
      set_specs.sets_entry_index    = elem_ind2;
      set_specs.sets_dist_index     = df_ind2;
      set_specs.sets_entry_list     = elem_list2;
      set_specs.sets_extra_list     = side_list2;
      set_specs.sets_dist_fact      = dist_fact2;
      error                         = ex_put_concat_sets(exoid, EX_SIDE_SET, &set_specs);
    }
    printf("after ex_put_concat_side_sets (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_SIDE_SET, 30, "COLOR", 100);
    printf("after ex_put_prop (%d), error = %d\n", n, error);

    error = ex_put_prop(exoidm[n], EX_SIDE_SET, 31, "COLOR", 101);
    printf("after ex_put_prop (%d), error = %d\n", n, error);
  }

  /* END COMMENTED OUT SECTION */

#endif /* TEST_INDIV_SIDESET */

  /* write QA records */

  num_qa_rec = 2;

  qa_record[0][0] = "TESTWTM";
  qa_record[0][1] = "testwtm";
  qa_record[0][2] = "07/07/93";
  qa_record[0][3] = "15:41:33";
  qa_record[1][0] = "FASTQ";
  qa_record[1][1] = "fastq";
  qa_record[1][2] = "07/07/93";
  qa_record[1][3] = "16:41:33";

  error = ex_put_qa(exoid, num_qa_rec, qa_record);
  printf("after ex_put_qa, error = %d\n", error);

  num_qa_rec2 = 2;

  qa_record2[0][0] = "TESTWTM";
  qa_record2[0][1] = "testwtm";
  qa_record2[0][2] = "07/07/93";
  qa_record2[0][3] = "15:41:33";
  qa_record2[1][0] = "FASTQ";
  qa_record2[1][1] = "fastq";
  qa_record2[1][2] = "07/07/93";
  qa_record2[1][3] = "16:41:33";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_qa(exoidm[n], num_qa_rec2, qa_record2);
    printf("after ex_put_qa (%d), error = %d\n", n, error);
  }

  /* write information records */

  num_info = 3;

  info[0] = "This is the first information record.";
  info[1] = "This is the second information record.";
  info[2] = "This is the third information record.";

  error = ex_put_info(exoid, num_info, info);
  printf("after ex_put_info, error = %d\n", error);

  num_info2 = 3;

  info2[0] = "This is the first information record.";
  info2[1] = "This is the second information record.";
  info2[2] = "This is the third information record.";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_info(exoidm[n], num_info2, info2);
    printf("after ex_put_info (%d), error = %d\n", n, error);
  }

  /* write results variables parameters and names */

  num_glo_vars = 1;

  var_names[0] = "glo_vars";

  error = ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);

  num_glo_vars2 = 1;

  var_names2[0] = "glo_vars";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_variable_param(exoidm[n], EX_GLOBAL, num_glo_vars2);
    printf("after ex_put_variable_param (%d), error = %d\n", n, error);
    error = ex_put_variable_names(exoidm[n], EX_GLOBAL, num_glo_vars2, var_names2);
    printf("after ex_put_variable_names (%d), error = %d\n", n, error);
  }

  num_nod_vars = 2;

  var_names[0] = "nod_var0";
  var_names[1] = "nod_var1";

  error = ex_put_variable_param(exoid, EX_NODAL, num_nod_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);

  num_nod_vars2 = 2;

  var_names2[0] = "nod_var0";
  var_names2[1] = "nod_var1";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_variable_param(exoidm[n], EX_NODAL, num_nod_vars2);
    printf("after ex_put_variable_param (%d), error = %d\n", n, error);
    error = ex_put_variable_names(exoidm[n], EX_NODAL, num_nod_vars2, var_names2);
    printf("after ex_put_variable_names (%d), error = %d\n", n, error);
  }

  num_ele_vars = 3;

  var_names[0] = "ele_var0";
  var_names[1] = "ele_var1";
  var_names[2] = "ele_var2";

  error = ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);

  num_ele_vars2 = 3;

  var_names2[0] = "ele_var20";
  var_names2[1] = "ele_var21";
  var_names2[2] = "ele_var22";

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_variable_param(exoidm[n], EX_ELEM_BLOCK, num_ele_vars2);
    printf("after ex_put_variable_param (%d), error = %d\n", n, error);
    error = ex_put_variable_names(exoidm[n], EX_ELEM_BLOCK, num_ele_vars, var_names);
    printf("after ex_put_variable_names (%d), error = %d\n", n, error);
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

  for (n = 0; n < nexofiles; n++) {
    error = ex_put_truth_table(exoidm[n], EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab);
    printf("after ex_put_elem_var_tab (%d), error = %d\n", n, error);
  }

  free(truth_tab);

  /* for each time step, write the analysis results;
   * the code below fills the arrays  glob_var_vals,
   * nodal_var_vals, and elem_var_vals with values for debugging purposes;
   * obviously the analysis code will populate these arrays
   */

  whole_time_step = 1;
  num_time_steps  = 10;

  glob_var_vals  = (float *)calloc(num_glo_vars, sizeof(float));
  nodal_var_vals = (float *)calloc(num_nodes, sizeof(float));
  elem_var_vals  = (float *)calloc(4, sizeof(float));

  for (i = 0; i < num_time_steps; i++) {
    time_value  = (float)(i + 1) / 100.;
    time_value2 = (float)(i + 1) / 100.;

    /* write time value */

    error = ex_put_time(exoid, whole_time_step, &time_value);
    printf("after ex_put_time, error = %d\n", error);

    for (n = 0; n < nexofiles; n++) {
      error = ex_put_time(exoidm[n], whole_time_step, &time_value2);
      printf("after ex_put_time (%d), error = %d\n", n, error);
    }

    /* write global variables */

    for (j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (float)(j + 2) * time_value;
    }

    error = ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals);
    printf("after ex_put_glob_vars, error = %d\n", error);

    for (n = 0; n < nexofiles; n++) {
      error = ex_put_var(exoidm[n], whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals);
      printf("after ex_put_glob_vars (%d), error = %d\n", n, error);
    }
    /* write nodal variables */

    for (k = 1; k <= num_nod_vars; k++) {
      for (j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (float)k + ((float)(j + 1) * time_value);
      }

      error = ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals);
      printf("after ex_put_nodal_var, error = %d\n", error);

      for (n = 0; n < nexofiles; n++) {
        error = ex_put_var(exoidm[n], whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals);
        printf("after ex_put_nodal_var (%d), error = %d\n", n, error);
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

        for (n = 0; n < nexofiles; n++) {
          error = ex_put_var(exoidm[n], whole_time_step, EX_ELEM_BLOCK, k, ebids[j],
                             num_elem_in_block[j], elem_var_vals);
          printf("after ex_put_elem_var (%d), error = %d\n", n, error);
        }
      }
    }

    whole_time_step++;

    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    error = ex_update(exoid);
    printf("after ex_update, error = %d\n", error);
    for (n = 0; n < nexofiles; n++) {
      error = ex_update(exoidm[n]);
      printf("after ex_update (%d), error = %d\n", n, error);
    }
  }
  free(glob_var_vals);
  free(nodal_var_vals);
  free(elem_var_vals);

  /* close the EXODUS files
   */
  error = ex_close(exoid);
  printf("after ex_close, error = %d\n", error);
  for (n = 0; n < nexofiles; n++) {
    error = ex_close(exoidm[n]);
    printf("after ex_close (%d), error = %d\n", n, error);
  }
  return 0;
}
