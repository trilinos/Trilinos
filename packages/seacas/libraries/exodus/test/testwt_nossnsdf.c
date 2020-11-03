/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwt1 - test write an ExodusII database file
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
  int num_elem_in_block[10], num_nodes_per_elem[10], num_attr[10];
  int num_node_sets, num_side_sets, error;
  int i, j, k, m, *elem_map, *connect, *node_map;
  int node_list[100], elem_list[100], side_list[100];
  int ebids[10], ids[10];
  int num_nodes_per_set[10], num_elem_per_set[10];
  int num_df_per_set[10];
  int df_ind[10], node_ind[10], elem_ind[10];
  int num_qa_rec, num_info;
  int num_glo_vars, num_nod_vars, num_ele_vars;
  int whole_time_step, num_time_steps;
  int CPU_word_size, IO_word_size;
  int prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float  time_value;
  float  x[100], y[100], z[100];
  float  attrib[100], dist_fact[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char * prop_names[2];

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* float or double */
  IO_word_size  = 0; /* use system default (4 bytes) */

  /* create EXODUS II file */

  exoid = ex_create("test.exo",     /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */

  num_dim       = 3;
  num_nodes     = 28;
  num_elem      = 8;
  num_elem_blk  = 7;
  num_node_sets = 2;
  num_side_sets = 5;
  /* num_side_sets = 6; Uncomment to test NULL side sets */

  error = ex_put_init(exoid, "This is testwt1", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

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

  /* Circle #1 */
  x[20] = 100.0;
  y[20] = 100.0;
  z[20] = 0.0;

  /* Sphere #1 */
  x[21] = 50.0;
  y[21] = 50.0;
  z[21] = 20.0;

  /* Wedge #1 */
  x[22] = 3.0;
  y[22] = 0.0;
  z[22] = 6.0;
  x[23] = 6.0;
  y[23] = 0.0;
  z[23] = 0.0;
  x[24] = 0.0;
  y[24] = 0.0;
  z[24] = 0.0;
  x[25] = 3.0;
  y[25] = 2.0;
  z[25] = 6.0;
  x[26] = 6.0;
  y[26] = 2.0;
  z[26] = 2.0;
  x[27] = 0.0;
  y[27] = 2.0;
  z[27] = 0.0;

  error = ex_put_coord(exoid, x, y, z);
  printf("after ex_put_coord, error = %d\n", error);

  coord_names[0] = "xcoor";
  coord_names[1] = "ycoor";
  coord_names[2] = "zcoor";

  error = ex_put_coord_names(exoid, coord_names);
  printf("after ex_put_coord_names, error = %d\n", error);

  /* write element order map */

  elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i;
  }

  error = ex_put_map(exoid, elem_map);
  printf("after ex_put_map, error = %d\n", error);

  free(elem_map);

  /* write element numbering map */

  elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i * 2;
  }

  error = ex_put_id_map(exoid, EX_ELEM_MAP, elem_map);
  printf("after ex_put_elem_num_map, error = %d\n", error);

  free(elem_map);

  /* write node numbering map */

  node_map = (int *)calloc(num_nodes, sizeof(int));

  for (i = 1; i <= num_nodes; i++) {
    node_map[i - 1] = i * 3;
  }

  error = ex_put_id_map(exoid, EX_NODE_MAP, node_map);
  printf("after ex_put_node_num_map, error = %d\n", error);

  free(node_map);

  /* write element block parameters */

  num_elem_in_block[0] = 1; /* element 1: Quad 1 */
  num_elem_in_block[1] = 2; /* elements 2, 3: Quad 1 & 2 */
  num_elem_in_block[2] = 1; /* element 4: Hex    */
  num_elem_in_block[3] = 1; /* element 5: Tetra  */
  num_elem_in_block[4] = 1; /* element 6: Circle */
  num_elem_in_block[5] = 1; /* element 7: Sphere */
  num_elem_in_block[6] = 1; /* element 8: Wedge  */

  num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
  num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
  num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
  num_nodes_per_elem[3] = 4; /* elements in block #3 are 4-node tetras */
  num_nodes_per_elem[4] = 1; /* elements in block #4 are 1-node circles */
  num_nodes_per_elem[5] = 1; /* elements in block #5 are 1-node spheres */
  num_nodes_per_elem[6] = 6; /* elements in block #6 are 6-node wedges */

  ebids[0] = 10;
  ebids[1] = 11;
  ebids[2] = 12;
  ebids[3] = 13;
  ebids[4] = 14;
  ebids[5] = 15;
  ebids[6] = 16;

  num_attr[0] = 3;
  num_attr[1] = 3;
  num_attr[2] = 3;
  num_attr[3] = 3;
  num_attr[4] = 3;
  num_attr[5] = 3;
  num_attr[6] = 3;

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "quad", num_elem_in_block[0],
                       num_nodes_per_elem[0], 0, 0, num_attr[0]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[1], "quad", num_elem_in_block[1],
                       num_nodes_per_elem[1], 0, 0, num_attr[1]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[2], "hex", num_elem_in_block[2],
                       num_nodes_per_elem[2], 0, 0, num_attr[2]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[3], "tetra", num_elem_in_block[3],
                       num_nodes_per_elem[3], 0, 0, num_attr[3]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[4], "circle", num_elem_in_block[4],
                       num_nodes_per_elem[4], 0, 0, num_attr[4]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[5], "sphere", num_elem_in_block[5],
                       num_nodes_per_elem[5], 0, 0, num_attr[5]);
  printf("after ex_put_elem_block, error = %d\n", error);

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[6], "wedge", num_elem_in_block[6],
                       num_nodes_per_elem[6], 0, 0, num_attr[6]);
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
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[5], "MATL", 60);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[6], "MATL", 70);
  printf("after ex_put_prop, error = %d\n", error);

  /* write element connectivity */

  connect    = (int *)calloc(8, sizeof(int));
  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;
  connect[4] = 5;
  connect[5] = 6;
  connect[6] = 7;
  connect[7] = 8;

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

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[4], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 22;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[5], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  connect[0] = 23;
  connect[1] = 24;
  connect[2] = 25;
  connect[3] = 26;
  connect[4] = 27;
  connect[5] = 28;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[6], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  free(connect);

  /* write element block attributes  (3 per block) */

  attrib[0]  = 1.0;
  attrib[1]  = 2.0;
  attrib[2]  = 3.0;
  attrib[3]  = 1.11;
  attrib[4]  = 2.11;
  attrib[5]  = 3.11;
  attrib[6]  = 1.12;
  attrib[7]  = 2.12;
  attrib[8]  = 3.12;
  attrib[9]  = 1.2;
  attrib[10] = 2.2;
  attrib[11] = 3.2;
  attrib[12] = 1.3;
  attrib[13] = 2.3;
  attrib[14] = 3.3;
  attrib[15] = 1.4;
  attrib[16] = 2.4;
  attrib[17] = 3.4;
  attrib[18] = 1.5;
  attrib[19] = 2.5;
  attrib[20] = 3.5;
  attrib[21] = 1.6;
  attrib[22] = 2.6;
  attrib[23] = 3.6;

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], &attrib[0]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[1], &attrib[3]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[2], &attrib[9]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[3], &attrib[12]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[4], &attrib[15]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[5], &attrib[18]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[6], &attrib[21]);
  printf("after ex_put_elem_attr, error = %d\n", error);

  /* write individual node sets */

  ids[0] = 20;
  ids[1] = 21;

  num_nodes_per_set[0] = 5;
  num_nodes_per_set[1] = 3;
  /* num_nodes_per_set[1] = 0; Uncomment to test NULL node sets */

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

  {
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = ids;
    set_specs.num_entries_per_set = num_nodes_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = node_ind;
    set_specs.sets_dist_index     = df_ind;
    set_specs.sets_entry_list     = node_list;
    set_specs.sets_extra_list     = NULL;
    set_specs.sets_dist_fact      = dist_fact;
    error                         = ex_put_concat_sets(exoid, EX_NODE_SET, &set_specs);
  }

  printf("after ex_put_concat_node_sets, error = %d\n", error);

  error = ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4);
  printf("after ex_put_prop, error = %d\n", error);
  error = ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5);
  printf("after ex_put_prop, error = %d\n", error);

  prop_array[0] = 1000;
  prop_array[1] = 2000;

  error = ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array);
  printf("after ex_put_prop_array, error = %d\n", error);

  ids[0] = 30;
  ids[1] = 31;
  ids[2] = 32;
  ids[3] = 33;
  ids[4] = 34;
  ids[5] = 35;

  /* side set #1  - quad */
  node_list[0] = 8;
  node_list[1] = 5;
  node_list[2] = 6;
  node_list[3] = 7;

  /* side set #2  - quad/hex, spanning 2 element types  */
  node_list[4] = 2;
  node_list[5] = 3;
  node_list[6] = 7;
  node_list[7] = 8;

  /* side set #3  - hex */
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

  /* side set #4  - tetras */
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

  /* side set #5  - circle and sphere */
  node_list[48] = 21;
  node_list[49] = 22;

  /* side set #6  - wedges */
  node_list[50] = 27;
  node_list[51] = 26;
  node_list[52] = 23;
  node_list[53] = 24;

  node_list[54] = 28;
  node_list[55] = 27;
  node_list[56] = 24;
  node_list[57] = 25;

  node_list[58] = 28;
  node_list[59] = 25;
  node_list[60] = 23;
  node_list[61] = 26;

  node_list[62] = 25;
  node_list[63] = 24;
  node_list[64] = 23;

  node_list[65] = 26;
  node_list[66] = 27;
  node_list[67] = 28;

  node_ind[0] = 0;
  node_ind[1] = 4;
  node_ind[2] = 8;
  node_ind[3] = 36;
  node_ind[4] = 47;
  node_ind[5] = 49;

  num_elem_per_set[0] = 2; /* two sides uses 2 elements */
  num_elem_per_set[1] = 2;
  num_elem_per_set[2] = 7;
  num_elem_per_set[3] = 4;
  num_elem_per_set[4] = 2;
  num_elem_per_set[5] = 5;
  /* num_elem_per_set[5] = 0; Uncomment to test NULL side sets */

  num_nodes_per_set[0] = 4;
  num_nodes_per_set[1] = 4;
  num_nodes_per_set[2] = 28;
  num_nodes_per_set[3] = 12;
  num_nodes_per_set[4] = 2;
  num_nodes_per_set[5] = 18;

  elem_ind[0] = 0;
  elem_ind[1] = 2;
  elem_ind[2] = 4;
  elem_ind[3] = 11;
  elem_ind[4] = 15;
  elem_ind[5] = 17;

  elem_list[0]  = 3;
  elem_list[1]  = 3; /* side set 1: Quad #2 */
  elem_list[2]  = 1;
  elem_list[3]  = 3; /* side set 2: Quad #1 & #2 */
  elem_list[4]  = 4;
  elem_list[5]  = 4; /* side set 3: Hex */
  elem_list[6]  = 4;
  elem_list[7]  = 4;
  elem_list[8]  = 4;
  elem_list[9]  = 4;
  elem_list[10] = 4;
  elem_list[11] = 5;
  elem_list[12] = 5; /* side set 4: Tetra */
  elem_list[13] = 5;
  elem_list[14] = 5;
  elem_list[15] = 6;
  elem_list[16] = 7; /* side set 5: Circle & Sphere */
  elem_list[17] = 8;
  elem_list[18] = 8; /* side set 6: Wedge  */
  elem_list[19] = 8;
  elem_list[20] = 8;
  elem_list[21] = 8;

  error = ex_cvt_nodes_to_sides(exoid, num_elem_per_set, num_nodes_per_set, elem_ind, node_ind,
                                elem_list, node_list, side_list);
  printf("after ex_cvt_nodes_to_sides, error = %d\n", error);

  num_df_per_set[0] = 4;
  num_df_per_set[1] = 4;
  num_df_per_set[2] = 0;
  num_df_per_set[3] = 0;
  num_df_per_set[4] = 0;
  num_df_per_set[5] = 0;

  df_ind[0] = 0;
  df_ind[1] = 4;

  /* side set #1 df */
  dist_fact[0] = 30.0;
  dist_fact[1] = 30.1;
  dist_fact[2] = 30.2;
  dist_fact[3] = 30.3;

  /* side set #2 df */
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

  /* write QA records */

  num_qa_rec = 2;

  qa_record[0][0] = "TESTWT1";
  qa_record[0][1] = "testwt1";
  qa_record[0][2] = "03/16/94";
  qa_record[0][3] = "15:41:33";
  qa_record[1][0] = "FASTQ";
  qa_record[1][1] = "fastq";
  qa_record[1][2] = "07/07/93";
  qa_record[1][3] = "16:41:33";

  error = ex_put_qa(exoid, num_qa_rec, qa_record);
  printf("after ex_put_qa, error = %d\n", error);

  /* write information records */

  num_info = 3;

  info[0] = "This is the first information record.";
  info[1] = "This is the second information record.";
  info[2] = "This is the third information record.";

  error = ex_put_info(exoid, num_info, info);
  printf("after ex_put_info, error = %d\n", error);

  /* write results variables parameters and names */

  num_glo_vars = 1;

  var_names[0] = "glo vars";

  error = ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_name(exoid, EX_GLOBAL, 1, var_names[0]);
  printf("after ex_put_variable_name, error = %d\n", error);

  num_nod_vars = 2;

  var_names[0] = "nod_var0";
  var_names[1] = "nod_var1";

  error = ex_put_variable_param(exoid, EX_NODAL, num_nod_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);

  num_ele_vars = 3;

  var_names[0] = "ele_var0";
  var_names[1] = "ele_var1";
  var_names[2] = "ele_var2";

  error = ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  error = ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);

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

    /* write global variables */

    for (j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (float)(j + 2) * time_value;
    }

    error = ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals);
    printf("after ex_put_glob_vars, error = %d\n", error);

    /* write nodal variables */

    for (k = 1; k <= num_nod_vars; k++) {
      for (j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (float)k + ((float)(j + 1) * time_value);
      }

      error = ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals);
      printf("after ex_put_nodal_var, error = %d\n", error);
    }

    /* write element variables */

    for (k = 1; k <= num_ele_vars; k++) {
      for (j = 0; j < num_elem_blk; j++) {
        for (m = 0; m < num_elem_in_block[j]; m++) {
          elem_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          /* printf("elem_var_vals[%d]: %f\n",m,elem_var_vals[m]); */
        }
        if (k == 1 && j == 2) {
          continue; /* skip element block 3, variable 1 */
        }

        error = ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, ebids[j], num_elem_in_block[j],
                           elem_var_vals);
        printf("after ex_put_elem_var, error = %d\n", error);
      }
    }

    whole_time_step++;

    /* update the data file; this should be done at the end of every time step
     * to ensure that no data is lost if the analysis dies
     */
    error = ex_update(exoid);
    printf("after ex_update, error = %d\n", error);
  }
  free(glob_var_vals);
  free(nodal_var_vals);
  free(elem_var_vals);

  /* close the EXODUS files
   */
  error = ex_close(exoid);
  printf("after ex_close, error = %d\n", error);
  return 0;
}
