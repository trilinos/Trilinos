/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_elem_in_block[10], num_nodes_per_elem[10];
  int  num_node_sets, num_side_sets, error;
  int  i, j, k, kk, m, *elem_map, *connect;
  int  node_list[100], elem_list[100], side_list[100];
  int  ebids[10], ssids[10], nsids[10], nattr[10];
  int  num_nodes_per_set[10], num_elem_per_set[10];
  int  num_df_per_set[10];
  int  num_qa_rec, num_info;
  int  num_glo_vars, num_nod_vars, num_ele_vars, num_nset_vars, num_sset_vars;
  int *truth_tab, *nset_tab, *sset_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;
  int  prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals, *nset_var_vals, *sset_var_vals;
  float  time_value;
  float  x[100], y[100], z[100];
  float  attrib[1], dist_fact[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[7];
  char * prop_names[2];
  char * eb_type[10];

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
  num_elem      = 7;
  num_elem_blk  = 7;
  num_node_sets = 2;
  num_side_sets = 5;

  error = ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

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

  /* write element order map */

  elem_map = (int *)calloc(num_elem, sizeof(int));

  for (i = 1; i <= num_elem; i++) {
    elem_map[i - 1] = i;
  }

  error = ex_put_map(exoid, elem_map);
  printf("after ex_put_map, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(elem_map);

  /* write element block parameters */

  num_elem_in_block[0] = 1;
  num_elem_in_block[1] = 1;
  num_elem_in_block[2] = 1;
  num_elem_in_block[3] = 1;
  num_elem_in_block[4] = 1;
  num_elem_in_block[5] = 1;
  num_elem_in_block[6] = 1;

  num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
  num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
  num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
  num_nodes_per_elem[3] = 4; /* elements in block #4 are 4-node tetras */
  num_nodes_per_elem[4] = 6; /* elements in block #5 are 6-node wedges */
  num_nodes_per_elem[5] = 8; /* elements in block #6 are 8-node tetras */
  num_nodes_per_elem[6] = 3; /* elements in block #7 are 3-node tris   */

  ebids[0] = 10;
  ebids[1] = 11;
  ebids[2] = 12;
  ebids[3] = 13;
  ebids[4] = 14;
  ebids[5] = 15;
  ebids[6] = 16;

  nattr[0] = nattr[1] = nattr[2] = nattr[3] = 1;
  nattr[4] = nattr[5] = nattr[6] = 1;

  eb_type[0] = "quad";
  eb_type[1] = "quad";
  eb_type[2] = "hex";
  eb_type[3] = "tetra";
  eb_type[4] = "wedge";
  eb_type[5] = "tetra";
  eb_type[6] = "tri";

  error = ex_put_concat_elem_block(exoid, ebids, eb_type, num_elem_in_block, num_nodes_per_elem,
                                   nattr, 0);
  printf("after ex_put_concat_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element block properties */

  prop_names[0] = "MATL";
  prop_names[1] = "DENSITY";
  error         = ex_put_prop_names(exoid, EX_ELEM_BLOCK, 2, prop_names);
  printf("after ex_put_prop_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], "MATL", 10);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], "MATL", 20);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], "MATL", 30);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[3], "MATL", 40);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[4], "MATL", 50);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[5], "MATL", 60);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[6], "MATL", 70);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element connectivity */

  connect    = (int *)calloc(8, sizeof(int));
  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 5;
  connect[1] = 6;
  connect[2] = 7;
  connect[3] = 8;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[1], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[3], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 21;
  connect[1] = 22;
  connect[2] = 23;
  connect[3] = 24;
  connect[4] = 25;
  connect[5] = 26;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[4], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;
  connect[4] = 27;
  connect[5] = 28;
  connect[6] = 30;
  connect[7] = 29;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[5], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 31;
  connect[1] = 32;
  connect[2] = 33;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[6], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(connect);

  /* write element block attributes */

  attrib[0] = 3.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  attrib[0] = 6.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[1], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[2], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[3], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[4], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[5], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[6], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write individual node sets */

  nsids[0]             = 20;
  nsids[1]             = 21;
  num_nodes_per_set[0] = 5;
  num_nodes_per_set[1] = 3;
  num_df_per_set[0]    = 5;
  num_df_per_set[1]    = 3;

  {
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = nsids;
    set_specs.num_entries_per_set = num_nodes_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = NULL;
    set_specs.sets_dist_index     = NULL;
    set_specs.sets_entry_list     = NULL;
    set_specs.sets_extra_list     = NULL;
    set_specs.sets_dist_fact      = NULL;

    error = ex_put_concat_sets(exoid, EX_NODE_SET, &set_specs);
  }
  printf("after ex_put_concat_node_sets, error = %d\n", error);

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

  /* Define the sideset params at one time, then write individually */
  ssids[0] = 30;
  ssids[1] = 31;
  ssids[2] = 32;
  ssids[3] = 33;
  ssids[4] = 34;

  num_elem_per_set[0] = 2;
  num_elem_per_set[1] = 2;
  num_elem_per_set[2] = 7;
  num_elem_per_set[3] = 8;
  num_elem_per_set[4] = 10;

  num_df_per_set[0] = 4;
  num_df_per_set[1] = 4;
  num_df_per_set[2] = 0;
  num_df_per_set[3] = 0;
  num_df_per_set[4] = 0;

  {
    struct ex_set_specs set_specs;

    set_specs.sets_ids            = ssids;
    set_specs.num_entries_per_set = num_elem_per_set;
    set_specs.num_dist_per_set    = num_df_per_set;
    set_specs.sets_entry_index    = NULL;
    set_specs.sets_dist_index     = NULL;
    set_specs.sets_entry_list     = NULL;
    set_specs.sets_extra_list     = NULL;
    set_specs.sets_dist_fact      = NULL;
    error                         = ex_put_concat_sets(exoid, EX_SIDE_SET, &set_specs);
  }
  printf("after ex_put_concat_side_sets, error = %d\n", error);

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

  error = ex_put_set(exoid, EX_SIDE_SET, 30, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_SIDE_SET, 30, dist_fact);
  printf("after ex_put_side_set_dist_fact, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* side set #2  - quad, spanning 2 elements  */

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

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_SIDE_SET, 31, dist_fact);
  printf("after ex_put_side_set_dist_fact, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set(exoid, EX_SIDE_SET, 32, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set(exoid, EX_SIDE_SET, 33, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set(exoid, EX_SIDE_SET, 34, elem_list, side_list);
  printf("after ex_put_side_set, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100);
  printf("after ex_put_prop, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101);
  printf("after ex_put_prop, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  ex_put_all_var_param(exoid, num_glo_vars, num_nod_vars, num_ele_vars, truth_tab, num_nset_vars,
                       nset_tab, num_sset_vars, sset_tab);
  printf("after ex_put_all_var_param, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(truth_tab);
  free(nset_tab);
  free(sset_tab);

  var_names[0] = "glo_vars";
  error        = ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /*              12345678901234567890123456789012 */
  var_names[0] = "node_variable_a_very_long_name_0";
  var_names[1] = "nod_var1";
  error        = ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  var_names[0] = "ele_var0";
  var_names[1] = "ele_var1";
  var_names[2] = "ele_var2";
  error        = ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  var_names[0] = "nset_var0";
  var_names[1] = "nset_var1";
  var_names[2] = "nset_var2";
  var_names[3] = "nset_var3";
  error        = ex_put_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  var_names[0] = "sset_var0";
  var_names[1] = "sset_var1";
  var_names[2] = "sset_var2";
  var_names[3] = "sset_var3";
  var_names[4] = "sset_var4";
  var_names[5] = "sset_var5";
  var_names[6] = "sset_var6";
  error        = ex_put_variable_names(exoid, EX_SIDE_SET, num_sset_vars, var_names);
  printf("after ex_put_variable_names, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

    error = ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 0, num_glo_vars, glob_var_vals);
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

    /* write nodeset variables */

    kk = 0;
    for (j = 0; j < num_node_sets; j++) {
      for (k = 0; k < num_nset_vars; k++) {
        if (kk++ % 2 == 0) {
          for (m = 0; m < num_nodes_per_set[j]; m++) {
            nset_var_vals[m] = (float)(k + 1) + (float)(j + 2) + ((float)(m + 1) * time_value);
          }
          error = ex_put_var(exoid, whole_time_step, EX_NODE_SET, k + 1, nsids[j],
                             num_nodes_per_set[j], nset_var_vals);
          printf("after ex_put_nset_var, error = %d\n", error);
          if (error) {
            ex_close(exoid);
            exit(-1);
          }
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
          error = ex_put_var(exoid, whole_time_step, EX_SIDE_SET, k + 1, ssids[j],
                             num_elem_per_set[j], sset_var_vals);
          printf("after ex_put_sset_var, error = %d\n", error);
          if (error) {
            ex_close(exoid);
            exit(-1);
          }
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
  free(nset_var_vals);
  free(sset_var_vals);

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
