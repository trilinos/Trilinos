/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testwt_ss - test write an ExodusII database file
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
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_elem_in_block[10], num_nodes_per_elem[10];
  int  num_node_sets, num_side_sets, error;
  int *connect;
  int  node_list[100], elem_list[100], side_list[100];
  int  ebids[10], ids[10];
  int  num_nodes_per_set[10], num_elem_per_set[10];
  int  num_df_per_set[10];
  int  df_ind[10], node_ind[10], elem_ind[10];
  int  num_qa_rec, num_info;
  int  CPU_word_size, IO_word_size;

  float x[100], y[100], z[100];
  float attrib[1], dist_fact[100];
  char *coord_names[3], *qa_record[2][4], *info[3];

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
  num_elem      = 8;
  num_elem_blk  = 8;
  num_node_sets = 2;
  num_side_sets = 9;

  error = ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  x[0] = 0.0f;
  y[0] = 0.0f;
  z[0] = 0.0f;
  x[1] = 1.0f;
  y[1] = 0.0f;
  z[1] = 0.0f;
  x[2] = 1.0f;
  y[2] = 1.0f;
  z[2] = 0.0f;
  x[3] = 0.0f;
  y[3] = 1.0f;
  z[3] = 0.0f;

  /* Quad #2 */
  x[4] = 1.0f;
  y[4] = 0.0f;
  z[4] = 0.0f;
  x[5] = 2.0f;
  y[5] = 0.0f;
  z[5] = 0.0f;
  x[6] = 2.0f;
  y[6] = 1.0f;
  z[6] = 0.0f;
  x[7] = 1.0f;
  y[7] = 1.0f;
  z[7] = 0.0f;

  /* Hex #1 */
  x[8]  = 0.0f;
  y[8]  = 0.0f;
  z[8]  = 0.0f;
  x[9]  = 10.0f;
  y[9]  = 0.0f;
  z[9]  = 0.0f;
  x[10] = 10.0f;
  y[10] = 0.0f;
  z[10] = -10.0f;
  x[11] = 1.0f;
  y[11] = 0.0f;
  z[11] = -10.0f;
  x[12] = 1.0f;
  y[12] = 10.0f;
  z[12] = 0.0f;
  x[13] = 10.0f;
  y[13] = 10.0f;
  z[13] = 0.0f;
  x[14] = 10.0f;
  y[14] = 10.0f;
  z[14] = -10.0f;
  x[15] = 1.0f;
  y[15] = 10.0f;
  z[15] = -10.0f;

  /* Tetra #1 */
  x[16] = 0.0f;
  y[16] = 0.0f;
  z[16] = 0.0f;
  x[17] = 1.0f;
  y[17] = 0.0f;
  z[17] = 5.0f;
  x[18] = 10.0f;
  y[18] = 0.0f;
  z[18] = 2.0f;
  x[19] = 7.0f;
  y[19] = 5.0f;
  z[19] = 3.0f;

  /* Wedge #1 */
  x[20] = 3.0f;
  y[20] = 0.0f;
  z[20] = 6.0f;
  x[21] = 6.0f;
  y[21] = 0.0f;
  z[21] = 0.0f;
  x[22] = 0.0f;
  y[22] = 0.0f;
  z[22] = 0.0f;
  x[23] = 3.0f;
  y[23] = 2.0f;
  z[23] = 6.0f;
  x[24] = 6.0f;
  y[24] = 2.0f;
  z[24] = 2.0f;
  x[25] = 0.0f;
  y[25] = 2.0f;
  z[25] = 0.0f;

  /* Tetra #2 */
  x[26] = 2.7f;
  y[26] = 1.7f;
  z[26] = 2.7f;
  x[27] = 6.0f;
  y[27] = 1.7f;
  z[27] = 3.3f;
  x[28] = 5.7f;
  y[28] = 1.7f;
  z[28] = 1.7f;
  x[29] = 3.7f;
  y[29] = 0.0f;
  z[29] = 2.3f;

  /* TriShell #1 */
  x[30] = 2.7f;
  y[30] = 1.7f;
  z[30] = 2.7f;
  x[31] = 6.0f;
  y[31] = 1.7f;
  z[31] = 3.3f;
  x[32] = 5.7f;
  y[32] = 1.7f;
  z[32] = 1.7f;

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

  /* write element block parameters */

  num_elem_in_block[0] = 1;
  num_elem_in_block[1] = 1;
  num_elem_in_block[2] = 1;
  num_elem_in_block[3] = 1;
  num_elem_in_block[4] = 1;
  num_elem_in_block[5] = 1;
  num_elem_in_block[6] = 1;
  num_elem_in_block[7] = 1;

  num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
  num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
  num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
  num_nodes_per_elem[3] = 4; /* elements in block #4 are 4-node tetras */
  num_nodes_per_elem[4] = 6; /* elements in block #5 are 6-node wedges */
  num_nodes_per_elem[5] = 8; /* elements in block #6 are 8-node tetras */
  num_nodes_per_elem[6] = 4; /* elements in block #7 are 4-node shells */
  num_nodes_per_elem[7] = 3; /* elements in block #8 are 3-node shells */

  ebids[0] = 10;
  ebids[1] = 11;
  ebids[2] = 12;
  ebids[3] = 13;
  ebids[4] = 14;
  ebids[5] = 15;
  ebids[6] = 16;
  ebids[7] = 17;

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "quad", num_elem_in_block[0],
                       num_nodes_per_elem[0], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[1], "quad", num_elem_in_block[1],
                       num_nodes_per_elem[1], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[2], "hex", num_elem_in_block[2],
                       num_nodes_per_elem[2], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[3], "tetra", num_elem_in_block[3],
                       num_nodes_per_elem[3], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[4], "wedge", num_elem_in_block[4],
                       num_nodes_per_elem[4], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[5], "tetra", num_elem_in_block[5],
                       num_nodes_per_elem[5], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[6], "shell", num_elem_in_block[6],
                       num_nodes_per_elem[6], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[7], "triangle", num_elem_in_block[7],
                       num_nodes_per_elem[7], 0, 0, 1);
  printf("after ex_put_elem_block, error = %d\n", error);

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

  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[6], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 30;
  connect[1] = 31;
  connect[2] = 32;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[7], connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(connect);

  /* write element block attributes */
  attrib[0] = 3.14159f;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[0], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  attrib[0] = 6.14159f;
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

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[7], attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write individual side sets */

  /* write concatenated side sets; this produces the same information as
   * the above code which writes individual side sets
   */

  ids[0] = 30;
  ids[1] = 31;
  ids[2] = 32;
  ids[3] = 33;
  ids[4] = 34;
  ids[5] = 35;
  ids[6] = 36;
  ids[7] = 37;
  ids[8] = 38;

  /* side set #1  - NULL side set */
  /* do nothing except set num_elem_per_set to 0 */

  /* side set #2  - NULL side set */
  /* do nothing except set num_elem_per_set to 0 */

  /* side set #3  - quad; 2 sides */

  node_list[0] = 8;
  node_list[1] = 5;
  elem_list[0] = 2;

  node_list[2] = 6;
  node_list[3] = 7;
  elem_list[1] = 2;

  /* side set #4  - quad; 2 sides spanning 2 elements  */

  node_list[4] = 2;
  node_list[5] = 3;
  elem_list[2] = 1;

  node_list[6] = 7;
  node_list[7] = 8;
  elem_list[3] = 2;

  /* side set #5  - hex; 7 sides */

  node_list[8]  = 9;
  node_list[9]  = 12;
  node_list[10] = 11;
  node_list[11] = 10;
  elem_list[4]  = 3;

  node_list[12] = 11;
  node_list[13] = 12;
  node_list[14] = 16;
  node_list[15] = 15;
  elem_list[5]  = 3;

  node_list[16] = 16;
  node_list[17] = 15;
  node_list[18] = 11;
  node_list[19] = 12;
  elem_list[6]  = 3;

  node_list[20] = 10;
  node_list[21] = 11;
  node_list[22] = 15;
  node_list[23] = 14;
  elem_list[7]  = 3;

  node_list[24] = 13;
  node_list[25] = 16;
  node_list[26] = 12;
  node_list[27] = 9;
  elem_list[8]  = 3;

  node_list[28] = 14;
  node_list[29] = 13;
  node_list[30] = 9;
  node_list[31] = 10;
  elem_list[9]  = 3;

  node_list[32] = 16;
  node_list[33] = 13;
  node_list[34] = 14;
  node_list[35] = 15;
  elem_list[10] = 3;

  /* side set #6  - 4-node tetras; 4 sides */

  node_list[36] = 17;
  node_list[37] = 18;
  node_list[38] = 20;
  elem_list[11] = 4;

  node_list[39] = 18;
  node_list[40] = 19;
  node_list[41] = 20;
  elem_list[12] = 4;

  node_list[42] = 17;
  node_list[43] = 20;
  node_list[44] = 19;
  elem_list[13] = 4;

  node_list[45] = 17;
  node_list[46] = 19;
  node_list[47] = 18;
  elem_list[14] = 4;

  /* side set #7  - shells; front and back faces */

  node_list[48] = 1;
  node_list[49] = 2;
  node_list[50] = 3;
  node_list[51] = 4;
  elem_list[15] = 7;

  node_list[52] = 4;
  node_list[53] = 3;
  node_list[54] = 2;
  node_list[55] = 1;
  elem_list[16] = 7;

  /* side set #8  - shells; 4 edges */

  node_list[56] = 1;
  node_list[57] = 2;
  elem_list[17] = 7;

  node_list[58] = 2;
  node_list[59] = 3;
  elem_list[18] = 7;

  node_list[60] = 3;
  node_list[61] = 4;
  elem_list[19] = 7;

  node_list[62] = 4;
  node_list[63] = 1;
  elem_list[20] = 7;

  /* side set #9 --  3-node shells -- front and back */

  node_list[64] = 30;
  node_list[65] = 31;
  node_list[66] = 32;
  elem_list[21] = 8;

  node_list[67] = 32;
  node_list[68] = 31;
  node_list[69] = 30;
  elem_list[22] = 8;

  /* set up indices */
  node_ind[0] = 0;
  node_ind[1] = 0;
  node_ind[2] = 0;
  node_ind[3] = 4;
  node_ind[4] = 8;
  node_ind[5] = 36;
  node_ind[6] = 48;
  node_ind[7] = 56;
  node_ind[8] = 64;

  num_elem_per_set[0] = 0;
  num_elem_per_set[1] = 0;
  num_elem_per_set[2] = 2;
  num_elem_per_set[3] = 2;
  num_elem_per_set[4] = 7;
  num_elem_per_set[5] = 4;
  num_elem_per_set[6] = 2;
  num_elem_per_set[7] = 4;
  num_elem_per_set[8] = 2;

  num_nodes_per_set[0] = 0;
  num_nodes_per_set[1] = 0;
  num_nodes_per_set[2] = 4;
  num_nodes_per_set[3] = 4;
  num_nodes_per_set[4] = 28;
  num_nodes_per_set[5] = 12;
  num_nodes_per_set[6] = 8;
  num_nodes_per_set[7] = 8;
  num_nodes_per_set[8] = 6;

  elem_ind[0] = 0;
  elem_ind[1] = 0;
  elem_ind[2] = 0;
  elem_ind[3] = 2;
  elem_ind[4] = 4;
  elem_ind[5] = 11;
  elem_ind[6] = 15;
  elem_ind[7] = 17;
  elem_ind[8] = 21;

  error = ex_cvt_nodes_to_sides(exoid, num_elem_per_set, num_nodes_per_set, elem_ind, node_ind,
                                elem_list, node_list, side_list);
  printf("after ex_cvt_nodes_to_sides, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  num_df_per_set[0] = 0;
  num_df_per_set[1] = 0;
  num_df_per_set[2] = 4;
  num_df_per_set[3] = 4;
  num_df_per_set[4] = 0;
  num_df_per_set[5] = 0;
  num_df_per_set[6] = 0;
  num_df_per_set[7] = 0;
  num_df_per_set[8] = 0;

  df_ind[0] = 0;
  df_ind[1] = 0;
  df_ind[2] = 0;
  df_ind[3] = 4;
  df_ind[4] = 0;
  df_ind[5] = 0;
  df_ind[6] = 0;
  df_ind[7] = 0;
  df_ind[8] = 0;

  dist_fact[0] = 30.0f;
  dist_fact[1] = 30.1f;
  dist_fact[2] = 30.2f;
  dist_fact[3] = 30.3f;

  dist_fact[4] = 31.0f;
  dist_fact[5] = 31.1f;
  dist_fact[6] = 31.2f;
  dist_fact[7] = 31.3f;

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

    error = ex_put_concat_sets(exoid, EX_SIDE_SET, &set_specs);
  }
  printf("after ex_put_concat_side_sets, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* THIS SECTION IS COMMENTED OUT
      END COMMENTED OUT SECTION */

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
