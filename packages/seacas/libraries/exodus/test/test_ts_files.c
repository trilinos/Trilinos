/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <exodusII.h>
#include <exodusII_int.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Create a separate exodus file on each thread (same file as written
 * by testwt.c)
 */

#define NUM_THREADS 8

typedef struct
{
  long threadid;
  int  exoid;
} param;

void *init_file(void *varg)
{
  param *arg = (param *)varg;

  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_face_in_sset[10], num_nodes_in_nset[10];
  int  num_node_sets, num_side_sets, error;
  int  i, j, k, m, *elem_map, *connect;
  int  node_list[100], elem_list[100], side_list[100];
  int  ssids[10], nsids[10];
  int  num_qa_rec, num_info;
  int  num_glo_vars, num_nod_vars, num_ele_vars, num_sset_vars, num_nset_vars;
  int *truth_tab;
  int  whole_time_step, num_time_steps;
  int  CPU_word_size, IO_word_size;
  int  prop_array[2];

  float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
  float *sset_var_vals, *nset_var_vals;
  float  time_value;
  float  x[100], y[100], z[100];
  float  attrib[1], dist_fact[100];
  char * coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char * block_names[10], *nset_names[10], *sset_names[10];
  char * prop_names[2], *attrib_names[2];
  char * title = "This is a test";

  struct ex_block blocks[10];

  char name[32];

  /* Name is "test_{thread}.exo" */
  sprintf(name, "test%ld.exo", arg->threadid);

  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */

  exoid = ex_create(name,           /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for %s, exoid = %d\n", name, exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */

  num_dim       = 3;
  num_nodes     = 33;
  num_elem      = 7;
  num_elem_blk  = 7;
  num_node_sets = 2;
  num_side_sets = 5;

  error = ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
                      num_side_sets);

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

  /* Add nodal attributes */
  error = ex_put_attr_param(exoid, EX_NODAL, 0, 2);
  printf("after ex_put_attr_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_one_attr(exoid, EX_NODAL, 0, 1, x);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_one_attr(exoid, EX_NODAL, 0, 2, y);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  {
    attrib_names[0] = "Node_attr_1";
    attrib_names[1] = "Node_attr_2";
    error           = ex_put_attr_names(exoid, EX_NODAL, 0, attrib_names);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
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
  for (i = 0; i < 10; i++) {
    blocks[i].type                = EX_ELEM_BLOCK;
    blocks[i].id                  = 0;
    blocks[i].num_entry           = 0;
    blocks[i].num_nodes_per_entry = 0;
    blocks[i].num_edges_per_entry = 0;
    blocks[i].num_faces_per_entry = 0;
    blocks[i].num_attribute       = 0;
  }

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

  error = ex_put_block_params(exoid, num_elem_blk, blocks);
  printf("after ex_put_block_params, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* Write element block names */
  for (i = 0; i < num_elem_blk; i++) {
    error = ex_put_name(exoid, EX_ELEM_BLOCK, blocks[i].id, block_names[i]);
    printf("after ex_put_names, error = %d\n", error);
  }

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* write element block properties */

  /*               12345678901234567890123456789012 */
  prop_names[0] = "MATERIAL_PROPERTY_LONG_NAME_32CH";
  prop_names[1] = "DENSITY";
  error         = ex_put_prop_names(exoid, EX_ELEM_BLOCK, 2, prop_names);
  printf("after ex_put_prop_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[0].id, prop_names[0], 10);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[1].id, prop_names[0], 20);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[2].id, prop_names[0], 30);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[3].id, prop_names[0], 40);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[4].id, prop_names[0], 50);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[5].id, prop_names[0], 60);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_ELEM_BLOCK, blocks[6].id, prop_names[0], 70);
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

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[0].id, connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 5;
  connect[1] = 6;
  connect[2] = 7;
  connect[3] = 8;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[1].id, connect, NULL, NULL);
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

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[2].id, connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 17;
  connect[1] = 18;
  connect[2] = 19;
  connect[3] = 20;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[3].id, connect, NULL, NULL);
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

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[4].id, connect, NULL, NULL);
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

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[5].id, connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  connect[0] = 31;
  connect[1] = 32;
  connect[2] = 33;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[6].id, connect, NULL, NULL);
  printf("after ex_put_elem_conn, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  free(connect);

  /* write element block attributes */
  attrib[0] = 3.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[0].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[0].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  attrib[0] = 6.14159;
  error     = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[1].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[2].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[3].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[4].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[5].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_ELEM_BLOCK, blocks[6].id, attrib);
  printf("after ex_put_elem_attr, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  attrib_names[0] = "THICKNESS";
  for (i = 0; i < num_elem_blk; i++) {
    error = ex_put_attr_names(exoid, EX_ELEM_BLOCK, blocks[i].id, attrib_names);
    printf("after ex_put_elem_attr_names, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
  }

  /* write individual node sets */

  num_nodes_in_nset[0] = 5;
  num_nodes_in_nset[1] = 3;

  nsids[0] = 20;
  nsids[1] = 21;

  error = ex_put_set_param(exoid, EX_NODE_SET, nsids[0], 5, 5);
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

  error = ex_put_set(exoid, EX_NODE_SET, nsids[0], node_list, NULL);
  printf("after ex_put_node_set, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[0], dist_fact);
  printf("after ex_put_node_set_dist_fact, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_param(exoid, EX_NODE_SET, nsids[1], 3, 3);
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

  error = ex_put_set(exoid, EX_NODE_SET, nsids[1], node_list, NULL);
  printf("after ex_put_node_set, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_set_dist_fact(exoid, EX_NODE_SET, nsids[1], dist_fact);
  printf("after ex_put_node_set_dist_fact, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  /* Write node set names */
  nset_names[0] = "nset_1";
  nset_names[1] = "nset_2";

  error = ex_put_names(exoid, EX_NODE_SET, nset_names);
  printf("after ex_put_names, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_NODE_SET, nsids[0], "FACE", 4);
  printf("after ex_put_prop, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_prop(exoid, EX_NODE_SET, nsids[1], "FACE", 5);
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

  /* Add nodeset attributes */
  error = ex_put_attr_param(exoid, EX_NODE_SET, nsids[0], 1);
  printf("after ex_put_attr_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  error = ex_put_attr(exoid, EX_NODE_SET, nsids[0], x);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  {
    attrib_names[0] = "Nodeset_attribute";
    error           = ex_put_attr_names(exoid, EX_NODE_SET, nsids[0], attrib_names);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
  }

  /* write individual side sets */
  num_face_in_sset[0] = 2;
  num_face_in_sset[1] = 2;
  num_face_in_sset[2] = 7;
  num_face_in_sset[3] = 8;
  num_face_in_sset[4] = 10;

  ssids[0] = 30;
  ssids[1] = 31;
  ssids[2] = 32;
  ssids[3] = 33;
  ssids[4] = 34;

  /* side set #1  - quad */

  error = ex_put_set_param(exoid, EX_SIDE_SET, ssids[0], 2, 4);
  printf("after ex_put_side_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set_param(exoid, EX_SIDE_SET, 31, 2, 4);
  printf("after ex_put_side_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set_param(exoid, EX_SIDE_SET, 32, 7, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set_param(exoid, EX_SIDE_SET, 33, 8, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  error = ex_put_set_param(exoid, EX_SIDE_SET, 34, 10, 0);
  printf("after ex_put_side_set_param, error = %d\n", error);

  if (error) {
    ex_close(exoid);
    exit(-1);
  }

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

  /* Write side set names */
  sset_names[0] = "sset_1";
  sset_names[1] = "sset_2";
  sset_names[2] = "sset_3";
  sset_names[3] = "sset_4";
  sset_names[4] = "sset_5";

  error = ex_put_names(exoid, EX_SIDE_SET, sset_names);
  printf("after ex_put_names, error = %d\n", error);

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

  /* write QA records; test empty and just blank-filled records */

  num_qa_rec = 2;

  qa_record[0][0] = "TESTWT";
  qa_record[0][1] = "testwt";
  qa_record[0][2] = "07/07/93";
  qa_record[0][3] = "15:41:33";
  qa_record[1][0] = "Thirty-Two character QA Record|";
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
  info[2] = "This info record is exactly 80 characters long.  last character should be pipe |";

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

  {
    num_nset_vars = 3;

    var_names[0] = "ns_var0";
    var_names[1] = "ns_var1";
    var_names[2] = "ns_var2";

    error = ex_put_variable_param(exoid, EX_NODE_SET, num_nset_vars);
    printf("after ex_put_variable_param, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    error = ex_put_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names);
    printf("after ex_put_variable_names, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
  }

  {
    num_sset_vars = 3;

    var_names[0] = "ss_var0";
    var_names[1] = "ss_var1";
    var_names[2] = "ss_var2";

    error = ex_put_variable_param(exoid, EX_SIDE_SET, num_sset_vars);
    printf("after ex_put_variable_param, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    error = ex_put_variable_names(exoid, EX_SIDE_SET, num_sset_vars, var_names);
    printf("after ex_put_variable_names, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
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
  sset_var_vals  = (float *)calloc(10, CPU_word_size);
  nset_var_vals  = (float *)calloc(10, CPU_word_size);

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

    /* write sideset variables */

    for (k = 1; k <= num_sset_vars; k++) {
      for (j = 0; j < num_side_sets; j++) {
        for (m = 0; m < num_face_in_sset[j]; m++) {
          sset_var_vals[m] = (float)(k + 2) + (float)(j + 3) + ((float)(m + 1) * time_value);
          /* printf("sset_var_vals[%d]: %f\n",m,sset_var_vals[m]); */
        }
        error = ex_put_var(exoid, whole_time_step, EX_SIDE_SET, k, ssids[j], num_face_in_sset[j],
                           sset_var_vals);
        printf("after ex_put_sset_var, error = %d\n", error);
        if (error) {
          ex_close(exoid);
          exit(-1);
        }
      }
    }

    /* write nodeset variables */

    for (k = 1; k <= num_nset_vars; k++) {
      for (j = 0; j < num_node_sets; j++) {
        for (m = 0; m < num_nodes_in_nset[j]; m++) {
          nset_var_vals[m] = (float)(k + 3) + (float)(j + 4) + ((float)(m + 1) * time_value);
          /* printf("nset_var_vals[%d]: %f\n",m,nset_var_vals[m]); */
        }
        error = ex_put_var(exoid, whole_time_step, EX_NODE_SET, k, nsids[j], num_nodes_in_nset[j],
                           nset_var_vals);
        printf("after ex_put_nset_var, error = %d\n", error);
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
  free(sset_var_vals);
  free(nset_var_vals);

  /* close the EXODUS files
   */
  error = ex_close(exoid);
  printf("after ex_close, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }
  arg->exoid = exoid;
  return arg;
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int       rc;
  long      t;

  param arg[NUM_THREADS];

  printf("Running on %d threads\n", NUM_THREADS);
  for (t = 0; t < NUM_THREADS; t++) {
    arg[t].threadid = t;
    rc              = pthread_create(&threads[t], NULL, init_file, (void *)(arg + t));
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    pthread_join(threads[t], NULL);
  }
}
