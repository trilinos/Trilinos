/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testrdwt - test reading from one ExodusII file and writing to another
 *            ExodusII file open concurrently
 *
 * author - Sandia National Laboratories
 *          Larry A. Schoof - Original
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
 *  database read and write routines. It tests reading from an open EXODUSII
 *  file and writing to another concurrently opened EXODUSII file.
 *
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
  int num_elem_in_block, num_nodes_per_elem, num_attr;
  int num_nodes_in_set;
  int num_sides_in_set, num_df_in_set;

  char title[MAX_LINE_LENGTH + 1], elem_type[MAX_STR_LENGTH + 1];

  /* Specify compute and i/o word size */

  int CPU_word_size = 0; /* sizeof(float) */
  int IO_word_size  = 4; /* float */

  /* open EXODUS II file for reading */

  ex_opts(EX_VERBOSE | EX_ABORT);

  float version;
  int   exoid = ex_open("test.exo",     /* filename path */
                        EX_READ,        /* access mode */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size,  /* I/O float word size in bytes */
                        &version);      /* returned version number */
  printf("after ex_open for test.exo\n");
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* create EXODUS II file for writing */

  int exoid2 = ex_create("test2.exo",    /* filename path */
                         EX_CLOBBER,     /* create mode */
                         &CPU_word_size, /* CPU float word size in bytes */
                         &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test2.exo, exoid = %d\n", exoid2);

  /* read initialization parameters */

  int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk,
                          &num_node_sets, &num_side_sets);

  printf("after ex_get_init, error = %d\n", error);

  /* write initialization parameters */

  error = ex_put_init(exoid2, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
                      num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

  /* read nodal coordinate values */

  float *x = (float *)calloc(num_nodes, sizeof(float));
  float *y = (float *)calloc(num_nodes, sizeof(float));
  float *z = (num_dim >= 3) ? (float *)calloc(num_nodes, sizeof(float)) : NULL;

  error = ex_get_coord(exoid, x, y, z);
  printf("\nafter ex_get_coord, error = %3d\n", error);

  /* write nodal coordinate values */

  error = ex_put_coord(exoid2, x, y, z);
  printf("after ex_put_coord, error = %d\n", error);

  free(x);
  free(y);
  if (num_dim >= 3) {
    free(z);
  }

  /* read nodal coordinate names */

  char *coord_names[3];
  for (int i = 0; i < num_dim; i++) {
    coord_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_coord_names(exoid, coord_names);
  printf("\nafter ex_get_coord_names, error = %3d\n", error);

  /* write nodal coordinate names */

  error = ex_put_coord_names(exoid2, coord_names);
  printf("after ex_put_coord_names, error = %d\n", error);

  for (int i = 0; i < num_dim; i++) {
    free(coord_names[i]);
  }

  /* read element order map */

  int *elem_map = (int *)calloc(num_elem, sizeof(int));

  error = ex_get_id_map(exoid, EX_ELEM_MAP, elem_map);
  printf("\nafter ex_get_id_map, error = %3d\n", error);

  /* write element order map */

  error = ex_put_id_map(exoid2, EX_ELEM_MAP, elem_map);
  printf("after ex_put_id_map, error = %d\n", error);

  free(elem_map);

  /* read and write element block parameters and element connectivity */

  int *ids = (int *)calloc(num_elem_blk, sizeof(int));
  error    = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
  printf("\nafter ex_get_elem_blk_ids, error = %3d\n", error);

  float attrib[] = {3.14159};
  for (int i = 0; i < num_elem_blk; i++) {
    error = ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &num_elem_in_block,
                         &num_nodes_per_elem, 0, 0, &num_attr);
    printf("\nafter ex_get_elem_block, error = %d\n", error);

    error = ex_put_block(exoid2, EX_ELEM_BLOCK, ids[i], elem_type, num_elem_in_block,
                         num_nodes_per_elem, 0, 0, num_attr);
    printf("after ex_put_elem_block, error = %d\n", error);

    int *connect = (int *)calloc((num_nodes_per_elem * num_elem_in_block), sizeof(int));

    error = ex_get_conn(exoid, EX_ELEM_BLOCK, ids[i], connect, NULL, NULL);
    printf("\nafter ex_get_elem_conn, error = %d\n", error);

    error = ex_put_conn(exoid2, EX_ELEM_BLOCK, ids[i], connect, NULL, NULL);
    printf("after ex_put_elem_conn, error = %d\n", error);

    /* write element block attributes */
    error = ex_put_attr(exoid2, EX_ELEM_BLOCK, ids[i], attrib);
    printf("after ex_put_elem_attr, error = %d\n", error);

    free(connect);
  }

  /* read and write element block properties */

  int num_props = ex_inquire_int(exoid, EX_INQ_EB_PROP);
  printf("\nafter ex_inquire, error = %d\n", error);

  char *prop_names[3];
  for (int i = 0; i < num_props; i++) {
    prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_prop_names(exoid, EX_ELEM_BLOCK, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  error = ex_put_prop_names(exoid2, EX_ELEM_BLOCK, num_props, prop_names);
  printf("after ex_put_prop_names, error = %d\n", error);

  for (int i = 0; i < num_props; i++) {
    for (int j = 0; j < num_elem_blk; j++) {
      int prop_value;
      error = ex_get_prop(exoid, EX_ELEM_BLOCK, ids[j], prop_names[i], &prop_value);
      printf("after ex_get_prop, error = %d\n", error);

      if (i > 0) { /* first property is the ID which is already stored */
        error = ex_put_prop(exoid2, EX_ELEM_BLOCK, ids[j], prop_names[i], prop_value);
        printf("after ex_put_prop, error = %d\n", error);
      }
    }
  }

  for (int i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }

  free(ids);

  /* read and write individual node sets */

  ids = (int *)calloc(num_node_sets, sizeof(int));

  error = ex_get_ids(exoid, EX_NODE_SET, ids);
  printf("\nafter ex_get_node_set_ids, error = %3d\n", error);

  for (int i = 0; i < num_node_sets; i++) {
    error = ex_get_set_param(exoid, EX_NODE_SET, ids[i], &num_nodes_in_set, &num_df_in_set);
    printf("\nafter ex_get_node_set_param, error = %3d\n", error);

    error = ex_put_set_param(exoid2, EX_NODE_SET, ids[i], num_nodes_in_set, num_df_in_set);
    printf("after ex_put_node_set_param, error = %d\n", error);

    int   *node_list = (int *)calloc(num_nodes_in_set, sizeof(int));
    float *dist_fact = (float *)calloc(num_nodes_in_set, sizeof(float));

    error = ex_get_set(exoid, EX_NODE_SET, ids[i], node_list, NULL);
    printf("\nafter ex_get_node_set, error = %3d\n", error);

    error = ex_put_set(exoid2, EX_NODE_SET, ids[i], node_list, NULL);
    printf("after ex_put_node_set, error = %d\n", error);

    if (num_df_in_set > 0) {
      error = ex_get_set_dist_fact(exoid, EX_NODE_SET, ids[i], dist_fact);
      printf("\nafter ex_get_node_set_dist_fact, error = %3d\n", error);

      error = ex_put_set_dist_fact(exoid2, EX_NODE_SET, ids[i], dist_fact);
      printf("after ex_put_node_set, error = %d\n", error);
    }

    free(node_list);
    free(dist_fact);
  }
  free(ids);

  /* read node set properties */
  num_props = ex_inquire_int(exoid, EX_INQ_NS_PROP);
  printf("\nafter ex_inquire, error = %d\n", error);

  for (int i = 0; i < num_props; i++) {
    prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }
  int *prop_values = (int *)calloc(num_node_sets, sizeof(int));

  error = ex_get_prop_names(exoid, EX_NODE_SET, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  error = ex_put_prop_names(exoid2, EX_NODE_SET, num_props, prop_names);
  printf("after ex_put_prop_names, error = %d\n", error);

  for (int i = 0; i < num_props; i++) {
    error = ex_get_prop_array(exoid, EX_NODE_SET, prop_names[i], prop_values);
    printf("after ex_get_prop_array, error = %d\n", error);

    error = ex_put_prop_array(exoid2, EX_NODE_SET, prop_names[i], prop_values);
    printf("after ex_put_prop_array, error = %d\n", error);
  }
  for (int i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }
  free(prop_values);

  /* read and write individual side sets */

  ids = (int *)calloc(num_side_sets, sizeof(int));

  error = ex_get_ids(exoid, EX_SIDE_SET, ids);
  printf("\nafter ex_get_side_set_ids, error = %3d\n", error);

  for (int i = 0; i < num_side_sets; i++) {
    error = ex_get_set_param(exoid, EX_SIDE_SET, ids[i], &num_sides_in_set, &num_df_in_set);
    printf("\nafter ex_get_side_set_param, error = %3d\n", error);

    error = ex_put_set_param(exoid2, EX_SIDE_SET, ids[i], num_sides_in_set, num_df_in_set);
    printf("after ex_put_side_set_param, error = %d\n", error);

    /* Note: The # of elements is same as # of sides!  */
    int    num_elem_in_set = num_sides_in_set;
    int   *elem_list       = (int *)calloc(num_elem_in_set, sizeof(int));
    int   *side_list       = (int *)calloc(num_sides_in_set, sizeof(int));
    int   *node_ctr_list   = (int *)calloc(num_elem_in_set, sizeof(int));
    int   *node_list       = (int *)calloc(num_elem_in_set * 21, sizeof(int));
    float *dist_fact       = (float *)calloc(num_df_in_set, sizeof(float));

    error = ex_get_set(exoid, EX_SIDE_SET, ids[i], elem_list, side_list);
    printf("\nafter ex_get_side_set, error = %3d\n", error);

    error = ex_put_set(exoid2, EX_SIDE_SET, ids[i], elem_list, side_list);
    printf("after ex_put_side_set, error = %d\n", error);

    error = ex_get_side_set_node_list(exoid, ids[i], node_ctr_list, node_list);
    printf("\nafter ex_get_side_set_node_list, error = %3d\n", error);

    if (num_df_in_set > 0) {
      error = ex_get_set_dist_fact(exoid, EX_SIDE_SET, ids[i], dist_fact);
      printf("\nafter ex_get_side_set_dist_fact, error = %3d\n", error);

      error = ex_put_set_dist_fact(exoid2, EX_SIDE_SET, ids[i], dist_fact);
      printf("after ex_put_side_set_dist_fact, error = %d\n", error);
    }

    free(elem_list);
    free(side_list);
    free(node_ctr_list);
    free(node_list);
    free(dist_fact);
  }

  /* read side set properties */
  num_props = ex_inquire_int(exoid, EX_INQ_SS_PROP);
  printf("\nafter ex_inquire, error = %d\n", error);

  for (int i = 0; i < num_props; i++) {
    prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_prop_names(exoid, EX_SIDE_SET, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  for (int i = 0; i < num_props; i++) {
    for (int j = 0; j < num_side_sets; j++) {
      int prop_value;
      error = ex_get_prop(exoid, EX_SIDE_SET, ids[j], prop_names[i], &prop_value);
      printf("after ex_get_prop, error = %d\n", error);

      if (i > 0) { /* first property is ID so it is already stored */
        error = ex_put_prop(exoid2, EX_SIDE_SET, ids[j], prop_names[i], prop_value);
        printf("after ex_put_prop, error = %d\n", error);
      }
    }
  }
  for (int i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }
  free(ids);

  /* read and write QA records */

  int num_qa_rec = ex_inquire_int(exoid, EX_INQ_QA);

  char *qa_record[2][4];
  for (int i = 0; i < num_qa_rec; i++) {
    for (int j = 0; j < 4; j++) {
      qa_record[i][j] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }
  }

  error = ex_get_qa(exoid, qa_record);
  printf("\nafter ex_get_qa, error = %3d\n", error);

  error = ex_put_qa(exoid2, num_qa_rec, qa_record);
  printf("after ex_put_qa, error = %d\n", error);

  for (int i = 0; i < num_qa_rec; i++) {
    for (int j = 0; j < 4; j++) {
      free(qa_record[i][j]);
    }
  }
  /* read and write information records */

  int num_info = ex_inquire_int(exoid, EX_INQ_INFO);
  printf("\nafter ex_inquire, error = %3d\n", error);

  char *info[3];
  for (int i = 0; i < num_info; i++) {
    info[i] = (char *)calloc((MAX_LINE_LENGTH + 1), sizeof(char));
  }

  error = ex_get_info(exoid, info);
  printf("\nafter ex_get_info, error = %3d\n", error);

  error = ex_put_info(exoid2, num_info, info);
  printf("after ex_put_info, error = %d\n", error);

  for (int i = 0; i < num_info; i++) {
    free(info[i]);
  }

  /* close the EXODUS files */

  error = ex_close(exoid);
  printf("after ex_close, error = %d\n", error);
  error = ex_close(exoid2);
  printf("after ex_close (2), error = %d\n", error);
  return 0;
}
