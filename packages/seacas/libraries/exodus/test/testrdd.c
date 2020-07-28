/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/******************************************************************************
 * testrdd - read exodus file test.exo created by testwt - double precision
 *
 * author - Sandia National Laboratories
 *          Larry A. Schoof - Original
 *
 *
 * environment - UNIX
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *
 * exit conditions -
 *
 * revision history -
 *
 *   Original L. A. Schoof
 *   04/05/93 V.R. Yarberry - revised so that output resembles Fortran version
 *   06/25/93 VRY - revised to match 2.00 API.
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

void *my_calloc(size_t length, size_t size)
{
  if (length == 0 || size == 0) {
    return NULL;
  }
  return calloc(length, size);
}

int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets;
  int  num_side_sets, error;
  int  i, j, k, node_ctr;
  int *elem_map, *connect, *node_list, *node_ctr_list, *elem_list, *side_list;
  int *ids;
  int *num_nodes_per_set, *num_elem_per_set;
  int *num_df_per_set;
  int *node_ind, *elem_ind, *df_ind, num_qa_rec, num_info;
  int  num_glo_vars, num_nod_vars, num_ele_vars;
  int *truth_tab;
  int  num_time_steps;
  int *num_elem_in_block, *num_nodes_per_elem, *num_attr;
  int  num_nodes_in_set, num_elem_in_set;
  int  num_sides_in_set, num_df_in_set;
  int  list_len      = 0;
  int  elem_list_len = 0;
  int  node_list_len = 0;
  int  df_list_len   = 0;
  int  node_num, time_step, var_index, beg_time, end_time, elem_num;
  int  CPU_word_size, IO_word_size;
  int  num_props, prop_value, *prop_values;

  double  time_value, *time_values, *var_values;
  double *x, *y, *z;
  double  attrib[1], *dist_fact;
  float   version, fdum;

  char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char  title[MAX_LINE_LENGTH + 1], elem_type[MAX_STR_LENGTH + 1];
  char *cdum = 0;
  char *prop_names[3];

  CPU_word_size = 8; /* sizeof(double) */
  IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */

  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  printf("\nafter ex_open\n");
  if (exoid < 0) {
    exit(1);
  }

  printf("test.exo is an EXODUSII file; version %4.2f\n", version);
  printf("         CPU word size %1d\n", CPU_word_size);
  printf("         I/O word size %1d\n", IO_word_size);

  /* read database parameters */

  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);

  printf("after ex_get_init, error = %3d\n", error);

  printf("database parameters:\n");
  printf("title =  '%s'\n", title);
  printf("num_dim = %3d\n", num_dim);
  printf("num_nodes = %3d\n", num_nodes);
  printf("num_elem = %3d\n", num_elem);
  printf("num_elem_blk = %3d\n", num_elem_blk);
  printf("num_node_sets = %3d\n", num_node_sets);
  printf("num_side_sets = %3d\n", num_side_sets);

  /* read nodal coordinates values and names from database */

  x = (double *)my_calloc(num_nodes, sizeof(double));
  y = (double *)my_calloc(num_nodes, sizeof(double));
  if (num_dim >= 3) {
    z = (double *)my_calloc(num_nodes, sizeof(double));
  }
  else {
    z = 0;
  }

  error = ex_get_coord(exoid, x, y, z);
  printf("\nafter ex_get_coord, error = %3d\n", error);

  printf("x coords = \n");
  for (i = 0; i < num_nodes; i++) {
    printf("%5.1f\n", x[i]);
  }

  printf("y coords = \n");
  for (i = 0; i < num_nodes; i++) {
    printf("%5.1f\n", y[i]);
  }

  if (num_dim >= 3) {
    printf("z coords = \n");
    for (i = 0; i < num_nodes; i++) {
      printf("%5.1f\n", z[i]);
    }
  }

  /*
     error = ex_get_1_coord (exoid, 2, x, y, z);
     printf ("\nafter ex_get_1_coord, error = %3d\n", error);

     printf ("x coord of node 2 = \n");
     printf ("%f \n", x[0]);

     printf ("y coord of node 2 = \n");
     printf ("%f \n", y[0]);
  */
  free(x);
  free(y);
  if (num_dim >= 3) {
    free(z);
  }

  for (i = 0; i < num_dim; i++) {
    coord_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_coord_names(exoid, coord_names);
  printf("\nafter ex_get_coord_names, error = %3d\n", error);

  printf("x coord name = '%s'\n", coord_names[0]);
  printf("y coord name = '%s'\n", coord_names[1]);

  for (i = 0; i < num_dim; i++) {
    free(coord_names[i]);
  }

  /* read element order map */

  elem_map = (int *)my_calloc(num_elem, sizeof(int));

  error = ex_get_map(exoid, elem_map);
  printf("\nafter ex_get_map, error = %3d\n", error);

  for (i = 0; i < num_elem; i++) {
    printf("elem_map(%d) = %d \n", i, elem_map[i]);
  }

  free(elem_map);

  /* read element block parameters */

  ids                = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_elem_in_block  = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_nodes_per_elem = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_attr           = (int *)my_calloc(num_elem_blk, sizeof(int));

  error = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
  printf("\nafter ex_get_elem_blk_ids, error = %3d\n", error);

  for (i = 0; i < num_elem_blk; i++) {
    error = ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &(num_elem_in_block[i]),
                         &(num_nodes_per_elem[i]), 0, 0, &(num_attr[i]));
    printf("\nafter ex_get_elem_block, error = %d\n", error);

    printf("element block id = %2d\n", ids[i]);
    printf("element type = '%s'\n", elem_type);
    printf("num_elem_in_block = %2d\n", num_elem_in_block[i]);
    printf("num_nodes_per_elem = %2d\n", num_nodes_per_elem[i]);
    printf("num_attr = %2d\n", num_attr[i]);
  }

  /* read element block properties */
  error = ex_inquire(exoid, EX_INQ_EB_PROP, &num_props, &fdum, cdum);
  printf("\nafter ex_inquire, error = %d\n", error);
  printf("\nThere are %2d properties for each element block\n", num_props);

  for (i = 0; i < num_props; i++) {
    prop_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_prop_names(exoid, EX_ELEM_BLOCK, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  for (i = 0; i < num_props; i++) {
    for (j = 0; j < num_elem_blk; j++) {
      error = ex_get_prop(exoid, EX_ELEM_BLOCK, ids[j], prop_names[i], &prop_value);
      if (error == 0) {
        printf("element block %2d, property(%2d): '%s'= %5d\n", j + 1, i + 1, prop_names[i],
               prop_value);
      }
      else {
        printf("after ex_get_prop, error = %d\n", error);
      }
    }
  }

  for (i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }

  /* read element connectivity */

  for (i = 0; i < num_elem_blk; i++) {
    connect = (int *)my_calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));

    error = ex_get_conn(exoid, EX_ELEM_BLOCK, ids[i], connect, NULL, NULL);
    printf("\nafter ex_get_elem_conn, error = %d\n", error);

    printf("connect array for elem block %2d\n", ids[i]);

    for (j = 0; j < num_nodes_per_elem[i]; j++) {
      printf("%3d\n", connect[j]);
    }
    free(connect);
  }

  /* read element block attributes */

  for (i = 0; i < num_elem_blk; i++) {
    error = ex_get_attr(exoid, EX_ELEM_BLOCK, ids[i], attrib);
    printf("\n after ex_get_elem_attr, error = %d\n", error);

    if (error == 0) {
      printf("element block %d attributes = %10.8f\n", ids[i], *attrib);
    }
  }

  free(ids);
  free(num_nodes_per_elem);
  free(num_attr);

  /* read individual node sets */

  ids = (int *)my_calloc(num_node_sets, sizeof(int));

  error = ex_get_ids(exoid, EX_NODE_SET, ids);
  printf("\nafter ex_get_node_set_ids, error = %3d\n", error);

  for (i = 0; i < num_node_sets; i++) {
    error = ex_get_set_param(exoid, EX_NODE_SET, ids[i], &num_nodes_in_set, &num_df_in_set);
    printf("\nafter ex_get_node_set_param, error = %3d\n", error);

    printf("\nnode set %2d parameters: \n", ids[i]);
    printf("num_nodes = %2d\n", num_nodes_in_set);

    node_list = (int *)my_calloc(num_nodes_in_set, sizeof(int));
    dist_fact = (double *)my_calloc(num_nodes_in_set, sizeof(double));

    error = ex_get_set(exoid, EX_NODE_SET, ids[i], node_list, NULL);
    printf("\nafter ex_get_node_set, error = %3d\n", error);

    if (num_df_in_set > 0) {
      error = ex_get_set_dist_fact(exoid, EX_NODE_SET, ids[i], dist_fact);
      printf("\nafter ex_get_node_set_dist_fact, error = %3d\n", error);
    }

    printf("\nnode list for node set %2d\n", ids[i]);

    for (j = 0; j < num_nodes_in_set; j++) {
      printf("%3d\n", node_list[j]);
    }

    if (num_df_in_set > 0) {
      printf("dist factors for node set %2d\n", ids[i]);

      for (j = 0; j < num_df_in_set; j++) {
        printf("%5.2f\n", dist_fact[j]);
      }
    }
    else {
      printf("no dist factors for node set %2d\n", ids[i]);
    }

    free(node_list);
    free(dist_fact);
  }
  free(ids);

  /* read node set properties */
  error = ex_inquire(exoid, EX_INQ_NS_PROP, &num_props, &fdum, cdum);
  printf("\nafter ex_inquire, error = %d\n", error);
  printf("\nThere are %2d properties for each node set\n", num_props);

  for (i = 0; i < num_props; i++) {
    prop_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }
  prop_values = (int *)my_calloc(num_node_sets, sizeof(int));

  error = ex_get_prop_names(exoid, EX_NODE_SET, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  for (i = 0; i < num_props; i++) {
    error = ex_get_prop_array(exoid, EX_NODE_SET, prop_names[i], prop_values);
    if (error == 0) {
      for (j = 0; j < num_node_sets; j++) {
        printf("node set %2d, property(%2d): '%s'= %5d\n", j + 1, i + 1, prop_names[i],
               prop_values[j]);
      }
    }
    else {
      printf("after ex_get_prop_array, error = %d\n", error);
    }
  }
  for (i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }
  free(prop_values);

  /* read concatenated node sets; this produces the same information as
   * the above code which reads individual node sets
   */

  error = ex_inquire(exoid, EX_INQ_NODE_SETS, &num_node_sets, &fdum, cdum);
  printf("\nafter ex_inquire, error = %3d\n", error);

  ids               = (int *)my_calloc(num_node_sets, sizeof(int));
  num_nodes_per_set = (int *)my_calloc(num_node_sets, sizeof(int));
  num_df_per_set    = (int *)my_calloc(num_node_sets, sizeof(int));
  node_ind          = (int *)my_calloc(num_node_sets, sizeof(int));
  df_ind            = (int *)my_calloc(num_node_sets, sizeof(int));

  error = ex_inquire(exoid, EX_INQ_NS_NODE_LEN, &list_len, &fdum, cdum);
  printf("\nafter ex_inquire: EX_INQ_NS_NODE_LEN = %d, error = %3d\n", list_len, error);

  node_list = (int *)my_calloc(list_len, sizeof(int));

  error = ex_inquire(exoid, EX_INQ_NS_DF_LEN, &list_len, &fdum, cdum);
  printf("\nafter ex_inquire: EX_INQ_NS_DF_LEN = %d, error = %3d\n", list_len, error);
  dist_fact = (double *)my_calloc(list_len, sizeof(double));

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

    error = ex_get_concat_sets(exoid, EX_NODE_SET, &set_specs);
  }
  printf("\nafter ex_get_concat_node_sets, error = %3d\n", error);

  printf("\nconcatenated node set info\n");

  printf("ids = \n");
  for (i = 0; i < num_node_sets; i++) {
    printf("%3d\n", ids[i]);
  }

  printf("num_nodes_per_set = \n");
  for (i = 0; i < num_node_sets; i++) {
    printf("%3d\n", num_nodes_per_set[i]);
  }

  printf("node_ind = \n");
  for (i = 0; i < num_node_sets; i++) {
    printf("%3d\n", node_ind[i]);
  }

  printf("node_list = \n");
  for (i = 0; i < list_len; i++) {
    printf("%3d\n", node_list[i]);
  }

  printf("dist_fact = \n");
  for (i = 0; i < list_len; i++) {
    printf("%5.3f\n", dist_fact[i]);
  }

  free(ids);
  free(num_nodes_per_set);
  free(df_ind);
  free(node_ind);
  free(num_df_per_set);
  free(node_list);
  free(dist_fact);

  /* read individual side sets */

  ids = (int *)my_calloc(num_side_sets, sizeof(int));

  error = ex_get_ids(exoid, EX_SIDE_SET, ids);
  printf("\nafter ex_get_side_set_ids, error = %3d\n", error);

  for (i = 0; i < num_side_sets; i++) {
    error = ex_get_set_param(exoid, EX_SIDE_SET, ids[i], &num_sides_in_set, &num_df_in_set);
    printf("\nafter ex_get_side_set_param, error = %3d\n", error);

    printf("side set %2d parameters:\n", ids[i]);
    printf("num_sides = %3d\n", num_sides_in_set);
    printf("num_dist_factors = %3d\n", num_df_in_set);

    /* Note: The # of elements is same as # of sides!  */
    num_elem_in_set = num_sides_in_set;
    elem_list       = (int *)my_calloc(num_elem_in_set, sizeof(int));
    side_list       = (int *)my_calloc(num_sides_in_set, sizeof(int));
    node_ctr_list   = (int *)my_calloc(num_elem_in_set, sizeof(int));
    node_list       = (int *)my_calloc(num_elem_in_set * 21, sizeof(int));
    dist_fact       = (double *)my_calloc(num_df_in_set, sizeof(double));

    error = ex_get_set(exoid, EX_SIDE_SET, ids[i], elem_list, side_list);
    printf("\nafter ex_get_side_set, error = %3d\n", error);

    error = ex_get_side_set_node_list(exoid, ids[i], node_ctr_list, node_list);
    printf("\nafter ex_get_side_set_node_list, error = %3d\n", error);

    if (num_df_in_set > 0) {
      error = ex_get_set_dist_fact(exoid, EX_SIDE_SET, ids[i], dist_fact);
      printf("\nafter ex_get_side_set_dist_fact, error = %3d\n", error);
    }

    printf("element list for side set %2d\n", ids[i]);

    for (j = 0; j < num_elem_in_set; j++) {
      printf("%3d\n", elem_list[j]);
    }

    printf("side list for side set %2d\n", ids[i]);
    for (j = 0; j < num_sides_in_set; j++) {
      printf("%3d\n", side_list[j]);
    }

    node_ctr = 0;
    printf("node list for side set %2d\n", ids[i]);
    for (k = 0; k < num_elem_in_set; k++) {
      for (j = 0; j < node_ctr_list[k]; j++) {
        printf("%3d\n", node_list[node_ctr + j]);
      }
      node_ctr += node_ctr_list[k];
    }

    if (num_df_in_set > 0) {
      printf("dist factors for side set %2d\n", ids[i]);

      for (j = 0; j < num_df_in_set; j++) {
        printf("%5.3f\n", dist_fact[j]);
      }
    }
    else {
      printf("no dist factors for side set %2d\n", ids[i]);
    }

    free(elem_list);
    free(side_list);
    free(node_ctr_list);
    free(node_list);
    free(dist_fact);
  }

  /* read side set properties */
  error = ex_inquire(exoid, EX_INQ_SS_PROP, &num_props, &fdum, cdum);
  printf("\nafter ex_inquire, error = %d\n", error);
  printf("\nThere are %2d properties for each side set\n", num_props);

  for (i = 0; i < num_props; i++) {
    prop_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_prop_names(exoid, EX_SIDE_SET, prop_names);
  printf("after ex_get_prop_names, error = %d\n", error);

  for (i = 0; i < num_props; i++) {
    for (j = 0; j < num_side_sets; j++) {
      error = ex_get_prop(exoid, EX_SIDE_SET, ids[j], prop_names[i], &prop_value);
      if (error == 0) {
        printf("side set %2d, property(%2d): '%s'= %5d\n", j + 1, i + 1, prop_names[i], prop_value);
      }
      else {
        printf("after ex_get_prop, error = %d\n", error);
      }
    }
  }
  for (i = 0; i < num_props; i++) {
    free(prop_names[i]);
  }
  free(ids);

  error = ex_inquire(exoid, EX_INQ_SIDE_SETS, &num_side_sets, &fdum, cdum);
  printf("\nafter ex_inquire: EX_INQ_SIDE_SETS = %d,  error = %d\n", num_side_sets, error);

  if (num_side_sets > 0) {
    error = ex_inquire(exoid, EX_INQ_SS_ELEM_LEN, &elem_list_len, &fdum, cdum);
    printf("\nafter ex_inquire: EX_INQ_SS_ELEM_LEN = %d,  error = %d\n", elem_list_len, error);

    error = ex_inquire(exoid, EX_INQ_SS_NODE_LEN, &node_list_len, &fdum, cdum);
    printf("\nafter ex_inquire: EX_INQ_SS_NODE_LEN = %d,  error = %d\n", node_list_len, error);

    error = ex_inquire(exoid, EX_INQ_SS_DF_LEN, &df_list_len, &fdum, cdum);
    printf("\nafter ex_inquire: EX_INQ_SS_DF_LEN = %d,  error = %d\n", df_list_len, error);
  }

  /* read concatenated side sets; this produces the same information as
   * the above code which reads individual side sets
   */

  /* concatenated side set read */

  ids              = (int *)my_calloc(num_side_sets, sizeof(int));
  num_elem_per_set = (int *)my_calloc(num_side_sets, sizeof(int));
  num_df_per_set   = (int *)my_calloc(num_side_sets, sizeof(int));
  elem_ind         = (int *)my_calloc(num_side_sets, sizeof(int));
  df_ind           = (int *)my_calloc(num_side_sets, sizeof(int));
  elem_list        = (int *)my_calloc(elem_list_len, sizeof(int));
  side_list        = (int *)my_calloc(elem_list_len, sizeof(int));
  dist_fact        = (double *)my_calloc(df_list_len, sizeof(double));

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

    error = ex_get_concat_sets(exoid, EX_SIDE_SET, &set_specs);
  }

  printf("\nafter ex_get_concat_side_sets, error = %3d\n", error);

  printf("concatenated side set info\n");

  printf("ids = \n");
  for (i = 0; i < num_side_sets; i++) {
    printf("%3d\n", ids[i]);
  }

  printf("num_elem_per_set = \n");
  for (i = 0; i < num_side_sets; i++) {
    printf("%3d\n", num_elem_per_set[i]);
  }

  printf("num_dist_per_set = \n");
  for (i = 0; i < num_side_sets; i++) {
    printf("%3d\n", num_df_per_set[i]);
  }

  printf("elem_ind = \n");
  for (i = 0; i < num_side_sets; i++) {
    printf("%3d\n", elem_ind[i]);
  }

  printf("dist_ind = \n");
  for (i = 0; i < num_side_sets; i++) {
    printf("%3d\n", df_ind[i]);
  }

  printf("elem_list = \n");
  for (i = 0; i < elem_list_len; i++) {
    printf("%3d\n", elem_list[i]);
  }

  printf("side_list = \n");
  for (i = 0; i < elem_list_len; i++) {
    printf("%3d\n", side_list[i]);
  }

  printf("dist_fact = \n");
  for (i = 0; i < df_list_len; i++) {
    printf("%5.3f\n", dist_fact[i]);
  }

  free(ids);
  free(num_elem_per_set);
  free(num_df_per_set);
  free(df_ind);
  free(elem_ind);
  free(elem_list);
  free(side_list);
  free(dist_fact);

  /* end of concatenated side set read */

  /* read QA records */

  ex_inquire(exoid, EX_INQ_QA, &num_qa_rec, &fdum, cdum);

  for (i = 0; i < num_qa_rec; i++) {
    for (j = 0; j < 4; j++) {
      qa_record[i][j] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }
  }

  error = ex_get_qa(exoid, qa_record);
  printf("\nafter ex_get_qa, error = %3d\n", error);

  printf("QA records = \n");
  for (i = 0; i < num_qa_rec; i++) {
    for (j = 0; j < 4; j++) {
      printf(" '%s'\n", qa_record[i][j]);
      free(qa_record[i][j]);
    }
  }

  /* read information records */

  error = ex_inquire(exoid, EX_INQ_INFO, &num_info, &fdum, cdum);
  printf("\nafter ex_inquire, error = %3d\n", error);

  for (i = 0; i < num_info; i++) {
    info[i] = (char *)my_calloc((MAX_LINE_LENGTH + 1), sizeof(char));
  }

  error = ex_get_info(exoid, info);
  printf("\nafter ex_get_info, error = %3d\n", error);

  printf("info records = \n");
  for (i = 0; i < num_info; i++) {
    printf(" '%s'\n", info[i]);
    free(info[i]);
  }

  /* read global variables parameters and names */

  error = ex_get_variable_param(exoid, EX_GLOBAL, &num_glo_vars);
  printf("\nafter ex_get_variable_param, error = %3d\n", error);

  for (i = 0; i < num_glo_vars; i++) {
    var_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names);
  printf("\nafter ex_get_variable_names, error = %3d\n", error);

  printf("There are %2d global variables; their names are :\n", num_glo_vars);
  for (i = 0; i < num_glo_vars; i++) {
    printf(" '%s'\n", var_names[i]);
    free(var_names[i]);
  }

  /* read nodal variables parameters and names */

  error = ex_get_variable_param(exoid, EX_NODAL, &num_nod_vars);
  printf("\nafter ex_get_variable_param, error = %3d\n", error);

  for (i = 0; i < num_nod_vars; i++) {
    var_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
  printf("\nafter ex_get_variable_names, error = %3d\n", error);

  printf("There are %2d nodal variables; their names are :\n", num_nod_vars);
  for (i = 0; i < num_nod_vars; i++) {
    printf(" '%s'\n", var_names[i]);
    free(var_names[i]);
  }

  /* read element variables parameters and names */

  error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_ele_vars);
  printf("\nafter ex_get_variable_param, error = %3d\n", error);

  for (i = 0; i < num_ele_vars; i++) {
    var_names[i] = (char *)my_calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
  printf("\nafter ex_get_variable_names, error = %3d\n", error);

  printf("There are %2d element variables; their names are :\n", num_ele_vars);
  for (i = 0; i < num_ele_vars; i++) {
    printf(" '%s'\n", var_names[i]);
    free(var_names[i]);
  }

  /* read element variable truth table */

  truth_tab = (int *)my_calloc((num_elem_blk * num_ele_vars), sizeof(int));

  error = ex_get_truth_table(exoid, EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab);
  printf("\nafter ex_get_elem_var_tab, error = %3d\n", error);

  printf("This is the element variable truth table:\n");

  k = 0;
  for (i = 0; i < num_elem_blk * num_ele_vars; i++) {
    printf("%2d\n", truth_tab[k++]);
  }
  free(truth_tab);

  /* determine how many time steps are stored */

  error = ex_inquire(exoid, EX_INQ_TIME, &num_time_steps, &fdum, cdum);
  printf("\nafter ex_inquire, error = %3d\n", error);
  printf("There are %2d time steps in the database.\n", num_time_steps);

  /* read time value at one time step */

  time_step = 3;
  error     = ex_get_time(exoid, time_step, &time_value);
  printf("\nafter ex_get_time, error = %3d\n", error);

  printf("time value at time step %2d = %5.3f\n", time_step, time_value);

  /* read time values at all time steps */

  time_values = (double *)my_calloc(num_time_steps, sizeof(double));

  error = ex_get_all_times(exoid, time_values);
  printf("\nafter ex_get_all_times, error = %3d\n", error);

  printf("time values at all time steps are:\n");
  for (i = 0; i < num_time_steps; i++) {
    printf("%5.3f\n", time_values[i]);
  }

  free(time_values);

  /* read all global variables at one time step */

  var_values = (double *)my_calloc(num_glo_vars, sizeof(double));

  error = ex_get_var(exoid, time_step, EX_GLOBAL, 1, 1, num_glo_vars, var_values);
  printf("\nafter ex_get_glob_vars, error = %3d\n", error);

  printf("global variable values at time step %2d\n", time_step);
  for (i = 0; i < num_glo_vars; i++) {
    printf("%5.3f\n", var_values[i]);
  }

  free(var_values);

  /* read a single global variable through time */

  var_index = 1;
  beg_time  = 1;
  end_time  = -1;

  var_values = (double *)my_calloc(num_time_steps, sizeof(double));

  error = ex_get_var_time(exoid, EX_GLOBAL, var_index, 1, beg_time, end_time, var_values);
  printf("\nafter ex_get_glob_var_time, error = %3d\n", error);

  printf("global variable %2d values through time:\n", var_index);
  for (i = 0; i < num_time_steps; i++) {
    printf("%5.3f\n", var_values[i]);
  }

  free(var_values);

  /* read a nodal variable at one time step */

  var_values = (double *)my_calloc(num_nodes, sizeof(double));

  error = ex_get_var(exoid, time_step, EX_NODAL, var_index, 1, num_nodes, var_values);
  printf("\nafter ex_get_nodal_var, error = %3d\n", error);

  printf("nodal variable %2d values at time step %2d\n", var_index, time_step);
  for (i = 0; i < num_nodes; i++) {
    printf("%5.3f\n", var_values[i]);
  }

  free(var_values);

  /* read a nodal variable through time */

  var_values = (double *)my_calloc(num_time_steps, sizeof(double));

  node_num = 1;
  error    = ex_get_var_time(exoid, EX_NODAL, var_index, node_num, beg_time, end_time, var_values);
  printf("\nafter ex_get_nodal_var_time, error = %3d\n", error);

  printf("nodal variable %2d values for node %2d through time:\n", var_index, node_num);
  for (i = 0; i < num_time_steps; i++) {
    printf("%5.3f\n", var_values[i]);
  }

  free(var_values);

  /* read an element variable at one time step */

  ids = (int *)my_calloc(num_elem_blk, sizeof(int));

  error = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
  printf("\n after ex_get_elem_blk_ids, error = %3d\n", error);

  for (i = 0; i < num_elem_blk; i++) {
    var_values = (double *)my_calloc(num_elem_in_block[i], sizeof(double));

    error = ex_get_var(exoid, time_step, EX_ELEM_BLOCK, var_index, ids[i], num_elem_in_block[i],
                       var_values);
    printf("\nafter ex_get_elem_var, error = %3d\n", error);

    if (!error) {
      printf("element variable %2d values of element block %2d at time step %2d\n", var_index,
             ids[i], time_step);
      for (j = 0; j < num_elem_in_block[i]; j++) {
        printf("%5.3f\n", var_values[j]);
      }
    }

    free(var_values);
  }
  free(num_elem_in_block);
  free(ids);

  /* read an element variable through time */

  var_values = (double *)my_calloc(num_time_steps, sizeof(double));

  var_index = 2;
  elem_num  = 2;
  error =
      ex_get_var_time(exoid, EX_ELEM_BLOCK, var_index, elem_num, beg_time, end_time, var_values);
  printf("\nafter ex_get_elem_var_time, error = %3d\n", error);

  printf("element variable %2d values for element %2d through time:\n", var_index, elem_num);
  for (i = 0; i < num_time_steps; i++) {
    printf("%5.3f\n", var_values[i]);
  }

  free(var_values);

  error = ex_close(exoid);
  printf("\nafter ex_close, error = %3d\n", error);
  return 0;
}
