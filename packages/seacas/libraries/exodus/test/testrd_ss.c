/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testrd_ss - read exodus file test.exo created by testwt_ss
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
  int *connect, *node_list, *node_ctr_list, *elem_list, *side_list;
  int *ids;
  int *num_elem_per_set;
  int *num_df_per_set;
  int *elem_ind, *df_ind;
  int *num_elem_in_block, *num_nodes_per_elem, *num_attr;
  int  num_elem_in_set;
  int  num_sides_in_set, num_df_in_set;
  int  elem_list_len = 0;
  int  node_list_len = 0;
  int  df_list_len   = 0;
  int  CPU_word_size, IO_word_size;
  int  idum;

  float *dist_fact;
  float  version, fdum;

  char  title[MAX_LINE_LENGTH + 1], elem_type[MAX_STR_LENGTH + 1];
  char *cdum = 0;

  CPU_word_size = 0; /* sizeof(float) */
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
  /*   printf ("         CPU word size %1d\n",CPU_word_size);  */
  printf("         I/O word size %1d\n", IO_word_size);
  ex_inquire(exoid, EX_INQ_API_VERS, &idum, &version, cdum);
  printf("EXODUSII API; version %4.2f\n", version);

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

  /* read element block parameters */

  ids                = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_elem_in_block  = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_nodes_per_elem = (int *)my_calloc(num_elem_blk, sizeof(int));
  num_attr           = (int *)my_calloc(num_elem_blk, sizeof(int));

  error = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
  printf("\nafter ex_get_elem_blk_ids, error = %3d\n", error);

  for (i = 0; i < num_elem_blk; i++) {
    error = ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &(num_elem_in_block[i]),
                         &(num_nodes_per_elem[i]), NULL, NULL, &(num_attr[i]));
    printf("\nafter ex_get_elem_block, error = %d\n", error);

    printf("element block id = %2d\n", ids[i]);
    printf("element type = '%s'\n", elem_type);
    printf("num_elem_in_block = %2d\n", num_elem_in_block[i]);
    printf("num_nodes_per_elem = %2d\n", num_nodes_per_elem[i]);
    printf("num_attr = %2d\n", num_attr[i]);
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
  free(ids);
  free(num_elem_in_block);
  free(num_nodes_per_elem);
  free(num_attr);

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
    dist_fact       = (float *)my_calloc(num_df_in_set, sizeof(float));

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
      printf("%3d nodes for side %3d\n", node_ctr_list[k], k);
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
  free(ids);

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
  dist_fact        = (float *)my_calloc(df_list_len, sizeof(float));

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

  error = ex_close(exoid);
  printf("\nafter ex_close, error = %3d\n", error);
  return 0;
}
