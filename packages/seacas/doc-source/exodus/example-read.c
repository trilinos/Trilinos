#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  int  num_nset_vars, num_sset_vars;
  int *truth_tab;
  int  num_time_steps;
  int *num_elem_in_block, *num_nodes_per_elem, *num_attr;
  int  num_nodes_in_set, num_elem_in_set;
  int  num_sides_in_set, num_df_in_set;
  int  list_len, elem_list_len, node_list_len, df_list_len;
  int  node_num, time_step, var_index, beg_time, end_time, elem_num;
  int  CPU_word_size, IO_word_size;
  int  num_props, prop_value, *prop_values;
  int  idum;

  ex_set_specs ss_specs, ns_specs;

  float  time_value, *time_values, *var_values;
  float *x, *y, *z;
  float *attrib, *dist_fact;
  float  version;

  char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char *block_names[10], *nset_names[10], *sset_names[10];
  char *attrib_names[10];
  char  name[MAX_STR_LENGTH + 1];
  char  title[MAX_LINE_LENGTH + 1], elem_type[MAX_STR_LENGTH + 1];
  char *cdum;
  char *prop_names[3];

  cdum = 0;

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  if (exoid < 0)
    exit(1);

  printf("test.exo is an EXODUSII file; version %4.2f\n", version);
  printf("         I/O word size %1d\n", IO_word_size);
  ex_inquire(exoid, EX_INQ_API_VERS, &idum, &version, cdum);
  printf("EXODUSII API; version %4.2f\n", version);

  ex_inquire(exoid, EX_INQ_LIB_VERS, &idum, &version, cdum);
  printf("EXODUSII Library API; version %4.2f (%d)\n", version, idum);

  /* read database parameters */
  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);

  printf("database parameters:\n");
  printf("title =  '%s'\n", title);
  printf("num_dim = %3d\n", num_dim);
  printf("num_nodes = %3d\n", num_nodes);
  printf("num_elem = %3d\n", num_elem);
  printf("num_elem_blk = %3d\n", num_elem_blk);
  printf("num_node_sets = %3d\n", num_node_sets);
  printf("num_side_sets = %3d\n", num_side_sets);

  /* read nodal coordinates values and names from database */
  x = (float *)calloc(num_nodes, sizeof(float));
  y = (float *)calloc(num_nodes, sizeof(float));
  if (num_dim >= 3)
    z = (float *)calloc(num_nodes, sizeof(float));
  else
    z = 0;

  error = ex_get_coord(exoid, x, y, z);
  free(x);
  free(y);
  if (num_dim >= 3)
    free(z);

  for (i = 0; i < num_dim; i++) {
    coord_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_coord_names(exoid, coord_names);

  for (i = 0; i < num_dim; i++)
    free(coord_names[i]);

  {
    int num_attrs = 0;
    error         = ex_get_attr_param(exoid, EX_NODAL, 0, &num_attrs);
    if (num_attrs > 0) {
      for (j = 0; j < num_attrs; j++) {
        attrib_names[j] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }
      error = ex_get_attr_names(exoid, EX_NODAL, 0, attrib_names);

      if (error == 0) {
        attrib = (float *)calloc(num_nodes, sizeof(float));
        for (j = 0; j < num_attrs; j++) {
          error = ex_get_one_attr(exoid, EX_NODAL, 0, j + 1, attrib);
          free(attrib_names[j]);
        }
        free(attrib);
      }
    }
  }

  /* read element order map */

  elem_map = (int *)calloc(num_elem, sizeof(int));

  error = ex_get_map(exoid, elem_map);

  free(elem_map);

  /* read element block parameters */
  if (num_elem_blk > 0) {
    ids                = (int *)calloc(num_elem_blk, sizeof(int));
    num_elem_in_block  = (int *)calloc(num_elem_blk, sizeof(int));
    num_nodes_per_elem = (int *)calloc(num_elem_blk, sizeof(int));
    num_attr           = (int *)calloc(num_elem_blk, sizeof(int));

    error = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);

    for (i = 0; i < num_elem_blk; i++) {
      block_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_names(exoid, EX_ELEM_BLOCK, block_names);

    for (i = 0; i < num_elem_blk; i++) {
      ex_get_name(exoid, EX_ELEM_BLOCK, ids[i], name);
      /* 'name' should equal 'block_names[i]' at this point */

      error = ex_get_elem_block(exoid, ids[i], elem_type, &(num_elem_in_block[i]),
                                &(num_nodes_per_elem[i]), &(num_attr[i]));
      free(block_names[i]);
    }

    /* read element block properties */
    num_props = ex_inquire_int(exoid, EX_INQ_EB_PROP);
    for (i = 0; i < num_props; i++) {
      prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_prop_names(exoid, EX_ELEM_BLOCK, prop_names);

    for (i = 1; i < num_props; i++) { /* Prop 1 is id; skip that here */
      for (j = 0; j < num_elem_blk; j++) {
        error = ex_get_prop(exoid, EX_ELEM_BLOCK, ids[j], prop_names[i], &prop_value);
      }
    }

    for (i = 0; i < num_props; i++)
      free(prop_names[i]);
  }

  /* read element connectivity */
  for (i = 0; i < num_elem_blk; i++) {
    if (num_elem_in_block[i] > 0) {
      connect = (int *)calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));

      error = ex_get_conn(exoid, EX_ELEM_BLOCK, ids[i], connect, 0, 0);
      free(connect);
    }
  }

  /* read element block attributes */
  for (i = 0; i < num_elem_blk; i++) {
    if (num_elem_in_block[i] > 0) {
      for (j            = 0; j < num_attr[i]; j++)
        attrib_names[j] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));

      attrib = (float *)calloc(num_attr[i] * num_elem_in_block[i], sizeof(float));
      error  = ex_get_attr(exoid, EX_ELEM_BLOCK, ids[i], attrib);

      if (error == 0) {
        error = ex_get_attr_names(exoid, EX_ELEM_BLOCK, ids[i], attrib_names);
      }
      free(attrib);
      for (j = 0; j < num_attr[i]; j++)
        free(attrib_names[j]);
    }
  }

  if (num_elem_blk > 0) {
    free(ids);
    free(num_nodes_per_elem);
    free(num_attr);
  }

  /* read individual node sets */
  if (num_node_sets > 0) {
    ids = (int *)calloc(num_node_sets, sizeof(int));

    error = ex_get_ids(exoid, EX_NODE_SET, ids);

    for (i = 0; i < num_node_sets; i++) {
      nset_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    /* Get all nodeset names in one call */
    error = ex_get_names(exoid, EX_NODE_SET, nset_names);

    for (i = 0; i < num_node_sets; i++) {
      /* Can also get the names one at a time... */
      ex_get_name(exoid, EX_NODE_SET, ids[i], name);
      /* 'name' should equal 'block_names[i]' at this point */

      error = ex_get_set_param(exoid, EX_NODE_SET, ids[i], &num_nodes_in_set, &num_df_in_set);
      free(nset_names[i]);
      node_list = (int *)calloc(num_nodes_in_set, sizeof(int));
      dist_fact = (float *)calloc(num_nodes_in_set, sizeof(float));

      error = ex_get_set(exoid, EX_NODE_SET, ids[i], node_list, 0);

      if (num_df_in_set > 0) {
        error = ex_get_set_dist_fact(exoid, EX_NODE_SET, ids[i], dist_fact);
      }

      free(node_list);
      free(dist_fact);

      {
        int num_attrs = 0;
        error         = ex_get_attr_param(exoid, EX_NODE_SET, ids[i], &num_attrs);
        if (num_attrs > 0) {
          for (j = 0; j < num_attrs; j++) {
            attrib_names[j] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
          }
          error = ex_get_attr_names(exoid, EX_NODE_SET, ids[i], attrib_names);

          if (error == 0) {
            attrib = (float *)calloc(num_nodes_in_set, sizeof(float));
            for (j = 0; j < num_attrs; j++) {
              error = ex_get_one_attr(exoid, EX_NODE_SET, ids[i], j + 1, attrib);
              free(attrib_names[j]);
            }
            free(attrib);
          }
        }
      }
    }
    free(ids);

    /* read node set properties */
    num_props = ex_inquire_int(exoid, EX_INQ_NS_PROP);

    for (i = 0; i < num_props; i++) {
      prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }
    prop_values = (int *)calloc(num_node_sets, sizeof(int));

    error = ex_get_prop_names(exoid, EX_NODE_SET, prop_names);

    for (i = 0; i < num_props; i++) {
      error = ex_get_prop_array(exoid, EX_NODE_SET, prop_names[i], prop_values);
    }
    for (i = 0; i < num_props; i++)
      free(prop_names[i]);
    free(prop_values);

    /* read concatenated node sets; this produces the same information as
     * the above code which reads individual node sets
     */
    {
      num_node_sets                = ex_inquire_int(exoid, EX_INQ_NODE_SETS);
      ns_specs.sets_ids            = (int *)calloc(num_node_sets, sizeof(int));
      ns_specs.num_entries_per_set = (int *)calloc(num_node_sets, sizeof(int));
      ns_specs.num_dist_per_set    = (int *)calloc(num_node_sets, sizeof(int));
      ns_specs.sets_entry_index    = (int *)calloc(num_node_sets, sizeof(int));
      ns_specs.sets_dist_index     = (int *)calloc(num_node_sets, sizeof(int));

      list_len                 = ex_inquire_int(exoid, EX_INQ_NS_NODE_LEN);
      ns_specs.sets_entry_list = (int *)calloc(list_len, sizeof(int));

      list_len                 = ex_inquire_int(exoid, EX_INQ_NS_DF_LEN);
      ns_specs.sets_dist_fact  = (float *)calloc(list_len, sizeof(float));
      ns_specs.sets_extra_list = NULL;

      error = ex_get_concat_sets(exoid, EX_NODE_SET, &ns_specs);
    }
  }

  /* read individual side sets */
  if (num_side_sets > 0) {
    ids = (int *)calloc(num_side_sets, sizeof(int));

    error = ex_get_ids(exoid, EX_SIDE_SET, ids);
    for (i = 0; i < num_side_sets; i++) {
      sset_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_names(exoid, EX_SIDE_SET, sset_names);

    for (i = 0; i < num_side_sets; i++) {
      ex_get_name(exoid, EX_SIDE_SET, ids[i], name);

      error = ex_get_set_param(exoid, EX_SIDE_SET, ids[i], &num_sides_in_set, &num_df_in_set);
      free(sset_names[i]);

      /* Note: The # of elements is same as # of sides!  */
      num_elem_in_set = num_sides_in_set;
      elem_list       = (int *)calloc(num_elem_in_set, sizeof(int));
      side_list       = (int *)calloc(num_sides_in_set, sizeof(int));
      node_ctr_list   = (int *)calloc(num_elem_in_set, sizeof(int));
      node_list       = (int *)calloc(num_elem_in_set * 21, sizeof(int));
      dist_fact       = (float *)calloc(num_df_in_set, sizeof(float));

      error = ex_get_set(exoid, EX_SIDE_SET, ids[i], elem_list, side_list);
      error = ex_get_side_set_node_list(exoid, ids[i], node_ctr_list, node_list);

      if (num_df_in_set > 0) {
        error = ex_get_set_dist_fact(exoid, EX_SIDE_SET, ids[i], dist_fact);
      }

      free(elem_list);
      free(side_list);
      free(node_ctr_list);
      free(node_list);
      free(dist_fact);
    }

    /* read side set properties */
    num_props = ex_inquire_int(exoid, EX_INQ_SS_PROP);

    for (i = 0; i < num_props; i++) {
      prop_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_prop_names(exoid, EX_SIDE_SET, prop_names);

    for (i = 0; i < num_props; i++) {
      for (j = 0; j < num_side_sets; j++) {
        error = ex_get_prop(exoid, EX_SIDE_SET, ids[j], prop_names[i], &prop_value);
      }
    }
    for (i = 0; i < num_props; i++)
      free(prop_names[i]);
    free(ids);

    num_side_sets = ex_inquire_int(exoid, EX_INQ_SIDE_SETS);

    if (num_side_sets > 0) {
      elem_list_len = ex_inquire_int(exoid, EX_INQ_SS_ELEM_LEN);
      node_list_len = ex_inquire_int(exoid, EX_INQ_SS_NODE_LEN);
      df_list_len   = ex_inquire_int(exoid, EX_INQ_SS_DF_LEN);
    }

    /* read concatenated side sets; this produces the same information as
     * the above code which reads individual side sets
     */

    /* concatenated side set read */

    if (num_side_sets > 0) {
      ss_specs.sets_ids            = (int *)calloc(num_side_sets, sizeof(int));
      ss_specs.num_entries_per_set = (int *)calloc(num_side_sets, sizeof(int));
      ss_specs.num_dist_per_set    = (int *)calloc(num_side_sets, sizeof(int));
      ss_specs.sets_entry_index    = (int *)calloc(num_side_sets, sizeof(int));
      ss_specs.sets_dist_index     = (int *)calloc(num_side_sets, sizeof(int));
      ss_specs.sets_entry_list     = (int *)calloc(elem_list_len, sizeof(int));
      ss_specs.sets_extra_list     = (int *)calloc(elem_list_len, sizeof(int));
      ss_specs.sets_dist_fact      = (float *)calloc(df_list_len, sizeof(float));

      error = ex_get_concat_sets(exoid, EX_SIDE_SET, &ss_specs);
    }
  }
  /* end of concatenated side set read */

  /* read QA records */
  num_qa_rec = ex_inquire_int(exoid, EX_INQ_QA);

  for (i = 0; i < num_qa_rec; i++) {
    for (j = 0; j < 4; j++) {
      qa_record[i][j] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }
  }

  error = ex_get_qa(exoid, qa_record);

  /* read information records */
  num_info = ex_inquire_int(exoid, EX_INQ_INFO);
  for (i = 0; i < num_info; i++) {
    info[i] = (char *)calloc((MAX_LINE_LENGTH + 1), sizeof(char));
  }

  error = ex_get_info(exoid, info);

  for (i = 0; i < num_info; i++) {
    free(info[i]);
  }

  /* read global variables parameters and names */
  error = ex_get_variable_param(exoid, EX_GLOBAL, &num_glo_vars);

  for (i = 0; i < num_glo_vars; i++) {
    var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }

  error = ex_get_variable_names(exoid, EX_GLOBAL, num_glo_vars, var_names);

  for (i = 0; i < num_glo_vars; i++) {
    free(var_names[i]);
  }

  /* read nodal variables parameters and names */
  num_nod_vars = 0;
  if (num_nodes > 0) {
    error = ex_get_variable_param(exoid, EX_NODAL, &num_nod_vars);

    for (i = 0; i < num_nod_vars; i++) {
      var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_variable_names(exoid, EX_NODAL, num_nod_vars, var_names);
    for (i = 0; i < num_nod_vars; i++) {
      free(var_names[i]);
    }
  }

  /* read element variables parameters and names */
  num_ele_vars = 0;
  if (num_elem > 0) {
    error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_ele_vars);

    for (i = 0; i < num_ele_vars; i++) {
      var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    error = ex_get_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
    for (i = 0; i < num_ele_vars; i++) {
      free(var_names[i]);
    }

    /* read element variable truth table */
    if (num_ele_vars > 0) {
      truth_tab = (int *)calloc((num_elem_blk * num_ele_vars), sizeof(int));

      error = ex_get_truth_table(exoid, EX_ELEM_BLOCK, num_elem_blk, num_ele_vars, truth_tab);
      free(truth_tab);
    }
  }

  /* read nodeset variables parameters and names */
  num_nset_vars = 0;
  if (num_node_sets > 0) {
    error = ex_get_variable_param(exoid, EX_NODE_SET, &num_nset_vars);

    if (num_nset_vars > 0) {
      for (i = 0; i < num_nset_vars; i++) {
        var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }

      error = ex_get_variable_names(exoid, EX_NODE_SET, num_nset_vars, var_names);

      for (i = 0; i < num_nset_vars; i++) {
        free(var_names[i]);
      }

      /* read nodeset variable truth table */
      if (num_nset_vars > 0) {
        truth_tab = (int *)calloc((num_node_sets * num_nset_vars), sizeof(int));

        error = ex_get_truth_table(exoid, EX_NODE_SET, num_node_sets, num_nset_vars, truth_tab);
        free(truth_tab);
      }
    }
  }

  /* read sideset variables parameters and names */
  num_sset_vars = 0;
  if (num_side_sets > 0) {
    error = ex_get_variable_param(exoid, EX_SIDE_SET, &num_sset_vars);

    if (num_sset_vars > 0) {
      for (i = 0; i < num_sset_vars; i++) {
        var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }

      error = ex_get_variable_names(exoid, EX_SIDE_SET, num_sset_vars, var_names);

      for (i = 0; i < num_sset_vars; i++) {
        free(var_names[i]);
      }

      /* read sideset variable truth table */
      if (num_sset_vars > 0) {
        truth_tab = (int *)calloc((num_side_sets * num_sset_vars), sizeof(int));

        error = ex_get_truth_table(exoid, EX_SIDE_SET, num_side_sets, num_sset_vars, truth_tab);
        free(truth_tab);
      }
    }
  }

  /* determine how many time steps are stored */
  num_time_steps = ex_inquire_int(exoid, EX_INQ_TIME);

  /* read time value at one time step */
  time_step = 3;
  error     = ex_get_time(exoid, time_step, &time_value);

  /* read time values at all time steps */
  time_values = (float *)calloc(num_time_steps, sizeof(float));

  error = ex_get_all_times(exoid, time_values);

  free(time_values);

  /* read all global variables at one time step */
  var_values = (float *)calloc(num_glo_vars, sizeof(float));

  error = ex_get_var(exoid, time_step, EX_GLOBAL, 1, 0, num_glo_vars, var_values);
  free(var_values);

  /* read a single global variable through time */
  var_index = 1;
  beg_time  = 1;
  end_time  = -1;

  var_values = (float *)calloc(num_time_steps, sizeof(float));

  error = ex_get_var_time(exoid, EX_GLOBAL, var_index, 1, beg_time, end_time, var_values);
  free(var_values);

  /* read a nodal variable at one time step */
  if (num_nodes > 0) {
    var_values = (float *)calloc(num_nodes, sizeof(float));

    error = ex_get_var(exoid, time_step, EX_NODAL, var_index, 1, num_nodes, var_values);
    free(var_values);

    /* read a nodal variable through time */
    var_values = (float *)calloc(num_time_steps, sizeof(float));

    node_num = 1;
    error = ex_get_var_time(exoid, EX_NODAL, var_index, node_num, beg_time, end_time, var_values);
    free(var_values);
  }

  /* read an element variable at one time step */
  if (num_elem_blk > 0) {
    ids = (int *)calloc(num_elem_blk, sizeof(int));

    error = ex_get_elem_blk_ids(exoid, ids);

    for (i = 0; i < num_elem_blk; i++) {
      if (num_elem_in_block[i] > 0) {
        var_values = (float *)calloc(num_elem_in_block[i], sizeof(float));

        error = ex_get_var(exoid, time_step, EX_ELEM_BLOCK, var_index, ids[i], num_elem_in_block[i],
                           var_values);

        free(var_values);
      }
    }
    free(num_elem_in_block);
    free(ids);
  }
  /* read an element variable through time */
  if (num_ele_vars > 0) {
    var_values = (float *)calloc(num_time_steps, sizeof(float));

    var_index = 2;
    elem_num  = 2;
    error =
        ex_get_var_time(exoid, EX_ELEM_BLOCK, var_index, elem_num, beg_time, end_time, var_values);

    free(var_values);
  }

  /* read a sideset variable at one time step */
  if (num_sset_vars > 0) {
    for (i = 0; i < num_side_sets; i++) {
      var_values = (float *)calloc(ss_specs.num_entries_per_set[i], sizeof(float));

      error = ex_get_var(exoid, time_step, EX_SIDE_SET, var_index, ss_specs.sets_ids[i],
                         ss_specs.num_entries_per_set[i], var_values);

      free(var_values);
    }
  }

  /* read a nodeset variable at one time step */
  if (num_nset_vars > 0) {
    for (i = 0; i < num_node_sets; i++) {
      var_values = (float *)calloc(ns_specs.num_entries_per_set[i], sizeof(float));

      error = ex_get_var(exoid, time_step, EX_NODE_SET, var_index, ns_specs.sets_ids[i],
                         ns_specs.num_entries_per_set[i], var_values);

      free(var_values);
    }
  }
  error = ex_close(exoid);
  return 0;
}
