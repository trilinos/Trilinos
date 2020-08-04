/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "exodusII.h"
#include "exodusII_int.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    if ((error = (funcall)) != NC_NOERR) {                                                         \
      fprintf(stderr, "ERROR Calling %s, error = %d\n", TOSTRING(funcall), error);                 \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main()
{
  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */

  int CPU_word_size = 8; /* sizeof(float) */
  int IO_word_size  = 8; /* (4 bytes) */
  int error;

  /* create EXODUS II file */

  int exoid = ex_create("test.exo",     /* filename path */
                        EX_CLOBBER,     /* create mode */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size); /* I/O float word size in bytes */

  /* initialize file with parameters */
  constexpr int num_dim       = 3;
  constexpr int num_nodes     = 64;
  constexpr int num_elem      = 9;
  constexpr int num_elem_blk  = 3;
  constexpr int num_node_sets = 0;
  constexpr int num_side_sets = 0;
  std::string   title         = "This is a test";

  EXCHECK(ex_put_init(exoid, title.c_str(), num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets));

  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  std::vector<double> x(num_nodes);
  std::vector<double> y(num_nodes);
  std::vector<double> z(num_nodes);

  x[0] = 0.00;
  x[1] = 0.75;
  x[2] = 1.50;
  x[3] = 2.25;
  x[4] = 3.00;

  x[5] = x[0];
  x[6] = 1.125;
  x[7] = 1.625;
  x[8] = 3.000;

  x[9]  = x[0];
  x[10] = 0.50;
  x[11] = 1.00;
  x[12] = x[2];
  x[13] = 2.00;
  x[14] = 2.50;
  x[15] = 3.00;

  x[16] = x[10];
  x[17] = x[11];
  x[18] = x[12];
  x[19] = x[13];
  x[20] = x[14];

  x[21] = x[11];
  x[22] = x[13];
  x[23] = 0.0;
  x[24] = x[6];
  x[25] = x[7];
  x[26] = 3.0;

  x[27] = x[0];
  x[28] = x[1];
  x[29] = x[2];
  x[30] = x[3];
  x[31] = x[4];

  y[0] = y[1] = y[2] = y[3] = y[4] = 0.00;
  y[5]                             = 0.75;
  y[6]                             = 0.50;
  y[7]                             = y[6];
  y[8]                             = y[5];

  y[9]  = 1.50;
  y[10] = 1.25;
  y[11] = y[12] = y[13] = 1.00;
  y[14]                 = y[10];
  y[15]                 = y[9];

  y[16] = 1.75;
  y[17] = y[18] = y[19] = 2.00;
  y[20]                 = y[16];

  y[21] = y[22] = y[9];
  y[23] = y[26] = 2.25;
  y[24] = y[25] = 2.50;

  y[27] = y[28] = y[29] = y[30] = y[31] = 3.0;

  for (int i = 0; i < num_nodes / 2; i++) {
    z[i] = 0.0;

    x[i + num_nodes / 2] = x[i];
    y[i + num_nodes / 2] = y[i];
    z[i + num_nodes / 2] = 1.0;
  }

  EXCHECK(ex_put_coord(exoid, x.data(), y.data(), z.data()));

  const char *coord_names[3];
  coord_names[0] = "xcoor";
  coord_names[1] = "ycoor";
  coord_names[2] = "zcoor";

  EXCHECK(ex_put_coord_names(exoid, (char **)coord_names));

  std::vector<ex_block> blocks(num_elem_blk);

  blocks[0].id                  = 10;
  blocks[0].type                = EX_ELEM_BLOCK;
  blocks[0].num_entry           = 4;
  blocks[0].num_nodes_per_entry = 16;

  blocks[1].id                  = 20;
  blocks[1].type                = EX_ELEM_BLOCK;
  blocks[1].num_entry           = 4;
  blocks[1].num_nodes_per_entry = 12;

  blocks[2].id                  = 30;
  blocks[2].type                = EX_ELEM_BLOCK;
  blocks[2].num_entry           = 1;
  blocks[2].num_nodes_per_entry = 16;

  ex_copy_string(blocks[0].topology, "hex", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[1].topology, "wedge", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[2].topology, "hex", MAX_STR_LENGTH + 1);

  const char *block_names[3];
  block_names[0] = "block_10";
  block_names[1] = "block_20";
  block_names[2] = "block_30";

  EXCHECK(ex_put_block_params(exoid, num_elem_blk, blocks.data()));

  /* Write element block names */
  for (int i = 0; i < num_elem_blk; i++) {
    EXCHECK(ex_put_name(exoid, EX_ELEM_BLOCK, blocks[i].id, block_names[i]));
  }

  /* write element connectivity */
  int connect0[64] = {1,  3,  12, 10, 2,  7,  11, 6,  3,  5,  16, 14, 4,  9,  15, 8,
                      10, 18, 30, 28, 17, 25, 29, 24, 20, 16, 32, 30, 21, 27, 31, 26};
  int connect1[48] = {3,  14, 12, 8,  13, 7,  14, 16, 20, 15, 21, 23,
                      18, 20, 30, 19, 26, 25, 10, 12, 18, 11, 22, 17};
  int connect2[16] = {12, 14, 20, 18, 13, 23, 19, 22};

  for (int i = 0; i < blocks[0].num_entry * blocks[0].num_nodes_per_entry / 2; i++) {
    connect0[i + 32] = connect0[i] + 32;
  }

  for (int i = 0; i < blocks[1].num_entry * blocks[1].num_nodes_per_entry / 2; i++) {
    connect1[i + 24] = connect1[i] + 32;
  }

  for (int i = 0; i < blocks[2].num_entry * blocks[2].num_nodes_per_entry / 2; i++) {
    connect2[i + 8] = connect2[i] + 32;
  }

  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[0].id, connect0, NULL, NULL));
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[1].id, connect1, NULL, NULL));
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[2].id, connect2, NULL, NULL));

  /* write information records; test empty and just blank-filled records */
  constexpr int num_info = 3;
  const char *  info[3];
  info[0] = "This is the first information record.";
  info[1] = "";
  info[2] = "                                     ";

  EXCHECK(ex_put_info(exoid, num_info, (char **)info));

  /* write results variables parameters and names */
  const int num_glo_vars = 1;

  const char *var_names[4];
  var_names[0] = "glo_vars";

  EXCHECK(ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, (char **)var_names));

  const int num_nod_vars = 2;
  /*              12345678901234567890123456789012 */
  var_names[0] = "node_variable_a_very_long_name_0";
  var_names[1] = "nod_var1";

  EXCHECK(ex_put_variable_param(exoid, EX_NODAL, num_nod_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, (char **)var_names));

  const int num_ele_vars = 3;
  /*              0        1         2         3   */
  /*              12345678901234567890123456789012 */
  var_names[0] = "this_variable_name_is_short";
  var_names[1] = "this_variable_name_is_just_right";
  var_names[2] = "this_variable_name_is_tooooo_long";

  EXCHECK(ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_ele_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_ele_vars, (char **)var_names));

  // for each time step, write the analysis results;
  // the code below fills the arrays glob_var_vals,
  // nodal_var_vals, and elem_var_vals with values for debugging purposes;
  int       whole_time_step = 1;
  const int num_time_steps  = 10;

  std::vector<double> glob_var_vals(num_glo_vars);
  std::vector<double> nodal_var_vals(num_nodes);
  std::vector<double> elem_var_vals(4);

  for (int i = 0; i < num_time_steps; i++) {
    double time_value = (double)(i + 1) / 100.;

    /* write time value */
    EXCHECK(ex_put_time(exoid, whole_time_step, &time_value));

    // write global variables
    for (int j = 0; j < num_glo_vars; j++) {
      glob_var_vals[j] = (double)(j + 2) * time_value;
    }

    EXCHECK(
        ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, glob_var_vals.data()));

    // write nodal variables
    for (int k = 1; k <= num_nod_vars; k++) {
      for (int j = 0; j < num_nodes; j++) {
        nodal_var_vals[j] = (double)k + ((double)(j + 1) * time_value);
      }
      EXCHECK(ex_put_var(exoid, whole_time_step, EX_NODAL, k, 1, num_nodes, nodal_var_vals.data()));
    }

    // write element variables
    for (int k = 1; k <= num_ele_vars; k++) {
      for (int j = 0; j < num_elem_blk; j++) {
        for (int m = 0; m < blocks[j].num_entry; m++) {
          elem_var_vals[m] = (double)(k + 1) + (double)(j + 2) + ((double)(m + 1) * time_value);
        }
        EXCHECK(ex_put_var(exoid, whole_time_step, EX_ELEM_BLOCK, k, blocks[j].id,
                           blocks[j].num_entry, elem_var_vals.data()));
      }
    }

    whole_time_step++;

    // update the data file; this should be done at the end of every time step
    // to ensure that no data is lost if the analysis dies
    EXCHECK(ex_update(exoid));
  }

  // close the EXODUS files
  EXCHECK(ex_close(exoid));
  return 0;
}
