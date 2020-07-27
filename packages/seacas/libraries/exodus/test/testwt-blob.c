/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exodusII.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    error = (funcall);                                                                             \
    printf("after %s, error = %d\n", TOSTRING(funcall), error);                                    \
    if (error != EX_NOERR) {                                                                       \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int num_node_sets, num_side_sets, num_assembly, num_blob;
  int error;
  int CPU_word_size, IO_word_size;

  char *title = "This is a test";

  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */

  CPU_word_size = 8;
  IO_word_size  = 8;

  /* create EXODUS II file */

  exoid = ex_create("test-blob.exo", /* filename path */
                    EX_CLOBBER,      /* create mode */
                    &CPU_word_size,  /* CPU double word size in bytes */
                    &IO_word_size);  /* I/O double word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  /* initialize file with parameters */
  {
    num_dim       = 3;
    num_nodes     = 0;
    num_elem      = 0;
    num_elem_blk  = 0;
    num_node_sets = 0;
    num_side_sets = 0;
    num_assembly  = 0;
    num_blob      = 3;

    ex_init_params par = {.num_dim       = num_dim,
                          .num_nodes     = num_nodes,
                          .num_elem      = num_elem,
                          .num_elem_blk  = num_elem_blk,
                          .num_node_sets = num_node_sets,
                          .num_side_sets = num_side_sets,
                          .num_assembly  = num_assembly,
                          .num_blob      = num_blob};

    ex_copy_string(par.title, title, MAX_LINE_LENGTH + 1);
    EXCHECK(ex_put_init_ext(exoid, &par));
  }

  /* ======================================================================== */
  /* Define and Output Blobs */
  /* Define blobs.
   */
  ex_blob blob[] = {{100, "Tempus", 10}, // ID, Name, Size
                    {200, "IOSS", 20},
                    {300, "Solver", 15}};

  for (int i = 0; i < num_blob; i++) {
    EXCHECK(ex_put_blob(exoid, blob[i]));
  }

  /* ======================================================================== */
  /* Add some arbitrary attributes to the blobs */
  {
    double scale     = 1.5;
    double offset[]  = {1.1, 2.2, 3.3};
    char * dimension = "length";
    int    units[]   = {1, 0, 0, -1};

    EXCHECK(ex_put_double_attribute(exoid, EX_BLOB, 100, "Scale", 1, &scale));
    EXCHECK(ex_put_double_attribute(exoid, EX_BLOB, 200, "Offset", 3, offset));
    EXCHECK(ex_put_text_attribute(exoid, EX_BLOB, 300, "Dimension", dimension));

    ex_attribute attribute = {.entity_type = EX_BLOB,
                              .entity_id   = 100,
                              .name        = {"Units"},
                              .type        = EX_INTEGER,
                              .value_count = 4,
                              .values      = units};
    EXCHECK(ex_put_attribute(exoid, attribute));

    ex_attribute attr_offset = {EX_BLOB, 300, "Offset", EX_DOUBLE, 3, offset};
    EXCHECK(ex_put_attribute(exoid, attr_offset));

    EXCHECK(ex_put_text_attribute(exoid, EX_GLOBAL, 0, "SOLID_MODEL", "STEP-X-43-1547836-Rev 0"));
  }

  /* ======================================================================== */
  /* Transient and Reduction Variables */
  int num_red_blob_vars = 4;
  int num_blob_vars     = 3;
  EXCHECK(ex_put_reduction_variable_param(exoid, EX_BLOB, num_red_blob_vars));
  EXCHECK(ex_put_variable_param(exoid, EX_BLOB, num_blob_vars));

  {
    char *red_var_names[] = {"Momentum_X", "Momentum_Y", "Momentum_Z", "Kinetic_Energy"};
    EXCHECK(ex_put_reduction_variable_names(exoid, EX_BLOB, num_red_blob_vars, red_var_names));
  }

  {
    char *var_names[] = {"X", "XDOT", "XDDOT"};
    EXCHECK(ex_put_variable_names(exoid, EX_BLOB, num_blob_vars, var_names));
  }

  { /* Output time steps ... */
    int64_t max_count = 0;
    for (int k = 0; k < num_blob; k++) {
      if (blob[k].num_entry > max_count) {
        max_count = blob[k].num_entry;
      }
    }

    double *var_vals = (double *)calloc(num_red_blob_vars, CPU_word_size);
    double *vals     = (double *)calloc(max_count, CPU_word_size);
    for (int ts = 0; ts < 4; ts++) {
      double time_val = (double)(ts + 1) / 100.0f;

      EXCHECK(ex_put_time(exoid, ts + 1, &time_val));

      /* write blob reduction variables */
      for (int k = 0; k < num_blob; k++) {
        for (int var_idx = 0; var_idx < num_red_blob_vars; var_idx++) {
          var_vals[var_idx] = (double)(var_idx + 2) * time_val + k;
        }
        EXCHECK(
            ex_put_reduction_vars(exoid, ts + 1, EX_BLOB, blob[k].id, num_red_blob_vars, var_vals));
      }

      /* write blob variables */
      for (int k = 0; k < num_blob; k++) {
        for (int var_idx = 0; var_idx < num_blob_vars; var_idx++) {
          for (int j = 0; j < blob[k].num_entry; j++) {
            vals[j] = (double)(var_idx + 2) * time_val + j + k / 100.0;
          }
          EXCHECK(
              ex_put_var(exoid, ts + 1, EX_BLOB, var_idx + 1, blob[k].id, blob[k].num_entry, vals));
        }
      }
    }
    free(var_vals);
    free(vals);
  }

  /* close the EXODUS files
   */
  EXCHECK(ex_close(exoid));
  return 0;
}

/*!  Reduction Variables:
 *   - Transient values applied to the object instead of each member of the object.
 *     - for example, a value for the element block; not for each element in the block.
 *   - Each object type has the same variables
 *   - All values for a specific object are output at same time -- Vals for "blob 100"
 *   - GLOBAL values can be considered reduction variables
 *     - for BW compat, can output either with:
 *       - ex_put_var(exoid, ts, EX_GLOBAL, ?, ?, #gvar, gvar), or
 *       - ex_put_reduction_var(exoid, ts, EX_GLOBAL, ?, #gvar, gvar)
 *     -
 *   - ex_put_reduc_var(exoid, ts, EX_BLOB, id, #var, var)
 *   - ex_get_reduc_var(exoid, ts, EX_BLOB, id, #var, &var)
 *
 *   - ex_put_reduc_variable_names()
 *   - ex_put_reduc_variable_param()

 */
/* int ex_put_reduction_var(exoid, time_step, obj_type, obj_id, #values, values); */
