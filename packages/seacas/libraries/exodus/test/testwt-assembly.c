/*
 * Copyright(C) 1999-2021, 2024 National Technology & Engineering Solutions
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
#define TOSTRING(x)  STRINGIFY(x)

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    int error = (funcall);                                                                         \
    printf("after %s, error = %d\n", TOSTRING(funcall), error);                                    \
    if (error != EX_NOERR) {                                                                       \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  ex_opts(EX_VERBOSE);

  /* Specify compute and i/o word size */
  int CPU_word_size = 8;
  int IO_word_size  = 8;

  /* create EXODUS II file */
  int exoid = ex_create("test-assembly.exo", /* filename path */
                        EX_CLOBBER,          /* create mode */
                        &CPU_word_size,      /* CPU double word size in bytes */
                        &IO_word_size);      /* I/O double word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  int num_elem_blk = 7;
  int num_assembly = 4;
  {
    /* initialize file with parameters */
    int            num_dim       = 3;
    int            num_nodes     = 1;
    int            num_elem      = 7;
    int            num_node_sets = 0;
    int            num_side_sets = 0;
    ex_init_params par           = {.num_dim       = num_dim,
                                    .num_nodes     = num_nodes,
                                    .num_elem      = num_elem,
                                    .num_elem_blk  = num_elem_blk,
                                    .num_node_sets = num_node_sets,
                                    .num_side_sets = num_side_sets,
                                    .num_assembly  = num_assembly};

    char *title = "This is a test";
    ex_copy_string(par.title, title, MAX_LINE_LENGTH + 1);
    EXCHECK(ex_put_init_ext(exoid, &par));
  }

  double coord[] = {0.0};
  EXCHECK(ex_put_coord(exoid, coord, coord, coord));

  /* ======================================================================== */
  /* write element block parameters */
  struct ex_block blocks[7];
  for (int i = 0; i < num_elem_blk; i++) {
    blocks[i] = (ex_block){.type = EX_ELEM_BLOCK, .num_entry = 1, .id = i + 10};
  }

  ex_copy_string(blocks[0].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[1].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[2].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[3].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[4].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[5].topology, "sphere", MAX_STR_LENGTH + 1);
  ex_copy_string(blocks[6].topology, "sphere", MAX_STR_LENGTH + 1);

  blocks[0].num_nodes_per_entry = 1;
  blocks[1].num_nodes_per_entry = 1;
  blocks[2].num_nodes_per_entry = 1;
  blocks[3].num_nodes_per_entry = 1;
  blocks[4].num_nodes_per_entry = 1;
  blocks[5].num_nodes_per_entry = 1;
  blocks[6].num_nodes_per_entry = 1;

  EXCHECK(ex_put_block_params(exoid, num_elem_blk, blocks));

  int connect[] = {1};
  for (int i = 0; i < num_elem_blk; i++) {
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[i].id, connect, NULL, NULL));
  }

  /* Write element block names */
  char *block_names[] = {"block_A", "block_B", "block_C", "block_D",
                         "block_E", "block_F", "block_G"};
  for (int i = 0; i < num_elem_blk; i++) {
    EXCHECK(ex_put_name(exoid, EX_ELEM_BLOCK, blocks[i].id, block_names[i]));
  }

  /* ======================================================================== */
  /* Define and Output Assemblies */
  int64_t list_100[] = {200, 300, 400};
  int64_t list_200[] = {10, 11, 12, 13};
  int64_t list_300[] = {14, 15, 16};
  int64_t list_400[] = {10, 16};

  /* Define assemblies.  Note that first and last do not have an entity list at time of definition.
   */
  ex_assembly assembly[] = {{100, "Root", EX_ASSEMBLY, 3, NULL},
                            {200, "Child2", EX_ELEM_BLOCK, 4, list_200},
                            {300, "Child3", EX_ELEM_BLOCK, 3, list_300},
                            {400, "Child4", EX_ELEM_BLOCK, 2, NULL}};

  for (int i = 0; i < num_assembly; i++) {
    EXCHECK(ex_put_assembly(exoid, assembly[i]));
  }

  /* Now verify that can put entity list at a later time... */
  assembly[0].entity_list = list_100;
  EXCHECK(ex_put_assembly(exoid, assembly[0]));

  assembly[3].entity_list = list_400;
  EXCHECK(ex_put_assembly(exoid, assembly[3]));

  /* ======================================================================== */
  /* Add some arbitrary attributes to the assemblies */
  {
    double scale     = 1.5;
    double offset[]  = {1.1, 2.2, 3.3};
    char  *dimension = "length";
    int    units[]   = {1, 0, 0, -1};

    EXCHECK(ex_put_double_attribute(exoid, EX_ASSEMBLY, 100, "Scale", 1, &scale));
    EXCHECK(ex_put_double_attribute(exoid, EX_ASSEMBLY, 200, "Offset", 3, offset));
    EXCHECK(ex_put_text_attribute(exoid, EX_ASSEMBLY, 300, "Dimension", dimension));
    EXCHECK(ex_put_integer_attribute(exoid, EX_ASSEMBLY, 400, "Units", 4, units));

    ex_attribute attribute = {.entity_type = EX_ASSEMBLY,
                              .entity_id   = 100,
                              .name        = {"Units"},
                              .type        = EX_INTEGER,
                              .value_count = 4,
                              .values      = units};
    EXCHECK(ex_put_attribute(exoid, attribute));

    ex_attribute attr_offset = {EX_ASSEMBLY, 300, "Offset", EX_DOUBLE, 3, offset};
    EXCHECK(ex_put_attribute(exoid, attr_offset));

    /* Make sure this works for non-assemblies also... */
    EXCHECK(ex_put_integer_attribute(exoid, EX_ELEM_BLOCK, 11, "Units", 4, units));
    EXCHECK(ex_put_text_attribute(exoid, EX_GLOBAL, 0, "SOLID_MODEL", "STEP-X-43-1547836-Rev 0"));
  }

  /* ======================================================================== */
  /* Transient Variables */
  int num_assem_vars = 4;
  EXCHECK(ex_put_reduction_variable_param(exoid, EX_ASSEMBLY, num_assem_vars));

  {
    char *var_names[] = {"Momentum_X", "Momentum_Y", "Momentum_Z", "Kinetic_Energy"};
    EXCHECK(ex_put_reduction_variable_names(exoid, EX_ASSEMBLY, num_assem_vars, var_names));
  }

  { /* Output time steps ... */
    double *var_vals = (double *)calloc(num_assem_vars, CPU_word_size);
    for (int ts = 0; ts < 10; ts++) {
      double time_val = (double)(ts + 1) / 100.0;

      EXCHECK(ex_put_time(exoid, ts + 1, &time_val));

      /* write assembly variables */
      for (int k = 0; k < num_assembly; k++) {
        for (int var_idx = 0; var_idx < num_assem_vars; var_idx++) {
          var_vals[var_idx] = (double)(var_idx + 2) * time_val + k;
        }
        EXCHECK(ex_put_reduction_vars(exoid, ts + 1, EX_ASSEMBLY, assembly[k].id, num_assem_vars,
                                      var_vals));
      }
    }
    free(var_vals);
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
 *   - All values for a specific object are output at same time -- Vals for "assembly 100"
 *   - GLOBAL values can be considered reduction variables
 *     - for BW compat, can output either with:
 *       - ex_put_var(exoid, ts, EX_GLOBAL, ?, ?, #gvar, gvar), or
 *       - ex_put_reduction_var(exoid, ts, EX_GLOBAL, ?, #gvar, gvar)
 *     -
 *   - ex_put_reduc_var(exoid, ts, EX_ASSEMBLY, id, #var, var)
 *   - ex_get_reduc_var(exoid, ts, EX_ASSEMBLY, id, #var, &var)
 *
 *   - ex_put_reduc_variable_names()
 *   - ex_put_reduc_variable_param()

 */
/* int ex_put_reduction_var(exoid, time_step, obj_type, obj_id, #values, values); */
