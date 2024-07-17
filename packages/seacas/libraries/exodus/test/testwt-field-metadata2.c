/*
 * Copyright(C) 2022, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <assert.h>
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
  int exoid = ex_create("test-field-metadata.exo", /* filename path */
                        EX_CLOBBER,                /* create mode */
                        &CPU_word_size,            /* CPU double word size in bytes */
                        &IO_word_size);            /* I/O double word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  int num_elem_blk = 3;
  {
    ex_init_params par = {.num_dim       = 3,
                          .num_nodes     = 1,
                          .num_elem      = 7,
                          .num_elem_blk  = num_elem_blk,
                          .num_node_sets = 0,
                          .num_side_sets = 0,
                          .num_assembly  = 0};

    char *title = "This is a test";
    ex_copy_string(par.title, title, MAX_LINE_LENGTH + 1);
    EXCHECK(ex_put_init_ext(exoid, &par));
  }

  double coord[] = {0.0};
  EXCHECK(ex_put_coord(exoid, coord, coord, coord));

  /* ======================================================================== */
  /* write element block parameters */
  struct ex_block blocks[num_elem_blk];
  for (int i = 0; i < num_elem_blk; i++) {
    blocks[i] = (ex_block){.type = EX_ELEM_BLOCK, .num_entry = 1, .id = i + 10};
    ex_copy_string(blocks[i].topology, "sphere", MAX_STR_LENGTH + 1);
    blocks[i].num_nodes_per_entry = 1;
  }

  EXCHECK(ex_put_block_params(exoid, num_elem_blk, blocks));

  int connect[] = {1};
  for (int i = 0; i < num_elem_blk; i++) {
    EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, blocks[i].id, connect, NULL, NULL));
  }

  /* Write element block names */
  for (int i = 0; i < num_elem_blk; i++) {
    char block_names[32];
    sprintf(block_names, "block_%c", i + 'A');
    EXCHECK(ex_put_name(exoid, EX_ELEM_BLOCK, blocks[i].id, block_names));
  }

  {
    int units[] = {1, 0, 0, -1};

    EXCHECK(ex_put_integer_attribute(exoid, EX_ELEM_BLOCK, 11, "Units", 4, units));
    EXCHECK(ex_put_text_attribute(exoid, EX_GLOBAL, 0, "SOLID_MODEL", "STEP-X-43-1547836-Rev 0"));
  }

  /* ======================================================================== */
  /* Transient Variables */
  char *var_names[] = {
      "DispX",        "DispY",        "DispZ",        "Velocity%X",   "Velocity%Y",
      "Velocity%Z",   "Gradient-X$1",
      "Gradient-Y$1", "Gradient-Z$1", "Gradient-X$2", "Gradient-Y$2", "Gradient-Z$2",
      "Gradient-X$3", "Gradient-Y$3", "Gradient-Z$3", "Gradient-X$4", "Gradient-Y$4",
      "Gradient-Z$4", "Gradient-X$5", "Gradient-Y$5", "Gradient-Z$5", "Gradient-X$6",
      "Gradient-Y$6", "Gradient-Z$6", "Gradient-X$7", "Gradient-Y$7", "Gradient-Z$7",
      "Gradient-X$8", "Gradient-Y$8", "Gradient-Z$8", "Gradient-X$9", "Gradient-Y$9", "Gradient-Z$9", "Curl@1",       "Curl@2",
      "Curl@3",       "Curl@4",       "Curl@5",       "Curl@6",       "Curl@7",
      "Curl@8",       "User_h2o",     "User_gas",     "User_ch4",     "User_methane"};
  int num_block_vars = sizeof(var_names) / sizeof(var_names[0]);
  int num_node_vars  = 6;

  EXCHECK(ex_put_variable_param(exoid, EX_ELEM_BLOCK, num_block_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_block_vars, var_names));
  EXCHECK(ex_put_variable_param(exoid, EX_NODAL, num_node_vars));
  EXCHECK(ex_put_variable_names(exoid, EX_NODAL, num_node_vars, var_names));

  int vname = 0;
  {
    struct ex_field field = (ex_field){.entity_type            = EX_ELEM_BLOCK,
                                       .entity_id              = blocks[0].id,
                                       .name                   = "Disp",
                                       .type                   = {EX_VECTOR_3D},
                                       .nesting                = 1,
                                       .component_separator[0] = 0};
    EXCHECK(ex_put_field_metadata(exoid, field));

    /* Put same field on the nodes... */
    field.entity_type = EX_NODAL;
    EXCHECK(ex_put_field_metadata(exoid, field));

    int cardinality =
        field.cardinality[0] != 0 ? field.cardinality[0] : ex_field_cardinality(field.type[0]);
    for (int i = 0; i < cardinality; i++) {
      const char *name = ex_component_field_name(&field, (int[]){i + 1});
      assert(strcmp(var_names[vname++], name) == 0);
    }
  }

  {
    struct ex_field field = (ex_field){.entity_type            = EX_ELEM_BLOCK,
                                       .entity_id              = blocks[0].id,
                                       .name                   = "Velocity",
                                       .type                   = {EX_VECTOR_3D},
                                       .nesting                = 1,
                                       .component_separator[0] = '%'};
    EXCHECK(ex_put_field_metadata(exoid, field));

    /* Put same field on the nodes... */
    field.entity_type = EX_NODAL;
    EXCHECK(ex_put_field_metadata(exoid, field));

    int cardinality =
        field.cardinality[0] != 0 ? field.cardinality[0] : ex_field_cardinality(field.type[0]);
    for (int i = 0; i < cardinality; i++) {
      const char *name = ex_component_field_name(&field, (int[]){i + 1});
      assert(strcmp(var_names[vname++], name) == 0);
    }
  }

  {
    struct ex_field field = (ex_field){.entity_type         = EX_ELEM_BLOCK,
                                       .entity_id           = blocks[1].id,
                                       .name                = "Gradient",
                                       .type                = {EX_VECTOR_3D, EX_FIELD_TYPE_SEQUENCE},
                                       .nesting             = 2,
                                       .cardinality         = {0,9},
                                       .component_separator = {'-', '$'}};
    EXCHECK(ex_put_field_metadata(exoid, field));

    struct ex_field field2 = (ex_field){.entity_type         = EX_ELEM_BLOCK,
                                        .entity_id           = blocks[1].id,
                                        .name                = "Curl",
                                        .type                = {EX_FIELD_TYPE_SEQUENCE},
                                        .nesting             = 1,
					.cardinality         = {9},
                                        .component_separator = {'@'}};
    EXCHECK(ex_put_field_metadata(exoid, field2));
  }

  {
    struct ex_field field = (ex_field){.entity_type         = EX_ELEM_BLOCK,
                                       .entity_id           = blocks[1].id,
                                       .name                = "User",
                                       .type                = {EX_FIELD_TYPE_USER_DEFINED},
                                       .nesting             = 1,
                                       .cardinality         = {4},
                                       .component_separator = {'_'}};
    EXCHECK(ex_put_field_metadata(exoid, field));
    EXCHECK(ex_put_field_suffices(exoid, field, "h2o,gas,ch4,methane"));
  }

  { /* Output time steps ... */
    for (int ts = 0; ts < 1; ts++) {
      double time_val = (double)(ts + 1) / 100.0;

      EXCHECK(ex_put_time(exoid, ts + 1, &time_val));

      /* write variables */
      for (int k = 0; k < num_elem_blk; k++) {
        double *var_vals = (double *)calloc(blocks[k].num_entry, CPU_word_size);
        for (int var_idx = 0; var_idx < num_block_vars; var_idx++) {
          for (int elem = 0; elem < blocks[k].num_entry; elem++) {
            var_vals[elem] = (double)(var_idx + 2) * time_val + elem;
          }
          EXCHECK(ex_put_var(exoid, ts + 1, EX_ELEM_BLOCK, var_idx + 1, blocks[k].id,
                             blocks[k].num_entry, var_vals));
        }
        free(var_vals);
      }
    }
  }

  /* close the EXODUS files
   */
  EXCHECK(ex_close(exoid));
  return 0;
}
