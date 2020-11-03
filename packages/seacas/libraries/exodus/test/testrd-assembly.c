/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testrd - read exodus file test.exo created by testwt
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
 *
 *****************************************************************************/

#include "exodusII.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    error = (funcall);                                                                             \
    printf("after %s, error = %d\n", TOSTRING(funcall), error);                                    \
    if (error != EX_NOERR && error != EX_WARN) {                                                   \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  int  exoid, num_elem_blk, num_assembly;
  int  num_assembly_vars;
  int  error;
  int  i;
  int *ids;
  int *num_elem_in_block  = NULL;
  int *num_nodes_per_elem = NULL;
  int *num_attr           = NULL;
  int  CPU_word_size, IO_word_size;
  int  idum;

  float version;
  float fdum;

  char *var_names[10];
  char *block_names[10];
  char  name[MAX_STR_LENGTH + 1];
  char  elem_type[MAX_STR_LENGTH + 1];
  char  title_chk[MAX_LINE_LENGTH + 1];
  char *cdum = 0;

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  exoid = ex_open("test-assembly.exo", /* filename path */
                  EX_READ,             /* access mode = READ */
                  &CPU_word_size,      /* CPU word size */
                  &IO_word_size,       /* IO word size */
                  &version);           /* ExodusII library version */

  printf("\nafter ex_open\n");
  if (exoid < 0) {
    exit(1);
  }

  printf("test-assembly.exo is an EXODUSII file; version %4.2f\n", version);
  /*   printf ("         CPU word size %1d\n",CPU_word_size);  */
  printf("         I/O word size %1d\n", IO_word_size);
  ex_inquire(exoid, EX_INQ_API_VERS, &idum, &version, cdum);
  printf("EXODUSII API; version %4.2f\n", version);

  ex_inquire(exoid, EX_INQ_LIB_VERS, &idum, &version, cdum);
  printf("EXODUSII Library API; version %4.2f (%d)\n", version, idum);

  /* read database parameters */
  {
    ex_init_params par;
    EXCHECK(ex_get_init_ext(exoid, &par));

    printf("database parameters:\n");
    printf("title =  '%s'\n", par.title);
    printf("num_dim = %" PRId64 "\n", par.num_dim);
    printf("num_assembly = %" PRId64 "\n", par.num_assembly);
    printf("num_nodes = %" PRId64 "\n", par.num_nodes);
    printf("num_edge = %" PRId64 "\n", par.num_edge);
    printf("num_face = %" PRId64 "\n", par.num_face);
    printf("num_elem = %" PRId64 "\n", par.num_elem);
    printf("num_elem_blk = %" PRId64 "\n", par.num_elem_blk);
    printf("num_node_sets = %" PRId64 "\n", par.num_node_sets);
    printf("num_side_sets = %" PRId64 "\n", par.num_side_sets);

    /* Check that ex_inquire gives same title */
    EXCHECK(ex_inquire(exoid, EX_INQ_TITLE, &idum, &fdum, title_chk));
    if (strcmp(par.title, title_chk) != 0) {
      printf("error in ex_inquire for EX_INQ_TITLE %s, vs %s\n", par.title, title_chk);
    }
    num_elem_blk = par.num_elem_blk;
    num_assembly = par.num_assembly;
  }

  /* read element block parameters */
  if (num_elem_blk > 0) {
    ids                = (int *)calloc(num_elem_blk, sizeof(int));
    num_elem_in_block  = (int *)calloc(num_elem_blk, sizeof(int));
    num_nodes_per_elem = (int *)calloc(num_elem_blk, sizeof(int));
    num_attr           = (int *)calloc(num_elem_blk, sizeof(int));

    EXCHECK(ex_get_ids(exoid, EX_ELEM_BLOCK, ids));

    for (i = 0; i < num_elem_blk; i++) {
      block_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    EXCHECK(ex_get_names(exoid, EX_ELEM_BLOCK, block_names));

    for (i = 0; i < num_elem_blk; i++) {
      EXCHECK(ex_get_name(exoid, EX_ELEM_BLOCK, ids[i], name));
      if (strcmp(name, block_names[i]) != 0) {
        printf("error in ex_get_name for block id %d\n", ids[i]);
      }
      EXCHECK(ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &(num_elem_in_block[i]),
                           &(num_nodes_per_elem[i]), 0, 0, &(num_attr[i])));
      printf("element block id = %2d\n", ids[i]);
      printf("element type = '%s'\n", elem_type);
      printf("num_elem_in_block = %2d\n", num_elem_in_block[i]);
      printf("num_nodes_per_elem = %2d\n", num_nodes_per_elem[i]);
      printf("num_attr = %2d\n", num_attr[i]);
      printf("name = '%s'\n", block_names[i]);
    }
  }

  /* Read assembly information */
  /* Verify ex_inquire_int gives same value for assembly count... */
  int chk_num_assembly = ex_inquire_int(exoid, EX_INQ_ASSEMBLY);
  if (chk_num_assembly != num_assembly) {
    printf("error in ex_inquire_int for EX_INQ_ASSEMBLY: %d vs %d\n", chk_num_assembly,
           num_assembly);
  }
  if (num_assembly > 0) {
    int assembly_ids[10];
    ex_get_ids(exoid, EX_ASSEMBLY, assembly_ids);

    char *assembly_names[10];
    for (i = 0; i < num_assembly; i++) {
      assembly_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    ex_assembly assemblies[10];
    int64_t     entity[10];
    for (i = 0; i < num_assembly; i++) {
      assemblies[i].id   = assembly_ids[i];
      assemblies[i].name = assembly_names[i];
      /* Clear out name to make sure still getting same name */
      assemblies[i].name[0] = '\0';

      assemblies[i].entity_list = (int64_t *)&entity;
      EXCHECK(ex_get_assembly(exoid, &assemblies[i]));
      printf("Assembly named '%s' has id %" PRId64 ". It contains %d entities of type '%s'\n\t",
             assemblies[i].name, assemblies[i].id, assemblies[i].entity_count,
             ex_name_of_object(assemblies[i].type));
      for (int j = 0; j < assemblies[i].entity_count; j++) {
        printf("%" PRId64 ", ", entity[j]);
      }
      printf("\n");
    }

    ex_assembly assmbly[10];
    for (i = 0; i < num_assembly; i++) {
      assmbly[i].name = NULL;
      assmbly[i].name = assembly_names[i];
      /* Clear out name to make sure still getting same name */
      assmbly[i].name[0]     = '\0';
      assmbly[i].entity_list = NULL;
    }
    EXCHECK(ex_get_assemblies(exoid, assmbly));
    for (i = 0; i < num_assembly; i++) {
      printf("Assembly named '%s' has id %" PRId64 ". It contains %d entities of type '%s'\n",
             assmbly[i].name, assmbly[i].id, assmbly[i].entity_count,
             ex_name_of_object(assmbly[i].type));
    }

    /* Read attributes... */
    ex_attribute attr[10];

    for (i = 0; i < num_assembly; i++) {
      memset(attr, 0, sizeof(ex_attribute) * 10);
      int att_count = ex_get_attribute_count(exoid, EX_ASSEMBLY, assmbly[i].id);
      printf("Assembly named '%s' with id %" PRId64 ". It contains %d attributes:\n",
             assmbly[i].name, assmbly[i].id, att_count);

      ex_get_attribute_param(exoid, EX_ASSEMBLY, assmbly[i].id, attr);
      ex_get_attributes(exoid, att_count, attr);

      for (int j = 0; j < att_count; j++) {
        printf("\tName: '%s', Type = %d, Value Count = %d\n\t", attr[j].name, attr[j].type,
               (int)attr[j].value_count);
        for (int k = 0; k < attr[j].value_count; k++) {
          if (attr[j].type == EX_INTEGER) {
            int *vals = attr[j].values;
            printf("\t%d", vals[k]);
          }
          else if (attr[j].type == EX_DOUBLE) {
            double *vals = attr[j].values;
            printf("\t%g", vals[k]);
          }
          else if (attr[j].type == EX_CHAR) {
            char *vals = attr[j].values;
            if (vals[k] != '\0') {
              printf("\t%c", vals[k]);
            }
          }
        }
        printf("\n");
        free(attr[j].values);
        attr[j].values = NULL;
      }
      printf("\n");
    }

    { /* Global attributes (includes exodus-internal attributes) */
      ex_attribute attr[10];
      int          att_count = ex_get_attribute_count(exoid, EX_GLOBAL, 0);
      printf("GLOBAL contains %d attributes:\n", att_count);
      for (int j = 0; j < att_count; j++) {
        ex_get_attribute_param(exoid, EX_GLOBAL, 0, attr);
        printf("\tName: '%s', Type = %d, Value Count = %d\n", attr[j].name, attr[j].type,
               (int)attr[j].value_count);
      }
    }

    EXCHECK(ex_get_reduction_variable_param(exoid, EX_ASSEMBLY, &num_assembly_vars));

    if (num_assembly_vars > 0) {
      for (i = 0; i < num_assembly_vars; i++) {
        var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }

      EXCHECK(ex_get_reduction_variable_names(exoid, EX_ASSEMBLY, num_assembly_vars, var_names));

      printf("There are %2d assembly reduction variables; their names are :\n", num_assembly_vars);
      for (i = 0; i < num_assembly_vars; i++) {
        printf(" '%s'\n", var_names[i]);
        free(var_names[i]);
      }

      /* determine how many time steps are stored */
      int num_time_steps = ex_inquire_int(exoid, EX_INQ_TIME);
      printf("There are %2d time steps in the database.\n", num_time_steps);

      float *var_values = (float *)calloc(num_assembly_vars, sizeof(float));
      for (i = 0; i < num_time_steps; i++) {
        float time_value;
        EXCHECK(ex_get_time(exoid, i + 1, &time_value));
        printf("Time at step %d is %f.\n", i + 1, time_value);

        for (int k = 0; k < num_assembly; k++) {
          EXCHECK(ex_get_reduction_vars(exoid, i + 1, EX_ASSEMBLY, assmbly[k].id, num_assembly_vars,
                                        var_values));
          printf("Values for Assembly %" PRId64 " at step %d: %f\t%f\t%f\t%f\n", assmbly[k].id,
                 i + 1, var_values[0], var_values[1], var_values[2], var_values[3]);
        }
      }
      free(var_values);
    }
  }
  /*  free(block_names[i]); */
  EXCHECK(ex_close(exoid));
  return 0;
}
