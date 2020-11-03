/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

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
  int exoid;
  int num_blob;
  int num_red_vars;
  int num_vars;
  int error;
  int i;
  int CPU_word_size;
  int IO_word_size;
  int idum;

  float version;

  char *var_names[10];
  char *cdum = 0;

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  exoid = ex_open("test-blob.exo", /* filename path */
                  EX_READ,         /* access mode = READ */
                  &CPU_word_size,  /* CPU word size */
                  &IO_word_size,   /* IO word size */
                  &version);       /* ExodusII library version */

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
    printf("num_blobs = %" PRId64 "\n", par.num_blob);
    printf("num_nodes = %" PRId64 "\n", par.num_nodes);
    printf("num_edge = %" PRId64 "\n", par.num_edge);
    printf("num_face = %" PRId64 "\n", par.num_face);
    printf("num_elem = %" PRId64 "\n", par.num_elem);
    printf("num_elem_blk = %" PRId64 "\n", par.num_elem_blk);
    printf("num_node_sets = %" PRId64 "\n", par.num_node_sets);
    printf("num_side_sets = %" PRId64 "\n", par.num_side_sets);

    num_blob = par.num_blob;
  }

  /* Read blob information */
  /* Verify ex_inquire_int gives same value for blob count... */
  int chk_num_blob = ex_inquire_int(exoid, EX_INQ_BLOB);
  if (chk_num_blob != num_blob) {
    printf("error in ex_inquire_int for EX_INQ_BLOB: %d vs %d\n", chk_num_blob, num_blob);
  }
  if (num_blob > 0) {
    int blob_ids[10];
    ex_get_ids(exoid, EX_BLOB, blob_ids);

    char *blob_names[10];
    for (i = 0; i < num_blob; i++) {
      blob_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }

    ex_blob blobs[10];
    for (i = 0; i < num_blob; i++) {
      blobs[i].id   = blob_ids[i];
      blobs[i].name = blob_names[i];
      /* Clear out name to make sure still getting same name */
      blobs[i].name[0] = '\0';

      EXCHECK(ex_get_blob(exoid, &blobs[i]));
      printf("Blob named '%s' has id %" PRId64 ". It contains %" PRId64 " entries.\n\t",
             blobs[i].name, blobs[i].id, blobs[i].num_entry);
      printf("\n");
    }

    ex_blob blb[10];
    for (i = 0; i < num_blob; i++) {
      blb[i].name = NULL;
      blb[i].name = blob_names[i];
      /* Clear out name to make sure still getting same name */
      blb[i].name[0] = '\0';
    }
    EXCHECK(ex_get_blobs(exoid, blb));
    for (i = 0; i < num_blob; i++) {
      printf("Blob named '%s' has id %" PRId64 ". It contains %" PRId64 " entries.\n", blb[i].name,
             blb[i].id, blb[i].num_entry);
    }

    /* Read attributes... */
    ex_attribute attr[10];

    for (i = 0; i < num_blob; i++) {
      memset(attr, 0, sizeof(ex_attribute) * 10);
      int att_count = ex_get_attribute_count(exoid, EX_BLOB, blb[i].id);
      printf("Blob named '%s' with id %" PRId64 ". It contains %d attributes:\n", blb[i].name,
             blb[i].id, att_count);

      ex_get_attribute_param(exoid, EX_BLOB, blb[i].id, attr);
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

    EXCHECK(ex_get_reduction_variable_param(exoid, EX_BLOB, &num_red_vars));
    EXCHECK(ex_get_variable_param(exoid, EX_BLOB, &num_vars));

    if (num_red_vars > 0) {
      for (i = 0; i < num_red_vars; i++) {
        var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }

      EXCHECK(ex_get_reduction_variable_names(exoid, EX_BLOB, num_red_vars, var_names));

      printf("There are %2d blob reduction variables; their names are :\n", num_red_vars);
      for (i = 0; i < num_red_vars; i++) {
        printf(" '%s'\n", var_names[i]);
        free(var_names[i]);
      }
    }

    if (num_vars > 0) {
      for (i = 0; i < num_vars; i++) {
        var_names[i] = (char *)calloc((MAX_STR_LENGTH + 1), sizeof(char));
      }

      EXCHECK(ex_get_variable_names(exoid, EX_BLOB, num_vars, var_names));

      printf("There are %2d blob variables; their names are :\n", num_vars);
      for (i = 0; i < num_vars; i++) {
        printf(" '%s'\n", var_names[i]);
        free(var_names[i]);
      }
    }

    /* determine how many time steps are stored */
    int num_time_steps = ex_inquire_int(exoid, EX_INQ_TIME);
    printf("There are %2d time steps in the database.\n", num_time_steps);

    int64_t max_count = 0;
    for (int k = 0; k < num_blob; k++) {
      if (blobs[k].num_entry > max_count) {
        max_count = blobs[k].num_entry;
      }
    }

    float *vals       = (float *)calloc(max_count, CPU_word_size);
    float *var_values = (float *)calloc(num_red_vars, sizeof(float));
    for (i = 0; i < num_time_steps; i++) {
      float time_value;
      EXCHECK(ex_get_time(exoid, i + 1, &time_value));
      printf("Time at step %d is %f.\n", i + 1, time_value);

      for (int k = 0; k < num_blob; k++) {
        EXCHECK(ex_get_reduction_vars(exoid, i + 1, EX_BLOB, blb[k].id, num_red_vars, var_values));
        printf("Values for Blob %" PRId64 " at step %d: %f\t%f\t%f\t%f\n", blb[k].id, i + 1,
               var_values[0], var_values[1], var_values[2], var_values[3]);

        for (int var_idx = 0; var_idx < num_vars; var_idx++) {
          EXCHECK(ex_get_var(exoid, i + 1, EX_BLOB, var_idx + 1, blobs[k].id, blobs[k].num_entry,
                             vals));
          for (int j = 0; j < blobs[k].num_entry; j++) {
            printf("%5.3f\n", vals[j]);
          }
        }
      }
    }
    free(var_values);
    free(vals);
  }
  EXCHECK(ex_close(exoid));
  return 0;
}
