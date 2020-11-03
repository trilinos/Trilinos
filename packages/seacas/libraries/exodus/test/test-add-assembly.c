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
    if (error != EX_NOERR) {                                                                       \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

int main(int argc, char **argv)
{
  int  exoid;
  int  num_assembly;
  int  error;
  int *ids;
  int  CPU_word_size;
  int  IO_word_size;
  int  idum;

  float version;
  float fdum;

  char  title_chk[MAX_LINE_LENGTH + 1];
  char *cdum = 0;

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  exoid = ex_open("test-assembly.exo", /* filename path */
                  EX_WRITE,            /* access mode = READ */
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

    num_assembly = par.num_assembly;

    /* Check that ex_inquire gives same title */
    EXCHECK(ex_inquire(exoid, EX_INQ_TITLE, &idum, &fdum, title_chk));
    if (strcmp(par.title, title_chk) != 0) {
      printf("error in ex_inquire for EX_INQ_TITLE %s, vs %s\n", par.title, title_chk);
    }
  }

  /* See if we can add a new assembly to an existing file... */
  ids = (int *)calloc(num_assembly, sizeof(int));
  EXCHECK(ex_get_ids(exoid, EX_ASSEMBLY, ids));
  {
    int64_t     list_222[] = {100, 200, 300, 400};
    ex_assembly assembly   = {222, "NewAssembly", EX_ASSEMBLY, 4, NULL};
    EXCHECK(ex_put_assembly(exoid, assembly));

    assembly.entity_list = list_222;
    EXCHECK(ex_put_assembly(exoid, assembly));
  }

  /*  free(block_names[i]); */
  EXCHECK(ex_close(exoid));
  return 0;
}
