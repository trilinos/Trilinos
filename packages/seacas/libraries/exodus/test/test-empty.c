/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <exodusII.h>
#include <stdio.h>

int main()
{
  float version = 0.0;

  ex_opts(EX_VERBOSE | EX_ABORT);
  int CPU_word_size = 0; /* sizeof(float) */
  int IO_word_size  = 4; /* (4 bytes) */

  /* ======================================== */
  /* Create an empty exodus file             */
  /* ====================================== */
  int exoid = ex_create("test.exo", EX_CLOBBER, &CPU_word_size, &IO_word_size);
  ex_close(exoid);

  /* ======================================== */
  /* Now try to open and read the empty file */
  /* ====================================== */

  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  printf("test.exo exoid = %d\n", exoid);

  {
    char title[MAX_LINE_LENGTH + 1];
    int  num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
    int  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk,
                            &num_node_sets, &num_side_sets);
    printf("after ex_get_init, error = %3d\n", error);
    if (error) {
      exit(-1);
    }
    else {
      exit(0);
    }
  }
}
