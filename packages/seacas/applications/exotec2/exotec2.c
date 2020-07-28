/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "add_to_log.h" // for add_to_log
#include <exodusII.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern void tec(int iexi, const char *filename);

static void banner()
{
  time_t     time_val;
  struct tm *time_structure;
  char *     asc_time;

  time_val       = time((time_t *)NULL);
  time_structure = localtime(&time_val);
  asc_time       = asctime(time_structure);

  printf("                           *** EXOTEC Version 2.03 ***\n");
  printf("                                Revised 2017/08/10\n\n");
  printf("                          EXODUS --> TECPLOT TRANSLATOR\n\n");
  printf("                         Run on %s\n\n", asc_time);
}

int main(int argc, char *argv[])
{
  int   exoid;
  char *filename;
  int   cpu_word_size = 8;
  int   io_word_size  = 0;
  float version       = 0.0;

  banner();

  if (argc != 3) {
    fprintf(stderr, "ERROR: Usage is exotec2 exo_in tec_out\n\n");
    exit(1);
  }

  /* Open the files... */
  filename = argv[1];
  exoid    = ex_open(filename, EX_READ, &cpu_word_size, &io_word_size, &version);

  if (exoid < 0) {
    fprintf(stderr, "Cannot open file '%s' - exiting.\n", filename);
    exit(1);
  }

  /* Write the tec file... */
  filename = argv[2];

  tec(exoid, filename);

  ex_close(exoid);

  add_to_log(argv[0], 0.0);
}
