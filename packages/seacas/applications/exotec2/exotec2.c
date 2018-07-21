/*
 * Copyright(C) 2008-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *
 * * Neither the name of NTESS nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  printf("                           *** EXOTEC Version 2.02 ***\n");
  printf("                                Revised 2011/07/01\n\n");
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
