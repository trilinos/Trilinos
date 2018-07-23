/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
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
 *
 */

#include "defs.h"  // for FALSE, TRUE
#include <stdio.h> // for printf, fclose, fprintf, etc

static int input_assign_normal(FILE *finassign, char *inassignname, int nvtxs, int *assignment);
static int input_assign_inv();

int input_assign(FILE *finassign, char *inassignname, int nvtxs, int *assignment)
{
  extern int IN_ASSIGN_INV; /* read assignment in inverse form? */
  int        i;             /* return value */

  if (IN_ASSIGN_INV) {
    i = input_assign_inv(finassign, inassignname, nvtxs, assignment);
  }
  else {
    i = input_assign_normal(finassign, inassignname, nvtxs, assignment);
  }

  return (i);
}

static int input_assign_normal(FILE *finassign, char *inassignname, int nvtxs, int *assignment)
/*
  FILE     *finassign;		input assignment file
  char     *inassignname;		name of input assignment file
  int       nvtxs;		number of vertices to output
  int    *assignment;		values to be printed
*/

{
  extern FILE *Output_File; /* output file or null */
  extern int   CHECK_INPUT; /* print warning messages? */
  extern int   DEBUG_TRACE; /* trace main execution path */
  int          flag;        /* logical conditional */
  int          end_flag;    /* return flag from read_int() */
  int          i, j;        /* loop counter */
  int          read_int();

  if (DEBUG_TRACE > 0) {
    printf("<Entering input_assign>\n");
  }

  /* Get the assignment vector one line at a time, checking as you go. */
  /* First read past any comments at top. */
  end_flag = 1;
  while (end_flag == 1) {
    assignment[0] = read_int(finassign, &end_flag);
  }

  if (assignment[0] < 0) {
    printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n", 1, inassignname,
           assignment[0]);
    return (1);
  }

  if (end_flag == -1) {
    printf("ERROR: No values found in assignment file `%s'\n", inassignname);
    return (1);
  }

  flag = 0;
  if (assignment[0] > nvtxs) {
    flag = assignment[1];
  }
  for (i = 1; i < nvtxs; i++) {
    j = fscanf(finassign, "%d", &(assignment[i]));
    if (j != 1) {
      printf("ERROR: Too few values in assignment file `%s'.\n", inassignname);
      return (1);
    }
    if (assignment[i] < 0) {
      printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n", i + 1, inassignname,
             assignment[i]);
      return (1);
    }
    if (assignment[i] > nvtxs) { /* warn since probably an error */
      if (assignment[i] > flag) {
        flag = assignment[i];
      }
    }
  }

  if (flag && CHECK_INPUT) {
    printf("WARNING: Possible error in assignment file `%s'\n", inassignname);
    printf("         More assignment sets (%d) than vertices (%d)\n", flag, nvtxs);
    if (Output_File != NULL) {
      fprintf(Output_File, "WARNING: Possible error in assignment file `%s'\n", inassignname);
      fprintf(Output_File, "         More assignment sets (%d) than vertices (%d)\n", flag, nvtxs);
    }
  }

  /* Check for spurious extra stuff in file. */
  flag     = FALSE;
  end_flag = 0;
  while (!flag && end_flag != -1) {
    read_int(finassign, &end_flag);
    if (!end_flag) {
      flag = TRUE;
    }
  }
  if (flag && CHECK_INPUT) {
    printf("WARNING: Possible error in assignment file `%s'\n", inassignname);
    printf("         Numerical data found after expected end of file\n");
    if (Output_File != NULL) {
      fprintf(Output_File, "WARNING: Possible error in assignment file `%s'\n", inassignname);
      fprintf(Output_File, "         Numerical data found after expected end of file\n");
    }
  }
  return (0);
}

static int input_assign_inv(FILE *finassign,    /* input assignment file */
                            char *inassignname, /* name of input assignment file */
                            int   nvtxs,        /* number of vertices to output */
                            int * assignment    /* values to be printed */
)
{
  extern int DEBUG_TRACE; /* trace main execution path */
  int        set;         /* set number being read */
  int        size;        /* number of vertices in set */
  int        total;       /* total number of vertices read */
  int        done;        /* have I hit end of file yet? */
  int        end_flag;    /* return flag from read_int() */
  int        i, j, k;     /* loop counter */
  int        read_int();

  if (DEBUG_TRACE > 0) {
    printf("<Entering input_assign_inv>\n");
  }

  /* Get the assignment vector one line at a time, checking as you go. */

  /* Initialize assignment to help error detection. */
  for (i = 0; i < nvtxs; i++) {
    assignment[i] = -1;
  }

  /* First read past any comments at top. */
  total    = 0;
  set      = 0;
  end_flag = 1;
  while (end_flag == 1) {
    size = read_int(finassign, &end_flag);
  }

  if (end_flag == -1) {
    printf("ERROR: In assignment file `%s'\n", inassignname);
    printf("       No values found\n");
    return (1);
  }

  if (size < 0) {
    printf("ERROR: In assignment file `%s'\n", inassignname);
    printf("       Size of set %d less than zero (%d)\n", set, size);
    return (1);
  }

  if (total + size > nvtxs) {
    printf("ERROR: In assignment file `%s'\n", inassignname);
    printf("       Total set sizes greater than nvtxs (%d)\n", nvtxs);
    return (1);
  }

  done = FALSE;
  while (!done && total < nvtxs) {
    for (i = 1; i <= size; i++) {
      j = fscanf(finassign, "%d", &k);
      if (j != 1) {
        printf("ERROR: Too few values in assignment file `%s'.\n", inassignname);
        return (1);
      }

      if (k <= 0 || k > nvtxs) {
        printf("ERROR: In assignment file `%s'\n", inassignname);
        printf("       Entry %d of set %d invalid (%d)\n", total + i, set, k);
        return (1);
      }

      if (assignment[k - 1] != -1) {
        printf("ERROR: In assignment file `%s'\n", inassignname);
        printf("       Vertex %d assigned to multiple sets\n", k);
        return (1);
      }

      assignment[k - 1] = set;
    }

    total += size;
    j = fscanf(finassign, "%d", &size);
    ++set;
    if (j != 1) {
      if (total != nvtxs) {
        printf("ERROR: Too few values in assignment file `%s'.\n", inassignname);
        return (1);
      }

      done = TRUE;
      size = 0;
    }

    if (size < 0) {
      printf("ERROR: In assignment file `%s'\n", inassignname);
      printf("       Size of set %d less than zero (%d)\n", set, size);
      return (1);
    }

    if (total + size > nvtxs) {
      printf("ERROR: In assignment file `%s'\n", inassignname);
      printf("       Total set sizes greater than nvtxs (%d)\n", nvtxs);
      return (1);
    }
  }
  return (0);
}
