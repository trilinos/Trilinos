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

#include "defs.h"
#include "params.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

static char line[LINE_LENGTH];       /* space to hold values */
static int  offset    = 0;           /* offset into line for next data */
static int  break_pnt = LINE_LENGTH; /* place in sequence to pause */
static int  save_pnt;                /* place in sequence to save */

static void flush_line();

double read_val(FILE *infile,  /* file to read value from */
                int * end_flag /* 0 => OK, 1 => EOL, -1 => EOF */
)
{
  double val;         /* return value */
  char * ptr;         /* ptr to next string to read */
  char * ptr2;        /* ptr to next string to read */
  int    length;      /* length of line to read */
  int    length_left; /* length of line still around */
  int    white_seen;  /* have I detected white space yet? */
  int    done;        /* checking for end of scan */
  int    i;           /* loop counter */
  double strtod();

  *end_flag = 0;

  if (offset == 0 || offset >= break_pnt) {
    if (offset >= break_pnt) { /* Copy rest of line back to beginning. */
      length_left = LINE_LENGTH - save_pnt - 1;
      ptr2        = line;
      ptr         = &line[save_pnt];
      for (i = length_left; i; i--) {
        *ptr2++ = *ptr++;
      }
      length = save_pnt + 1;
    }
    else {
      length      = LINE_LENGTH;
      length_left = 0;
    }

    line[LINE_LENGTH - 1] = ' ';
    line[LINE_LENGTH - 2] = ' ';
    /* Now read next line, or next segment of current one. */
    ptr2 = fgets(&line[length_left], length, infile);

    if (ptr2 == NULL) { /* We've hit end of file. */
      *end_flag = -1;
      return (0.0);
    }

    if (line[LINE_LENGTH - 1] == '\0' && line[LINE_LENGTH - 2] != '\0' &&
        line[LINE_LENGTH - 2] != '\n' && line[LINE_LENGTH - 2] != '\f') {
      /* Line too long.  Find last safe place in line. */
      break_pnt  = LINE_LENGTH - 1;
      save_pnt   = break_pnt;
      white_seen = FALSE;
      done       = FALSE;
      while (!done) {
        --break_pnt;
        if (line[break_pnt] != '\0') {
          if (isspace(line[break_pnt])) {
            if (!white_seen) {
              save_pnt   = break_pnt + 1;
              white_seen = TRUE;
            }
          }
          else if (white_seen) {
            done = TRUE;
          }
        }
      }
    }
    else {
      break_pnt = LINE_LENGTH;
    }

    offset = 0;
  }

  while (offset < LINE_LENGTH && isspace(line[offset])) {
    offset++;
  }
  if (offset == LINE_LENGTH || line[offset] == '%' || line[offset] == '#') {
    *end_flag = 1;
    if (break_pnt < LINE_LENGTH) {
      flush_line(infile);
    }
    return (0.0);
  }

  ptr = &(line[offset]);
  val = strtod(ptr, &ptr2);

  if (ptr2 == ptr) { /* End of input line. */
    offset    = 0;
    *end_flag = 1;
    return (0.0);
  }

  offset = (int)(ptr2 - line) / sizeof(char);

  return (val);
}

int read_int(FILE *infile,  /* file to read value from */
             int * end_flag /* 0 => OK, 1 => EOL, -1 => EOF */
)
{
  int   val;         /* return value */
  char *ptr;         /* ptr to next string to read */
  char *ptr2;        /* ptr to next string to read */
  int   length;      /* length of line to read */
  int   length_left; /* length of line still around */
  int   white_seen;  /* have I detected white space yet? */
  int   done;        /* checking for end of scan */
  int   i;           /* loop counter */

  *end_flag = 0;

  if (offset == 0 || offset >= break_pnt) {
    if (offset >= break_pnt) { /* Copy rest of line back to beginning. */
      length_left = LINE_LENGTH - save_pnt - 1;
      ptr2        = line;
      ptr         = &line[save_pnt];
      for (i = length_left; i; i--) {
        *ptr2++ = *ptr++;
      }
      length = save_pnt + 1;
    }
    else {
      length      = LINE_LENGTH;
      length_left = 0;
    }

    line[LINE_LENGTH - 1] = ' ';
    line[LINE_LENGTH - 2] = ' ';
    /* Now read next line, or next segment of current one. */
    ptr2 = fgets(&line[length_left], length, infile);

    if (ptr2 == NULL) { /* We've hit end of file. */
      *end_flag = -1;
      return (0);
    }

    if (line[LINE_LENGTH - 1] == '\0' && line[LINE_LENGTH - 2] != '\0' &&
        line[LINE_LENGTH - 2] != '\n' && line[LINE_LENGTH - 2] != '\f') {
      /* Line too long.  Find last safe place in line. */
      break_pnt  = LINE_LENGTH - 1;
      save_pnt   = break_pnt;
      white_seen = FALSE;
      done       = FALSE;
      while (!done) {
        --break_pnt;
        if (line[break_pnt] != '\0') {
          if (isspace(line[break_pnt])) {
            if (!white_seen) {
              save_pnt   = break_pnt + 1;
              white_seen = TRUE;
            }
          }
          else if (white_seen) {
            done = TRUE;
          }
        }
      }
    }
    else {
      break_pnt = LINE_LENGTH;
    }

    offset = 0;
  }

  while (offset < LINE_LENGTH && isspace(line[offset])) {
    offset++;
  }
  if (offset == LINE_LENGTH || line[offset] == '%' || line[offset] == '#') {
    *end_flag = 1;
    if (break_pnt < LINE_LENGTH) {
      flush_line(infile);
    }
    return (0);
  }

  ptr = &(line[offset]);
  val = (int)strtol(ptr, &ptr2, 10);

  if (ptr2 == ptr) { /* End of input line. */
    offset    = 0;
    *end_flag = 1;
    return (0);
  }

  offset = (int)(ptr2 - line) / sizeof(char);

  return (val);
}

static void flush_line(FILE *infile)
{
  int c; /* character being read */

  offset = 0;
  c      = getc(infile);
  while (c != '\n' && c != '\f') {
    c = getc(infile);
  }
}
