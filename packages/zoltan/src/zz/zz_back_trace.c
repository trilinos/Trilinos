// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdio.h>

/*****************************************************************************/
/*****************************************************************************/
/* Save the ZOLTAN_TRACE stack to be printed on ZOLTAN_PRINT_ERROR           */
/*****************************************************************************/

#define TRACE_SIZE 30
#define TRACE_STRING_LENGTH 128

static char trace[TRACE_SIZE][TRACE_STRING_LENGTH];
static int trace_base=-1;
static int trace_top=-1;

void Zoltan_add_back_trace(char *yo)
{
int i;
char *c;

  if (trace_top < 0){   /* nothing in trace at this point */
    trace_base = trace_top = 0;
  }
  else{
    trace_top++;
    if (trace_top == TRACE_SIZE) trace_top= 0;

    if (trace_top == trace_base){
      trace_base = trace_top + 1;
      if (trace_base == TRACE_SIZE) trace_base = 0;
    }
  }

  c = trace[trace_top];   /* snprintf is not portable */
  for (i=0; (i < TRACE_STRING_LENGTH-1) && *yo; i++){
    *c++ = *yo++;
  }
  *c = '\0';
}

void Zoltan_remove_back_trace()
{
  if (trace_top < 0) return;  /* trace is empty */

  if (trace_top == trace_base){
    trace_top = trace_base = -1;
  }
  else{
    trace_top--;
    if (trace_top < 0) trace_top = TRACE_SIZE - 1;
  }
}

void Zoltan_print_trace(int rank)
{
int t, i;

  if (trace_top < 0) return;

  t = trace_top;

  fprintf(stderr,"\n[%d] Trace:\n",rank);
  for (i=0; i < TRACE_SIZE; i++){

    fprintf(stderr,"[%d] (%d) %s\n", rank, t, trace[t]);
    if (t == trace_base) break;

    t--;
    if (t < 0) t = TRACE_SIZE-1;
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
