/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <exodusII.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Create a separate exodus file on each thread (same file as written
 * by testwt.c)
 */

#define NUM_THREADS 8

typedef struct
{
  long threadid;
  int  exoid;
} param;

void *error_test(void *varg)
{
  param *     arg = (param *)varg;
  int         i   = 0;
  char        name[32];
  const char *routine;
  const char *message;
  int         err_num;
  /*  ex_opts(EX_VERBOSE); */
  sprintf(name, "Thread%ld", arg->threadid);
  for (i = EX_MEMFAIL; i <= EX_INTERNAL; i++) {
    ex_err(name, "Testing thread-safe exodus", i);
    ex_get_err(&message, &routine, &err_num);
    if (err_num != i) {
      fprintf(stderr, "Thread %ld: ERROR: called error (%d) does not match stored value (%d)\n",
              arg->threadid, i, err_num);
    }
  }

  /* Netcdf error codes... (negative) */
  for (i = NC_EBADID; i >= NC4_LAST_ERROR; i--) {
    ex_err(name, "Testing thread-safe exodus", i);
    ex_get_err(&message, &routine, &err_num);
    if (err_num != i) {
      fprintf(stderr, "Thread %ld: ERROR: called error (%d) does not match stored value (%d)\n",
              arg->threadid, i, err_num);
    }
  }

  return arg;
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int       rc;
  long      t;

  param arg[NUM_THREADS];

  printf("Running on %d threads\n", NUM_THREADS);
  for (t = 0; t < NUM_THREADS; t++) {
    arg[t].threadid = t;
    rc              = pthread_create(&threads[t], NULL, error_test, (void *)(arg + t));
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    pthread_join(threads[t], NULL);
  }
}
