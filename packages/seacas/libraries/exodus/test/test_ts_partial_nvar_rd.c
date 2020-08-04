/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <exodusII.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Create an exodus file with NUM_NODES nodes and 0 elements.
 * Initial mesh is created on main thread.
 * There are
 * Nodal var data and names is written on multiple threads
 * -- Each thread handles a portion of the nodes for each nodal variable
 * -- Each thread handles a portion of the nodes for each coordinate variable
 */

#define NUM_THREADS 8
#define NUM_NODES 64 /* should be multiple of NUM_THREADS */
#define NUM_NODAL_VAR 4

static int ulpsDistance(float a, float b)
{
  int ia, ib;
  if (a == b)
    return 0;
  memcpy(&ia, &a, sizeof(a));
  memcpy(&ib, &b, sizeof(b));

  if ((ia < 0) != (ib < 0))
    return INT_MAX;
  int distance = ia - ib;
  if (distance < 0)
    distance = -distance;
  return distance;
}
static int approx_equal(float f1, float f2) { return ulpsDistance(f1, f2) <= 2; }

typedef struct
{
  long threadid;
  int  exoid;
  int  timestep;
} param;

void *input_nodal_var(void *varg)
{
  param *arg      = (param *)varg;
  int    num_node = ex_inquire_int(arg->exoid, EX_INQ_NODES);
  int    segment  = num_node / NUM_THREADS;
  int    begin    = arg->threadid * segment;

  float *data = malloc(segment * sizeof(float));
  int    i, var;

  if (arg->timestep == 1) {
    ex_get_partial_coord(arg->exoid, begin + 1, segment, data, NULL, NULL);

    for (i = 0; i < segment; i++) {
      if (!approx_equal(data[i], (float)(begin + i) / 100.0)) {
        fprintf(stderr,
                "ERROR: Thread %ld: X Coordinate mismatch at node %d: Got: %f, expected %f\n",
                arg->threadid, i + 1, data[i], (float)(begin + i) / 100.0);
      }
    }

    ex_get_partial_coord(arg->exoid, begin + 1, segment, NULL, data, NULL);

    for (i = 0; i < segment; i++) {
      if (!approx_equal(data[i], (float)(begin + i) / 100.0)) {
        fprintf(stderr,
                "ERROR: Thread %ld: Y Coordinate mismatch at node %d: Got: %f, expected %f\n",
                arg->threadid, i + 1, data[i], (float)(begin + i) / 100.0);
      }
    }

    ex_get_partial_coord(arg->exoid, begin + 1, segment, NULL, NULL, data);

    for (i = 0; i < segment; i++) {
      if (!approx_equal(data[i], (float)(begin + i) / 100.0)) {
        fprintf(stderr,
                "ERROR: Thread %ld: Z Coordinate mismatch at node %d: Got: %f, expected %f\n",
                arg->threadid, i + 1, data[i], (float)(begin + i) / 100.0);
      }
    }
  }

  if (arg->threadid < NUM_NODAL_VAR) {
    char db_name[33];
    char ex_name[33];

    sprintf(ex_name, "NodalVar%ld", arg->threadid + 1);
    ex_get_variable_name(arg->exoid, EX_NODAL, arg->threadid + 1, db_name);
    if (strcmp(db_name, ex_name) != 0) {
      fprintf(stderr,
              "ERROR: Thread %ld: Incorrect variable name for variable %ld: Got: %s, expected %s\n",
              arg->threadid, arg->threadid + 1, db_name, ex_name);
    }
  }

  for (var = 1; var <= NUM_NODAL_VAR; var++) {
    ex_get_partial_var(arg->exoid, arg->timestep, EX_NODAL, var, 1, begin + 1, segment, data);
    for (i = 0; i < segment; i++) {
      if (!approx_equal(data[i], (arg->timestep - 1) * 10 + var + (float)(begin + i) / 100.0)) {
        fprintf(stderr,
                "ERROR: Thread %ld: Nodal variable data mismatch in variable %d at node %d: "
                "Got: %f, expected %f\n",
                arg->threadid, var, i + 1, data[i],
                (arg->timestep - 1) * 10 + var + (float)(begin + i) / 100.0);
      }
    }
  }
  free(data);
  return NULL;
}

int init_file(int num_nodal_vars)
{
  int   exoid;
  float version;
  int   CPU_word_size = 0; /* sizeof(float) */
  int   IO_word_size  = 0; /* use what is stored in file */

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  if (exoid < 0)
    exit(1);

  if (!ex_inquire_int(exoid, EX_INQ_THREADSAFE)) {
    fprintf(stderr,
            "ERROR: This exodus library is not compiled to allow thread-safe operations.\n");
    exit(1);
  }
  return exoid;
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int       rc;
  long      t;

  int exoid = init_file(NUM_NODAL_VAR);

  param arg[NUM_THREADS];

  printf("Running on %d threads\n", NUM_THREADS);
  for (t = 0; t < NUM_THREADS; t++) {
    arg[t].exoid    = exoid;
    arg[t].threadid = t;
    arg[t].timestep = 1;

    rc = pthread_create(&threads[t], NULL, input_nodal_var, (void *)(arg + t));
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    pthread_join(threads[t], NULL);
  }

  ex_close(exoid);
}
