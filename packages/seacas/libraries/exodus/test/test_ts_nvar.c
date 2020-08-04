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

/*
 * Create an exodus file with NUM_NODES nodes and 0 elements.
 * Initial mesh is created on main thread.
 * Nodal var data and names is written on multiple threads
 * -- Each thread handles a single nodal variable
 * -- Coordinate data is written by first 3 threads
 */

#define NUM_THREADS 8
#define NUM_NODES 64

typedef struct
{
  long threadid;
  int  exoid;
  int  timestep;
} param;

void *output_nodal_var(void *varg)
{
  char   name[33];
  param *arg      = (param *)varg;
  int    num_node = ex_inquire_int(arg->exoid, EX_INQ_NODES);
  float *data     = malloc(num_node * sizeof(float));
  int    i;
  for (i = 0; i < num_node; i++) {
    data[i] = (arg->timestep - 1) * 10 + arg->threadid + 1 + (float)i / 100.0;
  }

  if (arg->timestep == 1) {
    switch (arg->threadid) {
    case 0: ex_put_coord(arg->exoid, data, NULL, NULL); break;
    case 1: ex_put_coord(arg->exoid, NULL, data, NULL); break;
    case 2: ex_put_coord(arg->exoid, NULL, NULL, data); break;
    default: break;
    }
  }

  sprintf(name, "NodalVar%ld", arg->threadid + 1);
  ex_put_variable_name(arg->exoid, EX_NODAL, arg->threadid + 1, name);

  ex_put_var(arg->exoid, arg->timestep, EX_NODAL, arg->threadid + 1, 1, num_node, data);
  free(data);
  return NULL;
}

int init_file(int num_nodal_vars)
{
  /* Specify compute and i/o word size */
  int CPU_word_size = 0; /* sizeof(float) */
  int IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */
  int exoid = ex_create("test.exo",     /* filename path */
                        EX_CLOBBER,     /* create mode */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size); /* I/O float word size in bytes */
  printf("after ex_create for test.exo, exoid = %d\n", exoid);
  printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

  if (!ex_inquire_int(exoid, EX_INQ_THREADSAFE)) {
    fprintf(stderr,
            "ERROR: This exodus library is not compiled to allow thread-safe operations.\n");
    exit(1);
  }

  /* initialize file with parameters */
  int num_dim       = 3;
  int num_nodes     = NUM_NODES;
  int num_elem      = 0;
  int num_elem_blk  = 0;
  int num_node_sets = 0;
  int num_side_sets = 0;

  int error = ex_put_init(exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk,
                          num_node_sets, num_side_sets);

  printf("after ex_put_init, error = %d\n", error);

  error = ex_put_variable_param(exoid, EX_NODAL, num_nodal_vars);
  printf("after ex_put_variable_param, error = %d\n", error);
  if (error) {
    ex_close(exoid);
    exit(-1);
  }

  float time_value = 0.0;
  error            = ex_put_time(exoid, 1, &time_value);
  printf("after ex_put_time, error = %d\n", error);

  return exoid;
}

int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int       rc;
  long      t;

  int exoid = init_file(NUM_THREADS);

  param arg[NUM_THREADS];

  printf("Running on %d threads\n", NUM_THREADS);
  for (t = 0; t < NUM_THREADS; t++) {
    arg[t].exoid    = exoid;
    arg[t].threadid = t;
    arg[t].timestep = 1;

    rc = pthread_create(&threads[t], NULL, output_nodal_var, (void *)(arg + t));
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    pthread_join(threads[t], NULL);
  }

  ex_close(exoid);

  /* Last thing that main() should do */
}
