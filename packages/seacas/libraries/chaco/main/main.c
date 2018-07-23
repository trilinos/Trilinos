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

#include "chaco.h"
#include "defs.h"
#include "params.h"
#include "smalloc.h"
#include <stdio.h>

int main(void)
{
  extern int    Using_Main;                 /* is main routine being called? */
  extern char * Graph_File_Name;            /* name of graph input file */
  extern char * Geometry_File_Name;         /* name of coordinate input file */
  extern char * Assign_In_File_Name;        /* name of assignment input input file */
  extern char * PARAMS_FILENAME;            /* name of file with parameter updates */
  extern double EIGEN_TOLERANCE;            /* tolerance for eigen calculations */
  extern int    OUTPUT_ASSIGN;              /* whether to write assignment to file */
  extern int    DEBUG_MEMORY;               /* debug memory allocation and freeing? */
  extern int    DEBUG_TRACE;                /* trace main execution path */
  extern int    DEBUG_PARAMS;               /* debug flag for reading parameters */
  extern long   RANDOM_SEED;                /* seed for random number generators */
  extern int    ECHO;                       /* controls amount of output */
  extern int    PROMPT;                     /* prompt for input or not? */
  extern int    PRINT_HEADERS;              /* print lines for output sections? */
  extern int    MATCH_TYPE;                 /* matching routine to call */
  extern double input_time;                 /* times data file input */
  extern double start_time;                 /* time partitioning starts */
  FILE *        fin;                        /* input file */
  FILE *        fingeom;                    /* geometry input file (for inertial method) */
  FILE *        finassign;                  /* assignment file if reading in */
  FILE *        params_file;                /* file with parameter value updates */
  double *      goal;                       /* desired set sizes */
  float *       x, *y, *z;                  /* coordinates for inertial method */
  int *         start;                      /* start of edge list for each vertex */
  int *         adjacency;                  /* edge list data */
  float *       ewgts;                      /* weights for all edges */
  int *         vwgts;                      /* weights for all vertices */
  int           global_method;              /* global partitioning method */
  int           local_method;               /* local partitioning method */
  int *         assignment;                 /* set number of each vtx (length nvtxs+1) */
  double        eigtol;                     /* tolerance in eigenvector calculation */
  int           nvtxs;                      /* number of vertices in graph */
  int           ndims;                      /* dimension of recursive partitioning */
  int           architecture;               /* 0 => hypercube, d => d-dimensional mesh */
  int           ndims_tot;                  /* total number of cube dimensions to divide */
  int           mesh_dims[3];               /* dimensions of mesh of processors */
  long          seed;                       /* for random graph mutations */
  int           rqi_flag;                   /* use RQI/Symmlq eigensolver? */
  int           vmax;                       /* if so, how many vertices to coarsen down to? */
  int           igeom;                      /* geometry dimension if inertial method */
  char          graphname[NAME_LENGTH];     /* name of graph input file */
  char          geomname[NAME_LENGTH];      /* name of geometry input file */
  char          inassignname[NAME_LENGTH];  /* assignment input file name */
  char          outassignname[NAME_LENGTH]; /* assignment output file name */
  char          outfilename[NAME_LENGTH];   /* name of output file */
  char *        outassignptr;               /* name or null pointer for output assignment */
  char *        outfileptr;                 /* name or null pointer for output file */
  int           another;                    /* run another problem? */
  double        time;                       /* timing marker */
  int           flag;                       /* return code from input routines */
  double        seconds();                  /* returns elapsed time in seconds */
  int           affirm();
  int           input_graph(), input_geom();
  void          input_queries(), read_params(), clear_timing();

  /*malloc_debug(2);*/

  if (DEBUG_TRACE > 0) {
    printf("<Entering main>\n");
  }

  if (PRINT_HEADERS) {
    printf("\n                    Chaco 2.0\n");
    printf("          Sandia National Laboratories\n\n");
  }

  Using_Main  = TRUE;
  another     = TRUE;
  params_file = fopen(PARAMS_FILENAME, "r");
  if (params_file == NULL && DEBUG_PARAMS > 1) {
    printf("Parameter file `%s' not found; using default parameters.\n", PARAMS_FILENAME);
  }

  while (another) {

    start_time = time = seconds();

    x = y = z  = NULL;
    goal       = NULL;
    assignment = NULL;

    read_params(params_file);

    input_queries(&fin, &fingeom, &finassign, graphname, geomname, inassignname, outassignname,
                  outfilename, &architecture, &ndims_tot, mesh_dims, &global_method, &local_method,
                  &rqi_flag, &vmax, &ndims);

    if (global_method == 7)
      Assign_In_File_Name = inassignname;
    else
      Assign_In_File_Name = NULL;

    if (OUTPUT_ASSIGN > 0)
      outassignptr = outassignname;
    else
      outassignptr = NULL;

    if (ECHO < 0)
      outfileptr = outfilename;
    else
      outfileptr = NULL;

    flag = input_graph(fin, graphname, &start, &adjacency, &nvtxs, &vwgts, &ewgts);
    if (flag) {
      sfree(ewgts);
      sfree(vwgts);
      sfree(adjacency);
      sfree(start);
      goto skip;
    }

    Graph_File_Name = graphname;

    assignment = smalloc(nvtxs * sizeof(int));
    if (global_method == 7) {
      flag = input_assign(finassign, inassignname, nvtxs, assignment);
      if (flag)
        goto skip;
      Assign_In_File_Name = inassignname;
    }

    if (global_method == 3 ||
        (MATCH_TYPE == 5 && (global_method == 1 || (global_method == 2 && rqi_flag)))) {
      /* Read in geometry data. */
      flag = input_geom(fingeom, geomname, nvtxs, &igeom, &x, &y, &z);
      if (flag)
        goto skip;
      Geometry_File_Name = geomname;
    }
    else {
      x = y = z = NULL;
    }

    input_time += seconds() - time;

    eigtol = EIGEN_TOLERANCE;
    seed   = RANDOM_SEED;

    interface(nvtxs, start, adjacency, vwgts, ewgts, x, y, z, outassignptr, outfileptr, assignment,
              architecture, ndims_tot, mesh_dims, goal, global_method, local_method, rqi_flag, vmax,
              ndims, eigtol, seed);

  skip:
    if (global_method == 3) {
      if (z != NULL)
        sfree(z);
      if (y != NULL)
        sfree(y);
      if (x != NULL)
        sfree(x);
    }
    sfree(assignment);

    if (DEBUG_MEMORY > 0) {
      printf("\n");
    }

    if (PROMPT) {
      another = affirm("\nRun Another Problem");
    }
    else {
      another = affirm("");
    }
    if (another) {
      clear_timing();
      printf("\n------------------------------------------------\n\n");
      fflush(stdout);
    }
  }
  if (params_file != NULL)
    fclose(params_file);

  if (DEBUG_TRACE > 1) {
    printf("<Leaving main>\n");
  }

  return (0);
}
