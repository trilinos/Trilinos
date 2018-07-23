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

#include <stdio.h> // for printf, scanf, NULL, fopen, etc

void input_queries(FILE **fin,           /* input file */
                   FILE **fingeom,       /* geometry input file (for inertial method) */
                   FILE **finassign,     /* assignment input file (if just simulating) */
                   char * inname,        /* name of graph input file */
                   char * geomname,      /* name of geometry input file */
                   char * inassignname,  /* name of assignment input file */
                   char * outassignname, /* name of assignment output file */
                   char * outfilename,   /* name of file for outputting run results */
                   int *  architecture,  /* 0=> hypercube, d=> d-dimensional mesh */
                   int *  ndims_tot,     /* target number of hypercube dimensions */
                   int    mesh_dims[3],  /* mesh dimensions */
                   int *  global_method, /* what global partitioning strategy to use? */
                   int *  local_method,  /* what local refinement strategy to use? */
                   int *  rqi_flag,      /* should I use multilevel eigensolver? */
                   int *  vmax,          /* if so, how far should I coarsen? */
                   int *  ndims          /* number of divisions at each stage */
)
{
  extern int SEQUENCE;      /* sequence instead of partition graph? */
  extern int ARCHITECTURE;  /* 0=> hypercube, d=> d-dimensional mesh */
  extern int OUTPUT_ASSIGN; /* write assignments to file? */
  extern int ECHO;          /* copy input to screen? results to file? */
  extern int DEBUG_TRACE;   /* trace main execution path */
  extern int PROMPT;        /* prompt for input? */
  extern int MATCH_TYPE;    /* max-matching routine to call */
  int        eigensolver;   /* which kind of eigensolver to use */
  int        nprocs;        /* number of processors being divided into */

  int input_int();

  if (DEBUG_TRACE > 0) {
    printf("<Entering input_queries>\n");
  }

  *architecture = ARCHITECTURE;

  /* Name and open input graph file. */
  *fin = NULL;
  while (*fin == NULL) {
    if (PROMPT) {
      printf("Graph input file: ");
    }
    (void)scanf("%s", inname);

    *fin = fopen(inname, "r");
    if (*fin == NULL) {
      printf("Graph file `%s' not found.\n", inname);
    }
  }

  /* Name output assignment file. */
  if (OUTPUT_ASSIGN && !SEQUENCE) {
    if (PROMPT) {
      printf("Assignment output file: ");
    }
    (void)scanf("%s", outassignname);
  }

  /* Name output results file. */
  if (ECHO < 0) {
    if (PROMPT) {
      printf("File name for saving run results: ");
    }
    (void)scanf("%s", outfilename);
  }

  /* Initialize the method flags */
  *rqi_flag      = 0;
  *global_method = 0;
  *fingeom       = NULL;
  *finassign     = NULL;

  /* Get global method, if any. */
  if (SEQUENCE) {
    *global_method = 2;
  }
  else {
    while (*global_method < 1 || *global_method > 7) {
      if (PROMPT) {
        printf("Global partitioning method:\n");
        printf("  (1) Multilevel-KL\n");
        printf("  (2) Spectral\n");
        printf("  (3) Inertial\n");
        printf("  (4) Linear\n");
        printf("  (5) Random\n");
        printf("  (6) Scattered\n");
        printf("  (7) Read-from-file\n");
      }
      *global_method = input_int();
    }
  }

  if (*global_method == 7) { /* Name and open input assignment file. */
    while (*finassign == NULL) {
      if (PROMPT) {
        printf("Assignment input file: ");
      }
      (void)scanf("%s", inassignname);

      *finassign = fopen(inassignname, "r");
      if (*finassign == NULL) {
        printf("Assignment file `%s' not found.\n", inassignname);
      }
    }
  }

  else if (*global_method == 3) {
    while (*fingeom == NULL) {
      if (PROMPT) {
        printf("Geometry input file name: ");
      }
      (void)scanf("%s", geomname);

      *fingeom = fopen(geomname, "r");
      if (*fingeom == NULL) {
        printf("Geometry input file `%s' not found.\n", geomname);
      }
    }
  }
  else if (*global_method == 2) {
    eigensolver = 0;
    while (eigensolver < 1 || eigensolver > 2) {
      if (PROMPT) {
        printf("Eigensolver:\n");
        printf("  (1) Multilevel RQI/Symmlq\n");
        printf("  (2) Lanczos\n");
      }
      eigensolver = input_int();
    }
    if (eigensolver == 1) {
      if (MATCH_TYPE == 5) { /* geometric matching */
        while (*fingeom == NULL) {
          if (PROMPT) {
            printf("Geometry input file name: ");
          }
          (void)scanf("%s", geomname);

          *fingeom = fopen(geomname, "r");
          if (*fingeom == NULL) {
            printf("Geometry input file `%s' not found.\n", geomname);
          }
        }
      }
      *rqi_flag = 1;
      *vmax     = 0;
      while (*vmax <= 0) {
        if (PROMPT) {
          printf("Number of vertices to coarsen down to: ");
        }
        *vmax = input_int();
      }
    }
  }
  else if (*global_method == 1) {
    if (MATCH_TYPE == 5) { /* geometric matching */
      while (*fingeom == NULL) {
        if (PROMPT) {
          printf("Geometry input file name: ");
        }
        (void)scanf("%s", geomname);

        *fingeom = fopen(geomname, "r");
        if (*fingeom == NULL) {
          printf("Geometry input file `%s' not found.\n", geomname);
        }
      }
    }
    *vmax = 0;
    while (*vmax <= 1) {
      if (PROMPT) {
        printf("Number of vertices to coarsen down to: ");
      }
      *vmax = input_int();
    }
  }

  if (SEQUENCE) {
    *local_method = 2;
    if (*architecture == 0) {
      *ndims_tot = 1;
    }
    else if (*architecture > 0) {
      mesh_dims[0] = 2;
      mesh_dims[1] = mesh_dims[2] = 1;
    }
    *ndims = 1;
    goto End_Label;
  }

  /* Get local method, if any */
  *local_method = 0;
  if (*global_method == 1) {
    *local_method = 1;
  }
  else {
    while (*local_method < 1 || *local_method > 2) {
      if (PROMPT) {
        printf("Local refinement method:\n");
        printf("  (1) Kernighan-Lin\n");
        printf("  (2) None\n");
      }
      *local_method = input_int();
    }
  }

  /* Now learn about the parallel architecture. */
  if (*architecture == 0) {
    /* Get total number of hypercube dimensions in which to partition. */
    *ndims_tot = 0;
    while (*ndims_tot < 1) {
      if (PROMPT) {
        printf("Total number of target hypercube dimensions: ");
      }
      *ndims_tot = input_int();
      if (*ndims_tot < 1) {
        printf(" Number of divisions must be at least 1\n");
      }
    }
    nprocs = 1 << (*ndims_tot);
  }

  else { /* Get dimensions of mesh. */
    mesh_dims[1] = mesh_dims[2] = 1;
    if (*architecture == 2) {
      if (PROMPT) {
        printf("X and Y extent of of 2-D mesh: ");
      }
      mesh_dims[0] = input_int();
      mesh_dims[1] = input_int();
    }
    else if (*architecture == 3) {
      if (PROMPT) {
        printf("X, Y and Z extent of 3-D mesh: ");
      }
      mesh_dims[0] = input_int();
      mesh_dims[1] = input_int();
      mesh_dims[2] = input_int();
    }
    else { /* Anything else => 1-D mesh */
      if (PROMPT) {
        printf("Size of 1-D mesh: ");
      }
      mesh_dims[0]  = input_int();
      *architecture = 1;
    }
    nprocs = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
  }

  /* Get number of dimensions in which to partition at each level. */
  *ndims = 0;
  if (nprocs <= 3) {
    *ndims = 1;
  }
  else if (nprocs <= 7) {
    if (PROMPT) {
      printf("Partitioning dimension: \n");
    }
    while (*ndims < 1 || *ndims > 2) {
      if (PROMPT) {
        printf("  (1) Bisection\n");
        printf("  (2) Quadrisection\n");
      }
      *ndims = input_int();
    }
  }
  else {
    if (PROMPT) {
      printf("Partitioning dimension: \n");
    }
    while (*ndims < 1 || *ndims > 3) {
      if (PROMPT) {
        printf("  (1) Bisection\n");
        printf("  (2) Quadrisection\n");
        printf("  (3) Octasection\n");
      }
      *ndims = input_int();
    }
  }
End_Label:

  if (*global_method == 1 || *rqi_flag) {
    if (*vmax < 2 * (1 << *ndims)) {
      *vmax = 2 * (1 << *ndims);
    }
  }
}
