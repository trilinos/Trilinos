/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef CHACO_H
#define CHACO_H
/**
 * Calling Chaco from other programs. Throughout this document we have
 * assumed that Chaco is being used as a stand{alone program. However,
 * this needn't be the case. We designed version 2.0 to allow for easy
 * interface with other codes written in either C or Fortran. The
 * mechanism for this interface is described below. Some familiarity
 * with the remainder of this document is assumed.  The interface()
 * routine and can be found in \code/main/interface.c".  This is the
 * routine that Chaco itself invokes after prompting the user for all
 * the necessary input. Consequently, no functionality is lost by
 * calling interface() yourself.
 *
 * The input can still be checked for consistency, all the output
 * options are still active and the ability to modify parameters at
 * run time by reading a file as described in x6.10 is
 * maintained. (The parameters ARCHITECTURE, EIGEN TOLERANCE and
 * RANDOM SEED are made obsolete by the arguments to interface() as
 * detailed below, and DEBUG INPUT and PROMPT become irrelevant, but
 * all other parameters remain active.) The ability to control the
 * goals argument described below actually gives you greater
 * functionality than you would have in standalone mode.
 *
 * The interface routine returns 0 if the partitioning is successful,
 * and 1 otherwise.  Typically, a return code of 1 indicates the
 * detection of some inconsistencies in the input arguments. The
 * arguments to interface() describe the graph, input and output files
 * and arrays, properties of the desired decomposition and the
 * requested partitioning algorithm. The arguments are described below
 * in the order in which they occur.
 *
 * A. Arguments describing the graph.
 *
 *      1. nvtxs. Type int. This is the number of vertices in the graph. Vertices are
 *         numbered from 1 to nvtxs.
 *      2. start. Type int *. Although Chaco internally uses a C structure to represent
 *         the graph, a simpler representation at the start allows for interface with Fortran
 *         programs. The start array is of size (nvtxs+1). It's values are indices into the
 *         adjacency array. The values in adjacency from start i ; 1] to start i];1 are the
 *         vertices adjacent to vertex i in the graph. (Note that C arrays begin at zero,
 *         so in Fortran, the relevant range would be start i] to start(i + 1) ; 1.)
 *      3. adjacency. Type int *. As indicated in the description of start, this array contains
 *         a list of edges for all vertices in the graph. Note that if the FREE GRAPH
 *         parameter from x6.8 is set to TRUE, then after converting to a new data structure,
 *         both start and adjacency are freed. If this is inappropriate for your application
 *         (e.g. you want to keep the graph, or you didn't dynamically allocate
 *         these arrays), then you should set FREE GRAPH to FALSE.
 *      4. vwgts. Type int *. This array of length nvtxs specifies weights for all the
 *         vertices. If you pass in a NULL pointer, then all vertices are given unit weight.
 *         Vertex weights should be positive.
 *
 *      5. ewgts. Type float * (Fortran type real*4). This array
 *         specifies weights for all the edges. It is of the same length
 *         as adjacency and is indexed in the same way. If you use
 *         Kernighan-Lin or the multilevel partitioner, these values
 *         will be rounded to the nearest integer. We suggest scaling
 *         them so they are neither very small nor very big. Edge
 *         weights should be positive.
 *
 *      6. x. Type float *. If you are using the inertial
 *         partitioner, you need to specify geometric coordinates for
 *         each vertex. This array of length nvtxs specifies the x
 *         coordinate for each vertex.
 *
 *      7. y. Type float *. This array specifies the y coordinate for
 *         each vertex. If it is NULL, the geometry is assumed to
 *         one-dimensional.
 *
 *      8. z. Type float *. This array specifies the z coordinate for
 *         he each vertex. If z is NULL and y is not NULL, the geometry
 *         is assumed to be two-dimensional.
 *
 * B. Output le names.
 *
 *      9. outassignname. Type char *. If you desire the final
 *         assignment to be written to a file, this argument gives the
 *         name of that file. If this argument is NULL or if the
 *         parameter OUTPUT ASSIGN is 0, then the assignment is not
 *         written to a file.
 *
 *     10. out filename. Type char *. This is the name of a file in
 *         which the results of the run are printed. If it is NULL or if
 *         the parameter ECHO is not negative, then no file output is
 *         performed.
 *
 * C. Assignment.
 *     11. assignment. Type int *. This is the only output argument to interface().
 *         It is an array of length nvtxs and returns the set number to which each vertex
 *         is assigned. The set number for vertex i is returned in assignment i ; 1] (or for
 *         Fortran, in assignment(i)). This can also be an input argument if global method,
 *         argument 16 below, is set to 7. A description of what functionality can be used
 *         with an input assignment can be found in x4.4
 *   NOTE: This argument was a short in the original implementation and documentation.
 *         Since this limits the processor decompositon to < 32,768 processors, it needed
 *         to be changed to an integer as were all other shorts in the library.
 *
 * D. Description of the target machine.
 *      12. architecture. Type int. This parameter designates the topology of the par-
 *         allel machine for which you are partitioning. Current capabilities include a
 *         hypercube (indicated by a value of 0), and a one-, two- or three-dimensional
 *         mesh (indicated by a value of 1, 2 or 3 respectively.) Note that this argument
 *         overrides the ARCHITECTURE parameter.
 *      13. ndims tot. Type int. If architecture is zero, indicating a hypercube, this
 *         value is the number of dimensions in the hypercube.
 *      14. mesh dims. Type int array of size 3. If architecture is 1, 2 or 3, indicating a
 *         mesh, the values in this array denote the size of the mesh in each dimension.
 *      15. goal. Type double *. This optional array speci es the desired sizes of the
 *         di erent sets. The total number of sets is implicit in the architectural speci -
 *         cations provided by the preceding three parameters. If a null value is passed
 *         for goal, the code will try to make each set have the same vertex weight sum.
 *         If it is not null, the goal array should be as long as the total number of sets.
 *         The value in goal i] (or, for Fortran, goal(i + 1)) should be the desired sum
 *         of vertex weights of vertices assigned to set i. Note that set numbers begin
 *         at zero. Chaco will try to get as close to this goal as possible, but may not
 *         succeed exactly. The sum of all the goals should equal the sum of all the vertex
 *         weights, and values should be nonnegative.
 *         Although the default is to make all set sizes equal, there are applications where
 *         this may be undesirable. One example would be if you are decomposing a
 *         computation among processors of di erent speeds. All the code in Chaco
 *         handles this more general case, and should work for any consistent values in
 *         goal.
 *
 * E. Partitioning options.
 *      16. global method. Type int. This argument speci es the global partitioning
 *         method and should have a value from 1 and 7. These values are the same
 *         as those on the "Global method" menu when running Chaco in standalone
 *         method, as reviewed in x5.4.
 *
 *      17. local method. Type int. This argument speci es the local
 *         partitioning method and should have a value of 1 or 2. These
 *         values are the same as those on the "Local method" menu when
 *         running Chaco in stand{alone method, as reviewed in x5.4.
 *
 *      18. rqi flag. Type int. If you requested spectral partitioning and wish to use the
 *          multilevel RQI/Symmlq eigensolver, this argument should be set to 1. If you
 *          wish instead to use Lanczos, it should be set to 0.
 *
 *      19. vmax. Type int. If you are using either the multilevel-KL partitioner, or
 *         the multilevel RQI/Symmlq eigensolver, you need to specify when the coarsest
 *         graph is small enough. When a coarse graph has no more than vmax vertices,
 *         the recursive coarsening is finished.
 *      20. ndims. Type int. This argument should have a value of 1, 2 or 3 indicating
 *         partitioning by bisection, quadrisection or octasection.
 *      21. eigtol. Type double. If you are using a spectral method or multilevel-KL,
 *         this argument specifies the tolerance you request for the eigensolver. A
 *         discussion of an appropriate choice can be found in the description of the
 *         EIGEN TOLERANCE parameter in x6.2. Note that this argument overrides the
 *         value of the EIGEN TOLERANCE parameter.
 *      22. seed. Type long. This is a seed for the random number generator "rand()".
 *         Note that it overrides the RANDOM SEED parameter.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef _MSC_VER
#define INTER_FACE interFace
#else
#define INTER_FACE interface
#endif

extern int INTER_FACE(int    nvtxs,                 /**< number of vertices in full graph */
                      int *  start,                 /**< start of edge list for each vertex */
                      int *  adjacency,             /**< edge list data */
                      int *  vwgts,                 /**< weights for all vertices */
                      float *ewgts,                 /**< weights for all edges */
                      float *x, float *y, float *z, /**< coordinates for inertial method */
                      char *  outassignname,        /**< name of assignment output file */
                      char *  outfilename,          /**< output file name */
                      int *   assignment,           /**< set number of each vtx (length n) */
                      int     architecture,         /**< 0 => hypercube, d => d-dimensional mesh */
                      int     ndims_tot,     /**< total number of cube dimensions to divide */
                      int     mesh_dims[3],  /**< dimensions of mesh of processors */
                      double *goal,          /**< desired set sizes for each set */
                      int     global_method, /**< global partitioning algorithm */
                      int     local_method,  /**< local partitioning algorithm */
                      int     rqi_flag,      /**< should I use RQI/Symmlq eigensolver? */
                      int     vmax,          /**< how many vertices to coarsen down to? */
                      int     ndims,         /**< number of eigenvectors (2^d sets) */
                      double  eigtol,        /**< tolerance on eigenvectors */
                      long    seed);            /**< for random graph mutations */

/* Chaco interface to read assignment vector from file */
extern int input_assign(FILE *, char *, int, int *);

#define CHACO_VERSION_MAJOR 3
#define CHACO_VERSION_MINOR 0
#define CHACO_VERSION_PATCH 0

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif

#endif
