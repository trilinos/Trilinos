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

#include "params.h" // for NAME_LENGTH

#define TRUE 1
#define FALSE 0

/* Input and output control parameters */

int CHECK_INPUT    = TRUE;  /* Check input for consistency? (TRUE/FALSE) */
int ECHO           = 2;     /* Print input/param options? to file? (-2..2) */
int OUTPUT_METRICS = 2;     /* Controls displaying of results (-2..2) */
int OUTPUT_TIME    = 2;     /* At what level to display timings (0..2) */
int OUTPUT_ASSIGN  = FALSE; /* Write assignments to file? (TRUE/FALSE) */
int OUT_ASSIGN_INV = FALSE; /* If so, use inverse form? (TRUE/FALSE) */
int IN_ASSIGN_INV  = FALSE; /* Input file in inverse form? (TRUE/FALSE) */
int PROMPT         = TRUE;  /* Prompt for input? (TRUE/FALSE) */
int PRINT_HEADERS  = TRUE;  /* Print pretty output headers (TRUE/FALSE) */

/* Eigenvector calculation parameters */

int LANCZOS_TYPE = 3;                   /* type of Lanczos to use */
                                        /* 1 => full orthog, 2 => full inverse operator */
                                        /* 3 =>  selective orthogonalization */
double EIGEN_TOLERANCE          = 1e-3; /* Numerical eigen-tolerance */
double SRESTOL                  = -1.;  /* Rel resid tol on T evec; autoset if <= 0 */
int    LANCZOS_SO_INTERVAL      = 10;   /* Itns. between SO orthog checks; set >= 2 */
int    LANCZOS_MAXITNS          = -1;   /* Max Lanczos its; autoset if <= 0 */
double BISECTION_SAFETY         = 10;   /* Divides Lanczos bisection tol */
int    LANCZOS_CONVERGENCE_MODE = 0;    /* Lanczos convergence test type: */
                                        /* 0=> residual,  1=> partition */
int RQI_CONVERGENCE_MODE = 1;           /* RQI convergence test type: */
                                        /* 0=> residual,  1=> partition */
int    LANCZOS_SO_PRECISION = 2;        /* 2 => double Lanczos, 1 => float */
int    WARNING_EVECS        = 2;        /* Warnings in eigenvector generation (0..3) */
double WARNING_ORTHTOL      = 2;        /* Warning if Ares and bjitol have this ratio */
double WARNING_MISTOL       = 100;      /* Warning if Ares and bjitol have this ratio */
int    LANCZOS_TIME         = FALSE;    /* Detailed Lanczos times? (TRUE/FALSE) */
int    TIME_KERNELS         = FALSE;    /* Time numerical kernels? (TRUE/FALSE) */

/* Other parameters for spectral methods */

int    MAKE_CONNECTED = TRUE;   /* Connect graph if using spectral method? */
int    PERTURB        = TRUE;   /* Randomly perturb matrix in spectral method? */
int    NPERTURB       = 2;      /* If so, how many edges to modify? */
double PERTURB_MAX    = 3.0e-3; /* Largest value for perturbation */
int    MAPPING_TYPE   = 1;      /* How to map from eigenvectors to partition */
                                /* 0 => cut at origin, 1 => min-cost assign */
int COARSE_NLEVEL_RQI = 2;      /* # levels between RQI calls in uncoarsening */
int OPT3D_NTRIES      = 5;      /* # local opts to look for global min in opt3d */

/* Kernighan--Lin/Fiduccia--Mattheyses parameters */

int    KL_METRIC     = 2;    /* KL interset cost: 1=>cuts, 2=>hops */
int    KL_RANDOM     = TRUE; /* Use randomness in Kernighan-Lin? (TRUE/FALSE)*/
int    KL_BAD_MOVES  = 20;   /* Number of unhelpful moves in a row allowed */
int    KL_NTRIES_BAD = 1;    /* # unhelpful passes before quitting KL */
int    KL_UNDO_LIST  = TRUE; /* Only resort changed vtxs? (TRUE/FALSE) */
double KL_IMBALANCE  = 0.0;  /* Fractional imbalance allowed by KL */

/* Coarsening parameters */

double COARSEN_RATIO_MIN = .7;    /* Min vtx reduction each coarsen stage */
int    COARSE_NLEVEL_KL  = 2;     /* # levels between KL calls in uncoarsening */
int    MATCH_TYPE        = 1;     /* Type of contraction matching (1..5) */
int    HEAVY_MATCH       = FALSE; /* Encourage heavy match edges? (TRUE/FALSE) */
int    COARSE_KL_BOTTOM  = TRUE;  /* Force KL at lowest level (TRUE/FALSE) */
int    COARSEN_VWGTS     = TRUE;  /* Sum vtx weights in coarsening? (TRUE/FALSE) */
int    COARSEN_EWGTS     = TRUE;  /* Sum edge weights in coarsening? (TRUE/FALSE) */
int    KL_ONLY_BNDY      = TRUE;  /* Start moving vtxs on boundary? (TRUE/FALSE) */

/* Parameters for post-processing options */

int REFINE_PARTITION  = FALSE; /* Postprocess to improve cuts? */
int INTERNAL_VERTICES = FALSE; /* ... to up internal vtxs? (TRUE/FALSE) */
int CONNECTED_DOMAINS = FALSE; /* ... to force connectivity? (T/F) */
int REFINE_MAP        = FALSE; /* ... to reduce hops? (T/F) */

/* Architecture parameters */

int ARCHITECTURE = 0; /* 0=> hypercube, d=> d-dimensional mesh (0..3)*/

/* Miscellaneous parameters */

int    TERM_PROP                 = FALSE;          /* Invoke terminal propagation? (TRUE/FALSE) */
double CUT_TO_HOP_COST           = 1;              /* ..if so, relative importance of cuts/hops */
int    SEQUENCE                  = FALSE;          /* Only do spectral ordering? (TRUE/FALSE) */
char   SEQ_FILENAME[NAME_LENGTH] = "Sequence.out"; /* If so, file name */
long   RANDOM_SEED               = 7654321L;       /* Seed for random number generator */
int    NSQRTS                    = 1000;           /* # square roots to precompute if coarsening */
int    MAKE_VWGTS                = FALSE;          /* Make vtx weights degrees+1? (TRUE/FALSE) */
int    FREE_GRAPH                = TRUE;           /* Free input graph data? (TRUE/FALSE) */
char * PARAMS_FILENAME           = "User_Params";  /* File of parameter changes */

/* Parameters that control debugging output */

int DEBUG_EVECS       = 0; /* Debug flag for eigenvector generation (0..5) */
int DEBUG_KL          = 0; /* Debug flag for Kernighan-Lin (0..3) */
int DEBUG_INERTIAL    = 0; /* Debug flag for inertial method (0..1) */
int DEBUG_CONNECTED   = 0; /* Debug flag for connected components (0..1) */
int DEBUG_PERTURB     = 0; /* Debug flag for matrix perturbation (0..1) */
int DEBUG_ASSIGN      = 0; /* Debug flag for assignment to sets (0..1) */
int DEBUG_OPTIMIZE    = 0; /* Debug flag for optimization/rotation (0..2) */
int DEBUG_BPMATCH     = 0; /* Debug flag for bipartite matching code (0..2) */
int DEBUG_COARSEN     = 0; /* Debug flag for coarsening/uncoarsening (0..1) */
int DEBUG_MEMORY      = 0; /* Debug flag for smalloc/sfree (0..3) */
int DEBUG_INPUT       = 0; /* Debug flag for having read input files (0..1) */
int DEBUG_PARAMS      = 2; /* Debug flag for reading parameter file (0..2) */
int DEBUG_INTERNAL    = 0; /* Debug flag for internal vertices (0..2) */
int DEBUG_REFINE_PART = 0; /* Debug flag for refine partition (0..1) */
int DEBUG_REFINE_MAP  = 0; /* Debug flag for refining mapping (0..1) */
int DEBUG_SIMULATOR   = 0; /* Debug flag for comm simulator (0..2) */
int DEBUG_TRACE       = 0; /* Trace main execution path (0..2) */
int DEBUG_MACH_PARAMS = 0; /* Print computed machine params? (0..1) */
