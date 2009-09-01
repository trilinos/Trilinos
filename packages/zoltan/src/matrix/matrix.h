/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/*****************************************************************************
 * This module wants to be an "abstract" view of user data before being
 * "specialized" into real Zoltan datastructure like HyperGraph or Graph.
 *
 * Basic pre-computations defined here will thus be available for all Zoltan
 * internal data structures:
 *   - Redistribution
 *   - Symmetrization
 *   - Filtering
 */
#ifndef __MATRIX_H
#define __MATRIX_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg_comm.h" /* Not useful in the future ? */
#include "zoltan_dd.h"

/* Would be a lot better in C++ ! */
/*******************/
/* Matrix is a structure to represent user data before specializing it
 * into a graph or an hypergraph.  Currently, X are the objects and Y
 * are the hyperedges for the hypergraph.  For a graph point of view,
 * X and Y are the same entities and pins are the edges.
 ************/

typedef enum {ADD_WEIGHT=0, MAX_WEIGHT, CMP_WEIGHT} WgtOp;

/* This structure defines how the matrix will be constructed */
typedef struct Zoltan_matrix_options_ {
  int enforceSquare;           /* Want to build a graph */
                               /* Cols & rows will have a consistent numbering */
  WgtOp pinwgtop;              /* How to deal with duplicate arcs */
  int randomize;               /* Do we have to randomize input ? */
  int pinwgt;                  /* Do we want pinwgt ? */
  int local;                   /* If only local edges have to be kept */
  int final_output;            /* final_output flag, not used yet */
  int symmetrize;              /* What kind of symmetry we have to apply, not used yet */
  int keep_distribution;       /* Conserve the original distribution, cannot work with randomize */
} Zoltan_matrix_options;

/* This structure is a CS view of a part of the matrix/hypergraph */
typedef struct Zoltan_matrix_ {
  Zoltan_matrix_options opts;  /* How to build the matrix */
  int           redist;        /* HG queries have been used or matrix distribution has changed*/
  int           completed;     /* Matrix is ready to be specialized in HG or G */
  int           globalX;       /* Overall number on X dimension */
  int           globalY;       /* Overall number on Y dimension */
  int           offsetY;       /* Used for bipartite graph: GNO >= offsetY are edges */
  int           nY;            /* Local number in Y dimension */
  int           nY_ori;        /* nY in the initial (user ?) distribution */
  int           ywgtdim;       /* Wgt dimensions for Y */
  int           nPins;         /* Local number of Pins */
  int           pinwgtdim;     /* Wgt dimensions for pins */
  int          *yGNO;          /* Local edges gnos */
  int          *ystart;        /* Indirection array to describe a column */
  int          *yend;          /* end of local pins, usually ystart+1
			         (and is ystart+1 after matrix complete) */
  int          *pinGNO;        /* array of gno of other extremtiy */
  float        *pinwgt;        /* Wgt for pins */

  /* These fields are used only before matrix_complete */
  /* Allow us to move only pins and CSR structure without having to worry
   * about vertex and edge data. */
  struct Zoltan_DD_Struct *ddX; /* Map xGNO -> xGID, xwgt, Input_Parts */
  struct Zoltan_DD_Struct *ddY; /* Map yGNO -> yGID, ywgt */

  /* These fields are used after matrix_complete */
  float        *ywgt;           /* Wgt for local Y */
  ZOLTAN_ID_PTR yGID;           /* Local Y GID */
} Zoltan_matrix;

  /* Overstructure to handle distribution */
typedef struct Zoltan_matrix_2d_ {
  Zoltan_matrix   mtx;          /* The "matrix" */
  PHGComm         *comm;        /* How data are distributed */
  int             *dist_x;      /* Distribution on x axis */
  int             *dist_y;      /* Distribution on y axis */
} Zoltan_matrix_2d;

/* Auxiliary struct used internaly */
typedef struct Zoltan_Arc_ {
    int yGNO;
    int pinGNO;
    int offset;
} Zoltan_Arc;

/*--------------
 * FUNCTIONS
 */

/* Build a matrix from the user data get by the queries.
 * The matrix is only a translation of user data from GID space to GNO space.
 * The matrix cannot be directly used as a graph or hypergraph, you have
 * to make a distribution of this matrix. */
int
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix_options *opt, Zoltan_matrix* matrix);

/* Free a matrix object */
void
Zoltan_Matrix_Free(Zoltan_matrix *m);

/* Free a matrix2d object */
void
Zoltan_Matrix2d_Free(Zoltan_matrix_2d *m);

/* This function compute the indices of the diagonal terms.
   This function needs that diagonal terms are declared at most
   1 time locally.
 */
int
Zoltan_Matrix_Mark_Diag(ZZ* zz, const Zoltan_matrix* const m,
			int *n_nnz, int **nnz);

 /* This function removes nnz which are listed as arguments (list of indexes in
    pin* arrays.
    nnz array have to be sorted.
 */
int
Zoltan_Matrix_Delete_nnz(ZZ* zz, Zoltan_matrix* m,
			 const int n_nnz, const int* const nnz);

/* Performs a permutation of the matrix, perm_y A perm_y^t.
 * At this time we only do symmetric permutations (don't know xGNO !).
 * This call "uncomplete" the matrix and must be followed by
 * Zoltan_Matrix2d_Distribute.
 */
int
Zoltan_Matrix_Permute(ZZ* zz, Zoltan_matrix *m, const int* const perm_y);

/* Distribute the matrix in the 2D layout defined by user in outmat
 * if !copy, inmat is not usable after this call */
int
Zoltan_Matrix2d_Distribute (ZZ* zz, Zoltan_matrix inmat,
			    Zoltan_matrix_2d *outmat, int copy);

/* Compute a 2D datalayout that fit the constraints given.
 */
int
Zoltan_Distribute_layout (ZZ *zz, const PHGComm * const inlayout,
			  int hiRank, int loRank,
			  int reqx, int reqy,
			  PHGComm *outlayout);

/* ShortCuts to compute a layout. */
/* Distribute in a square way, typically for a hypergraph */
int Zoltan_Distribute_Square (ZZ * zz, PHGComm *layout) ;
/* Distribute in a linear 1D way, typically for a graph */
int Zoltan_Distribute_LinearY (ZZ * zz, PHGComm *layout) ;

/* Compute a symmertrization of the matrix.
 * if bipartite == 0, A+At transformation is done (matrix has to be square).
 * else, a bipartite graph of the matrix is built (matrix can have any shape).
 */
int
Zoltan_Matrix_Sym(ZZ* zz, Zoltan_matrix *matrix, int bipartite);


int
Zoltan_Matrix_Remove_DupArcs(ZZ *zz, int size, Zoltan_Arc *arcs, float* pinwgt,
			     Zoltan_matrix *outmat);

/* Function that group duplicate nnz */
int
Zoltan_Matrix_Remove_Duplicates(ZZ *zz, Zoltan_matrix inmat, Zoltan_matrix *outmat);

/* This code has to be called just before specializing the matrix into
 * a graph or an hypergraph.
 */
int
Zoltan_Matrix_Complete(ZZ* zz, Zoltan_matrix* m);

/* This code is used to fill the adjproc array which is used in some
 * place in Zoltan.
 * Perhaps not usefull to do this to call Scotch ?
 */
int
Zoltan_Matrix2d_adjproc (ZZ* zz, const Zoltan_matrix_2d * const mat, int **adjproc);

/* Declare the matrix as empty */
/* Internal use only */
void
Zoltan_Matrix_Reset(Zoltan_matrix* m);


#ifdef __cplusplus
}
#endif

#endif
