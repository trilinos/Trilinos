/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
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
#ifndef __ZOLTAN_MATRIX_H
#define __ZOLTAN_MATRIX_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifdef KDDKDD_DEBUG
#define KDDKDDKDD(me, s) if (me==0){printf("KDDKDD %s  ", s);fflush(stdout);system("date");}
#else
#define KDDKDDKDD(me, s) {}
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
typedef enum {MATRIX_FULL_DD=0, MATRIX_FAST, MATRIX_NO_REDIST} SpeedOpt;


/* Hash function used to describe how to distribute data */
/* (Y, X, data, &part_y) */
typedef int distFnct(ZOLTAN_GNO_TYPE, ZOLTAN_GNO_TYPE, void *, int*);
int Zoltan_Distribute_Origin(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y);
int Zoltan_Distribute_Linear(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y);
int Zoltan_Distribute_Cyclic(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y);
int Zoltan_Distribute_Partition(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y);

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
  SpeedOpt speed;
  int fast_build_base;         /* smallest GID (typically 0 or 1); for fast
                                  builds only.  User-specified. */
} Zoltan_matrix_options;


/* Confusing memory allocations: Some fields are allocated, and then their
 *  pointers are put in another structure that frees the memory when done.
 *  We need to keep track of which fields get freed and which don't.
 */

#define FIELD_DIST_Y  0
#define FIELD_YSTART  1
#define FIELD_PINGNO  2
#define FIELD_YGID    3
#define FIELD_PINWGT 4
#define FIELD_NUMBER_OF_FIELDS  5

#define FIELD_FREE_WHEN_DONE(flag, x)  (flag |= (1 << x))
#define FIELD_DO_NOT_FREE_WHEN_DONE(flag, x)  (flag &= ~(1 << x))
#define FIELD_QUERY_DO_FREE(flag, x)  (flag & (1 << x))

/* This structure is a CS view of a part of the matrix/hypergraph */
typedef struct Zoltan_matrix_ {
  Zoltan_matrix_options opts;  /* How to build the matrix */
  int           redist;        /* HG queries have been used or matrix distribution has changed*/
  int           completed;     /* Matrix is ready to be specialized in HG or G */
  int           bipartite;
  ZOLTAN_GNO_TYPE globalX;       /* Overall number on X dimension */
  ZOLTAN_GNO_TYPE globalY;       /* Overall number on Y dimension */
  int           nY;            /* Local number in Y dimension */
  int           nY_ori;        /* nY in the initial (user ?) distribution */
  int           ywgtdim;       /* Wgt dimensions for Y  TODO: where are the weights*/
  int           nPins;         /* Local number of Pins */
  int           pinwgtdim;     /* Wgt dimensions for pins */
  ZOLTAN_GNO_TYPE  *yGNO;       /* Local edges gnos */
  int          *ystart;        /* Indirection array to describe a column */
  int          *yend;          /* end of local pins, usually ystart+1
			         (and is ystart+1 after matrix complete) */
  ZOLTAN_GNO_TYPE  *pinGNO;     /* array of gno of other extremtiy */
  float        *pinwgt;        /* Wgt for pins */

  /* These fields are used only before matrix_complete */
  /* Allow us to move only pins and CSR structure without having to worry
   * about vertex and edge data. */
  struct Zoltan_DD_Struct *ddX; /* Map xGNO -> xGID, xpid */
  struct Zoltan_DD_Struct *ddY; /* Map yGNO -> yGID, ypid */

  /* These fields are used after matrix_complete */
  ZOLTAN_ID_PTR yGID;           /* Local Y GID */

  int *ypid;           /* Initial processor */
  int *ybipart;
} Zoltan_matrix;

  /* Overstructure to handle distribution */
typedef struct Zoltan_matrix_2d_ {
  Zoltan_matrix   mtx;          /* The "matrix" */
  PHGComm         *comm;        /* How data are distributed */
  ZOLTAN_GNO_TYPE *dist_x;      /* Distribution on x axis */
  ZOLTAN_GNO_TYPE *dist_y;      /* Distribution on y axis */
  distFnct        *hashDistFct; /* How to distribute nnz */
  void            *hashDistData;/* Used by hashDist */
  int             delete_flag;  /* 0x001: dist_y, 0x010: ystart, 
                                 * 0x100: pinGNO, 0x1000: yGID */

} Zoltan_matrix_2d;

/* Auxiliary struct used internaly */
typedef struct Zoltan_Arc_ {
    ZOLTAN_GNO_TYPE GNO[2];
    int part_y;
} Zoltan_Arc;

/*--------------
 * FUNCTIONS
 */

/* Build a matrix from the user data get by the queries.
 * The matrix is only a translation of user data from GID space to GNO space.
 * The matrix cannot be directly used as a graph or hypergraph, you have
 * to make a distribution of this matrix. */
int
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix_options *opt, Zoltan_matrix* matrix,
  int, int, ZOLTAN_ID_PTR, ZOLTAN_GNO_TYPE *);

/* Free a matrix object */
void
Zoltan_Matrix_Free(Zoltan_matrix *m, int delete_flag);

/* Free a matrix2d object */
void
Zoltan_Matrix2d_Free(Zoltan_matrix_2d *m);

void
Zoltan_Matrix2d_Init(Zoltan_matrix_2d *m);

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
Zoltan_Matrix_Permute(ZZ* zz, Zoltan_matrix *m, ZOLTAN_GNO_TYPE * perm_y);

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

int Zoltan_Distribute_Set(Zoltan_matrix_2d* mat,
			  distFnct *hashDistFct, void * hashDistData);

void* Zoltan_Distribute_Partition_Register(ZZ* zz, int size, ZOLTAN_GNO_TYPE * yGNO, int *part, int nProc, int nPart);
void Zoltan_Distribute_Partition_Free(void** dist);


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

int
Zoltan_Matrix_Construct_CSR(ZZ *zz, int size, Zoltan_Arc *arcs, float* pinwgt,
			    Zoltan_matrix *outmat, int offset);


/* This code has to be called just before specializing the matrix into
 * a graph or an hypergraph.
 * Warning: Matrix cannot be modified afterwards.
 */
int
Zoltan_Matrix_Complete(ZZ* zz, Zoltan_matrix* m);

/* Return an array of locally owned GID */
ZOLTAN_ID_PTR Zoltan_Matrix_Get_GID(ZZ* zz, Zoltan_matrix* m);

int
Zoltan_Matrix_Vertex_Info(ZZ* zz, const Zoltan_matrix * const m,
			  ZOLTAN_ID_PTR lid,
			  float *wwgt, int *input_part);

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
