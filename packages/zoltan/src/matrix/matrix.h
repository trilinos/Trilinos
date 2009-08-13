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

#ifndef __MATRIX_H
#define __MATRIX_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg_comm.h"
#include "zoltan_dd.h"

  /* Would be a lot better in C++ ! */


/* This structure is a CS view of a part of the matrix/hypergraph */
typedef struct Zoltan_matrix_ {
  int           transpose;   /* Need to transpose to have a CSC view ? */
  int           enforceSquare;
  int           completed;
  int           globalX;   /* Overall number of objects */
  int           globalY;
  int           offsetY;     /* Used for bipartite graph: GNO >= offsetY are edges */
  int           nY;
  int           ywgtdim;
  int           nPins;
  int          *yGNO;       /* Local edges gnos */
  int          *ystart;     /* Indirection array to describe a column */
  int          *yend;       /* end of local pins, usually ystart+1 
			       (and is ystart+1 after matrix complete) */
  int          *pinGNO;     /* array of gno of other extremtiy */
  float        *pinwgt;

  /* These fields are used only before matrix_complete */
  /* Allow us to move only pins and CSR structure without having to worry
   * about vertex and edge data. */
  struct Zoltan_DD_Struct *ddX; /* Map xGNO -> xGID, xwgt, Input_Parts */
  struct Zoltan_DD_Struct *ddY; /* Map yGNO -> yGID, ywgt */

  /* These fields are used after matrix_complete */
  float        *ywgt;
  ZOLTAN_ID_PTR yGID;
} Zoltan_matrix;

typedef struct Zoltan_matrix_2d_ {
  Zoltan_matrix   mtx;
  PHGComm         *comm;
  int             *dist_x;
  int             *dist_y;
} Zoltan_matrix_2d;

int
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix* matrix);

void
Zoltan_Matrix_Free(ZZ *zz, Zoltan_matrix *m);

void
Zoltan_Matrix2d_Free(ZZ *zz, Zoltan_matrix_2d *m);

void
Zoltan_Matrix_Reset(Zoltan_matrix* m);


/* if !copy, inmat is not usable after this call */
int
Zoltan_Matrix2d_Distribute (ZZ* zz, const Zoltan_matrix inmat,
			    Zoltan_matrix_2d *outmat, int copy);

int
Zoltan_Distribute_layout (ZZ *zz, const PHGComm * const inlayout,
			  int hiRank, int loRank,
			  int reqx, int reqy,
			  PHGComm *outlayout);

int Zoltan_Distribute_Square (ZZ * zz, PHGComm *layout) ;
int Zoltan_Distribute_LinearY (ZZ * zz, PHGComm *layout) ;

int
Zoltan_Matrix_Bipart(ZZ* zz, Zoltan_matrix *matrix, int nProc, int myProc);

int
Zoltan_Matrix_Complete(ZZ* zz,Zoltan_matrix* m);

int
Zoltan_Matrix2d_adjproc (ZZ* zz, const Zoltan_matrix_2d * const mat, int **adjproc);

#ifdef __cplusplus
}
#endif

#endif
