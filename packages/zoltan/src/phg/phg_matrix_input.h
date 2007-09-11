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

#ifndef __PHG_MATRIX_INPUT_H
#define __PHG_MATRIX_INPUT_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "params_const.h" /* for MAX_PARAM_STRING_LEN */

/********************************************************************
 * Data structures used when converting an input sparse matrix
 * to a PHG problem, running the PHG problem and retaining results
 * for user query.  
 *
 * For performance reasons, we'll use ints for matrix indices and
 * sizes, but for really large matrices these may someday be longs.
 * So typedef the row and column indices and matrix dimensions.
 *
 * If IJTYPE is changed, change the parameters of the
 * CSC and CSR query functions in zoltan.h and the User's Guide.
 *
 * This should be moved out of phg to somewhere more general.  The
 * Zoltan algorithm used for the sparse matrix need not be PHG.
 ********************************************************************/

typedef unsigned int IJTYPE; /* matrix row ID, column ID or matrix size */

enum ObjectType {ROW_TYPE = 1, COL_TYPE = 2, NZ_TYPE = 3};

enum PartitionType {MP_ROW_TYPE=1, MP_COLUMN_TYPE=2, MP_GENERAL_TYPE=3};

struct obj_node{
  IJTYPE i;
  IJTYPE j;
  IJTYPE objLID;
  struct obj_node *next;
};
typedef struct _obj_lookup{
  struct obj_node *htTop;
  struct obj_node **ht;
  IJTYPE table_size;
  int    key_size;
} obj_lookup;

struct Zoltan_MP_Data_Struct{
  struct Zoltan_Struct *zzLib;    /* Problem created by Zoltan_Matrix_Partition() */

  /* Parameters */
  int approach;       /* a PartitionType, the LB_APPROACH parameter */
  char method[MAX_PARAM_STRING_LEN]; /* string, partitioning method */

  /* The local portion of sparse matrix returned by the query function */
  int input_type;    /* a ObjectType, how they supply the matrix (CSC or CSR) */
  IJTYPE numRC;      /* number of rows or columns */
  IJTYPE *rcGID;     /* row or column GIDs   */
  IJTYPE *pinIndex;  /* index into pinGIDs array, last is num pins */
  IJTYPE *pinGID;    /* non-zeroes column or row GIDs */

  /* Mirror specification of sparse matrix: if input was CSR, create CSC, 
   * or in input was CSC, create CSR */

  IJTYPE numCR;  
  IJTYPE *crGID;
  IJTYPE *mirrorPinIndex;
  IJTYPE *mirrorPinGID; 

  /* Global values filled out by process_matrix_input().                  */
  IJTYPE rowBaseID;   /* most likely zero or one */
  IJTYPE colBaseID;   /* most likely zero or one */
  IJTYPE nRows;
  IJTYPE nCols;
  IJTYPE nNonZeros;

  /* Hypergraph generated from sparse matrix (if not obvious 
   * from sparse matrix representation) */

  IJTYPE nMyVtx;    /* my number of vertices in hypergraph */
  IJTYPE *vtxGID;   /* vertex GIDs */
  double *vtxWgt;   /* weight for each vertex (1 double) */
  IJTYPE nHedges;   /* number of hyperedges */
  IJTYPE *hindex;   /* index into list of pins for each h.e., last is npins */
  IJTYPE *hvertex;  /* vtx GID of pins in my hyperedges */

  /* Results, to supply data to query functions */
  int *rowproc;
  int *rowpart;
  obj_lookup *row_lookup;
  int *colproc;
  int *colpart;
  obj_lookup *col_lookup;
  int *pinproc;
  int *pinpart;
  obj_lookup *pin_lookup;
};

typedef struct Zoltan_MP_Data_Struct ZOLTAN_MP_DATA;

int Zoltan_MP_Get_NonZero_Assignment(struct Zoltan_Struct *zz, int nNZ,
        IJTYPE *rowIDs, IJTYPE *colIDs, int *nzProcs, int *nzParts);
int Zoltan_MP_Get_Column_Assignment(struct Zoltan_Struct *zz, int nCols, IJTYPE *colIDs,
        int *colProcs, int *colParts);
int Zoltan_MP_Get_Row_Assignment(struct Zoltan_Struct *zz, int nRows, IJTYPE *rowIDs,
        int *rowProcs, int *rowParts);
int Zoltan_Lookup_Obj(obj_lookup *lu, IJTYPE I, IJTYPE J);

void Zoltan_MP_Debug_Partitioning(struct Zoltan_Struct *zz);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __PHG_MATRIX_INPUT_H */
