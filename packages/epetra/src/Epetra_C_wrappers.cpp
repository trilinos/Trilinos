
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */



#include <iostream.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#ifdef EPETRA_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

//#include "paz_aztec.h"
#include "Epetra_Epetra.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_LocalBlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_DVBR_Matrix.h"
#include "Epetra_DCRS_Matrix.h"
#include "Epetra_c_wrappers.h"

// catch memory errors
static void freeStoreException()
{
    cout << endl; // flush standard output
    Epetra_error("Epetra: Out of memory", 0);
}

/////////////////////////////////////////////////////
//                  Epetra_DVBR_Matrix         //
///////////////////////////////////////////////////


EPETRA_MATRIX_POINTER MANGLE(petra_rdp_dvbr_matrix_create)
                   (EPETRA_MAP rowmap)
{
  Epetra_BlockMap& rowmap_ = *(Epetra_BlockMap *) rowmap;
  Epetra_DVBR_Matrix *B = new Epetra_DVBR_Matrix(rowmap_);
  return((EPETRA_MATRIX_POINTER) B);
}

int MANGLE(petra_rdp_dvbr_matrix_allocate)
          (EPETRA_MATRIX A, int* numNzBlks, int* blkColInds)
{
  Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
  return(B->allocate(numNzBlks, blkColInds));
}
int MANGLE(petra_rdp_dvbr_matrix_putblockrow)
          (EPETRA_MATRIX A, EPETRA_INT blk_row, EPETRA_INT num_nz_blocks, 
				  double* vals, int* blk_col_inds)
{
  Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
  return(B->putBlockRow( DEREF_ blk_row, DEREF_ num_nz_blocks, vals, 
                         blk_col_inds));
}

int MANGLE(petra_rdp_dvbr_matrix_fillcomplete)(EPETRA_MATRIX A)
{
  Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
  return(B->fillComplete());
}

int MANGLE(petra_rdp_dvbr_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x,
                                             EPETRA_VECTOR y)
{
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    const Epetra_Vector& x1 = *(Epetra_Vector *) x;
    Epetra_Vector& y1 = *(Epetra_Vector *) y;
    return(B->matvec(x1, y1));
}

int MANGLE(petra_rdp_dvbr_matrix_matmultivec)(EPETRA_MATRIX A,
                                          EPETRA_MULTIVECTOR x,
                                          EPETRA_MULTIVECTOR y)
{
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    const Epetra_MultiVector& x1 = *(Epetra_MultiVector *) x;
    Epetra_MultiVector& y1 = *(Epetra_MultiVector *) y;
    return(B->matvec(x1, y1));
}

void MANGLE(petra_rdp_dvbr_matrix_destroy)(EPETRA_MATRIX A)
{
    delete (Epetra_DVBR_Matrix *) A;
}

/////////////////////////////////////////////////////
//                  Epetra_DCRS_Matrix         //
///////////////////////////////////////////////////


EPETRA_MATRIX_POINTER MANGLE(petra_rdp_dcrs_matrix_create) (EPETRA_MAP rowmap)
{
  Epetra_Map& rowmap_ = *(Epetra_Map *) rowmap;
  Epetra_DCRS_Matrix *B = new Epetra_DCRS_Matrix(rowmap_);
  return((EPETRA_MATRIX_POINTER) B);
}

int MANGLE(petra_rdp_dcrs_matrix_allocate)
          (EPETRA_MATRIX A, int* rowLengths)
{
  Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
  return(B->allocate(rowLengths));
}
int MANGLE(petra_rdp_dcrs_matrix_putrow)(EPETRA_MATRIX A, EPETRA_INT row,
EPETRA_INT num_nz, 
				  double* vals, int* col_inds)
{
  Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
  return(B->putRow( DEREF_ row, DEREF_ num_nz, vals, col_inds));
}

int MANGLE(petra_rdp_dcrs_matrix_sumintodiagonal)
          (EPETRA_MATRIX A, double* diagonal)
{
  Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
  return(B->sumIntoDiagonal( diagonal));
}

int MANGLE(petra_rdp_dcrs_matrix_fillcomplete)(EPETRA_MATRIX A)
{
  Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
  return(B->fillComplete());
}

int MANGLE(petra_rdp_dcrs_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x,
                                             EPETRA_VECTOR y)
{
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    const Epetra_Vector& x1 = *(Epetra_Vector *) x;
    Epetra_Vector& y1 = *(Epetra_Vector *) y;
    return(B->matvec(x1, y1));
}

int MANGLE(petra_rdp_dcrs_matrix_matmultivec)(EPETRA_MATRIX A,
                                          EPETRA_MULTIVECTOR x,
                                          EPETRA_MULTIVECTOR y)
{
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    const Epetra_MultiVector& x1 = *(Epetra_MultiVector *) x;
    Epetra_MultiVector& y1 = *(Epetra_MultiVector *) y;
    return(B->matvec(x1, y1));
}


void MANGLE(petra_rdp_dcrs_matrix_destroy)(EPETRA_MATRIX A)
{
    delete (Epetra_DCRS_Matrix *) A;
}

/////////////////////////////////////////////////////
//                  Epetra_Vector                  //
///////////////////////////////////////////////////

EPETRA_VECTOR_POINTER MANGLE(petra_rdp_vector_create)(EPETRA_MAP map)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_Vector *vector = new Epetra_Vector(map_);
  return((EPETRA_VECTOR_POINTER ) vector);
}

int MANGLE(petra_rdp_vector_putvector)(EPETRA_VECTOR x, double *vector)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  const double * t = (const double *) vector;
  return(x_->putVector(t));
}

int MANGLE(petra_rdp_vector_putscalar)(EPETRA_VECTOR x, EPETRA_DOUBLE scalar)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  return(x_->putScalar(DEREF_ scalar));
}

int MANGLE(petra_rdp_vector_norm1)(EPETRA_VECTOR x, double *scalar)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  return(x_->norm1(scalar));
}

int MANGLE(petra_rdp_vector_norm2)(EPETRA_VECTOR x, double *scalar)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  return(x_->norm2(scalar));
}

int MANGLE(petra_rdp_vector_random)(EPETRA_VECTOR x)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  return(x_->random());
}

int MANGLE(petra_rdp_vector_lincomb)(EPETRA_VECTOR x, EPETRA_VECTOR b, EPETRA_DOUBLE
				 scalar, EPETRA_VECTOR c)
{
  Epetra_Vector *x_ = (Epetra_Vector *) x;
  Epetra_Vector& b_ = *(Epetra_Vector *) b;
  Epetra_Vector& c_ = *(Epetra_Vector *) c;
  return(x_->linComb(b_, DEREF_ scalar, c_));
}

void MANGLE(petra_rdp_vector_destroy)(EPETRA_VECTOR x)
{
    delete (Epetra_Vector *) x;
}

/////////////////////////////////////////////////////
//                  Epetra_MultiVector         //
///////////////////////////////////////////////////

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create)()
{
  Epetra_MultiVector *vector = new Epetra_MultiVector();
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create1)
                         (EPETRA_MAP map, EPETRA_INT numVectors)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_MultiVector *vector = new Epetra_MultiVector(map_, DEREF_ numVectors);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create2)(EPETRA_MAP map, 
                                double *A, EPETRA_INT lda, EPETRA_INT numVectors)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_MultiVector *vector = new Epetra_MultiVector(map_, A, DEREF_ lda,
          DEREF_ numVectors);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create3)(EPETRA_MAP map, 
                                double **in_multiVector, EPETRA_INT numVectors)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_MultiVector *vector = new Epetra_MultiVector(map_, in_multiVector,
                                                 DEREF_ numVectors);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create4)
                         (EPETRA_MULTIVECTOR in_multiVector)
{
  Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
  Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create5)(EPETRA_MULTIVECTOR 
                        in_multiVector, EPETRA_INT numVectors, int *vecIndices)
{
  Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
  Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_,
                                   DEREF_ numVectors, vecIndices);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create6)(EPETRA_MULTIVECTOR
                        in_multiVector, EPETRA_INT startindex, EPETRA_INT numvectors)
{
  Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
  Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_, DEREF_ startindex,
                                   DEREF_ numvectors);
  return((EPETRA_MULTIVECTOR_POINTER) vector);
}

int MANGLE(petra_rdp_multivector_putmultivector)
          (EPETRA_MULTIVECTOR multiVector, 
                                       double **in_multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  const double ** t = (const double **) in_multiVector;
  return(multiVector_->putMultiVector(t));
}

int MANGLE(petra_rdp_multivector_allocate)(EPETRA_MULTIVECTOR multiVector, 
                                EPETRA_MAP map, EPETRA_INT numVectors)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->allocate(map_, DEREF_ numVectors));
}

int MANGLE(petra_rdp_multivector_putscalar)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->putScalar(DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_scale)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->scale(DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_scalecopy)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
           EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->scaleCopy(multiVector_in_, DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_dotprod)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
           double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->dotProd(multiVector_in_, scalar));
}

int MANGLE(petra_rdp_multivector_addvec)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar, 
           EPETRA_MULTIVECTOR multiVector_in)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->addVec(DEREF_ scalar, multiVector_in_));
}

int MANGLE(petra_rdp_multivector_norm1)
          (EPETRA_MULTIVECTOR multiVector, double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->norm1(scalar));
}

int MANGLE(petra_rdp_multivector_norm2)
          (EPETRA_MULTIVECTOR multiVector, double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->norm2(scalar));
}

int MANGLE(petra_rdp_multivector_lincomb)(EPETRA_MULTIVECTOR multiVector,
                               EPETRA_MULTIVECTOR b, 
                               EPETRA_DOUBLE scalar, EPETRA_MULTIVECTOR c)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& b_ = *(Epetra_MultiVector *) b;
  Epetra_MultiVector& c_ = *(Epetra_MultiVector *) c;
  return(multiVector_->linComb(b_,DEREF_ scalar,c_));
}

int MANGLE(petra_rdp_multivector_random)
          (EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->random());
}

int MANGLE(petra_rdp_multivector_reduce)(EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->reduce());
}

int MANGLE(petra_rdp_multivector_numvectors)(EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->numVectors());
}

int MANGLE(petra_rdp_multivector_gemm)(EPETRA_MULTIVECTOR multiVector,
                   EPETRA_INT transa, EPETRA_INT transb, EPETRA_DOUBLE alpha,
                   EPETRA_MULTIVECTOR A, EPETRA_MULTIVECTOR B,
                   EPETRA_DOUBLE beta )
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& A_ = *(Epetra_MultiVector *) A;
  Epetra_MultiVector& B_ = *(Epetra_MultiVector *) B;
  bool transa_ = !(DEREF_ transa==0);
  bool transb_ = !(DEREF_ transb==0);
  return(multiVector_->GEMM(transa_, transb_, DEREF_ alpha, A_, B_, DEREF_ beta));
}

void MANGLE(petra_rdp_multivector_destroy)(EPETRA_MULTIVECTOR multiVector)
{
    delete (Epetra_MultiVector *) multiVector;
}

/////////////////////////////////////////////////////
//                  Epetra_Map                     //
///////////////////////////////////////////////////

EPETRA_MAP_POINTER MANGLE(petra_map_create1)(EPETRA_INT numGlobalEquations,
                                            EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_Map *map = new Epetra_Map(DEREF_ numGlobalEquations, comm_);
  return((EPETRA_MAP_POINTER ) map);
}


EPETRA_MAP_POINTER MANGLE(petra_map_create2)(EPETRA_INT numGlobalEquations,
					   EPETRA_INT numLocalEquations,
					   int *updateList, 
					   EPETRA_INT indexBase,
					   EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_Map *map = new Epetra_Map(DEREF_ numGlobalEquations, DEREF_ numLocalEquations,
                                updateList, DEREF_ indexBase, comm_);
  return((EPETRA_MAP_POINTER ) map);
}

int MANGLE(petra_map_numlocalequations)(EPETRA_MAP map)
{
  Epetra_Map * map_ = (Epetra_Map *) map;
  return(map_->numLocalEquations());
}

#ifndef EPETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
int * MANGLE(petra_map_getupdatelist)(EPETRA_MAP map)
{
  Epetra_Map * map_ = (Epetra_Map *) map;
  return(map_->getUpdateList());
}
#endif
EPETRA_COMM_POINTER MANGLE(petra_map_getcommunicator)(EPETRA_MAP map)
{
  Epetra_Map * map_ = (Epetra_Map *) map;
  return((EPETRA_COMM_POINTER) &(map_->getCommunicator()));
}

void MANGLE(petra_map_destroy)(EPETRA_MAP map)
{
    delete (Epetra_Map *) map;
}

/////////////////////////////////////////////////////
//                  Epetra_BlockMap                //
///////////////////////////////////////////////////

EPETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create1)(
		      EPETRA_INT numGlobalEquations, EPETRA_INT numLocalEquations, int *updateList,
		      EPETRA_INT numGlobalBlocks, EPETRA_INT numLocalBlocks, 
		      int *blockUpdateList,
		      int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_BlockMap *blockmap = new Epetra_BlockMap(DEREF_ numGlobalEquations,
                           DEREF_ numLocalEquations, updateList,
                           DEREF_ numGlobalBlocks, DEREF_ numLocalBlocks,
                           blockUpdateList,
                           blockSizes, DEREF_ indexBase, comm_);
  return((EPETRA_BLOCKMAP_POINTER ) blockmap);
}

EPETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create2)(
		      EPETRA_INT numGlobalBlocks, EPETRA_INT numLocalBlocks, 
		      int *blockUpdateList,
		      int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_BlockMap *blockmap = new Epetra_BlockMap(
                           DEREF_ numGlobalBlocks, DEREF_ numLocalBlocks, 
			   blockUpdateList,
                           blockSizes, DEREF_ indexBase, comm_);
  return((EPETRA_BLOCKMAP_POINTER ) blockmap);
}

void MANGLE(petra_blockmap_destroy)(EPETRA_BLOCKMAP blockmap)
{
    delete (Epetra_BlockMap *) blockmap;
}

/////////////////////////////////////////////////////
//                  Epetra_LocalMap                //
///////////////////////////////////////////////////

EPETRA_LOCALMAP_POINTER MANGLE(petra_localmap_create)(EPETRA_INT numLocalEquations,
                               EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalMap *localmap = new Epetra_LocalMap(DEREF_ numLocalEquations,
                    DEREF_ indexBase, comm_);
  return((EPETRA_LOCALMAP_POINTER ) localmap);
}

void MANGLE(petra_localmap_destroy)(EPETRA_LOCALMAP localmap)
{
    delete (Epetra_LocalMap *) localmap;
}

/////////////////////////////////////////////////////
//                  Epetra_LocalBlockMap           //
///////////////////////////////////////////////////

EPETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create1)(
		      EPETRA_INT numLocalEquations,
		      EPETRA_INT numLocalBlocks,
		      int* blockSizes,
		      EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalBlockMap *localblockmap = new
                      Epetra_LocalBlockMap(DEREF_ numLocalEquations,
                                  DEREF_ numLocalBlocks,
                                  blockSizes,
                                  DEREF_ indexBase, comm_);
  return((EPETRA_LOCALBLOCKMAP_POINTER ) localblockmap);
}

EPETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create2)(
		      EPETRA_INT numLocalBlocks,
		      int* blockSizes,
		      EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalBlockMap *localblockmap = new
                      Epetra_LocalBlockMap(DEREF_ numLocalBlocks,
                                  blockSizes,
                                  DEREF_ indexBase, comm_);
  return((EPETRA_LOCALBLOCKMAP_POINTER ) localblockmap);
}

void MANGLE(petra_localblockmap_destroy)(EPETRA_LOCALBLOCKMAP localblockmap)
{
    delete (Epetra_LocalBlockMap *) localblockmap;
}

/////////////////////////////////////////////////////
//                  Epetra_Comm                    //
///////////////////////////////////////////////////

#ifdef EPETRA_MPI
EPETRA_COMM_POINTER MANGLE(petra_comm_create)(MPI_Comm MPI_communicator)
{
  Epetra_Comm *comm = new Epetra_Comm(MPI_communicator);
  return((EPETRA_COMM_POINTER ) comm);
}
#else
EPETRA_COMM_POINTER MANGLE(petra_comm_create_serial)()
{
  Epetra_Comm *comm = new Epetra_Comm();
  return((EPETRA_COMM_POINTER ) comm);
}
#endif

int MANGLE(petra_comm_getmypid)(EPETRA_COMM communicator)
{
  Epetra_Comm *communicator_ = (Epetra_Comm *) communicator;
  return(communicator_->getMyPID());
  
}
int MANGLE(petra_comm_getnumproc)(EPETRA_COMM communicator)
{
  Epetra_Comm* communicator_ = (Epetra_Comm *) communicator;
  return(communicator_->getNumProc());
  
}

void MANGLE(petra_comm_barrier)(EPETRA_COMM communicator)
{
  Epetra_Comm* communicator_ = (Epetra_Comm *) communicator;
  communicator_->barrier();
  
}

void MANGLE(petra_comm_destroy)(EPETRA_COMM communicator)
{
    delete (Epetra_Comm *) communicator;

}


#ifdef __cplusplus
}
#endif
