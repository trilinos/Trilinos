

#include <iostream.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#ifdef PETRA_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

//#include "paz_aztec.h"
#include "Petra_Petra.h"
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_LocalMap.h"
#include "Petra_BlockMap.h"
#include "Petra_LocalBlockMap.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_DVBR_Matrix.h"
#include "Petra_RDP_DCRS_Matrix.h"
#include "Petra_c_wrappers.h"

// catch memory errors
static void freeStoreException()
{
    cout << endl; // flush standard output
    Petra_error("Petra: Out of memory", 0);
}

/////////////////////////////////////////////////////
//                  Petra_RDP_DVBR_Matrix         //
///////////////////////////////////////////////////


PETRA_MATRIX_POINTER MANGLE(petra_rdp_dvbr_matrix_create)
                   (PETRA_MAP rowmap)
{
  Petra_BlockMap& rowmap_ = *(Petra_BlockMap *) rowmap;
  Petra_RDP_DVBR_Matrix *B = new Petra_RDP_DVBR_Matrix(rowmap_);
  return((PETRA_MATRIX_POINTER) B);
}

int MANGLE(petra_rdp_dvbr_matrix_allocate)
          (PETRA_MATRIX A, int* numNzBlks, int* blkColInds)
{
  Petra_RDP_DVBR_Matrix *B = (Petra_RDP_DVBR_Matrix *) A;
  return(B->allocate(numNzBlks, blkColInds));
}
int MANGLE(petra_rdp_dvbr_matrix_putblockrow)
          (PETRA_MATRIX A, PETRA_INT blk_row, PETRA_INT num_nz_blocks, 
				  double* vals, int* blk_col_inds)
{
  Petra_RDP_DVBR_Matrix *B = (Petra_RDP_DVBR_Matrix *) A;
  return(B->putBlockRow( DEREF_ blk_row, DEREF_ num_nz_blocks, vals, 
                         blk_col_inds));
}

int MANGLE(petra_rdp_dvbr_matrix_fillcomplete)(PETRA_MATRIX A)
{
  Petra_RDP_DVBR_Matrix *B = (Petra_RDP_DVBR_Matrix *) A;
  return(B->fillComplete());
}

int MANGLE(petra_rdp_dvbr_matrix_matvec)(PETRA_MATRIX A, PETRA_VECTOR x,
                                             PETRA_VECTOR y)
{
    Petra_RDP_DVBR_Matrix *B = (Petra_RDP_DVBR_Matrix *) A;
    const Petra_RDP_Vector& x1 = *(Petra_RDP_Vector *) x;
    Petra_RDP_Vector& y1 = *(Petra_RDP_Vector *) y;
    return(B->matvec(x1, y1));
}

int MANGLE(petra_rdp_dvbr_matrix_matmultivec)(PETRA_MATRIX A,
                                          PETRA_MULTIVECTOR x,
                                          PETRA_MULTIVECTOR y)
{
    Petra_RDP_DVBR_Matrix *B = (Petra_RDP_DVBR_Matrix *) A;
    const Petra_RDP_MultiVector& x1 = *(Petra_RDP_MultiVector *) x;
    Petra_RDP_MultiVector& y1 = *(Petra_RDP_MultiVector *) y;
    return(B->matvec(x1, y1));
}

void MANGLE(petra_rdp_dvbr_matrix_destroy)(PETRA_MATRIX A)
{
    delete (Petra_RDP_DVBR_Matrix *) A;
}

/////////////////////////////////////////////////////
//                  Petra_RDP_DCRS_Matrix         //
///////////////////////////////////////////////////


PETRA_MATRIX_POINTER MANGLE(petra_rdp_dcrs_matrix_create) (PETRA_MAP rowmap)
{
  Petra_Map& rowmap_ = *(Petra_Map *) rowmap;
  Petra_RDP_DCRS_Matrix *B = new Petra_RDP_DCRS_Matrix(rowmap_);
  return((PETRA_MATRIX_POINTER) B);
}

int MANGLE(petra_rdp_dcrs_matrix_allocate)
          (PETRA_MATRIX A, int* rowLengths)
{
  Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
  return(B->allocate(rowLengths));
}
int MANGLE(petra_rdp_dcrs_matrix_putrow)(PETRA_MATRIX A, PETRA_INT row,
PETRA_INT num_nz, 
				  double* vals, int* col_inds)
{
  Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
  return(B->putRow( DEREF_ row, DEREF_ num_nz, vals, col_inds));
}

int MANGLE(petra_rdp_dcrs_matrix_sumintodiagonal)
          (PETRA_MATRIX A, double* diagonal)
{
  Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
  return(B->sumIntoDiagonal( diagonal));
}

int MANGLE(petra_rdp_dcrs_matrix_fillcomplete)(PETRA_MATRIX A)
{
  Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
  return(B->fillComplete());
}

int MANGLE(petra_rdp_dcrs_matrix_matvec)(PETRA_MATRIX A, PETRA_VECTOR x,
                                             PETRA_VECTOR y)
{
    Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
    const Petra_RDP_Vector& x1 = *(Petra_RDP_Vector *) x;
    Petra_RDP_Vector& y1 = *(Petra_RDP_Vector *) y;
    return(B->matvec(x1, y1));
}

int MANGLE(petra_rdp_dcrs_matrix_matmultivec)(PETRA_MATRIX A,
                                          PETRA_MULTIVECTOR x,
                                          PETRA_MULTIVECTOR y)
{
    Petra_RDP_DCRS_Matrix *B = (Petra_RDP_DCRS_Matrix *) A;
    const Petra_RDP_MultiVector& x1 = *(Petra_RDP_MultiVector *) x;
    Petra_RDP_MultiVector& y1 = *(Petra_RDP_MultiVector *) y;
    return(B->matvec(x1, y1));
}


void MANGLE(petra_rdp_dcrs_matrix_destroy)(PETRA_MATRIX A)
{
    delete (Petra_RDP_DCRS_Matrix *) A;
}

/////////////////////////////////////////////////////
//                  Petra_RDP_Vector                  //
///////////////////////////////////////////////////

PETRA_VECTOR_POINTER MANGLE(petra_rdp_vector_create)(PETRA_MAP map)
{
  Petra_Map& map_ = *(Petra_Map *) map;
  Petra_RDP_Vector *vector = new Petra_RDP_Vector(map_);
  return((PETRA_VECTOR_POINTER ) vector);
}

int MANGLE(petra_rdp_vector_putvector)(PETRA_VECTOR x, double *vector)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  const double * t = (const double *) vector;
  return(x_->putVector(t));
}

int MANGLE(petra_rdp_vector_putscalar)(PETRA_VECTOR x, PETRA_DOUBLE scalar)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  return(x_->putScalar(DEREF_ scalar));
}

int MANGLE(petra_rdp_vector_norm1)(PETRA_VECTOR x, double *scalar)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  return(x_->norm1(scalar));
}

int MANGLE(petra_rdp_vector_norm2)(PETRA_VECTOR x, double *scalar)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  return(x_->norm2(scalar));
}

int MANGLE(petra_rdp_vector_random)(PETRA_VECTOR x)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  return(x_->random());
}

int MANGLE(petra_rdp_vector_lincomb)(PETRA_VECTOR x, PETRA_VECTOR b, PETRA_DOUBLE
				 scalar, PETRA_VECTOR c)
{
  Petra_RDP_Vector *x_ = (Petra_RDP_Vector *) x;
  Petra_RDP_Vector& b_ = *(Petra_RDP_Vector *) b;
  Petra_RDP_Vector& c_ = *(Petra_RDP_Vector *) c;
  return(x_->linComb(b_, DEREF_ scalar, c_));
}

void MANGLE(petra_rdp_vector_destroy)(PETRA_VECTOR x)
{
    delete (Petra_RDP_Vector *) x;
}

/////////////////////////////////////////////////////
//                  Petra_RDP_MultiVector         //
///////////////////////////////////////////////////

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create)()
{
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector();
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create1)
                         (PETRA_MAP map, PETRA_INT numVectors)
{
  Petra_Map& map_ = *(Petra_Map *) map;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(map_, DEREF_ numVectors);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create2)(PETRA_MAP map, 
                                double *A, PETRA_INT lda, PETRA_INT numVectors)
{
  Petra_Map& map_ = *(Petra_Map *) map;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(map_, A, DEREF_ lda,
          DEREF_ numVectors);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create3)(PETRA_MAP map, 
                                double **in_multiVector, PETRA_INT numVectors)
{
  Petra_Map& map_ = *(Petra_Map *) map;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(map_, in_multiVector,
                                                 DEREF_ numVectors);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create4)
                         (PETRA_MULTIVECTOR in_multiVector)
{
  Petra_RDP_MultiVector & in_multiVector_ = *(Petra_RDP_MultiVector *) in_multiVector;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(in_multiVector_);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create5)(PETRA_MULTIVECTOR 
                        in_multiVector, PETRA_INT numVectors, int *vecIndices)
{
  Petra_RDP_MultiVector & in_multiVector_ = *(Petra_RDP_MultiVector *) in_multiVector;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(in_multiVector_,
                                   DEREF_ numVectors, vecIndices);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create6)(PETRA_MULTIVECTOR
                        in_multiVector, PETRA_INT startindex, PETRA_INT numvectors)
{
  Petra_RDP_MultiVector & in_multiVector_ = *(Petra_RDP_MultiVector *) in_multiVector;
  Petra_RDP_MultiVector *vector = new Petra_RDP_MultiVector(in_multiVector_, DEREF_ startindex,
                                   DEREF_ numvectors);
  return((PETRA_MULTIVECTOR_POINTER) vector);
}

int MANGLE(petra_rdp_multivector_putmultivector)
          (PETRA_MULTIVECTOR multiVector, 
                                       double **in_multiVector)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  const double ** t = (const double **) in_multiVector;
  return(multiVector_->putMultiVector(t));
}

int MANGLE(petra_rdp_multivector_allocate)(PETRA_MULTIVECTOR multiVector, 
                                PETRA_MAP map, PETRA_INT numVectors)
{
  Petra_Map& map_ = *(Petra_Map *) map;
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->allocate(map_, DEREF_ numVectors));
}

int MANGLE(petra_rdp_multivector_putscalar)
          (PETRA_MULTIVECTOR multiVector, PETRA_DOUBLE scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->putScalar(DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_scale)
          (PETRA_MULTIVECTOR multiVector, PETRA_DOUBLE scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->scale(DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_scalecopy)
          (PETRA_MULTIVECTOR multiVector, PETRA_MULTIVECTOR multiVector_in,
           PETRA_DOUBLE scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  Petra_RDP_MultiVector& multiVector_in_ = *(Petra_RDP_MultiVector *) multiVector_in;
  return(multiVector_->scaleCopy(multiVector_in_, DEREF_ scalar));
}

int MANGLE(petra_rdp_multivector_dotprod)
          (PETRA_MULTIVECTOR multiVector, PETRA_MULTIVECTOR multiVector_in,
           double *scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  Petra_RDP_MultiVector& multiVector_in_ = *(Petra_RDP_MultiVector *) multiVector_in;
  return(multiVector_->dotProd(multiVector_in_, scalar));
}

int MANGLE(petra_rdp_multivector_addvec)
          (PETRA_MULTIVECTOR multiVector, PETRA_DOUBLE scalar, 
           PETRA_MULTIVECTOR multiVector_in)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  Petra_RDP_MultiVector& multiVector_in_ = *(Petra_RDP_MultiVector *) multiVector_in;
  return(multiVector_->addVec(DEREF_ scalar, multiVector_in_));
}

int MANGLE(petra_rdp_multivector_norm1)
          (PETRA_MULTIVECTOR multiVector, double *scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->norm1(scalar));
}

int MANGLE(petra_rdp_multivector_norm2)
          (PETRA_MULTIVECTOR multiVector, double *scalar)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->norm2(scalar));
}

int MANGLE(petra_rdp_multivector_lincomb)(PETRA_MULTIVECTOR multiVector,
                               PETRA_MULTIVECTOR b, 
                               PETRA_DOUBLE scalar, PETRA_MULTIVECTOR c)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  Petra_RDP_MultiVector& b_ = *(Petra_RDP_MultiVector *) b;
  Petra_RDP_MultiVector& c_ = *(Petra_RDP_MultiVector *) c;
  return(multiVector_->linComb(b_,DEREF_ scalar,c_));
}

int MANGLE(petra_rdp_multivector_random)
          (PETRA_MULTIVECTOR multiVector)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->random());
}

int MANGLE(petra_rdp_multivector_reduce)(PETRA_MULTIVECTOR multiVector)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->reduce());
}

int MANGLE(petra_rdp_multivector_numvectors)(PETRA_MULTIVECTOR multiVector)
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  return(multiVector_->numVectors());
}

int MANGLE(petra_rdp_multivector_gemm)(PETRA_MULTIVECTOR multiVector,
                   PETRA_INT transa, PETRA_INT transb, PETRA_DOUBLE alpha,
                   PETRA_MULTIVECTOR A, PETRA_MULTIVECTOR B,
                   PETRA_DOUBLE beta )
{
  Petra_RDP_MultiVector *multiVector_ = (Petra_RDP_MultiVector *) multiVector;
  Petra_RDP_MultiVector& A_ = *(Petra_RDP_MultiVector *) A;
  Petra_RDP_MultiVector& B_ = *(Petra_RDP_MultiVector *) B;
  bool transa_ = !(DEREF_ transa==0);
  bool transb_ = !(DEREF_ transb==0);
  return(multiVector_->GEMM(transa_, transb_, DEREF_ alpha, A_, B_, DEREF_ beta));
}

void MANGLE(petra_rdp_multivector_destroy)(PETRA_MULTIVECTOR multiVector)
{
    delete (Petra_RDP_MultiVector *) multiVector;
}

/////////////////////////////////////////////////////
//                  Petra_Map                     //
///////////////////////////////////////////////////

PETRA_MAP_POINTER MANGLE(petra_map_create1)(PETRA_INT numGlobalEquations,
                                            PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_Map *map = new Petra_Map(DEREF_ numGlobalEquations, comm_);
  return((PETRA_MAP_POINTER ) map);
}


PETRA_MAP_POINTER MANGLE(petra_map_create2)(PETRA_INT numGlobalEquations,
					   PETRA_INT numLocalEquations,
					   int *updateList, 
					   PETRA_INT indexBase,
					   PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_Map *map = new Petra_Map(DEREF_ numGlobalEquations, DEREF_ numLocalEquations,
                                updateList, DEREF_ indexBase, comm_);
  return((PETRA_MAP_POINTER ) map);
}

int MANGLE(petra_map_numlocalequations)(PETRA_MAP map)
{
  Petra_Map * map_ = (Petra_Map *) map;
  return(map_->numLocalEquations());
}

#ifndef PETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
int * MANGLE(petra_map_getupdatelist)(PETRA_MAP map)
{
  Petra_Map * map_ = (Petra_Map *) map;
  return(map_->getUpdateList());
}
#endif
PETRA_COMM_POINTER MANGLE(petra_map_getcommunicator)(PETRA_MAP map)
{
  Petra_Map * map_ = (Petra_Map *) map;
  return((PETRA_COMM_POINTER) &(map_->getCommunicator()));
}

void MANGLE(petra_map_destroy)(PETRA_MAP map)
{
    delete (Petra_Map *) map;
}

/////////////////////////////////////////////////////
//                  Petra_BlockMap                //
///////////////////////////////////////////////////

PETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create1)(
		      PETRA_INT numGlobalEquations, PETRA_INT numLocalEquations, int *updateList,
		      PETRA_INT numGlobalBlocks, PETRA_INT numLocalBlocks, 
		      int *blockUpdateList,
		      int* blockSizes, PETRA_INT indexBase, PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_BlockMap *blockmap = new Petra_BlockMap(DEREF_ numGlobalEquations,
                           DEREF_ numLocalEquations, updateList,
                           DEREF_ numGlobalBlocks, DEREF_ numLocalBlocks,
                           blockUpdateList,
                           blockSizes, DEREF_ indexBase, comm_);
  return((PETRA_BLOCKMAP_POINTER ) blockmap);
}

PETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create2)(
		      PETRA_INT numGlobalBlocks, PETRA_INT numLocalBlocks, 
		      int *blockUpdateList,
		      int* blockSizes, PETRA_INT indexBase, PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_BlockMap *blockmap = new Petra_BlockMap(
                           DEREF_ numGlobalBlocks, DEREF_ numLocalBlocks, 
			   blockUpdateList,
                           blockSizes, DEREF_ indexBase, comm_);
  return((PETRA_BLOCKMAP_POINTER ) blockmap);
}

void MANGLE(petra_blockmap_destroy)(PETRA_BLOCKMAP blockmap)
{
    delete (Petra_BlockMap *) blockmap;
}

/////////////////////////////////////////////////////
//                  Petra_LocalMap                //
///////////////////////////////////////////////////

PETRA_LOCALMAP_POINTER MANGLE(petra_localmap_create)(PETRA_INT numLocalEquations,
                               PETRA_INT indexBase, PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_LocalMap *localmap = new Petra_LocalMap(DEREF_ numLocalEquations,
                    DEREF_ indexBase, comm_);
  return((PETRA_LOCALMAP_POINTER ) localmap);
}

void MANGLE(petra_localmap_destroy)(PETRA_LOCALMAP localmap)
{
    delete (Petra_LocalMap *) localmap;
}

/////////////////////////////////////////////////////
//                  Petra_LocalBlockMap           //
///////////////////////////////////////////////////

PETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create1)(
		      PETRA_INT numLocalEquations,
		      PETRA_INT numLocalBlocks,
		      int* blockSizes,
		      PETRA_INT indexBase, PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_LocalBlockMap *localblockmap = new
                      Petra_LocalBlockMap(DEREF_ numLocalEquations,
                                  DEREF_ numLocalBlocks,
                                  blockSizes,
                                  DEREF_ indexBase, comm_);
  return((PETRA_LOCALBLOCKMAP_POINTER ) localblockmap);
}

PETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create2)(
		      PETRA_INT numLocalBlocks,
		      int* blockSizes,
		      PETRA_INT indexBase, PETRA_COMM comm)
{
  Petra_Comm& comm_ = *(Petra_Comm *) comm;
  Petra_LocalBlockMap *localblockmap = new
                      Petra_LocalBlockMap(DEREF_ numLocalBlocks,
                                  blockSizes,
                                  DEREF_ indexBase, comm_);
  return((PETRA_LOCALBLOCKMAP_POINTER ) localblockmap);
}

void MANGLE(petra_localblockmap_destroy)(PETRA_LOCALBLOCKMAP localblockmap)
{
    delete (Petra_LocalBlockMap *) localblockmap;
}

/////////////////////////////////////////////////////
//                  Petra_Comm                    //
///////////////////////////////////////////////////

#ifdef PETRA_MPI
PETRA_COMM_POINTER MANGLE(petra_comm_create)(MPI_Comm MPI_communicator)
{
  Petra_Comm *comm = new Petra_Comm(MPI_communicator);
  return((PETRA_COMM_POINTER ) comm);
}
#else
PETRA_COMM_POINTER MANGLE(petra_comm_create_serial)()
{
  Petra_Comm *comm = new Petra_Comm();
  return((PETRA_COMM_POINTER ) comm);
}
#endif

int MANGLE(petra_comm_getmypid)(PETRA_COMM communicator)
{
  Petra_Comm *communicator_ = (Petra_Comm *) communicator;
  return(communicator_->getMyPID());
  
}
int MANGLE(petra_comm_getnumproc)(PETRA_COMM communicator)
{
  Petra_Comm* communicator_ = (Petra_Comm *) communicator;
  return(communicator_->getNumProc());
  
}

void MANGLE(petra_comm_barrier)(PETRA_COMM communicator)
{
  Petra_Comm* communicator_ = (Petra_Comm *) communicator;
  communicator_->barrier();
  
}

void MANGLE(petra_comm_destroy)(PETRA_COMM communicator)
{
    delete (Petra_Comm *) communicator;

}


#ifdef __cplusplus
}
#endif
