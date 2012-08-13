/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include "Epetra_ConfigDefs.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#endif

#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_C_wrappers.h"
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#endif

  /////////////////////////////////////////////////////
//                  Epetra_Comm                    //
///////////////////////////////////////////////////

#ifdef EPETRA_MPI
  EPETRA_OBJECT_PTR MANGLE(epetra_mpicomm_create1)() {
    Epetra_Comm *comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    return((EPETRA_OBJECT_PTR ) comm_);
  }
  EPETRA_OBJECT_PTR MANGLE(epetra_mpicomm_create2)(MPI_Comm * comm) {
    Epetra_Comm *comm_ = new Epetra_MpiComm(*comm);
    return((EPETRA_OBJECT_PTR ) comm_);
  }
#endif

  EPETRA_OBJECT_PTR MANGLE(epetra_serialcomm_create)() {
    Epetra_Comm *comm = new Epetra_SerialComm();
    return((EPETRA_OBJECT_PTR ) comm);
  }

  int MANGLE(epetra_comm_mypid)(EPETRA_OBJECT_REF comm) {
    Epetra_Comm *comm_ = (Epetra_Comm *) comm;
    return(comm_->MyPID());
  
  }
  int MANGLE(epetra_comm_numproc)(EPETRA_OBJECT_REF comm) {
    Epetra_Comm* comm_ = (Epetra_Comm *) comm;
    return(comm_->NumProc());
  
  }

  void MANGLE(epetra_comm_barrier)(EPETRA_OBJECT_REF comm) {
    Epetra_Comm* comm_ = (Epetra_Comm *) comm;
    comm_->Barrier();
  
  }

  void MANGLE(epetra_comm_destroy)(EPETRA_OBJECT_REF comm) {
    delete (Epetra_Comm *) comm;
  }

  /////////////////////////////////////////////////////
  //                  Epetra_Map                     //
  ///////////////////////////////////////////////////

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  EPETRA_OBJECT_PTR MANGLE(epetra_map_create1)(EPETRA_INT numGlobalElements,
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create2)(EPETRA_INT numGlobalElements,
                 EPETRA_INT numMyElements,
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(numMyElements), 
             EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create3)(EPETRA_INT numGlobalElements,
                 EPETRA_INT numLocalElements,
                 int *updateList, 
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(numLocalElements),
             updateList, EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  EPETRA_OBJECT_PTR MANGLE(epetra_map_create1_64)(EPETRA_LONG_LONG numGlobalElements,
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create2_64)(EPETRA_LONG_LONG numGlobalElements,
                 EPETRA_INT numMyElements,
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(numMyElements), 
             EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create3_64)(EPETRA_LONG_LONG numGlobalElements,
                 EPETRA_INT numLocalElements,
                 long long *updateList, 
                 EPETRA_INT indexBase,
                 EPETRA_OBJECT_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(numLocalElements),
             updateList, EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }
#endif

  int MANGLE(epetra_map_nummyelements)(EPETRA_OBJECT_REF map) {
    Epetra_Map * map_ = (Epetra_Map *) map;
    return(map_->NumMyElements());
  }

  long long MANGLE(epetra_map_numglobalelements)(EPETRA_OBJECT_REF map) {
    Epetra_Map * map_ = (Epetra_Map *) map;
    return(map_->NumGlobalElements64());
  }
#ifndef EPETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int * MANGLE(epetra_map_myglobalelements)(EPETRA_OBJECT_REF map) {
    Epetra_Map * map_ = (Epetra_Map *) map;
    return(map_->MyGlobalElements());
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  long long * MANGLE(epetra_map_myglobalelements_64)(EPETRA_OBJECT_REF map) {
    Epetra_Map * map_ = (Epetra_Map *) map;
    return(map_->MyGlobalElements64());
  }
#endif
#endif

  EPETRA_OBJECT_PTR MANGLE(epetra_map_comm)(EPETRA_OBJECT_REF map) {
    Epetra_Map * map_ = (Epetra_Map *) map;
    return((EPETRA_OBJECT_PTR) &(map_->Comm()));
  }

  void MANGLE(epetra_map_destroy)(EPETRA_OBJECT_REF map)
  {
    delete (Epetra_Map *) map;
  }

  /////////////////////////////////////////////////////
  //                  Epetra_Vector                  //
  ///////////////////////////////////////////////////

  EPETRA_OBJECT_PTR MANGLE(epetra_vector_create1)(EPETRA_OBJECT_REF map) {
    Epetra_Map& map_ = *(Epetra_Map *) map;
    Epetra_Vector *vector = new Epetra_Vector(map_);
    return((EPETRA_OBJECT_PTR ) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_vector_create2)(EPETRA_INT CopyValues, EPETRA_OBJECT_REF map,
              double * V) {
    Epetra_DataAccess CV = View;
    if (EPETRA_DEREF(CopyValues)==1) CV = Copy;
    Epetra_Map& map_ = *(Epetra_Map *) map;
    Epetra_Vector *vector = new Epetra_Vector(CV, map_, V);
    return((EPETRA_OBJECT_PTR ) vector);
  }

  int MANGLE(epetra_vector_putscalar)(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalar) {
    Epetra_Vector *x_ = (Epetra_Vector *) x;
    return(x_->PutScalar(EPETRA_DEREF(scalar)));
  }

  int MANGLE(epetra_vector_norm1)(EPETRA_OBJECT_REF x, double *scalar) {
    Epetra_Vector *x_ = (Epetra_Vector *) x;
    return(x_->Norm1(scalar));
  }

  int MANGLE(epetra_vector_norm2)(EPETRA_OBJECT_REF x, double *scalar) {
    Epetra_Vector *x_ = (Epetra_Vector *) x;
    return(x_->Norm2(scalar));
  }

  int MANGLE(epetra_vector_random)(EPETRA_OBJECT_REF x) {
    Epetra_Vector *x_ = (Epetra_Vector *) x;
    return(x_->Random());
  }

  int MANGLE(epetra_vector_update)(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalara, EPETRA_OBJECT_REF a, 
           EPETRA_DOUBLE scalarb, EPETRA_OBJECT_REF b, EPETRA_DOUBLE scalarx) {
    Epetra_Vector *x_ = (Epetra_Vector *) x;
    Epetra_Vector& a_ = *(Epetra_Vector *) a;
    Epetra_Vector& b_ = *(Epetra_Vector *) b;
    return(x_->Update(EPETRA_DEREF(scalara), a_, EPETRA_DEREF(scalarb), b_, EPETRA_DEREF(scalarx)));
  }

  void MANGLE(epetra_vector_print)(EPETRA_OBJECT_REF x) {
    cout << *(Epetra_Vector *) x;
  }

  void MANGLE(epetra_vector_destroy)(EPETRA_OBJECT_REF x) {
    delete (Epetra_Vector *) x;
  }

  /////////////////////////////////////////////////////
#ifdef SKIP4NOW /* Comment this out for now */
  /////////////////////////////////////////////////////
  //                  Epetra_DVBR_Matrix         //
  ///////////////////////////////////////////////////


  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_dvbr_matrix_create)
    (EPETRA_MAP rowmap)
  {
    Epetra_BlockMap& rowmap_ = *(Epetra_BlockMap *) rowmap;
    Epetra_DVBR_Matrix *B = new Epetra_DVBR_Matrix(rowmap_);
    return((EPETRA_OBJECT_PTR) B);
  }

  int MANGLE(epetra_rdp_dvbr_matrix_allocate)
    (EPETRA_MATRIX A, int* numNzBlks, int* blkColInds)
  {
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    return(B->allocate(numNzBlks, blkColInds));
  }
  int MANGLE(epetra_rdp_dvbr_matrix_putblockrow)
    (EPETRA_MATRIX A, EPETRA_INT blk_row, EPETRA_INT num_nz_blocks, 
     double* vals, int* blk_col_inds)
  {
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    return(B->putBlockRow( EPETRA_DEREF(blk_row), EPETRA_DEREF(num_nz_blocks), vals, 
         blk_col_inds));
  }

  int MANGLE(epetra_rdp_dvbr_matrix_fillcomplete)(EPETRA_MATRIX A)
  {
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    return(B->fillComplete());
  }

  int MANGLE(epetra_rdp_dvbr_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x,
              EPETRA_VECTOR y)
  {
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    const Epetra_Vector& x1 = *(Epetra_Vector *) x;
    Epetra_Vector& y1 = *(Epetra_Vector *) y;
    return(B->matvec(x1, y1));
  }

  int MANGLE(epetra_rdp_dvbr_matrix_matmultivec)(EPETRA_MATRIX A,
             EPETRA_MULTIVECTOR x,
             EPETRA_MULTIVECTOR y)
  {
    Epetra_DVBR_Matrix *B = (Epetra_DVBR_Matrix *) A;
    const Epetra_MultiVector& x1 = *(Epetra_MultiVector *) x;
    Epetra_MultiVector& y1 = *(Epetra_MultiVector *) y;
    return(B->matvec(x1, y1));
  }

  void MANGLE(epetra_rdp_dvbr_matrix_destroy)(EPETRA_MATRIX A)
  {
    delete (Epetra_DVBR_Matrix *) A;
  }

  /////////////////////////////////////////////////////
  //                  Epetra_DCRS_Matrix         //
  ///////////////////////////////////////////////////


  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_dcrs_matrix_create) (EPETRA_MAP rowmap)
  {
    Epetra_Map& rowmap_ = *(Epetra_Map *) rowmap;
    Epetra_DCRS_Matrix *B = new Epetra_DCRS_Matrix(rowmap_);
    return((EPETRA_OBJECT_PTR) B);
  }

  int MANGLE(epetra_rdp_dcrs_matrix_allocate)
    (EPETRA_MATRIX A, int* rowLengths)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    return(B->allocate(rowLengths));
  }
  int MANGLE(epetra_rdp_dcrs_matrix_putrow)(EPETRA_MATRIX A, EPETRA_INT row,
              EPETRA_INT num_nz, 
              double* vals, int* col_inds)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    return(B->putRow( EPETRA_DEREF(row), EPETRA_DEREF(num_nz), vals, col_inds));
  }

  int MANGLE(epetra_rdp_dcrs_matrix_sumintodiagonal)
    (EPETRA_MATRIX A, double* diagonal)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    return(B->sumIntoDiagonal( diagonal));
  }

  int MANGLE(epetra_rdp_dcrs_matrix_fillcomplete)(EPETRA_MATRIX A)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    return(B->fillComplete());
  }

  int MANGLE(epetra_rdp_dcrs_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x,
              EPETRA_VECTOR y)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    const Epetra_Vector& x1 = *(Epetra_Vector *) x;
    Epetra_Vector& y1 = *(Epetra_Vector *) y;
    return(B->matvec(x1, y1));
  }

  int MANGLE(epetra_rdp_dcrs_matrix_matmultivec)(EPETRA_MATRIX A,
             EPETRA_MULTIVECTOR x,
             EPETRA_MULTIVECTOR y)
  {
    Epetra_DCRS_Matrix *B = (Epetra_DCRS_Matrix *) A;
    const Epetra_MultiVector& x1 = *(Epetra_MultiVector *) x;
    Epetra_MultiVector& y1 = *(Epetra_MultiVector *) y;
    return(B->matvec(x1, y1));
  }


  void MANGLE(epetra_rdp_dcrs_matrix_destroy)(EPETRA_MATRIX A)
  {
    delete (Epetra_DCRS_Matrix *) A;
  }

  //                  Epetra_MultiVector         //
  ///////////////////////////////////////////////////

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create)()
  {
    Epetra_MultiVector *vector = new Epetra_MultiVector();
    return((EPETRA_OBJECT_PTR) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create1)
    (EPETRA_MAP map, EPETRA_INT numVectors)
  {
    Epetra_Map& map_ = *(Epetra_Map *) map;
    Epetra_MultiVector *vector = new Epetra_MultiVector(map_, EPETRA_DEREF(numVectors));
    return((EPETRA_OBJECT_PTR) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create2)(EPETRA_MAP map, 
                    double *A, EPETRA_INT lda, EPETRA_INT numVectors)
  {
    Epetra_Map& map_ = *(Epetra_Map *) map;
    Epetra_MultiVector *vector = new Epetra_MultiVector(map_, A, EPETRA_DEREF(lda),
              EPETRA_DEREF(numVectors));
    return((EPETRA_OBJECT_PTR) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create3)(EPETRA_MAP map, 
                    double **in_multiVector, EPETRA_INT numVectors)
  {
    Epetra_Map& map_ = *(Epetra_Map *) map;
    Epetra_MultiVector *vector = new Epetra_MultiVector(map_, in_multiVector,
              EPETRA_DEREF(numVectors));
    return((EPETRA_OBJECT_PTR) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create4)
    (EPETRA_MULTIVECTOR in_multiVector)
  {
    Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
    Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_);
    return((EPETRA_OBJECT_PTR) vector);
  }

  EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create5)(EPETRA_MULTIVECTOR 
                    in_multiVector, EPETRA_INT numVectors, int *vecIndices)
  {
    Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
    Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_,
              EPETRA_DEREF(numVectors), vecIndices));
  return((EPETRA_OBJECT_PTR) vector);
}

EPETRA_OBJECT_PTR MANGLE(epetra_rdp_multivector_create6)(EPETRA_MULTIVECTOR
                  in_multiVector, EPETRA_INT startindex, EPETRA_INT numvectors)
{
  Epetra_MultiVector & in_multiVector_ = *(Epetra_MultiVector *) in_multiVector;
  Epetra_MultiVector *vector = new Epetra_MultiVector(in_multiVector_, EPETRA_DEREF(startindex),
                  EPETRA_DEREF(numvectors));
  return((EPETRA_OBJECT_PTR) vector);
}

int MANGLE(epetra_rdp_multivector_putmultivector)
  (EPETRA_MULTIVECTOR multiVector, 
   double **in_multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  const double ** t = (const double **) in_multiVector;
  return(multiVector_->putMultiVector(t));
}

int MANGLE(epetra_rdp_multivector_allocate)(EPETRA_MULTIVECTOR multiVector, 
              EPETRA_MAP map, EPETRA_INT numVectors)
{
  Epetra_Map& map_ = *(Epetra_Map *) map;
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->allocate(map_, EPETRA_DEREF(numVectors)));
}

int MANGLE(epetra_rdp_multivector_putscalar)
  (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->putScalar(EPETRA_DEREF(scalar)));
}

int MANGLE(epetra_rdp_multivector_scale)
  (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->scale(EPETRA_DEREF(scalar)));
}

int MANGLE(epetra_rdp_multivector_scalecopy)
  (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
   EPETRA_DOUBLE scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->scaleCopy(multiVector_in_, EPETRA_DEREF(scalar)));
}

int MANGLE(epetra_rdp_multivector_dotprod)
  (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
   double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->dotProd(multiVector_in_, scalar));
}

int MANGLE(epetra_rdp_multivector_addvec)
  (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar, 
   EPETRA_MULTIVECTOR multiVector_in)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& multiVector_in_ = *(Epetra_MultiVector *) multiVector_in;
  return(multiVector_->addVec(EPETRA_DEREF(scalar), multiVector_in_)));
}

int MANGLE(epetra_rdp_multivector_norm1)
  (EPETRA_MULTIVECTOR multiVector, double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->norm1(scalar));
}

int MANGLE(epetra_rdp_multivector_norm2)
  (EPETRA_MULTIVECTOR multiVector, double *scalar)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->norm2(scalar));
}

int MANGLE(epetra_rdp_multivector_lincomb)(EPETRA_MULTIVECTOR multiVector,
             EPETRA_MULTIVECTOR b, 
             EPETRA_DOUBLE scalar, EPETRA_MULTIVECTOR c)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& b_ = *(Epetra_MultiVector *) b;
  Epetra_MultiVector& c_ = *(Epetra_MultiVector *) c;
  return(multiVector_->linComb(b_,EPETRA_DEREF(scalar,c_)));
}

int MANGLE(epetra_rdp_multivector_random)
  (EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->random());
}

int MANGLE(epetra_rdp_multivector_reduce)(EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->reduce());
}

int MANGLE(epetra_rdp_multivector_numvectors)(EPETRA_MULTIVECTOR multiVector)
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  return(multiVector_->numVectors());
}

int MANGLE(epetra_rdp_multivector_gemm)(EPETRA_MULTIVECTOR multiVector,
          EPETRA_INT transa, EPETRA_INT transb, EPETRA_DOUBLE alpha,
          EPETRA_MULTIVECTOR A, EPETRA_MULTIVECTOR B,
          EPETRA_DOUBLE beta )
{
  Epetra_MultiVector *multiVector_ = (Epetra_MultiVector *) multiVector;
  Epetra_MultiVector& A_ = *(Epetra_MultiVector *) A;
  Epetra_MultiVector& B_ = *(Epetra_MultiVector *) B;
  bool transa_ = !(EPETRA_DEREF(transa==0));
  bool transb_ = !(EPETRA_DEREF(transb==0));
  return(multiVector_->GEMM(transa_, transb_, EPETRA_DEREF(alpha), A_, B_, EPETRA_DEREF(beta)));
}

void MANGLE(epetra_rdp_multivector_destroy)(EPETRA_MULTIVECTOR multiVector)
{
  delete (Epetra_MultiVector *) multiVector;
}

/////////////////////////////////////////////////////
//                  Epetra_BlockMap                //
///////////////////////////////////////////////////

EPETRA_OBJECT_PTR MANGLE(epetra_blockmap_create1)(
              EPETRA_INT numGlobalElements, EPETRA_INT numLocalElements, int *updateList,
              EPETRA_INT numGlobalBlocks, EPETRA_INT numLocalBlocks, 
              int *blockUpdateList,
              int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_BlockMap *blockmap = new Epetra_BlockMap(EPETRA_DEREF(numGlobalElements),
              EPETRA_DEREF(numLocalElements), updateList,
              EPETRA_DEREF(numGlobalBlocks), EPETRA_DEREF(numLocalBlocks),
              blockUpdateList,
              blockSizes, EPETRA_DEREF(indexBase), comm_);
  return((EPETRA_OBJECT_PTR ) blockmap);
}

EPETRA_OBJECT_PTR MANGLE(epetra_blockmap_create2)(
              EPETRA_INT numGlobalBlocks, EPETRA_INT numLocalBlocks, 
              int *blockUpdateList,
              int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_BlockMap *blockmap = new Epetra_BlockMap(
              EPETRA_DEREF(numGlobalBlocks), EPETRA_DEREF(numLocalBlocks), 
              blockUpdateList,
              blockSizes, EPETRA_DEREF(indexBase), comm_);
  return((EPETRA_OBJECT_PTR ) blockmap);
}

void MANGLE(epetra_blockmap_destroy)(EPETRA_BLOCKMAP blockmap)
{
  delete (Epetra_BlockMap *) blockmap;
}

/////////////////////////////////////////////////////
//                  Epetra_LocalMap                //
///////////////////////////////////////////////////

EPETRA_OBJECT_PTR MANGLE(epetra_localmap_create)(EPETRA_INT numLocalElements,
                   EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalMap *localmap = new Epetra_LocalMap(EPETRA_DEREF(numLocalElements),
              EPETRA_DEREF(indexBase), comm_);
  return((EPETRA_OBJECT_PTR ) localmap);
}

void MANGLE(epetra_localmap_destroy)(EPETRA_LOCALMAP localmap)
{
  delete (Epetra_LocalMap *) localmap;
}

/////////////////////////////////////////////////////
//                  Epetra_LocalBlockMap           //
///////////////////////////////////////////////////

EPETRA_OBJECT_PTR MANGLE(epetra_localblockmap_create1)(
                  EPETRA_INT numLocalElements,
                  EPETRA_INT numLocalBlocks,
                  int* blockSizes,
                  EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalBlockMap *localblockmap = new
    Epetra_LocalBlockMap(EPETRA_DEREF(numLocalElements),
       EPETRA_DEREF(numLocalBlocks),
       blockSizes,
       EPETRA_DEREF(indexBase), comm_);
  return((EPETRA_OBJECT_PTR ) localblockmap);
}

EPETRA_OBJECT_PTR MANGLE(epetra_localblockmap_create2)(
                  EPETRA_INT numLocalBlocks,
                  int* blockSizes,
                  EPETRA_INT indexBase, EPETRA_COMM comm)
{
  Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
  Epetra_LocalBlockMap *localblockmap = new
    Epetra_LocalBlockMap(EPETRA_DEREF(numLocalBlocks),
       blockSizes,
       EPETRA_DEREF(indexBase), comm_);
  return((EPETRA_OBJECT_PTR ) localblockmap);
}

void MANGLE(epetra_localblockmap_destroy)(EPETRA_LOCALBLOCKMAP localblockmap)
{
  delete (Epetra_LocalBlockMap *) localblockmap;
}


#endif /* 0 */
