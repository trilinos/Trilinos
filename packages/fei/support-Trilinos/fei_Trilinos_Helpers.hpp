/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Trilinos_Helpers_hpp_
#define _fei_Trilinos_Helpers_hpp_

#include "fei_trilinos_macros.hpp"
#include "fei_fwd.hpp"

#include <fei_Include_Trilinos.hpp>

#include <fei_mpi.h>
#include <fei_SharedPtr.hpp>

#include <fei_LinearProblemManager.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_MatrixGraph.hpp>

namespace Trilinos_Helpers {

#ifdef HAVE_FEI_EPETRA

  /** Epetra_Map objects are light-weight wrt copying, since they employ a
      memory-model that uses a reference-counted pointer to an internal data
      object. Thus, it is reasonable for this function to return a pass-by-value
      Epetra_Map object.
  */
  Epetra_Map create_Epetra_Map(MPI_Comm comm,
                               const std::vector<int>& local_eqns);

  /** Epetra_BlockMap objects are light-weight wrt copying, since they employ a
      memory-model that uses a reference-counted pointer to an internal data
      object. Thus, it is reasonable for this function to return a pass-by-value
      Epetra_BlockMap object.
  */
  Epetra_BlockMap
    create_Epetra_BlockMap(const fei::SharedPtr<fei::VectorSpace>& vecspace);

  Epetra_CrsGraph
    create_Epetra_CrsGraph(const fei::SharedPtr<fei::MatrixGraph>& matgraph,
                           bool blockEntries,
                           bool orderRowsWithLocalColsFirst=false);

  fei::SharedPtr<fei::Matrix>
    create_from_Epetra_Matrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              bool blockEntryMatrix,
                              fei::SharedPtr<fei::Reducer> reducer,
                              bool orderRowsWithLocalColsFirst=false);

  fei::SharedPtr<fei::Matrix>
    create_from_LPM_EpetraBasic(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                 bool blockEntryMatrix,
                                 fei::SharedPtr<fei::Reducer> reducer,
                                 fei::SharedPtr<fei::LinearProblemManager>
                                   lpm_epetrabasic);
#endif

  /** Copies parameters from fei::ParameterSet to Teuchos::ParameterList.
    Does not clear any pre-existing contents from the Teuchos:ParameterList.
  */
  void copy_parameterset(const fei::ParameterSet& paramset,
                         Teuchos::ParameterList& paramlist);

  /** Copies parameters from Teuchos::ParameterList to fei::ParameterSet.
    Does not clear any pre-existing contents from the fei:ParameterSet.
  */
  void copy_parameterlist(const Teuchos::ParameterList& paramlist,
                          fei::ParameterSet& paramset);

#ifdef HAVE_FEI_EPETRA
  /** Extracts a pointer to a Epetra_MultiVector from a fei::Vector. Throws
    an exception if unsuccessful.
  */
  Epetra_MultiVector*
    get_Epetra_MultiVector(fei::Vector* feivec, bool soln_vec);

  /** Extracts a pointer to a Epetra_VbrMatrix from a fei::Matrix. Throws
    an exception if unsuccessful.
  */
  Epetra_VbrMatrix* get_Epetra_VbrMatrix(fei::Matrix* feimat);

  /** Extracts a pointer to a Epetra_CrsMatrix from a fei::Matrix. Throws
    an exception if unsuccessful.
  */
  Epetra_CrsMatrix* get_Epetra_CrsMatrix(fei::Matrix* feimat);

  /** Extracts pointers to epetra objects from fei container-objects.
    If epetra objects can't be obtained, output arguments are set
    to NULL.
  */
  void get_Epetra_pointers(fei::SharedPtr<fei::Matrix> feiA,
                           fei::SharedPtr<fei::Vector> feix,
                           fei::SharedPtr<fei::Vector> feib,
                           Epetra_CrsMatrix*& crsA,
                           Epetra_Operator*& opA,
                           Epetra_MultiVector*& x,
                           Epetra_MultiVector*& b);

  /** Fill the VbrMatrix with zeros. */
  int zero_Epetra_VbrMatrix(Epetra_VbrMatrix* mat);

  /** Return the pointer to the first local coefficient in the matrix.
    Only works for Epetra_CrsMatrix.
  */
  inline
  double* getBeginPointer(const fei::SharedPtr<fei::Matrix>& feiA)
  {
    Epetra_CrsMatrix* A = get_Epetra_CrsMatrix(feiA.get());
    return (*A)[0];
  }

  /** Get the offset of the specified matrix location from the
    beginning of the local matrix storage.
    global_row_index and global_col_index are Epetra global indices.
  */
  inline
  int getOffsetG(const fei::SharedPtr<fei::Matrix>& feiA,
                 int global_row_index, int global_col_index)
  {
    Epetra_CrsMatrix* A = get_Epetra_CrsMatrix(feiA.get());
    const Epetra_Map& erowmap = A->RowMap();
    const Epetra_Map& ecolmap = A->ColMap();
    int local_row = erowmap.LID(global_row_index);
    int local_col = ecolmap.LID(global_col_index);
  
    int* rowOffsets;
    int* colIndices;
    double* coefs;
    A->ExtractCrsDataPointers(rowOffsets, colIndices, coefs);
  
    int* row_ptr = &colIndices[rowOffsets[local_row]];
    int* end_row = &colIndices[rowOffsets[local_row+1]];
  
    int col_offset = 0;
    for(; row_ptr != end_row; ++row_ptr) {
      if (*row_ptr == local_col) break;
      ++col_offset;
    }
  
    return rowOffsets[local_row] + col_offset;
  }

  /** Get the offset of the specified matrix location from the
    beginning of the local matrix storage.
    local_row_index and local_col_index are Epetra local indices.
  */
  inline
  int getOffsetL(const fei::SharedPtr<fei::Matrix>& feiA,
                 int local_row_index, int local_col_index)
  {
    Epetra_CrsMatrix* A = get_Epetra_CrsMatrix(feiA.get());
  
    int* rowOffsets;
    int* colIndices;
    double* coefs;
    A->ExtractCrsDataPointers(rowOffsets, colIndices, coefs);
  
    int* row_ptr = &colIndices[rowOffsets[local_row_index]];
    int* end_row = &colIndices[rowOffsets[local_row_index+1]];
  
    int col_offset = 0;
    for(; row_ptr != end_row; ++row_ptr) {
      if (*row_ptr == local_col_index) break;
      ++col_offset;
    }
  
    return rowOffsets[local_row_index] + col_offset;
  }

#endif // HAVE_FEI_EPETRA

}//namespace Trilinos_Helpers

#endif // _Trilinos_Helpers_hpp_

