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


#ifndef _fei_MatrixTraits_Epetra_h_
#define _fei_MatrixTraits_Epetra_h_

#include <fei_trilinos_macros.hpp>

#ifdef HAVE_FEI_EPETRA

//
//IMPORTANT NOTE: Make sure that wherever this file is included from, it
//appears BEFORE any include of fei_Vector_Impl.hpp or fei_Matrix_Impl.hpp !!!
//
#include <fei_MatrixTraits.hpp>
#include <snl_fei_BlockMatrixTraits.hpp>
#include <fei_VectorTraits_Epetra.hpp>
#include <fei_Include_Trilinos.hpp>
#include <fei_Vector_Impl.hpp>

namespace fei {
  /** Declare an Epetra_CrsMatrix specialization of the
      fei::MatrixTraits struct.

      This allows Epetra_CrsMatrix to be used as the template parameter
      of the fei::Matrix class.
  */
  template<>
  struct MatrixTraits<Epetra_CrsMatrix> {
    static const char* typeName()
      { return("Epetra_CrsMatrix"); }

    static double* getBeginPointer(Epetra_CrsMatrix* mat)
      {
        return (*mat)[0];
      }

    static int getOffset(Epetra_CrsMatrix* A, int row, int col)
      {
        const Epetra_Map& erowmap = A->RowMap();
        const Epetra_Map& ecolmap = A->ColMap();
        int local_row = erowmap.LID(row);
        int local_col = ecolmap.LID(col);
    
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

    static int setValues(Epetra_CrsMatrix* mat, double scalar)
      {
        return( mat->PutScalar(scalar) );
      }

    static int getNumLocalRows(Epetra_CrsMatrix* mat, int& numRows)
    {
      numRows = mat->NumMyRows();
      return(0);
    }

    static int getRowLength(Epetra_CrsMatrix* mat, int row, int& length)
      {
        length = mat->NumGlobalEntries(row);
        if (length < 0) return(-1);
        return( 0 );
      }

    static int copyOutRow(Epetra_CrsMatrix* mat,
                      int row, int len, double* coefs, int* indices)
      {
        int dummy;
        return(mat->ExtractGlobalRowCopy(row, len, dummy, coefs, indices));
      }

    static int putValuesIn(Epetra_CrsMatrix* mat,
                     int numRows, const int* rows,
                     int numCols, const int* cols,
                     const double* const* values,
                           bool sum_into)
      {
        //!!! STATIC DATA NOT THREAD-SAFE !!!!!
        static std::vector<int> idx;
        idx.resize(numRows+numCols);
        int* idx_row = &idx[0];
        int* idx_col = idx_row+numRows;
        for(int i=0; i<numRows; ++i) {
          idx_row[i] = mat->RowMap().LID(rows[i]);
        }
        for(int i=0; i<numCols; ++i) {
          idx_col[i] = mat->ColMap().LID(cols[i]);
        }
        if (sum_into) {
          for(int i=0; i<numRows; ++i) {
            int err = mat->SumIntoMyValues(idx_row[i], numCols,
                                               (double*)values[i],
                                               idx_col);
            if (err != 0) {
              return(err);
            }
          }
        }
        else {
          for(int i=0; i<numRows; ++i) {
            int err = mat->ReplaceMyValues(idx_row[i], numCols,
                                               (double*)values[i],
                                               idx_col);
            if (err != 0) {
              return(err);
            }
          }
        }
        return(0);
      }

    static int globalAssemble(Epetra_CrsMatrix* mat)
    {
      if (!mat->Filled()) {
        int err = mat->FillComplete();
        if (err != 0) {
          fei::console_out() << "MatrixTraits<Epetra_CrsMatrix>::globalAssemble"
                   << " ERROR in mat->FillComplete" << FEI_ENDL;
          return(-1);
        }
      }

      if (!mat->StorageOptimized()) {
        mat->OptimizeStorage();
      }

      return( 0 );
    }

    static int matvec(Epetra_CrsMatrix* mat,
                      fei::Vector* x,
                      fei::Vector* y)
    {
      fei::Vector_Impl<Epetra_MultiVector>* evx =
        dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>* >(x);
      fei::Vector_Impl<Epetra_MultiVector>* evy =
        dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>* >(y);

      if (evx == NULL || evy == NULL) {
        return(-1);
      }

      Epetra_MultiVector* ex = evx->getUnderlyingVector();
      Epetra_MultiVector* ey = evy->getUnderlyingVector();

      return( mat->Multiply(false, *ex, *ey) );
    }

  };//struct MatrixTraits<Epetra_CrsMatrix>
}//namespace fei

namespace snl_fei {
  /** Declare an Epetra_VbrMatrix specialization of the
      snl_fei::BlockMatrixTraits struct.

      This allows Epetra_VbrMatrix to be used as the template parameter
      for the fei::Matrix class.
  */
  template<>
  struct BlockMatrixTraits<Epetra_VbrMatrix> {
    static const char* typeName()
      { return("Epetra_VbrMatrix"); }

    static int putScalar(Epetra_VbrMatrix* mat, double scalar)
      {
        return( mat->PutScalar(scalar) );
      }

    static int getRowLength(Epetra_VbrMatrix* mat, int row, int& length)
      {
        length = mat->NumGlobalBlockEntries(row);
        return(0);
      }

    static int getPointRowLength(Epetra_VbrMatrix* mat, int row, int& length)
    {
      const Epetra_Map& map = mat->RowMatrixRowMap();
      int minLocalRow = map.MinMyGID();
      int localRow = row - minLocalRow;
      int error = mat->NumMyRowEntries(localRow, length);
      return(error);
    }

    static int copyOutRow(Epetra_VbrMatrix* mat,
                          int row, int numBlkCols,
                          int rowDim,
                          int* blkCols,
                          int* colDims,
                          double* coefs,
                          int coefsLen,
                          int& blkRowLength)
      {
        int checkRowDim;
        int error = mat->BeginExtractGlobalBlockRowCopy(row, numBlkCols,
                                                        checkRowDim,
                                                        blkRowLength,
                                                        blkCols, colDims);
        if (error != 0 || checkRowDim != rowDim || blkRowLength != numBlkCols) {
          return(error);
        }

        int offset = 0;
        for(int i=0; i<numBlkCols; ++i) {
          if (offset >= coefsLen) {
            cerr << "BlockMatrixTraits::copyOutRow ran off end of coefs array."
                 << endl;
            return(-2);
          }
          int numValues = rowDim*colDims[i];
          error = mat->ExtractEntryCopy(numValues, &(coefs[offset]),
                                        rowDim, false);
          if (error != 0) {
            return(error);
          }
          offset += numValues;
        }

        return(0);
      }

    static int copyOutPointRow(Epetra_VbrMatrix* mat,
                               int firstLocalOffset,
                               int row,
                               int len,
                               double* coefs,
                               int* indices,
                               int& rowLength)
      {
        int error = mat->ExtractMyRowCopy(row-firstLocalOffset,
                                          len, rowLength,
                                          coefs, indices);

        const Epetra_Map& colmap = mat->RowMatrixColMap();
        for(int i=0; i<len; ++i) {
          indices[i] = colmap.GID(indices[i]);
        }

        return(error);
      }

    static int sumIn(Epetra_VbrMatrix* mat,
                     int blockRow,
                     int rowDim,
                     int numBlockCols,
                     const int* blockCols,
                     const int* colDims,
                     int LDA,
                     const double* values)
    {
      int err, voffset = 0;
      for(int j=0; j<numBlockCols; ++j) {
        err = mat->DirectSubmitBlockEntry(blockRow, blockCols[j],
                                          &(values[voffset]), LDA,
                                         rowDim, colDims[j], true/*sum_into*/);
        if (err != 0) return(err);

        voffset += colDims[j]*LDA;
      }

      return(0);
    }

    static int copyIn(Epetra_VbrMatrix* mat,
                      int blockRow,
                      int rowDim,
                      int numBlockCols,
                      const int* blockCols,
                      const int* colDims,
                      int LDA,
                      const double* values)
    {
      int err, voffset = 0;
      for(int j=0; j<numBlockCols; ++j) {
        err = mat->DirectSubmitBlockEntry(blockRow, blockCols[j],
                                         &(values[voffset]), LDA,
                                    rowDim, colDims[j], false/*replace*/);
        if (err != 0) return(err);

        voffset += colDims[j]*LDA;
      }

      return(0);
    }

    static int sumIn(Epetra_VbrMatrix* mat,
                     int row,
                     int rowDim,
                     int numCols,
                     const int* cols,
                     const int* LDAs,
                     const int* colDims,
                     const double* const* values)
      {
        int err = 0;
        for(int i=0; i<numCols; ++i) {
          err = mat->DirectSubmitBlockEntry(row, cols[i], values[i],
                                            LDAs[i], rowDim, colDims[i], true);
          if (err != 0) return(err);
        }

        return(err);
      }

    static int copyIn(Epetra_VbrMatrix* mat,
                      int row,
                      int rowDim,
                      int numCols,
                      const int* cols,
                      const int* LDAs,
                      const int* colDims,
                      const double* const* values)
      {
        int err = 0;
        for(int i=0; i<numCols; ++i) {
          err = mat->DirectSubmitBlockEntry(row, cols[i], values[i],
                                          LDAs[i], rowDim, colDims[i], false);
          if (err != 0) return(err);
        }

        return(err);
      }

    static int globalAssemble(Epetra_VbrMatrix* mat)
    {
      return( mat->FillComplete() );
    }
  };//struct BlockMatrixTraits<Epetra_VbrMatrix>
}//namespace snl_fei
#endif //HAVE_FEI_EPETRA
#endif // _fei_MatrixTraits_Epetra_hpp_
