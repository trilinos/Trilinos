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

#ifndef _AztecDMSR_Matrix_h_
#define _AztecDMSR_Matrix_h_

#ifdef HAVE_FEI_AZTECOO


//
// This class is a wrapper for the Aztec DMSR matrix data structure.
//
// Important usage notes:
//
// * The 'oneBased' argument to the constructor indicates whether the
//   matrix should use 1-based indices (row numbers and column indices) in
//   the input and output arguments to its interfaces (e.g., getRow),
//   with the exception of the update list -- keep reading.
//   'oneBased' should be 1 for 1-based indices, 0 for 0-based.
//   Here's the confusing part -- the update list should contain 0-based
//   indices, regardless of the value of 'oneBased'.  That's because the
//   update list gets used internally by Aztec functions that only work
//   in 0-based numbers.
//
// * The 'rowLengths' array, input argument to the configure function,
//   must contain the lengths of each row, *NOT* including the
//   coefficient on the diagonal.
//
#include <az_aztec.h>
#include <fei_SharedPtr.hpp>
#include <fei_Aztec_Map.hpp>
#include "fei_iostream.hpp"
#include "fei_fstream.hpp"
#include "fei_sstream.hpp"

namespace fei_trilinos {

class Aztec_LSVector;

class AztecDMSR_Matrix {
    
  public:
    // Constructor.
    AztecDMSR_Matrix(fei::SharedPtr<Aztec_Map> map);

    //Copy constructor
    AztecDMSR_Matrix(const AztecDMSR_Matrix& src);

    virtual ~AztecDMSR_Matrix ();

    // Mathematical functions.
    void matvec(const Aztec_LSVector& x, Aztec_LSVector& y) const;

    void put(double s);
    void getDiagonal(Aztec_LSVector& diagVector) const;

    fei::SharedPtr<Aztec_Map> getAztec_Map() const {return(amap_);};

    int rowLength(int row) const;
    
    // ... to read matrix.
    void getRow(int row, int& length, double *coefs, int *colInd) const;
    void getRow(int row, int& length, double *coefs) const;
    void getRow(int row, int& length, int *colInd) const;

    /** write access to the diagonal entry for the specified row. */
    int setDiagEntry(int row, double value);

    /** Read-only access to the diagonal entry for the specified row. */
    double getDiagEntry(int row) const;

    // ... to write matrix.
    int putRow(int row, int len, const double *coefs, 
                       const int *colInd);

    int sumIntoRow(int numRows, const int* rows,
                 int numCols, const int* cols,
                 const double* const* coefs);

    int sumIntoRow(int row, int len, const double *coefs, 
                           const int *colInd);

    int addScaledMatrix(double scalar, const AztecDMSR_Matrix& source);

    void scale(double scalar);

    /** Special direct-access pointer function.
     */
    int getOffDiagRowPointers(int row, int*& colIndices, double*& coefs,
			      int& offDiagRowLength);

    void allocate(int *rowLengths);

    //inform about structure, including column-indices, so that val and bindx
    //can be allocated *and* so that bindx can be populated.
    void allocate(int *rowLengths, const int* const* colIndices);

    //inform that data fill is complete, so AZ_transform can be called.
    void fillComplete();

    bool isFilled() const {return(isFilled_);};
    void setFilled(bool flag) {isFilled_ = flag;};
    bool isAllocated() const {return(isAllocated_);};
    void setAllocated(bool flag) {isAllocated_ = flag;};

    void copyStructure(AztecDMSR_Matrix& source);

    bool readFromFile(const char *filename);
    bool writeToFile(const char *fileName) const;
    bool rowMax() const {return true;};
    double rowMax(int row) const;
 
    int getNumNonZeros() {return(nnzeros_);};

    double* getBeginPointer() { return val; }

    int getOffset(int row, int col)
    {
      int localRow;
      if (!amap_->inUpdate(row,localRow)){
        std::ostringstream oss;
        oss << "row "<<row<<" not found";
        std::string str = oss.str();
        throw std::runtime_error(str.c_str());
      }

      if (row == col) return localRow;

      int* row_ptr = &bindx[bindx[localRow]];
      int* end_row = &bindx[bindx[localRow+1]];

      int col_offset = 0;
      for(; row_ptr != end_row; ++row_ptr) {
        if (amap_->getTransformedEqn(*row_ptr) == col) break;
        ++col_offset;
      }
      if (row_ptr == end_row){
        FEI_OSTRINGSTREAM osstr;
        osstr << "Col "<<col << " not found for row "<<row;
        throw std::runtime_error(osstr.str());
      }
      return bindx[localRow] + col_offset;
    }

    //Aztec-specific functions:

    AZ_MATRIX* getAZ_MATRIX_PTR() {return(Amat_);};

  private:
    void messageAbort(const char* mesg);
    int insert(int item, int offset, int* list, int& len, int allocLen);
    int insert(double item, int offset, double* list, int& len, int allocLen);
    void expand_array(int*& array, int& arraylen, int newlen);
    void expand_array(double*& array, int& arraylen, int newlen);

    bool isFilled_;
    bool isAllocated_;

    int localOffset_;
    int localSize_;

    fei::SharedPtr<Aztec_Map> amap_;

    AZ_MATRIX* Amat_;

    bool arraysAllocated_;
    double *val;
    int *bindx;
    int *rowLengths_;
    int nnzeros_; //val and bindx are of length nnzeros_+1

    int N_update_;

    int* tmp_array_;
    int tmp_array_len_;
    double* dtmp_array_;
    int dtmp_array_len_;

    bool azTransformed_;
};

}//namespace fei_trilinos

#endif //HAVE_FEI_AZTECOO

#endif
