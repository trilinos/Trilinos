//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#ifndef EpetraExt_CRSMATRIX_TRANSPOSE_H
#define EpetraExt_CRSMATRIX_TRANSPOSE_H

#include <EpetraExt_Transform.h>

class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_Map;
class Epetra_Export;

namespace EpetraExt {

//! Transform to form the explicit transpose of a Epetra_RowMatrix
class RowMatrix_Transpose : public SameTypeTransform<Epetra_RowMatrix>
{

 public:

  //! Destructor
  ~RowMatrix_Transpose();

  //! Constructor
  /*!
    \param In
    TransposeRowMap - Map to be used for row mapping of transpose matrix
    \param In
    IgnoreNonLocalCols - Whether to ignore non-local columns for the transpose
   */
  RowMatrix_Transpose( Epetra_Map * TransposeRowMap = 0,
		       bool IgnoreNonLocalCols = false)
  : TransposeMatrix_(0),
    TransposeRowMap_(TransposeRowMap),
    IgnoreNonLocalCols_(IgnoreNonLocalCols),
    NumMyRows_(0),
    NumMyCols_(0),
    MaxNumEntries_(0),
    Indices_(0),
    Values_(0),
    OrigMatrixIsCrsMatrix_(false)
  {}

  //! Transpose Tranform Operator
  NewTypeRef operator()( OriginalTypeRef orig );

  //! Foward Data Migration
  bool fwd();

  //! Reverse Data Migration
  bool rvs();

  //! Release the pointer to TransposeMatrix_ (so you can take the matrix out w/o worring about deallocation)
  void ReleaseTranspose() {TransposeMatrix_=0;}

 private:
  Epetra_CrsMatrix * BuildTempTrans();

  Epetra_CrsMatrix * TransposeMatrix_;

  Epetra_Map * TransposeRowMap_;

  bool MakeDataContiguous_;
  bool IgnoreNonLocalCols_;

  int NumMyRows_;
  int NumMyCols_;
  int MaxNumEntries_;
  int * Indices_;
  double * Values_;

  bool OrigMatrixIsCrsMatrix_;

};

} //namespace EpetraExt

#endif //EpetraExt_ROWMATRIX_TRANSPOSE_H
