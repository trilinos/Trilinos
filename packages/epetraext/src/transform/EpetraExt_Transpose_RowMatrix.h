// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef EpetraExt_CRSMATRIX_TRANSPOSE_H
#define EpetraExt_CRSMATRIX_TRANSPOSE_H

#include <EpetraExt_Transform.h>

class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_Map;
class Epetra_Export;

namespace EpetraExt {

class RowMatrix_Transpose : public SameTypeTransform<Epetra_RowMatrix> {

 public:

  ~RowMatrix_Transpose();

  RowMatrix_Transpose( bool IgnoreNonLocalCols = false,
                       bool MakeDataContiguous = false )
  : IgnoreNonLocalCols_(IgnoreNonLocalCols),
    MakeDataContiguous_(MakeDataContiguous),
    TransposeMatrix_(0),
    TransposeExporter_(0),
    TransposeRowMap_(0),
    NumMyRows_(0),
    NumMyCols_(0),
    MaxNumEntries_(0),
    Indices_(0),
    Values_(0),
    TransNumNz_(0),
    TransIndices_(0),
    TransValues_(0),
    TransMyGlobalEquations_(0),
    OrigMatrixIsCrsMatrix_(false)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();

  bool rvs();

 private:

  Epetra_CrsMatrix * TransposeMatrix_;

  Epetra_Export * TransposeExporter_;

  Epetra_Map * TransposeRowMap_;

  bool MakeDataContiguous_;
  bool IgnoreNonLocalCols_;

  int NumMyRows_;
  int NumMyCols_;
  int MaxNumEntries_;
  int * Indices_;
  double * Values_;

  int * TransNumNz_;
  int ** TransIndices_;
  double ** TransValues_;
  int * TransMyGlobalEquations_;

  bool OrigMatrixIsCrsMatrix_;

};

} //namespace EpetraExt

#endif //EpetraExt_ROWMATRIX_TRANSPOSE_H
