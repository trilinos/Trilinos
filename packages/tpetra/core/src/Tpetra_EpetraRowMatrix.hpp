/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#ifndef TPETRA_EPETRAROWMATRIX_HPP
#define TPETRA_EPETRAROWMATRIX_HPP

#include <Tpetra_ConfigDefs.hpp>

#if defined(HAVE_TPETRA_EPETRA)

#include <Epetra_Comm.h>
#include <Epetra_BasicRowMatrix.h>
#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra
{

//! A class for wrapping a Tpetra::RowMatrix object in the Epetra_RowMatrix interface.
template<class TpetraMatrixType>
class EpetraRowMatrix : public Epetra_BasicRowMatrix {
public:
  EpetraRowMatrix(const Teuchos::RCP<TpetraMatrixType> &mat, const Epetra_Comm &comm);
  virtual ~EpetraRowMatrix();

  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  //not implemented
  int ExtractMyEntryView(int CurEntry, double * & Value, int & RowIndex, int & ColIndex);

  //not implemented
  int ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const;

  int NumMyRowEntries(int MyRow, int & NumEntries) const;

private:
  Teuchos::RCP<TpetraMatrixType> tpetra_matrix_;
};//class EpetraRowMatrix

template<class TpetraMatrixType>
EpetraRowMatrix<TpetraMatrixType>::EpetraRowMatrix(
  const Teuchos::RCP<TpetraMatrixType> &mat, const Epetra_Comm &comm
  )
 : Epetra_BasicRowMatrix(comm),
   tpetra_matrix_(mat)
{
  typedef typename TpetraMatrixType::global_ordinal_type GO;
  GO globalNumRows = tpetra_matrix_->getRowMap()->getGlobalNumElements();
  GO globalNumCols = tpetra_matrix_->getColMap()->getGlobalNumElements();
  Teuchos::ArrayView<const GO> row_elem_list = tpetra_matrix_->getRowMap()->getNodeElementList();
  Teuchos::ArrayView<const GO> col_elem_list = tpetra_matrix_->getColMap()->getNodeElementList();
  Epetra_Map rowmap(globalNumRows, row_elem_list.size(), row_elem_list.getRawPtr(), 0, comm);
  Epetra_Map colmap(globalNumCols, col_elem_list.size(), col_elem_list.getRawPtr(), 0, comm);
  SetMaps(rowmap, colmap);
}

template<class TpetraMatrixType>
EpetraRowMatrix<TpetraMatrixType>::~EpetraRowMatrix()
{
}

template<class TpetraMatrixType>
int EpetraRowMatrix<TpetraMatrixType>::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  Teuchos::ArrayView<int> inds(Indices, Length);
  Teuchos::ArrayView<double> vals(Values, Length);
  size_t num_entries = NumEntries;
  tpetra_matrix_->getLocalRowCopy(MyRow, inds, vals, num_entries);
  NumEntries = num_entries;
  return 0;
}

template<class TpetraMatrixType>
int EpetraRowMatrix<TpetraMatrixType>::ExtractMyEntryView(int CurEntry, double * & Value, int & RowIndex, int & ColIndex)
{
  //not implemented
  return -1;
}

template<class TpetraMatrixType>
int EpetraRowMatrix<TpetraMatrixType>::ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const
{
  //not implemented
  return -1;
}

template<class TpetraMatrixType>
int EpetraRowMatrix<TpetraMatrixType>::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  NumEntries = tpetra_matrix_->getNumEntriesInLocalRow(MyRow);
  return 0;
}

}//namespace Tpetra

#endif // defined(HAVE_TPETRA_EPETRA)

//here is the include-guard #endif:

#endif
