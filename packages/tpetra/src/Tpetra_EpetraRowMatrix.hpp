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
  int globalNumRows = tpetra_matrix_->getRowMap()->getGlobalNumElements();
  int globalNumCols = tpetra_matrix_->getColMap()->getGlobalNumElements();
  Teuchos::ArrayView<const int> row_elem_list = tpetra_matrix_->getRowMap()->getNodeElementList();
  Teuchos::ArrayView<const int> col_elem_list = tpetra_matrix_->getColMap()->getNodeElementList();
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
