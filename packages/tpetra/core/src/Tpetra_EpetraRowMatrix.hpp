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

#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_EPETRA)

#include <Epetra_Comm.h>
#include <Epetra_BasicRowMatrix.h>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TestForException.hpp>
#include <memory> // std::shared_ptr
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace Tpetra {
namespace Details {

// Epetra_MpiComm actually has reference-counting value semantics,
// just like std::shared_ptr.  We only return
// std::shared_ptr<Epetra_Comm> because Epetra_Comm is an abstract
// base class, so we must return it by pointer.
std::shared_ptr<Epetra_Comm>
makeEpetraCommFromTeuchosComm (const Teuchos::Comm<int>& teuchosComm);

} // namespace Details
} // namespace Tpetra

namespace { // (anonymous)

template<class EpetraGlobalOrdinalType, class TpetraMapType>
Epetra_Map
tpetraToEpetraMapTmpl (const TpetraMapType& tpetraMap)
{
  using Teuchos::ArrayView;
  typedef typename TpetraMapType::global_ordinal_type TGO;
  typedef typename TpetraMapType::local_ordinal_type LO;
  typedef EpetraGlobalOrdinalType EGO;

  const TGO gblNumInds = static_cast<TGO> (tpetraMap.getGlobalNumElements ());
  const LO lclNumInds = static_cast<LO> (tpetraMap.getNodeNumElements ());
  ArrayView<const TGO> global_index_list = tpetraMap.getNodeElementList ();

  std::vector<EGO> global_index_list_epetra;
  const EGO* global_index_list_epetra_ptr = NULL;
  if (std::is_same<TGO, EGO>::value) {
    global_index_list_epetra_ptr =
      reinterpret_cast<const EGO*> (global_index_list.getRawPtr ());
  }
  else {
    global_index_list_epetra.resize (lclNumInds);
    for (LO k = 0; k < lclNumInds; ++k) {
      // TODO (mfh 11 Oct 2017) Detect overflow from TGO to EGO.
      global_index_list_epetra[k] = static_cast<EGO> (global_index_list[k]);
    }
    global_index_list_epetra_ptr = global_index_list_epetra.data ();
  }
  const EGO indexBase = tpetraMap.getIndexBase ();
  std::shared_ptr<Epetra_Comm> epetraComm =
    Tpetra::Details::makeEpetraCommFromTeuchosComm (* (tpetraMap.getComm ()));
  // It's OK for the shared_ptr to fall out of scope.  Subclasses of
  // Epetra_Comm have reference-counted value semantics, so passing
  // the Epetra_Comm by const reference into Epetra_Map's constructor
  // is safe.
  //
  // TODO (mfh 11 Oct 2017) Detect overflow from TGO to EGO, or from
  // LO to int.
  return Epetra_Map (static_cast<EGO> (gblNumInds),
                     static_cast<int> (lclNumInds),
                     global_index_list_epetra_ptr, indexBase, *epetraComm);
}

} // namespace (anonymous)

namespace Tpetra {

//! A class for wrapping a Tpetra::RowMatrix object in the Epetra_RowMatrix interface.
template<class TpetraMatrixType>
class EpetraRowMatrix : public Epetra_BasicRowMatrix {
public:
  EpetraRowMatrix(const Teuchos::RCP<TpetraMatrixType> &mat, const Epetra_Comm &comm);
  virtual ~EpetraRowMatrix() {};

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
  using Teuchos::RCP;
  typedef typename TpetraMatrixType::map_type tpetra_map_type;
#if ! defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  // Prefer using Epetra64 if it is enabled.
  typedef long long EGO;
#elif ! defined(EPETRA_NO_32BIT_GLOBAL_INDICES)
  // We don't have Epetra64, but we do have 32-bit indices.
  typedef int EGO;
#else
#  error "Epetra was not configured correctly.  Neither 64-bit nor 32-bit indices are enabled."
#endif
  const char tfecfFuncName[] = "EpetraRowMatrix: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (mat.is_null (), std::invalid_argument,
     "The input Tpetra matrix is null.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (mat->getRowMap ().is_null (), std::invalid_argument,
     "The input Tpetra matrix's row Map is null.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (mat->getColMap ().is_null (), std::invalid_argument,
     "The input Tpetra matrix's column Map is null.");

  RCP<const tpetra_map_type> tpetraRowMap = mat->getRowMap ();
  Epetra_Map epetraRowMap =
    tpetraToEpetraMapTmpl<EGO, tpetra_map_type> (*tpetraRowMap);
  RCP<const tpetra_map_type> tpetraColMap = mat->getColMap ();
  Epetra_Map epetraColMap =
    tpetraToEpetraMapTmpl<EGO, tpetra_map_type> (*tpetraColMap);
  this->SetMaps (epetraRowMap, epetraColMap);
}

template<class TpetraMatrixType>
int EpetraRowMatrix<TpetraMatrixType>::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  static_assert (std::is_same<typename TpetraMatrixType::scalar_type, double>::value,
                 "This code assumes that Tpetra::CrsMatrix's scalar_type is int.");
  static_assert (std::is_same<typename TpetraMatrixType::local_ordinal_type, int>::value,
                 "This code assumes that Tpetra::CrsMatrix's local_ordinal_type is int.");
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
