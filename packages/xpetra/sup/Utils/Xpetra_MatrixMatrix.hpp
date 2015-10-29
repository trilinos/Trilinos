// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_XPETRA_TPETRA

namespace Xpetra {

/*!
  @class Helpers
  @brief Xpetra utility class containing transformation routines  between Xpetra::Matrix and Epetra/Tpetra objects

  Note: this class is not in the Xpetra_UseShortNames.hpp
 */
template <class Scalar,
class LocalOrdinal  = int,
class GlobalOrdinal = LocalOrdinal,
class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
class Helpers {

private:
#include "Xpetra_UseShortNames.hpp"

public:

#ifdef HAVE_XPETRA_EPETRA
  static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    RCP<const Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Op);
    if (crsOp == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
    return A;
  } //Op2EpetraCrs

  static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    RCP<Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Op);
    if (crsOp == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstEpetraCrs

  static const Epetra_CrsMatrix & Op2EpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op) {
    RCP<const Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    try {
      const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & >(Op);
      RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp.getCrsMatrix();
      const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      return *A;
    } catch(...) {
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    }
  } //Op2EpetraCrs

  static Epetra_CrsMatrix & Op2NonConstEpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op) {
    RCP<Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    try {
      const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & >(Op);
      RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp.getCrsMatrix();
      const RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal> > &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal> >(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
      A = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(tmp_ECrsMtx->getEpetra_CrsMatrix());
      return *A;
    } catch(...) {
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    }
  } //Op2NonConstEpetraCrs
#endif // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
  static RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    // Get the underlying Tpetra Mtx
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Op);
    if (crsOp == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrix();
    return A;
  } //Op2TpetraCrs

  static RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    // Get the underlying Tpetra Mtx
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Op);
    if (crsOp == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstTpetraCrs

  static const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op2TpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op) {
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    // Get the underlying Tpetra Mtx
    try{
      const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & crsOp =dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & >(Op);
      RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp.getCrsMatrix();
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tmp_CrsMtx);
      if (tmp_TCrsMtx == Teuchos::null)
        throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
      A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      return *A;
    } catch (...) {
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    }
  } //Op2TpetraCrs

  static Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op2NonConstTpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Op) {
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    // Get the underlying Tpetra Mtx
    try{
      const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> & >(Op);
      RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_CrsMtx = crsOp.getCrsMatrix();
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tmp_CrsMtx);
      if (tmp_TCrsMtx == Teuchos::null)
        throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
      A = Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tmp_TCrsMtx->getTpetra_CrsMatrix());
      return *A;
    } catch (...) {
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    }
  } //Op2NonConstTpetraCrs

#endif // HAVE_XPETRA_TPETRA

};

template <class Scalar,
class LocalOrdinal  = int,
class GlobalOrdinal = LocalOrdinal,
class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
class MatrixMatrix {

private:
#undef XPETRA_MATRIXMATRIX_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
  In a parallel setting, A and B need not have matching distributions,
  but C needs to have the same row-map as A (if transposeA is false).
  At this time C=AT*B and C=A*BT are known to not work. However,
  C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
       then C's graph must already contain all nonzero locations that
       will be produced when forming the product A*B. On exit,
       C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

   */
  static void Multiply(
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
      bool transposeA,
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
      bool transposeB,
      Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
      bool call_FillComplete_on_result = true,
      bool doOptimizeStorage = true,
      const std::string & label = std::string()) {

    if(transposeA == false && C.getRowMap()->isSameAs(*A.getRowMap()) == false) {
      std::string msg = "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as row map of A";
      throw(Xpetra::Exceptions::RuntimeError(msg));
    }
    else if(transposeA == true && C.getRowMap()->isSameAs(*A.getDomainMap()) == false) {
      std::string msg = "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as domain map of A";
      throw(Xpetra::Exceptions::RuntimeError(msg));
    }


    if (!A.isFillComplete())
      throw(Xpetra::Exceptions::RuntimeError("A is not fill-completed"));
    if (!B.isFillComplete())
      throw(Xpetra::Exceptions::RuntimeError("B is not fill-completed"));

    bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

    if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#       ifndef HAVE_XPETRA_EPETRAEXT
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply requires EpetraExt to be compiled."));
#else
      Epetra_CrsMatrix & epA = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstEpetraCrs(A);
      Epetra_CrsMatrix & epB = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstEpetraCrs(B);
      Epetra_CrsMatrix & epC = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstEpetraCrs(C);

      int i = EpetraExt::MatrixMatrix::Multiply(epA,transposeA,epB,transposeB,epC,haveMultiplyDoFillComplete);
      if (haveMultiplyDoFillComplete) {
        // Due to Epetra wrapper intricacies, we need to explicitly call
        // fillComplete on Xpetra matrix here. Specifically, EpetraCrsMatrix
        // only keeps an internal variable to check whether we are in resumed
        // state or not, but never touches the underlying Epetra object. As
        // such, we need to explicitly update the state of Xpetra matrix to
        // that of Epetra one afterwords
        C.fillComplete();
      }

      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }

#endif
    } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & tpA = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & tpB = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &       tpC = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstTpetraCrs(C);

      //18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
      //Previously, Tpetra's matrix matrix multiply did not support fillComplete.
      Tpetra::MatrixMatrix::Multiply(tpA,transposeA,tpB,transposeB,tpC,haveMultiplyDoFillComplete,label);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
    }

    if(call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
      params->set("Optimize Storage",doOptimizeStorage);
      C.fillComplete((transposeB) ? B.getRangeMap() : B.getDomainMap(),
          (transposeA) ? A.getDomainMap() : A.getRangeMap(),
              params);
    }

    // transfer striding information
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(A));
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB); // TODO use references instead of RCPs
  } // end Multiply

  /**
   @brief Helper function to do matrix-matrix multiply

   Given CrsMatrix objects A, B and C, form the product C = A*B.
   In a parallel setting, A and B need not have matching distributions,
   but C needs to have the same row-map as A (if transposeA is false).
   At this time C=AT*B and C=A*BT are known to not work. However,
   C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. If Teuchos::null, a new CrsMatrix is created with optimal number of nnz per row.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

   */
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Multiply(const Matrix& A, bool transposeA,
      const Matrix& B, bool transposeB,
      RCP<Matrix> C_in,
      Teuchos::FancyOStream& fos,
      bool doFillComplete           = true,
      bool doOptimizeStorage        = true,
      const std::string & label     = std::string()) {

    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    // Default case: Xpetra Multiply
    RCP<Matrix> C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow = Teuchos::as<double>(0);

      if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
        // For now, follow what ML and Epetra do.
        GO numRowsA = A.getGlobalNumRows();
        GO numRowsB = B.getGlobalNumRows();
        nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
            sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
        nnzPerRow *=  nnzPerRow;
        double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
        double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
        if (totalNnz < minNnz)
          totalNnz = minNnz;
        nnzPerRow = totalNnz / A.getGlobalNumRows();

        fos << "Matrix product nnz per row estimate = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
      }

      if (transposeA) C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else            C = MatrixFactory::Build(A.getRowMap(),    Teuchos::as<LO>(nnzPerRow));

    } else {
      C->resumeFill(); // why this is not done inside of Tpetra MxM?
      fos << "Reuse C pattern" << std::endl;
    }

    Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage,label); // call Multiply routine from above

    return C;
  }

  /*! @brief Helper function to do matrix-matrix multiply

  Returns C = AB.

  @param A left matrix
  @param transposeA if true, use the transpose of A
  @param B right matrix
  @param transposeB if true, use the transpose of B
  @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
   */
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Multiply(const Matrix & A,
      bool transposeA,
      const Matrix & B,
      bool transposeB,
      Teuchos::FancyOStream &fos,
      bool callFillCompleteOnResult = true,
      bool doOptimizeStorage        = true,
      const std::string & label     = std::string()){
    return Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage,label);
  }

#ifdef HAVE_XPETRA_EPETRAEXT
  // Michael Gee's MLMultiply
  static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
      const Epetra_CrsMatrix& epB,
      Teuchos::FancyOStream& fos) {
    throw(Xpetra::Exceptions::RuntimeError("MLTwoMatrixMultiply only available for SC=double and GO=LO=int."));
    return Teuchos::null;
  }
#endif //ifdef HAVE_XPETRA_EPETRAEXT

  /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
   */
  static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(BlockedCrsMatrix& A, bool transposeA,
      BlockedCrsMatrix& B, bool transposeB,
      Teuchos::FancyOStream& fos,
      bool doFillComplete    = true,
      bool doOptimizeStorage = true) {
    if (transposeA || transposeB)
      throw Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

    // Preconditions
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgmapextractor = A.getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domapextractor = B.getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));

    for (size_t i = 0; i < A.Rows(); ++i) { // loop over all block rows of A
      for (size_t j = 0; j < B.Cols(); ++j) { // loop over all block columns of B
        RCP<Matrix> Cij;

        for (size_t l = 0; l < B.Rows(); ++l) { // loop for calculating entry C_{ij}
          RCP<CrsMatrix> crmat1 = A.getMatrix(i,l);
          RCP<CrsMatrix> crmat2 = B.getMatrix(l,j);

          if (crmat1.is_null() || crmat2.is_null()) {
            continue;
          }

          RCP<CrsMatrixWrap> crop1 = rcp(new CrsMatrixWrap(crmat1));
          RCP<CrsMatrixWrap> crop2 = rcp(new CrsMatrixWrap(crmat2));

          RCP<Matrix> temp = Multiply (*crop1, false, *crop2, false, fos);

          if (Cij.is_null ())
            Cij = temp;
          else
            Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd (*temp, false, 1.0, *Cij, 1.0);
        }

        if (!Cij.is_null())  {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));

          RCP<CrsMatrixWrap> crsCij = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Cij);
          TEUCHOS_TEST_FOR_EXCEPTION(Cij.is_null(), Xpetra::Exceptions::BadCast,
              "MatrixFactory failed in generating a CrsMatrixWrap." );

          RCP<CrsMatrix> crsMatCij = crsCij->getCrsMatrix();
          C->setMatrix(i, j, crsMatCij);

        } else {
          C->setMatrix(i, j, Teuchos::null);
        }
      }
    }

    if (doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  } // TwoMatrixMultiplyBlock

  /*! @brief Helper function to calculate B = alpha*A + beta*B.

  @param A      left matrix operand
  @param transposeA indicate whether to use transpose of A
  @param alpha  scalar multiplier for A
  @param B      right matrix operand
  @param beta   scalar multiplier for B

  @return sum in B.

  Note that B does not have to be fill-completed.
   */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
      throw Exceptions::RuntimeError("TwoMatrixAdd for Epetra matrices needs <double,int,int> for Scalar, LocalOrdinal and GlobalOrdinal.");
    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpA = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2TpetraCrs(A);
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpB = Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
    }
  } //MatrixMatrix::TwoMatrixAdd()


  /*! @brief Helper function to calculate C = alpha*A + beta*B.

  @param A          left matrix operand
  @param transposeA indicate whether to use transpose of A
  @param alpha      scalar multiplier for A, defaults to 1.0
  @param B          right matrix operand
  @param transposeB indicate whether to use transpose of B
  @param beta       scalar multiplier for B, defaults to 1.0
  @param C          resulting sum

  It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
   */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
      const Matrix& B, bool transposeB, const SC& beta,
      RCP<Matrix>& C,  Teuchos::FancyOStream &fos, bool AHasFixedNnzPerRow = false) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (C == Teuchos::null) {
      if (!A.isFillComplete() || !B.isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Global statistics are not available for estimates.");

      size_t maxNzInA     = A.getGlobalMaxNumRowEntries();
      size_t maxNzInB     = B.getGlobalMaxNumRowEntries();
      size_t numLocalRows = A.getNodeNumRows();

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        // first check if either A or B has at most 1 nonzero per row
        // the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

        if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

        } else {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }

        fos << "MatrixMatrix::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
            << ", using static profiling" << std::endl;
        C = rcp(new CrsMatrixWrap(A.getRowMap(), exactNnzPerRow, Xpetra::StaticProfile));

      } else {
        // general case
        double nnzPerRowInA = Teuchos::as<double>(A.getGlobalNumEntries()) / A.getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B.getGlobalNumEntries()) / B.getGlobalNumRows();
        LO    nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        fos << "MatrixMatrix::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
            << ", max possible nnz per row in sum = " << maxPossible
            << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
            << std::endl;
        C = rcp(new CrsMatrixWrap(A.getRowMap(), nnzToAllocate, pft));
      }
      if (transposeB)
        fos << "MatrixMatrix::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw Exceptions::RuntimeError("MatrixMatrix::Add for Epetra only available with Scalar = double, LO = GO = int.");
    } else if (C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpA =
          Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpB =
          Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >  tpC =
          Xpetra::Helpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, transposeB, beta, tpC);
#else
      throw Exceptions::RuntimeError("Xpetra must be compile with Tpetra.");
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
    if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
    ///////////////////////// EXPERIMENTAL

  } //MatrixMatrix::TwoMatrixAdd()


}; // class MatrixMatrix


// specialization MatrixMatrix for SC=double, LO=GO=int
template<>
class MatrixMatrix<double,int,int> {
  typedef double                                                     SC;
  typedef int                                                        LO;
  typedef int                                                        GO;
  typedef KokkosClassic::DefaultNode::DefaultNodeType                NO;
  typedef Xpetra::Map<int,int,NO>                                    Map;
  typedef Xpetra::Matrix<double,int,int,NO>                          Matrix;
  typedef Xpetra::MultiVector<double,int,int,NO>                     MultiVector;

public:

  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
  In a parallel setting, A and B need not have matching distributions,
  but C needs to have the same row-map as A (if transposeA is false).
  At this time C=AT*B and C=A*BT are known to not work. However,
  C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
       then C's graph must already contain all nonzero locations that
       will be produced when forming the product A*B. On exit,
       C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

   */
  static void Multiply(
      const Xpetra::Matrix<SC, LO, GO, NO>& A,
      bool transposeA,
      const Xpetra::Matrix<SC, LO, GO, NO>& B,
      bool transposeB,
      Xpetra::Matrix<SC, LO, GO, NO>& C,
      bool call_FillComplete_on_result = true,
      bool doOptimizeStorage = true,
      const std::string & label = std::string()) {
    if(transposeA == false && C.getRowMap()->isSameAs(*A.getRowMap()) == false) {
      std::string msg = "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as row map of A";
      throw(Xpetra::Exceptions::RuntimeError(msg));
    }
    else if(transposeA == true && C.getRowMap()->isSameAs(*A.getDomainMap()) == false) {
      std::string msg = "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as domain map of A";
      throw(Xpetra::Exceptions::RuntimeError(msg));
    }


    if (!A.isFillComplete())
      throw(Xpetra::Exceptions::RuntimeError("A is not fill-completed"));
    if (!B.isFillComplete())
      throw(Xpetra::Exceptions::RuntimeError("B is not fill-completed"));

    bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

    if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#       ifndef HAVE_XPETRA_EPETRAEXT
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply requires EpetraExt to be compiled."));
#else
      Epetra_CrsMatrix & epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(A);
      Epetra_CrsMatrix & epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(B);
      Epetra_CrsMatrix & epC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(C);

      int i = EpetraExt::MatrixMatrix::Multiply(epA,transposeA,epB,transposeB,epC,haveMultiplyDoFillComplete);
      if (haveMultiplyDoFillComplete) {
        // Due to Epetra wrapper intricacies, we need to explicitly call
        // fillComplete on Xpetra matrix here. Specifically, EpetraCrsMatrix
        // only keeps an internal variable to check whether we are in resumed
        // state or not, but never touches the underlying Epetra object. As
        // such, we need to explicitly update the state of Xpetra matrix to
        // that of Epetra one afterwords
        C.fillComplete();
      }

      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }

#endif
    } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT
      const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<SC,LO,GO,NO> &       tpC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(C);

      //18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
      //Previously, Tpetra's matrix matrix multiply did not support fillComplete.
      Tpetra::MatrixMatrix::Multiply(tpA,transposeA,tpB,transposeB,tpC,haveMultiplyDoFillComplete,label);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=int enabled."));
#endif
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
    }

    if(call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
      params->set("Optimize Storage",doOptimizeStorage);
      C.fillComplete((transposeB) ? B.getRangeMap() : B.getDomainMap(),
          (transposeA) ? A.getDomainMap() : A.getRangeMap(),
              params);
    }

    // transfer striding information
    RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(A));
    RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB); // TODO use references instead of RCPs
  } // end Multiply

  /**
   @brief Helper function to do matrix-matrix multiply

   Given CrsMatrix objects A, B and C, form the product C = A*B.
   In a parallel setting, A and B need not have matching distributions,
   but C needs to have the same row-map as A (if transposeA is false).
   At this time C=AT*B and C=A*BT are known to not work. However,
   C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. If Teuchos::null, a new CrsMatrix is created with optimal number of nnz per row.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

   */
  static RCP<Xpetra::Matrix<SC,LO,GO,NO> > Multiply(const Matrix& A, bool transposeA,
      const Matrix& B, bool transposeB,
      RCP<Matrix> C_in,
      Teuchos::FancyOStream& fos,
      bool doFillComplete           = true,
      bool doOptimizeStorage        = true,
      const std::string & label     = std::string()) {

    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    // Optimization using ML Multiply when available and requested
    // This feature is currently not supported. We would have to introduce the HAVE_XPETRA_ML_MMM flag
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT) && defined(HAVE_XPETRA_ML_MMM)
    if (B.getDomainMap()->lib() == Xpetra::UseEpetra && !transposeA && !transposeB) {
      RCP<const Epetra_CrsMatrix> epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2EpetraCrs(rcpFromRef(A));
      RCP<const Epetra_CrsMatrix> epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2EpetraCrs(rcpFromRef(B));
      RCP<Epetra_CrsMatrix>       epC = MLTwoMatrixMultiply(*epA, *epB, fos);

      RCP<Matrix> C = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<SC,LO,GO,NO> (epC);
      if (doFillComplete) {
        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage", doOptimizeStorage);
        C->fillComplete(B.getDomainMap(), A.getRangeMap(), params);
      }

      // Fill strided maps information
      // This is necessary since the ML matrix matrix multiplication routine has no handling for this
      // TODO: move this call to MLMultiply...
      C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);

      return C;
    }
#endif // EPETRA + EPETRAEXT + ML

    // Default case: Xpetra Multiply
    RCP<Matrix> C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow = Teuchos::as<double>(0);

      if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
        // For now, follow what ML and Epetra do.
        GO numRowsA = A.getGlobalNumRows();
        GO numRowsB = B.getGlobalNumRows();
        nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
            sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
        nnzPerRow *=  nnzPerRow;
        double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
        double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
        if (totalNnz < minNnz)
          totalNnz = minNnz;
        nnzPerRow = totalNnz / A.getGlobalNumRows();

        fos << "Matrix product nnz per row estimate = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
      }

      if (transposeA) C = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else            C = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(A.getRowMap(),    Teuchos::as<LO>(nnzPerRow));

    } else {
      C->resumeFill(); // why this is not done inside of Tpetra MxM?
      fos << "Reuse C pattern" << std::endl;
    }

    Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage,label); // call Multiply routine from above

    return C;
  }

  /*! @brief Helper function to do matrix-matrix multiply

  Returns C = AB.

  @param A left matrix
  @param transposeA if true, use the transpose of A
  @param B right matrix
  @param transposeB if true, use the transpose of B
  @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
   */
  static RCP<Xpetra::Matrix<SC,LO,GO,NO> > Multiply(const Matrix & A,
      bool transposeA,
      const Matrix & B,
      bool transposeB,
      Teuchos::FancyOStream &fos,
      bool callFillCompleteOnResult = true,
      bool doOptimizeStorage        = true,
      const std::string & label     = std::string()){
    return Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage,label);
  }

#ifdef HAVE_XPETRA_EPETRAEXT
  // Michael Gee's MLMultiply
  static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
      const Epetra_CrsMatrix& epB,
      Teuchos::FancyOStream& fos) {
#if defined(HAVE_XPETRA_ML_MMM)  // Note: this is currently not supported
    ML_Comm* comm;
    ML_Comm_Create(&comm);
    fos << "****** USING ML's MATRIX MATRIX MULTIPLY ******" << std::endl;
#ifdef HAVE_MPI
    // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
    const Epetra_MpiComm * Mcomm = dynamic_cast<const Epetra_MpiComm*>(&(epA.Comm()));
    if (Mcomm)
      ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
# endif
    //in order to use ML, there must be no indices missing from the matrix column maps.
    EpetraExt::CrsMatrix_SolverMap Atransform;
    EpetraExt::CrsMatrix_SolverMap Btransform;
    const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(epA));
    const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(epB));

    if (!A.Filled())    throw Exceptions::RuntimeError("A has to be FillCompleted");
    if (!B.Filled())    throw Exceptions::RuntimeError("B has to be FillCompleted");

    // create ML operators from EpetraCrsMatrix
    ML_Operator* ml_As = ML_Operator_Create(comm);
    ML_Operator* ml_Bs = ML_Operator_Create(comm);
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As); // Should we test if the lightweight wrapper is actually used or if WrapEpetraCrsMatrix fall backs to the heavy one?
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&B),ml_Bs);
    ML_Operator* ml_AtimesB = ML_Operator_Create(comm);
    {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ML_2matmult kernel"));
      ML_2matmult(ml_As,ml_Bs,ml_AtimesB,ML_CSR_MATRIX); // do NOT use ML_EpetraCRS_MATRIX!!!
    }
    ML_Operator_Destroy(&ml_As);
    ML_Operator_Destroy(&ml_Bs);

    // For ml_AtimesB we have to reconstruct the column map in global indexing,
    // The following is going down to the salt-mines of ML ...
    // note: we use integers, since ML only works for Epetra...
    int N_local = ml_AtimesB->invec_leng;
    ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
    if (!getrow_comm)   throw(Exceptions::RuntimeError("ML_Operator does not have a CommInfo"));
    ML_Comm* comm_AB = ml_AtimesB->comm;   // get comm object
    if (N_local != B.DomainMap().NumMyElements())
      throw(Exceptions::RuntimeError("Mismatch in local row dimension between ML and Epetra"));
    int N_rcvd = 0;
    int N_send = 0;
    int flag   = 0;
    for (int i = 0; i < getrow_comm->N_neighbors; i++) {
      N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
      N_send += (getrow_comm->neighbors)[i].N_send;
      if ( ((getrow_comm->neighbors)[i].N_rcv    != 0) &&
          ((getrow_comm->neighbors)[i].rcv_list != NULL) ) flag = 1;
    }
    // For some unknown reason, ML likes to have stuff one larger than
    // neccessary...
    std::vector<double> dtemp(N_local + N_rcvd + 1); // "double" vector for comm function
    std::vector<int>    cmap (N_local + N_rcvd + 1); // vector for gids
    for (int i = 0; i < N_local; ++i) {
      cmap[i]  = B.DomainMap().GID(i);
      dtemp[i] = (double) cmap[i];
    }
    ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm_AB); // do communication
    if (flag) { // process received data
      int count = N_local;
      const int neighbors = getrow_comm->N_neighbors;
      for (int i = 0; i < neighbors; i++) {
        const int nrcv = getrow_comm->neighbors[i].N_rcv;
        for (int j = 0; j < nrcv; j++)
          cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
      }
    } else {
      for (int i = 0; i < N_local+N_rcvd; ++i)
        cmap[i] = (int)dtemp[i];
    }
    dtemp.clear();  // free double array

    // we can now determine a matching column map for the result
    Epetra_Map gcmap(-1, N_local+N_rcvd, &cmap[0], B.ColMap().IndexBase(), A.Comm());

    int     allocated = 0;
    int     rowlength;
    double* val       = NULL;
    int*    bindx     = NULL;

    const int myrowlength    = A.RowMap().NumMyElements();
    const Epetra_Map& rowmap = A.RowMap();

    // Determine the maximum bandwith for the result matrix.
    // replaces the old, very(!) memory-consuming guess:
    // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
    int educatedguess = 0;
    for (int i = 0; i < myrowlength; ++i) {
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (rowlength>educatedguess)
        educatedguess = rowlength;
    }

    // allocate our result matrix and fill it
    RCP<Epetra_CrsMatrix> result = rcp(new Epetra_CrsMatrix(::Copy, A.RangeMap(), gcmap, educatedguess, false));

    std::vector<int> gcid(educatedguess);
    for (int i = 0; i < myrowlength; ++i) {
      const int grid = rowmap.GID(i);
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (!rowlength) continue;
      if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
      for (int j = 0; j < rowlength; ++j) {
        gcid[j] = gcmap.GID(bindx[j]);
        if (gcid[j] < 0)
          throw Exceptions::RuntimeError("Error: cannot find gcid!");
      }
      int err = result->InsertGlobalValues(grid, rowlength, val, &gcid[0]);
      if (err != 0 && err != 1) {
        std::ostringstream errStr;
        errStr << "Epetra_CrsMatrix::InsertGlobalValues returned err=" << err;
        throw Exceptions::RuntimeError(errStr.str());
      }
    }
    // free memory
    if (bindx) ML_free(bindx);
    if (val)   ML_free(val);
    ML_Operator_Destroy(&ml_AtimesB);
    ML_Comm_Destroy(&comm);

    return result;
#else // no MUELU_ML
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "No ML multiplication available. This feature is currently not supported by Xpetra.");
    return Teuchos::null;
#endif
  }
#endif //ifdef HAVE_XPETRA_EPETRAEXT

  /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
   */
  static RCP<Xpetra::BlockedCrsMatrix<SC,LO,GO,NO> > TwoMatrixMultiplyBlock(Xpetra::BlockedCrsMatrix<SC,LO,GO,NO>& A, bool transposeA,
      Xpetra::BlockedCrsMatrix<SC,LO,GO,NO>& B, bool transposeB,
      Teuchos::FancyOStream& fos,
      bool doFillComplete    = true,
      bool doOptimizeStorage = true) {
    if (transposeA || transposeB)
      throw Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

    // Preconditions
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    RCP<const Xpetra::MapExtractor<SC,LO,GO,NO> > rgmapextractor = A.getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<SC,LO,GO,NO> > domapextractor = B.getDomainMapExtractor();

    RCP<Xpetra::BlockedCrsMatrix<SC,LO,GO,NO> > C =
        rcp(new Xpetra::BlockedCrsMatrix<SC,LO,GO,NO>(rgmapextractor, domapextractor, 33 /* TODO fix me */));

    for (size_t i = 0; i < A.Rows(); ++i) { // loop over all block rows of A
      for (size_t j = 0; j < B.Cols(); ++j) { // loop over all block columns of B
        RCP<Matrix> Cij;

        for (size_t l = 0; l < B.Rows(); ++l) { // loop for calculating entry C_{ij}
          RCP<Xpetra::CrsMatrix<SC,LO,GO,NO> > crmat1 = A.getMatrix(i,l);
          RCP<Xpetra::CrsMatrix<SC,LO,GO,NO> > crmat2 = B.getMatrix(l,j);

          if (crmat1.is_null() || crmat2.is_null()) {
            continue;
          }

          RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crop1 = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO>(crmat1));
          RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crop2 = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO>(crmat2));

          RCP<Xpetra::Matrix<SC,LO,GO,NO> > temp = Multiply (*crop1, false, *crop2, false, fos);

          if (Cij.is_null ())
            Cij = temp;
          else
            Xpetra::MatrixMatrix<SC,LO,GO,NO>::TwoMatrixAdd (*temp, false, 1.0, *Cij, 1.0);
        }

        if (!Cij.is_null())  {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));

          RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crsCij = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(Cij);
          TEUCHOS_TEST_FOR_EXCEPTION(Cij.is_null(), Xpetra::Exceptions::BadCast,
              "MatrixFactory failed in generating a CrsMatrixWrap." );

          RCP<Xpetra::CrsMatrix<SC,LO,GO,NO> > crsMatCij = crsCij->getCrsMatrix();
          C->setMatrix(i, j, crsMatCij);

        } else {
          C->setMatrix(i, j, Teuchos::null);
        }
      }
    }

    if (doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  } // TwoMatrixMultiplyBlock

  /*! @brief Helper function to calculate B = alpha*A + beta*B.

  @param A      left matrix operand
  @param transposeA indicate whether to use transpose of A
  @param alpha  scalar multiplier for A
  @param B      right matrix operand
  @param beta   scalar multiplier for B

  @return sum in B.

  Note that B does not have to be fill-completed.
   */
  static void TwoMatrixAdd(const Xpetra::Matrix<SC,LO,GO,NO>& A, bool transposeA, SC alpha, Xpetra::Matrix<SC,LO,GO,NO>& B, SC beta) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2EpetraCrs(A);
      Epetra_CrsMatrix&       epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(B);

      //FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, beta);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value " + Teuchos::toString(rv));
      std::ostringstream buf;
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with EpetraExt.");
#endif
    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT
      const Tpetra::CrsMatrix<SC,LO,GO,NO>& tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
      Tpetra::CrsMatrix<SC,LO,GO,NO>& tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=int enabled."));
#endif
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
    }
  } //MatrixMatrix::TwoMatrixAdd()


  /*! @brief Helper function to calculate C = alpha*A + beta*B.

  @param A          left matrix operand
  @param transposeA indicate whether to use transpose of A
  @param alpha      scalar multiplier for A, defaults to 1.0
  @param B          right matrix operand
  @param transposeB indicate whether to use transpose of B
  @param beta       scalar multiplier for B, defaults to 1.0
  @param C          resulting sum

  It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
   */
  static void TwoMatrixAdd(const Xpetra::Matrix<SC,LO,GO,NO>& A, bool transposeA, const SC& alpha,
      const Xpetra::Matrix<SC,LO,GO,NO>& B, bool transposeB, const SC& beta,
      RCP<Xpetra::Matrix<SC,LO,GO,NO> >& C,  Teuchos::FancyOStream &fos, bool AHasFixedNnzPerRow = false) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (C == Teuchos::null) {
      if (!A.isFillComplete() || !B.isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Global statistics are not available for estimates.");

      size_t maxNzInA     = A.getGlobalMaxNumRowEntries();
      size_t maxNzInB     = B.getGlobalMaxNumRowEntries();
      size_t numLocalRows = A.getNodeNumRows();

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        // first check if either A or B has at most 1 nonzero per row
        // the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

        if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

        } else {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }

        fos << "MatrixMatrix::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
            << ", using static profiling" << std::endl;
        C = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO>(A.getRowMap(), exactNnzPerRow, Xpetra::StaticProfile));

      } else {
        // general case
        double nnzPerRowInA = Teuchos::as<double>(A.getGlobalNumEntries()) / A.getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B.getGlobalNumEntries()) / B.getGlobalNumRows();
        LO    nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        fos << "MatrixMatrix::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
            << ", max possible nnz per row in sum = " << maxPossible
            << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
            << std::endl;
        C = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO>(A.getRowMap(), nnzToAllocate, pft));
      }
      if (transposeB)
        fos << "MatrixMatrix::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2EpetraCrs(A);
      const Epetra_CrsMatrix& epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2EpetraCrs(B);
      RCP<Epetra_CrsMatrix>   epC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(C);
      Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

      //FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, transposeB, beta, ref2epC);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value of " + Teuchos::toString(rv));
#else
      throw Exceptions::RuntimeError("MueLu must be compile with EpetraExt.");
#endif
    } else if (C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT
      const Tpetra::CrsMatrix<SC,LO,GO,NO>& tpA =
          Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC,LO,GO,NO>& tpB =
          Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> >  tpC =
          Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, transposeB, beta, tpC);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=int enabled."));
#endif
#else
      throw Exceptions::RuntimeError("Xpetra must be compile with Tpetra.");
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
    if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
    ///////////////////////// EXPERIMENTAL

  } //MatrixMatrix::TwoMatrixAdd()
};

} // end namespace Xpetra

#define XPETRA_MATRIXMATRIX_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_HPP_ */
