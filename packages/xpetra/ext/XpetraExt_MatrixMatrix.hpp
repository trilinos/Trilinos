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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * XpetraExt_MatrixMatrix_decl.hpp
 *
 *  Created on: Oct 10, 2012
 *      Author: wiesner
 */

#ifndef XPETRAEXT_MATRIXMATRIX_DECL_HPP_
#define XPETRAEXT_MATRIXMATRIX_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "XpetraExt_MatrixMatrix.hpp"

#include "Xpetra_Matrix.hpp"

#ifdef HAVE_XPETRA_EPETRA
namespace Xpetra {
  class EpetraCrsMatrix; // TODO: replace by include of _fwd.hpp
  //  class
}
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
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_XPETRA_TPETRA

/*! \file Xpetra_MatrixMatrix.hpp

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Xpetra {

namespace MatrixMatrix {

  // TODO: code factorization with MueLu_Utilities for conversion functions.
  // TODO: _decl and _def

#ifdef HAVE_XPETRA_EPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op) {
  RCP<const Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Op);
  if (crsOp == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp->getCrsMatrix();
  const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
  if (tmp_ECrsMtx == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
  A = tmp_ECrsMtx->getEpetra_CrsMatrix();
  return A;
} //Op2EpetraCrs


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op) {
  RCP<Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Op);
  if (crsOp == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp->getCrsMatrix();
  const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
  if (tmp_ECrsMtx == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
  A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
  return A;
} //Op2NonConstEpetraCrs

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Epetra_CrsMatrix & Op2EpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op) {
  RCP<const Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  try {
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & >(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp.getCrsMatrix();
    const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
    return *A;
  } catch(...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
} //Op2EpetraCrs


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Epetra_CrsMatrix & Op2NonConstEpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op) {
  RCP<Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  try {
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & >(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp.getCrsMatrix();
    const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(tmp_CrsMtx);
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
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op2TpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op) {
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
  // Get the underlying Tpetra Mtx
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Op);
  if (crsOp == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp->getCrsMatrix();
  const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(tmp_CrsMtx);
  if (tmp_ECrsMtx == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
  A = tmp_ECrsMtx->getTpetra_CrsMatrix();
  return A;
} //Op2TpetraCrs


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op) {
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
  // Get the underlying Tpetra Mtx
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > crsOp = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Op);
  if (crsOp == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp->getCrsMatrix();
  const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(tmp_CrsMtx);
  if (tmp_ECrsMtx == Teuchos::null)
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
  A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
  return A;
} //Op2NonConstTpetraCrs

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op2TpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op) {
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
  // Get the underlying Tpetra Mtx
  try{
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & crsOp =dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & >(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp.getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(tmp_CrsMtx);
    if (tmp_TCrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_TCrsMtx->getTpetra_CrsMatrix();
    return *A;
  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
} //Op2TpetraCrs

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op2NonConstTpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & Op) {
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
  // Get the underlying Tpetra Mtx
  try{
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & >(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tmp_CrsMtx = crsOp.getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(tmp_CrsMtx);
    if (tmp_TCrsMtx == Teuchos::null)
      throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(tmp_TCrsMtx->getTpetra_CrsMatrix());
    return *A;
  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
} //Op2NonConstTpetraCrs

#endif // HAVE_XPETRA_TPETRA

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
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void Multiply(
  const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& A,
  bool transposeA,
  const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& B,
  bool transposeB,
  Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& C,
  bool call_FillComplete_on_result = true,
  bool doOptimizeStorage = true) {

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
      Epetra_CrsMatrix & epA = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(A);
      Epetra_CrsMatrix & epB = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(B);
      Epetra_CrsMatrix & epC = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(C);


      int i = EpetraExt::MatrixMatrix::Multiply(epA,transposeA,epB,transposeB,epC,haveMultiplyDoFillComplete);
      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }

#endif
    } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & tpA = Xpetra::MatrixMatrix::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & tpB = Xpetra::MatrixMatrix::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> &       tpC = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(C);

      //18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
      //Previously, Tpetra's matrix matrix multiply did not support fillComplete.
      Tpetra::MatrixMatrix::Multiply(tpA,transposeA,tpB,transposeB,tpC,haveMultiplyDoFillComplete);
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
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Teuchos::rcpFromRef(A));
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB); // TODO use references instead of RCPs
} // end Multiply

#if TO_BE_FIXED
/** Given CrsMatrix objects A and B, form the sum B = a*A + b*B
 * Currently not functional.

@param A Input, must already have had 'FillComplete()' called.
@param transposeA Input, whether to use transpose of matrix A.
@param scalarA Input, scalar multiplier for matrix A.
@param B Result. On entry to this method, it doesn't matter whether
         FillComplete() has already been called on B or not. If it has,
   then B's graph must already contain all nonzero locations that
   will be produced when forming the sum.
@param scalarB Input, scalar multiplier for matrix B.

 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void Add(
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& A,
    bool transposeA,
    Scalar scalarA,
    Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& B,
    Scalar scalarB ) {

  if (!(A.getRowMap()->isSameAs(*(B.getRowMap())))) {
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Add: matrix row maps are not the same."));
  }

  /*RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(A));
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(B));*/

  if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_XPETRA_EPETRAEXT
    const Epetra_CrsMatrix & epA = Xpetra::MatrixMatrix::Op2EpetraCrs(A);
    Epetra_CrsMatrix & epB = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(B);

    //FIXME is there a bug if beta=0?
    int i = EpetraExt::MatrixMatrix::Add(epA,transposeA,scalarA,epB,scalarB);
    if (i != 0) {
      std::ostringstream buf;
      buf << i;
      std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
      throw(Xpetra::Exceptions::RuntimeError(msg));
    }
#else
    throw(Exceptions::RuntimeError("Xpetra must be compiled with EpetraExt."));
#endif
  } else if(A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    //RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA = Xpetra::MatrixMatrix::Op2TpetraCrs(rcpA);
    //RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpB = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(rcpB);

    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & tpA = Xpetra::MatrixMatrix::Op2TpetraCrs(A);
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> &       tpB = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(B);

    //Tpetra::MatrixMatrix::Add(*tpA, transposeA, scalarA, *tpB, scalarB);
    Tpetra::MatrixMatrix::Add(tpA, transposeA, scalarA, tpB, scalarB);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }
} // end Add

/** Given CrsMatrix objects A and B, form the sum C = a*A + b*B

@param A Input, must already have had 'FillComplete()' called.
@param transposeA Input, whether to use transpose of matrix A.
@param scalarA Input, scalar multiplier for matrix A.
@param B Input, must already have had 'FillComplete()' called.
@param transposeB Input, whether to use transpose of matrix B.
@param scalarB Input, scalar multiplier for matrix B.
@param C Result. On entry to this method, C can be NULL or a pointer
         to an unfilled or filled matrix. If C is NULL then a new
         object is allocated and must be deleted by the user.
         If C is not NULL and FillComplete has already been called then the sparsity pattern is assumed to be fixed and compatible  with the sparsity of A+B. If FillComplete has not been called then the sum is completed and the function
         returns without calling FillComplete on C.

 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void Add(
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& A,
    bool transposeA,
    Scalar scalarA,
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& B,
    bool transposeB,
    Scalar scalarB,
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > C) {

  if ( !(A.getRowMap()->isSameAs(*(B.getRowMap()))) ) {
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Add: matrix row maps are not the same."));
  }
  if (C==Teuchos::null)
    //FIXME 5 is a complete guess as to the #nonzeros per row
    C = rcp( new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(A.getRowMap(), 5) );

  if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_XPETRA_EPETRAEXT
      const Epetra_CrsMatrix& epA = Xpetra::MatrixMatrix::Op2EpetraCrs(A);
      const Epetra_CrsMatrix& epB = Xpetra::MatrixMatrix::Op2EpetraCrs(B);
      RCP<Epetra_CrsMatrix>       epC = Xpetra::MatrixMatrix::Op2NonConstEpetraCrs(C);
      Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

      //FIXME is there a bug if beta=0?
      int i = EpetraExt::MatrixMatrix::Add(epA,transposeA,scalarA,epB,transposeB,scalarB,ref2epC);

      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
        throw(Xpetra::Exceptions::RuntimeError(msg));
      }
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compile with EpetraExt."));
#endif
  } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & tpA = Xpetra::MatrixMatrix::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & tpB = Xpetra::MatrixMatrix::Op2TpetraCrs(B);
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >       tpC = Xpetra::MatrixMatrix::Op2NonConstTpetraCrs(C);

    Tpetra::MatrixMatrix::Add(tpA, transposeA, scalarA, tpB, transposeB, scalarB, tpC);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compile with Tpetra."));
#endif
  }

  ///////////////////////// EXPERIMENTAL
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Teuchos::rcpFromRef(A));
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(Teuchos::rcpFromRef(B));
  if(A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpA);
  if(B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpB);
  ///////////////////////// EXPERIMENTAL
}
#endif // TO_BE_FIXED


} // end namespace MatrixMatrix

} // end namespace Xpetra


#endif /* XPETRAEXT_MATRIXMATRIX_DECL_HPP_ */
