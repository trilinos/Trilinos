// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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

#ifndef MUELU_TPETRAOPERATOR_DEF_HPP 
#define MUELU_TPETRAOPERATOR_DEF_HPP 

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_TPETRA

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include "MueLu_TpetraOperator_decl.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

// ------------- getDomainMap -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const {

  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> XMatrix;

  RCP<MueLu::Level>  L0 = Hierarchy_->GetLevel(0);
  RCP<XMatrix> A = L0->Get< RCP<XMatrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(A);
  if(tpbA != Teuchos::null)
    return Xpetra::toTpetra(tpbA->getDomainMap());

  RCP< Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(A);
  return tpA->getDomainMap();
}

// ------------- getRangeMap -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const {

  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> XMatrix;

  RCP<MueLu::Level>  L0 = Hierarchy_->GetLevel(0);
  RCP<XMatrix> A = L0->Get< RCP<XMatrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(A);
  if(tpbA != Teuchos::null)
    return Xpetra::toTpetra(tpbA->getRangeMap());

  RCP< Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(A);
  return tpA->getRangeMap();
}

// ------------- apply -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                                                               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                                                               Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {

  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV; 
  typedef Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> XTMV;

  try {

    TMV & temp_x = const_cast<TMV &>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV       tY(rcpFromRef(Y));

    tY.putScalar(0.0);
    Hierarchy_->Iterate(tX, 1, tY, true);
  }
  
  catch(std::exception& e) {
    //FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::TpetraOperator::ApplyInverse():" << std::endl
        << e.what() << std::endl;
  }
}

// ------------- apply -----------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const {
  return false;
}

} // namespace
#endif //ifdef HAVE_MUELU_TPETRA

#endif //ifdef MUELU_TPETRAOPERATOR_DEF_HPP 
