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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef MUELU_TPETRAOPERATOR_DEF_HPP
#define MUELU_TPETRAOPERATOR_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include "MueLu_TpetraOperator_decl.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;

  RCP<const Map> domainMap;
  if (!Hierarchy_.is_null())
    domainMap = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A")->getDomainMap();
  else
    domainMap = Operator_->getDomainMap();

  RCP<const BlockedMap> bDomainMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(domainMap);
  if (bDomainMap.is_null() == false) {
    return Xpetra::toTpetraNonZero(bDomainMap->getFullMap());
  }
  return Xpetra::toTpetraNonZero(domainMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;

  RCP<const Map> rangeMap;
  if (!Hierarchy_.is_null())
    rangeMap = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A")->getRangeMap();
  else
    rangeMap = Operator_->getRangeMap();

  RCP<const BlockedMap> bRangeMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rangeMap);
  if (bRangeMap.is_null() == false) {
    return Xpetra::toTpetraNonZero(bRangeMap->getFullMap());
  }
  return Xpetra::toTpetraNonZero(rangeMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                                                                      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                                                                      Teuchos::ETransp mode, Scalar /* alpha */, Scalar /* beta */) const {
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TMV;
  typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTMV;

  TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS, std::logic_error, "MueLu::TpetraOperator does not support applying the adjoint operator");

  try {
    TMV& temp_x = const_cast<TMV&>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV tY(rcpFromRef(Y));

    if (!Hierarchy_.is_null())
      Hierarchy_->Iterate(tX, tY, 1, true);
    else
      Operator_->apply(tX, tY);

  } catch (std::exception& e) {
    std::cerr << "MueLu::TpetraOperator::apply : detected an exception" << std::endl
              << e.what() << std::endl;
    throw;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasTransposeApply() const {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetHierarchy() const {
  return Hierarchy_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetOperator() const {
  return Operator_;
}

}  // namespace MueLu

#endif  // ifdef MUELU_TPETRAOPERATOR_DEF_HPP
