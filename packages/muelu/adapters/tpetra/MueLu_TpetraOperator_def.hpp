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

#ifdef HAVE_MUELU_TPETRA

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include "MueLu_TpetraOperator_decl.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

  RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
  return Xpetra::toTpetraNonZero(A->getDomainMap());
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

  RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
  return Xpetra::toTpetraNonZero(A->getRangeMap());
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                                                               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                                                               Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>       TMV;
  typedef Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> XTMV;

  try {
    TMV& temp_x = const_cast<TMV &>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV       tY(rcpFromRef(Y));

    tY.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    Hierarchy_->Iterate(tX, tY, 1, true);

  } catch (std::exception& e) {
    //FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::TpetraOperator::ApplyInverse():" << std::endl
        << e.what() << std::endl;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasTransposeApply() const {
  return false;
}

} // namespace
#endif //ifdef HAVE_MUELU_TPETRA

#endif //ifdef MUELU_TPETRAOPERATOR_DEF_HPP
