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
#ifndef XPETRA_TPETRAHALFPRECISIONOPEARTOR_HPP
#define XPETRA_TPETRAHALFPRECISIONOPEARTOR_HPP

#include "Xpetra_ConfigDefs.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

namespace Xpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<typename Teuchos::ScalarTraits<Scalar>::halfPrecision, LocalOrdinal, GlobalOrdinal, Node> >
  convertToHalfPrecision(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A) {
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
    using HalfScalar = typename Teuchos::ScalarTraits<Scalar>::halfPrecision;
    typedef typename Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpCrs;

    RCP<XpCrs> xpCrs = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(A, true)->getCrsMatrix();
    auto tpCrs = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(xpCrs, true)->getTpetra_CrsMatrix();
    auto newTpCrs = tpCrs->template convert<HalfScalar>();
    auto newXpCrs = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>(newTpCrs));
    auto newA = Teuchos::rcp(new Xpetra::CrsMatrixWrap<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>(Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>>(newXpCrs)));

    return newA;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
      "Xpetra::convertToHalfPrecision only available for Tpetra with SC=double and SC=float enabled");
    TEUCHOS_UNREACHABLE_RETURN(false);
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::halfPrecision, LocalOrdinal, GlobalOrdinal, Node> >
  convertToHalfPrecision(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X) {
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
    using HalfScalar = typename Teuchos::ScalarTraits<Scalar>::halfPrecision;

    auto tpX = toTpetra(*X);
    auto newTpX = tpX.template convert<HalfScalar>();
    auto newX = rcp(new Xpetra::TpetraMultiVector<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>(newTpX));

    return newX;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
      "Xpetra::convertToHalfPrecision only available for Tpetra with SC=double and SC=float enabled");
    TEUCHOS_UNREACHABLE_RETURN(false);
#endif
  }

/*!  @brief Wraps an existing halfer precision Xpetra::Operator as a Xpetra::Operator.
*/
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class TpetraHalfPrecisionOperator : public Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::halfPrecision HalfScalar;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    TpetraHalfPrecisionOperator(const RCP<Xpetra::Operator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> >& op) : Op_(op) {
      Allocate(1);
    }

    void Allocate(int numVecs) {
      X_ = Xpetra::MultiVectorFactory<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Op_->getDomainMap(),numVecs);
      Y_ = Xpetra::MultiVectorFactory<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Op_->getRangeMap(),numVecs);
    }

    //! Destructor.
    virtual ~TpetraHalfPrecisionOperator() { }

    //@}

    //! Returns the Tpetra::Map object associated with the domain of this TpetraOperator.
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      return Op_->getDomainMap();
    }

    //! Returns the Tpetra::Map object associated with the range of this TpetraOperator.
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      return Op_->getRangeMap();
    }

    //! Returns in Y the result of a Xpetra::TpetraOperator applied to a Xpetra::MultiVector X.
    /*!
      \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta = Teuchos::ScalarTraits<Scalar>::one()) const{
      typedef Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>     tMV;
      typedef Xpetra::TpetraMultiVector<HalfScalar,LocalOrdinal,GlobalOrdinal,Node> tMVHalf;
      Tpetra::deep_copy(*Teuchos::rcp_dynamic_cast<tMVHalf>(X_)->getTpetra_MultiVector(),
                        toTpetra(X));
      Op_->apply(*X_, *Y_, mode, Teuchos::as<HalfScalar>(alpha), Teuchos::as<HalfScalar>(beta));
      Tpetra::deep_copy(toTpetra(Y),
                        *Teuchos::rcp_dynamic_cast<tMVHalf>(Y_)->getTpetra_MultiVector());
    }

    //! Indicates whether this TpetraOperator supports applying the adjoint TpetraOperator.
    bool hasTransposeApply() const { return false; }

    //! Compute a residual R = B - (*this) * X
    void residual(const Xpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & X,
                  const Xpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & B,
                  Xpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & R) const {
      using STS = Teuchos::ScalarTraits<Scalar>;
      R.update(STS::one(),B,STS::zero());
      this->apply (X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
    }

    //! @name Xpetra specific
    //@{

    //! Direct access to the underlying TpetraOperator.
    RCP<Xpetra::Operator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> > GetHalfPrecisionOperator() const { return Op_; }

    void SetHalfPrecisionOperator(const RCP<Xpetra::Operator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> >& op) { Op_ = op; };

    //@}

  private:
    RCP<Xpetra::Operator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> > Op_;
    RCP<Xpetra::MultiVector<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> > X_, Y_;
  };

} // namespace

#endif // XPETRA_TPETRAHALFPRECISIONOPEARTOR_HPP
