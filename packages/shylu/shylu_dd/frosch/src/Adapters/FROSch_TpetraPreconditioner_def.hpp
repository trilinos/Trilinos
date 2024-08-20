// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TPETRAPRECONDITIONER_DEF_HPP
#define _FROSCH_TPETRAPRECONDITIONER_DEF_HPP

#if defined(HAVE_XPETRA_TPETRA)
#include <FROSch_TpetraPreconditioner_decl.hpp>

namespace FROSch {

    using namespace Teuchos;

    template <class SC,class LO,class GO,class NO>
    TpetraPreconditioner<SC,LO,GO,NO>::TpetraPreconditioner()
    {
    }

    template <class SC,class LO,class GO,class NO>
    TpetraPreconditioner<SC,LO,GO,NO>::TpetraPreconditioner(Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> > preconditioner) :
    preconditioner_(preconditioner)
    {
    }

    template <class SC,class LO,class GO,class NO>
    TpetraPreconditioner<SC,LO,GO,NO>::~TpetraPreconditioner()
    {
    }

    template <class SC,class LO,class GO,class NO>
    int TpetraPreconditioner<SC,LO,GO,NO>::
    initialize(bool useDefaultParameters)
    {
      return preconditioner_->initialize(useDefaultParameters);
    }

    template <class SC,class LO,class GO,class NO>
    int TpetraPreconditioner<SC,LO,GO,NO>::
    compute()
    {
      return preconditioner_->compute();
    }

    // Y = alpha * A^mode * X + beta * Y
    template <class SC,class LO,class GO,class NO>
    void TpetraPreconditioner<SC,LO,GO,NO>::
    apply(const Tpetra::MultiVector<SC,LO,GO,NO> &X,
                Tpetra::MultiVector<SC,LO,GO,NO> &Y,
          ETransp mode,
          SC alpha,
          SC beta) const
    {
      using XMultiVector = Xpetra::TpetraMultiVector<SC,LO,GO,NO>;
      auto xTmp = const_cast<Tpetra::MultiVector<SC,LO,GO,NO>&>(X);
      auto xRef = Xpetra::rcpFromRef(xTmp);
      const XMultiVector xX(xRef);
            XMultiVector xY(Xpetra::rcpFromRef(Y));
      preconditioner_->apply(xX, xY, mode, alpha, beta);
    }

    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > TpetraPreconditioner<SC,LO,GO,NO>::
    getDomainMap() const
    {
      return Xpetra::toTpetraNonZero(preconditioner_->getDomainMap());
    }

    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > TpetraPreconditioner<SC,LO,GO,NO>::
    getRangeMap() const
    {
      return Xpetra::toTpetraNonZero(preconditioner_->getRangeMap());
    }

    template <class SC,class LO,class GO,class NO>
    void TpetraPreconditioner<SC,LO,GO,NO>::
    describe(FancyOStream &out,
             const EVerbosityLevel verbLevel) const
    {
      return preconditioner_->describe(out, verbLevel);
    }

    template <class SC,class LO,class GO,class NO>
    string TpetraPreconditioner<SC,LO,GO,NO>::
    description() const
    {
      return preconditioner_->description();
    }

    template <class SC,class LO,class GO,class NO>
    bool TpetraPreconditioner<SC,LO,GO,NO>::isInitialized() const
    {
      return preconditioner_->isInitialized();
    }

    template <class SC,class LO,class GO,class NO>
    bool TpetraPreconditioner<SC,LO,GO,NO>::isComputed() const
    {
      return preconditioner_->isComputed();
    }

    template <class SC,class LO,class GO,class NO>
    void TpetraPreconditioner<SC,LO,GO,NO>::residual(const TMultiVector & X,
                                                     const TMultiVector & B,
                                                           TMultiVector& R) const
    {
      SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one;
      apply(X,R);
      R.update(one,B,negone);
    }

    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> >
    TpetraPreconditioner<SC,LO,GO,NO>::getSchwarzPreconditioner()
    {
      return preconditioner_;
    }
}
#endif

#endif
