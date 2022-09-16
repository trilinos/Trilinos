//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

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
    TpetraPreconditioner<SC,LO,GO,NO>::TpetraPreconditioner(Teuchos::RCP<SchwarzPreconditioner<SC,LO,GO,NO> > preconditioner) :
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
    Teuchos::RCP<SchwarzPreconditioner<SC,LO,GO,NO> >
    TpetraPreconditioner<SC,LO,GO,NO>::getSchwarzPreconditioner()
    {
      return preconditioner_;
    }
}
#endif

#endif
