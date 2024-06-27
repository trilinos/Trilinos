// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_TPETRA_HPP
#define BELOSSOLVERFACTORY_TPETRA_HPP

#include "Belos_Details_Tpetra_registerSolverFactory.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOperatorTraits_Tpetra.hpp"

namespace Belos {

template<class SC, class MV, class OP>
class TpetraSolverFactory : public Impl::SolverFactoryParent<SC, MV, OP>
{
  public:
    TpetraSolverFactory() {
      Details::Tpetra::registerSolverFactory();
    };
};

namespace Impl {

template<class SC, class LO, class GO, class NT>
class SolverFactorySelector<SC,Tpetra::MultiVector<SC, LO, GO, NT>,Tpetra::Operator<SC, LO, GO, NT>> {
  public:
    typedef TpetraSolverFactory<SC,Tpetra::MultiVector<SC, LO, GO, NT>,Tpetra::Operator<SC, LO, GO, NT>> type;
};

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_TPETRA_HPP
