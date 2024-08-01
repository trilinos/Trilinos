// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_XPETRA_HPP
#define BELOSSOLVERFACTORY_XPETRA_HPP

#include "Belos_Details_Xpetra_registerSolverFactory.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosXpetraAdapterMultiVector.hpp"
#include "BelosOperatorT.hpp"

namespace Belos {

template<class Scalar, class MV, class OP>
class XpetraSolverFactory : public Impl::SolverFactoryParent<Scalar, MV, OP>
{
  public:
    XpetraSolverFactory() {
      Details::Xpetra::registerSolverFactory();
    };
};

namespace Impl {

template<class SC, class LO, class GO, class NT>
class SolverFactorySelector<SC,::Xpetra::MultiVector<SC, LO, GO, NT>,::Belos::OperatorT<Xpetra::MultiVector<SC, LO, GO, NT>>> {
  public:
    typedef XpetraSolverFactory<SC,::Xpetra::MultiVector<SC, LO, GO, NT>,::Belos::OperatorT<Xpetra::MultiVector<SC, LO, GO, NT>>> type;
};

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_XPETRA_HPP
