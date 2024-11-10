// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_EPETRA_HPP
#define BELOSSOLVERFACTORY_EPETRA_HPP

#include "Belos_Details_Epetra_registerSolverFactory.hpp"
#include "BelosSolverFactory.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

namespace Belos {

class EpetraSolverFactory : public Impl::SolverFactoryParent<double, Epetra_MultiVector, Epetra_Operator>
{
  public:
    EpetraSolverFactory() {
      Details::Epetra::registerSolverFactory();
    };
};

namespace Impl {

template<>
class SolverFactorySelector<double, Epetra_MultiVector, Epetra_Operator> {
  public:
    typedef EpetraSolverFactory type;
};

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_EPETRA_HPP
