// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STRATIMIKOS_BELOS_PREC_TPETRA_HELPERS_HPP
#define STRATIMIKOS_BELOS_PREC_TPETRA_HELPERS_HPP

#include "Stratimikos_LinearSolverBuilder.hpp"
#include "Thyra_BelosTpetraPreconditionerFactory_decl.hpp"
#include "Thyra_BelosTpetraPreconditionerFactory_def.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>

namespace Stratimikos {

  template <typename MatrixType>
  void enableBelosPrecTpetra(LinearSolverBuilder<typename MatrixType::scalar_type>& builder, const std::string& stratName = "BelosPrecTpetra")
  {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableBelosPrecTpetra cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef typename MatrixType::scalar_type scalar_type;
    typedef Thyra::PreconditionerFactoryBase<scalar_type>       Base;
    typedef Thyra::BelosTpetraPreconditionerFactory<MatrixType> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }

} // namespace Stratimikos

#endif
