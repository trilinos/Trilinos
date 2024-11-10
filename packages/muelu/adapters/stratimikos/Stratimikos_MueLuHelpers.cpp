// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stratimikos_MueLuHelpers.hpp"

#include "MueLu_ConfigDefs.hpp"

#include "Thyra_MueLuPreconditionerFactory.hpp"
#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#include "Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

namespace Stratimikos {

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#if 0
  void enableMueLuTpetraQ2Q1(DefaultLinearSolverBuilder &builder, const std::string &stratName) {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableMueLuTpetraQ2Q1 cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef Thyra::PreconditionerFactoryBase<double>                      Base;
    typedef Thyra::MueLuTpetraQ2Q1PreconditionerFactory<double, int, int> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }
#endif
#endif

}  // namespace Stratimikos
