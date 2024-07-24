// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _STRATIMIKOS_FROSCH_DEF_HPP
#define _STRATIMIKOS_FROSCH_DEF_HPP

#include "Stratimikos_FROSch_decl.hpp"

#include "Thyra_FROSchFactory_def.hpp"


namespace Stratimikos {

  using namespace std;
  using namespace Teuchos;
  using namespace Thyra;

  template <typename LO,typename GO,typename NO>
  void enableFROSch (LinearSolverBuilder<double>& builder,
                     const string& stratName)
  {
    const RCP<const ParameterList> precValidParams = sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), logic_error,
                               "Stratimikos::enableFROSch cannot add \"" + stratName +"\" because it is already included in builder!");

    using Base = PreconditionerFactoryBase<double>;
    if (!stratName.compare("FROSch")) {
      using Impl = FROSchFactory<double, LO, GO, NO>;
      builder.setPreconditioningStrategyFactory(abstractFactoryStd<Base, Impl>(), stratName);
    }
  }

  template <typename SC, typename LO,typename GO,typename NO>
  void enableFROSch (LinearSolverBuilder<SC>& builder,
                     const string& stratName)
  {
    const RCP<const ParameterList> precValidParams = sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), logic_error,
                               "Stratimikos::enableFROSch cannot add \"" + stratName +"\" because it is already included in builder!");

    using Base = PreconditionerFactoryBase<SC>;
    if (!stratName.compare("FROSch")) {
      using Impl = FROSchFactory<SC, LO, GO, NO>;
      builder.setPreconditioningStrategyFactory(abstractFactoryStd<Base, Impl>(), stratName);
    }
  }

} // namespace Stratimikos

#endif
