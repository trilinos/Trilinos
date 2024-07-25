// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _STRATIMIKOS_FROSCH_DECL_HPP
#define _STRATIMIKOS_FROSCH_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#ifdef HAVE_SHYLU_DDFROSCH_STRATIMIKOS
#include "Stratimikos_LinearSolverBuilder.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"


namespace Thyra {

    template <class SC,
              class LO,
              class GO,
              class NO>
    class FROSchFactory;
}

namespace Stratimikos {

    using namespace std;
    using namespace Teuchos;
    using namespace Thyra;

    template <typename LO,
              typename GO,
              typename NO>
    void enableFROSch(LinearSolverBuilder<double>& builder,
                      const string& stratName = "FROSch");

    template <typename SC = double,
              typename LO = int,
              typename GO = int,
              typename NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    void enableFROSch(LinearSolverBuilder<SC>& builder,
                      const string& stratName = "FROSch");
} // namespace Stratimikos

#endif

#endif
