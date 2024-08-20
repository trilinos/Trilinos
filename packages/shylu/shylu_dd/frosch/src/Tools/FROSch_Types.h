// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TYPES_H
#define _FROSCH_TYPES_H

#include <ShyLU_DDFROSch_config.h>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    #if defined HAVE_XPETRA_EPETRA || defined HAVE_TPETRA_INT_INT
    typedef int DefaultGlobalOrdinal;
    #elif !defined HAVE_TPETRA_INT_LONG_LONG
    typedef long DefaultGlobalOrdinal;
    #else
    typedef long long DefaultGlobalOrdinal;
    #endif

    enum DofOrdering {NodeWise=0,DimensionWise=1,Custom=2};

    enum class NullSpaceType
    {
      Laplace = 0,
      Elasticity = 1
    };

    enum Verbosity {None=0,All=1};

}

#endif
