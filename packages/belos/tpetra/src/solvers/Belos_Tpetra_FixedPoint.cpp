// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosSolverFactory.hpp"
#include "BelosSolverFactory_Tpetra.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_FixedPoint_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;
  using solver_type = ::Belos::FixedPointSolMgr<SC, MV, OP>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra FixedPointSolMgr<"
              << TypeNameTraits<SC>::name () << ", "
              << TypeNameTraits<LO>::name () << ", "
              << TypeNameTraits<GO>::name () << ", "
              << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "FIXED POINT";
  registerSolverSubclassForTypes<solver_type, SC, MV, OP> (solverName);
}

void register_FixedPoint (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_FIXEDPOINT
#  undef BELOS_TPETRA_REGISTER_FIXEDPOINT
#endif // BELOS_TPETRA_REGISTER_FIXEDPOINT
#define BELOS_TPETRA_REGISTER_FIXEDPOINT( SC, LO, GO, NT ) register_FixedPoint_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_FIXEDPOINT )

#undef BELOS_TPETRA_REGISTER_FIXEDPOINT
}

} // namespace Impl
} // namespace BelosTpetra
