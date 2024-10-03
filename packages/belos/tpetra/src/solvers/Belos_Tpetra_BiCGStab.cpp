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
#include "BelosBiCGStabSolMgr.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_BiCGStab_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;
  using solver_type = ::Belos::BiCGStabSolMgr<SC, MV, OP>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra BiCGStabSolMgr<"
              << TypeNameTraits<SC>::name () << ", "
              << TypeNameTraits<LO>::name () << ", "
              << TypeNameTraits<GO>::name () << ", "
              << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "BICGSTAB";
  registerSolverSubclassForTypes<solver_type, SC, MV, OP> (solverName);
}

void register_BiCGStab (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_BICGSTAB
#  undef BELOS_TPETRA_REGISTER_BICGSTAB
#endif // BELOS_TPETRA_REGISTER_BICGSTAB
#define BELOS_TPETRA_REGISTER_BICGSTAB( SC, LO, GO, NT ) register_BiCGStab_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_BICGSTAB )

#undef BELOS_TPETRA_REGISTER_BICGSTAB
}

} // namespace Impl
} // namespace BelosTpetra
