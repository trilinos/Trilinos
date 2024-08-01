// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosSolverFactory.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOperatorTraits_Tpetra.hpp"
#include "Belos_Tpetra_Cg.hpp"
#include "Belos_Tpetra_SolverManager.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_Cg_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;
  using solver_type = CgSolverManager<SC, MV, OP>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra CgSolverManager<"
	      << TypeNameTraits<SC>::name () << ", "
	      << TypeNameTraits<LO>::name () << ", "
	      << TypeNameTraits<GO>::name () << ", "
	      << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "TPETRA CG";
  registerSolverSubclassForTypes<solver_type, SC, MV, OP> (solverName);
}

void register_Cg (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_CG
#  undef BELOS_TPETRA_REGISTER_CG
#endif // BELOS_TPETRA_REGISTER_CG
#define BELOS_TPETRA_REGISTER_CG( SC, LO, GO, NT ) register_Cg_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_CG )

#undef BELOS_TPETRA_REGISTER_CG
}

} // namespace Impl
} // namespace BelosTpetra
