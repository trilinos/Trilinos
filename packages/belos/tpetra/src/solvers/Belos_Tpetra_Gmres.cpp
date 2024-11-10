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
#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_SolverManager.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_Gmres_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;
  using solver_type = GmresSolverManager<SC, MV, OP>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra GmresSolverManager<"
	      << TypeNameTraits<SC>::name () << ", "
	      << TypeNameTraits<LO>::name () << ", "
	      << TypeNameTraits<GO>::name () << ", "
	      << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "TPETRA GMRES";
  registerSolverSubclassForTypes<solver_type, SC, MV, OP> (solverName);
}

void register_Gmres (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_GMRES
#  undef BELOS_TPETRA_REGISTER_GMRES
#endif // BELOS_TPETRA_REGISTER_GMRES
#define BELOS_TPETRA_REGISTER_GMRES( SC, LO, GO, NT ) register_Gmres_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_GMRES )

#undef BELOS_TPETRA_REGISTER_GMRES
}

} // namespace Impl
} // namespace BelosTpetra
