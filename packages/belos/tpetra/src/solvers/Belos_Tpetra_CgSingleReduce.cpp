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
#include "Belos_Tpetra_CgSingleReduce.hpp"
#include "Belos_Tpetra_SolverManager.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_CgSingleReduce_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;
  using solver_type = CgSingleReduceSolverManager<SC, MV, OP>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra CgSingleReduceSolverManager<"
	      << TypeNameTraits<SC>::name () << ", "
	      << TypeNameTraits<LO>::name () << ", "
	      << TypeNameTraits<GO>::name () << ", "
	      << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "TPETRA CG SINGLE REDUCE";
  registerSolverSubclassForTypes<solver_type, SC, MV, OP> (solverName);
}

void register_CgSingleReduce (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE
#  undef BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE
#endif // BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE
#define BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE( SC, LO, GO, NT ) register_CgSingleReduce_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE )

#undef BELOS_TPETRA_REGISTER_CG_SINGLE_REDUCE
}

} // namespace Impl
} // namespace BelosTpetra
