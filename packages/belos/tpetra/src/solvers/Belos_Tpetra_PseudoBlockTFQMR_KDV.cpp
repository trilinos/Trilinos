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
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_PseudoBlockTFQMR_KDV_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra PseudoBlockTFQMRSolMgr<"
              << TypeNameTraits<SC>::name () << ", "
              << TypeNameTraits<LO>::name () << ", "
              << TypeNameTraits<GO>::name () << ", "
              << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "PSEUDOBLOCK TFQMR";
  {
    using DM = typename MV::wrapped_dual_view_type::DVT;
    using solver_type = ::Belos::PseudoBlockTFQMRSolMgr<SC, MV, OP, DM>;
    registerSolverSubclassForTypes<solver_type, SC, MV, OP, DM> (solverName);
  }
}

void register_PseudoBlockTFQMR_KDV (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV
#  undef BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV
#endif // BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV
#define BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV( SC, LO, GO, NT ) register_PseudoBlockTFQMR_KDV_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV )

#undef BELOS_TPETRA_REGISTER_PSEUDOBLOCKTFQMR_KDV
}

} // namespace Impl
} // namespace BelosTpetra
