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
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>

namespace BelosTpetra {
namespace Impl {

template<class SC, class LO, class GO, class NT>
void register_PseudoBlockCG_tmpl (const bool verbose)
{
  using ::Belos::Impl::registerSolverSubclassForTypes;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = ::Tpetra::Operator<SC, LO, GO, NT>;

  if (verbose) {
    using Teuchos::TypeNameTraits;
    std::cout << "Registering BelosTpetra PseudoBlockCGSolMgr<"
              << TypeNameTraits<SC>::name () << ", "
              << TypeNameTraits<LO>::name () << ", "
              << TypeNameTraits<GO>::name () << ", "
              << TypeNameTraits<NT>::name () << ">" << std::endl;
  }
  const char solverName[] = "PSEUDOBLOCK CG";
  {
    using DM = ::Teuchos::SerialDenseMatrix<int, SC>;
    using solver_type = ::Belos::PseudoBlockCGSolMgr<SC, MV, OP, DM>;
    registerSolverSubclassForTypes<solver_type, SC, MV, OP, DM> (solverName);
  }

  {
    using DM = ::Kokkos::DualView<typename KokkosKernels::ArithTraits<SC>::val_type **, Kokkos::LayoutLeft>;
    using solver_type = ::Belos::PseudoBlockCGSolMgr<SC, MV, OP, DM>;
    registerSolverSubclassForTypes<solver_type, SC, MV, OP, DM> (solverName);
  }
}

void register_PseudoBlockCG (const bool verbose)
{
  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG
#  undef BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG
#endif // BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG
#define BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG( SC, LO, GO, NT ) register_PseudoBlockCG_tmpl<SC, LO, GO, NT> (verbose);

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG )

#undef BELOS_TPETRA_REGISTER_PSEUDOBLOCKCG
}

} // namespace Impl
} // namespace BelosTpetra
