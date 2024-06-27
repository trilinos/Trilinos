// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/// \file Belos_Details_Tpetra_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Tpetra

#include "BelosSolverFactory_Tpetra.hpp"

namespace BelosTpetra {
namespace Impl {

extern void register_BiCGStab (const bool verbose);
extern void register_BlockCG (const bool verbose);
extern void register_BlockGmres (const bool verbose);
extern void register_Cg (const bool verbose);
extern void register_CgPipeline (const bool verbose);
extern void register_CgSingleReduce (const bool verbose);
extern void register_FixedPoint (const bool verbose);
extern void register_GCRODR (const bool verbose);
extern void register_Gmres (const bool verbose);
extern void register_GmresPipeline (const bool verbose);
extern void register_GmresPoly (const bool verbose);
extern void register_GmresSingleReduce (const bool verbose);
extern void register_GmresSstep (const bool verbose);
extern void register_LSQR (const bool verbose);
extern void register_Minres (const bool verbose);
extern void register_PCPG (const bool verbose);
extern void register_PseudoBlockCG (const bool verbose);
extern void register_PseudoBlockGmres (const bool verbose);
extern void register_PseudoBlockTFQMR (const bool verbose);
extern void register_TFQMR (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#include "TpetraCore_ETIHelperMacros.h"
TPETRA_ETI_MANGLING_TYPEDEFS()

namespace Belos {
namespace Details {
namespace Tpetra {

void registerSolverFactory() {
  ::BelosTpetra::Impl::register_BiCGStab (false);
  ::BelosTpetra::Impl::register_BlockCG (false);
  ::BelosTpetra::Impl::register_BlockGmres (false);
  ::BelosTpetra::Impl::register_Cg (false);
  ::BelosTpetra::Impl::register_CgPipeline (false);
  ::BelosTpetra::Impl::register_CgSingleReduce (false);
  ::BelosTpetra::Impl::register_FixedPoint (false);
  ::BelosTpetra::Impl::register_GCRODR (false);
  ::BelosTpetra::Impl::register_Gmres (false);
  ::BelosTpetra::Impl::register_GmresPipeline (false);
  ::BelosTpetra::Impl::register_GmresPoly (false);
  ::BelosTpetra::Impl::register_GmresSingleReduce (false);
  ::BelosTpetra::Impl::register_GmresSstep (false);
  ::BelosTpetra::Impl::register_LSQR (false);
  ::BelosTpetra::Impl::register_Minres (false);
  ::BelosTpetra::Impl::register_PCPG (false);
  ::BelosTpetra::Impl::register_PseudoBlockCG (false);
  ::BelosTpetra::Impl::register_PseudoBlockGmres (false);
  ::BelosTpetra::Impl::register_PseudoBlockTFQMR (false);
  ::BelosTpetra::Impl::register_TFQMR (false);
}

} // namespace Tpetra
} // namespace Details
} // namespace Belos

