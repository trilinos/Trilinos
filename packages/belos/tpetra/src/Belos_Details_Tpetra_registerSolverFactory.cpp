//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER


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

