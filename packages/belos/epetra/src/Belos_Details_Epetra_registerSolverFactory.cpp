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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER


/// \file Belos_Details_Epetra_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Epetra

#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"

#include "BelosBiCGStabSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosMinresSolMgr.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"

namespace Belos {
namespace Details {
namespace Epetra {

void registerSolverFactory() {
  typedef double ST;
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  Impl::registerSolverSubclassForTypes<BiCGStabSolMgr<ST,MV,OP>, ST, MV, OP> ("BICGSTAB");
  Impl::registerSolverSubclassForTypes<BlockCGSolMgr<ST,MV,OP>, ST, MV, OP> ("BLOCK CG");
  Impl::registerSolverSubclassForTypes<BlockGmresSolMgr<ST,MV,OP>, ST, MV, OP> ("BLOCK GMRES");
  Impl::registerSolverSubclassForTypes<FixedPointSolMgr<ST,MV,OP>, ST, MV, OP> ("FIXED POINT");
  Impl::registerSolverSubclassForTypes<GCRODRSolMgr<ST,MV,OP>, ST, MV, OP> ("GCRODR");
  Impl::registerSolverSubclassForTypes<LSQRSolMgr<ST,MV,OP>, ST, MV, OP> ("LSQR");
  Impl::registerSolverSubclassForTypes<MinresSolMgr<ST,MV,OP>, ST, MV, OP> ("MINRES");
  Impl::registerSolverSubclassForTypes<PCPGSolMgr<ST,MV,OP>, ST, MV, OP> ("PCPG");
  Impl::registerSolverSubclassForTypes<PseudoBlockCGSolMgr<ST,MV,OP>, ST, MV, OP> ("PSEUDOBLOCK CG");
  Impl::registerSolverSubclassForTypes<PseudoBlockGmresSolMgr<ST,MV,OP>, ST, MV, OP> ("PSEUDOBLOCK GMRES");
  Impl::registerSolverSubclassForTypes<PseudoBlockTFQMRSolMgr<ST,MV,OP>, ST, MV, OP> ("PSEUDOBLOCK TFQMR");
  Impl::registerSolverSubclassForTypes<RCGSolMgr<ST,MV,OP>, ST, MV, OP> ("RCG");
  Impl::registerSolverSubclassForTypes<TFQMRSolMgr<ST,MV,OP>, ST, MV, OP> ("TFQMR");
}

} // namespace Epetra
} // namespace Details
} // namespace Belos

