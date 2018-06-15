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

/// \file Belos_Details_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Belos

#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosSolverFactory.hpp"

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

#ifdef HAVE_TEUCHOS_COMPLEX
#define BELOS_DEFINE_REGISTER_SOLVER_MANAGER(manager,name)                           \
  Impl::registerSolverSubclassForTypes<manager<ST,MV,OP>, ST, MV, OP> (name);        \
  Impl::registerSolverSubclassForTypes<manager<cST,cMV,cOP>, cST, cMV, cOP> (name);
#else // HAVE_TEUCHOS_COMPLEX
#define BELOS_DEFINE_REGISTER_SOLVER_MANAGER(manager,name)            \
  Impl::registerSolverSubclassForTypes<manager<ST,MV,OP>, ST, MV, OP> (name);
#endif // HAVE_TEUCHOS_COMPLEX

void registerSolverFactory () {
  typedef double ST;
  typedef MultiVec<ST> MV;
  typedef Operator<ST> OP;

#ifdef HAVE_TEUCHOS_COMPLEX
  typedef std::complex<double> cST;
  typedef MultiVec<cST> cMV;
  typedef Operator<cST> cOP;
#endif // HAVE_TEUCHOS_COMPLEX

  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BiCGStabSolMgr, "BICGSTAB")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BlockCGSolMgr, "BLOCK CG")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BlockGmresSolMgr, "BLOCK GMRES")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(FixedPointSolMgr, "FIXED POINT")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(GCRODRSolMgr, "GCRODR")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(LSQRSolMgr, "LSQR")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(MinresSolMgr, "MINRES")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(PCPGSolMgr, "PCPG")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(PseudoBlockCGSolMgr, "PSEUDOBLOCK CG")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(PseudoBlockGmresSolMgr, "PSEUDOBLOCK GMRES")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(PseudoBlockTFQMRSolMgr, "PSEUDOBLOCK TFQMR")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(RCGSolMgr, "RCG")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(TFQMRSolMgr, "TFQMR")
}

} // namespace Details
} // namespace Belos

