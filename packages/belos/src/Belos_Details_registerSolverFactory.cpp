// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "BelosGmresPolySolMgr.hpp"
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
  Impl::registerSolverSubclassForTypes<manager<fST,fMV,fOP>, fST, fMV, fOP> (name);        \
  Impl::registerSolverSubclassForTypes<manager<dST,dMV,dOP>, dST, dMV, dOP> (name);        \
  Impl::registerSolverSubclassForTypes<manager<cST,cMV,cOP>, cST, cMV, cOP> (name);  \
  Impl::registerSolverSubclassForTypes<manager<cfST,cfMV,cfOP>, cfST, cfMV, cfOP> (name);
#else // HAVE_TEUCHOS_COMPLEX
#define BELOS_DEFINE_REGISTER_SOLVER_MANAGER(manager,name)            \
  Impl::registerSolverSubclassForTypes<manager<fST,fMV,fOP>, fST, fMV, fOP> (name);  \
  Impl::registerSolverSubclassForTypes<manager<dST,dMV,dOP>, dST, dMV, dOP> (name);
#endif // HAVE_TEUCHOS_COMPLEX

void registerSolverFactory () {
  typedef double dST;
  typedef MultiVec<dST> dMV;
  typedef Operator<dST> dOP;

  typedef float fST;
  typedef MultiVec<fST> fMV;
  typedef Operator<fST> fOP;

#ifdef HAVE_TEUCHOS_COMPLEX
  typedef std::complex<double> cST;
  typedef MultiVec<cST> cMV;
  typedef Operator<cST> cOP;

  typedef std::complex<float> cfST;
  typedef MultiVec<cfST> cfMV;
  typedef Operator<cfST> cfOP;
#endif // HAVE_TEUCHOS_COMPLEX

  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BiCGStabSolMgr, "BICGSTAB")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BlockCGSolMgr, "BLOCK CG")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(BlockGmresSolMgr, "BLOCK GMRES")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(FixedPointSolMgr, "FIXED POINT")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(GCRODRSolMgr, "GCRODR")
  BELOS_DEFINE_REGISTER_SOLVER_MANAGER(GmresPolySolMgr, "HYBRID BLOCK GMRES")
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

