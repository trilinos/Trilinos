// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/// \file Belos_Details_Epetra_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Epetra

#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"

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
  Impl::registerSolverSubclassForTypes<GmresPolySolMgr<ST,MV,OP>, ST, MV, OP> ("HYBRID BLOCK GMRES");
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

