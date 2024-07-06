// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_GENERIC_HPP
#define BELOSSOLVERFACTORY_GENERIC_HPP

// The GenericSolverFactory is different from the others because it is
// manually included and registers all the solvers on construction. It exists
// for cases which are not covered by the other solver types.
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

template<class SC, class MV, class OP>
class GenericSolverFactory : public Impl::SolverFactoryParent<SC, MV, OP>
{
  public:
    GenericSolverFactory() {
      Impl::registerSolverSubclassForTypes<BiCGStabSolMgr<SC,MV,OP>, SC, MV, OP> ("BICGSTAB");
      Impl::registerSolverSubclassForTypes<BlockCGSolMgr<SC,MV,OP>, SC, MV, OP> ("BLOCK CG");
      Impl::registerSolverSubclassForTypes<BlockGmresSolMgr<SC,MV,OP>, SC, MV, OP> ("BLOCK GMRES");
      Impl::registerSolverSubclassForTypes<FixedPointSolMgr<SC,MV,OP>, SC, MV, OP> ("FIXED POINT");
      Impl::registerSolverSubclassForTypes<GCRODRSolMgr<SC,MV,OP>, SC, MV, OP> ("GCRODR");
      Impl::registerSolverSubclassForTypes<LSQRSolMgr<SC,MV,OP>, SC, MV, OP> ("LSQR");
      Impl::registerSolverSubclassForTypes<MinresSolMgr<SC,MV,OP>, SC, MV, OP> ("MINRES");
      Impl::registerSolverSubclassForTypes<PCPGSolMgr<SC,MV,OP>, SC, MV, OP> ("PCPG");
      Impl::registerSolverSubclassForTypes<PseudoBlockCGSolMgr<SC,MV,OP>, SC, MV, OP> ("PSEUDOBLOCK CG");
      Impl::registerSolverSubclassForTypes<PseudoBlockGmresSolMgr<SC,MV,OP>, SC, MV, OP> ("PSEUDOBLOCK GMRES");
      Impl::registerSolverSubclassForTypes<PseudoBlockTFQMRSolMgr<SC,MV,OP>, SC, MV, OP> ("PSEUDOBLOCK TFQMR");
      Impl::registerSolverSubclassForTypes<RCGSolMgr<SC,MV,OP>, SC, MV, OP> ("RCG");
      Impl::registerSolverSubclassForTypes<TFQMRSolMgr<SC,MV,OP>, SC, MV, OP> ("TFQMR");
    };
};

} // namespace Belos

#endif // BELOSSOLVERFACTORY_GENERIC_HPP
