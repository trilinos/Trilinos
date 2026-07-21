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
  Impl::registerSolverSubclassForTypes<manager<fST,fMVt,fOPt,fSDM>, fST, fMVt, fOPt, fSDM> (name);        \
  Impl::registerSolverSubclassForTypes<manager<dST,dMVt,dOPt,dSDM>, dST, dMVt, dOPt, dSDM> (name);        \
  Impl::registerSolverSubclassForTypes<manager<cST,cMVt,cOPt,cSDM>, cST, cMVt, cOPt, cSDM> (name);  \
  Impl::registerSolverSubclassForTypes<manager<cfST,cfMVt,cfOPt,cfSDM>, cfST, cfMVt, cfOPt, cfSDM> (name); \
  Impl::registerSolverSubclassForTypes<manager<fST,fMVk,fOPk,fKDV>, fST, fMVk, fOPk, fKDV> (name);        \
  Impl::registerSolverSubclassForTypes<manager<dST,dMVk,dOPk,dKDV>, dST, dMVk, dOPk, dKDV> (name);        \
  Impl::registerSolverSubclassForTypes<manager<cST,cMVk,cOPk,cKDV>, cST, cMVk, cOPk, cKDV> (name);  \
  Impl::registerSolverSubclassForTypes<manager<cfST,cfMVk,cfOPk,cfKDV>, cfST, cfMVk, cfOPk, cfKDV> (name);
#else // HAVE_TEUCHOS_COMPLEX
#define BELOS_DEFINE_REGISTER_SOLVER_MANAGER(manager,name)            \
  Impl::registerSolverSubclassForTypes<manager<fST,fMVt,fOPt,fSDM>, fST, fMVt, fOPt, fSDM> (name);  \
  Impl::registerSolverSubclassForTypes<manager<dST,dMVt,dOPt,dSDM>, dST, dMVt, dOPt, dSDM> (name);  \
  Impl::registerSolverSubclassForTypes<manager<fST,fMVk,fOPk,fKDV>, fST, fMVk, fOPk, fKDV> (name);  \
  Impl::registerSolverSubclassForTypes<manager<dST,dMVk,dOPk,dKDV>, dST, dMVk, dOPk, dKDV> (name);
#endif // HAVE_TEUCHOS_COMPLEX

void registerSolverFactory () {
  typedef double dST;
  typedef Teuchos::SerialDenseMatrix<int, dST> dSDM;
  typedef Kokkos::DualView<typename KokkosKernels::ArithTraits<dST>::val_type **, Kokkos::LayoutLeft> dKDV;
  typedef MultiVec<dST, dSDM> dMVt;
  typedef Operator<dST, dSDM> dOPt;
  typedef MultiVec<dST, dKDV> dMVk;
  typedef Operator<dST, dKDV> dOPk;

  typedef float fST;
  typedef Teuchos::SerialDenseMatrix<int, fST> fSDM;
  typedef Kokkos::DualView<typename KokkosKernels::ArithTraits<fST>::val_type **, Kokkos::LayoutLeft> fKDV;
  typedef MultiVec<fST, fSDM> fMVt;
  typedef Operator<fST, fSDM> fOPt;
  typedef MultiVec<fST, fKDV> fMVk;
  typedef Operator<fST, fKDV> fOPk;

#ifdef HAVE_TEUCHOS_COMPLEX
  typedef std::complex<double> cST;
  typedef Teuchos::SerialDenseMatrix<int, cST> cSDM;
  typedef Kokkos::DualView<typename KokkosKernels::ArithTraits<cST>::val_type **, Kokkos::LayoutLeft> cKDV;
  typedef MultiVec<cST, cSDM> cMVt;
  typedef Operator<cST, cSDM> cOPt;
  typedef MultiVec<cST, cKDV> cMVk;
  typedef Operator<cST, cKDV> cOPk;

  typedef std::complex<float> cfST;
  typedef Teuchos::SerialDenseMatrix<int, cfST> cfSDM;
  typedef Kokkos::DualView<typename KokkosKernels::ArithTraits<cfST>::val_type **, Kokkos::LayoutLeft> cfKDV;
  typedef MultiVec<cfST, cfSDM> cfMVt;
  typedef Operator<cfST, cfSDM> cfOPt;
  typedef MultiVec<cfST, cfKDV> cfMVk;
  typedef Operator<cfST, cfKDV> cfOPk;
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

