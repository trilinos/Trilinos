// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Belos_Details_Xpetra_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Xpetra

#include "BelosSolverFactory_Xpetra.hpp"

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

#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_TPETRA
# include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()
#endif

namespace Belos {
namespace Details {
namespace Xpetra {

void registerSolverFactory() {

// Note the formatting here was taken from the original Xpetra macros but then
// changed things so we wouldn't have any MueLu references. HAVE_MUELU_EPETRA
// becomes HAVE_XPETRA_EPETRA for example. That means we have quite a bit of
// duplication here but would like to test this and submit the proposed design
// before going further with decisions regarding how to eliminate the duplication.

#if   (defined(HAVE_XPETRA_EPETRA) &&  defined(EPETRA_HAVE_OMP) && (!defined(HAVE_XPETRA_TPETRA) || !defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with OpenMP node, but Tpetra is a) not enabled, or b) is not instantiated on OpenMP, or c) is not instantiated on OpenMP with <double,int,int>
  typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode EpetraNode;
#elif (defined(HAVE_XPETRA_EPETRA) && !defined(EPETRA_HAVE_OMP) && (!defined(HAVE_XPETRA_TPETRA) || !defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with Serial node, but Tpetra is a) not enabled, or b) is not instantiated on Serial, or c) is not instantiated on Serial with <double,int,int>
  typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode EpetraNode;
#endif


// Epetra = on, Tpetra = off
#if defined(HAVE_XPETRA_EPETRA) && !defined(HAVE_XPETRA_TPETRA)
  #define BELOS_XPETRA_CALL(INSTMACRO) INSTMACRO(double, int, int, EpetraNode)
#endif

// Epetra = on, Tpetra = on
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_TPETRA)
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
  #define BELOS_XPETRA_CALL(INSTMACRO) INSTMACRO(double, int, int, EpetraNode) TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)
# else
  #define BELOS_XPETRA_CALL(INSTMACRO) TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)
#endif
#endif

// Epetra = off, Tpetra = on
#if !defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_TPETRA)
  #define BELOS_XPETRA_CALL(INSTMACRO) TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)
#endif

  #define BELOS_LCL_CALL_FOR_MANAGER(manager,name,SC, LO, GO, NT)              \
    Impl::registerSolverSubclassForTypes<manager<                              \
              SC,                                                              \
              ::Xpetra::MultiVector<SC, LO, GO, NT>,                           \
              ::Belos::OperatorT<::Xpetra::MultiVector<SC, LO, GO, NT>>>,      \
              SC,                                                              \
              ::Xpetra::MultiVector<SC, LO, GO, NT>,                           \
              ::Belos::OperatorT<::Xpetra::MultiVector<SC, LO, GO, NT>>> (name);

  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(BiCGStabSolMgr, "BICGSTAB", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(BlockCGSolMgr, "BLOCK CG", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(BlockGmresSolMgr, "BLOCK GMRES", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(FixedPointSolMgr, "FIXED POINT", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(GCRODRSolMgr, "GCRODR", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(LSQRSolMgr, "LSQR", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(MinresSolMgr, "MINRES", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(PCPGSolMgr, "PCPG", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(PseudoBlockCGSolMgr, "PSEUDOBLOCK CG", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(PseudoBlockGmresSolMgr, "PSEUDOBLOCK GMRES", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(PseudoBlockTFQMRSolMgr, "PSEUDOBLOCK TFQMR", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(RCGSolMgr, "RCG", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )

  #undef LCL_CALL
  #define LCL_CALL( SC, LO, GO, NT ) BELOS_LCL_CALL_FOR_MANAGER(TFQMRSolMgr, "TFQMR", SC, LO, GO, NT)
  BELOS_XPETRA_CALL( LCL_CALL )
}

} // namespace Xpetra
} // namespace Details
} // namespace Belos

