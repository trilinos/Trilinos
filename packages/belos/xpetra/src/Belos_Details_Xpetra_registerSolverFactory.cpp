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

# include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()

namespace Belos {
namespace Details {
namespace Xpetra {

void registerSolverFactory() {

// Note the formatting here was taken from the original Xpetra macros but then
// changed things so we wouldn't have any MueLu references. HAVE_MUELU_EPETRA
// becomes HAVE_XPETRA_EPETRA for example. That means we have quite a bit of
// duplication here but would like to test this and submit the proposed design
// before going further with decisions regarding how to eliminate the duplication.

#if   (defined(HAVE_XPETRA_EPETRA) &&  defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with OpenMP node, but Tpetra is a) not enabled, or b) is not instantiated on OpenMP, or c) is not instantiated on OpenMP with <double,int,int>
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode EpetraNode;
#elif (defined(HAVE_XPETRA_EPETRA) && !defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with Serial node, but Tpetra is a) not enabled, or b) is not instantiated on Serial, or c) is not instantiated on Serial with <double,int,int>
  typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode;
#endif


// Epetra = on, Tpetra = off

// Epetra = on, Tpetra = on
#if defined(HAVE_XPETRA_EPETRA)
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
  #define BELOS_XPETRA_CALL(INSTMACRO) INSTMACRO(double, int, int, EpetraNode) TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)
# else
  #define BELOS_XPETRA_CALL(INSTMACRO) TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)
#endif
#endif

// Epetra = off, Tpetra = on
#if !defined(HAVE_XPETRA_EPETRA)
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

