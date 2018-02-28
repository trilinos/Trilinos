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


/// \file Belos_Details_Tpetra_registerSolverFactory
/// \brief Implement Injection and Inversion (DII) for Tpetra

namespace Belos {
namespace Details {
namespace Tpetra {

void registerGCRODRSolMgr();
void registerPseudoBlockGmresSolMgr();
void registerPseudoBlockCGSolMgr();
void registerBlockGmresSolMgr();
void registerBlockCGSolMgr();
void registerFixedPointSolMgr();
void registerLSQRSolMgr();
void registerLSQRSolMgr();
void registerPCPGSolMgr();
void registerRCGSolMgr();
void registerBiCGStabSolMgr();
void registerMinresSolMgr();
void registerTFQMRSolMgr();
void registerPseudoBlockTFQMRSolMgr();

void registerSolverFactory() {
  registerGCRODRSolMgr();
  registerPseudoBlockGmresSolMgr();
  registerPseudoBlockCGSolMgr();
  registerBlockGmresSolMgr();
  registerBlockCGSolMgr();
  registerFixedPointSolMgr();
  registerLSQRSolMgr();
  registerPCPGSolMgr();
  registerRCGSolMgr();
  registerBiCGStabSolMgr();
  registerMinresSolMgr();
  registerTFQMRSolMgr();
  registerPseudoBlockTFQMRSolMgr();
}

} // namespace Tpetra
} // namespace Details
} // namespace Belos

namespace { // (anonymous)
  class Register_Belos_Details_Tpetra_SolverFactory {
  public:
    Register_Belos_Details_Tpetra_SolverFactory () {
      Belos::Details::Tpetra::registerSolverFactory();
    }
  };

  // ensure that the function actually gets called as premain
  // if it doesn't SolverFactoryParent constructor will manually call
  // this registerSolverFactory() at construction.

  // This would add premain registration.
  // However two tests were seg faulting on clang (not linux) for parallel:
  // Belos_Tpetra_MultipleSolves_MPI_4 & Ifpack2_AdditiveSchwarz_RILUK_MPI_4
  // Note there are three places where I commented this out.
  // TODO: Investigate that and decide if we want premain solve.

  // Register_Belos_Details_Tpetra_SolverFactory
  //   register_belos_details_tpetra_solverFactory;

} // namespace (anonymous)
