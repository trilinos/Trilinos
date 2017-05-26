// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_FACTORY_HPP
#define ANASAZI_FACTORY_HPP

/*! \file AnasaziFactory.hpp
 *  \brief The Anasazi::Factory provides a factory to produce solver managers.
*/

#include "AnasaziConfigDefs.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziRTRSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziTraceMinDavidsonSolMgr.hpp"
#include "AnasaziTraceMinSolMgr.hpp"

namespace Anasazi {

class Factory {
public:

  template<class ScalarType, class MV, class OP>
  static
  Teuchos::RCP<SolverManager<ScalarType,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > problem,
           Teuchos::ParameterList &pl ) {
    using Teuchos::rcp;

    if (solverType == "BLOCK_DAVIDSON")
      return rcp(new BlockDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "BLOCK_KRYLOV_SCHUR")
      return rcp(new BlockKrylovSchurSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "LOBPCG")
      return rcp(new LOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "RTR")
      return rcp(new RTRSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "SIMPLE_LOBPCG")
      return rcp(new SimpleLOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else
      TEUCHOS_TEST_FOR_EXCEPTION(
                                 true, std::invalid_argument, "Anasazi::Factory::create: "
                                 "Invalid solver type \"" << solverType << "\".");
  }

  // Some of the solver managers cannot be instantiated on std::complex.
  template<class MV, class OP>
  static
  Teuchos::RCP<SolverManager<double,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<Eigenproblem<double,MV,OP> > problem,
           Teuchos::ParameterList &pl ) {
    using Teuchos::rcp;
    using ScalarType = double;

    if (solverType == "BLOCK_DAVIDSON")
      return rcp(new BlockDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "BLOCK_KRYLOV_SCHUR")
      return rcp(new BlockKrylovSchurSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "GENERALIZED_DAVIDSON")
      return rcp(new GeneralizedDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "LOBPCG")
      return rcp(new LOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "RTR")
      return rcp(new RTRSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "SIMPLE_LOBPCG")
      return rcp(new SimpleLOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "TRACE_MIN")
      return rcp(new Experimental::TraceMinSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (solverType == "TRACE_MIN_DAVIDSON")
      return rcp(new Experimental::TraceMinDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else
      TEUCHOS_TEST_FOR_EXCEPTION(
                                 true, std::invalid_argument, "Anasazi::Factory::create: "
                                 "Invalid solverType type \"" << solverType << "\".");
  }

  template<class ScalarType, class MV, class OP>
  static
  Teuchos::RCP<SolverManager<ScalarType,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<BasicEigenproblem<ScalarType,MV,OP> > &problem,
           Teuchos::ParameterList &pl ) {
    auto eproblem = Teuchos::rcp_static_cast<Eigenproblem<ScalarType,MV,OP>>(problem);
    return create(solverType, eproblem, pl);
  }
};


} // end Anasazi namespace

#endif /* ANASAZI_FACTORY_HPP */

