// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <algorithm>
#include <string>

/*! \class Anasazi::Factory
  \brief This provides a factory to build Anasazi solvers using parameter lists.
*/


namespace Anasazi {

class Factory {
public:

  /** \brief Create an instance of Anasazi::SolverManager given the string
   * name of the solver type.
   *
   * \param solverType [in] Name of solver type to be created.
   * \param problem [in] Anasazi eigenproblem used to define the solver
   *
   * Throw an exception if the solver with that input name does not exist.
   * Otherwise, return a newly created solver object.
   *
   * \note Some of the solver managers cannot be used with std::complex.
   */
  template<class ScalarType, class MV, class OP>
  static
  Teuchos::RCP<SolverManager<ScalarType,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > problem,
           Teuchos::ParameterList &pl ) {
    using Teuchos::rcp;

    std::string type = solverType;
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    if (type == "block_davidson" || type == "block davidson")
      return rcp(new BlockDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "block_krylov_schur" || type == "block krylov schur")
      return rcp(new BlockKrylovSchurSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "lobpcg")
      return rcp(new LOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "rtr")
      return rcp(new RTRSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "simple_lobpcg" || type == "simple lobpcg")
      return rcp(new SimpleLOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Anasazi::Factory::create: Invalid solver type \"" << solverType << "\".");

    return Teuchos::null;
  }

  template<class MV, class OP>
  static
  Teuchos::RCP<SolverManager<double,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<Eigenproblem<double,MV,OP> > problem,
           Teuchos::ParameterList &pl ) {
    using Teuchos::rcp;
    using ScalarType = double;

    std::string type = solverType;
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    if (type == "block_davidson" || type == "block davidson")
      return rcp(new BlockDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "block_krylov_schur" || type == "block krylov schur")
      return rcp(new BlockKrylovSchurSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "generalized_davidson" || type == "generalized davidson")
      return rcp(new GeneralizedDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "lobpcg")
      return rcp(new LOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "rtr")
      return rcp(new RTRSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "simple_lobpcg" || type == "simple lobpcg")
      return rcp(new SimpleLOBPCGSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "TRACE_MIN" || type == "trace min")
      return rcp(new Experimental::TraceMinSolMgr<ScalarType,MV,OP>(problem, pl));

    else if (type == "trace_min_davidson" || type == "trace min davidson")
      return rcp(new Experimental::TraceMinDavidsonSolMgr<ScalarType,MV,OP>(problem, pl));

    else
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
        "Anasazi::Factory::create: Invalid solverType type \"" << solverType << "\".");

    return Teuchos::null;
  }

  //! Specialize create for BasicEigenproblem type.
  template<class ScalarType, class MV, class OP>
  static
  Teuchos::RCP<SolverManager<ScalarType,MV,OP> >
  create ( const std::string& solverType,
           const Teuchos::RCP<BasicEigenproblem<ScalarType,MV,OP> > &problem,
           Teuchos::ParameterList &pl ) {
    Teuchos::RCP<Eigenproblem<ScalarType,MV,OP>> eproblem =
            Teuchos::rcp_static_cast<Eigenproblem<ScalarType,MV,OP>>(problem);
    return create(solverType, eproblem, pl);
  }
};


} // end Anasazi namespace

#endif /* ANASAZI_FACTORY_HPP */

