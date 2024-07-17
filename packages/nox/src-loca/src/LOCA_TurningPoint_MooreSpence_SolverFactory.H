// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H
#define LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace TurningPoint {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
}

namespace LOCA {

  namespace TurningPoint {

    namespace MooreSpence {

      /*!
       * \brief Factory for creating solver objects for solving Moore-Spence
       * turning point equations
       */
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Solver Method"
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Salinger Bordering" (LOCA::TurningPoint::MooreSpence::SalingerBordering) [Default]
       *   <li> "Phipps Bordering" (LOCA::TurningPoint::MooreSpence::PhippsBordering)
       *   </ul>
       * </ul>
       */
      class SolverFactory {

      public:

    //! Constructor
    SolverFactory(
          const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! Destructor
    virtual ~SolverFactory();

    //! Create solver strategy
    /*!
     * \param topParams [in] Parsed top-level parameter list.
     * \param solverParams [in] Solver parameters as described above
     */
    Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>
    create(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //! Return strategy name given by \c solverParams
    const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

      private:

    //! Private to prohibit copying
    SolverFactory(const SolverFactory&);

    //! Private to prohibit copying
    SolverFactory& operator = (const SolverFactory&);

      protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

      }; // Class SolverFactory

    } // Namespace MooreSpence

  } // Namespace TurningPoint

} // Namespace LOCA

#endif // LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H
