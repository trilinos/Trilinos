// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H
#define LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H

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
  namespace Hopf {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*!
       * \brief Factory for creating solver objects for solving Moore-Spence
       * Hopf equations
       */
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Solver Method"
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Salinger Bordering" (LOCA::Hopf::MooreSpence::SalingerBordering) [Default]
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
    Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>
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

  } // Namespace Hopf

} // Namespace LOCA

#endif // LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H
