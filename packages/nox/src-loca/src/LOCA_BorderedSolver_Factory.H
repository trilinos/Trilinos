// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_FACTORY_H
#define LOCA_BORDEREDSOLVER_FACTORY_H

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
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace BorderedSolver {

    //! Factory for creating %BorderedSolver strategy objects
    /*!
     * The parameters passed to the create() through the \c solverParams
     * argument method should specify the "Bordered Solver Method"
     * as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Bordered Solver Method" - Name of the method. Valid choices are
     *   <ul>
     *   <li> "Bordering" (LOCA::BorderedSolver::Bordering) [Default]
     *   <li> "Nested" (LOCA::BorderedSolver::Nested)
     *   </ul>
     * </ul>
     *
     * There are also %Epetra and %LAPACK specific strategies that can
     * be instantiated by the LOCA::Epetra::Factory and LOCA::LAPACK::Factory.
     * See LOCA::BorderedSolver::LAPACKDirectSolve,
     * LOCA::BorderedSolver::EpetraHouseholder and
     * LOCA::BorderedSolver::Epetra::Augmented.
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create bordered system solver strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param solverParams [in] Solver parameters as described above
       */
      Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>
      create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Return strategy name given by \c solverParams
      const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_FACTORY_H
