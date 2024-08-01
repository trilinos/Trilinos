// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TPETRA_FACTORY_H
#define LOCA_TPETRA_FACTORY_H

#include "LOCA_Abstract_Factory.H"    // base class
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

// Forward declarations
namespace LOCA{
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Tpetra {

    //! Implementation of the LOCA::Abstract::Factory for Tpetra groups.
    class Factory : public LOCA::Abstract::Factory {

    public:

      //! Constructor
      Factory();

      //! Destructor
      virtual ~Factory();

      //! Initialize factory
      virtual void
      init(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      /*!
       * @name Strategy create methods
       */
      //@{

       //! Create bordered system solver strategy
      virtual bool
      createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy);

      //@}

    private:

      //! Private to prohibit copying
      Factory(const Factory& fac);

      //! Private to prohibit copying
      Factory& operator = (const Factory& fac);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace Tpetra

} // Namespace LOCA

#endif
