// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_ANASAZIOPERATOR_FACTORY_H
#define LOCA_ANASAZIOPERATOR_FACTORY_H

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
  namespace AnasaziOperator {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace AnasaziOperator {

    //! Factory for creating Anasazi operator strategy objects
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the
     * "Operator" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Operator" - Name of the Anasazi operator. Valid choices are
     *   <ul>
     *   <li> "Jacobian Inverse" (LOCA::AnasaziOperator::JacobianInverse)
     *        [Default]
     *   <li> "Shift-Invert" (LOCA::AnasaziOperator::ShiftInvert)
     *   <li> "Cayley" (LOCA::AnasaziOperator::Cayley)
     *   </ul>
     * </ul>
     *
     * There is also an %Epetra specific strategy that can
     * be instantiated by the LOCA::Epetra::Factory.  See
     * LOCA::Epetra::AnasaziOperator::Floquet.
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create Anasazi operator strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param eigenParams [in] %Eigensolver parameters as described above
       * \param solverParams [in] Linear solver parameters
       * \param grp [in] Group representing Jacobian/mass matrices
       */
      Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>
      create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
    const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
    const Teuchos::RCP<NOX::Abstract::Group>& grp);

      //! Return strategy name given by \c eigenParams
      const std::string& strategyName(Teuchos::ParameterList& eigenParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace AnasaziOperator

} // Namespace LOCA

#endif // LOCA_ANASAZIOPERATOR_FACTORY_H
