// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENSOLVER_FACTORY_H
#define LOCA_EIGENSOLVER_FACTORY_H

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
  namespace Eigensolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace Eigensolver {

    //! Factory for creating %Eigensolver strategy objects
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the
     * "Method" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Method" - Name of the eigensolver method. Valid choices are
     *   <ul>
     *   <li> "Default" (LOCA::Eigensolver::DefaultStrategy) [Default]
     *   <li> "Anasazi" (LOCA::Eigensolver::AnasaziStrategy)
     *   </ul>
     * </ul>
     *
     * There is also a %LAPACK specific strategy that can
     * be instantiated by the LOCA::LAPACK::Factory.  See
     * LOCA::Eigensolver::DGGEVStrategy.
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create eigensolver strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param eigenParams [in] %Eigensolver parameters as described above
       */
      Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>
      create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

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

  } // Namespace Eigensolver

} // Namespace LOCA

#endif // LOCA_EIGENSOLVER_FACTORY_H
