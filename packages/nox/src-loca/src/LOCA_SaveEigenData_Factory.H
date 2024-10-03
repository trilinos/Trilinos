// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SAVEEIGENDATA_FACTORY_H
#define LOCA_SAVEEIGENDATA_FACTORY_H

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
  namespace SaveEigenData {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace SaveEigenData {

    //! Factory for creating strategy objects to save eigenvectors/values
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the
     * "Method" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Save Eigen Data Method" - Name of the method. Valid choices are
     *   <ul>
     *   <li> "Default" (LOCA::SaveEigenData::DefaultStrategy) [Default]
     *   <li> "User-Defined" - User defined strategy
     *   </ul>
     * </ul>
     * User-defined strategies are defined by supplying the parameter
     * "User-Defined Save Eigen Data Name" which is the std::string name of the
     * strategy, and then a parameter with this name that is of the type
     * Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>.  This
     * allows the user to supply a user-defined strategy without providing
     * a factory to instantiate it.  By supplying the name of the parameter
     * storing the strategy, the user can provide multiple strategies in the
     * parameter list and select among them by setting
     * "User-Defined Save Eigen Data Name" to be the name of the strategy.
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create strategy to save eigenvectors/values
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param eigenParams [in] %Eigensolver parameters as described above
       */
      Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>
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

  } // Namespace SaveEigenData

} // Namespace LOCA

#endif // LOCA_SAVEEIGENDATA_FACTORY_H
