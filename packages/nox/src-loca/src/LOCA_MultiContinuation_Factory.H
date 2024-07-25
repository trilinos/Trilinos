// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_FACTORY_H
#define LOCA_MULTICONTINUATION_FACTORY_H

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
  namespace MultiContinuation {
    class AbstractStrategy;
    class AbstractGroup;
  }
  namespace MultiPredictor {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    //! Factory for creating continuation strategy objects
    /*!
     * The parameters passed to the create() through the \c stepperParams
     * argument method should specify the "Continuation Method" as described
     * below, as well as any additional parameters for the particular strategy.
     * <ul>
     * <li> "Continuation Method" - Name of the continuation method. Valid
     *      choices are
     *   <ul>
     *   <li> "Arc Length" (LOCA::MultiContinuation::ArcLengthGroup) [Default]
     *   <li> "Natural" (LOCA::MultiContinuation::NaturalGroup)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create continuation strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param stepperParams [in] Stepper parameters as described above
       * \param grp [in] Underlying group
       * \param pred [in] Predictor strategy
       * \param paramIDs [in] Indicies of continuation parameters
       */
      Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>
      create(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs);

      //! Return strategy name given by \c stepperParams
      const std::string& strategyName(Teuchos::ParameterList& stepperParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace MultiContinuation

} // Namespace LOCA

#endif // LOCA_MULTICONTINUATION_FACTORY_H
