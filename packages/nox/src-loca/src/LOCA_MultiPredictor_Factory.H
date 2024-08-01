// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTIPREDICTOR_FACTORY_H
#define LOCA_MULTIPREDICTOR_FACTORY_H

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
  namespace MultiPredictor {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace MultiPredictor {

    //! Factory for creating %Predictor strategy objects
    /*!
     * The parameters passed to the create() through the \c predictorParams
     * argument method should specify the
     * "Method" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Method" - Name of the predictor method. Valid choices are
     *   <ul>
     *   <li> "Constant" (LOCA::MultiPredictor::Constant)
     *   <li> "Tangent" (LOCA::MultiPredictor::Tangent)
     *   <li> "Secant" (LOCA::MultiPredictor::Secant) [Default]
     *   <li> "Random" (LOCA::MultiPredictor::Random)
     *   <li> "Restart" (LOCA::MultiPredictor::Restart)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create predictor strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param predictorParams [in] %Predictor parameters as described above
       */
      Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
      create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& predictorParams);

      //! Return strategy name given by \c predictorParams
      const std::string& strategyName(Teuchos::ParameterList& predictorParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace MultiPredictor

} // Namespace LOCA

#endif // LOCA_MULTIPREDICTOR_FACTORY_H
