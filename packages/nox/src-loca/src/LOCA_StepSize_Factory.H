// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STEPSIZE_FACTORY_H
#define LOCA_STEPSIZE_FACTORY_H

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
  namespace StepSize {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace StepSize {

    //! Factory for creating step size control strategy objects
    /*!
     * The parameters passed to the create() through the \c stepsizeParams
     * argument method should specify the
     * "Method" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Method" - Name of the step size control method. Valid choices are
     *   <ul>
     *   <li> "Constant" (LOCA::StepSize::Constant)
     *   <li> "Adaptive" LOCA::StepSize::Adaptive) [Default]
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create step size control strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param stepsizeParams [in] Step size parameters as described above
       */
      Teuchos::RCP<LOCA::StepSize::AbstractStrategy>
      create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams);

      //! Return strategy name given by \c stepsizeParams
      const std::string& strategyName(Teuchos::ParameterList& stepsizeParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace StepSize

} // Namespace LOCA

#endif // LOCA_STEPSIZE_FACTORY_H
