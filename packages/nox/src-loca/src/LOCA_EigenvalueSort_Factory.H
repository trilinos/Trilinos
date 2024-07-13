// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENVALUESORT_FACTORY_H
#define LOCA_EIGENVALUESORT_FACTORY_H

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
  namespace EigenvalueSort {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace EigenvalueSort {

    //! Factory for creating %EigenvalueSort strategy objects
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the sorting method "Sorting Method"
     * as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Sorting Order" - Name of the sorting method. Valid choices are
     *   <ul>
     *   <li> "LM" (LOCA::EigenvalueSort::LargestMagnitude) [Default]
     *   <li> "LR" (LOCA::EigenvalueSort::LargestReal)
     *   <li> "LI" (LOCA::EigenvalueSort::LargestImaginary)
     *   <li> "SM" (LOCA::EigenvalueSort::SmallestMagnitude)
     *   <li> "SR" (LOCA::EigenvalueSort::SmallestReal)
     *   <li> "SI" (LOCA::EigenvalueSort::SmallestImaginary)
     *   <li> "CA" (LOCA::EigenvalueSort::LargestRealInverseCayley)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create sorting strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param eigenParams [in] %Eigensolver parameters as described above
       */
      Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>
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

  } // Namespace EigenvalueSort

} // Namespace LOCA

#endif
