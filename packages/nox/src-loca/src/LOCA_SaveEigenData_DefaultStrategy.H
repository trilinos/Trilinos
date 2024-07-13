// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SAVEEIGENDATA_DEFAULTSTRATEGY_H
#define LOCA_SAVEEIGENDATA_DEFAULTSTRATEGY_H

#include "LOCA_SaveEigenData_AbstractStrategy.H" // base class

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace SaveEigenData {

    //! Default strategy for saving eigenvector/value data
    /*!
     * This class implements a default strategy for saving eigenvectors
     * and eigenvalues that does nothing and exists so the LOCA::Stepper
     * always has an object to pass eigen data to.
     */
    class DefaultStrategy : public LOCA::SaveEigenData::AbstractStrategy {

    public:

      //! Constructor
      DefaultStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      virtual ~DefaultStrategy();

      //! Save eigenvalues/eigenvectors
      /*!
       * The implementation here does nothing and always returns
       * NOX::Abstract::Group::Ok.
       */
      virtual NOX::Abstract::Group::ReturnType
      save(Teuchos::RCP< std::vector<double> >& evals_r,
       Teuchos::RCP< std::vector<double> >& evals_i,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      DefaultStrategy(const DefaultStrategy&);

      //! Private to prohibit copying
      DefaultStrategy& operator = (const DefaultStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class DefaultStrategy

  } // Namespace SaveEigenData

} // Namespace LOCA

#endif
