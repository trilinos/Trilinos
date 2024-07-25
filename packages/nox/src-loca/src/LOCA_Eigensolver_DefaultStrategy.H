// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENSOLVER_DEFAULTSTRATEGY_H
#define LOCA_EIGENSOLVER_DEFAULTSTRATEGY_H

#include "LOCA_Eigensolver_AbstractStrategy.H" // base class

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

  namespace Eigensolver {

    //! Default eigensolver strategy
    /*!
     * This class implements a default eigensolver strategy that does
     * not compute any eigenvalues.
     */
    class DefaultStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      DefaultStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      virtual ~DefaultStrategy();

      //! Compute eigenvalues/eigenvectors
      /*!
       * The implementation here does nothing and always returns
       * NOX::Abstract::Group::Ok.  Note that this implies the returned
       * ref-count pointers are null.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
          NOX::Abstract::Group& group,
          Teuchos::RCP< std::vector<double> >& evals_r,
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

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
