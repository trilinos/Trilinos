// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENSOLVER_DGGEVSTRATEGY_H
#define LOCA_EIGENSOLVER_DGGEVSTRATEGY_H

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

    //! LAPACK DGGEV eigensolver strategy
    /*!
     * This class implements an eigensolver strategy using LAPACK routines
     * DGGEV and DGEEV
     */
    class DGGEVStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      DGGEVStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      virtual ~DGGEVStrategy();

      //! Compute eigenvalues/eigenvectors
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
          NOX::Abstract::Group& group,
          Teuchos::RCP< std::vector<double> >& evals_r,
          Teuchos::RCP< std::vector<double> >& evals_i,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      DGGEVStrategy(const DGGEVStrategy&);

      //! Private to prohibit copying
      DGGEVStrategy& operator = (const DGGEVStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parsed top-level parameter list
      Teuchos::RCP<LOCA::Parameter::SublistParser> topParams;

      //! Eigensolver parameter list
      Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! The number of requested eigenvalues
      int nev;

      //! Which eigenvalues are of interest
      std::string which;

    }; // Class DGGEVStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
