// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENSOLVER_ANASAZISTRATEGY_H
#define LOCA_EIGENSOLVER_ANASAZISTRATEGY_H

#include "LOCA_Eigensolver_AbstractStrategy.H" // base class

#include "Teuchos_ParameterList.hpp"           // class element
#include "AnasaziOutputManager.hpp"            // class element
#include "Anasazi_LOCA_Sort.H"                 // class element

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

    //! Anasazi eigensolver strategy
    /*!
     * This class implements an eigensolver strategy using the generic
     * Trilinos eigensolver package Anasazi.  In particular, this strategy
     * uses the Anasazi::BlockKrylovSchur solver.  Since Anasazi is a generic
     * solver, this strategy will work with any group implementation.  This
     * strategy references the following parameters passed through the
     * \c eigenParams argument to the constructor (this list is passed
     * directly to the Anasazi::BlockKrylovSchulSolMgr solver manager):
     * <ul>
     * <li> "Operator" -- [string] (default: "Jacobian inverse") Operator
     *      to compute eigenvalues of.
     * <li> "Block Size" -- [int] (default: 1) Block size
     * <li> "Num Blocks" -- [int] (default: 30) Maximum number of blocks
     *       (equals the maximum length of the Arnoldi factorization
     * <li> "Num Eigenvalues" -- [int] (default: 4) Number of requested
     *       eigenvalues
     * <li> "Convergence Tolerance" -- [double] (default: 1.0e-7) Tolerance
     *       for the converged eigenvalues
     * <li> "Step Size" -- [int] (default: 1) Checks convergence
     *       every so many steps
     * <li> "Maximum Restarts" -- [int] (default: 1) Number of restarts allowed
     * <li> "Symmetric" -- [bool] (default:  false) Is the operator symmetric
     * <li> "Verbosity" -- [Anasazi::MsgType]
     *       (default: Anasazi::Errors+Anasazi::Warnings+Anasazi::FinalSummary)
     *       Verbosity level
     * <li> "Sorting Order" -- [string" (default: "LM") Sorting order of
     *       printed eigenvalues
     * </ul>
     */
    class AnasaziStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.  Used to
       * obtain linear-solver parameters and to instantiate sorting strategy.
       * \param eigenParams [in] %Eigensolver parameters as described above.
       * solver.
       */
      AnasaziStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      virtual ~AnasaziStrategy();

      //! Compute eigenvalues/eigenvectors
      /*!
       * The implementation here the sets up and calls the Anasazi
       * BlockKrylovSchur solver for computing eigenvalues.
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
      AnasaziStrategy(const AnasaziStrategy&);

      //! Private to prohibit copying
      AnasaziStrategy& operator = (const AnasaziStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parsed parameter list
      Teuchos::RCP<LOCA::Parameter::SublistParser> topParams;

      //! Typedef to simplify templating
      typedef NOX::Abstract::MultiVector MV;

      //! Typedef to simplify templating
      typedef LOCA::AnasaziOperator::AbstractStrategy OP;

      //! Eigensolver parameter list
      Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! Linear solver parameter list
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! The block size
      int blksz;

      //! The number of requested eigenvalues
      int nev;

      //! Is problem symmetric
      bool isSymmetric;
    }; // Class AnasaziStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
