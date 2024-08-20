// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H
#define LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H

#include <vector>
#include "Teuchos_RCP.hpp"
#include "NOX_Abstract_Group.H" // for ReturnType

namespace LOCA {

  //! Collection of eigensolver strategies
  namespace Eigensolver {

    //! Abstract interface class for Eigensolver strategies
    /*!
     * AbstractStrategy defines an abstract interface for eigensolver
     * strategies. It is used by LOCA::Stepper to compute eigenvalues of
     * the steady-state solution after each continuation step.
     *
     * The interface currently defines one pure virtual method,
     * computeEigenvalues(), to compute the eigenvalues.  Derived classes
     * should implement this method for a particular eigensolver strategy.
     * Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c eigenParams
     * is a parameter list of eigensolver parameters.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Compute eigenvalues/eigenvectors in group \c group
      /*!
       * \param group [in] NOX Group to compute eigenvalues of
       * \param evals_r [out] Real eigenvalues
       * \param evals_i [out] Imaginary eigenvalues
       * \param evecs_r [out] Real eigenvectors
       * \param evecs_i [out] Imaginary eigenvectors
       * \returns ReturnType code indicating success or failure
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
          NOX::Abstract::Group& group,
          Teuchos::RCP< std::vector<double> >& evals_r,
          Teuchos::RCP< std::vector<double> >& evals_i,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i) = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif // LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H
