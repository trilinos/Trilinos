// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H
#define LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H

#include <vector>
#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H" // for ReturnType

namespace LOCA {

  //! Collection strategies to save eigenvector/value data
  namespace SaveEigenData {

    //! Abstract interface class strategies to save eigenvector/value data
    /*!
     * AbstractStrategy defines an abstract interface for saving eigenvectors
     * and eigenvalues that are computed at each continuation step.  This
     * is important because this data is often useful for restarting
     * continuations near bifurcation points and gives the user flexibility
     * in how this data is stored.
     *
     * The interface currently defines one pure virtual method,
     * save(), to save any eigenvectors or values as specified by the user.
     * Derived classes should implement this method for a particular strategy
     * to save this data, which is usually highly application code dependent.
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
     * is a parameter list of eigensolver parameters.  This list should also
     * specify which and how many eigenvectors/values to save as defined
     * by the strategy.
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

      //! Save eigenvalues/eigenvectors
      /*!
       * \param evals_r [out] Real eigenvalues
       * \param evals_i [out] Imaginary eigenvalues
       * \param evecs_r [out] Real eigenvectors
       * \param evecs_i [out] Imaginary eigenvectors
       * \returns ReturnType code indicating success or failure
       */
      virtual NOX::Abstract::Group::ReturnType
      save(Teuchos::RCP< std::vector<double> >& evals_r,
       Teuchos::RCP< std::vector<double> >& evals_i,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i) = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace SaveEigenData

} // Namespace LOCA

#endif // LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H
