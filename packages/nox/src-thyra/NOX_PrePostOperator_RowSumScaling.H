// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H
#define NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_PrePostOperator.H"
#include "Thyra_VectorBase.hpp"

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

  /** \brief Updates an inverse row sum scaling vector at the beginning of a solve.

  */
  class RowSumScaling : public NOX::Abstract::PrePostOperator {

  public:

    enum ENOX_WhenToUpdateScaling {
      UpdateInvRowSumVectorAtBeginningOfSolve,
      UpdateInvRowSumVectorAtBeginningOfIteration
    };

    /** <tt>sv</tt> is the scaling vector */
    RowSumScaling(const Teuchos::RCP< ::Thyra::VectorBase<double> >& inv_row_sum_vec,
          ENOX_WhenToUpdateScaling s);

    void runPreIterate(const NOX::Solver::Generic& solver);

    void runPreSolve(const NOX::Solver::Generic& solver);

    Teuchos::RCP<const ::Thyra::VectorBase<double> >
    getInvRowSumScalingVector() const;

  private:

    void computeScaling(const NOX::Solver::Generic& solver);

    //! Inverse row sum vector
    Teuchos::RCP< ::Thyra::VectorBase<double> > inv_row_sum_vec_;

    ENOX_WhenToUpdateScaling when_to_update;
  };
}

#endif
