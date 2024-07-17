// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STATUS_MAXITERS_H
#define LOCA_STATUS_MAXITERS_H

#include "LOCA_StatusTest_Abstract.H"

#include "Teuchos_RCP.hpp"
#include "LOCA_GlobalData.H"

namespace LOCA {

namespace StatusTest {

  // TODO Update documentation.

//! Failure test based on the maximum number of continuation steps.
/*!
  Let \f$k\f$ denote the current number of iterations (accessed via \c
  LOCA::Stepper::getNumTotalSteps) and \f$k_{\max}\f$ denote the
  tolerance set in the constructor of this status test.
  This test returns LOCA::StatusTest::Failed if
  \f$
  k \geq k_{\rm max}.
  \f$
  Otherwise, it returns LOCA::StatusTest::NotFinished.

  \note If checkStatus is called with the type set to
  LOCA::StatusTest::None, it then the status is set to to
  LOCA::Status::Unevaluated and returned. (Also #niters is set to -1.)
*/
class MaxIters : public Abstract {

public:

  //! Constructor. Specify the maximum number of nonlinear solver iterations, \f$k_{\max}\f$ ands optinally an error stream for printing errors.
  MaxIters(int maxIterations,
           bool return_failed_on_max_steps = true,
           const Teuchos::RCP<const LOCA::GlobalData> globalData = Teuchos::null);

  //! Destructor.
  virtual ~MaxIters();

  virtual LOCA::StatusTest::StatusType
//  checkStatus(const LOCA::Stepper& stepper,
  checkStatus(const LOCA::Abstract::Iterator& stepper,
              LOCA::StatusTest::CheckType status);

  virtual LOCA::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  //! Returns the Maximum number of iterations set in the constructor.
  virtual int getMaxIters() const;

  /*!
    \brief Returns the current number of iterations taken by the solver.

    Returns -1 if the status of this test is NOX::StatusTest::Unevaluated.
   */
  virtual int getNumIters() const;

private:

  //! Maximum number of iterations
  int maxiters;

  //! Determines if LOCA::StatusTest::Failed or LOCA::StatusTest::Finished
  //! if returned when the \c maxiters has been reached.
  bool return_failed_on_max_steps_;

  //! Current number of iterations (if known)
  int niters;

  //! %Status
  LOCA::StatusTest::StatusType status;

  //! Ostream used to print errors
  Teuchos::RCP<const LOCA::GlobalData> globalDataPtr_;

};

} // namespace Status
} // namespace LOCA

#endif
