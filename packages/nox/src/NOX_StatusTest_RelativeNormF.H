// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUSTEST_RELATIVENORMF_H
#define NOX_STATUSTEST_RELATIVENORMF_H

#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Abstract_Vector.H"        // for NormType
#include "NOX_Utils.H"                  // for std::cerr

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

/** \brief Requires reduction in residual norm compared to initial norm.

    At iteration zero, this computes \f$||F_0||\f$ and is unconverged.

    Returns a converged status if \f$||F_k|| < ||F_0|| *
    \mbox{tolerance}\f$, where \f$k\f$ is the nonlinear iteration number.

*/
class RelativeNormF : public Generic {

public:

  RelativeNormF(double tolerance,
                bool scale_norms_by_vector_length = false,
                const NOX::Utils* u = NULL,
                NOX::Abstract::Vector::NormType ntype = NOX::Abstract::Vector::TwoNorm);

  // derived
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  // derived
  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

private:

  NOX::StatusTest::StatusType status;

  //! Tolerance required for convergence
  double tolerance;

  //! Initial Norm of F
  double normF_0;

  //! Current Norm of F
  double normF;

  //! If set to true, the norms are scaled by the size of the vector
  bool scale_norms_by_vector_length;

  //! Ostream used to print errors
  NOX::Utils utils;

  //! Type of norm to use
  NOX::Abstract::Vector::NormType norm_type;
};

} // namespace Status
} // namespace NOX

#endif
