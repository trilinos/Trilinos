// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LAPACK_INTERFACE_H
#define NOX_LAPACK_INTERFACE_H

// Forward declarations
namespace NOX {
  namespace LAPACK {
    class Vector;
    template <typename T> class Matrix;
  }
}

namespace NOX {
  namespace LAPACK {

    /*! \brief Virtual interface for NOX::BLAS::Group to provide a
      link to the external code for RHS and Jacobian fills. */
    class Interface {

    public:

      //! Constructor
      Interface() {};

      //! Destructor
      virtual ~Interface() {};

      //! Get initial guess
      virtual const NOX::LAPACK::Vector& getInitialGuess() = 0;

      //! Compute RHS given the specified input vector x.  Returns true if computation was successful.
      virtual bool computeF(NOX::LAPACK::Vector& rhs, const NOX::LAPACK::Vector &x) = 0;

      //! Compute Jacobian given the specified input vector x.  Returns true if computation was successful.
      virtual bool computeJacobian(NOX::LAPACK::Matrix<double>& J, const NOX::LAPACK::Vector & x) = 0;

    };
  } // namespace LAPACK
} // namespace NOX

#endif
