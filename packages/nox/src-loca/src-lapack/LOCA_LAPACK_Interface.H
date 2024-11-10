// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_LAPACK_INTERFACE_H
#define LOCA_LAPACK_INTERFACE_H

#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Interface.H"

// Forward declarations
namespace LOCA {
  class ParameterVector;
}

namespace LOCA {
  namespace LAPACK {

    /*!
      \brief Virtual interface for LOCA::LAPACK::Group to provide a
      link to the external code for RHS and Jacobian fills.
    */
    /*!
      This interface is derived from the NOX::LAPACK::Interface and
      and additionally provides a method for setting problem parameters.
    */
    class Interface : public NOX::LAPACK::Interface {

    public:

      //! Constructor
      Interface() {}

      //! Destructor
      virtual ~Interface() {}

      //! Set the problem parameters
      virtual void setParams(const ParameterVector& p) = 0;

      //! Call user's own print routine for vector-parameter pair
      virtual void printSolution(const NOX::LAPACK::Vector& /* x_ */,
                                 const double /* conParam */) {}

      /*!
       * \brief Compute shifted matrix alpha*J + beta*M where J is the
       * Jacobian matrix and M is the mass matrix.  Returns true if
       * computation was successful.
       */
      /*!
       * The default implementation assumes the mass matrix is the
       * identity matrix.
       */
      virtual bool computeShiftedMatrix(
                   double alpha, double beta,
                   const NOX::LAPACK::Vector& x,
                   NOX::LAPACK::Matrix<double>& A)
      {
    bool res = true;
    if (alpha != 0.0) {
      res = computeJacobian(A, x);
      A.scale(alpha);
    }
    else
      A.scale(0.0);
    if (beta != 0.0)
      for (int i=0; i<A.numRows(); i++)
        if (i < A.numCols())
          A(i,i) += beta;
    return res;
      }

      //! Projects solution to a few scalars for multiparameter continuation
      /*! Default implementation is the max norm
       */
      virtual void projectToDraw(const NOX::LAPACK::Vector& x,
                 double *px) const {
    px[0] = x.norm(NOX::Abstract::Vector::MaxNorm);
      }

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const { return 1; }

    };
  } // namespace LAPACK
} // namespace LOCA

#endif
