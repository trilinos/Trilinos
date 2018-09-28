// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
