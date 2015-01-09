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

#ifndef PITCHFORKPROBLEMINTERFACE_H
#define PITCHFORKPROBLEMINTERFACE_H

#include "LOCA_LAPACK_Group.H"

/*!
  \brief LOCA-LAPACK problem interface for a 1D finite difference discretization of a simple PDE exhibiting a pitchfork bifurcation.
*/
/*!
  %PitchforkProblemInterface implements the LOCA::LAPACK::Interface for a
  1D finite-difference discretization of the PDE:
  \f[
      \frac{d^2 T}{d x^2} + \lambda T - \alpha T^2 + \beta T^3 = 0
  \f]
   subject to the boundar conditions \f$T(-1) = T(1) = \beta\f$.  The
   parameters are \f$\alpha\f$, \f$\beta\f$, \f$\lambda\f$, and \f$n\f$, the
   size of the discretization.
*/
class PitchforkProblemInterface : public LOCA::LAPACK::Interface {

public:

  //! Constructor.  Also generates an initial guess
  PitchforkProblemInterface(int N, double a, double b, double l);

  //! Destructor
  ~PitchforkProblemInterface() {};

  //! Returns initial guess
  const NOX::LAPACK::Vector& getInitialGuess();

  //! Evaluates residual
  bool
  computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector &x);

  //! Evalues jacobian
  bool
  computeJacobian(NOX::LAPACK::Matrix<double>& J,
          const NOX::LAPACK::Vector & x);

  //! Sets parameters
  void setParams(const LOCA::ParameterVector& p);

  //! Prints solution after successful step
  void printSolution(const NOX::LAPACK::Vector &x, const double conParam);


private:

  //! Computes the source term
  double source_term(double x);

  //! Computes the derivative of the source term
  double source_deriv(double x);

private:

  //! Stores the initial guess
  NOX::LAPACK::Vector initialGuess;

  //! Stores \f$\alpha\f$
  double alpha;

  //! Stores \f$\beta\f$
  double beta;

  //! Stores \f$\lambda\f$
  double lambda;

  //! Stores the size of the discretization
  int n;

};

#endif
