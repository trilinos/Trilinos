// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef CHANPROBLEMINTERFACE_H
#define CHANPROBLEMINTERFACE_H

#include "LOCA_LAPACK_Group.H"

/*!
  \brief LOCA-LAPACK problem interface for the Chan problem.
*/
/*!
  %ChanProblemInterface implements the LOCA::LAPACK::Interface for a 
  1D finite-difference discretization of the Chan problem:
  \f[
      \frac{d^2 T}{d x^2} + \alpha s \left(1 + \frac{x + 0.5 x^2}{1 + 0.01 x^2}\right) = 0
  \f]
   subject to the boundar conditions \f$T(0) = T(1) = \beta\f$.  The parameters
   are \f$\alpha\f$, \f$\beta\f$, \f$s\f$, and \f$n\f$, the size of the 
   discretization.  The scaling factor \f$s\f$ is used to test continuation
   at different scales of \f$\alpha\f$.
*/
class ChanProblemInterface : public LOCA::LAPACK::Interface {

public:
 
  //! Constructor.  Also generates an initial guess
  ChanProblemInterface(int N, double a, double b, double s);

  //! Destructor
  ~ChanProblemInterface() {};

  //! Returns initial guess
  const NOX::LAPACK::Vector& getInitialGuess();

  //! Evaluates residual
  bool
  computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector &x);
  
  //! Evalues jacobian
  bool
  computeJacobian(NOX::LAPACK::Matrix& J, const NOX::LAPACK::Vector & x);

  //! Sets parameters
  void setParams(const LOCA::ParameterVector& p);

private:

  //! Computes the source term
  double source_term(double x);

  //! Computes the derivative of the source term
  double source_deriv(double x);

  //! Computes the scaled source parameter
  double source_param(double a, double s);

private:

  //! Stores the initial guess
  NOX::LAPACK::Vector initialGuess;

  //! Stores \f$\alpha\f$
  double alpha;

  //! Stores \f$\beta\f$
  double beta;

  //! Stores the scale factor \f$s\f$
  double scale;

  //! Stores the size of the discretization
  int n;

};

#endif
