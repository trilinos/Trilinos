// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
