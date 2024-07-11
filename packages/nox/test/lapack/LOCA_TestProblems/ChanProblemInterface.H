// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef CHANPROBLEMINTERFACE_H
#define CHANPROBLEMINTERFACE_H

#include "Teuchos_RCP.hpp"

#include "LOCA_LAPACK_Interface.H" // base class
#include "NOX_LAPACK_Vector.H" // class data element
#include <fstream>

namespace LOCA {
  class GlobalData;
}

/*!
  \brief LOCA-LAPACK problem interface for the Chan problem.
*/
/*!
  %ChanProblemInterface implements the LOCA::LAPACK::Interface for a
  1D finite-difference discretization of the Chan problem:
  \f[
      \frac{d^2 T}{d x^2} + \alpha s \left(1 + \frac{x + 0.5 x^2}{1 + 0.01 x^2}\right) = 0
  \f]
   subject to the boundary conditions \f$T(0) = T(1) = \beta\f$.  The parameters
   are \f$\alpha\f$, \f$\beta\f$, \f$s\f$, and \f$n\f$, the size of the
   discretization.  The scaling factor \f$s\f$ is used to test continuation
   at different scales of \f$\alpha\f$.
*/
class ChanProblemInterface :
  public LOCA::LAPACK::Interface {

public:

  //! Constructor.  Also generates an initial guess
  ChanProblemInterface(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int N, double a, double b, double s);

  //! Constructor.  Also generates an initial guess
  ChanProblemInterface(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int N, double a, double b, double s, std::ofstream& file);

  //! Destructor
  ~ChanProblemInterface() {};

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

  //! Computes the initial guess
  void init();

  //! Computes the source term
  double source_term(double x);

  //! Computes the derivative of the source term
  double source_deriv(double x);

  //! Computes the scaled source parameter
  double source_param(double a, double s);

  //! Private to prohibit copying
  ChanProblemInterface(const ChanProblemInterface&);

  //! Private to prohibit copying
  ChanProblemInterface& operator=(const ChanProblemInterface&);

private:

  //! Global data
  Teuchos::RCP<LOCA::GlobalData> globalData;

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

  //! Stores pointer to output file stream
  std::ofstream* outputFilePtr;

};

#endif
