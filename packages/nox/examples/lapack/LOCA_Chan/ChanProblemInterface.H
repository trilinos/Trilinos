// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef CHANPROBLEMINTERFACE_H
#define CHANPROBLEMINTERFACE_H

#include "Teuchos_RefCountPtr.hpp"

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
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		  int N, double a, double b, double s);
 
  //! Constructor.  Also generates an initial guess
  ChanProblemInterface(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    int N, double a, double b, double s, ofstream& file);

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
  Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

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
  ofstream* outputFilePtr;

};

#endif
