// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#ifndef NOX_STATUS_MAXITERS_H
#define NOX_STATUS_MAXITERS_H

#include "NOX_StatusTest_Generic.H"	// base class
#include "NOX_Utils.H"                  // class data element

namespace NOX {

namespace StatusTest {

//! Failure test based on the maximum number of nonlinear solver iterations.
/*!
  Let \f$k\f$ denote the current number of iterations (accessed via
  NOX::Solver::getNumIterations) and \f$k_{\max}\f$ denote the
  tolerance set in the constructor of this status test.
  This test returns NOX::StatusTest::Failed if
  \f$
  k \geq k_{\rm max}.
  \f$
  Otherwise, it returns NOX::StatusTest::Unconverged.

  \note If checkStatus is called with the type set to
  NOX::StatusTest::None, it then the status is set to to
  NOX::Status::Unevaluated and returned. (Also #niters is set to -1.)
*/
class MaxIters : public Generic {

public:

  //! Constructor. Specify the maximum number of nonlinear solver iterations, \f$k_{\max}\f$ ands optinally an error stream for printing errors.
  MaxIters(int maxIterations, const NOX::Utils* u = NULL);

  //! Destructor.
  virtual ~MaxIters();

  virtual NOX::StatusTest::StatusType 
  checkStatus(const NOX::Solver::Generic& problem, 
	      NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

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

  //! Current number of iterations (if known)
  int niters;

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Ostream used to print errors
  NOX::Utils utils;

};

} // namespace Status
} // namespace NOX

#endif
