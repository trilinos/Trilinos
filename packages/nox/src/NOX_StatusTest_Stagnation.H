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

#ifndef NOX_STATUS_STAGNATION_H
#define NOX_STATUS_STAGNATION_H

#include "NOX_StatusTest_Generic.H"	// base class

// Forward declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

namespace StatusTest {

//! Failure test based on the convergence rate between nonlinear iterations.
/*!
  This status test returns a Failed status if the convergence rate is worse 
  than a specified tolerance for \f$ n \f$ consecutive iterations.

  In many codes if the Jacobian is not accurate, the solver can stagnate.  
  In this case, further nonlinear iterations make negligible progress in 
  reducing the norm of F.  To prevent unnecessary churning, this status test 
  will identify when a code stagnates and ends the solve.

  This status test returns a Failed status if \f$ n \f$ nonlinear steps 
  have consecutively failed the following test:

  At nonlinear iteration \f$ k \f$, the convergence rate, \f$ \eta \f$ is 
  computed by:
  
  \f[ \eta = \frac{\| F_k \|}{\| F_{k-1} \|} \f]

  If \f$ \eta \f$ is greater than or equal to the specified tolerance, then 
  a counter is incremented.  If the counter hits the user specified number of 
  iterations, the status test returns a Failed status.  NOTE: For the counter 
  to increment, the steps have to fail CONSECUTIVELY. The counter is reset 
  every time a convergence rate passes the requested tolerance.

  Based on experience the following values are recomended:

  <ul>
  <li> For Newton solves:  n = 50, tolerance = 1.0<br>
  <li> For Newton solves with a line search: n = 15, tolerance = 0.99
  <\ul>
*/
class Stagnation : public Generic {

public:

  //! Constructor.
  /*! n is the number of consecutive nonlinear iterations with a 
    convergence rate failure for this test to return a Failed status. 
   */
  Stagnation(int n = 50, double tolerance = 1.0);

  //! Destructor.
  virtual ~Stagnation();

  virtual NOX::StatusTest::StatusType checkStatus(const NOX::Solver::Generic& 
						  problem);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

  //! Returns the used specified number of steps that can consecutively fail the tolerance test before the test returns a failed status.
  virtual int getMaxNumSteps() const;
  
  //! Returns the current number of steps that have consecutively failed the tolerance test.
  virtual int getCurrentNumSteps() const;

  //! Returns the user specified tolerance.
  virtual double getTolerance() const;

  //! Returns the current convergence rate.
  virtual double getConvRate() const;

private:

  //! User supplied value of n. 
  int maxSteps;

  //! Current number of consecutive nonlinear iterations that have failed the specified tolerance.
  int numSteps;

  //! The last nonlinear iteration index. 
  /*! This is used to prevent counting a step multiple times if by chance 
    the status test is called multiple times between iterations.
  */ 
  int lastIteration;

  //! User specified tolerance.
  double tolerance;

  //! Currently computed convergence rate.
  double convRate;
  
  //! %Status
  NOX::StatusTest::StatusType status;


};

} // namespace Status
} // namespace NOX

#endif
