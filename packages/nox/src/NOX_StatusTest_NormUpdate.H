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

#ifndef NOX_STATUSTEST_NORM_UPDATE_H
#define NOX_STATUSTEST_NORM_UPDATE_H

#include "NOX_StatusTest_Generic.H"	// base class
#include "NOX_Abstract_Vector.H"        // for NormType

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

//! Various convergence tests based on the norm of the change in the solution vector, \f$ x \f$, between outer iterations.
/*! 

  Return \c StatusTest::Converged if \f$\alpha < \beta\f$ where \f$\alpha\f$
  represents the norm of \f$ \Delta x \f$ and \f$\beta\f$ represents the tolerance.  We define:

\f[ \Delta x = x_k - x_{k-1} \f]

where \f$ x_k \f$ is the solution vector of the \f$ k \f$-th outer iterate.

  We define \f$\gamma\f$ to be the optional scale factor. We define it to be

  - \f$\gamma = \frac{1}{n}\f$ if ScaleType is \c Scaled, and

  - \f$\gamma = 1\f$ if ScaleType is \c Unscaled, and

  The norm of \f$ \Delta x \f$, \f$\alpha\f$, is defined as follows:

  - If the Abstract::Vector::NormType is \c TWO, then
    \f[ \alpha = \sqrt{ \gamma \sum_{i=1}^n \Delta x_i^2 } \f]

  - If the Abstract::Vector::NormType is \c ONE, then
    \f[ \alpha = \gamma \sum_{i=1}^n | \Delta x_i | \f]

  - If the Abstract::Vector::NormType is \c INF, then
    \f[ \alpha = \gamma \max_{i} | \Delta x_i |  \f]

  $\beta$ is a scalar set by the user in the constructor. 

  \f[  \beta  = \mbox{tolerance} \f]

  NOTE:  If no iterations have been performed, then there is no information for\f$ x_{k-1} \f$.  If the test is called while the number of iterations is zero (i.e. before an iteration takes place), the code returns the status test as Unconverged and sets the value to -1.0.  Once an iteration has taken place the Norm of the update is calculated as above.  

*/
class NormUpdate : public Generic {

public:

  //! Type that determines whether to scale the norm by the problem size.
  enum ScaleType {
    //! No norm scaling
    Unscaled, 
    //! Scale the norm by the length of the vector
    Scaled
  };

  //! Constructor for absolute norm. 
  /*! This constructor defaults to the \c Absolute tolerance type. */
  NormUpdate(double tolerance, NOX::Abstract::Vector::NormType ntype, ScaleType stype = Scaled);

  //! Constructor for absolute norm
  /*! This constructor defaults to the \c Absolute ToleranceType and \c TWO NormType. */
  NormUpdate(double tolerance, ScaleType stype = Scaled);

  //! Destructor.
  virtual ~NormUpdate();

  virtual NOX::StatusTest::StatusType checkStatus(const NOX::Solver::Generic& problem);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

  /* @name Accessor Functions
     Used to query current values of variables in the status test.
  */
  //@{

  //! Returns the value of the Update-norm computed in the last call to checkStatus.
  virtual double getNormUpdate() const;
   
  //! Returns the true tolerance.
  virtual double getTolerance() const;

  //@}

private:

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Vector containing the update for the current outer iteration 
  NOX::Abstract::Vector* updateVector;

  //! Type of norm to use
  NOX::Abstract::Vector::NormType normType;

  //! Scaling to use
  ScaleType scaleType;

  //! Tolerance value specified by user in constructor
  double tolerance;

  //! Norm of the update to be compared to trueTolerance
  double normUpdate;

};

} // namespace Status
} // namespace NOX

#endif
