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

#ifndef NOX_STATUSTEST_NORMF_H
#define NOX_STATUSTEST_NORMF_H

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

//! Various convergence tests based on the norm of the residual.
/*! 

  Return \c StatusTest::Converged if \f$\alpha < \beta\f$ where \f$\alpha\f$
  represents the norm of \f$F(x)\f$ and \f$\beta\f$ represents the
  tolerance.

  We define \f$\gamma\f$ to be the optional scale factor. We define it to be

  - \f$\gamma = \frac{1}{n}\f$ if ScaleType is \c Scaled, and

  - \f$\gamma = 1\f$ if ScaleType is \c Unscaled, and

  The norm of \f$F(x)\f$, \f$\alpha\f$, is defined as follows:

  - If the Abstract::Vector::NormType is \c TWO, then
    \f[ \alpha = \sqrt{ \gamma \sum_{i=1}^n F_i^2 } \f]

  - If the Abstract::Vector::NormType is \c ONE, then
    \f[ \alpha = \gamma \sum_{i=1}^n | F_i | \f]

  - If the Abstract::Vector::NormType is \c INF, then
    \f[ \alpha = \gamma \max_{i} | F_i |  \f]

  We set $\beta$ as follows. 

  - If an initial guess is provided, we use a relative tolerance
  defined by \f[ \beta = \frac{\mbox{tolerance}}{\alpha_0} \f] Here
  \f$\alpha_0\f$ is the \f$\alpha\f$ as defined above associated with
  the initial guess.

  - Otherwise, we use an absolte tolerance defined by
  \f[  \beta  = \mbox{tolerance} \f]

*/
class NormF : public Generic {

public:

  //! Type that determines whether to scale the norm by the problem size.
  enum ScaleType {
    //! No norm scaling
    Unscaled, 
    //! Scale the norm by the length of the vector
    Scaled
  };

  //! Type that determines whether the norm is absolute or relative to the intial guess
  enum ToleranceType {
    //! Relative to starting guess
    Relative, 
    //! Absolute
    Absolute
  };

  //! Constructor for absolute norm. 
  /*! This constructor defaults to the \c Absolute tolerance type. */
  NormF(double tolerance, NOX::Abstract::Vector::NormType ntype, ScaleType stype = Scaled);

  //! Constructor for absolute norm
  /*! This constructor defaults to the \c Absolute ToleranceType and \c TWO NormType. */
  NormF(double tolerance, ScaleType stype = Scaled);

  //! Constructor with initial guess (for relative norms)
  /*! This constructor defaults to the \c Relative tolerance type. */
  NormF(NOX::Abstract::Group& initialGuess, double tolerance, 
	NOX::Abstract::Vector::NormType ntype, ScaleType stype = Scaled);

  //! Constructor with initial guess (for relative norms)
  /*! This constructor defaults to the \c Relative ToleranceType and \c TWO NormType. */
  NormF(NOX::Abstract::Group& initialGuess, double tolerance, ScaleType stype = Scaled);

  //! Destructor.
  virtual ~NormF();

  virtual NOX::StatusTest::StatusType checkStatus(const NOX::Solver::Generic& problem);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

  /* @name Accessor Functions
     Used to query current values of variables in the status test.
  */
  //@{

  //! Returns the value of the F-norm computed in the last call to checkStatus.
  virtual double getNormF() const;
   
  //! Returns the true tolerance.
  virtual double getTrueTolerance() const;

  //! Returns the specified tolerance set in the constructor.
  virtual double getSpecifiedTolerance() const;

  //! Returns the initial tolerance.
  virtual double getInitialTolerance() const;

  //@}


private:

  /*! \brief Calculate the norm of F for the given group according to
    the scaling type, norm type, and tolerance type. */
  double computeNorm(const NOX::Abstract::Group& grp);

  /*! In the case of a relative norm calculation, initializes 
    \c trueTolerance based on the F-value at the initial guess.*/
  void relativeSetup(NOX::Abstract::Group& initialGuess);

private:

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Type of norm to use
  NOX::Abstract::Vector::NormType normType;

  //! Scaling to use
  ScaleType scaleType;

  //! Tolerance type (i.e., relative or absolute)
  ToleranceType toleranceType;

  //! Tolerance required for convergence.
  double specifiedTolerance;

  //! Initial tolerance
  double initialTolerance;

  //! True tolerance value, i.e., specifiedTolerance / initialTolerance
  double trueTolerance;

  //! Norm of F to be compared to trueTolerance
  double normF;

};

} // namespace Status
} // namespace NOX

#endif
