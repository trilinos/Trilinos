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

#ifndef NOX_PARAMETER_USERNORM_H
#define NOX_PARAMETER_USERNORM_H

#include "NOX_Parameter_Arbitrary.H" // base class

#include "NOX_Common.H"  // for ostream

// Forward declaration
namespace NOX {
  namespace Abstract {
    class Vector;
  }
}

namespace NOX {

namespace Parameter {

//! Derived from NOX::Parameter::Arbitrary class to support a user defined norm that can be passed in throught the parameter list.
/*!  

  This class allows the user to write their own norm and pass it into
  the NOX solvers.  Check the documentation for each class (line
  searches and Directions) to see whether it supports a UserNorm.
  Currently it is used by NOX::Direction::Newton for computing
  spcialized norms for the adjustable forcing term and in the
  NOX::LineSearch::Polynomial and NOX::LineSearch::MoreThuente2 line
  searches.
  
To create and use a user defined norm:

<ol>

<li> Create a norm that derives from
NOX::Parameter::UserNorm. For example, the norm \c Foo might be
defined as shown below.

\code
class Foo : public NOX::Parameter::UserNorm {
// Insert class definition here
}
\endcode

<li> Create the appropriate entries in the parameter list, as follows.

\code
Foo foo();
params.sublist("Line Search").sublist("Polynomial").setParameter("User Defined Norm", foo);
\endcode

</ol>

*/
class UserNorm : public virtual NOX::Parameter::Arbitrary {

 public:			

  //! Default Constructor.
  UserNorm(){};

  //! Destructor.
  ~UserNorm(){};

  //! Computes a user defined norm based on the input vector.  We assume the norm is of the type <a,b>^1/2.
  virtual double norm(const NOX::Abstract::Vector& input) const = 0;

  //@{ \name Auxiliary Methods.

  /*! In certain algorithms, the algorithm itself may depend on the type of norm calculted.  For example, in the TrustRegionBased solver, the dogleg computation is based on the l-2 norm.  If the user wishes to define their own norm, the trust region algorithm will fail to compute the correct direction because it calls dot() due to the l-2 norm assumption.  Therefore, we allow the user to supply auxiliary functions that may be required if the norm definition is changed.  
   */  
  //! Computes the dot product \f$ a^Tb \f$ consistent with the the norm calculation above.
  virtual double dot (const NOX::Abstract::Vector& a,
		      const NOX::Abstract::Vector& b) const = 0;
  
  //@}

};
} // namespace Parameter
} // namespace NOX

#endif
