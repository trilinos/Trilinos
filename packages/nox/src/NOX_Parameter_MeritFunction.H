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

#ifndef NOX_PARAMETER_MERITFUNCTION_H
#define NOX_PARAMETER_MERITFUNCTION_H

#include "NOX_Parameter_Arbitrary.H" // base class

#include "NOX_Common.H"  // for ostream

// Forward declaration
namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

namespace Parameter {

//! Derived from NOX::Parameter::Arbitrary class to support a user defined merit function that can be passed in throught the parameter list.
/*!
  This class allows the user to define their own merit function for use in a line search.  Each line search type will specify in it's input parameter list if it supports this functionality.
*/
class MeritFunction : public virtual NOX::Parameter::Arbitrary {

 public:			

  //! Default Constructor.
  MeritFunction(){};

  //! Destructor.
  ~MeritFunction(){};

  //! Computes the merit function \f$ f \f$ to be used in a line search.
  virtual double computef(const NOX::Abstract::Group& grp) const = 0;

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient.
  */
  virtual double computeSlope(const NOX::Abstract::Vector& dir,
			      const NOX::Abstract::Group& grp) const = 0;

};
} // namespace Parameter
} // namespace NOX

#endif
