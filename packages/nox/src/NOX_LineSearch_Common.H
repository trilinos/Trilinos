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

#ifndef NOX_LINESEARCH_COMMON_H
#define NOX_LINESEARCH_COMMON_H

#include "NOX_LineSearch_Generic.H"
#include "NOX_Common.H"

// Forward declaration
namespace NOX {
class Utils;
}

namespace NOX {

namespace LineSearch {

//! %Common line search utilites
/*!  
  This is an abstract class that other line searches can derive
  from to get access to some common functionality.
*/
class Common : public Generic {

public:

  //! Default constructor
  Common(const NOX::Utils& u);

  //! Destructor
  virtual ~Common();

protected:

  //! Print out step information
  /*!
    NEED TO INSERT EXMAPLE OUTPUT HERE
  */
  void printStep(int n, double step, double oldf, double newf, const string s = "") const;

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient associated with the given group.
  */
  double computeSlope(const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp);

protected:

  //! Printing utilities
  const NOX::Utils& utils;

private:

  //! A vector that *may* be used in computeSlope()
  NOX::Abstract::Vector* vecPtr;

};
} // namespace LineSearch
} // namespace NOX

#endif
