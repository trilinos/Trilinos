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

#ifdef WITH_PRERELEASE

#ifndef SOLVER_PICARD_H
#define SOLVER_PICARD_H

#include "NOX_Direction_Generic.H"	    // base class

// Forward declaration
namespace NOX {
class Utils;
}

namespace NOX {
namespace Direction {

//! Calculates a search direction using the simplest of methods, Picard.
/*!

Calculates the direction
\f[
s = -  F(x)
\f]

  To use this direction, specify that the "Method" is "Picard" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  Currently, no options may be specified in the "Picard" sublist
  of the "Direction" sublist of the solver parameters.

*/

class Picard : public Generic {

public:

  //! Constructor. 
  Picard(const NOX::Utils& u, Parameter::List& params);

  //! Destructor.
  virtual ~Picard();

  //! derived
  virtual bool reset(Parameter::List& p);

  //! derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
                       const Solver::Generic& solver);

private:
  
  //! Printing Utils
  const NOX::Utils& utils;

  //! Input parameters.
  Parameter::List* paramsPtr;	

  //! Number of nonlinear iterations.
  int niter;                    

};
} // namespace Direction
} // namespace NOX

#endif
#endif
