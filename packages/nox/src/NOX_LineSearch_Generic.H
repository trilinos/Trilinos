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

#ifndef NOX_LINESEARCH_GENERIC_H
#define NOX_LINESEARCH_GENERIC_H

namespace NOX {
namespace Abstract {
class Vector;
class Group;
}
namespace Parameter {
class List;
}
}

namespace NOX {

//! %NOX linesearches
namespace LineSearch {

//! %Generic line search interface
/*!
  Every line search should respect the following Parameter:

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)
*/
class Generic {

public:

  //! Default constructor
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset parameters
  virtual bool reset(Parameter::List& params) = 0;
  
  //! Perform a line search from oldgrp along dirction dir. 
  //! Update newgrp = oldgrp + step * dir.
  virtual bool operator()(Abstract::Group& newgrp, double& step,
			  const Abstract::Group& oldgrp, 
			  const Abstract::Vector& dir) = 0;

};
} // namespace LineSearch
} // namespace NOX

#endif
