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

#ifndef NOX_DIRECTION_GENERIC_H
#define NOX_DIRECTION_GENERIC_H

// Forward declarations
namespace NOX {
namespace Abstract {
class Vector;
class Group;
}
namespace Solver {
class Generic;
}
namespace Parameter {
class List;
}
}

namespace NOX {

//! Search direction strategies
namespace Direction {

//! %Generic direction interface
/*!
  %Generic interface to calculating search directions.
 */
class Generic {

public:

  //! Default base constructor. 
  /*! 
    Constructors of derived objects should look like reset.
   */
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset parameters
  virtual bool reset(const Parameter::List& params) = 0;
  
  /*!
    \brief Compute the direction vector for a specific method given
    the current group.

     Access to the solver is also allowed for getting additional
     status-type information like iteration number, etc.
  */
  virtual bool operator()(Abstract::Vector& dir, Abstract::Group& grp, 
			  const Solver::Generic& solver) = 0;

};
} // namespace Direction
} // namespace NOX

#endif
