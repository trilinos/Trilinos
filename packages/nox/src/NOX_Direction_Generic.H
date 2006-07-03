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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_GENERIC_H
#define NOX_DIRECTION_GENERIC_H

#include "Teuchos_RefCountPtr.hpp"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Solver {
    class Generic;
    class LineSearchBased;
  }
}

namespace NOX {

//! Search direction strategies
/*!
  The Generic class spells out the interface. The Manager class can
  instantiate any other concrete implementation. The remaining classes
  are concrete implementations.
*/
namespace Direction {

//! %Generic direction interface
/*!
  %Generic interface for calculating a search direction, \f$d\f$, to
  be used in updating the iterate.
 */
class Generic {

public:

  //! Constructor. 
  /*! 
    Constructors of derived objects should look like reset().
   */
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset direction based on possibly new parameters.
  virtual bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
		     Teuchos::ParameterList& params) = 0;
  
  /*!
    \brief Compute the direction vector, \c dir, for a specific method given
    the current group, \c grp. 

    The \c grp is not const so that we can compute the F vector, the
    Jacobian matrix, the %Newton vector, and so on.

    Const access to the solver is used for getting additional
    information such as the past solution, the iteration number, and
    so on.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver) = 0;

  /*!

    \brief Same as compute(NOX::Abstract::Vector&,
    NOX::Abstract::Group&, const NOX::Solver::Generic&)

    Enables direct support for line search based solvers for the
    purpose of efficiency since the LineSearchBased object has a
    getStep() function that some directions require.

    If it is not redefined in the derived class, it will just call the
    compute with the NOX::Solver::Generic argument.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::LineSearchBased& solver);
};
} // namespace Direction
} // namespace NOX

#endif
