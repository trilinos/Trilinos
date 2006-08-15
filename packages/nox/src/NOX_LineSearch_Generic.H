// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_GENERIC_H
#define NOX_LINESEARCH_GENERIC_H

#include "Teuchos_RefCountPtr.hpp"

// Forward Declarations
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
  virtual bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
		     Teuchos::ParameterList& params) = 0;
  
  //! Perform a line search
  /*!
    Let
    
    - \f$x_{\rm new}\f$ denotes the new solution to be calculated (corresponding to \c grp)
    - \f$\lambda\f$ denotes the step to be calculated (\c step), 
    - \f$d\f$ denotes the search direction (\c dir).
    - \f$x_{\rm old}\f$ denotes the previous solution (i.e., the result of s.getPreviousSolutionGroup().getX())

    In the end, we should have computed \f$ \lambda  \f$ and updated \c grp so that
    \f[
    x_{\rm new} = x_{\rm old} + \lambda d.
    \f]

    Ideally, \f$ \|F(x_{\rm new})\| < \|F(x_{\rm old})\| \f$.
   */
  virtual bool compute(NOX::Abstract::Group& grp, 
		       double& step,
		       const NOX::Abstract::Vector& dir, 
		       const NOX::Solver::Generic& s) = 0;

};
} // namespace LineSearch
} // namespace NOX

#endif
