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

#ifndef NOX_LINESEARCH_UTILS_SLOPE_H
#define NOX_LINESEARCH_UTILS_SLOPE_H

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "Teuchos_RefCountPtr.hpp"

// Forward declaration
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

namespace LineSearch {

namespace Utils {

//! %Common line search utilites for computing the slope of a function.
/*!  
  This class provides routines for computing the slope of a give function.  There are two methods, one that uses a Jacobian and the other that estimates the action of the Jacobian by directional derivatives. 
*/

class Slope {

public:

  //! Default constructor
  Slope(const Teuchos::RefCountPtr<NOX::GlobalData>& gd);

  //! Destructor
  virtual ~Slope();

  //! Reset method.
  void reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd);

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient associated with the given group.
  */
  double computeSlope(const NOX::Abstract::Vector& dir, 
		      const NOX::Abstract::Group& grp);

  //! This is a variant of the computeSlope() method above optimized to work with out having to compute an explicit Jacobian.  
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x) = d^TJ^TF
    \f]

    Here \f$d\f$ represents the input parameter \c dir \f$\nabla
    f(x)\f$ is the gradient associated with the given group (for nonlinear solves this equates to \f$ J^TF \f$ where \f$ J \f$ is the Jacobian and \f$ F \f$ is the original nonlinear function).

    We can rewrite this equation as:

    \f[ d^TJ^TF = F^TJd \f]
    
    which allows us to use directional derivatives to estimate \f$ J^TF \f$:

    \f[ F^TJd = F^T \frac{F(x + \eta d) - F(x)}{\eta} \f]

    This may allow for faster computations of the slope if the Jacobian is expensive to evaluate.

where \f$\eta\f$ is a scalar perturbation calculated by:

\f[ \eta = \lambda * (\lambda + \frac{\| x\|}{\| d\|} ) \f]

\f$ \lambda \f$ is a constant fixed at 1.0e-6.

  */
  double computeSlopeWithOutJac(const NOX::Abstract::Vector& dir, 
				const NOX::Abstract::Group& grp);

private:
  //! Disallow default constructor.
  Slope() {};

private:

  //! Printing object
  NOX::Utils utils;

  //! A vector that *may* be used in computeSlope()
  Teuchos::RefCountPtr<NOX::Abstract::Vector> vecPtr;

  //! An extra group that will only be allocated if computeSlopeWithOutJac() is called.
  Teuchos::RefCountPtr<NOX::Abstract::Group> grpPtr;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
