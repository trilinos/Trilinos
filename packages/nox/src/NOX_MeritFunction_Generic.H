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

#ifndef NOX_MERITFUNCTION_GENERIC_H
#define NOX_MERITFUNCTION_GENERIC_H

#include "NOX_Common.H"  // for ostream

// Forward declaration
namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

namespace MeritFunction {

//! Base class to support a user defined merit function that can be passed to line searches and directions through the parameter list.
/*!
  This class allows the user to define their own merit function for use in a line search.  Each line search type will specify in it's input parameter list if it supports this functionality.

To create and use a user defined merit function:

<ol>

<li> Create a merit function that derives from
NOX::Parameter::MeritFunction. For example, the merit function \c Foo might be
defined as shown below.

\code
class Foo : public NOX::Parameter::MeritFunction {
// Insert class definition here
}
\endcode

<li> Create the appropriate entries in the parameter list, as follows.

\code
Foo foo();
params.sublist("Solver Options").set("User Defined Merit Function", foo);
\endcode

</ol>

*/
class Generic {

 public:			

  //! Default Constructor.
  Generic(){};

  //! Destructor.
  virtual ~Generic(){};

  //! Computes the merit function, \f$ f(x) \f$.
  virtual double computef(const NOX::Abstract::Group& grp) const = 0;

  //! Computes the gradient of the merit function, \f$ \nabla f \f$, and returns the result in the \c result vector.  
  virtual void computeGradient(const NOX::Abstract::Group& group,
			       NOX::Abstract::Vector& result) const = 0;

  //! Computes the inner product of the given direction and the gradient associated with the merit function.  Returns the steepest descent direction in the \c result vector.
  /*!
    Calculates and returns \f$ \zeta \f$:
    \f[
      \zeta = \nabla f(x)^T d
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient of the merit function.
  */
  virtual double computeSlope(const NOX::Abstract::Vector& dir,
			      const NOX::Abstract::Group& grp) const = 0;
  
  //! Compute the quadratic model,\f$ m(d) \f$, for the given merit function.
  /*! Computes and returns \f$ m(d) \f$:
    \f[
      m(d) = f(x) + \nabla f(x)^T d + d^T \nabla^2 f(x) d + d^T \mathbf{B} d 
    \f]

    Here \f$d\f$ represents the input parameter \c dir. \f$ B \f$ is
    the Hessian of the merit function,\f$\nabla^2 f(x)\f$, but can be
    approximated with the restriction that it is a symmetric and has
    uniform boundedness in the iterate sequence (see J. Nocedal and
    S. J. Wright, "Numerical Optimization", Springer, 1999. Chapters 4
    and 6).
  */
  virtual double computeQuadraticModel(const NOX::Abstract::Vector& dir,
				   const NOX::Abstract::Group& grp) const = 0;

  //! Computes the vector in the steepest descent direction that minimizes the quadratic model.
  /*! The quadratic model is defined as: 
    \f[ 
      m(d) = f(x) + \nabla f(x)^T d + d^T \nabla^2 f(x) d + d^T \mathbf{B} d 
    \f]
    
    where \f$ B \f$ is ideally the Hessian of the merit
    function,\f$\nabla^2 f(x)\f$, but can be approximated with the
    restriction that it is a symmetric and has uniform boundedness in
    the iterate sequence (see J. Nocedal and S. J. Wright, "Numerical
    Optimization", Springer, 1999. Chapters 4 and 6).

    The \c result vector should be computed as:
    \f[
      result = -\frac{\nabla f^T \nabla f}{\nabla f^T B \nabla f} \nabla f
    \f]
  */
  virtual void computeQuadraticMinimizer(const NOX::Abstract::Group& grp,
				     NOX::Abstract::Vector& result) const = 0;

  //! Returns the name of the merit function.
  virtual const string& name() const = 0;

};
} // namespace MeritFunction
} // namespace NOX

#endif
