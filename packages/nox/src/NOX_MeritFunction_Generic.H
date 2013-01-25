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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "NOX_Common.H"  // for std::ostream
#include <string>

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
  virtual const std::string& name() const = 0;

};
} // namespace MeritFunction
} // namespace NOX

#endif
