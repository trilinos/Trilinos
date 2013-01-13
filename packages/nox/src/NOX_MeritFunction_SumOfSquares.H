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

#ifndef NOX_MERITFUNCTION_SUMOFSQUARES_H
#define NOX_MERITFUNCTION_SUMOFSQUARES_H

#include "NOX_Common.H"  // for std::ostream
#include "Teuchos_RCP.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_LineSearch_Utils_Slope.H"

namespace NOX {

namespace MeritFunction {

  //! Sum of squares merit function.
  /*!
    A basic merit function used in many nonlinear equation solvers:
    
    \f[
      f = \frac{1}{2} \| F(x) \| ^2
    \f]

    Where the norm is the 2-Norm using the NOX::Abstract::Vector's
    inner product.

    This is the default merit function used in nox.

    This merit function is taken from: J. E. Dennis Jr. and Robert B.
    Schnabel, "Numerical Methods for Unconstrained Optimization and
    Nonlinear Equations," Prentice Hall, 1983
  */
  class SumOfSquares : public virtual NOX::MeritFunction::Generic {

  public:			

    //! Constructor.
    SumOfSquares(const Teuchos::RCP<NOX::Utils>& u);
    
    //! Destructor.
    virtual ~SumOfSquares();

    //! Computes the merit function, \f$ f(x) = \frac{1}{2}\| F(x) \|^2 \f$. 
    virtual double computef(const NOX::Abstract::Group& grp) const;
    
    //! Computes the gradient, \f$ g = \nabla f(x) = J(x)^T F(x) \f$.
    virtual void computeGradient(const NOX::Abstract::Group& group,
				 NOX::Abstract::Vector& result) const;
    
    //! Computes the slope, \f$ s(x,d) = d^T \nabla f(x) = d^T J(x)^T F(x) \f$.
    /*! If the Jacobian is not computed in the \c grp object, then the
      slope can be approximated using directional derivatives.  More
      information can be found in the method computeSlopeWithoutJac.
     */
    virtual double computeSlope(const NOX::Abstract::Vector& dir,
				const NOX::Abstract::Group& grp) const;
    
    //! Computes the quadratic model, \f$ m(x,d) = f(x) + \nabla f(x)^T d + d^T \nabla^2 f(x) d \f$.
    /*!
      We approximate \f$ \nabla^2f(x) \approx J^TJ \f$: 

      \f[
        m(d) = f(x) + (J(x)^T F)^T d + \frac{1}{2} d^T B d
      \f]
    */
    virtual double computeQuadraticModel(const NOX::Abstract::Vector& dir,
				   const NOX::Abstract::Group& grp) const;
    
    //! Computes the vector  in the steepest descent direction that minimizes, the quadratic model.
    /*!
      Computes the vector \c result:
    \f[
      result = \frac{\nabla f^T \nabla f}{\nabla f^T B \nabla f} \nabla f = -\frac{(J^T F)^T (J^T F)}{(J J^T F)^T (J J^T F)} J^T F
    \f]
    */
    virtual void computeQuadraticMinimizer(const NOX::Abstract::Group& grp,
				     NOX::Abstract::Vector& result) const;

    virtual const std::string& name() const;

  private:
    
    //! Disallow default ctor.
    SumOfSquares() {};

    //! This is a variant of the computeSlope() method above optimized to work with out having to compute an explicit Jacobian.  
    /*!
      Calculates and returns 
      \f[
      \zeta = d^T \nabla f(x) = d^TJ^TF
      \f]

    Here \f$d\f$ represents the input parameter \c dir \f$\nabla
    f(x)\f$ is the gradient associated with the given group (for
    nonlinear solves this equates to \f$ J^TF \f$ where \f$ J \f$ is
    the Jacobian and \f$ F \f$ is the original nonlinear function).

    We can rewrite this equation as:

    \f[ d^TJ^TF = F^TJd \f]
    
    which allows us to use directional derivatives to estimate \f$ J^TF \f$:

    \f[ F^TJd = F^T \frac{F(x + \eta d) - F(x)}{\eta} \f]

    This may allow for faster computations of the slope if the
    Jacobian is expensive to evaluate.

    where \f$\eta\f$ is a scalar perturbation calculated by:

    \f[ \eta = \lambda * (\lambda + \frac{\| x\|}{\| d\|} ) \f]

    \f$ \lambda \f$ is a constant fixed at 1.0e-6.

    */
    virtual double 
    computeSlopeWithoutJacobian(const NOX::Abstract::Vector& dir,
				const NOX::Abstract::Group& grp) const;

    //! This is a variant of the computeSlope() method above that works when the Jacobian operator is available but doesn't support transpose
    /*!
      Calculates and returns 
      \f[
      \zeta = d^T \nabla f(x) = (Jd)^TF
      \f]

    Variables are as defined above.  This alternative form for the slope
    is provided to support Jacobian operators which do not
    support the transpose operation.
    */
    virtual double 
    computeSlopeWithoutJacobianTranspose(const NOX::Abstract::Vector& dir,
				         const NOX::Abstract::Group& grp) const;

  private:

    //!Printing utilities.
    Teuchos::RCP<Utils> utils;

    //! Temporary vector for computations.
    mutable Teuchos::RCP<NOX::Abstract::Vector> tmpVecPtr;

    //! Temporary vector for computations.
    /*! Only allocated if the method computeJacobianWithOutJac is called. */
    mutable Teuchos::RCP<NOX::Abstract::Group> tmpGrpPtr;

    //! Name of this function
    std::string meritFunctionName;

  };
} // namespace MeritFunction
} // namespace NOX

#endif
