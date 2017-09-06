// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINTFROMOBJECTIVE_H
#define ROL_CONSTRAINTFROMOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::ConstraintFromObjective
    \brief Creates a constraint from an objective function and a target value

    Example:  Suppose we have an objective function f(x) and we wish to impose
              a condition f(x) >= target, then this class creates the scalar
              constraint c(x) = f(x)-target 
*/


namespace ROL {

template<class Real> 
class ConstraintFromObjective : public Constraint<Real> {

    template <typename T> using RCP = Teuchos::RCP<T>;

  using V    = ROL::Vector<Real>;
  using Fill = ROL::Elementwise::Fill<Real>;
  using Sum  = ROL::Elementwise::ReductionSum<Real>;

private:

  RCP<ROL::Objective<Real>>  obj_;
  Real                       target_;
  Sum                        sum_;

public:

  ConstraintFromObjective( const RCP<Objective<Real>>& obj, Real target ) :
    obj_(obj), target_(target) {}

  /** \brief Update constraint
  */
  void update( const V& x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
  }
  
  /** \brief Evaluate constraint c(x) = f(x)-target
  */
  void value( V& c, const V& x, Real& tol ) {
    Real f = obj->value(x,tol);
    c.applyUnary(Fill(f-target_));
  }

  /** \brief Apply the constraint Jacobian at \f$x\f$, \f$c'(x) \in L(\mathcal{X}, \mathcal{C})\f$,
             to vector \f$v\f$.
      \f[ 
          c'(x)v = \lim_{t\rightarrow 0} \frac{d}{dt} f(x+tv) = \langle \nabla f(x),v\rangle 
      \f]
  */

  void applyJacobian( V& jv, const V& v, const V& x, Real& tol ) {
    jv.applyUnary(Fill(obj->dirDeriv(x,v,tol)));
  }

  /** \brief Apply the adjoint of the the constraint Jacobian at 
             \f$x\f$, \f$c'(x)^\ast \in L(\mathcal{C}^\ast, \mathcal{X}^\ast)\f$,
             to vector \f$v\f$.

      \f[ c'(x)^\ast v = v\nabla f(x) \f]
  */
  void applyAdjointJacobian( V& ajv, const V& v, const V& x, Real& tol ) {
    obj->gradient(ajv,x,tol);
    ajv.scale(v.reduce( sum_ ) );
  }

  /** \brief Apply the derivative of the adjoint of the constraint Jacobian at \f$x\f$
             to vector \f$u\f$ in direction \f$v\f$,
             according to \f$ v \mapsto c''(x)(v,\cdot)^\ast u \f$.

     \f[ c"(x)(v,\cdot)^\ast = u \frac{d}{dt} \nabla f(x+tv)\big|_{t=0}
  */
  void applyAdjointHessian( V& ahuv, const V& u, const V& v, const V& x, tol ) {
    obj->hessVec(ahuv,v,x,tol);
    ahuv.scale(u.reduce( sum_ ) );
  }
};

}; // ConstraintFromObjective

} // namespace ROL

#endif // ROL_CONSTRAINTFROMOBJECTIVE_H
