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

#ifndef ROL_CONSTRAINTS_H
#define ROL_CONSTRAINTS_H

#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::Constraints
    \brief Provides the interface to evaluate constraint functions.
*/


namespace ROL {

template <class Real>
class Constraints {
private:
  bool activated_;

  std::vector<Teuchos::RCP<BoundConstraint<Real> > > ic_;
  std::vector<Teuchos::RCP<EqualityConstraint<Real> > >   ec_;

public:

  virtual ~Constraints() {}

  Constraints(void) : activated_(true) {}

  Constraints(std::vector<Teuchos::RCP<BoundConstraint<Real> > > & ic) : activated_(true), ic_(ic) {}

  Constraints(Teuchos::RCP<BoundConstraint<Real> > & ic) : activated_(true) {
    ic_.clear();
    ic_.push_back(ic);
  }

  Constraints(std::vector<Teuchos::RCP<EqualityConstraint<Real> > > & ec) : activated_(true), ec_(ec) {}

  Constraints(Teuchos::RCP<EqualityConstraint<Real> > & ec) : activated_(true) {
    ec_.clear();
    ec_.push_back(ec);
  }

  Constraints(std::vector<Teuchos::RCP<BoundConstraint<Real> > > & ic, std::vector<Teuchos::RCP<EqualityConstraint<Real> > > & ec) : activated_(true), ic_(ic), ec_(ec) {}

  Constraints(Teuchos::RCP<BoundConstraint<Real> > &ic, Teuchos::RCP<EqualityConstraint<Real> > & ec) : activated_(true) {
    ic_.clear();
    ic_.push_back(ic);
    ec_.clear();
    ec_.push_back(ec);
  }

  Constraints(Teuchos::RCP<BoundConstraint<Real> > &ic, std::vector<Teuchos::RCP<EqualityConstraint<Real> > > & ec) : activated_(true), ec_(ec) {
    ic_.clear();
    ic_.push_back(ic);
  }

  Constraints(std::vector<Teuchos::RCP<BoundConstraint<Real> > > &ic, Teuchos::RCP<EqualityConstraint<Real> > & ec) : activated_(true), ic_(ic) {
    ec_.clear();
    ec_.push_back(ec);
  }

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if x is changed,
                iter is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->update(x, flag, iter);
      }
    }
  }

  /** \brief Project optimization variables onto constraint set.
                x is the optimization variable
  */
  void project( Vector<Real> &x ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->project(x);
      }
    }
  }
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneUpperActive(v, g, x, eps);
      }
    }
  }
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneLowerActive(v, g, x, eps);
      }
    }
  }
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneUpperActive(v, x, eps);
      }
    }
  }
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneLowerActive(v, x, eps);
      }
    }
  }
  void setVectorToUpperBound( Vector<Real> &u ) {
    if (this->ic_[0]->isActivated()) {
      this->ic_[0]->setVectorToUpperBound(u);
    }
  }
  void setVectorToLowerBound( Vector<Real> &l ) {
    if (this->ic_[0]->isActivated()) {
      this->ic_[0]->setVectorToLowerBound(l);
    }
  }

  /** \brief Remove active set variables that are also in the binding set.
                v is the vector to be pruned 
                g is the gradient of the objective function at x
                x is the optimization variable
                eps is the active set tolerance
  */
  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneActive(v, g, x, eps);
      }
    }
  }

  /** \brief Remove active set variables.
                v is the vector to be pruned 
                x is the optimization variable
                eps is the active set tolerance
  */
  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        this->ic_[i]->pruneActive(v, x, eps);
      }
    }
  }

  /** \brief Check if the vector, v, is feasible
  */
  bool isFeasible( const Vector<Real> &v ) {
    bool iFeas = true;
    for (unsigned i=0; i<this->ic_.size(); i++) {
      if (this->ic_[i]->isActivated()) {
        iFeas = iFeas && (this->ic_[i]->isFeasible(v));
      }
    }
    return iFeas;
  }


  /** \brief Evaluate constraints:
             either for all activated equality constraints or
             for all activated inequality constraints.
  */
  void value(std::vector<Teuchos::RCP<Vector<Real> > > &c, const Vector<Real> &x, const EConstraint & ctype, std::vector<Real> & tol) {
    if (ctype == CONSTRAINT_EQUALITY) {
      for (unsigned i=0; i<size(this->ec_.size()); i++) {
        if (this->ec_[i]->isActivated()) {
          this->ec_[i]->value(*(c[i]), x, tol[i]);
        }
      }
    }
    else if (ctype == CONSTRAINT_INEQUALITY) {
      for (unsigned i=0; i<size(this->ic_.size()); i++) {
        if (this->ic_[i]->isActivated()) {
          this->ic_[i]->value(*(c[i]), x, tol[i]);
        }
      }
    }
  }


  /** \brief Apply constraint Jacobians or their adjoints:
             either for all activated equality constraints or
             for all activated inequality constraints.
  */
  void applyJacobian(std::vector<Teuchos::RCP<Vector<Real> > > &jv,
                     std::vector<Teuchos::RCP<const Vector<Real> > > &v,
                     const Vector<Real> &x,
                     const bool &adj,
                     const EConstraint &ctype,
                     std::vector<Real> &tol) {
    if (ctype == CONSTRAINT_EQUALITY) {
      for (unsigned i=0; i<size(this->ec_.size()); i++) {
        if (this->ec_[i]->isActivated()) {
          this->ec_[i]->applyJacobian(*(jv[i]), *(v[i]), x, adj, tol[i]);
        }
      }
    }
    else if (ctype == CONSTRAINT_INEQUALITY) {
      for (unsigned i=0; i<size(this->ic_.size()); i++) {
        if (this->ic_[i]->isActivated()) {
          this->ic_[i]->applyJacobian(*(jv[i]), *(v[i]), x, adj, tol[i]);
        }
      }
    }
  }


  /** \brief Compute the action of the operators W that are onto
             the null space (kernel) of the contraint Jacobians:
             either for all activated equality constraints or
             for all activated inequality constraints.
  */
  void maptoJacobianKernel(std::vector<Teuchos::RCP<Vector<Real> > > &wv,
                           std::vector<Teuchos::RCP<const Vector<Real> > > &v,
                           const Vector<Real> &x,
       	       	           const EConstraint &ctype,
                           std::vector<Real> &tol) {
    if (ctype == CONSTRAINT_EQUALITY) {
      for (unsigned i=0; i<size(this->ec_.size()); i++) {
        if (this->ec_[i]->isActivated()) {
          this->ec_[i]->maptoJacobianKernel(*(wv[i]), *(v[i]), x, tol[i]);
        }
      }
    }
    else if (ctype == CONSTRAINT_INEQUALITY) {
      for (unsigned i=0; i<size(this->ic_.size()); i++) {
        if (this->ic_[i]->isActivated()) {
          this->ic_[i]->maptoJacobianKernel(*(wv[i]), *(v[i]), x, tol[i]);
        }
      }
    }
  }


  /** \brief Apply constraint Hessians to (v,u), c''(x)(v,u) = (c''(x)u)v:
             either for all activated equality constraints or
             for all activated inequality constraints.
  */
  void applyHessian(std::vector<Teuchos::RCP<Vector<Real> > > &huv,
                    std::vector<Teuchos::RCP<const Vector<Real> > > &u,
                    std::vector<Teuchos::RCP<const Vector<Real> > > &v,
                    const Vector<Real> &x,
       	       	    const EConstraint &ctype,
                    std::vector<Real> &tol) {
    if (ctype == CONSTRAINT_EQUALITY) {
      for (unsigned i=0; i<size(this->ec_.size()); i++) {
        if (this->ec_[i]->isActivated()) {
          this->ec_[i]->applyHessian(*(huv[i]), *(u[i]), *(v[i]), x, tol[i]);
        }
      }
    }
    else if (ctype == CONSTRAINT_INEQUALITY) {
      for (unsigned i=0; i<size(this->ic_.size()); i++) {
        if (this->ic_[i]->isActivated()) {
          this->ic_[i]->applyHessian(*(huv[i]), *(u[i]), *(v[i]), x, tol[i]);
        }
      }
    }
  }


  /** \brief Remove the inactive set variables that are not in the binding set.
                v is the vector to be pruned 
                g is the gradient of the objective function at x
                x is the optimization variable
                eps is the active set tolerance
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) { 
    Teuchos::RCP<Vector<Real> > tmp = x.clone(); 
    tmp->set(v);
    this->pruneActive(*tmp,g,x,eps);
    v.axpy(-1.0,*tmp);
  }

  /** \brief Remove the inactive set variables.
                v is the vector to be pruned 
                x is the optimization variable
                eps is the active set tolerance
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) { 
    Teuchos::RCP<Vector<Real> > tmp = x.clone(); 
    tmp->set(v);
    this->pruneActive(*tmp,x,eps);
    v.axpy(-1.0,*tmp);
  }
 
  /** \brief Compute projected gradient.
  *             g is the gradient of the objective function at x
  *             x is the optimization variable
  */
  void computeProjectedGradient( Vector<Real> &g, const Vector<Real> &x ) {
    Teuchos::RCP<Vector<Real> > tmp = g.clone();
    tmp->set(g);
    this->pruneActive(g,*tmp,x);
  }
 
  /** \brief Compute projected step P(x+v)-x.
               v is the step vector
               x is the optimization variables
  */
  void computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) { 
    v.plus(x);
    this->project(v);
    v.axpy(-1.0,x);
  }

  /** \brief Turn on constraints 
  */
  void activate(void)    { this->activated_ = true;  }

  /** \brief Turn off constraints
  */
  void deactivate(void)  { this->activated_ = false; }

  /** \brief Check if constraints are on
  */
  bool isActivated(void) { return this->activated_;  }

}; // class Constraints

} // namespace ROL

#endif
