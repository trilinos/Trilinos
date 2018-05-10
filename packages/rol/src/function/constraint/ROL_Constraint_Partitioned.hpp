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

#ifndef ROL_CONSTRAINT_PARTITIONED_H
#define ROL_CONSTRAINT_PARTITIONED_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

/** @ingroup func_group
 *  \class ROL::Constraint_Partitioned
 *  \brief Has both inequality and equality constraints.
 *        Treat inequality constraint as equality with slack variable
 */

template<class Real>
class Constraint_Partitioned : public Constraint<Real> {
private:
  std::vector<ROL::Ptr<Constraint<Real> > > cvec_;
  std::vector<bool> isInequality_;      // Label whether cvec_[i] is inequality
  ROL::Ptr<Vector<Real> > scratch_; // Scratch vector for intermediate computation
  int  ncval_;                          // Number of constraint evaluations
  bool initialized_;                    // Is scratch vector initialized?

  Vector<Real>& getOpt( Vector<Real> &xs ) {
    try {
      return *dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  const Vector<Real>& getOpt( const Vector<Real> &xs ) {
    try {
      return *dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  Vector<Real>& getSlack( Vector<Real> &xs, const int ind ) {
    return *dynamic_cast<PartitionedVector<Real>&>(xs).get(ind);
  }

  const Vector<Real>& getSlack( const Vector<Real> &xs, const int ind ) {
    return *dynamic_cast<const PartitionedVector<Real>&>(xs).get(ind);
  }
  

public:
  Constraint_Partitioned(const std::vector<ROL::Ptr<Constraint<Real> > > &cvec,
                         bool isInequality = false)
   : cvec_(cvec),
     scratch_(ROL::nullPtr), ncval_(0), initialized_(false) {
    isInequality_.clear(); isInequality_.resize(cvec.size(),isInequality);
  }

  Constraint_Partitioned(const std::vector<ROL::Ptr<Constraint<Real> > > &cvec,
                         const std::vector<bool>                             &isInequality)
   : cvec_(cvec), isInequality_(isInequality),
     scratch_(ROL::nullPtr), ncval_(0), initialized_(false) {}

  int getNumberConstraintEvaluations(void) const {
    return ncval_;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const int ncon = static_cast<int>(cvec_.size());
    for (int i = 0; i < ncon; ++i) {
      cvec_[i]->update(getOpt(x),flag,iter);
    }
  }

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    PartitionedVector<Real> &cpv
      = dynamic_cast<PartitionedVector<Real>&>(c);

    const int ncon = static_cast<int>(cvec_.size());
    int cnt = 1;
    for (int i = 0; i < ncon; ++i) {
      cvec_[i]->value(*cpv.get(i), getOpt(x), tol);
      if (isInequality_[i]) {
        cpv.get(i)->axpy(static_cast<Real>(-1),getSlack(x,cnt));
        cnt++;
      }
    }
    ++ncval_;
  }

  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) {
    PartitionedVector<Real> &jvpv
      = dynamic_cast<PartitionedVector<Real>&>(jv);

    const int ncon = static_cast<int>(cvec_.size());
    int cnt = 1;
    for (int i = 0; i < ncon; ++i) {
      cvec_[i]->applyJacobian(*jvpv.get(i), getOpt(v), getOpt(x), tol);
      if (isInequality_[i]) {
        jvpv.get(i)->axpy(static_cast<Real>(-1),getSlack(v,cnt));
        cnt++;
      }
    }
  }

  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) {
    if (!initialized_) {
      scratch_ = getOpt(ajv).clone();
      initialized_ = true;
    }

    const PartitionedVector<Real> &vpv
      = dynamic_cast<const PartitionedVector<Real>&>(v);

    const int ncon = static_cast<int>(cvec_.size());
    int cnt = 1;
    getOpt(ajv).zero();
    for (int i = 0; i < ncon; ++i) {
      scratch_->zero();
      cvec_[i]->applyAdjointJacobian(*scratch_, *vpv.get(i), getOpt(x), tol);
      getOpt(ajv).plus(*scratch_);
      if (isInequality_[i]) {
        getSlack(ajv,cnt).set(*vpv.get(i));
        getSlack(ajv,cnt).scale(static_cast<Real>(-1));
        cnt++;
      }
    }
  }

  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) {
    if (!initialized_) {
      scratch_ = getOpt(ahuv).clone();
      initialized_ = true;
    }

    const PartitionedVector<Real> &upv
      = dynamic_cast<const PartitionedVector<Real>&>(u);

    const int ncon = static_cast<int>(cvec_.size());
    int cnt = 1;
    getOpt(ahuv).zero();
    for (int i = 0; i < ncon; ++i) {
      ROL::Ptr<const Vector<Real> > ui = upv.get(i);
      scratch_->zero();
      cvec_[i]->applyAdjointHessian(*scratch_, *ui, getOpt(v), getOpt(x), tol);
      getOpt(ahuv).plus(*scratch_);
      if (isInequality_[i]) {
        getSlack(ahuv,cnt).zero();
        cnt++;
      }
    }
  }

  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) {
    PartitionedVector<Real> &pvpv
      = dynamic_cast<PartitionedVector<Real>&>(pv);
    const PartitionedVector<Real> &vpv
      = dynamic_cast<const PartitionedVector<Real>&>(v);
    
    const int ncon = static_cast<int>(cvec_.size());
    for (int i = 0; i < ncon; ++i) {
      cvec_[i]->applyPreconditioner(*pvpv.get(i), *vpv.get(i), getOpt(x), getOpt(g), tol);
    }
  }

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) {
    Constraint<Real>::setParameter(param);
    const int ncon = static_cast<int>(cvec_.size());
    for (int i = 0; i < ncon; ++i) {
      cvec_[i]->setParameter(param);
    }
  }
}; // class Constraint_Partitioned

} // namespace ROL

#endif
