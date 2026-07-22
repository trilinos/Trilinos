// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SIMULATED_BOUND_CONSTRAINT_H
#define ROL_SIMULATED_BOUND_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_SimulatedVector.hpp"

/** @ingroup func_group
    \class ROL::SimulatedBoundConstraint
    \brief A BoundConstraint formed from a single bound constraint
           replacated according to a SampleGenerator.

*/
namespace ROL {

template <class Real>
class SimulatedBoundConstraint : public BoundConstraint<Real> {
private:
  const Ptr<SampleGenerator<Real>> sampler_;
  const Ptr<BoundConstraint<Real>> bnd_;
  Ptr<Vector<Real>> l_;
  Ptr<Vector<Real>> u_;

  const Vector<Real>& getVector(const Vector<Real> &x, int k) const {
    try {
      return *(dynamic_cast<const SimulatedVector<Real>&>(x).get(k));
    }
    catch (const std::bad_cast &e) {
      return x;
    }
  }

  Vector<Real>& getVector(Vector<Real> &x, int k) const {
    try {
      return *(dynamic_cast<SimulatedVector<Real>&>(x).get(k));
    }
    catch (const std::bad_cast &e) {
      return x;
    }
  }

public:
  ~SimulatedBoundConstraint() {}

  SimulatedBoundConstraint(const Ptr<SampleGenerator<Real>> &sampler,
                           const Ptr<BoundConstraint<Real>> &bnd )
    : sampler_(sampler), bnd_(bnd) {
    int nsamp = sampler_->numMySamples();
    std::vector<Ptr<Vector<Real>>> lvec(nsamp), uvec(nsamp);
    for ( int k=0; k<sampler_->numMySamples(); ++k) {
      lvec[k] = bnd_->getLowerBound()->clone();
      lvec[k]->set(*bnd_->getLowerBound());
      uvec[k] = bnd_->getUpperBound()->clone();
      uvec[k]->set(*bnd_->getUpperBound());
    }
    l_ = makePtr<SimulatedVector<Real>>(lvec,sampler_->getBatchManager());
    u_ = makePtr<SimulatedVector<Real>>(uvec,sampler_->getBatchManager());
  }

  void project( Vector<Real> &x ) {
    for( int k=0; k<sampler_->numMySamples(); ++k ) {
      if( bnd_->isActivated() ) {
        bnd_->project(getVector(x,k));
      }
    }
  }

  void projectInterior( Vector<Real> &x ) {
    for( int k=0; k<sampler_->numMySamples(); ++k ) {
      if( bnd_->isActivated() ) {
        bnd_->projectInterior(getVector(x,k));
      }
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneUpperActive(getVector(v,k),getVector(x,k),eps);
      }
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneUpperActive(getVector(v,k),getVector(g,k),getVector(x,k),xeps,geps);
      }
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
   if( bnd_->isActivated() ) {
     for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneLowerActive(getVector(v,k),getVector(x,k),eps);
      }
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneLowerActive(getVector(v,k),getVector(g,k),getVector(x,k),xeps,geps);
      }
    }
  }

  const Ptr<const Vector<Real>> getLowerBound( void ) const {
    return l_;
  }

  const Ptr<const Vector<Real>> getUpperBound( void ) const {
    return u_;
  }

  bool isFeasible( const Vector<Real> &v ) {
    bool feasible = true;
    if(bnd_->isActivated()) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        feasible = feasible && bnd_->isFeasible(getVector(v,k));
      }
    }
    return feasible;
  }

  void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->applyInverseScalingFunction(getVector(dv,k),getVector(v,k),getVector(x,k),getVector(g,k));
      }
    }
  }

  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->applyScalingFunctionJacobian(getVector(dv,k),getVector(v,k),getVector(x,k),getVector(g,k));
      }
    }
  }
}; // class SimulatedBoundConstraint
} // namespace ROL

#endif
