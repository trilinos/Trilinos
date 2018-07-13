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
  const ROL::Ptr<SampleGenerator<Real> > sampler_;
  const ROL::Ptr<BoundConstraint<Real> > bnd_;
  ROL::Ptr<Vector<Real> > l_;
  ROL::Ptr<Vector<Real> > u_;

  const Vector<Real>& getVector(const Vector<Real> &x, int k) const {
    try {
      return *(dynamic_cast<const SimulatedVector<Real>&>(x).get(k));
    }
    catch (const std::bad_cast &e) {
      return x;
    }
  }

  Vector<Real>& getVector(Vector<Real> &x, int k) {
    try {
      return *(dynamic_cast<SimulatedVector<Real>&>(x).get(k));
    }
    catch (const std::bad_cast &e) {
      return x;
    }
  }
 
public:
  ~SimulatedBoundConstraint() {}

  SimulatedBoundConstraint(const ROL::Ptr<SampleGenerator<Real> > &sampler,
                           const ROL::Ptr<BoundConstraint<Real> > &bnd )
    : sampler_(sampler), bnd_(bnd) {
    int nsamp = sampler_->numMySamples();
    std::vector<ROL::Ptr<Vector<Real> > > lvec(nsamp), uvec(nsamp);
    for ( int k=0; k<sampler_->numMySamples(); ++k) {
      lvec[k] = bnd_->getLowerBound()->clone();
      lvec[k]->set(*bnd_->getLowerBound());
      uvec[k] = bnd_->getUpperBound()->clone();
      uvec[k]->set(*bnd_->getUpperBound());
    }
    l_ = ROL::makePtr<SimulatedVector<Real>>(lvec,sampler_->getBatchManager());
    u_ = ROL::makePtr<SimulatedVector<Real>>(uvec,sampler_->getBatchManager());
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

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneUpperActive(getVector(v,k),getVector(x,k),eps);
      }
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneUpperActive(getVector(v,k),getVector(g,k),getVector(x,k),eps);
      }
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
   if( bnd_->isActivated() ) {
     for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneLowerActive(getVector(v,k),getVector(x,k),eps);
      }
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if( bnd_->isActivated() ) {
      for( int k=0; k<sampler_->numMySamples(); ++k ) {
        bnd_->pruneLowerActive(getVector(v,k),getVector(g,k),getVector(x,k),eps);
      }
    }
  }
 
  const ROL::Ptr<const Vector<Real> > getLowerBound( void ) const {
    return l_;
  }

  const ROL::Ptr<const Vector<Real> > getUpperBound( void ) const {
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
}; // class SimulatedBoundConstraint
} // namespace ROL

#endif
