// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUND_CONSTRAINT_PARTITIONED_H
#define ROL_BOUND_CONSTRAINT_PARTITIONED_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::BoundConstraint_Partitioned
    \brief A composite composite BoundConstraint formed
           from bound constraints on subvectors of a PartitionedVector

*/
namespace ROL {

template<typename Real>
class BoundConstraint_Partitioned : public BoundConstraint<Real> {

  typedef Vector<Real>                          V;
  typedef PartitionedVector<Real>               PV;
  typedef typename std::vector<Real>::size_type uint;

private:
  std::vector<Ptr<BoundConstraint<Real>>> bnd_;

  using BoundConstraint<Real>::lower_;
  using BoundConstraint<Real>::upper_;
  Ptr<V> l_;
  Ptr<V> u_;

  uint dim_;

  bool hasLvec_;
  bool hasUvec_;

public:
  ~BoundConstraint_Partitioned() {}

  BoundConstraint_Partitioned(const std::vector<Ptr<BoundConstraint<Real>>> &bnd,
                              const std::vector<Ptr<Vector<Real>>> &x)
    : bnd_(bnd), dim_(bnd.size()), hasLvec_(true), hasUvec_(true) {
    BoundConstraint<Real>::deactivate();
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        BoundConstraint<Real>::activate();
        break;
      }
    }
    std::vector<Ptr<Vector<Real>>> lp(dim_);
    std::vector<Ptr<Vector<Real>>> up(dim_);
    for( uint k=0; k<dim_; ++k ) {
      try {
        lp[k] = x[k]->clone();
        if (bnd_[k]->isLowerActivated()) {
          lp[k]->set(*bnd_[k]->getLowerBound());
        }
        else {
          lp[k]->setScalar(-BoundConstraint<Real>::computeInf(*x[k]));
        }
      }
      catch (std::exception &e1) {
        try {
          lp[k] = x[k]->clone();
          lp[k]->setScalar(-BoundConstraint<Real>::computeInf(*x[k]));
        }
        catch (std::exception &e2) {
          lp[k] = nullPtr;
          hasLvec_ = false;
        }
      }
      try {
        up[k] = x[k]->clone();
        if (bnd_[k]->isUpperActivated()) {
          up[k]->set(*bnd_[k]->getUpperBound());
        }
        else {
          up[k]->setScalar( BoundConstraint<Real>::computeInf(*x[k]));
        }
      }
      catch (std::exception &e1) {
        try {
          up[k] = x[k]->clone();
          up[k]->setScalar( BoundConstraint<Real>::computeInf(*x[k]));
        }
        catch (std::exception &e2) {
          up[k] = nullPtr;
          hasUvec_ = false;
        }
      }
    }
    if (hasLvec_) {
      lower_ = makePtr<PV>(lp);
    }
    if (hasUvec_) {
      upper_ = makePtr<PV>(up);
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
  }

  void project( Vector<Real> &x ) {
    PV &xpv = dynamic_cast<PV&>(x);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->project(*xpv.get(k));
      }
    }
  }

  void projectInterior( Vector<Real> &x ) {
    PV &xpv = dynamic_cast<PV&>(x);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->projectInterior(*xpv.get(k));
      }
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
          PV &vpv = dynamic_cast<PV&>(v);
    const PV &xpv = dynamic_cast<const PV&>(x);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->pruneUpperActive(*(vpv.get(k)),*(xpv.get(k)),eps);
      }
    }
 }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
          PV &vpv = dynamic_cast<PV&>(v);
    const PV &gpv = dynamic_cast<const PV&>(g);
    const PV &xpv = dynamic_cast<const PV&>(x);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->pruneUpperActive(*(vpv.get(k)),*(gpv.get(k)),*(xpv.get(k)),xeps,geps);
      }
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
          PV &vpv = dynamic_cast<PV&>(v);
    const PV &xpv = dynamic_cast<const PV&>(x);
   for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->pruneLowerActive(*(vpv.get(k)),*(xpv.get(k)),eps);
      }
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
          PV &vpv = dynamic_cast<PV&>(v);
    const PV &gpv = dynamic_cast<const PV&>(g);
    const PV &xpv = dynamic_cast<const PV&>(x);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->pruneLowerActive(*(vpv.get(k)),*(gpv.get(k)),*(xpv.get(k)),xeps,geps);
      }
    }
  }

  bool isFeasible( const Vector<Real> &v ) {
    bool feasible = true;
    const PV &vs = dynamic_cast<const PV&>(v);
    for( uint k=0; k<dim_; ++k ) {
      if(bnd_[k]->isActivated()) {
        feasible = feasible && bnd_[k]->isFeasible(*(vs.get(k)));
      }
    }
    return feasible;
  }

  void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
          PV &dvpv = dynamic_cast<PV&>(dv);
    const PV &vpv  = dynamic_cast<const PV&>(v);
    const PV &xpv  = dynamic_cast<const PV&>(x);
    const PV &gpv  = dynamic_cast<const PV&>(g);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->applyInverseScalingFunction(*(dvpv.get(k)),*(vpv.get(k)),*(xpv.get(k)),*(gpv.get(k)));
      }
    }
  }

  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
          PV &dvpv = dynamic_cast<PV&>(dv);
    const PV &vpv  = dynamic_cast<const PV&>(v);
    const PV &xpv  = dynamic_cast<const PV&>(x);
    const PV &gpv  = dynamic_cast<const PV&>(g);
    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        bnd_[k]->applyScalingFunctionJacobian(*(dvpv.get(k)),*(vpv.get(k)),*(xpv.get(k)),*(gpv.get(k)));
      }
    }
  }
}; // class BoundConstraint_Partitioned



template<typename Real>
Ptr<BoundConstraint<Real>>
CreateBoundConstraint_Partitioned( const Ptr<BoundConstraint<Real>> &bnd1,
                                   const Ptr<BoundConstraint<Real>> &bnd2 ) {


  typedef BoundConstraint<Real>             BND;
  typedef BoundConstraint_Partitioned<Real> BNDP;
  Ptr<BND> temp[] = {bnd1, bnd2};
  return makePtr<BNDP>( std::vector<Ptr<BND>>(temp,temp+2) );
}


} // namespace ROL

#endif
