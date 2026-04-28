// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_MEANVALUEOBJECTIVE_HPP
#define ROL_OED_MEANVALUEOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_OED_DesignVector.hpp"
#include "ROL_OED_ObjectiveArray.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {
namespace OED {

template<typename Real>

class MeanValueObjective : public Objective<Real> {
    private: // empty for now
        const Ptr<ObjectiveArray<Real>> objVec_;
        const Ptr<Vector<Real>> dwa_;
        // std::vector<Real> weights_;

    public: 
        // MeanValueObjective(const Ptr<ObjectiveArray<Real>> &objVec, const std::vector<Real> & weights)
        // : objVec_(objVec), dwa_(((objVec->buildDesignVector())->dual()).clone()), weights_(weights) {}

        MeanValueObjective(const Ptr<ObjectiveArray<Real>> &objVec)
        : objVec_(objVec), dwa_(((objVec->buildDesignVector())->dual()).clone()) {
            // weights_.resize(objVec_->numObjectives(), 1.0/static_cast<Real>(objVec_->numObjectives())); //cast numObjectives to real number 
        }

    void update(const Vector<Real> &x, UpdateType type, int iter = -1 ) override {
      for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) 
        objVec_->getObjective(i)->update(x, type, iter);
    }

    void setParameter(const std::vector<Real> &param) {
      Objective<Real>::setParameter(param);
      for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) 
        objVec_->getObjective(i)->setParameter(param);
    }
 
    Real value( const Vector<Real> &x, Real &tol ) override {
      Real val(0);
      for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) 
        val+= objVec_->getWeight(i)*objVec_->getObjective(i)->value(x,tol);
      return val; //one line in for doesn't need scoping braces
    }

    //gradient -- returns vector
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override {
      g.zero(); //initialize at 0
      for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) {
        objVec_->getObjective(i)->gradient(*dwa_,x,tol);
        g.axpy(objVec_->getWeight(i), *dwa_);
      }
    }

    //hessvec looks like gradient but with modified inputs
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override {
      hv.zero(); //initialize at 0
      for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) {
        objVec_->getObjective(i)->hessVec(*dwa_,v,x,tol);
        hv.axpy(objVec_->getWeight(i), *dwa_);
      }
    }
    
}; //class MeanValueObjective

} // end ROL namespace
} // end OED namespace

#endif //ending the header search space command
