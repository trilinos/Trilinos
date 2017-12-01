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

#ifndef ROL_CONVEXCOMBINATIONRISKMEASURE_HPP
#define ROL_CONVEXCOMBINATIONRISKMEASURE_HPP

#include "ROL_RiskMeasureFactory.hpp"

/** @ingroup risk_group
    \class ROL::ConvexCombinationRiskMeasure
    \brief Provides an interface for a convex combination of
           risk measures.

    This function provides the capability to produce a convex combination
    of risk measure, i.e.,
    \f[
       \mathcal{R}(X) = \sum_{k=1}^n \lambda_k \mathcal{R}_k(X)
    \f]
    where \f$\mathcal{R}_k\f$ are risk measures and \f$\lambda_k \ge 0\f$
    with \f$\lambda_1 + \ldots + \lambda_n = 1\f$.  In general,
    \f$\mathcal{R}\f$ is not law-invariant or coherent unless each
    \f$\mathcal{R}_k\f$ is.
*/

namespace ROL {

template<class Real>
class ConvexCombinationRiskMeasure : public RiskMeasure<Real> {
private:
  typedef typename std::vector<Real>::size_type uint;

  std::vector<Real> lambda_;
  std::vector<ROL::SharedPointer<ROL::ParameterList> > parlist_;
  std::vector<ROL::SharedPointer<RiskMeasure<Real> > > risk_;
  uint size_;

  ROL::SharedPointer<Vector<Real> > dualVector0_;
  bool firstReset_;

  void checkInputs(void) const {
    uint lSize = lambda_.size(), rSize = risk_.size();
    TEUCHOS_TEST_FOR_EXCEPTION((lSize!=rSize),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Convex combination parameter and risk measure arrays have different sizes!");
    Real sum(0), zero(0), one(1);
    for (uint i = 0; i < lSize; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION((lambda_[i]>one || lambda_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Element of convex combination parameter array out of range!");
      TEUCHOS_TEST_FOR_EXCEPTION(risk_[i] == ROL::nullPointer, std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Risk measure pointer is null!");
      sum += lambda_[i];
    }
    TEUCHOS_TEST_FOR_EXCEPTION((std::abs(sum-one) > std::sqrt(ROL_EPSILON<Real>())),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Coefficients do not sum to one!");
  }

public:
  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Convex Combination Risk Measure" and
      within the "Convex Combination Risk Measure" sublist should have the following parameters
      \li "Convex Combination Parameters" (greater than 0 and sum to 1)
      \li Sublists labeled 1 to n with risk measure definitions.
  */
  ConvexCombinationRiskMeasure(ROL::ParameterList &parlist)
    : RiskMeasure<Real>(), size_(0), firstReset_(true) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Convex Combination Risk Measure");
    // Get convex combination parameters
    lambda_ = ROL::getArrayFromStringParameter<Real>(list,"Convex Combination Parameters");
    size_ = lambda_.size();
    // Build risk measures
    risk_.clear(); risk_.resize(size_,ROL::nullPointer);
    parlist_.clear(); parlist_.resize(size_,ROL::nullPointer);
    for (uint i = 0; i < size_; ++i) {
      std::ostringstream convert;
      convert << i;
      std::string si = convert.str();
      ROL::ParameterList &ilist = list.sublist(si);

      std::string name = ilist.get<std::string>("Name");
      parlist_[i] = ROL::makeShared<ROL::ParameterList>();
      parlist_[i]->sublist("SOL").sublist("Risk Measure").set("Name",name);
      parlist_[i]->sublist("SOL").sublist("Risk Measure").sublist(name) = ilist;

      risk_[i] = RiskMeasureFactory<Real>(*parlist_[i]);
    }
    // Check inputs
    checkInputs();
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x) {
    ROL::SharedPointer<std::vector<Real> > stati;
    int N = 0, Ni = 0;
    // Must make x a risk vector with appropriate statistic
    const RiskVector<Real> &xr = dynamic_cast<const RiskVector<Real>&>(x);
    ROL::SharedPointer<const Vector<Real> > xptr = xr.getVector();
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    ROL::SharedPointer<const std::vector<Real> > stat
      = xr.getStatistic(comp,index);
    x0 = ROL::constPointerCast<Vector<Real> >(xptr);
    for (uint i = 0; i < size_; ++i) {
      // Build temporary risk vector
      RiskVector<Real> xri(parlist_[i],x0);
      // Set statistic from original risk vector
      stati = xri.getStatistic(0);
      if (stati != ROL::nullPointer) {
        Ni = stati->size();
        for (int j = 0; j < Ni; ++j) {
          (*stati)[j] = (*stat)[N+j];
        }
        xri.setStatistic(*stati,0);
      }
      else {
        Ni = 0;
      }
      N += Ni;
      // Reset current risk measure
      risk_[i]->reset(x0,xri);
    }
    if (firstReset_) {
      dualVector0_ = x0->dual().clone();
      firstReset_ = false;
    }
    dualVector0_->zero();
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x,
             ROL::SharedPointer<Vector<Real> > &v0, const Vector<Real> &v) {
    ConvexCombinationRiskMeasure<Real>::reset(x0,x);
    ROL::SharedPointer<std::vector<Real> > xstati, vstati;
    int N = 0, Ni = 0;
    // Must make x and v risk vectors with appropriate statistics
    const RiskVector<Real> &xr = dynamic_cast<const RiskVector<Real>&>(x);
    const RiskVector<Real> &vr = dynamic_cast<const RiskVector<Real>&>(v);
    ROL::SharedPointer<const Vector<Real> > xptr = xr.getVector();
    ROL::SharedPointer<const Vector<Real> > vptr = vr.getVector();
    x0 = ROL::constPointerCast<Vector<Real> >(xptr);
    v0 = ROL::constPointerCast<Vector<Real> >(vptr);
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    ROL::SharedPointer<const std::vector<Real> > xstat
      = xr.getStatistic(comp,index);
    ROL::SharedPointer<const std::vector<Real> > vstat
      = vr.getStatistic(comp,index);
    for (uint i = 0; i < size_; ++i) {
      // Build temporary risk vector
      RiskVector<Real> xri(parlist_[i],x0), vri(parlist_[i],v0);
      // Set statistic from original risk vector
      xstati = xri.getStatistic(0);
      vstati = vri.getStatistic(0);
      if (xstati != ROL::nullPointer) {
        Ni = xstati->size();
        for (int j = 0; j < Ni; ++j) {
          (*xstati)[j] = (*xstat)[N+j];
          (*vstati)[j] = (*vstat)[N+j];
        }
        xri.setStatistic(*xstati,0);
        vri.setStatistic(*vstati,0);
      }
      else {
        Ni = 0;
      }
      N += Ni;
      // Reset current risk measure
      risk_[i]->reset(x0,xri,v0,vri);
    }
    if (firstReset_) {
      dualVector0_ = x0->dual().clone();
      firstReset_ = false;
    }
    dualVector0_->zero();
  }

  void update(const Real val, const Real weight) {
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->update(val,weight);
    }
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val(0);
    for (uint i = 0; i < size_; ++i) {
      val += lambda_[i]*risk_[i]->getValue(sampler);
    }
    return val;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->update(val,g,weight);
    }
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    g.zero();
    // g does not have the correct dimension if it is a risk vector
    RiskVector<Real> &gr = dynamic_cast<RiskVector<Real>&>(g);
    ROL::SharedPointer<std::vector<Real> > stat, stati;
    stat = ROL::makeShared<std::vector<Real>>(0);
    for (uint i = 0; i < size_; ++i) {
      RiskVector<Real> gri(parlist_[i],dualVector0_);
      risk_[i]->getGradient(gri,sampler);
      (gr.getVector())->axpy(lambda_[i],*dualVector0_);
      stati = gri.getStatistic(0);
      if (stati != ROL::nullPointer) {
        for (uint j = 0; j < stati->size(); ++j) {
          stat->push_back(lambda_[i]*(*stati)[j]);
        }
      }
    }
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    gr.setStatistic(*stat,comp,index);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->update(val,g,gv,hv,weight);
    }
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    hv.zero();
    // hv does not have the correct dimension if it is a risk vector
    RiskVector<Real> &hvr = dynamic_cast<RiskVector<Real>&>(hv);
    ROL::SharedPointer<std::vector<Real> > stat, stati;
    stat = ROL::makeShared<std::vector<Real>>(0);
    for (uint i = 0; i < size_; ++i) {
      RiskVector<Real> hvri(parlist_[i],dualVector0_);
      risk_[i]->getHessVec(hvri,sampler);
      (hvr.getVector())->axpy(lambda_[i],*dualVector0_);
      stati = hvri.getStatistic(0);
      if (stati != ROL::nullPointer) {
        for (uint j = 0; j < stati->size(); ++j) {
          stat->push_back(lambda_[i]*(*stati)[j]);
        }
      }
    }
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    hvr.setStatistic(*stat,comp,index);
  }
};

}

#endif
