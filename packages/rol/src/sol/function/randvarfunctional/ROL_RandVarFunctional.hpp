// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RANDVARFUNCTIONAL_HPP
#define ROL_RANDVARFUNCTIONAL_HPP

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_ScalarController.hpp"
#include "ROL_VectorController.hpp"

/** @ingroup stochastic_group
    \class ROL::RandVarFunctional
    \brief Provides the interface to implement any functional that maps a
           random variable to a (extended) real number.

    Let \f$(\Omega,\mathcal{F},\mathbb{P})\f$ be a probability space.
    Here, \f$\Omega\f$ is the set of outcomes,
    \f$\mathcal{F}\subseteq 2^\Omega\f$ is a \f$\sigma\f$-algebra of events and
    \f$\mathbb{P}:\mathcal{F}\to[0,1]\f$ is a probability measure.  Moreover,
    let \f$\mathcal{X}\f$ be a class of random variables.
    A ``random variable functional'' is an extended real-valued functional that
    associates numerical values to random variables, i.e.,
    \f$\mathcal{R}:\mathcal{X}\to\mathbb{R}\cup\{+\infty\}\f$.  In most cases,
    \f$\mathcal{X} = L^p(\Omega,\mathcal{F},\mathbb{P})\f$ for some
    \f$1\le p\le \infty\f$.

    ROL's random variable functional base class is written in a way to exploit
    parallel sampling.  General risk measures may depend on global information
    such as the expected value of a random variable, \f$\mathbb{E}[X]\f$.  Thus,
    ROL::RandVarFunctional contains functions to update intermediate information
    and to compute desired quantities such as values, gradients and
    Hessians applied to vectors.
*/

namespace ROL {

template<class Real>
class RandVarFunctional {
private:
  bool storage_;
  bool storage_hessvec_;
  Ptr<ScalarController<Real>> value_storage_;
  Ptr<VectorController<Real>> gradient_storage_;
  Ptr<ScalarController<Real>> gradvec_storage_;
  Ptr<VectorController<Real>> hessvec_storage_;

protected:
  Real val_;
  Real gv_;
  Ptr<Vector<Real> > g_;
  Ptr<Vector<Real> > hv_;
  Ptr<Vector<Real> > dualVector_;
  bool firstReset_;

  std::vector<Real> point_;
  Real weight_;

  // Evaluate objective function at current parameter
  Real computeValue(Objective<Real> &obj, const Vector<Real> &x,
                    Real &tol) {
    Real val(0);
    bool isComputed = false;
    if (storage_) {
      isComputed = value_storage_->get(val,point_);
    }
    if (!isComputed || !storage_) {
      obj.setParameter(point_);
      val = obj.value(x,tol);
      if (storage_) {
        value_storage_->set(val,point_);
      }
    }
    return val;
  }
 
  // Evaluate gradient of objective function at current parameter
  void computeGradient(Vector<Real> &g, Objective<Real> &obj,
                       const Vector<Real> &x, Real &tol) {
    bool isComputed = false;
    if (storage_) {
      isComputed = gradient_storage_->get(g,point_);
    }
    if (!isComputed || !storage_) {
      obj.setParameter(point_);
      obj.gradient(g,x,tol);
      if ( storage_ ) {
        gradient_storage_->set(g,point_);
      }
    }
  }

  // Evaluate Gradient-times-a-vector at current parameter
  Real computeGradVec(Vector<Real> &g, Objective<Real> &obj,
                      const Vector<Real> &v, const Vector<Real> &x,
                      Real &tol) {
    Real gv(0);
    computeGradient(g,obj,x,tol);
    bool isComputed = false;
    if (storage_hessvec_) {
      isComputed = gradvec_storage_->get(gv,point_);
    }
    if (!isComputed || !storage_hessvec_) {
      //gv = g.dot(v.dual());
      gv = g.apply(v);
      if (storage_hessvec_) {
        gradvec_storage_->set(gv,point_);
      }
    }
    return gv;
  }

  // Evaluate Hessian-times-a-vector at current parameter
  void computeHessVec(Vector<Real> &hv, Objective<Real> &obj,
                      const Vector<Real> &v, const Vector<Real> &x,
                      Real &tol) {
    bool isComputed = false;
    if (storage_hessvec_) {
      isComputed = hessvec_storage_->get(hv,point_);
    }
    if (!isComputed || !storage_hessvec_) {
      obj.setParameter(point_);
      obj.hessVec(hv,v,x,tol);
      if (storage_hessvec_) {
        hessvec_storage_->set(hv,point_);
      }
    }
  }

public:
  virtual ~RandVarFunctional() {}

  RandVarFunctional(void) : storage_(false), storage_hessvec_(false),
                            value_storage_(nullPtr),
                            gradient_storage_(nullPtr),
                            gradvec_storage_(nullPtr),
                            hessvec_storage_(nullPtr),
                            val_(0), gv_(0), firstReset_(true),
                            point_({}), weight_(0) {}

  void useStorage(bool storage) {
    storage_ = storage;
    if (storage) {
      if (value_storage_ == nullPtr) {
        value_storage_    = makePtr<ScalarController<Real>>();
      }
      if (gradient_storage_ == nullPtr) {
        gradient_storage_ = makePtr<VectorController<Real>>();
      }
    }
  }

  void useHessVecStorage(bool storage) {
    storage_hessvec_ = storage;
    if (storage) {
      useStorage(storage);
      if (gradvec_storage_ == nullPtr) {
        gradvec_storage_ = makePtr<ScalarController<Real>>();
      }
      if (hessvec_storage_ == nullPtr) {
        hessvec_storage_ = makePtr<VectorController<Real>>();
      }
    }
  }

  virtual void setStorage(const Ptr<ScalarController<Real>> &value_storage,
                          const Ptr<VectorController<Real>> &gradient_storage) {
    value_storage_    = value_storage;
    gradient_storage_ = gradient_storage;
    useStorage(true);
  }

  virtual void setHessVecStorage(const Ptr<ScalarController<Real>> &gradvec_storage,
                                 const Ptr<VectorController<Real>> &hessvec_storage) {
    gradvec_storage_ = gradvec_storage;
    hessvec_storage_ = hessvec_storage;
    useHessVecStorage(true);
  }

  /** \brief Reset internal storage.

             @param[in]   x  is a vector used for initializing storage
  */
  virtual void resetStorage(bool flag = true) {
    if (storage_) {
      value_storage_->objectiveUpdate();
      if (flag) {
        gradient_storage_->objectiveUpdate();
        if (storage_hessvec_) {
          gradvec_storage_->objectiveUpdate();
          hessvec_storage_->objectiveUpdate();
        }
      }
    }
  }
  virtual void resetStorage(UpdateType type) {
    if (storage_) {
      value_storage_->objectiveUpdate(type);
      gradient_storage_->objectiveUpdate(type);
      if (storage_hessvec_) {
        gradvec_storage_->objectiveUpdate(type);
        hessvec_storage_->objectiveUpdate(type);
      }
    }
  }

  /** \brief Initialize temporary variables.

             @param[in]   x  is a vector used for initializing storage
  */
  virtual void initialize(const Vector<Real> &x) {
    // Create memory for class members
    if ( firstReset_ ) {
      g_          = x.dual().clone();
      hv_         = x.dual().clone();
      dualVector_ = x.dual().clone();
      firstReset_ = false;
    }
    // Zero member variables
    const Real zero(0);
    val_ = zero; gv_ = zero;
    g_->zero(); hv_->zero(); dualVector_->zero();
    if (storage_hessvec_) {
      gradvec_storage_->reset();
      hessvec_storage_->reset();
    }
  }

  virtual void setSample(const std::vector<Real> &point, const Real weight) {
    point_.assign(point.begin(),point.end());
    weight_ = weight;
  }

  /** \brief Compute statistic.

      @param[in]    xstat   is a ROL::Ptr to a std::vector containing the
                            statistic vector
  */
  virtual Real computeStatistic(const Ptr<const std::vector<Real>> &xstat) const {
    Real stat(0);
    if (xstat != nullPtr && !xstat->empty()) {
      stat = (*xstat)[0];
    }
    return stat;
  }

  /** \brief Update internal storage for value computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void updateValue(Objective<Real>         &obj,
                           const Vector<Real>      &x,
                           const std::vector<Real> &xstat,
                           Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_ += weight_ * val;
  }

  /** \brief Update internal risk measure storage for gradient computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    g        is the gradient of the random variable objective
                             function at the current sample point
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void updateGradient(Objective<Real>         &obj,
                              const Vector<Real>      &x,
                              const std::vector<Real> &xstat,
                              Real                    &tol) {
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_,*dualVector_);
  }

  /** \brief Update internal risk measure storage for Hessian-time-a-vector computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    g        is the gradient of the random variable objective
                             function at the current sample point
      @param[in]    gv       is the gradient of the random variable objective
                             function at the current sample point applied to
                             the vector v0
      @param[in]    hv       is the Hessian of the random variable objective
                             function at the current sample point applied to
                             the vector v0
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void updateHessVec(Objective<Real>         &obj,
                             const Vector<Real>      &v,
                             const std::vector<Real> &vstat,
                             const Vector<Real>      &x,
                             const std::vector<Real> &xstat,
                             Real                    &tol) {
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_,*dualVector_);
  }

  /** \brief Return risk measure value.

      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getValue returns \f$\mathcal{R}(f(x_0))\f$ where \f$f(x_0)\f$
      denotes the random variable objective function evaluated at \f$x_0\f$.
  */
  virtual Real getValue(const Vector<Real>      &x,
                        const std::vector<Real> &xstat,
                        SampleGenerator<Real>   &sampler) {
    Real val(0);
    sampler.sumAll(&val_,&val,1);
    return val;
  }

  /** \brief Return risk measure (sub)gradient.

      @param[out]   g       is the (sub)gradient of the risk measure
      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getGradient returns \f$\theta\in\partial\mathcal{R}(f(x_0))\f$
      where \f$f(x_0)\f$ denotes the random variable objective function evaluated
      at \f$x_0\f$ and \f$\partial\mathcal{R}(X)\f$ denotes the subdifferential
      of \f$\mathcal{R}\f$ at \f$X\f$.
  */
  virtual void getGradient(Vector<Real>            &g,
                           std::vector<Real>       &gstat,
                           const Vector<Real>      &x,
                           const std::vector<Real> &xstat,
                           SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*g_,g);
  }

  /** \brief Return risk measure Hessian-times-a-vector.

      @param[out]   hv      is the Hessian-times-a-vector of the risk measure
      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getHessVec returns \f$\nabla^2 \mathcal{R}(f(x_0))v_0\f$
      (if available)
      where \f$f(x_0)\f$ denotes the random variable objective function evaluated
      at \f$x_0\f$.
  */
  virtual void getHessVec(Vector<Real>            &hv,
                          std::vector<Real>       &hvstat,
                          const Vector<Real>      &v,
                          const std::vector<Real> &vstat,
                          const Vector<Real>      &x,
                          const std::vector<Real> &xstat,
                          SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
