// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEARREGRESSION_H
#define ROL_LINEARREGRESSION_H

#include "ROL_RegressionError.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_StochasticObjective.hpp"
#include "ROL_ErrorMeasureFactory.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_Problem.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_RiskBoundConstraint.hpp"

/** @ingroup algo_group
    \class ROL::LinearRegression
    \brief Provides the interface to construct linear regression problem.

    ---
*/


namespace ROL {

template <class Real>
class LinearRegression {
private:
  const Ptr<RegressionError<Real>> error_;
  const Ptr<SampleGenerator<Real>> data_;

  Ptr<RandVarFunctional<Real>>   em_;
  Ptr<StochasticObjective<Real>> obj_;
  Ptr<std::vector<Real>>         cdata_;
  Ptr<RiskVector<Real>>          c_;

  Ptr<std::vector<Real>>         lower_;
  Ptr<std::vector<Real>>         upper_;
  Ptr<BoundConstraint<Real>>     bnd_;
  Ptr<RiskBoundConstraint<Real>> rbnd_;

  bool initialized_;

public:
  LinearRegression(const Ptr<SampleGenerator<Real>> &data)
    : error_(makePtr<RegressionError<Real>>()), data_(data),
      lower_(nullPtr), upper_(nullPtr), bnd_(nullPtr), rbnd_(nullPtr),
      initialized_(false) {
    int dim = data_->getMyPoint(0).size();
    cdata_  = makePtr<std::vector<Real>>(dim,0);
    c_      = makePtr<RiskVector<Real>>(makePtr<StdVector<Real>>(cdata_));
  }

  void setErrorMeasure(ROL::ParameterList &parlist, bool reset = false) {
    if (!initialized_ || reset) {
      em_ = ErrorMeasureFactory<Real>(parlist);
      obj_ = makePtr<StochasticObjective<Real>>(error_,em_,data_);
      initialized_ = true;
    }
  }

  void setLowerBound(const std::vector<Real> &lower) {
    lower_ = makePtr<std::vector<Real>>(lower);
  }

  void setUpperBound(const std::vector<Real> &upper) {
    upper_ = makePtr<std::vector<Real>>(upper);
  }

  void reset(void) {
    c_->zero();
    initialized_ = false;
  }

  const Ptr<OptimizationProblem<Real>> getOptimizationProblem(void) {
    if (!initialized_) {
      throw Exception::NotImplemented("ROL::LinearRegression::getOptimizationProblem : setErrorMeasure was not called!");
    }
    if (lower_ != nullPtr && upper_ == nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*lower_,true);
    }
    if (lower_ == nullPtr && upper_ != nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*upper_,false);
    }
    if (lower_ != nullPtr && upper_ != nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*lower_,*upper_);
    }
    if (bnd_ != nullPtr) {
      rbnd_ = makePtr<RiskBoundConstraint<Real>>(bnd_);
      return makePtr<OptimizationProblem<Real>>(obj_,c_,rbnd_);
    }
    return makePtr<OptimizationProblem<Real>>(obj_,c_);
  }

  const Ptr<Problem<Real>> getProblem(void) {
    if (!initialized_) {
      throw Exception::NotImplemented("ROL::LinearRegression::getProblem : setErrorMeasure was not called!");
    }
    Ptr<Problem<Real>> prob
      = makePtr<Problem<Real>>(obj_,c_);
    if (lower_ != nullPtr && upper_ == nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*lower_,true);
    }
    if (lower_ == nullPtr && upper_ != nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*upper_,false);
    }
    if (lower_ != nullPtr && upper_ != nullPtr) {
      bnd_ = makePtr<StdBoundConstraint<Real>>(*lower_,*upper_);
    }
    if (bnd_ != nullPtr) {
      rbnd_ = makePtr<RiskBoundConstraint<Real>>(bnd_);
      prob->addBoundConstraint(rbnd_);
    }
    return prob;
  }

  const Ptr<std::vector<Real>> getCoefficients(void) const {
    return cdata_;
  }

  void print(std::ostream &out = std::cout, const std::string delim = "  ") const {
    int dim = cdata_->size();
    out << std::endl;
    for (int i = 0; i < dim; ++i) {
      out << delim << (*cdata_)[i];
    }
    out << std::endl << std::endl;
  }
}; // class LinearRegression

} // namespace ROL

#endif
