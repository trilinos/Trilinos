// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_EXTRACTQP_H
#define ROL_PDEOPT_EXTRACTQP_H

#include "ROL_HelperFunctions.hpp"
#include "ROL_Teuchos_LinearOperator.hpp"
#include "ROL_QuadraticObjective.hpp"
#include "ROL_PEBBL_IntegerProblem.hpp"
#include "../../TOOLS/pdevector.hpp"

template<class Real>
class extractQP {
private:
  const ROL::Ptr<ROL::Objective<Real>>       input_obj_;
  const ROL::Ptr<ROL::Vector<Real>>          input_x_;
  const ROL::Ptr<ROL::BoundConstraint<Real>> input_bnd_;

  ROL::Ptr<ROL::LinearOperator<Real>> H_;
  ROL::Ptr<ROL::Vector<Real>> g_;
  Real c_;
  ROL::Ptr<ROL::Objective<Real>> obj_;

  ROL::Ptr<ROL::Vector<Real>> x_;

  ROL::Ptr<ROL::Vector<Real>> l_;
  ROL::Ptr<ROL::Vector<Real>> u_;
  ROL::Ptr<ROL::BoundConstraint<Real>> bnd_;

  ROL::Ptr<ROL::Vector<Real>> dual_;
  ROL::Ptr<ROL::Vector<Real>> zero_;

  // Get hessian and wrap as TeuchosLinearOperator
  void computeHessian(void) {
    ROL::Ptr<Teuchos::SerialDenseMatrix<int,Real>> H
      = ROL::makePtr<Teuchos::SerialDenseMatrix<int,Real>>();
    *H = ROL::computeDenseHessian(*input_obj_,*zero_);
    H_ = ROL::makePtr<ROL::TeuchosLinearOperator<int,Real>>(H);
  }

  // Get gradient and wrap as TeuchosVector
  void computeGradient(void) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    input_obj_->gradient(*dual_,*zero_,tol);
    ROL::Ptr<std::vector<Real>> gs
      = ROL::dynamicPtrCast<PDE_OptVector<Real>>(dual_)->getParameter()->getVector();
    int size = gs->size();
    ROL::Ptr<Teuchos::SerialDenseVector<int,Real>> gt
      = ROL::makePtr<Teuchos::SerialDenseVector<int,Real>>(size);
    for (int i = 0; i < size; ++i) {
      (*gt)(i) = (*gs)[i];
    }
    g_ = ROL::makePtr<ROL::TeuchosVector<int,Real>>(gt);
  }

  // Get value at z = 0
  void computeValue(void) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    c_ = input_obj_->value(*zero_,tol);
  }

  // Build QuadraticObjective
  void buildObjective(void) {
    computeHessian();
    computeGradient();
    computeValue();
    obj_ = ROL::makePtr<ROL::QuadraticObjective<Real>>(H_,g_,c_);
  }

  // Wrap lower bound as TeuchosVectr
  void computeLowerBound(void) {
    ROL::Ptr<const ROL::Vector<Real>> l = input_bnd_->getLowerBound();
    ROL::Ptr<const std::vector<Real>> ls
      = ROL::dynamicPtrCast<const PDE_OptVector<Real>>(l)->getParameter()->getVector();
    int size = ls->size();
    ROL::Ptr<Teuchos::SerialDenseVector<int,Real>> lt
      = ROL::makePtr<Teuchos::SerialDenseVector<int,Real>>(size);
    for (int i = 0; i < size; ++i) {
      (*lt)(i) = (*ls)[i];
    }
    l_ = ROL::makePtr<ROL::TeuchosVector<int,Real>>(lt);
  }

  // Wrap upper bound as TeuchosVectr
  void computeUpperBound(void) {
    ROL::Ptr<const ROL::Vector<Real>> u = input_bnd_->getUpperBound();
    ROL::Ptr<const std::vector<Real>> us
      = ROL::dynamicPtrCast<const PDE_OptVector<Real>>(u)->getParameter()->getVector();
    int size = us->size();
    ROL::Ptr<Teuchos::SerialDenseVector<int,Real>> ut
      = ROL::makePtr<Teuchos::SerialDenseVector<int,Real>>(size);
    for (int i = 0; i < size; ++i) {
      (*ut)(i) = (*us)[i];
    }
    u_ = ROL::makePtr<ROL::TeuchosVector<int,Real>>(ut);
  }

  // Build bound constraint
  void buildBound(void) {
    computeLowerBound();
    computeUpperBound();
    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(l_,u_);
  }

  void buildSolution(void) {
    ROL::Ptr<const std::vector<Real>> xs
      = ROL::dynamicPtrCast<const PDE_OptVector<Real>>(input_x_)->getParameter()->getVector();
    int size = xs->size();
    ROL::Ptr<Teuchos::SerialDenseVector<int,Real>> xt
      = ROL::makePtr<Teuchos::SerialDenseVector<int,Real>>(size);
    for (int i = 0; i < size; ++i) {
      (*xt)(i) = (*xs)[i];
    }
    x_ = ROL::makePtr<ROL::TeuchosVector<int,Real>>(xt);
  }

public:
  extractQP(const ROL::Ptr<ROL::Objective<Real>>       &obj,
            const ROL::Ptr<ROL::Vector<Real>>          &x,
            const ROL::Ptr<ROL::BoundConstraint<Real>> &bnd)
    : input_obj_(obj), input_x_(x), input_bnd_(bnd) {
    dual_ = input_x_->dual().clone();
    zero_ = input_x_->clone(); zero_->zero();
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> operator()(void) {
    buildObjective();
    buildSolution();
    buildBound();
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj_,x_);
    problem->addBoundConstraint(bnd_);
    return problem;
  }
};

#endif
