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
#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_PEBBL_IntegerProblem.hpp"
#include "../../TOOLS/pdevector.hpp"

template<class Real>
class extractQP {
private:
  const ROL::Ptr<ROL::Objective<Real>>       input_obj_;
  const ROL::Ptr<ROL::Vector<Real>>          input_x_;
  const ROL::Ptr<ROL::BoundConstraint<Real>> input_bnd_;
  const ROL::Ptr<ROL::Constraint<Real>>      input_icon_;
  const ROL::Ptr<ROL::Vector<Real>>          input_imul_;
  const ROL::Ptr<ROL::BoundConstraint<Real>> input_ibnd_;

  ROL::Ptr<ROL::LinearOperator<Real>> H_;
  ROL::Ptr<ROL::Vector<Real>> g_;
  Real c_;
  ROL::Ptr<ROL::Objective<Real>> obj_;

  ROL::Ptr<ROL::Vector<Real>> x_;

  ROL::Ptr<ROL::Vector<Real>> J_;
  Real b_;
  ROL::Ptr<ROL::Constraint<Real>> icon_;

  ROL::Ptr<ROL::Vector<Real>> imul_;

  ROL::Ptr<ROL::Vector<Real>> il_;
  ROL::Ptr<ROL::Vector<Real>> iu_;
  ROL::Ptr<ROL::BoundConstraint<Real>> ibnd_;

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

  // Get jacobian of knapsack constraint and wrap as TeuchosVector
  void computeJacobian(void) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    ROL::Ptr<ROL::Vector<Real>> one
      = ROL::makePtr<ROL::StdVector<Real>>(1,1.0);
    input_icon_->applyAdjointJacobian(*dual_,*one,*zero_,tol);
    ROL::Ptr<std::vector<Real>> es
      = ROL::dynamicPtrCast<PDE_OptVector<Real>>(dual_)->getParameter()->getVector();
    int size = es->size();
    ROL::Ptr<Teuchos::SerialDenseVector<int,Real>> et
      = ROL::makePtr<Teuchos::SerialDenseVector<int,Real>>(size);
    for (int i = 0; i < size; ++i) {
      (*et)(i) = (*es)[i];
    }
    J_ = ROL::makePtr<ROL::TeuchosVector<int,Real>>(et);
  }

  // Get shift of knapsack constraint
  void computeBudget(void) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    ROL::Ptr<ROL::Vector<Real>> b = input_imul_->dual().clone();
    input_icon_->value(*b,*zero_,tol);
    b_ = -(*ROL::dynamicPtrCast<ROL::StdVector<Real>>(b)->getVector())[0];
  }

  // Build LinearConstraint
  void buildConstraint(void) {
    computeJacobian();
    computeBudget();
    icon_ = ROL::makePtr<ROL::ScalarLinearConstraint<Real>>(J_,b_);
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

  void buildMultiplier(void) {
    ROL::Ptr<std::vector<Real>> ms
      = ROL::dynamicPtrCast<ROL::StdVector<Real>>(input_imul_)->getVector();
    imul_ = ROL::makePtr<ROL::SingletonVector<Real>>((*ms)[0]);
  }

  // Wrap lower bound as TeuchosVectr
  void computeIneqLowerBound(void) {
    ROL::Ptr<const ROL::Vector<Real>> l = input_ibnd_->getLowerBound();
    ROL::Ptr<const std::vector<Real>> ls
      = ROL::dynamicPtrCast<const ROL::StdVector<Real>>(l)->getVector();
    il_ = ROL::makePtr<ROL::SingletonVector<Real>>((*ls)[0]);
  }

  // Wrap upper bound as TeuchosVectr
  void computeIneqUpperBound(void) {
    ROL::Ptr<const ROL::Vector<Real>> u = input_ibnd_->getUpperBound();
    ROL::Ptr<const std::vector<Real>> us
      = ROL::dynamicPtrCast<const ROL::StdVector<Real>>(u)->getVector();
    iu_ = ROL::makePtr<ROL::SingletonVector<Real>>((*us)[0]);
  }

  // Build bound constraint
  void buildIneqBound(void) {
    computeIneqLowerBound();
    computeIneqUpperBound();
    ibnd_ = ROL::makePtr<ROL::Bounds<Real>>(il_,iu_);
  }

public:
  extractQP(const ROL::Ptr<ROL::Objective<Real>>       &obj,
            const ROL::Ptr<ROL::Vector<Real>>          &x,
            const ROL::Ptr<ROL::BoundConstraint<Real>> &bnd,
            const ROL::Ptr<ROL::Constraint<Real>>      &icon,
            const ROL::Ptr<ROL::Vector<Real>>          &imul,
            const ROL::Ptr<ROL::BoundConstraint<Real>> &ibnd = ROL::nullPtr)
    : input_obj_(obj), input_x_(x), input_bnd_(bnd), input_icon_(icon), input_imul_(imul), input_ibnd_(ibnd) {
    dual_ = input_x_->dual().clone();
    zero_ = input_x_->clone(); zero_->zero();
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> operator()(void) {
    buildObjective();
    buildSolution();
    buildBound();
    buildConstraint();
    buildMultiplier();
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj_,x_);
    problem->addBoundConstraint(bnd_);
    if (input_ibnd_ != ROL::nullPtr) {
      buildIneqBound();
      problem->addLinearConstraint("Linear",icon_,imul_,ibnd_);
    }
    else {
      problem->addLinearConstraint("Linear",icon_,imul_);
    }
    return problem;
  }
};

#endif
