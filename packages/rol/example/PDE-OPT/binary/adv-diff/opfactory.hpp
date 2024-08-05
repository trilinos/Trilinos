// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef OPFACTORY_BINARY_ADVDIFF_HPP
#define OPFACTORY_BINARY_ADVDIFF_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "pde_adv_diff.hpp"
#include "qoi_adv_diff.hpp"
#include "mesh_adv_diff.hpp"
#include "extractQP.hpp"

template<class Real>
class BinaryAdvDiffFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  int dim_;

  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<MeshManager_adv_diff<Real>>  mesh_;
  ROL::Ptr<PDE_adv_diff<Real>>          pde_;
  ROL::Ptr<Linear_PDE_Constraint<Real>> con_;
  ROL::Ptr<Assembler<Real>>             assembler_;
  ROL::Ptr<ROL::Vector<Real>>           u_, z_, p_;

public:
  BinaryAdvDiffFactory(ROL::ParameterList                 &pl,
                 const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                 const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    int nx = pl.sublist("Problem").get("Number Controls - X", 3);
    int ny = pl.sublist("Problem").get("Number Controls - Y", 3);
    dim_ = nx*ny;
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    update();
    ROL::Ptr<ROL::Objective<Real>>       obj  = buildObjective();
    ROL::Ptr<ROL::Vector<Real>>          x    = buildSolutionVector();
    ROL::Ptr<ROL::BoundConstraint<Real>> bnd  = buildBoundConstraint();
    ROL::Ptr<ROL::Constraint<Real>>      icon = buildConstraint();
    ROL::Ptr<ROL::Vector<Real>>          imul = buildMultiplier();
    ROL::Ptr<ROL::BoundConstraint<Real>> ibnd = buildSlackBoundConstraint();
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj,x);
    problem->addBoundConstraint(bnd);
    if (ibnd==ROL::nullPtr) problem->addLinearConstraint("Budget",icon,imul);
    else                    problem->addLinearConstraint("Budget",icon,imul,ibnd);
    return problem;
  }

  void update(void) {
    mesh_      = ROL::makePtr<MeshManager_adv_diff<Real>>(pl_);
    pde_       = ROL::makePtr<PDE_adv_diff<Real>>(pl_);
    con_       = ROL::makePtr<Linear_PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
    assembler_ = con_->getAssembler();

    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr;
    u_ptr = assembler_->createStateVector();
    p_ptr = assembler_->createStateVector();
    r_ptr = assembler_->createResidualVector();
    ROL::Ptr<ROL::Vector<Real>> up, pp, rp, zp;
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(u_ptr,pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(p_ptr,pde_,assembler_,pl_);
    z_ = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_));
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    std::vector<ROL::Ptr<QoI<Real>>> qoi_vec(1,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<Real>>(pde_->getFE());
    Real stateCost = pl_.sublist("Problem").get("State Cost",4.0);
    std::vector<Real> wts = {stateCost};
    ROL::Ptr<ROL::Objective_SimOpt<Real>> obj
      = ROL::makePtr<PDE_Objective<Real>>(qoi_vec,wts,assembler_);
    bool storage = pl_.sublist("Problem").get("Use Storage",true);
    return ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj, con_, u_, z_, p_, storage, false);
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return z_;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    ROL::Ptr<ROL::Vector<Real>> zlop
      = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_,0.0));
    ROL::Ptr<ROL::Vector<Real>> zhip
      = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_,1.0));
    return ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
  }

  ROL::Ptr<ROL::Constraint<Real>> buildConstraint(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    Real budget(0);
    if (!useIneq) budget = pl_.sublist("Problem").get("Control Cost", 4.0);
    ROL::Ptr<QoI<Real>> qoi
      = ROL::makePtr<QoI_Control_Cost_adv_diff<Real>>(budget);
    return ROL::makePtr<IntegralOptConstraint<Real>>(qoi,assembler_);
  }

  ROL::Ptr<ROL::Vector<Real>> buildMultiplier(void) {
    return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildSlackBoundConstraint(void) {
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      Real budget = pl_.sublist("Problem").get("Control Cost", 4.0);
      ROL::Ptr<ROL::Vector<Real>> klop, khip;
      klop = ROL::makePtr<ROL::StdVector<Real>>(1,static_cast<Real>(0));
      khip = ROL::makePtr<ROL::StdVector<Real>>(1,budget);
      return ROL::makePtr<ROL::Bounds<Real>>(klop,khip);
    }
    return ROL::nullPtr;
  }
};

template<class Real>
class BinaryAdvDiffQPFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<BinaryAdvDiffFactory<Real>> factory_;
  ROL::Ptr<extractQP<Real>> extract_;

public:
  BinaryAdvDiffQPFactory(ROL::ParameterList                 &pl,
                   const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                   const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    factory_ = ROL::makePtr<BinaryAdvDiffFactory<Real>>(pl_,comm_,os_);
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    factory_->update();
    extract_ = ROL::makePtr<extractQP<Real>>(factory_->buildObjective(),
                                             factory_->buildSolutionVector(),
                                             factory_->buildBoundConstraint(),
                                             factory_->buildConstraint(),
                                             factory_->buildMultiplier(),
                                             factory_->buildSlackBoundConstraint());
    return (*extract_)();
  }
};
  

#endif
