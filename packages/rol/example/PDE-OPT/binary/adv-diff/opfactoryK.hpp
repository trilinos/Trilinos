// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef OPFACTORY_BINARY_ADVDIFFK_HPP
#define OPFACTORY_BINARY_ADVDIFFK_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

#include "../../TOOLS/linearpdeconstraintK.hpp"
#include "../../TOOLS/pdeconstraintK.hpp"
#include "../../TOOLS/pdeobjectiveK.hpp"
#include "../../TOOLS/pdevectorK.hpp"
#include "../../TOOLS/integralconstraintK.hpp"
#include "pde_adv_diffK.hpp"
#include "qoi_adv_diffK.hpp"
#include "mesh_adv_diffK.hpp"
#include "extractQPK.hpp"

template<class Real, class DeviceType>
class BinaryAdvDiffFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  int dim_;

  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<MeshManager_adv_diff<Real,DeviceType>>  mesh_;
  ROL::Ptr<PDE_adv_diff<Real,DeviceType>>          pde_;
  ROL::Ptr<Linear_PDE_Constraint<Real,DeviceType>> con_;
  ROL::Ptr<Assembler<Real,DeviceType>>             assembler_;
  ROL::Ptr<ROL::Vector<Real>>                      u_, z_, p_;

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
    mesh_      = ROL::makePtr<MeshManager_adv_diff<Real,DeviceType>>(pl_);
    pde_       = ROL::makePtr<PDE_adv_diff<Real,DeviceType>>(pl_);
    con_       = ROL::makePtr<Linear_PDE_Constraint<Real,DeviceType>>(pde_,mesh_,comm_,pl_,*os_);
    assembler_ = con_->getAssembler();

    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr;
    u_ptr = assembler_->createStateVector();
    p_ptr = assembler_->createStateVector();
    r_ptr = assembler_->createResidualVector();
    ROL::Ptr<ROL::Vector<Real>> up, pp, rp, zp;
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real,DeviceType>>(u_ptr,pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real,DeviceType>>(p_ptr,pde_,assembler_,pl_);
    z_ = ROL::makePtr<ROL::StdVector<Real>>(dim_);
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    std::vector<ROL::Ptr<QoI<Real,DeviceType>>> qoi_vec(1,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<Real,DeviceType>>(pde_->getFE());
    Real stateCost = pl_.sublist("Problem").get("State Cost",4.0);
    std::vector<Real> wts = {stateCost};
    ROL::Ptr<ROL::Objective_SimOpt<Real>> obj
      = ROL::makePtr<PDE_Objective<Real,DeviceType>>(qoi_vec,wts,assembler_);
    bool storage = pl_.sublist("Problem").get("Use Storage",true);
    return ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj, con_, u_, z_, p_, storage, false);
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return z_;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    auto zlop = ROL::makePtr<ROL::StdVector<Real>>(dim_,0.0);
    auto zhip = ROL::makePtr<ROL::StdVector<Real>>(dim_,1.0);
    return ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
  }

  ROL::Ptr<ROL::Constraint<Real>> buildConstraint(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    Real budget(0);
    if (!useIneq) budget = pl_.sublist("Problem").get("Control Cost", 4.0);
    auto qoi = ROL::makePtr<QoI_Control_Cost_adv_diff<Real,DeviceType>>(budget);
    return ROL::makePtr<IntegralOptConstraint<Real,DeviceType>>(qoi,assembler_);
  }

  ROL::Ptr<ROL::Vector<Real>> buildMultiplier(void) {
    return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildSlackBoundConstraint(void) {
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      Real budget = pl_.sublist("Problem").get("Control Cost", 4.0);
      auto klop = ROL::makePtr<ROL::StdVector<Real>>(1,static_cast<Real>(0));
      auto khip = ROL::makePtr<ROL::StdVector<Real>>(1,budget);
      return ROL::makePtr<ROL::Bounds<Real>>(klop,khip);
    }
    return ROL::nullPtr;
  }
};

template<class Real, class DeviceType>
class BinaryAdvDiffQPFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<BinaryAdvDiffFactory<Real,DeviceType>> factory_;
  ROL::Ptr<extractQP<Real>> extract_;

public:
  BinaryAdvDiffQPFactory(ROL::ParameterList                 &pl,
                   const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                   const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    factory_ = ROL::makePtr<BinaryAdvDiffFactory<Real,DeviceType>>(pl_,comm_,os_);
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
