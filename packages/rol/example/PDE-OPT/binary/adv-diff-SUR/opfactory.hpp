// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef OPFACTORY_BINARY_ADVDIFF_SUR_HPP
#define OPFACTORY_BINARY_ADVDIFF_SUR_HPP

#include "ROL_Bounds.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "pde_adv_diff.hpp"
#include "qoi_adv_diff.hpp"
#include "mesh_adv_diff.hpp"

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
  ROL::Ptr<ROL::Objective_SimOpt<Real>> obj_;
  ROL::Ptr<ROL::Objective<Real>>        robj_;
  ROL::Ptr<Assembler<Real>>             assembler_;
  ROL::Ptr<ROL::Vector<Real>>           u_, z_, p_;

public:
  BinaryAdvDiffFactory(ROL::ParameterList                 &pl,
                 const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                 const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {}

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    update();
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(buildObjective(),buildSolutionVector());
    problem->addBoundConstraint(buildBoundConstraint());
    return problem;
  }

  void update(void) {
    mesh_      = ROL::makePtr<MeshManager_adv_diff<Real>>(pl_);
    pde_       = ROL::makePtr<PDE_adv_diff<Real>>(pl_);
    con_       = ROL::makePtr<Linear_PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
    assembler_ = con_->getAssembler();

    bool usePC = pl_.sublist("Problem").get("Piecewise Constant Controls", true);
    int  order = pl_.sublist("Problem").get("Hilbert Curve Order", 2);
    int      n = std::pow(2,order);
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr, z_ptr;
    u_ptr = assembler_->createStateVector();
    p_ptr = assembler_->createStateVector();
    ROL::Ptr<ROL::Vector<Real>> up, pp, zp;
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(u_ptr,pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(p_ptr,pde_,assembler_,pl_);
    if (usePC) {
      Real XL  = pl_.sublist("Geometry").get("X0", 0.0);
      Real YL  = pl_.sublist("Geometry").get("Y0", 0.0);
      Real XU  = XL + pl_.sublist("Geometry").get("Width",  1.0);
      Real YU  = YL + pl_.sublist("Geometry").get("Height", 1.0);
      Real vol = (XU-XL)/static_cast<Real>(n) * (YU-YL)/static_cast<Real>(n);
      ROL::Ptr<std::vector<Real>> xvec, svec;
      xvec = ROL::makePtr<std::vector<Real>>(n*n,0);
      svec = ROL::makePtr<std::vector<Real>>(n*n,vol);
      ROL::Ptr<ROL::StdVector<Real>> xstd = ROL::makePtr<ROL::PrimalScaledStdVector<Real>>(xvec,svec);
      z_ = ROL::makePtr<PDE_OptVector<Real>>(xstd);
    }
    else {
      z_ptr = assembler_->createControlVector();
      z_ = ROL::makePtr<PDE_PrimalOptVector<Real>>(z_ptr,pde_,assembler_,pl_);
    }

    bool useL2       = pl_.sublist("Problem").get("Use L2", false);
    bool storage     = pl_.sublist("Problem").get("Use Storage", true);
    Real stateCost   = pl_.sublist("Problem").get("State Cost",   1.0);
    Real controlCost = pl_.sublist("Problem").get("Control Cost", 1.0);
    std::vector<Real> wts = {stateCost, controlCost};
    std::vector<ROL::Ptr<QoI<Real>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<Real>>(pde_->getFE());
    if (useL2) {
      qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_L2_adv_diff<Real>>(pde_->getFE(),pl_);
    }
    else {
      qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_adv_diff<Real>>(pde_->getFE(),pl_);
    }
    obj_  = ROL::makePtr<PDE_Objective<Real>>(qoi_vec, wts, assembler_);
    robj_ = ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj_, con_, u_, z_, p_, storage, false);
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    return robj_;
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return z_;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    ROL::Ptr<ROL::Vector<Real>> zlop = z_->clone(), zhip = z_->clone();
    zlop->setScalar(static_cast<Real>(0));
    zhip->setScalar(static_cast<Real>(1));
    return ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
    //return ROL::nullPtr;
  }

  void getState(ROL::Ptr<ROL::Vector<Real>> &u, const ROL::Ptr<ROL::Vector<Real>> &z) const {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    u = u_->clone();
    ROL::Ptr<ROL::Vector<Real>> r = u_->dual().clone();
    con_->solve(*r, *u, *z, tol);
  }

  ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  void print(std::ostream &stream = std::cout) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    ROL::Ptr<ROL::Vector<Real>> r = u_->dual().clone();
    con_->solve(*r,*u_,*z_,tol);
    assembler_->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(u_)->getVector(),"state.txt");
    bool usePC = pl_.sublist("Problem").get("Piecewise Constant Controls", true);
    if (usePC) {
      std::ofstream zfile;
      zfile.open("control.txt");
      ROL::Ptr<std::vector<Real>> zvec
        = ROL::dynamicPtrCast<PDE_OptVector<Real>>(z_)->getParameter()->getVector();
      for (unsigned i = 0; i < zvec->size(); ++i) {
        zfile << (*zvec)[i] << std::endl;
      }
      zfile.close();
      pde_->print();
    }
    else {
      assembler_->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(z_)->getVector(),"control.txt");
    }
    assembler_->printMeshData(stream);
  }

  void check(std::ostream &stream = std::cout) {
    update();
    ROL::Ptr<ROL::Vector<Real>> r1 = u_->dual().clone(); r1->randomize();
    ROL::Ptr<ROL::Vector<Real>> u1 = u_->clone();        u1->randomize();
    ROL::Ptr<ROL::Vector<Real>> u2 = u_->clone();        u2->randomize();
    ROL::Ptr<ROL::Vector<Real>> z1 = z_->clone();        z1->randomize();
    ROL::Ptr<ROL::Vector<Real>> z2 = z_->clone();        z2->randomize();
    con_->checkSolve(*u1,*z1,*r1,true,stream);
    con_->checkApplyJacobian_1(*u1,*z1,*u2,*r1,true,stream);
    con_->checkApplyJacobian_2(*u1,*z1,*z2,*r1,true,stream);
    con_->checkApplyAdjointHessian_11(*u1,*z1,*u1,*u2,*r1,true,stream);
    con_->checkApplyAdjointHessian_12(*u1,*z1,*u1,*u2,*z1,true,stream);
    con_->checkApplyAdjointHessian_21(*u1,*z1,*u1,*z2,*r1,true,stream);
    con_->checkApplyAdjointHessian_22(*u1,*z1,*u1,*z2,*z1,true,stream);
    obj_->checkGradient_1(*u1,*z1,*u2,true,stream);
    obj_->checkGradient_2(*u1,*z1,*z2,true,stream);
    obj_->checkHessVec_11(*u1,*z1,*u2,true,stream);
    obj_->checkHessVec_12(*u1,*z1,*z2,true,stream);
    obj_->checkHessVec_21(*u1,*z1,*u2,true,stream);
    obj_->checkHessVec_22(*u1,*z1,*z2,true,stream);
    robj_->checkGradient(*z1,*z2,true,stream);
    robj_->checkHessVec(*z1,*z2,true,stream);
  }
};

#endif
