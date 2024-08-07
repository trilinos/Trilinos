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

#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "objsum.hpp"
#include "tv_2d.hpp"

template<class Real>
class BinaryAdvDiffFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  int dim_;

  mutable ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<FEMdata<Real>>              fem_;
  ROL::Ptr<Assembler<Real>>            assembler_;
  ROL::Ptr<ROL::Objective<Real>>       penalty_, binary_;
  ROL::Ptr<ROL::Vector<Real>>          z_;
  ROL::Ptr<ROL::BoundConstraint<Real>> bnd_;  

public:
  BinaryAdvDiffFactory(ROL::ParameterList                 &pl,
                 const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                 const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    fem_        = ROL::makePtr<FEMdata<Real>>(comm,pl_,*os_);
    assembler_  = fem_->getAssembler();

    std::string costType = pl_.sublist("Problem").get("Control Cost Type", "TV");
    std::vector<ROL::Ptr<QoI<Real>>> qoi_pen(2,ROL::nullPtr);
    if (costType=="TV")
      qoi_pen[0] = ROL::makePtr<QoI_TVControl_Cost_adv_diff<Real>>(fem_->getControlFE(),pl_);
    else if (costType=="L1")
      qoi_pen[0] = ROL::makePtr<QoI_Control_Cost_adv_diff<Real>>(fem_->getControlFE(),pl_);
    else
      qoi_pen[0] = ROL::makePtr<QoI_Control_Cost_L2_adv_diff<Real>>(fem_->getControlFE(),pl_);
    qoi_pen[1] = ROL::makePtr<QoI_IntegralityControl_Cost_adv_diff<Real>>(fem_->getControlFE(),pl_);
    int ctrlOrder = pl_.sublist("Problem").get("Order of Control Discretization",0);
    if (ctrlOrder==0 && costType=="TV") {
      penalty_ = ROL::makePtr<Objective_TV_2D_C0<Real>>(pl_);
    }
    else {
      penalty_ = ROL::makePtr<IntegralOptObjective<Real>>(qoi_pen[0],assembler_);
    }
    binary_  = ROL::makePtr<IntegralOptObjective<Real>>(qoi_pen[1],assembler_);
    // Create template control vector
    bool usePC = pl_.sublist("Problem").get("Piecewise Constant Controls", true);
    if (usePC) {
      int nx = pl_.sublist("Problem").get("Number of X-Cells", 4);
      int ny = pl_.sublist("Problem").get("Number of Y-Cells", 2);
      Real width  = pl_.sublist("Geometry").get("Width",  2.0);
      Real height = pl_.sublist("Geometry").get("Height", 1.0);
      Real vol = width/static_cast<Real>(nx) * height/static_cast<Real>(ny);
      ROL::Ptr<std::vector<Real>> xvec, svec;
      xvec = ROL::makePtr<std::vector<Real>>(nx*ny,0);
      svec = ROL::makePtr<std::vector<Real>>(nx*ny,vol);
      std::ifstream file;
      std::stringstream name;
      name << "control_" << nx << "_" << ny << ".txt";
      file.open(name.str());
      if (file.is_open()) {
        for (int i = 0; i < nx; ++i) {
          for (int j = 0; j < ny; ++j) {
            file >> (*xvec)[i+j*nx];
          }
        }
      }
      file.close();
      z_ = ROL::makePtr<ROL::PrimalScaledStdVector<Real>>(xvec,svec);
    }
    else {
      z_ = fem_->createControlVector(pl_);
    }
    // Create bound constraint
    ROL::Ptr<ROL::Vector<Real>> zlop = z_->clone(), zhip = z_->clone();
    zlop->setScalar(static_cast<Real>(0));
    zhip->setScalar(static_cast<Real>(1));
    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    int bbout = pl_.sublist("Problem").get("BB Output Level", 0);
    ROL::Ptr<ROL::Objective<Real>>
      obj = ROL::makePtr<Sum_Objective<Real>>(fem_,penalty_,binary_,pl_,os_,bbout>2);
    ROL::Ptr<ROL::Vector<Real>> z = z_->clone(); z->set(*z_);
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj,z);
    problem->addBoundConstraint(bnd_);
    return problem;
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    ROL::Ptr<ROL::Vector<Real>> z = z_->clone(); z->set(*z_);
    return z;
  }

  bool controlType(void) const {
    bool usePC = pl_.sublist("Problem").get("Piecewise Constant Controls", true);
    return !usePC;
  }

  void getState(ROL::Ptr<ROL::Vector<Real>> &u, const ROL::Ptr<ROL::Vector<Real>> &z) const {
    u = fem_->createStateVector(pl_);
    Misfit_Objective<Real> obj(fem_,pl_);
    obj.solvePDE(*u,*z);
  }

  ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  void print(std::ostream &stream = std::cout) {
    assembler_->printMeshData(stream);
  }

  void check(std::ostream &stream = std::cout) {
    ROL::Ptr<ROL::Vector<Real>> z1 = z_->clone(); z1->randomize();
    ROL::Ptr<ROL::Vector<Real>> z2 = z_->clone(); z2->randomize();
    penalty_->checkGradient(*z1,*z2,true,stream);
    penalty_->checkHessVec(*z1,*z2,true,stream);
    binary_->checkGradient(*z1,*z2,true,stream);
    binary_->checkHessVec(*z1,*z2,true,stream);
    Misfit_Objective<Real> obj(fem_,pl_);
    obj.checkGradient(*z1,*z2,true,stream);
    obj.checkHessVec(*z1,*z2,true,stream);
  }
};

#endif
