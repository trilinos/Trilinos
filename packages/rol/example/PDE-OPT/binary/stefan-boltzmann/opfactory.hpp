// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef OPFACTORY_BINARY_STEFAN_BOLTZMANN_HPP
#define OPFACTORY_BINARY_STEFAN_BOLTZMANN_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
//#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "pde_stefan_boltzmann.hpp"
#include "obj_stefan_boltzmann.hpp"
#include "con_stefan_boltzmann.hpp"

template<class Real>
class BinaryStefanBoltzmannFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  mutable ROL::ParameterList   pl_;
  ROL::Ptr<Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream>       os_;

  ROL::Ptr<MeshManager<Real>>              mesh_;
  ROL::Ptr<BinaryStefanBoltzmannPDE<Real>> pde_;
  ROL::Ptr<QoI_StateCost<Real>>            qoi_;
  ROL::Ptr<Assembler<Real>>                assembler_;
  ROL::Ptr<ROL::Objective_SimOpt<Real>>    obj_;
  ROL::Ptr<ROL::Constraint_SimOpt<Real>>   con_;
  ROL::Ptr<ROL::Vector<Real>>              u_, p_, z_;
  ROL::Ptr<ROL::BoundConstraint<Real>>     bnd_;  
  ROL::Ptr<BudgetConstraint<Real>>         budget_con_;
  ROL::Ptr<ROL::Vector<Real>>              budget_mul_;
  ROL::Ptr<ROL::BoundConstraint<Real>>     budget_bnd_;

  bool useParam_, useBudget_;

public:
  BinaryStefanBoltzmannFactory(ROL::ParameterList                 &pl,
                         const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), os_(os) {
#ifndef HAVE_MPI
    comm_ = ROL::makePtr<Teuchos::SerialComm<int>>();
    initialize();
#endif
  }

#ifdef HAVE_MPI
  void setCommunicator(const ROL::Ptr<const MPI_Comm> &comm) {
    comm_ = ROL::makePtr<Teuchos::MpiComm<int>>(*comm);
    initialize();
  }
#endif
//  BinaryStefanBoltzmannFactory(ROL::ParameterList                 &pl,
//                         const ROL::Ptr<const Teuchos::Comm<int>> &comm,
//                         const ROL::Ptr<std::ostream>             &os)
//    : pl_(pl), comm_(comm), os_(os) {
//    // Create PDE Constraint
//    int probDim = pl_.sublist("Problem").get("Problem Dimension",2);
//    int nProcs  = comm_->getSize();
//    TEUCHOS_TEST_FOR_EXCEPTION(probDim<2||probDim>3, std::invalid_argument,
//      ">>> PDE-OPT/binary/stefan-boltzmann/example_01.cpp: Problem dim is not 2 or 3!");
//    //if (probDim == 2)      mesh_ = ROL::makePtr<MeshManager_Rectangle<Real>>(pl_);
//    //else if (probDim == 3)
//    mesh_ = ROL::makePtr<MeshReader<Real>>(pl_,nProcs);
//    pde_ = ROL::makePtr<BinaryStefanBoltzmannPDE<Real>>(pl_);
//    con_ = ROL::makePtr<PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
//    con_->setSolveParameters(pl_);
//    assembler_ = ROL::dynamicPtrCast<PDE_Constraint<Real>>(con_)->getAssembler();
//    // Create template vectors
//    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,pl_);
//    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,pl_);
//    useParam_ = pl.sublist("Problem").get("Use Parametric Control", true);
//    bool init = pl_.sublist("Problem").get("Input Control",false);
//    std::string inCtrlName = pl_.sublist("Problem").get("Input Control Name","control.txt");
//    if (useParam_) {
//      int nx = pl.sublist("Problem").get("Number X Control Patches", 4);
//      int ny = pl.sublist("Problem").get("Number Y Control Patches", 4);
//      z_ = ROL::makePtr<ROL::StdVector<Real>>(nx*ny);
//      if (init) {
//        ROL::Ptr<std::vector<Real>> zdata = ROL::dynamicPtrCast<ROL::StdVector<Real>>(z_)->getVector();
//        std::ifstream file;
//        file.open(inCtrlName);
//        for (int i = 0; i < ny; ++i) {
//          for (int j = 0; j < nx; ++j) {
//            file >> (*zdata)[j+i*nx];
//          }
//        }
//        file.close();
//      }
//    }
//    else {
//      z_ = ROL::makePtr<PDE_PrimalOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,pl_);
//      if (init) {
//        ROL::Ptr<Tpetra::MultiVector<>> zdata = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(z_)->getVector();
//        assembler_->inputTpetraVector(zdata, inCtrlName);
//      }
//    }
//    if (!init) z_->setScalar(static_cast<Real>(1));
//    u_->setScalar(static_cast<Real>(1));
//    // Create objective function
//    qoi_ = ROL::makePtr<QoI_StateCost<Real>>(pde_->getVolFE(),pl_);
//    obj_ = ROL::makePtr<PDE_Objective<Real>>(qoi_,assembler_);
//    // Create bound constraint
//    ROL::Ptr<ROL::Vector<Real>> zlop = z_->clone(), zhip = z_->clone();
//    zlop->setScalar(static_cast<Real>(0));
//    zhip->setScalar(static_cast<Real>(1));
//    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
//    // Create budget constraint
//    useBudget_ = pl.sublist("Problem").get("Use Budget Constraint",true);
//    if (useParam_ && useBudget_) {
//      budget_con_ = ROL::makePtr<BudgetConstraint<Real>>(pl);
//      budget_mul_ = budget_con_->createMultiplier();
//      budget_bnd_ = budget_con_->createBounds();
//    }
//  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    ROL::Ptr<ROL::Vector<Real>> z = z_->clone(); z->set(*z_);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<Real>>
      obj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj_, con_, u_, z_, p_, true, false);
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj,z);
    problem->addBoundConstraint(bnd_);
    if (useParam_ && useBudget_) {
      problem->addLinearConstraint("Budget",budget_con_,budget_mul_,budget_bnd_);
      problem->setProjectionAlgorithm(pl_);
    }
    return problem;
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    ROL::Ptr<ROL::Vector<Real>> z = z_->clone(); z->set(*z_);
    return z;
  }

  void getState(ROL::Ptr<ROL::Vector<Real>> &u, const ROL::Ptr<ROL::Vector<Real>> &z) const {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    u = u_->clone();
    ROL::Ptr<ROL::Vector<Real>> r = u->dual().clone();
    if (useParam_ && useBudget_)
      con_->solve(*r,*u,*ROL::dynamicPtrCast<ROL::PartitionedVector<Real>>(z)->get(0),tol);
    else
      con_->solve(*r,*u,*z,tol);
  }

  ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  void print(std::ostream &stream = std::cout) {
    assembler_->printMeshData(stream);
  }

  void printTpetraVector(const ROL::Ptr<ROL::Vector<Real>> &x, std::string name) const {
    assembler_->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(x)->getVector(),name);
  }

  void printControl(const ROL::Ptr<ROL::Vector<Real>> &x, std::string name) const {
    if (useParam_) {
      ROL::Ptr<std::vector<Real>> data;
      if (useBudget_)
        data = ROL::dynamicPtrCast<ROL::StdVector<Real>>(ROL::dynamicPtrCast<ROL::PartitionedVector<Real>>(x)->get(0))->getVector();
      else
        data = ROL::dynamicPtrCast<ROL::StdVector<Real>>(x)->getVector();  
      int dim = data->size();
      std::ofstream file; file.open(name);
      for (int i = 0; i < dim; ++i) {
        file << std::scientific << std::setprecision(16);
        file << (*data)[i] << std::endl;
      }
      file.close();
    }
    else {
      printTpetraVector(x,name);
    }
  }

  void check(std::ostream &stream = std::cout) {
    ROL::Ptr<ROL::Vector<Real>> z1 = z_->clone(); z1->randomize();
    ROL::Ptr<ROL::Vector<Real>> z2 = z_->clone(); z2->randomize();
    ROL::Ptr<ROL::Vector<Real>> u1 = u_->clone(); u1->randomize(500.0,2000.0);
    ROL::Ptr<ROL::Vector<Real>> u2 = u_->clone(); u2->randomize();
    ROL::Ptr<ROL::Vector<Real>> p  = p_->clone(); p->randomize();
    ROL::Ptr<ROL::Vector<Real>> du = u_->dual().clone(); du->randomize();
    ROL::Ptr<ROL::Vector<Real>> dz = z_->dual().clone(); dz->randomize();
    con_->checkSolve(*u1,*z1,*du,true,stream);
    con_->checkAdjointConsistencyJacobian_1(*p,*u2,*u1,*z1,true,stream);
    con_->checkAdjointConsistencyJacobian_2(*p,*z2,*u1,*z1,true,stream);
    con_->checkInverseJacobian_1(*du,*u2,*u1,*z1,true,stream);
    con_->checkInverseAdjointJacobian_1(*p,*du,*u1,*z1,true,stream);
    con_->checkApplyJacobian_1(*u1,*z1,*u2,*du,true,stream);
    con_->checkApplyJacobian_2(*u1,*z1,*z2,*du,true,stream);
    con_->checkApplyAdjointHessian_11(*u1,*z1,*p,*u2,*du,true,stream);
    con_->checkApplyAdjointHessian_21(*u1,*z1,*p,*z2,*du,true,stream);
    con_->checkApplyAdjointHessian_12(*u1,*z1,*p,*u2,*dz,true,stream);
    con_->checkApplyAdjointHessian_22(*u1,*z1,*p,*z2,*dz,true,stream);
    obj_->checkGradient_1(*u1,*z1,*u2,true,stream);
    obj_->checkGradient_2(*u1,*z1,*z2,true,stream);
    obj_->checkHessVec_11(*u1,*z1,*u2,true,stream);
    obj_->checkHessVec_21(*u1,*z1,*u2,true,stream);
    obj_->checkHessVec_12(*u1,*z1,*z2,true,stream);
    obj_->checkHessVec_22(*u1,*z1,*z2,true,stream);
  }

private:

  void initialize(void) {
    // Create PDE Constraint
    int probDim = pl_.sublist("Problem").get("Problem Dimension",2);
    int nProcs  = comm_->getSize();
    TEUCHOS_TEST_FOR_EXCEPTION(probDim<2||probDim>3, std::invalid_argument,
      ">>> PDE-OPT/binary/stefan-boltzmann/example_01.cpp: Problem dim is not 2 or 3!");
    //if (probDim == 2)      mesh_ = ROL::makePtr<MeshManager_Rectangle<Real>>(pl_);
    //else if (probDim == 3)
    mesh_ = ROL::makePtr<MeshReader<Real>>(pl_,nProcs);
    pde_ = ROL::makePtr<BinaryStefanBoltzmannPDE<Real>>(pl_);
    con_ = ROL::makePtr<PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
    con_->setSolveParameters(pl_);
    assembler_ = ROL::dynamicPtrCast<PDE_Constraint<Real>>(con_)->getAssembler();
    // Create template vectors
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,pl_);
    useParam_ = pl_.sublist("Problem").get("Use Parametric Control", true);
    bool init = pl_.sublist("Problem").get("Input Control",false);
    std::string inCtrlName = pl_.sublist("Problem").get("Input Control Name","control.txt");
    if (useParam_) {
      int nx = pl_.sublist("Problem").get("Number X Control Patches", 4);
      int ny = pl_.sublist("Problem").get("Number Y Control Patches", 4);
      z_ = ROL::makePtr<ROL::StdVector<Real>>(nx*ny);
      if (init) {
        ROL::Ptr<std::vector<Real>> zdata = ROL::dynamicPtrCast<ROL::StdVector<Real>>(z_)->getVector();
        std::ifstream file;
        file.open(inCtrlName);
        for (int i = 0; i < ny; ++i) {
          for (int j = 0; j < nx; ++j) {
            file >> (*zdata)[j+i*nx];
          }
        }
        file.close();
      }
    }
    else {
      z_ = ROL::makePtr<PDE_PrimalOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,pl_);
      if (init) {
        ROL::Ptr<Tpetra::MultiVector<>> zdata = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(z_)->getVector();
        assembler_->inputTpetraVector(zdata, inCtrlName);
      }
    }
    if (!init) z_->setScalar(static_cast<Real>(1));
    u_->setScalar(static_cast<Real>(1));
    // Create objective function
    qoi_ = ROL::makePtr<QoI_StateCost<Real>>(pde_->getVolFE(),pl_);
    obj_ = ROL::makePtr<PDE_Objective<Real>>(qoi_,assembler_);
    // Create bound constraint
    ROL::Ptr<ROL::Vector<Real>> zlop = z_->clone(), zhip = z_->clone();
    zlop->setScalar(static_cast<Real>(0));
    zhip->setScalar(static_cast<Real>(1));
    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
    // Create budget constraint
    useBudget_ = pl_.sublist("Problem").get("Use Budget Constraint",true);
    if (useParam_ && useBudget_) {
      budget_con_ = ROL::makePtr<BudgetConstraint<Real>>(pl_);
      budget_mul_ = budget_con_->createMultiplier();
      budget_bnd_ = budget_con_->createBounds();
    }
  }
};

#endif
