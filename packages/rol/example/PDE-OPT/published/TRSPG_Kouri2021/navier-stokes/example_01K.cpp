// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Solver.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_ConstraintFromObjective.hpp"

#include "../../../TOOLS/meshreaderK.hpp"
#include "../../../TOOLS/pdeconstraintK.hpp"
#include "../../../TOOLS/pdeobjectiveK.hpp"
#include "../../../TOOLS/pdevectorK.hpp"
#include "../../../TOOLS/integralconstraintK.hpp"

#include "pde_navier-stokesK.hpp"
#include "obj_navier-stokesK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  auto comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0))
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    bool output = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    output = (iprint > 0) && (myRank==0) && output;
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",output);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshReader<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    auto pde = ROL::makePtr<PDE_NavierStokes<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    auto u_ptr = assembler->createStateVector();    u_ptr->randomize();
    auto p_ptr = assembler->createStateVector();    p_ptr->randomize();
    auto z_ptr = assembler->createControlVector();  z_ptr->randomize();
    auto r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    auto xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize quadratic objective function.
    auto qoi = ROL::makePtr<QoI_Power_NavierStokes<RealT,DeviceT>>(*parlist,
                                                                   pde->getVelocityFE(),
                                                                   pde->getPressureFE(),
                                                                   pde->getControlFE(),
                                                                   pde->getStateFieldInfo(),
                                                                   pde->getControlFieldInfo());
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi,assembler);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    // Initialize volume constraint.
    RealT volFraction = parlist->sublist("Problem").get("Volume Fraction",0.5);
    auto qoi_vol = ROL::makePtr<QoI_Volume_NavierStokes<RealT,DeviceT>>(*parlist,
                                                                        pde->getControlFE(),
                                                                        pde->getControlFieldInfo());
    auto obj_vol = ROL::makePtr<IntegralOptObjective<RealT,DeviceT>>(qoi_vol,assembler);
    auto icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(obj_vol);
    auto imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0.0);

    // Build bound constraint
    auto lp = zp->clone(); lp->setScalar(0.0);
    auto hp = zp->clone(); hp->setScalar(1.0);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(lp, hp);
    // Build optimization problem
    auto optProb = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    optProb->addBoundConstraint(bnd);
    optProb->addLinearConstraint("Volume",icon,imul);
    optProb->setProjectionAlgorithm(*parlist);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      auto dup = up->clone(); dup->randomize();
      auto dzp = zp->clone(); dzp->randomize();
      con->checkApplyJacobian_1(*up,*zp,*dup,*up,true,*outStream);
      con->checkApplyJacobian_2(*up,*zp,*dzp,*up,true,*outStream);
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*up,true,*outStream);
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*up,true,*outStream);
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*zp,true,*outStream);
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*zp,true,*outStream);
      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);
      icon->checkApplyJacobian(*zp,*dzp,*imul,true,*outStream);
      icon->checkApplyAdjointHessian(*zp,*imul,*dzp,*zp,true,*outStream);
      robj->checkGradient(*zp,*dzp,true,*outStream);
      robj->checkHessVec(*zp,*dzp,true,*outStream);
      //optProb->check(*outStream);
    }

    // Solve optimization problem
    ROL::Ptr<ROL::Solver<RealT>> solver;
    std::vector<RealT> vals(6);

    // SPG
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    optProb->finalize(false,true,*outStream);
    parlist->sublist("Step").set("Type","Spectral Gradient");
    Teuchos::Time spgTimer("SPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    spgTimer.stop();
    *outStream << "Total optimization time = " << spgTimer.totalElapsedTime() << " seconds.\n";
    vals[0] = solver->getAlgorithmState()->value;

    // PQN
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    parlist->sublist("Step").set("Type","Line Search");
    parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Quasi-Newton Method");
    Teuchos::Time pqnTimer("PQN Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    pqnTimer.stop();
    *outStream << "Total optimization time = " << pqnTimer.totalElapsedTime() << " seconds.\n";
    vals[1] = solver->getAlgorithmState()->value;

    // LMTR
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time lmtrTimer("LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    lmtrTimer.stop();
    *outStream << "Total optimization time = " << lmtrTimer.totalElapsedTime() << " seconds.\n";
    vals[2] = solver->getAlgorithmState()->value;

    // TRSPG
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time trspgTimer("TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    trspgTimer.stop();
    *outStream << "Total optimization time = " << trspgTimer.totalElapsedTime() << " seconds.\n";
    vals[3] = solver->getAlgorithmState()->value;

    // AL-LMTR
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    optProb->edit();
    optProb->finalize(true,true,*outStream);
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time allmtrTimer("AL-LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    allmtrTimer.stop();
    *outStream << "Total optimization time = " << allmtrTimer.totalElapsedTime() << " seconds.\n";
    vals[4] = solver->getAlgorithmState()->value;

    // AL-LMTR
    up->zero(); pp->zero(); imul->zero();
    zp->setScalar(volFraction);
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time altrspgTimer("AL-TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(optProb,*parlist);
    solver->solve(*outStream);
    altrspgTimer.stop();
    *outStream << "Total optimization time = " << altrspgTimer.totalElapsedTime() << " seconds.\n";
    vals[5] = solver->getAlgorithmState()->value;

    RealT minval = std::min(vals[0],std::min(vals[1],std::min(vals[2],vals[3])));
    *outStream << std::scientific << std::setprecision(16);
    *outStream << "SPG      " << (vals[0]-minval)/minval << std::endl;
    *outStream << "PQN      " << (vals[1]-minval)/minval << std::endl;
    *outStream << "LMTR     " << (vals[2]-minval)/minval << std::endl;
    *outStream << "TRSPG    " << (vals[3]-minval)/minval << std::endl;
    *outStream << "AL-LMTR  " << (vals[4]-minval)/minval << std::endl;
    *outStream << "AL-TRSPG " << (vals[5]-minval)/minval << std::endl;
    *outStream << std::endl;

    // Output
    con->outputTpetraVector(z_ptr,"control.txt");
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //con->outputTpetraData();

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
