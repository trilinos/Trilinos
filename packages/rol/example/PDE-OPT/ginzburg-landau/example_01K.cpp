// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the optimal control of Ginzburg-Landau problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "../TOOLS/meshmanagerK.hpp"

#include "pde_ginzburg-landau_ex01K.hpp"
#include "obj_ginzburg-landau_ex01K.hpp"

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
    RealT tol(1e-8);// one(1);

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    parlist->sublist("Problem").set("Current Loading",static_cast<RealT>(1));
    parlist->sublist("Problem").set("State Scaling",  static_cast<RealT>(2e-6));
    parlist->sublist("Problem").set("Control Scaling",static_cast<RealT>(1e-1));

    parlist->sublist("Geometry").set("Width",1.0);
    int NX = parlist->sublist("Geometry").get("NX",16);
    parlist->sublist("Geometry").set("NY",NX);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    auto pde = ROL::makePtr<PDE_GinzburgLandau_ex01<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    // Create vectors.
    auto u_ptr = assembler->createStateVector();    u_ptr->randomize();
    auto p_ptr = assembler->createStateVector();    p_ptr->randomize();
    auto z_ptr = assembler->createControlVector();  z_ptr->randomize();
    auto r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);

    // Initialize compliance objective function.
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_GinzburgLandau_StateTracking_ex01<RealT,DeviceT>>(pde->getFE(),
                                                                            pde->getFieldHelper(),
                                                                            *parlist);
    qoi_vec[1] = ROL::makePtr<QoI_GinzburgLandau_ControlPenalty<RealT,DeviceT>>(pde->getFE(),
                                                                        pde->getBdryFE(),
                                                                        pde->getBdryCellLocIds(),
                                                                        pde->getFieldHelper(),
                                                                        *parlist);
    auto obj  = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,assembler);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      auto dup = up->clone(); dup->randomize(-1.0,1.0);
      auto dzp = zp->clone(); dzp->randomize(-1.0,1.0);
      auto rzp = zp->clone(); rzp->randomize(-1.0,1.0);
      auto dualup = up->dual().clone();
      auto dualzp = zp->dual().clone();
      ROL::Vector_SimOpt<RealT> x(up,zp), d(dup,dzp);

      *outStream << "\n\nCheck Opt Vector\n";
      zp->checkVector(*dzp,*rzp,true,*outStream);

      std::vector<ROL::Ptr<ROL::Objective_SimOpt<RealT>>> obj_vec(2,ROL::nullPtr);
      obj_vec[0] = ROL::makePtr<IntegralObjective<RealT,DeviceT>>(qoi_vec[0],assembler);
      obj_vec[1] = ROL::makePtr<IntegralObjective<RealT,DeviceT>>(qoi_vec[1],assembler);

      *outStream << "\n\nCheck Gradient of State Objective Function\n";
      obj_vec[0]->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of State Objective Function\n";
      obj_vec[0]->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of State Objective Function\n";
      obj_vec[0]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of State Objective Function\n";
      obj_vec[0]->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of State Objective Function\n";
      obj_vec[0]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj_vec[0]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj_vec[0]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of State Objective Function\n";
      obj_vec[0]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Control Objective Function\n";
      obj_vec[1]->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of Control Objective Function\n";
      obj_vec[1]->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Control Objective Function\n";
      obj_vec[1]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Control Objective Function\n";
      obj_vec[1]->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Control Objective Function\n";
      obj_vec[1]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj_vec[1]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj_vec[1]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Control Objective Function\n";
      obj_vec[1]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of Full Objective Function\n";
      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Full Objective Function\n";
      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Full Objective Function\n";
      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Full Objective Function\n";
      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(x,*pp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of PDE Constraint\n";
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of PDE Constraint\n";
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of PDE Constraint\n";
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dualzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of PDE Constraint\n";
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dualzp,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

      zp->zero();
      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    // Output uncontrolled state.
    up->zero(); zp->zero();
    con->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state_uncontrolled.txt");
    z_ptr->putScalar(static_cast<RealT>(1));

    auto problem = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    problem->finalize(false,true,*outStream);
    ROL::Solver<RealT> solver(problem,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    con->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state.txt");
    con->outputTpetraVector(z_ptr,"control.txt");

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
