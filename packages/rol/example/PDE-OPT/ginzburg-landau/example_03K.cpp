// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the stuctural topology optimization problem.
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

#include "pde_ginzburg-landau_ex02K.hpp"
#include "obj_ginzburg-landau_ex02K.hpp"

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
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex03.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    parlist->sublist("Problem").set("Current Loading",static_cast<RealT>(1));
    parlist->sublist("Problem").set("State Scaling",  static_cast<RealT>(2e-3));
    parlist->sublist("Problem").set("Control Scaling",static_cast<RealT>(1e-2));

    parlist->sublist("Geometry").set("Width",2.0);
    int NY = parlist->sublist("Geometry").get("NY",16);
    parlist->sublist("Geometry").set("NX",2*NY);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    auto pde = ROL::makePtr<PDE_GinzburgLandau_ex02<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    // Create vector.
    auto u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    auto p_ptr = assembler->createStateVector();    p_ptr->putScalar(0.0);
    auto z_ptr = assembler->createControlVector();  z_ptr->putScalar(0.0);
    auto r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT>> x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize compliance objective function.
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_GinzburgLandau_StateTracking_ex02<RealT,DeviceT>>(pde->getFE(),
                                                                            pde->getFieldHelper(),
                                                                            *parlist);
    qoi_vec[1] = ROL::makePtr<QoI_GinzburgLandau_ControlPenalty<RealT,DeviceT>>(pde->getFE(),
                                                                        pde->getBdryFE(),
                                                                        pde->getBdryCellLocIds(),
                                                                        pde->getFieldHelper(),
                                                                        *parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,assembler);
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Output uncontrolled state.
    up->zero(); zp->zero();
    con->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state_uncontrolled.txt");

    // Set up optimization problem
    ROL::Ptr<ROL::Problem<RealT>> problem;
    bool useReducedSpace = parlist->sublist("Problem").get("Use Reduced Space",true);
    if (useReducedSpace) {
      problem = ROL::makePtr<ROL::Problem<RealT>>(robj,zp);
    }
    else {
      problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
      problem->addConstraint("PDE",con,pp);
    }
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if (checkDeriv) problem->check(true,*outStream);

    // Solve optimization problem
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
