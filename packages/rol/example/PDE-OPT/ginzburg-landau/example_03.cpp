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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/meshmanager.hpp"

#include "pde_ginzburg-landau_ex02.hpp"
#include "obj_ginzburg-landau_ex02.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex03.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    parlist->sublist("Problem").set("Current Loading",static_cast<RealT>(1));
    parlist->sublist("Problem").set("State Scaling",  static_cast<RealT>(2e-3));
    parlist->sublist("Problem").set("Control Scaling",static_cast<RealT>(1e-2));

    parlist->sublist("Geometry").set("Width",2.0);
    int NY = parlist->sublist("Geometry").get("NY",16);
    parlist->sublist("Geometry").set("NX",2*NY);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_GinzburgLandau<RealT>>
      pde = ROL::makePtr<PDE_GinzburgLandau_ex02<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>>
      con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT>>
      pdecon = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    ROL::Ptr<Assembler<RealT>> assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create vector.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, z_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();    p_ptr->putScalar(0.0);
    z_ptr = assembler->createControlVector();  z_ptr->putScalar(0.0);
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, zp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT>> x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize compliance objective function.
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_GinzburgLandau_StateTracking_ex02<RealT>>(pde->getFE(),
                                                                            pde->getFieldHelper(),
                                                                            *parlist);
    qoi_vec[1] = ROL::makePtr<QoI_GinzburgLandau_ControlPenalty<RealT>>(pde->getFE(),
                                                                        pde->getBdryFE(),
                                                                        pde->getBdryCellLocIds(),
                                                                        pde->getFieldHelper(),
                                                                        *parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>
      obj = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,assembler);
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>>
      robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Output uncontrolled state.
    up->zero(); zp->zero();
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state_uncontrolled.txt");

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
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");

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
