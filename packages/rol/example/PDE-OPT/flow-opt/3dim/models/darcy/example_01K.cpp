// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Darcy porosity optimization problem.
*/

#include "Teuchos_Comm.hpp"
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

#include "../../../../TOOLS/meshreaderK.hpp"
#include "../../../../TOOLS/pdeconstraintK.hpp"
#include "../../../../TOOLS/pdeobjectiveK.hpp"
#include "../../../../TOOLS/pdevectorK.hpp"
#include "../../../../TOOLS/integralconstraintK.hpp"

#include "pde_darcyK.hpp"
#include "obj_darcyK.hpp"

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
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
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
    // Initialize PDE describing Darcy equations.
    auto pde = ROL::makePtr<PDE_Darcy<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vectors
    auto u_ptr = assembler->createStateVector(); u_ptr->randomize();
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto p_ptr = assembler->createStateVector(); p_ptr->randomize();
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);

    // Create control vectors
    bool useParamVar = parlist->sublist("Problem").get("Use Optimal Constant Velocity",false);
    int dim = 2;
    auto z_ptr = assembler->createControlVector(); z_ptr->randomize();
    ROL::Ptr<std::vector<RealT>>            z0_ptr;
    ROL::Ptr<ROL::StdVector<RealT>>         z0p;
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> z1p;
    ROL::Ptr<ROL::Vector<RealT>>            zp;
    if (useParamVar) {
      z0_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      z0p = ROL::makePtr<ROL::StdVector<RealT>>(z0_ptr);
      z1p = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
      zp = ROL::makePtr<PDE_OptVector<RealT>>(z1p,z0p,myRank);
    }
    else {
      zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    }

    // Create combined vector
    auto xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize quadratic objective function.
    auto qoi  = ROL::makePtr<QoI_VelocityTracking_Darcy<RealT,DeviceT>>(*parlist,
                                                                        pde->getPressureFE(),
                                                                        pde->getControlFE(),
                                                                        pde->getPermeability());
    auto obj  = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi,assembler);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    // Build bound constraint
    ROL::Ptr<ROL::Vector<RealT>> lp, hp;
    if (useParamVar) {
      auto l0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_NINF<RealT>());
      auto h0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_INF<RealT>());
      auto l1_ptr = assembler->createControlVector();
      auto h1_ptr = assembler->createControlVector();
      auto l1p = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(l1_ptr,pde,assembler,*parlist);
      auto h1p = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(h1_ptr,pde,assembler,*parlist);
      l1p->setScalar(0.0);
      h1p->setScalar(1.0);
      lp = ROL::makePtr<PDE_OptVector<RealT>>(l1p,l0p,myRank);
      hp = ROL::makePtr<PDE_OptVector<RealT>>(h1p,h0p,myRank);
    }
    else {
      lp = zp->clone(); lp->setScalar(0.0);
      hp = zp->clone(); hp->setScalar(1.0);
    }
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(lp, hp);
    // Build optimization problem
    auto optProb = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    optProb->addBoundConstraint(bnd);
    optProb->finalize(false,true,*outStream);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      auto rup = up->clone(); rup->randomize(-1.0,1.0);
      auto rzp = zp->clone(); rzp->randomize( 0.0,1.0);
      auto rpp = pp->clone(); rpp->randomize(-1.0,1.0);
      auto dup = up->clone(); dup->randomize(-1.0,1.0);
      auto dzp = zp->clone(); dzp->randomize( 0.0,1.0);
      con->checkApplyJacobian_1(*rup,*rzp,*dup,*rup,true,*outStream);
      con->checkApplyJacobian_2(*rup,*rzp,*dzp,*rup,true,*outStream);
      con->checkInverseJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      con->checkApplyAdjointHessian_11(*rup,*rzp,*rpp,*dup,*rup,true,*outStream);
      con->checkApplyAdjointHessian_21(*rup,*rzp,*rpp,*dzp,*rup,true,*outStream);
      con->checkApplyAdjointHessian_12(*rup,*rzp,*rpp,*dup,*rzp,true,*outStream);
      con->checkApplyAdjointHessian_22(*rup,*rzp,*rpp,*dzp,*rzp,true,*outStream);
      obj->checkGradient_1(*rup,*rzp,*dup,true,*outStream);
      obj->checkGradient_2(*rup,*rzp,*dzp,true,*outStream);
      obj->checkHessVec_11(*rup,*rzp,*dup,true,*outStream);
      obj->checkHessVec_12(*rup,*rzp,*dzp,true,*outStream);
      obj->checkHessVec_21(*rup,*rzp,*dup,true,*outStream);
      obj->checkHessVec_22(*rup,*rzp,*dzp,true,*outStream);
      robj->checkGradient(*rzp,*dzp,true,*outStream);
      robj->checkHessVec(*rzp,*dzp,true,*outStream);
      //optProb->check(*outStream);
    }

    // Solve optimization problem
    zp->setScalar(0.5);
    up->zero(); pp->zero();
    bool opt = parlist->sublist("Problem").get("Solve Optimization Problem",true);
    if (opt) {
      std::ifstream infile("control.txt");
      if (infile.good())
        assembler->inputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ifstream infile0; infile0.open("target.txt");
        if (infile0.good()) {
          for (int i = 0; i < dim; ++i) infile0 >> (*z0_ptr)[i];
          infile0.close();
        }
      }
      ROL::Solver<RealT> optSolver(optProb, *parlist);
      optSolver.solve(*outStream);
      con->outputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ofstream outfile; outfile.open("target.txt");
        for (const auto x : *z0_ptr) {
          outfile << std::scientific << std::setprecision(16);
          outfile << x << std::endl;
        }
        outfile.close();
      }
    }

    // Output
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    assembler->printCellAveragesPDE(pde,u_ptr,z_ptr);
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
