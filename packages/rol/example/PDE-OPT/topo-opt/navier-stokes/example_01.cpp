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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

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

#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"

#include "pde_navier-stokes.hpp"
#include "obj_navier-stokes.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
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

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    bool output = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    output = (iprint > 0) && (myRank==0) && output;
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",output);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>>
      meshMgr  = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<PDE_NavierStokes<RealT>>
      pde      = ROL::makePtr<PDE_NavierStokes<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>>
      con      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT>>
      pdecon   = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    ROL::Ptr<Assembler<RealT>>
      assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, z_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->randomize();
    p_ptr = assembler->createStateVector();    p_ptr->randomize();
    z_ptr = assembler->createControlVector();  z_ptr->randomize();
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, zp, rp, xp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize quadratic objective function.
    ROL::Ptr<QoI<RealT>> qoi
      = ROL::makePtr<QoI_Power_NavierStokes<RealT>>(*parlist,
                                                    pde->getVelocityFE(),
                                                    pde->getPressureFE(),
                                                    pde->getControlFE(),
                                                    pde->getStateFieldInfo(),
                                                    pde->getControlFieldInfo());
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>
      obj  = ROL::makePtr<PDE_Objective<RealT>>(qoi,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>>
      robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    // Initialize volume constraint.
    RealT volFraction = parlist->sublist("Problem").get("Volume Fraction",0.5);
    ROL::Ptr<QoI<RealT>>
      qoi_vol = ROL::makePtr<QoI_Volume_NavierStokes<RealT>>(*parlist,
                                                             pde->getControlFE(),
                                                             pde->getControlFieldInfo());
    ROL::Ptr<ROL::Objective<RealT>>
      obj_vol = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_vol,assembler);
    ROL::Ptr<ROL::Constraint<RealT>>
      icon    = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(obj_vol);
    ROL::Ptr<ROL::Vector<RealT>>
      iup     = ROL::makePtr<ROL::SingletonVector<RealT>>(0.0);
    ROL::Ptr<ROL::Vector<RealT>>
      imul    = ROL::makePtr<ROL::SingletonVector<RealT>>(0.0);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      ibnd    = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);

    // Build bound constraint
    ROL::Ptr<ROL::Vector<RealT>> lp = zp->clone(); lp->setScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> hp = zp->clone(); hp->setScalar(1.0);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
    bnd = ROL::makePtr<ROL::Bounds<RealT>>(lp, hp);
    // Build optimization problem
    ROL::Ptr<ROL::Problem<RealT>> optProb
      = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    optProb->addBoundConstraint(bnd);
    optProb->addLinearConstraint("Volume",icon,imul); //, ibnd);
    optProb->setProjectionAlgorithm(*parlist);
    bool useProj = parlist->sublist("Problem").get("Project Linear Constraints",true);
    optProb->finalize(!useProj,true,*outStream);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<ROL::Vector<RealT>> dup = up->clone(); dup->randomize();
      ROL::Ptr<ROL::Vector<RealT>> dzp = zp->clone(); dzp->randomize();
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
    up->zero(); pp->zero();
    bool opt = parlist->sublist("Problem").get("Solve Optimization Problem",true);
    if (opt) {
      std::ifstream infile("control.txt");
      if (infile.good()) assembler->inputTpetraVector(z_ptr,"control.txt");
      else               zp->setScalar(volFraction);
      Teuchos::Time algoTimer("Algorithm Time", true);
      ROL::Solver<RealT> optSolver(optProb,*parlist);
      optSolver.solve(*outStream);
      algoTimer.stop();
      *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";
      pdecon->outputTpetraVector(z_ptr,"control.txt");
    }
    else {
      zp->setScalar(1.0);
    }

    // Output
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //pdecon->outputTpetraData();

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
