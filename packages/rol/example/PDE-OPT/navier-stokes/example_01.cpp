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
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Solver.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
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

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<PDE_NavierStokes<RealT>>
      pde = ROL::makePtr<PDE_NavierStokes<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>>
      con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT>>
      pdecon = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    ROL::Ptr<Assembler<RealT>> assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr = assembler->createStateVector();    u_ptr->randomize();
    ROL::Ptr<Tpetra::MultiVector<>> p_ptr = assembler->createStateVector();    p_ptr->randomize();
    ROL::Ptr<Tpetra::MultiVector<>> z_ptr = assembler->createControlVector();  z_ptr->randomize();
    ROL::Ptr<Tpetra::MultiVector<>> r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, zp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);

    // Initialize quadratic objective function.
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_NavierStokes<RealT>>(*parlist,
                                                             pde->getVelocityFE(),
                                                             pde->getPressureFE(),
                                                             pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_NavierStokes<RealT>>(pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getVelocityBdryFE(),
                                                                 pde->getBdryCellLocIds(),
                                                                 pde->getFieldHelper());
    ROL::Ptr<StdObjective_NavierStokes<RealT>>
      std_obj = ROL::makePtr<StdObjective_NavierStokes<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>
      obj = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>>
      robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<ROL::Vector<RealT>> dup = up->clone(); dup->randomize(-1.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> dzp = zp->clone(); dzp->randomize(-1.0,1.0);
      ROL::Vector_SimOpt<RealT> d(dup,dzp);
      obj->checkGradient(x,d,true,*outStream);
      obj->checkHessVec(x,d,true,*outStream);
      con->checkApplyJacobian_1(*up,*zp,*dup,*up,true,*outStream);
      con->checkApplyJacobian_2(*up,*zp,*dzp,*up,true,*outStream);
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      robj->checkGradient(*zp,*dzp,true,*outStream);
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    up->zero(); zp->zero();
    ROL::Ptr<ROL::Problem<RealT>> optProb;
    bool useFullSpace = parlist->sublist("Problem").get("Full space",false);
    if ( useFullSpace ) {
      std::string step = parlist->sublist("Step").get("Type", "Composite Step");
      ROL::EStep els = ROL::StringToEStep(step);
      switch( els ) {
        case ROL::STEP_COMPOSITESTEP:
        case ROL::STEP_FLETCHER: {
          RealT tol(1.e-8);
          bool initSolve = parlist->sublist("Problem").get("Solve state for full space",true);
          if (initSolve) {
            con->solve(*rp,*up,*zp,tol);
            pdecon->outputTpetraVector(u_ptr,"state_uncontrolled.txt");
          }
          optProb = ROL::makePtr<ROL::Problem<RealT>>(obj, makePtrFromRef(x));
          optProb->addConstraint("PDE", con, rp);
          break;
        }
        default: {
          *outStream << "ERROR: Unsupported step." << std::endl;      
          errorFlag = 1; 
        }
      }
    }
    else {
      parlist->sublist("Step").set("Type","Trust Region");
      optProb = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    }
    optProb->finalize(false,true,*outStream);
    ROL::Solver<RealT> optSolver(optProb, *parlist);
    optSolver.solve(*outStream);

    // Output.
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //pdecon->outputTpetraData();
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
