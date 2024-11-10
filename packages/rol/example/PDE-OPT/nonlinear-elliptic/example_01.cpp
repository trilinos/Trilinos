// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Poisson-Boltzmann problem.
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
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"

#include "pde_nonlinear_elliptic.hpp"
#include "obj_nonlinear_elliptic.hpp"

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

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    RealT controlPenalty = parlist->sublist("Problem").get("Control penalty parameter",static_cast<RealT>(1.e-4));

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Nonlinear_Elliptic<RealT> > pde
      = ROL::makePtr<PDE_Nonlinear_Elliptic<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
    u_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > p_ptr = assembler->createStateVector();
    p_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > pp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    // Create control vector and set to ones
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = assembler->createControlVector();
    z_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);
    // Create residual vector and set to zeros
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    // Create state direction vector and set to random
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    du_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler);
    // Create control direction vector and set to random
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    dz_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > dzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize bound constraints.
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = assembler->createControlVector();
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > hip
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateTracking_Nonlinear_Elliptic<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_ControlPenalty_Nonlinear_Elliptic<RealT>>(pde->getFE());
    std::vector<RealT> weights = {static_cast<RealT>(1), controlPenalty};
    ROL::Ptr<PDE_Objective<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,weights,assembler);
    ROL::Ptr<ROL::Objective<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    // Run derivative checks
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);

    con->checkApplyJacobian(x,d,*rp,true,*outStream);
    con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
    con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);
    algo.run(*zp,*robj,*bnd,true,*outStream);

    // Output.
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");

    Teuchos::Array<RealT> res(1,0);
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
