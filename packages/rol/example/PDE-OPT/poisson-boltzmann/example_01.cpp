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

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_poisson_boltzmann.hpp"
#include "obj_poisson_boltzmann.hpp"

#include "ROL_Solver.hpp"

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
      meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Poisson_Boltzmann<RealT>>
      pde = ROL::makePtr<PDE_Poisson_Boltzmann<RealT>>(*parlist);
    ROL::Ptr<PDE_Constraint<RealT>>
      con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_L2Tracking_Poisson_Boltzmann<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_Poisson_Boltzmann<RealT>>(pde->getFE());
    ROL::Ptr<StdObjective_Poisson_Boltzmann<RealT>>
      std_obj = ROL::makePtr<StdObjective_Poisson_Boltzmann<RealT>>(*parlist);
    ROL::Ptr<PDE_Objective<RealT>>
      obj = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,con->getAssembler());

    // Create vectors
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr = con->getAssembler()->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> z_ptr = con->getAssembler()->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<>> p_ptr = con->getAssembler()->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> r_ptr = con->getAssembler()->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<>> du_ptr = con->getAssembler()->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> dz_ptr = con->getAssembler()->createControlVector();
    u_ptr->randomize();
    z_ptr->putScalar(1.0);
    p_ptr->putScalar(0.0);
    r_ptr->putScalar(0.0);
    du_ptr->randomize();
    dz_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> up, zp, pp, rp, dup, dzp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,con->getAssembler());
    zp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,con->getAssembler());
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,con->getAssembler());
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,con->getAssembler());
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,con->getAssembler());
    dzp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,con->getAssembler());
    // Create ROL SimOpt vectors
    ROL::Ptr<ROL::Vector_SimOpt<RealT>> x, d;
    x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);
    d = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(dup,dzp);

    // Run derivative checks
    obj->checkGradient(*x,*d,true,*outStream);
    obj->checkHessVec(*x,*d,true,*outStream);
    con->checkApplyJacobian(*x,*d,*up,true,*outStream);
    con->checkApplyAdjointHessian(*x,*dup,*d,*x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,*d,*x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    ROL::Ptr<ROL::Problem<RealT>>
      problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addConstraint("PDE",con,pp);
    problem->finalize(false,true,*outStream);
    ROL::Solver<RealT> solver(problem,*parlist);
    solver.solve(*outStream);

    // Output.
    con->getAssembler()->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state.txt");
    con->outputTpetraVector(z_ptr,"control.txt");

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
