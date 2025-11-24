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
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "../TOOLS/linearpdeconstraintK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_stokesK.hpp"
#include "mesh_stokesK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Stokes<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    auto pde = ROL::makePtr<PDE_Stokes<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<Linear_PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);
    con->outputTpetraData();

    // Create state vector and set to zeroes
    auto u_ptr = assembler->createStateVector();    u_ptr->randomize();
    auto r_ptr = assembler->createResidualVector(); r_ptr->randomize();
    auto z_ptr = assembler->createControlVector();  z_ptr->putScalar(0);
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      auto pp  = up->clone(); pp->randomize(-1.0,1.0);
      auto dup = up->clone(); dup->randomize(-1.0,1.0);

      *outStream << std::endl << "Check Jacobian_1 of Constraint" << std::endl;
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << std::endl << "Check Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian_1(*pp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Constraint Solve" << std::endl;
      con->checkSolve(*up,*zp,*rp,true,*outStream);
      *outStream << std::endl << "Check Inverse Jacobian_1 of Constraint" << std::endl;
      con->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Inverse Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkInverseAdjointJacobian_1(*rp,*pp,*up,*zp,true,*outStream);
    }

    RealT tol(1.e-8);
    up->zero();
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraData();
    con->outputTpetraVector(u_ptr,"state.txt");
    assembler->printMeshData(*outStream);
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
