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
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_poisson_boltzmannK.hpp"
#include "obj_poisson_boltzmannK.hpp"

#include "ROL_Solver.hpp"

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
  ROL::Ptr<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
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
    auto parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT, DeviceT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    auto pde = ROL::makePtr<PDE_Poisson_Boltzmann<RealT, DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT, DeviceT>>(pde, meshMgr, comm, *parlist, *outStream);
    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT, DeviceT>>> qoi_vec(2, ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_L2Tracking_Poisson_Boltzmann<RealT, DeviceT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_Poisson_Boltzmann<RealT, DeviceT>>(pde->getFE());
    auto std_obj = ROL::makePtr<StdObjective_Poisson_Boltzmann<RealT>>(*parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT, DeviceT>>(qoi_vec, std_obj, con->getAssembler());

    // Create vectors
    auto u_ptr = con->getAssembler()->createStateVector();
    auto z_ptr = con->getAssembler()->createControlVector();
    auto p_ptr = con->getAssembler()->createStateVector();
    auto r_ptr = con->getAssembler()->createResidualVector();
    auto du_ptr = con->getAssembler()->createStateVector();
    auto dz_ptr = con->getAssembler()->createControlVector();
    u_ptr->randomize();
    z_ptr->putScalar(1.0);
    p_ptr->putScalar(0.0);
    r_ptr->putScalar(0.0);
    du_ptr->randomize();
    dz_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> up, zp, pp, rp, dup, dzp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(u_ptr,pde,con->getAssembler());
    zp  = ROL::makePtr<PDE_PrimalOptVector<RealT, DeviceT>>(z_ptr,pde,con->getAssembler());
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(p_ptr,pde,con->getAssembler());
    rp  = ROL::makePtr<PDE_DualSimVector<RealT, DeviceT>>(r_ptr,pde,con->getAssembler());
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(du_ptr,pde,con->getAssembler());
    dzp = ROL::makePtr<PDE_PrimalOptVector<RealT, DeviceT>>(dz_ptr,pde,con->getAssembler());
    // Create ROL SimOpt vectors
    ROL::Ptr<ROL::Vector_SimOpt<RealT>> x, d;
    x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up, zp);
    d = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(dup, dzp);

    // Run derivative checks
    obj->checkGradient(*x,*d,true,*outStream);
    obj->checkHessVec(*x,*d,true,*outStream);
    con->checkApplyJacobian(*x,*d,*up,true,*outStream);
    con->checkApplyAdjointHessian(*x,*dup,*d,*x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,*d,*x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

    RealT tol(1.e-8);
    con->solve(*rp, *up, *zp, tol);
    auto problem = ROL::makePtr<ROL::Problem<RealT>>(obj, x);
    problem->addConstraint("PDE", con, pp);
    problem->finalize(false, true, *outStream);
    Teuchos::Time algoTimer("Algorithm Time", true);
    ROL::Solver<RealT> solver(problem, *parlist);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    con->getAssembler()->printMeshData(*outStream);
    con->solve(*rp, *up, *zp, tol);
    con->outputTpetraVector(u_ptr, "state.txt");
    con->outputTpetraVector(z_ptr, "control.txt");

    Teuchos::Array<RealT> res(1, 0);
    con->value(*rp, *up, *zp, tol);
    r_ptr->norm2(res.view(0, 1));
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
