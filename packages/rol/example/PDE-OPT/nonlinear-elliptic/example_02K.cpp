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
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"

#include "pde_nonlinear_ellipticK.hpp"
#include "obj_nonlinear_ellipticK.hpp"

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
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int> > serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();
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

    RealT controlPenalty = parlist->sublist("Problem").get("Control penalty parameter",static_cast<RealT>(1.e-4));

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    auto pde = ROL::makePtr<PDE_Nonlinear_Elliptic<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    auto u_ptr = assembler->createStateVector(); u_ptr->randomize();
    auto up    = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    // Create state vector and set to zeroes
    auto p_ptr = assembler->createStateVector(); p_ptr->randomize();
    auto pp    = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler);
    // Create control vector and set to ones
    auto z_ptr = assembler->createControlVector(); z_ptr->randomize();
    auto zp    = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler);
    // Create residual vector and set to zeros
    auto r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto rp    = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    // Create state direction vector and set to random
    auto du_ptr = assembler->createStateVector(); du_ptr->randomize();
    auto dup    = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(du_ptr,pde,assembler);
    // Create control direction vector and set to random
    auto dz_ptr = assembler->createControlVector(); dz_ptr->randomize();
    auto dzp    = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(dz_ptr,pde,assembler);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize bound constraints.
    auto lo_ptr = assembler->createControlVector(); lo_ptr->putScalar(0.0);
    auto hi_ptr = assembler->createControlVector(); hi_ptr->putScalar(1.0);
    auto lop = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(lo_ptr,pde,assembler);
    auto hip = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(hi_ptr,pde,assembler);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateTracking_Nonlinear_Elliptic<RealT,DeviceT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_ControlPenalty_Nonlinear_Elliptic<RealT,DeviceT>>(pde->getFE());
    std::vector<RealT> weights = {static_cast<RealT>(1), controlPenalty};
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,weights,assembler);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int stochDim = 20;
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    std::vector<RealT> tmp = {static_cast<RealT>(-1),static_cast<RealT>(1)};
    std::vector<std::vector<RealT>> bounds(stochDim,tmp);
    auto bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    zp->zero();
    ROL::OptimizationProblem<RealT> opt(robj,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic", static_cast<RealT>(1));
    opt.setStochasticObjective(*parlist,sampler);

    parlist->sublist("Step").set("Type","Trust Region");
    ROL::OptimizationSolver<RealT> solver(opt,*parlist);
    solver.solve(*outStream);

    // Output.
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
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
