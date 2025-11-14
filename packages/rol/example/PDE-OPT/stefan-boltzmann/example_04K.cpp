// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_04.cpp
    \brief Shows how to solve the Stefan-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>
//
#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_stoch_stefan_boltzmannK.hpp"
#include "obj_stoch_stefan_boltzmannK.hpp"
#include "mesh_stoch_stefan_boltzmannK.hpp"

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
    std::string filename = "input_ex04.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    // Problem dimensions
    const int controlDim = 1;
    RealT tol(1e-8);

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    auto pde = ROL::makePtr<StochasticStefanBoltzmannPDE<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    auto u_ptr   = assembler->createStateVector();    u_ptr->putScalar(0.0);
    auto p_ptr   = assembler->createStateVector();    p_ptr->putScalar(0.0);
    auto r_ptr   = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    auto zbc_ptr = assembler->createControlVector();  zbc_ptr->putScalar(280.0);
    auto zp_ptr  = ROL::makePtr<std::vector<RealT>>(controlDim,0.0);
    auto up      = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto pp      = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler);
    auto rp      = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    auto zbc     = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zbc_ptr,pde,assembler);
    auto zparam  = ROL::makePtr<ROL::StdVector<RealT>>(zp_ptr);
    auto zp      = ROL::makePtr<PDE_OptVector<RealT>>(zbc, zparam);
    auto xp      = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(3,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateCost<RealT,DeviceT>>(pde->getVolFE(),*parlist);
    qoi_vec[1] = ROL::makePtr<QoI_ControlCost<RealT,DeviceT>>(
      pde->getVolFE(),pde->getBdryFE(0),pde->getBdryCellLocIds(0),*parlist);
    qoi_vec[2] = ROL::makePtr<QoI_AdvectionCost<RealT,DeviceT>>();
    auto std_obj = ROL::makePtr<StochasticStefanBoltzmannStdObjective3<RealT>>(*parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,std_obj,assembler);
    auto objReduced = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Bounds for boundary control
    RealT lower_bc = parlist->sublist("Problem").get("Lower Control Bound", 280.0);
    RealT upper_bc = parlist->sublist("Problem").get("Upper Control Bound", 370.0);
    auto zlo_bc_ptr = assembler->createControlVector();
    auto zhi_bc_ptr = assembler->createControlVector();
    zlo_bc_ptr->putScalar(static_cast<RealT>(lower_bc));
    zhi_bc_ptr->putScalar(static_cast<RealT>(upper_bc));
    auto zlo_bc = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zlo_bc_ptr,pde,assembler);
    auto zhi_bc = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zhi_bc_ptr,pde,assembler);
    // Bounds for advection control
    RealT lower = parlist->sublist("Problem").get("Lower Advection Bound",-100.0);
    RealT upper = parlist->sublist("Problem").get("Upper Advection Bound", 100.0);
    auto zlo_param_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,lower);
    auto zhi_param_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,upper);
    auto zlo_adv = ROL::makePtr<ROL::StdVector<RealT>>(zlo_param_ptr);
    auto zhi_adv = ROL::makePtr<ROL::StdVector<RealT>>(zhi_param_ptr);
    // Combined bounds
    auto zlop = ROL::makePtr<PDE_OptVector<RealT>>(zlo_bc, zlo_adv);
    auto zhip = ROL::makePtr<PDE_OptVector<RealT>>(zhi_bc, zhi_adv);
    auto zbnd = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    auto ulo_ptr = assembler->createStateVector();
    auto uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    auto ulop = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(ulo_ptr,pde,assembler,*parlist);
    auto uhip = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(uhi_ptr,pde,assembler,*parlist);
    auto ubnd = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    auto bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT> >(ubnd,zbnd);

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Ptr<ROL::Problem<RealT>> problem;
    bool useReducedSpace = parlist->sublist("Problem").get("Use Reduced Space",true);
    if (useReducedSpace) {
      problem = ROL::makePtr<ROL::Problem<RealT>>(objReduced,zp);
      problem->addBoundConstraint(zbnd);
    }
    else {
      problem = ROL::makePtr<ROL::Problem<RealT>>(obj,xp);
      problem->addBoundConstraint(bnd);
      problem->addConstraint("PDE",con,pp);
    }
    problem->finalize(false,true,*outStream);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) problem->check(true,*outStream);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Solver<RealT> solver(problem,*parlist);
    (*zp_ptr)[0] = parlist->sublist("Problem").get("Advection Magnitude",0.0);
    u_ptr->putScalar(450.0);

    bool solveFS = parlist->sublist("Problem").get("Initial Solve for Full Space",true); 
    if (solveFS) con->solve(*rp,*up,*zp,tol);

    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    con->outputTpetraVector(zbc_ptr,"control.txt");
    con->outputTpetraVector(u_ptr,"state.txt");
    *outStream << std::endl << "Advection value: " << (*zp_ptr)[0] << std::endl;
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));

    /*************************************************************************/
    /***************** CHECK RESIDUAL NORM ***********************************/
    /*************************************************************************/
    *outStream << "Residual Norm: " << res[0] << std::endl << std::endl;
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
