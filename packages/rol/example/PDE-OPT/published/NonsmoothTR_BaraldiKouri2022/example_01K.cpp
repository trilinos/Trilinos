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
#include "ROL_GlobalMPISession.hpp"
#include "ROL_ParameterList.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_LinearCombinationObjective.hpp"

#include "../../TOOLS/pdevectorK.hpp"
#include "../../TOOLS/pdeobjectiveK.hpp"
#include "../../TOOLS/meshreaderK.hpp"

#include "src/pde_elasticityK.hpp"
#include "src/pde_filterK.hpp"
#include "src/filtered_compliance_robjK.hpp"
#include "src/volume_conK.hpp"
#include "src/obj_volumeK.hpp"

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
  auto comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0))
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  // *** Example body.
  try {
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    // Retrieve parameters.
    int probDim                = 3;
    const RealT volFraction    = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT initDens       = volFraction;
    const std::string minType  = "Compliance";
    const std::string hessAppr = "None";
    parlist->sublist("Problem").set("Problem Dimension",probDim);
    parlist->sublist("Mesh").set("File Name","meshfiles/truss3d.txt");
    parlist->sublist("Problem").set("Use Inexact Linear Solves",true);
    parlist->sublist("General").set("Inexact Objective Function",true);
    parlist->sublist("General").set("Inexact Gradient",true);

    /*** Initialize main data structure. ***/
    int nProcs = comm->getSize();
    if (nProcs == 1) nProcs = 0;
    auto meshMgr = ROL::makePtr<MeshReader<RealT,DeviceT>>(*parlist,nProcs);
    // Initialize PDE describing elasticity equations.
    auto pde = ROL::makePtr<PDE_Elasticity<RealT,DeviceT>>(*parlist);

    // Initialize the filter PDE.
    auto pdeFilter = ROL::makePtr<PDE_Filter<RealT,DeviceT>>(*parlist);
    pde->setDensityFields(pdeFilter->getFields());

    // Initialize reduced compliance objective function.
    auto robj_com = ROL::makePtr<TopOptFilteredComplianceObjective<RealT,DeviceT>>(
                    pde,pdeFilter,meshMgr,comm,*parlist,*outStream);
    auto assembler = robj_com->getAssembler();
    auto assemblerFilter = robj_com->getFilterAssembler();

    // Create vectors.
    auto u_ptr = assembler->createStateVector();         u_ptr->putScalar(0.0);
    auto p_ptr = assembler->createStateVector();         p_ptr->putScalar(0.0);
    auto r_ptr = assembler->createResidualVector();      r_ptr->putScalar(0.0);
    auto z_ptr = assemblerFilter->createControlVector(); z_ptr->putScalar(1.0);
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pdeFilter,assemblerFilter,*parlist);
    auto imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);

    // Build volume objective function.
    auto qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT,DeviceT>>(pdeFilter->getDensityFE(),volFraction);
    auto con_vol = ROL::makePtr<TopOptVolumeConstraint<RealT,DeviceT>>(qoi_vol,assemblerFilter,zp);

    // Normalize compliance objective function
    zp->setScalar(initDens);
    RealT cs = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT,DeviceT>>(robj_com)->normalize(*zp,tol);

    // Output problem details
    *outStream << std::endl;
    *outStream << "Problem Data"          << std::endl;
    *outStream << "  Dimension:         " << probDim << std::endl;
    *outStream << "  SIMP Power:        "
               << parlist->sublist("Problem").get("SIMP Power",3.0) << std::endl;
    *outStream << "  Volume Fraction:   " << volFraction << std::endl;
    *outStream << "  Initial Density:   " << initDens    << std::endl;
    *outStream << "  Compliance Scale:  " << cs  << std::endl;
    *outStream << std::endl;

    // Create objective, constraint, multiplier and bounds

    // Initialize bound constraints.
    RealT lval = 0.0, uval = 1.0;
    auto lop = zp->clone(); lop->setScalar(lval);
    auto hip = zp->clone(); hip->setScalar(uval);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Set up optimization problem.
    auto prob = ROL::makePtr<ROL::Problem<RealT>>(robj_com,zp);
    prob->addBoundConstraint(bnd);
    prob->addLinearConstraint("Volume",con_vol,imul);
    prob->setProjectionAlgorithm(*parlist);
    prob->finalize(false,true,*outStream);

    // Check derivatives.
    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) prob->check(true,*outStream);

    // Solve optimization problem.
    ROL::Ptr<ROL::Solver<RealT>> solver;

    // TRSPG
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time trspgTimer("TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    trspgTimer.stop();
    *outStream << "Total optimization time = " << trspgTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    robj_com->printToFile(*zp,*outStream);

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
