// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Poisson topology optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"

#include "../../TOOLS/pdevector.hpp"

#include "pde_poisson_topOpt.hpp"
#include "mesh_poisson_topOpt.hpp"
#include "../src/filtered_compliance_robj.hpp"
#include "../src/volume_con.hpp"
#include "../src/obj_volume.hpp"
#include "../src/pde_filter.hpp"

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
    RealT tol(1.e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT initDens     = parlist->sublist("Problem").get("Initial Density",volFraction);
    const bool useFilter     = parlist->sublist("Problem").get("Use Filter", true);
    const bool normalizeObj  = parlist->sublist("Problem").get("Normalize Compliance", true);
    const bool volEq         = parlist->sublist("Problem").get("Use Volume Equality Constraint",true);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshManager_Poisson_TopOpt<RealT>>(*parlist);

    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Poisson_TopOpt<RealT>>
      pde = ROL::makePtr<PDE_Poisson_TopOpt<RealT>>(*parlist);

    // Initialize the filter PDE.
    ROL::Ptr<PDE<RealT>> pdeFilter = pde;

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Objective<RealT>> robj;
    ROL::Ptr<Assembler<RealT>> assembler, assemblerFilter;
    if (useFilter) {
      pdeFilter = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
      pde->setDensityFields(pdeFilter->getFields());
      robj = ROL::makePtr<TopOptFilteredComplianceObjective<RealT>>(
             pde,pdeFilter,meshMgr,comm,*parlist,*outStream);
      assembler = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj)->getAssembler();
      assemblerFilter = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj)->getFilterAssembler();
    }
    else {
      robj = ROL::makePtr<TopOptComplianceObjective<RealT>>(
             pde,meshMgr,comm,*parlist,*outStream);
      assembler = ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj)->getAssembler();
      assemblerFilter = assembler;
    }

    // Create vectors
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr, z_ptr;
    u_ptr = assembler->createStateVector();         u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();         p_ptr->putScalar(0.0);
    r_ptr = assembler->createResidualVector();      r_ptr->putScalar(0.0);
    z_ptr = assemblerFilter->createControlVector(); z_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, rp, zp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pdeFilter,assemblerFilter,*parlist);

    // Build volume objective function.
    ROL::Ptr<QoI<RealT>> qoi_vol;
    if (useFilter)
      qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(ROL::dynamicPtrCast<PDE_Filter<RealT>>(pdeFilter)->getDensityFE(),
                                                        volFraction);
    else
      qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getDensityFE(), volFraction);
    ROL::Ptr<TopOptVolumeConstraint<RealT>>
      con_vol = ROL::makePtr<TopOptVolumeConstraint<RealT>>(qoi_vol,assemblerFilter,zp);
    ROL::Ptr<ROL::Vector<RealT>> imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    ROL::Ptr<ROL::Vector<RealT>> ilp  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    ROL::Ptr<ROL::Vector<RealT>> iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    zp->zero(); con_vol->value(*ilp,*zp,tol);
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilp,iup);

    // Define bound constraint
    ROL::Ptr<ROL::Vector<RealT>> zlo = zp->clone(); zlo->setScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT>> zhi = zp->clone(); zhi->setScalar(1.0);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd = ROL::makePtr<ROL::Bounds<RealT>>(zlo,zhi);

    // Normalize compliance objective function
    zp->setScalar(initDens);
    RealT cs(1);
    if (normalizeObj) {
      if (useFilter)
        cs = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj)->normalize(*zp,tol);
      else
        cs = ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj)->normalize(*zp,tol);
    }
    *outStream << std::endl;
    *outStream << "Problem Data"          << std::endl;
    *outStream << "  SIMP Power:        "
               << parlist->sublist("Problem").get("SIMP Power",3.0) << std::endl;
    *outStream << "  Use Filter:        " << useFilter << std::endl;
    *outStream << "  Volume Fraction:   " << volFraction << std::endl;
    *outStream << "  Initial Density:   " << initDens    << std::endl;
    *outStream << "  Use Equality:      " << volEq       << std::endl;
    *outStream << "  Compliance Scale:  " << cs          << std::endl;
    *outStream << std::endl;

    // Initialize optimization problem.
    ROL::Ptr<ROL::Problem<RealT>>
      opt = ROL::makePtr<ROL::Problem<RealT>>(robj,zp);
    opt->addBoundConstraint(bnd);
    if (volEq) opt->addLinearConstraint("Volume",con_vol,imul);
    else       opt->addLinearConstraint("Volume",con_vol,imul,ibnd);
    opt->setProjectionAlgorithm(*parlist);
    opt->finalize(false,true,*outStream);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) opt->check(true,*outStream);

    // Solve optimization problem
    ROL::Solver<RealT> solver(opt,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj)->summarize(*outStream);
    con_vol->summarize(*outStream);
    ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj)->printToFile(*zp,*outStream);

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
