// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the obstacle problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_obstacleK.hpp"
#include "energy_objectiveK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

template <class Real>
Real evaluateUpperBound(const std::vector<Real> & coord) {
  Real distance(0), one(1);
  if (coord[1] < coord[0]) {
    if (coord[1] < one - coord[0]) distance = coord[1];
    else                           distance = one - coord[0];
  }
  else {
    if (coord[1] < one - coord[0]) distance = coord[0];
    else                           distance = one - coord[1];
  }
  return distance;
}


template <class Real,class DeviceType>
void computeUpperBound(const ROL::Ptr<Tpetra::MultiVector<>> & ubVec,
                       const ROL::Ptr<const FE<Real,DeviceType>> & fe,
                       const Kokkos::DynRankView<Real,DeviceType> cellNodes,
                       const Kokkos::DynRankView<int,DeviceType> cellDofs,
                       const Teuchos::Array<typename Tpetra::Map<>::global_ordinal_type> & cellIds) {
  int c = fe->gradN().extent_int(0);
  int f = fe->gradN().extent_int(1);
  int d = fe->gradN().extent_int(3);
  Kokkos::DynRankView<Real,DeviceType> dofPoints("dofPoints",c,f,d);
  fe->computeDofCoords(dofPoints, cellNodes);
  
  std::vector<Real> coord(d);
  for (int i = 0; i < c; ++i) {
    for (int j = 0; j < f; ++j) {
      int fidx = cellDofs(cellIds[i],j);
      if (ubVec->getMap()->isNodeGlobalElement(fidx)) {
        for (int k = 0; k < d; ++k)
          coord[k] = dofPoints(i,j,k);
        ubVec->replaceGlobalValue(fidx,
                                  0,
                                  evaluateUpperBound<Real>(coord));
      }
    }
  }
} 


int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
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

    /*** Initialize PDE describing the obstacle problem ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    auto pde = ROL::makePtr<PDE_Obstacle<RealT,DeviceT>>(*parlist);
    auto obj = ROL::makePtr<EnergyObjective<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    auto assembler = obj->getAssembler();

    // Create state vectors
    auto u_ptr = assembler->createStateVector();
    auto du_ptr = assembler->createStateVector();
    u_ptr->randomize(); du_ptr->randomize();
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(du_ptr,pde,assembler);

    // Build bound constraints
    auto lo_ptr = assembler->createStateVector();
    auto hi_ptr = assembler->createStateVector();
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    computeUpperBound<RealT,DeviceT>(hi_ptr,pde->getFE(),
                             pde->getCellNodes(),
                             assembler->getDofManager()->getCellDofs(),
                             assembler->getCellIds());
    auto lop = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(lo_ptr,pde,assembler);
    auto hip = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(hi_ptr,pde,assembler);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Run derivative checks
    obj->checkGradient(*up,*dup,true,*outStream);
    obj->checkHessVec(*up,*dup,true,*outStream);

    // Build optimization problem
    du_ptr->putScalar(0.4);
    auto opt = ROL::makePtr<ROL::Problem<RealT>>(obj,up);
    opt->addBoundConstraint(bnd);
    opt->finalize(false,true,*outStream);

    up->set(*dup);
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::Solver<RealT> solverTR(opt,*parlist);
    std::clock_t timerTR = std::clock();
    solverTR.solve(*outStream);
    *outStream << "Trust Region Time: "
               << static_cast<RealT>(std::clock()-timerTR)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    up->set(*dup);
    parlist->sublist("Step").set("Type","Primal-Dual Active Set");
    ROL::Solver<RealT> solverPDAS(opt,*parlist);
    std::clock_t timerPDAS = std::clock();
    solverPDAS.solve(*outStream);
    *outStream << "Primal Dual Active Set Time: "
               << static_cast<RealT>(std::clock()-timerPDAS)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    // Output.
    assembler->printMeshData(*outStream);
    assembler->outputTpetraVector(u_ptr,"state.txt");
    assembler->outputTpetraVector(hi_ptr,"upperBound.txt");

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
