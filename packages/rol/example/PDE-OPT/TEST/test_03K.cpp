// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
    \brief Unit test for the mesh manager and the degree-of-freedom manager.
           Mesh type: BRICK with HEX CELLS and HCURL SPACE.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include <iostream>
#include <algorithm>

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/dofmanagerK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;
using basis_ptr = Intrepid2::BasisPtr<DeviceT, RealT, RealT>;
using scalar_view = typename MeshManager<RealT,DeviceT>::scalar_view;
using int_view = typename MeshManager<RealT,DeviceT>::int_view;

int main(int argc, char *argv[]) {

  ROL::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::ScopeGuard kokkosScope(argc, argv);
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_03.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshManager_Brick<RealT, DeviceT> meshmgr(*parlist);
    scalar_view nodes = meshmgr.getNodes();
    int_view cellToNodeMap = meshmgr.getCellToNodeMap();
    int_view cellToEdgeMap = meshmgr.getCellToEdgeMap();
    ROL::Ptr<std::vector<std::vector<std::vector<int> > > > sideSetsPtr = meshmgr.getSideSets();

    std::vector<std::vector<std::vector<int> > >  &sideSets = *sideSetsPtr;
    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl;  // << nodes;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl;  // << cellToNodeMap;
    *outStream << "Number of edges = " << meshmgr.getNumEdges() << std::endl;  // << cellToEdgeMap;
    // Print mesh info to file.
    std::ofstream meshfile;
    meshfile.open("cells.txt");
    for (int i=0; i<cellToNodeMap.extent_int(0); ++i) {
      meshfile << cellToNodeMap(i,0) << "  "
               << cellToNodeMap(i,1) << "  "
               << cellToNodeMap(i,2) << "  "
               << cellToNodeMap(i,3) << "  "
               << cellToNodeMap(i,4) << "  "
               << cellToNodeMap(i,5) << "  "
               << cellToNodeMap(i,6) << "  "
               << cellToNodeMap(i,7) << std::endl;
    }
    meshfile.close();
    meshfile.open("edges.txt");
    for (int i=0; i<cellToEdgeMap.extent_int(0); ++i) {
      meshfile << cellToEdgeMap(i,0)  << "  "
               << cellToEdgeMap(i,1)  << "  "
               << cellToEdgeMap(i,2)  << "  "
               << cellToEdgeMap(i,3)  << "  "
               << cellToEdgeMap(i,4)  << "  "
               << cellToEdgeMap(i,5)  << "  "
               << cellToEdgeMap(i,6)  << "  "
               << cellToEdgeMap(i,7)  << "  "
               << cellToEdgeMap(i,8)  << "  "
               << cellToEdgeMap(i,9)  << "  "
               << cellToEdgeMap(i,10) << "  "
               << cellToEdgeMap(i,11) << std::endl;
    }
    meshfile.close();
    meshfile.open("nodes.txt");
    for (int i=0; i<nodes.extent_int(0); ++i) {
      meshfile << nodes(i,0) << "  "
               << nodes(i,1) << "  "
               << nodes(i,2) << std::endl;
    }
    meshfile.close();
    meshfile.open("sideset.txt");
    for (int i=0; i<static_cast<int>(sideSets.size()); ++i) {
      for (int j=0; j<static_cast<int>(sideSets[i].size()); ++j) {
        if (sideSets[i][j].size() > 0) {
          for (int k=0; k<static_cast<int>(sideSets[i][j].size()); ++k) {
            meshfile << sideSets[i][j][k] << std::endl;
          }
        }
        meshfile << std::endl << std::endl;
      }
    }
    meshfile.close();

    basis_ptr basisPtrNedelec1 = ROL::makePtr<Intrepid2::Basis_HCURL_HEX_I1_FEM<DeviceT, RealT, RealT>>();

    std::vector<basis_ptr> basisPtrs(1, ROL::nullPtr);
    basisPtrs[0] = basisPtrNedelec1;

    ROL::Ptr<MeshManager<RealT, DeviceT>> meshmgrPtr = ROL::makePtrFromRef(meshmgr);

    DofManager<RealT, DeviceT> dofmgr(meshmgrPtr, basisPtrs);

    *outStream << "Number of node dofs = " << dofmgr.getNumNodeDofs() << std::endl;  // << *(dofmgr.getNodeDofs());
    *outStream << "Number of edge dofs = " << dofmgr.getNumEdgeDofs() << std::endl;  // << *(dofmgr.getEdgeDofs());
    *outStream << "Number of face dofs = " << dofmgr.getNumFaceDofs() << std::endl;  // << *(dofmgr.getFaceDofs());
    *outStream << "Number of void dofs = " << dofmgr.getNumVoidDofs() << std::endl;  // << *(dofmgr.getVoidDofs());
    *outStream << "Total number of dofs = " << dofmgr.getNumDofs() << std::endl;     // << *(dofmgr.getCellDofs());

    std::vector<std::vector<int> > fieldPattern = dofmgr.getFieldPattern();
    for (int i=0; i<dofmgr.getNumFields(); ++i) {
      *outStream << "\nField " << i << " pattern:   ";
      for (int j=0; j<dofmgr.getLocalFieldSize(i); ++j) {
        *outStream << fieldPattern[i][j] << " ";
      }
      *outStream << std::endl;
    }

    for (int i=0; i<dofmgr.getNumFields(); ++i) {
      *outStream << "\nField  " << i << std::endl;
      // *outStream << dofmgr.getFieldDofs(i);
    }

    bool correct = true;
    static const int checkDofs[] = {38, 42, 45, 41, 67, 71, 74, 70, 52, 53, 57, 56};
    for (int i=0; i<dofmgr.getLocalFieldSize(0); ++i) {
      correct = correct && ( (dofmgr.getFieldDofs(0))(11,i) == checkDofs[i] );
    }
    if (!correct) {
      errorFlag = -1;
    }

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
