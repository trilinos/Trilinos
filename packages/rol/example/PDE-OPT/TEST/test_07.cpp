// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_07.cpp
    \brief Unit test for the MeshReader mesh manager.
*/

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

#include "../TOOLS/meshreader.hpp"
#include "../TOOLS/dofmanager.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
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
    std::string filename = "input_07.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshReader<RealT> meshmgr(*parlist);

    ROL::Ptr<Intrepid::FieldContainer<RealT> > nodesPtr = meshmgr.getNodes();
    ROL::Ptr<Intrepid::FieldContainer<int> >   cellToNodeMapPtr = meshmgr.getCellToNodeMap();
    ROL::Ptr<Intrepid::FieldContainer<int> >   cellToEdgeMapPtr = meshmgr.getCellToEdgeMap();
    ROL::Ptr<Intrepid::FieldContainer<int> >   cellToFaceMapPtr = meshmgr.getCellToFaceMap();
    ROL::Ptr<std::vector<std::vector<std::vector<int> > > > sideSetsPtr = meshmgr.getSideSets(); 

    Intrepid::FieldContainer<RealT> &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    Intrepid::FieldContainer<int>   &cellToEdgeMap = *cellToEdgeMapPtr;
    Intrepid::FieldContainer<int>   &cellToFaceMap = *cellToFaceMapPtr;
    std::vector<std::vector<std::vector<int> > >  &sideSets = *sideSetsPtr;
    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl << nodes;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl << cellToNodeMap;
    *outStream << "Number of edges = " << meshmgr.getNumEdges() << std::endl << cellToEdgeMap;
    *outStream << "Number of faces = " << meshmgr.getNumFaces() << std::endl << cellToFaceMap;
    // Print sideset info to file.
    std::ofstream meshfile;
    meshfile.open("sideset.txt");
    for (int i=0; i<static_cast<int>(sideSets.size()); ++i) {
      for (int j=0; j<static_cast<int>(sideSets[i].size()); ++j) {
        if (sideSets[i][j].size() > 0) {
          for (int k=0; k<static_cast<int>(sideSets[i][j].size()); ++k) {
            meshfile << i << " " << j << " " << k << " " << sideSets[i][j][k] << std::endl;
          }
        }
        meshfile << std::endl;
      }
      meshfile << std::endl << std::endl;
    }
    meshfile.close();

    ROL::Ptr<Intrepid::Basis_HGRAD_HEX_C1_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrQ1 =
      ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<RealT, Intrepid::FieldContainer<RealT> >>();

    ROL::Ptr<Intrepid::Basis_HGRAD_HEX_C2_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrQ2 =
      ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<RealT, Intrepid::FieldContainer<RealT> >>();

    std::vector<ROL::Ptr<Intrepid::Basis<RealT, Intrepid::FieldContainer<RealT> > > > basisPtrs(3, ROL::nullPtr);
    basisPtrs[0] = basisPtrQ2;
    basisPtrs[1] = basisPtrQ1;
    basisPtrs[2] = basisPtrQ2;

    ROL::Ptr<MeshManager<RealT> > meshmgrPtr = ROL::makePtrFromRef(meshmgr);

    DofManager<RealT> dofmgr(meshmgrPtr, basisPtrs);

    *outStream << "Number of node dofs = " << dofmgr.getNumNodeDofs() << std::endl << *(dofmgr.getNodeDofs());
    *outStream << "Number of edge dofs = " << dofmgr.getNumEdgeDofs() << std::endl << *(dofmgr.getEdgeDofs());
    *outStream << "Number of face dofs = " << dofmgr.getNumFaceDofs() << std::endl << *(dofmgr.getFaceDofs());
    *outStream << "Number of void dofs = " << dofmgr.getNumVoidDofs() << std::endl << *(dofmgr.getVoidDofs());
    *outStream << "Total number of dofs = " << dofmgr.getNumDofs() << std::endl << *(dofmgr.getCellDofs());

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
      *outStream << *(dofmgr.getFieldDofs(i));
    }

    bool correct = true;
    static const int checkDofs[] = {11, 8, 38, 41, 23, 20, 44, 47, 53, 89, 91, 93, 61, 95, 97, 99, 71, 69, 101, 103, 109, 127, 129, 131, 133, 135, 141};
    for (int i=0; i<dofmgr.getLocalFieldSize(2); ++i) {
      correct = correct && ( (*(dofmgr.getFieldDofs(2)))(2,i) == checkDofs[i] );
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
