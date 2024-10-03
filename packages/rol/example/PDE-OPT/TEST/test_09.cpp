// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_09.cpp
    \brief Unit test for the mesh manager and the degree-of-freedom manager.
           Mesh type: INTERVAL with LINE CELLS and HGRAD SPACE.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"

#include <iostream>
#include <algorithm>

#include "../TOOLS/meshmanager.hpp"
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
    std::string filename = "input_09.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshManager_Interval<RealT> meshmgr(*parlist);
    ROL::Ptr<Intrepid::FieldContainer<RealT> > nodesPtr = meshmgr.getNodes();
    ROL::Ptr<Intrepid::FieldContainer<int> >   cellToNodeMapPtr = meshmgr.getCellToNodeMap();
    ROL::Ptr<std::vector<std::vector<std::vector<int> > > > sideSetsPtr = meshmgr.getSideSets(); 

    Intrepid::FieldContainer<RealT> &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    std::vector<std::vector<std::vector<int> > >  &sideSets = *sideSetsPtr;
    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl << nodes;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl << cellToNodeMap;
    // Print mesh info to file.
    std::ofstream meshfile;
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

    ROL::Ptr<Intrepid::Basis_HGRAD_LINE_Cn_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrL1 =
      ROL::makePtr<Intrepid::Basis_HGRAD_LINE_Cn_FEM<RealT, Intrepid::FieldContainer<RealT> >>(1, Intrepid::POINTTYPE_EQUISPACED);

    ROL::Ptr<Intrepid::Basis_HGRAD_LINE_Cn_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrL2 =
      ROL::makePtr<Intrepid::Basis_HGRAD_LINE_Cn_FEM<RealT, Intrepid::FieldContainer<RealT> >>(2, Intrepid::POINTTYPE_EQUISPACED);

    std::vector<ROL::Ptr<Intrepid::Basis<RealT, Intrepid::FieldContainer<RealT> > > > basisPtrs(3, ROL::nullPtr);
    basisPtrs[0] = basisPtrL2;
    basisPtrs[1] = basisPtrL1;
    basisPtrs[2] = basisPtrL2;

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

    /*bool correct = true;
    static const int checkDofs[] = {20, 23, 35, 32, 55, 63, 69, 61, 81};
    for (int i=0; i<dofmgr.getLocalFieldSize(2); ++i) {
      correct = correct && ( (*(dofmgr.getFieldDofs(2)))(5,i) == checkDofs[i] );
    }
    if (!correct) {
      errorFlag = -1;
    }*/

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
