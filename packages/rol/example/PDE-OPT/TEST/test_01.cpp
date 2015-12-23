// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Unit tests for the mesh manager and the degree-of-freedom manager.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"

#include <iostream>
#include <algorithm>

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/dofmanager.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshManager_Rectangle<RealT> meshmgr(*parlist);
    Intrepid::FieldContainer<RealT> nodes;
    Intrepid::FieldContainer<int>   cellToNodeMap;
    Intrepid::FieldContainer<int>   cellToEdgeMap;
    meshmgr.getNodes(nodes);
    meshmgr.getCellToNodeMap(cellToNodeMap);
    meshmgr.getCellToEdgeMap(cellToEdgeMap);
    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl << nodes;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl << cellToNodeMap;
    *outStream << "Number of edges = " << meshmgr.getNumEdges() << std::endl << cellToEdgeMap;
    // Print mesh to file.
    std::ofstream meshfile;
    meshfile.open("mesh.txt");
    for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
      meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,3), 0) << "  " << nodes(cellToNodeMap(i,3), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
    }
    meshfile.close();

    Teuchos::RCP<Intrepid::Basis_HGRAD_QUAD_C1_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrQ1 =
      Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<RealT, Intrepid::FieldContainer<RealT> >);

    Teuchos::RCP<Intrepid::Basis_HGRAD_QUAD_C2_FEM<RealT, Intrepid::FieldContainer<RealT> > > basisPtrQ2 =
      Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<RealT, Intrepid::FieldContainer<RealT> >);

    std::vector<Teuchos::RCP<Intrepid::Basis<RealT, Intrepid::FieldContainer<RealT> > > > basisPtrs(3, Teuchos::null);
    basisPtrs[0] = basisPtrQ2;
    basisPtrs[1] = basisPtrQ2;
    basisPtrs[2] = basisPtrQ1;

    Teuchos::RCP<MeshManager<RealT> > meshmgrPtr = Teuchos::rcpFromRef(meshmgr);

    DofManager<RealT> dofmgr(meshmgrPtr, basisPtrs);

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
