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

/*! \file  test_04.cpp
    \brief Unit test for the FE class on quads.
*/

#include "ROL_Types.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/dofmanager.hpp"
#include "../TOOLS/fe.hpp"

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
    std::string filename = "input_04.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshManager_Rectangle<RealT> meshmgr(*parlist);
    Teuchos::RCP<Intrepid::FieldContainer<RealT> > nodesPtr = meshmgr.getNodes();
    Teuchos::RCP<Intrepid::FieldContainer<int> >   cellToNodeMapPtr = meshmgr.getCellToNodeMap();

    Intrepid::FieldContainer<RealT> &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl;
    *outStream << "Number of edges = " << meshmgr.getNumEdges() << std::endl;

    Teuchos::RCP<MeshManager<RealT> > meshmgrPtr = Teuchos::rcpFromRef(meshmgr);

    // Basis.
    Teuchos::RCP<Intrepid::Basis<RealT, Intrepid::FieldContainer<RealT> > >    basis =
      Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<RealT, Intrepid::FieldContainer<RealT> >);
    // Cubature.
    Teuchos::RCP<Intrepid::Cubature<RealT, Intrepid::FieldContainer<RealT> > > cubature;
    shards::CellTopology cellType = basis->getBaseCellTopology();                        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<RealT> cubFactory;                                  // create cubature factory
    int cubDegree = 4;                                                                   // set cubature degree, e.g., 2
    cubature = cubFactory.create(cellType, cubDegree);                                   // create default cubature
    // Cell nodes.
    Intrepid::FieldContainer<RealT> cellNodes(meshmgr.getNumCells(), cellType.getNodeCount(), cellType.getDimension());
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      for (unsigned j=0; j<cellType.getNodeCount(); ++j) {
        for (unsigned k=0; k<cellType.getDimension(); ++k) {
          cellNodes(i,j,k) = nodes(cellToNodeMap(i,j), k);
        }
      }
    }
    // FE object.
    FE<RealT> fe(Teuchos::rcpFromRef(cellNodes), basis, cubature);

    // Check integration.
    Intrepid::FieldContainer<RealT> feVals1(meshmgr.getNumCells(), cubature->getNumPoints());
    Intrepid::FieldContainer<RealT> feVals2(meshmgr.getNumCells(), cubature->getNumPoints());
    Intrepid::FieldContainer<RealT> feGrads1(meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    Intrepid::FieldContainer<RealT> feGrads2(meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    Intrepid::FieldContainer<RealT> feCoeffs1(meshmgr.getNumCells(), basis->getCardinality());
    Intrepid::FieldContainer<RealT> feCoeffs2(meshmgr.getNumCells(), basis->getCardinality());
    Intrepid::FieldContainer<RealT> feIntegral(meshmgr.getNumCells());
    // values
    feCoeffs1.initialize(static_cast<RealT>(1));
    feCoeffs2.initialize(static_cast<RealT>(2));
    fe.evaluateValue(Teuchos::rcpFromRef(feVals1), Teuchos::rcpFromRef(feCoeffs1));
    fe.evaluateValue(Teuchos::rcpFromRef(feVals2), Teuchos::rcpFromRef(feCoeffs2));
    fe.computeIntegral(Teuchos::rcpFromRef(feIntegral), Teuchos::rcpFromRef(feVals1), Teuchos::rcpFromRef(feVals2));
    RealT valval(0);
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      valval += feIntegral(i);
    }
    if (std::abs(valval - static_cast<RealT>(1000)) > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag = -1;
    }
    // gradients
    Intrepid::FieldContainer<RealT> dofCoords(meshmgr.getNumCells(), basis->getCardinality(), cubature->getDimension());
    fe.computeDofCoords(Teuchos::rcpFromRef(dofCoords), Teuchos::rcpFromRef(cellNodes));
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      for (int j=0; j<basis->getCardinality(); ++j) {
        RealT x = dofCoords(i,j,0);
        RealT y = dofCoords(i,j,1);
        feCoeffs1(i,j) = x + y;
        feCoeffs2(i,j) = static_cast<RealT>(2)*x + static_cast<RealT>(3)*y;
      }
    }
    fe.evaluateGradient(Teuchos::rcpFromRef(feGrads1), Teuchos::rcpFromRef(feCoeffs1));
    fe.evaluateGradient(Teuchos::rcpFromRef(feGrads2), Teuchos::rcpFromRef(feCoeffs2));
    fe.computeIntegral(Teuchos::rcpFromRef(feIntegral), Teuchos::rcpFromRef(feGrads1), Teuchos::rcpFromRef(feGrads2));
    RealT gradgrad(0);
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      gradgrad += feIntegral(i);
    }
    if (std::abs(gradgrad - static_cast<RealT>(2500)) > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag = -1;
    }
    *outStream << std::setprecision(14) << "\nvalvalIntegral=" << valval << "  " << "gradgradIntegral=" << gradgrad << std::endl;

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
