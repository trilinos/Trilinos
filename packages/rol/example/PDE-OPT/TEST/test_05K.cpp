// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_04.cpp
    \brief Unit test for the FE_CURL class.
*/

#include "ROL_Types.hpp"

#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/dofmanagerK.hpp"
#include "../TOOLS/feK_curl.hpp"

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
    std::string filename = "input_05.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize mesh / degree-of-freedom manager. ***/
    MeshManager_Brick<RealT, DeviceT> meshmgr(*parlist);
    scalar_view nodes = meshmgr.getNodes();
    int_view cellToNodeMap = meshmgr.getCellToNodeMap();

    *outStream << "Number of nodes = " << meshmgr.getNumNodes() << std::endl;
    *outStream << "Number of cells = " << meshmgr.getNumCells() << std::endl;
    *outStream << "Number of edges = " << meshmgr.getNumEdges() << std::endl;

    ROL::Ptr<MeshManager<RealT, DeviceT>> meshmgrPtr = ROL::makePtrFromRef(meshmgr);

    // Basis.
    basis_ptr basis = ROL::makePtr<Intrepid2::Basis_HCURL_HEX_I1_FEM<DeviceT, RealT, RealT>>();
    // Cubature.
    ROL::Ptr<Intrepid2::Cubature<DeviceT, RealT, RealT>> cubature;
    shards::CellTopology cellType = basis->getBaseCellTopology();                        // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                        // create cubature factory
    int cubDegree = 4;                                                                   // set cubature degree, e.g., 2
    cubature = cubFactory.create<DeviceT, RealT, RealT>(cellType, cubDegree);            // create default cubature
    // Cell nodes.
    scalar_view cellNodes("cellNodes", meshmgr.getNumCells(), cellType.getNodeCount(), cellType.getDimension());
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      for (unsigned j=0; j<cellType.getNodeCount(); ++j) {
        for (unsigned k=0; k<cellType.getDimension(); ++k) {
          cellNodes(i,j,k) = nodes(cellToNodeMap(i,j), k);
        }
      }
    }
    // FE object.
    FE_CURL<RealT, DeviceT> fe(cellNodes, basis, cubature);

    // Check integration.
    scalar_view feVals1("feVals1", meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    scalar_view feVals2("feVals2", meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    scalar_view feCurls1("feCurls1", meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    scalar_view feCurls2("feCurls2", meshmgr.getNumCells(), cubature->getNumPoints(), cubature->getDimension());
    scalar_view feCoeffs1("feCoeffs1", meshmgr.getNumCells(), basis->getCardinality());
    scalar_view feCoeffs2("feCoeffs2", meshmgr.getNumCells(), basis->getCardinality());
    scalar_view feIntegral("feIntegral", meshmgr.getNumCells());
    // values
    Kokkos::deep_copy(feCoeffs1,static_cast<RealT>(0));
    Kokkos::deep_copy(feCoeffs2,static_cast<RealT>(0));
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      feCoeffs1(i, 0) = static_cast<RealT>(1);
      feCoeffs2(i, 0) = static_cast<RealT>(1);
    }
    fe.evaluateValue(feVals1, feCoeffs1);
    fe.evaluateValue(feVals2, feCoeffs2);
    fe.computeIntegral(feIntegral, feVals1, feVals2);
    RealT valval(0);
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      valval += feIntegral(i);
    }
    if (std::abs(valval - static_cast<RealT>(4)/static_cast<RealT>(9)) > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag = -1;
    }
    // curls
    fe.evaluateCurl(feCurls1, feCoeffs1);
    fe.evaluateCurl(feCurls2, feCoeffs2);
    fe.computeIntegral(feIntegral, feCurls1, feCurls2);
    RealT curlcurl(0);
    for (int i=0; i<meshmgr.getNumCells(); ++i) {
      curlcurl += feIntegral(i);
    }
    if (std::abs(curlcurl - static_cast<RealT>(8)/static_cast<RealT>(3)) > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag = -1;
    }
    *outStream << std::setprecision(14) << "\nvalvalIntegral=" << valval << "  "
                                        << "curlcurlIntegral=" << curlcurl << std::endl;

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
