// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
//
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
  \brief local L2 projection of a function onto continuous finite element space

  We consider HGrad, HCurl, HDiv and HVol spaces and Hexahedral, Tetrahedral, Quadrilateral and Triangular shapes.

  The mesh topology and global numbering of the Degrees of Freedom (DoFs) are computed using Panzer Dof Manager package

  We project polynomials of degree at most 2 that belong to the finite element spaces and
  we verify that the L2 error is order of round off.
  We also verify that degrees of freedoms shared by different elements are consistent across elements.

  We also project the same polynomials on the boundary sides, using the finte element spaces restricted to the boundary,
  an we verify that we obtain the same basis coefficients as the volume porjection at the boundary.

  Command line arguments:
  --space: the finite element spaces: "HGrad", "HCurl", "HDiv" and "HVOL"
  --shape: the elemet shape: "Hexahedron", "Tetrahedron", "Quadrilateral" and "Triangle".
  --nx, --ny, --nz:  the number of owned 1d elements in the axis directions. The total number of owned hexas is given by nx \times ny \times \nz
  --px, --py, --pz:  the number of processes in the axis directions. The total numebr of processes is px \times py \times pz.
  --verbose: whether to print out info to screen; default: 1
  --timigs-file: the file where to print timings; default: "", output redirected to screen

  \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_Kernels.hpp"

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "../cartesian_topology/CartesianConnManager.hpp"

#include <Tpetra_Export.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_FECrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_Assembly_Helpers.hpp>
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <array>
#include <set>
#include <random>
#include <algorithm>

namespace Discretization {

namespace Example {

// ************************************ Define Analytic functions **************************************

struct Fun {
  double
  KOKKOS_INLINE_FUNCTION
  operator()(const double& x, const double& y, const double& z, const int& comp) const {

    double f0 = y;
    double f1 = 0;//std::pow(y, degree-1);
    double f2 = 1;//std::pow(x+z, degree-1);

    switch (comp) {
    case 0:
      return f0 + a*x + (a1*z-a2*y);
    case 1:
      return f1 + a*y + (a2*x-a0*z);
    case 2:
      return f2 + a*z + (a0*y-a1*x);
    default:
      return 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  Fun(const Intrepid2::EFunctionSpace& space) {
    a=a0=a1=a2=0;
    if(space != Intrepid2::FUNCTION_SPACE_HCURL)
      a=1;
    if(space != Intrepid2::FUNCTION_SPACE_HDIV){
      a0 = 2; a1 = -1; a2 = 3;
    }
  }

  double a,a0,a1,a2;
};

enum ElemShape {HEX, TET, QUAD, TRI};


template<typename ValueType, typename DeviceSpaceType>
int feProjection(int argc, char *argv[]) {

  using HostSpaceType = typename
      Kokkos::Impl::HostMirror<DeviceSpaceType>::Space::execution_space;
  using DynRankView = Kokkos::DynRankView<ValueType,DeviceSpaceType>;
  using DynRankViewHost = Kokkos::DynRankView<ValueType,HostSpaceType>;

  using map_t = Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal>;

  using local_ordinal_t = typename map_t::local_ordinal_type;
  using global_ordinal_t = typename map_t::global_ordinal_type;
  using scalar_t = ValueType;

  using DynRankViewGId = Kokkos::DynRankView<global_ordinal_t,DeviceSpaceType>;

  using ct = Intrepid2::CellTools<DeviceSpaceType>;
  using ots = Intrepid2::OrientationTools<DeviceSpaceType>;
  using pts = Intrepid2::ProjectionTools<DeviceSpaceType>;
  using ProjectionStruct = Intrepid2::ProjectionStruct<DeviceSpaceType,scalar_t>;
  using rst = Intrepid2::RealSpaceTools<DeviceSpaceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceSpaceType>;

  using CellTopology = shards::CellTopology;

  int errorFlag = 0;


  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);

  //output stream/file
  Teuchos::RCP<Teuchos::FancyOStream> outStream;
  std::string timingsFile = "";

  try {


    // ************************************ GET INPUTS **************************************
    local_ordinal_t dim = 3;

    // parse command line arguments
    int degree = 2;
    local_ordinal_t nx = 2;
    local_ordinal_t ny = nx;
    local_ordinal_t nz = (dim == 3) ? nx :1;
    int np   = comm.getSize(); // number of processors
    int px = (dim == 2) ? std::sqrt(np) : std::cbrt(np); while(np%px!=0) --px;
    int py = (dim == 2) ? np/px : std::sqrt(np/px); while(np%py!=0) --py;
    int pz = np/(px*py);
    constexpr int bx=1, by=1, bz=1;  //blocks on each process. Here we assume there is only one block per process
    int verbose = 1;
    std::string shape("Hexahedron"), space("HGrad");

    //side of the brick domain containing the elements and its normal
    int sideFlag = 3;
    double sideNormal[3] = {-1.0,0.0,0.0};


    Teuchos::CommandLineProcessor clp;
    clp.setOption("shape",&shape);
    clp.setOption("space",&space);
    clp.setOption("nx",&nx);
    clp.setOption("ny",&ny);
    clp.setOption("nz",&nz);
    clp.setOption("px",&px);
    clp.setOption("py",&py);
    clp.setOption("pz",&pz);
    clp.setOption("basis-degree",&degree);
    clp.setOption("timings-file",&timingsFile);
    clp.setOption("verbose", &verbose);
    auto cmdResult = clp.parse(argc,argv);
    if(cmdResult!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      clp.printHelpMessage(argv[0], std::cout);
      errorFlag++;
    }

    ElemShape eShape;
    Teuchos::RCP<CellTopology> cellTopoPtr = Teuchos::rcp(new CellTopology);

    //brickBasis: used to compute the coordinates of the mesh vertices
    Teuchos::RCP< Intrepid2::Basis<HostSpaceType, scalar_t,scalar_t> > brickBasis;

    //linearBasis: used to compute the coordinates of the quadrature points for computing the L2 error
    Teuchos::RCP< Intrepid2::Basis<DeviceSpaceType, scalar_t,scalar_t> > linearBasis;

    if(shape == "Hexahedron") {
      cellTopoPtr = Teuchos::rcp(new CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >()));
      brickBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<HostSpaceType,scalar_t,scalar_t>);
      linearBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,scalar_t,scalar_t>);
      eShape = HEX;
      dim = 3;
    } else if (shape == "Tetrahedron") {
      cellTopoPtr = Teuchos::rcp(new CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >()));
      brickBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<HostSpaceType,scalar_t,scalar_t>);
      linearBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,scalar_t,scalar_t>);
      eShape = TET;
      dim = 3;
    } else if (shape == "Quadrilateral") {
      cellTopoPtr = Teuchos::rcp(new CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >()));
      brickBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<HostSpaceType,scalar_t,scalar_t>);
      linearBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,scalar_t,scalar_t>);
      eShape = QUAD;
      dim = 2;
    } else if (shape == "Triangle") {
      cellTopoPtr = Teuchos::rcp(new CellTopology(shards::getCellTopologyData<shards::Triangle<3> >()));
      brickBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<HostSpaceType,scalar_t,scalar_t>);
      linearBasis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceSpaceType,scalar_t,scalar_t>);
      eShape = TRI;
      dim = 2;
    } else {
      TEUCHOS_TEST_FOR_TERMINATION(true,
          "test_fe_projection.hpp: element shape " << shape << " is not supported.\n" <<
          "Supported shapes are: Hexahedron, Tetrahedron, Quadrilateral and Triangle." << std::endl);
    }
    int basisDimension;
    Teuchos::RCP< Intrepid2::Basis<DeviceSpaceType, scalar_t,scalar_t> > basis;
    using CG_DNBasis = Intrepid2::NodalBasisFamily<DeviceSpaceType,scalar_t,scalar_t>;
    Intrepid2::EFunctionSpace functionSpace;
    if(space == "HGrad") {
      functionSpace = Intrepid2::FUNCTION_SPACE_HGRAD;
      switch (eShape) {
      case HEX:
        basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_HEX(degree));
        break;
      case TET:
        basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_TET(degree));
        break;
      case QUAD:
        basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_QUAD(degree));
        break;
      case TRI:
        basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_TRI(degree));
      }
      basisDimension = 1;
    } else if (space == "HCurl") {
      functionSpace = Intrepid2::FUNCTION_SPACE_HCURL;
      switch (eShape) {
      case HEX:
        basis = Teuchos::rcp(new typename CG_DNBasis::HCURL_HEX(degree));
        break;
      case TET:
        basis = Teuchos::rcp(new typename CG_DNBasis::HCURL_TET(degree));
        break;
      case QUAD:
        basis = Teuchos::rcp(new typename CG_DNBasis::HCURL_QUAD(degree));
        break;
      case TRI:
        basis = Teuchos::rcp(new typename CG_DNBasis::HCURL_TRI(degree));
      }
      basisDimension = dim;
    } else if (space == "HDiv") {
      functionSpace = Intrepid2::FUNCTION_SPACE_HDIV;
      switch (eShape) {
      case HEX:
        basis = Teuchos::rcp(new typename CG_DNBasis::HDIV_HEX(degree));
        break;
      case TET:
        basis = Teuchos::rcp(new typename CG_DNBasis::HDIV_TET(degree));
        break;
      case QUAD:
        basis = Teuchos::rcp(new typename CG_DNBasis::HDIV_QUAD(degree));
        break;
      case TRI:
        basis = Teuchos::rcp(new typename CG_DNBasis::HDIV_TRI(degree));
      }
      basisDimension = dim;
    } else if (space == "HVol") {
      functionSpace = Intrepid2::FUNCTION_SPACE_HVOL;
      switch (eShape) {
      case HEX:
        basis = Teuchos::rcp(new typename CG_DNBasis::HVOL_HEX(degree));
        break;
      case TET:
        basis = Teuchos::rcp(new typename CG_DNBasis::HVOL_TET(degree));
        break;
      case QUAD:
        basis = Teuchos::rcp(new typename CG_DNBasis::HVOL_QUAD(degree));
        break;
      case TRI:
        basis = Teuchos::rcp(new typename CG_DNBasis::HVOL_TRI(degree));
      }
      basisDimension = 1;
    } else {
      TEUCHOS_TEST_FOR_TERMINATION(true,
          "test_fe_projection.hpp: function space " << space <<  " is not supported.\n" <<
          "Supported spaces are: HGrad, HCurl, HDiv and HVol." << std::endl);
    }

    if(dim == 2) {
      nz = 1;
      px = std::sqrt(np); while(np%px!=0) --px;
      py = np/px; while(np%py!=0) --py;
      pz = 1;
    }

    cmdResult = clp.parse(argc,argv);
    int targetCubDegree = 2; //order of exactness the quadrature rule to integrate the forcing term

    if(cmdResult!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      clp.printHelpMessage(argv[0], std::cout);
      errorFlag++;
    }

    outStream = ((comm.getRank () == 0) && verbose) ?
        getFancyOStream(Teuchos::rcpFromRef (std::cout)) :
        getFancyOStream(Teuchos::rcp (new Teuchos::oblackholestream ()));

    *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
    *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
    *outStream << "\n";


    auto meshTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Mesh Generation")));

    // *********************************** MESH TOPOLOGY **********************************

    // Get dimensions
    int numNodesPerElem = cellTopoPtr->getNodeCount();

    // build the topology
    auto connManager = Teuchos::rcp(new panzer::unit_test::CartesianConnManager);
    if(dim == 3) {
      connManager->initialize(comm,
          global_ordinal_t(nx*px),
          global_ordinal_t(ny*py),
          global_ordinal_t(nz*pz),
          px,py,pz,bx,by,bz,
          *cellTopoPtr);
    } else {
      connManager->initialize(comm,
          global_ordinal_t(nx*px),
          global_ordinal_t(ny*py),
          px,py,bx,by,
          *cellTopoPtr);
    }


    // *********************************** COMPUTE GLOBAL IDs OF VERTICES AND DOFs  ************************************

    // build the dof manager, and assocaite with the topology
    auto dofManager = Teuchos::rcp(new panzer::DOFManager);
    dofManager->setConnManager(connManager,*comm.getRawMpiComm());

    dofManager->setOrientationsRequired(basis->requireOrientation());

    auto basisCardinality = basis->getCardinality();
    Teuchos::RCP<panzer::Intrepid2FieldPattern> fePattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
    std::string blockId = (dim == 2) ? "block-0_0" : "block-0_0_0";
    dofManager->addField(blockId,fePattern);

    // try to get them all synced up
    comm.barrier();

    dofManager->buildGlobalUnknowns();

    local_ordinal_t numOwnedBricks = (dim == 2) ? nx*ny : nx*ny*nz;
    local_ordinal_t numOwnedElems = connManager->numSubElemsPerBrickElement()*numOwnedBricks;

    *outStream << "Number of cells on each processor nx X ny X nz: " << numOwnedElems << "\n";
    *outStream << "    nx" << "   ny" << "   nz\n";
    *outStream << std::setw(5) << nx <<
        std::setw(5) << ny <<
        std::setw(5) << nz << "\n\n";
    *outStream << "Number of processors px X py X pz: " << px*py*pz << "\n";
    *outStream << "    px" << "   py" << "   pz\n";
    *outStream << std::setw(5) << px <<
        std::setw(5) << py <<
        std::setw(5) << pz << "\n\n";

    global_ordinal_t totalNumCells = numOwnedElems*px*py*pz;
    *outStream << "Total number of cells: " << totalNumCells << ", number of DoFs per element: " << basisCardinality << "\n";

    // Print mesh information

    // Cube
    scalar_t leftX = -1.0, rightX = 1.0;
    scalar_t leftY = -1.0, rightY = 1.0;
    scalar_t leftZ = -1.0, rightZ = 1.0;

    // Mesh spacing
    scalar_t hx = (rightX-leftX)/((scalar_t)(nx*px*bx));
    scalar_t hy = (rightY-leftY)/((scalar_t)(ny*py*by));
    scalar_t hz = (rightZ-leftZ)/((scalar_t)(nz*pz*bz));


    // *********************************** COMPUTE COORDINATES OF PHYSICAL VERTICES  ************************************

    // Get coordinates of physical vertices
    DynRankView physVertexes("physVertexes", numOwnedElems, numNodesPerElem, dim);
    {
      auto physVertexesHost = Kokkos::create_mirror_view(physVertexes);

      //Intrepid2::Basis_HGRAD_QUAD_C1_FEM<HostSpaceType,scalar_t,scalar_t> brickBasis;
      int numNodesPerBrick = brickBasis->getCardinality();
      DynRankViewHost refVerticesHexa("refVerticesHexa", numNodesPerBrick, dim);
      DynRankViewHost physVerticesHexa("physVerticesHexa", numNodesPerBrick, dim);
      brickBasis->getDofCoords(refVerticesHexa);

      auto elemTriplet = connManager->getMyBrickElementsTriplet();
      double h[3] = {hx, hy, hz};

      for(int i=0; i<numOwnedBricks; ++i) {
        elemTriplet =  connManager->computeLocalBrickElementGlobalTriplet(i,connManager->getMyBrickElementsTriplet(),connManager->getMyBrickOffsetTriplet());
        double offset[3] = {leftX + elemTriplet.x*hx+hx/2, leftY +elemTriplet.y*hy+hy/2, leftX +elemTriplet.z*hz+hz/2};
        for(int j=0; j<numNodesPerBrick; ++j) {
          for(int k=0; k<dim; ++k)
            physVerticesHexa(j,k) = offset[k]+h[k]/2.0*refVerticesHexa(j,k);
        }
        for(int subElement=0; subElement<connManager->numSubElemsPerBrickElement(); ++subElement) {
          int cell = i*connManager->numSubElemsPerBrickElement()+subElement;
          for(int n=0; n<numNodesPerElem; ++n) {
            int j = connManager->getLocalBrickNodeFromSubElemNode(subElement,n);
            for(int k=0; k<dim; ++k) {
              physVertexesHost(cell,n,k) = physVerticesHexa(j,k);
            }
          }
        }
      }
      Kokkos::deep_copy(physVertexes, physVertexesHost);
    }

    meshTimer = Teuchos::null;


    // *********************************** COMPUTE ELEMENTS' ORIENTATION BASED ON GLOBAL IDs  ************************************

    auto orientationTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Building Orientations")));

    Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpaceType> elemOrts("elemOrts", numOwnedElems);

    if(basis->requireOrientation()) {
      //compute global ids of element vertices
      DynRankViewGId elemNodesGID("elemNodesGID", numOwnedElems, numNodesPerElem);

      {
        auto elemNodesGIDHost = Kokkos::create_mirror_view(elemNodesGID);
        for(int i=0; i<numOwnedElems; ++i) {
          const auto GIDs = connManager->getConnectivity(i);
          for(int j=0; j<numNodesPerElem; ++j) {
            elemNodesGIDHost(i,j) = GIDs[j];
          }
        }
        Kokkos::deep_copy(elemNodesGID, elemNodesGIDHost);

      }

      // compute orientations for cells (one time computation)
      ots::getOrientation(elemOrts, elemNodesGID, *cellTopoPtr);
    }

    orientationTimer = Teuchos::null;


    // ************************************ L2 Projection **************************************

    auto projectionTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Computing L2 Projections")));

    auto globalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(dofManager);
    std::vector<global_ordinal_t> ownedIndices, ownedAndGhostedIndices;
    globalIndexer->getOwnedIndices(ownedIndices);
    Teuchos::RCP<const map_t> ownedMap = Teuchos::rcp(new map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedIndices,0,Teuchos::rcpFromRef(comm)));

    *outStream << "Total number of DoFs: " << ownedMap->getGlobalNumElements() << ", number of owned DoFs: " << ownedMap->getLocalNumElements() << "\n";

    Teuchos::Array<global_ordinal_t> globalIdsInRow(basisCardinality);
    blockId = (dim == 2) ? "eblock-0_0" : "eblock-0_0_0";
    auto elmtOffsetKokkos = dofManager->getGIDFieldOffsetsKokkos(blockId,0);

    DynRankView basisCoeffsL2Proj("basisCoeffsL2Proj", numOwnedElems, basisCardinality);
    {
      ProjectionStruct projStruct;
      projStruct.createL2ProjectionStruct(basis.get(), targetCubDegree);

      auto evaluationPoints = projStruct.getAllEvalPoints();
      int numPoints     = evaluationPoints.extent(0);

      DynRankView refTargetAtEvalPoints, physTargetAtEvalPoints;
      if(functionSpace == Intrepid2::FUNCTION_SPACE_HCURL || functionSpace == Intrepid2::FUNCTION_SPACE_HDIV) {
        refTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints, dim);
        physTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints, dim);
      } else {
        refTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints);
        physTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints);
      }

      DynRankView physEvalPoints("physEvalPoints", numOwnedElems, numPoints, dim);
      ct::mapToPhysicalFrame(physEvalPoints,evaluationPoints,physVertexes,basis->getBaseCellTopology());

      //transform the target function and its derivative to the reference element (inverse of pullback operator)
      DynRankView jacobian("jacobian", numOwnedElems, numPoints, dim, dim);
      DynRankView jacobian_det("jacobian_det", numOwnedElems, numPoints);
      DynRankView jacobian_inv("jacobian_inv", numOwnedElems, numPoints, dim, dim);
      ct::setJacobian(jacobian, evaluationPoints, physVertexes, *cellTopoPtr);
      ct::setJacobianDet (jacobian_det, jacobian);
      ct::setJacobianInv (jacobian_inv, jacobian);


      Kokkos::parallel_for(Kokkos::RangePolicy<typename DeviceSpaceType::execution_space>(0,numOwnedElems),
          KOKKOS_LAMBDA (const int &ic) {
        Fun fun(functionSpace);
        for(int i=0;i<numPoints;i++) {
          auto x = physEvalPoints(ic,i,0), y = physEvalPoints(ic,i,1);
          auto z = (dim==2) ? 0.0 : physEvalPoints(ic,i,2);
          for(int d=0;d<basisDimension;d++)
            physTargetAtEvalPoints(ic,i,d) = fun(x,y,z,d);
        }
      });

      switch (functionSpace) {
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        fst::mapHGradDataFromPhysToRef(refTargetAtEvalPoints,physTargetAtEvalPoints);
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
        fst::mapHCurlDataFromPhysToRef(refTargetAtEvalPoints,jacobian,physTargetAtEvalPoints);
        break;
      case Intrepid2::FUNCTION_SPACE_HDIV:
        fst::mapHDivDataFromPhysToRef(refTargetAtEvalPoints,jacobian_inv,jacobian_det,physTargetAtEvalPoints);
        break;
      case Intrepid2::FUNCTION_SPACE_HVOL:
        fst::mapHVolDataFromPhysToRef(refTargetAtEvalPoints,jacobian_det,physTargetAtEvalPoints);
        break;
      default: {}
      }

      pts::getL2BasisCoeffs(basisCoeffsL2Proj,
          refTargetAtEvalPoints,
          elemOrts,
          basis.get(),
          &projStruct);
    }


    projectionTimer = Teuchos::null;


    // ************************************ CODE VERIFICATION **************************************

    // Make sure that DoFs are consistent on shared vertices/edges/faces
    {
      std::map<global_ordinal_t,scalar_t> mapL2Proj;
      std::map<global_ordinal_t,int> mapCell;
      std::map<global_ordinal_t, typename Intrepid2::Basis<HostSpaceType, scalar_t,scalar_t>::OrdinalTypeArrayStride1DHost> mapDofTag;
      std::vector<global_ordinal_t> elementGIDs(basisCardinality);
      Teuchos::TimeMonitor vTimer1 =  *Teuchos::TimeMonitor::getNewTimer("Verification, assemble solution");
      auto basisCoeffsL2ProjHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2Proj);
      auto elmtOffsetKokkosHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), elmtOffsetKokkos);

      for(int elemId=0; elemId<numOwnedElems; elemId++)
      {
        dofManager->getElementGIDs(elemId, elementGIDs);
        for(int nodeId=0; nodeId<basisCardinality; nodeId++)
        {
          global_ordinal_t gid = elementGIDs[elmtOffsetKokkosHost(nodeId)];

          if(ownedMap->isNodeGlobalElement(gid)) {
            auto it = mapL2Proj.find(gid);
            if (it==mapL2Proj.end()) {
              mapL2Proj.insert(std::make_pair(gid,basisCoeffsL2ProjHost(elemId, nodeId)));
              mapDofTag.insert(std::make_pair(gid,basis->getDofTag(nodeId)));
              mapCell.insert(std::make_pair(gid,elemId));
            }
            else if(std::abs(it->second-basisCoeffsL2ProjHost(elemId, nodeId))>1e-10) {
              auto dofTagNew = basis->getDofTag(nodeId);
              auto dofTagStored = mapDofTag.find(gid)->second;

              std::cout << "ERROR: DoFs shared by cells are not consistent \n"
                  "basisL2Proj(" << gid << "):" << it->second << " " << basisCoeffsL2ProjHost(elemId, nodeId) <<"\n"
                  "cell " << mapCell.find(gid)->second << ": " << dofTagStored(0) << ", " << dofTagStored(1) << ", " << dofTagStored(2) <<
                  " | cell " << elemId << ": " << dofTagNew(0) << ", " << dofTagNew(1) << ", " << dofTagNew(2)
                  << std::endl;
              errorFlag++;
            }
          }
        }
      }
    }

    // Compute the L2 norm of the diff between the function and its projection on the FE space
    Intrepid2::DefaultCubatureFactory cub_factory;
    int cubDegreeL2 = 2*std::max(targetCubDegree,degree);
    auto cell_cub = cub_factory.create<DeviceSpaceType, scalar_t, scalar_t>(cellTopoPtr->getBaseKey(), cubDegreeL2);
    int numRefCoords = cell_cub->getNumPoints();
    DynRankView refPoints("refPoints", numRefCoords, dim);
    DynRankView weights("weights", numRefCoords);
    cell_cub->getCubature(refPoints, weights);

    //Compute physical Dof Coordinates and Reference coordinates
    DynRankView physRefCoords("physRefCoords", numOwnedElems, numRefCoords, dim);
    {
      DynRankView linearBasisValuesAtRefCoords("linearBasisValuesAtRefCoords", numNodesPerElem, numRefCoords);
      linearBasis->getValues(linearBasisValuesAtRefCoords, refPoints);
      Kokkos::fence();
      Kokkos::MDRangePolicy<typename DeviceSpaceType::execution_space,Kokkos::Rank<3>> rangePolicy({0,0,0},{numOwnedElems, numRefCoords, dim});
      Kokkos::parallel_for(rangePolicy, KOKKOS_LAMBDA (const int &i, const int &j, const int &d) {
          for(int k=0; k<numNodesPerElem; ++k)
            physRefCoords(i,j,d) += physVertexes(i,k,d)*linearBasisValuesAtRefCoords(k,j);
      });
      Kokkos::fence();
    }

    // Compute basis function and compute FE solution in physical space
    DynRankView basisValuesAtRefCoordsOriented;
    DynRankView basisValuesAtRefCoordsCells;
    DynRankView basisValuesAtRefCoords;

    if(functionSpace == Intrepid2::FUNCTION_SPACE_HCURL || functionSpace == Intrepid2::FUNCTION_SPACE_HDIV) {
      basisValuesAtRefCoordsOriented = DynRankView("basisValuesAtRefCoordsOriented", numOwnedElems, basisCardinality, numRefCoords, dim);
      basisValuesAtRefCoordsCells = DynRankView("inValues", numOwnedElems, basisCardinality, numRefCoords, dim);
      basisValuesAtRefCoords = DynRankView("basisValuesAtRefCoords", basisCardinality, numRefCoords, dim);
    } else {
      basisValuesAtRefCoordsOriented = DynRankView("basisValuesAtRefCoordsOriented", numOwnedElems, basisCardinality, numRefCoords);
      basisValuesAtRefCoordsCells = DynRankView("inValues", numOwnedElems, basisCardinality, numRefCoords);
      basisValuesAtRefCoords = DynRankView("basisValuesAtRefCoords", basisCardinality, numRefCoords);
    }

    basis->getValues(basisValuesAtRefCoords, refPoints);
    rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

    // modify basis values to account for orientations
    ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
        basisValuesAtRefCoordsCells,
        elemOrts,
        basis.getRawPtr());

    // transform basis values to the reference element (pullback)
    DynRankView jacobianAtRefCoords("jacobianAtRefCoords", numOwnedElems, numRefCoords, dim, dim);
    ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, *cellTopoPtr);
    DynRankView jacobianAtRefCoords_det("jacobianAtRefCoords_det", numOwnedElems, numRefCoords);
    ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);

    DynRankView transformedBasisValuesAtRefCoordsOriented("transformedBasisValuesAtRefCoordsOriented", numOwnedElems, basisCardinality, numRefCoords, basisDimension);
    switch (functionSpace) {
    case Intrepid2::FUNCTION_SPACE_HGRAD:
    {
      fst::HGRADtransformVALUE(
          Kokkos::subview(transformedBasisValuesAtRefCoordsOriented,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0),
          basisValuesAtRefCoordsOriented);
    } break;
    case Intrepid2::FUNCTION_SPACE_HCURL:
    {
      DynRankView jacobianInvAtRefCoords("jacobianInvAtRefCoords", numOwnedElems, numRefCoords, dim, dim);
      ct::setJacobianInv (jacobianInvAtRefCoords, jacobianAtRefCoords);
      fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
          jacobianInvAtRefCoords,
          basisValuesAtRefCoordsOriented);
    } break;
    case Intrepid2::FUNCTION_SPACE_HDIV:
    {
      fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
          jacobianAtRefCoords,
          jacobianAtRefCoords_det,
          basisValuesAtRefCoordsOriented);
    } break;
    case Intrepid2::FUNCTION_SPACE_HVOL:
    {
      fst::HVOLtransformVALUE(
          Kokkos::subview(transformedBasisValuesAtRefCoordsOriented,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0),
          jacobianAtRefCoords_det,
          basisValuesAtRefCoordsOriented);
    } break;
    default:
    {
      TEUCHOS_TEST_FOR_TERMINATION(true,
          "test_fe_projection.hpp: function space not supported" << std::endl);
    }
    }


    // compute L2 error
    ValueType norm2(0);
    DynRankView projectedFunAtRefCoords("projectedFunAtRefCoords", numOwnedElems, numRefCoords, basisDimension);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<typename DeviceSpaceType::execution_space>(0,numOwnedElems),
        KOKKOS_LAMBDA (const int &ic, double &norm2Update) {
      Fun fun(functionSpace);
      for(int j=0; j<numRefCoords; ++j) {
        for(int k=0; k<basisCardinality; ++k) {
          for(int d=0; d<basisDimension; ++d)
            projectedFunAtRefCoords(ic,j,d) += basisCoeffsL2Proj(ic,k)*transformedBasisValuesAtRefCoordsOriented(ic,k,j,d);
        }

        auto x = physRefCoords(ic,j,0), y = physRefCoords(ic,j,1);
        auto z = (dim==2) ? 0.0 : physRefCoords(ic,j,2);
        for(int d=0; d<basisDimension; ++d) {
          scalar_t funAtRefCoords = fun(x,y,z,d);
          norm2Update += (funAtRefCoords - projectedFunAtRefCoords(ic,j,d))*
              (funAtRefCoords - projectedFunAtRefCoords(ic,j,d))*
              weights(j)*jacobianAtRefCoords_det(ic,j);
        }
      }
    },norm2);

    double totalNorm2=0;
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&norm2,&totalNorm2);

    *outStream << "L2 error: " << std::sqrt(totalNorm2) << std::endl;
    if(std::sqrt(totalNorm2)>1e-12) {
      std::cout << "ERROR: The projection error should be zero up to round off errors" <<std::endl;
      errorFlag++;
    }




    // ********************************** Side L2 Projection ***********************************
    if(functionSpace != Intrepid2::FUNCTION_SPACE_HVOL) {

      int sideDim = dim-1;

      //note this would not work with elements like wedges that have different side basis
      auto sideBasisPtr = basis->getSubCellRefBasis(sideDim,0);
      auto sideFunctionSpace = sideBasisPtr->getFunctionSpace();
      int numNodesPerSide = sideBasisPtr->getBaseCellTopology().getNodeCount();

      // build boundary sides informations and get the boundary sides map (sideParentElement, sideOrdinal)
      connManager->buildLocalBoundarySides();
      const auto boundarySides = connManager->getBoundarySides(sideFlag);
      int numBoundarySides = boundarySides.size();

      constexpr int maxNumNodesPerSide = 4;
      DynRankViewGId sideNodesGID("sideNodesGID", numBoundarySides, maxNumNodesPerSide);
      Kokkos::DynRankView<int,DeviceSpaceType> sideOrdinals("sideOrdinals", numBoundarySides);
      Kokkos::DynRankView<int,DeviceSpaceType> sideParentElem("sideParentElem", numBoundarySides);
      Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpaceType> sideOrts("sideOrts", numBoundarySides);

      //compute Kokkos maps for side ordinals and side parent cells.
      //compute side orientations when needed
      {
        auto sideNodesGIDHost = Kokkos::create_mirror_view(sideNodesGID);
        auto sideOrdinalsHost = Kokkos::create_mirror_view(sideOrdinals);
        auto sideParentElemHost = Kokkos::create_mirror_view(sideParentElem);
        for(int side=0; side<numBoundarySides; ++side) {
          int elem = boundarySides[side].first;
          sideParentElemHost(side) = elem;
          const auto GIDs = connManager->getConnectivity(elem);
          int elemSide = boundarySides[side].second;
          sideOrdinalsHost(side) = elemSide;
          for(size_t sideNode=0; sideNode < cellTopoPtr->getNodeCount(sideDim,elemSide); ++sideNode) {
            int node = cellTopoPtr->getNodeMap(sideDim, elemSide, sideNode);
            sideNodesGIDHost(side,sideNode) = GIDs[node];
          }
        }
        Kokkos::deep_copy(sideOrdinals, sideOrdinalsHost);
        Kokkos::deep_copy(sideParentElem, sideParentElemHost);

        // note, side orientations are needed when the parent basis needs orientations
        if(basis->requireOrientation()) {
          Kokkos::deep_copy(sideNodesGID, sideNodesGIDHost);
          ots::getOrientation(sideOrts, sideNodesGID, sideBasisPtr->getBaseCellTopology(), true);
        }
      }

      // compute side projections
      auto sideBasisCardinality = sideBasisPtr->getCardinality();
      DynRankView sideBasisCoeffsL2Proj("sideBasisCoeffsL2Proj", numBoundarySides, sideBasisCardinality);
      {
        ProjectionStruct sideProjStruct;
        sideProjStruct.createL2ProjectionStruct(sideBasisPtr.get(), targetCubDegree);

        auto sideEvaluationPoints = sideProjStruct.getAllEvalPoints();
        int numSidePoints = sideEvaluationPoints.extent(0);

        DynRankView physTargetAtSideEvalPoints,refTargetAtSideEvalPoints;
        if(sideFunctionSpace == Intrepid2::FUNCTION_SPACE_HCURL) {
          physTargetAtSideEvalPoints = DynRankView("physTargetAtSideEvalPoints", numBoundarySides, numSidePoints, dim);
          refTargetAtSideEvalPoints = DynRankView("refTargetAtSideEvalPoints", numBoundarySides, numSidePoints, sideDim);
        }
        else {
          physTargetAtSideEvalPoints = DynRankView("physTargetAtSideEvalPoints", numBoundarySides, numSidePoints);
          refTargetAtSideEvalPoints = DynRankView("refTargetAtSideEvalPoints", numBoundarySides, numSidePoints);
        }

        // maps reference points into physical space
        DynRankView physSideEvalPoints("physSideEvalPoints", numBoundarySides, numSidePoints, dim);
        DynRankView physVerticesOnSide("physVerticesOnSide", numBoundarySides, numNodesPerElem, dim);
        DynRankView sideEvaluationPoints3d("sideEvaluationPoints3d", numBoundarySides, numSidePoints, dim);
        {
          const auto subcellParametrization = Intrepid2::RefSubcellParametrization<DeviceSpaceType>::get(sideDim, cellTopoPtr->getKey());
          ct::mapToReferenceSubcell(sideEvaluationPoints3d, sideEvaluationPoints, subcellParametrization, sideOrdinals);

          Kokkos::DynRankView<int,DeviceSpaceType> sideNodeMap("sideNodeMap", cellTopoPtr->getSideCount(), maxNumNodesPerSide);
          auto sideNodeMapHost = Kokkos::create_mirror_view(sideNodeMap);
          for (size_t is=0; is < sideNodeMapHost.extent(0); is++) {
            for (size_t node=0; node < cellTopoPtr->getNodeCount(sideDim,is); node++)
              sideNodeMapHost(is, node) = cellTopoPtr->getNodeMap(sideDim, is, node);
          }
          Kokkos::deep_copy(sideNodeMap,sideNodeMapHost);

          DynRankView sidePhysVertices("sidePhysVertices", numBoundarySides, numNodesPerSide, dim);
          Kokkos::MDRangePolicy<typename DeviceSpaceType::execution_space,Kokkos::Rank<2>> rangePolicy({0,0},{numBoundarySides, dim});
          Kokkos::parallel_for(rangePolicy, KOKKOS_LAMBDA (const int &is, const int &d) {
            auto elem = sideParentElem(is);
            auto elemSide = sideOrdinals(is);

            for(int node=0; node < numNodesPerElem; ++node)
              physVerticesOnSide(is, node, d) = physVertexes(elem, node, d);

            for(int sideNode=0; sideNode < numNodesPerSide; ++sideNode) {
              int node = sideNodeMap(elemSide, sideNode);
              sidePhysVertices(is, sideNode, d) = physVertexes(elem, node, d);
            }
          });
          Kokkos::fence();

          ct::mapToPhysicalFrame(physSideEvalPoints,sideEvaluationPoints,sidePhysVertices,sideBasisPtr->getBaseCellTopology());
        }

        //evaluate the boundary term in the physical space and map it back to the reference element
        {
          DynRankView jacobian("jacobian",numBoundarySides, numSidePoints, dim, dim);
          DynRankView tangents("tangents", numBoundarySides, numSidePoints, dim, dim-1);
          DynRankView metricTensor("metricTensor", numBoundarySides, numSidePoints, dim-1, dim-1);
          DynRankView metricTensor_inv("metricTensor_inv", numBoundarySides, numSidePoints, dim-1, dim-1);
          DynRankView metricTensor_det("metricTensor_det", numBoundarySides, numSidePoints);

          ct::setJacobian(jacobian, sideEvaluationPoints3d, physVerticesOnSide, *cellTopoPtr);

          //compute  metric information
          if ((sideFunctionSpace == Intrepid2::FUNCTION_SPACE_HCURL) || (sideFunctionSpace == Intrepid2::FUNCTION_SPACE_HVOL)) {
            auto tangents0 = Kokkos::subview(tangents, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), 0);
            if(dim == 3) {
              auto tangents1 = Kokkos::subview(tangents, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), 1);
              ct::getPhysicalFaceTangents(tangents0, tangents1,jacobian,sideOrdinals,*cellTopoPtr);
            } else { //dim == 2
              ct::getPhysicalEdgeTangents(tangents0, jacobian,sideOrdinals,*cellTopoPtr);
            }
            rst::AtA(metricTensor,tangents);
            rst::det(metricTensor_det,metricTensor);  //note, this is equivalent to |normal|^2 = |tangents0 \times tangents1|^2
            if (sideFunctionSpace == Intrepid2::FUNCTION_SPACE_HCURL)
              rst::inverse(metricTensor_inv,metricTensor);
          }

          //evaluate the boundary terms in physical space
          double n_x(sideNormal[0]), n_y(sideNormal[1]), n_z(sideNormal[2]);
          Kokkos::parallel_for(Kokkos::RangePolicy<typename DeviceSpaceType::execution_space>(0,numBoundarySides),
              KOKKOS_LAMBDA (const int &is) {
            Fun fun(functionSpace);
            for(int pt=0;pt<numSidePoints;pt++) {

              auto x = physSideEvalPoints(is,pt,0), y = physSideEvalPoints(is,pt,1);
              auto z = (dim==2) ? 0.0 : physSideEvalPoints(is,pt,2);

              switch (functionSpace) {
              case Intrepid2::FUNCTION_SPACE_HGRAD:
                //compute f
                physTargetAtSideEvalPoints(is,pt) = fun(x,y,z,0);
                break;
              case Intrepid2::FUNCTION_SPACE_HDIV: {
                //compute f \dot n
                double fx(fun(x,y,z,0)), fy(fun(x,y,z,1)), fz(fun(x,y,z,2));
                physTargetAtSideEvalPoints(is,pt) = (dim == 3) ? n_x*fx+n_y*fy+n_z*fz : n_x*fx+n_y*fy;
              }
              break;
              case Intrepid2::FUNCTION_SPACE_HCURL:
              {
                //compute f \times n
                double fx(fun(x,y,z,0)), fy(fun(x,y,z,1)), fz(fun(x,y,z,2));
                if (dim == 3) {
                  physTargetAtSideEvalPoints(is,pt,0) = fy*n_z-fz*n_y;
                  physTargetAtSideEvalPoints(is,pt,1) = fz*n_x-fx*n_z;
                  physTargetAtSideEvalPoints(is,pt,2) = fx*n_y-fy*n_x;
                } else {
                  physTargetAtSideEvalPoints(is,pt) = fx*n_y-fy*n_x;
                }
              }
              break;
              default: {}
              }
            }
          });

          //map it back to reference space
          switch (functionSpace) {
          case Intrepid2::FUNCTION_SPACE_HGRAD:
            fst::mapHGradDataFromPhysSideToRefSide(refTargetAtSideEvalPoints,physTargetAtSideEvalPoints);
            break;
          case Intrepid2::FUNCTION_SPACE_HCURL:
            if(dim == 3)
              fst::mapHCurlDataCrossNormalFromPhysSideToRefSide(refTargetAtSideEvalPoints,tangents,metricTensor_inv,metricTensor_det,physTargetAtSideEvalPoints);
            else
              fst::mapHCurlDataCrossNormalFromPhysSideToRefSide(refTargetAtSideEvalPoints,metricTensor_det,physTargetAtSideEvalPoints);
            break;
          case Intrepid2::FUNCTION_SPACE_HDIV:
            fst::mapHDivDataDotNormalFromPhysSideToRefSide(refTargetAtSideEvalPoints,metricTensor_det,physTargetAtSideEvalPoints);
            break;
          default: {}
          }

          //perform the projections and get the basis coeffiecients
          pts::getL2BasisCoeffs(sideBasisCoeffsL2Proj,
              refTargetAtSideEvalPoints,
              sideOrts,
              sideBasisPtr.get(),
              &sideProjStruct);
        }

        // *************** CODE VERIFICATION FOR SIDE PROJECTIONS *********************
        {
          auto basisCoeffsL2ProjHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2Proj);
          auto sideBasisCoeffsL2ProjHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sideBasisCoeffsL2Proj);
          for(int side=0; side<numBoundarySides; ++side) {
            int elem = boundarySides[side].first;
            int elemSide = boundarySides[side].second;
            for(size_t sideNode=0; sideNode < cellTopoPtr->getNodeCount(sideDim,elemSide); ++sideNode) {
              int node = cellTopoPtr->getNodeMap(sideDim, elemSide, sideNode);
              if(sideBasisPtr->getDofCount(0,sideNode)!=0){
                int elemDof = basis->getDofOrdinal(0,node,0);
                int sideDof = sideBasisPtr->getDofOrdinal(0,sideNode,0);
                if(std::abs(basisCoeffsL2ProjHost(elem,elemDof)-sideBasisCoeffsL2ProjHost(side,sideDof))>1e-10) {
                  std::cout << "ERROR: Vertex DoFs don't match \n";
                  std::cout << "Point: " << node << ": " << basisCoeffsL2ProjHost(elem,elemDof) << " " << sideBasisCoeffsL2ProjHost(side,sideDof) << " " << basisCoeffsL2ProjHost(elem,elemDof) - sideBasisCoeffsL2ProjHost(side,sideDof) <<std::endl;
                  errorFlag++;
                }
              }
            }

            if(dim == 3) {
              for(size_t sideEdge=0; sideEdge < cellTopoPtr->getEdgeCount(2,elemSide); ++sideEdge) {
                int edge = Intrepid2::Orientation::getEdgeOrdinalOfFace(sideEdge,elemSide,*cellTopoPtr);
                for(int l=0; l<sideBasisPtr->getDofCount(1,sideEdge); l++) {
                  int elemDof = basis->getDofOrdinal(1,edge,l);
                  int sideDof = sideBasisPtr->getDofOrdinal(1,sideEdge,l);
                  if(std::abs(basisCoeffsL2ProjHost(elem,elemDof)-sideBasisCoeffsL2ProjHost(side,sideDof))>1e-10) {
                    std::cout << "ERROR: Edges DoFs don't match \n";
                    std::cout << "Edge: " << edge << ": " << basisCoeffsL2ProjHost(elem,elemDof) << " " << sideBasisCoeffsL2ProjHost(side,sideDof) << " " << basisCoeffsL2ProjHost(elem,elemDof) - sideBasisCoeffsL2ProjHost(side,sideDof) <<std::endl;
                    errorFlag++;
                  }
                }
              }
            }

            for(int l=0; l<sideBasisPtr->getDofCount(sideDim,0); l++) {
              int elemDof = basis->getDofOrdinal(sideDim,elemSide,l);
              int sideDof = sideBasisPtr->getDofOrdinal(sideDim,0,l);
              double diff = std::abs(basisCoeffsL2ProjHost(elem,elemDof)-sideBasisCoeffsL2ProjHost(side,sideDof));
              if(diff>1e-10) {
                std::cout << "ERROR: Side DoFs don't match \n";
                std::cout << "Side: " << elemSide << ": " << basisCoeffsL2ProjHost(elem,elemDof) << " " << sideBasisCoeffsL2ProjHost(side,sideDof) << ", " << diff <<std::endl;
                errorFlag++;
              }
            }
          }
        }
      }
    }

  } catch (const std::exception & err) {
    *outStream << " Exception\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }
 
 
  Teuchos::RCP<Teuchos::ParameterList> reportParams = parameterList(* (Teuchos::TimeMonitor::getValidReportParameters()));
  reportParams->set("Report format", "YAML");
  reportParams->set("YAML style", "spacious");
  if ( timingsFile != "" ){
    std::ofstream fout(timingsFile.c_str());
    Teuchos::TimeMonitor::report(Teuchos::rcpFromRef(comm).ptr(), fout, reportParams);
  } else {
    Teuchos::TimeMonitor::report(Teuchos::rcpFromRef(comm).ptr(), *outStream);
  }

  if (errorFlag != 0)
    *outStream << "End Result: TEST FAILED = " << errorFlag << "\n";
  else
    *outStream << "End Result: TEST PASSED\n";

  return errorFlag;
}
}
}
