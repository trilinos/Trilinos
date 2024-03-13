// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
  \brief  Assembly of high-order finite element matrix associated to Poisson problem -\nabla u + u = f with homogeneous Neumann BCs

  The mesh topology and global numbering of the Degrees of Freedom (DoFs) are computed using Panzer Dof Manager package
  Local matrices and right-hand-side are computed using high-order HGRAD elements on Hexahedra provided by the package Intrepid2
  Local matrices and right-hand-side are assembled into global matrix A and right-hand-side b using the
  FECrsMatrix and FEMultivector provided by Tpetra package

  The code is then verified by checking the norm of the residual r = b - A x.
  When the solution belong to the chosen finite element space, the residual norm should be close to machine eps.

  Command line arguments:
  --nx, --ny, --nz:  the number of owned 1d elements in the axis directions. The total number of owned hexas is given by nx \times ny \times \nz
  --px, --py, --pz:  the number of processes in the axis directions. The total numebr of processes is px \times py \times pz.
  --basis-degree: the degree of the HGRAD finite element basis; default: 4
  --quadrature-degree: the quadrature (cubature) degree, i.e. the maximum degree of polynomial that is integrated exactl; default: 2*basis_degree
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
#include "Intrepid2_LagrangianInterpolation.hpp"

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
    operator()(const double& x, const double& y, const double& z) const {
      return x*x*(x-1)*(x-1) + y*y*(y-1)*(y-1) + z*z*(z-1)*(z-1);
      //const double pi = 3.14159265358979323846;
      //return cos(pi*x)*cos(pi*y)*cos(pi*z);
    }
  };

  struct FunLapl {
    double
    KOKKOS_INLINE_FUNCTION
    operator()(const double& x, const double& y, const double& z) const {
      return  2*(x-1)*(x-1)+8*x*(x-1)+2*x*x + 2*(y-1)*(y-1)+8*y*(y-1)+2*y*y +2*(z-1)*(z-1)+8*z*(z-1)+2*z*z;
      //const double pi = 3.14159265358979323846;
      //  return -3 * pi * pi *cos(pi*x)*cos(pi*y)*cos(pi*z);
    }
  };

  struct FunRhs {
    double
    KOKKOS_INLINE_FUNCTION
    operator()(const double& x, const double& y, const double& z) const {
      return  fun(x,y,z) - funLapl(x,y,z);
    }
    Fun fun;
    FunLapl funLapl;
  };

  template<typename DynRankViewType>
  class EvalRhsFunctor {
  public:
    DynRankViewType funAtPoints;
    DynRankViewType points;
    Fun fun;
    FunLapl funLapl;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int elem) const
    {
      for(int i=0;i<static_cast<int>(points.extent(1));i++) {
        auto x = points(elem,i,0), y = points(elem,i,1), z= points(elem,i,2);
        funAtPoints(elem,i) =fun(x,y,z) - funLapl(x,y,z);
      }
    }
  };

  template<typename DynRankViewType>
  class EvalSolFunctor {
  public:
    DynRankViewType funAtPoints;
    DynRankViewType points;
    Fun fun;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int elem) const
    {
      for(int i=0;i<static_cast<int>(points.extent(1));i++) {
        auto x = points(elem,i,0), y = points(elem,i,1), z= points(elem,i,2);
        funAtPoints(elem,i) = fun(x,y,z);
      }
    }
  };

template<typename ValueType, typename DeviceSpaceType>
int feAssemblyHex(int argc, char *argv[]) {

  // host_memory/execution/mirror_space deprecated for kokkos@3.7.00, removed after release
  // see https://github.com/kokkos/kokkos/pull/3973
  using exec_space = typename Kokkos::is_space<DeviceSpaceType>::execution_space;
  using mem_space = typename Kokkos::is_space<DeviceSpaceType>::memory_space;
  using do_not_use_host_memory_space = std::conditional_t<
      std::is_same<mem_space, Kokkos::HostSpace>::value
#if defined(KOKKOS_ENABLE_CUDA)
          || std::is_same<mem_space, Kokkos::CudaUVMSpace>::value ||
          std::is_same<mem_space, Kokkos::CudaHostPinnedSpace>::value
#elif defined(KOKKOS_ENABLE_HIP)
          || std::is_same<mem_space,
                          Kokkos::HIPHostPinnedSpace>::value ||
          std::is_same<mem_space,
                       Kokkos::HIPManagedSpace>::value
#elif defined(KOKKOS_ENABLE_SYCL)
          || std::is_same<mem_space,
                          Kokkos::Experimental::SYCLSharedUSMSpace>::value ||
          std::is_same<mem_space,
                       Kokkos::Experimental::SYCLHostUSMSpace>::value
#endif
      ,
      mem_space, Kokkos::HostSpace>;

  using do_not_use_host_execution_space = std::conditional_t<
#if defined(KOKKOS_ENABLE_CUDA)
      std::is_same<exec_space, Kokkos::Cuda>::value ||
#elif defined(KOKKOS_ENABLE_HIP)
      std::is_same<exec_space, Kokkos::HIP>::value ||
#elif defined(KOKKOS_ENABLE_SYCL)
      std::is_same<exec_space, Kokkos::Experimental::SYCL>::value ||
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
      std::is_same<exec_space,
                   Kokkos::Experimental::OpenMPTarget>::value ||
#endif
          false,
      Kokkos::DefaultHostExecutionSpace, exec_space>;

  using host_mirror_space = std::conditional_t<
      std::is_same<exec_space, do_not_use_host_execution_space>::value &&
          std::is_same<mem_space, do_not_use_host_memory_space>::value,
      DeviceSpaceType,
      Kokkos::Device<do_not_use_host_execution_space,
                     do_not_use_host_memory_space>>;

  using HostSpaceType = typename host_mirror_space::execution_space;

  using DynRankView = Kokkos::DynRankView<ValueType,DeviceSpaceType>;

  using map_t = Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal>;

  using local_ordinal_t = typename map_t::local_ordinal_type;
  using global_ordinal_t = typename map_t::global_ordinal_type;
  using node_t = typename map_t::node_type;
  using fe_graph_t = Tpetra::FECrsGraph<local_ordinal_t,global_ordinal_t,node_t>;
  using scalar_t = ValueType;
  using fe_matrix_t = Tpetra::FECrsMatrix<scalar_t, local_ordinal_t, global_ordinal_t>;
  using fe_multivector_t = Tpetra::FEMultiVector<scalar_t, local_ordinal_t, global_ordinal_t>;
  using vector_t = Tpetra::Vector<scalar_t, local_ordinal_t, global_ordinal_t>;

  using DynRankViewGId = Kokkos::DynRankView<global_ordinal_t,DeviceSpaceType>;

  using ct = Intrepid2::CellTools<DeviceSpaceType>;
  using ots = Intrepid2::OrientationTools<DeviceSpaceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceSpaceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceSpaceType>;
  using li = Intrepid2::LagrangianInterpolation<DeviceSpaceType>;

  int errorFlag = 0;


#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);

  //output stream/file
  Teuchos::RCP<Teuchos::FancyOStream> outStream;
  std::string timingsFile = "";

  try {


    // ************************************ GET INPUTS **************************************
    constexpr local_ordinal_t dim = 3;
    int degree = 4;
    local_ordinal_t nx = 2;
    local_ordinal_t ny            = nx;
    local_ordinal_t nz            = nx;
    int np   = comm.getSize(); // number of processors
    int px = std::cbrt(np); while(np%px!=0) --px;
    int py = std::sqrt(np/px); while(np%py!=0) --py;
    int pz = np/(px*py);
    constexpr int bx=1, by=1, bz=1;  //blocks on each process. Here we assume there is only one block per process
    int verbose = 1;

    // parse command line arguments
    Teuchos::CommandLineProcessor clp;
    clp.setOption("nx",&nx);
    clp.setOption("ny",&ny);
    clp.setOption("nz",&nz);
    clp.setOption("px",&px);
    clp.setOption("py",&py);
    clp.setOption("pz",&pz);
    clp.setOption("basis-degree",&degree);
    int cubDegree = -1;
    clp.setOption("quadrature-degree",&cubDegree);
    clp.setOption("timings-file",&timingsFile);
    clp.setOption("verbose", &verbose);
    auto cmdResult = clp.parse(argc,argv);
    cubDegree = (cubDegree == -1) ? 2*degree : cubDegree;

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

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >() );

    // Get dimensions
    int numNodesPerElem = hexa.getNodeCount();

    // build the topology
    auto connManager = Teuchos::rcp(new panzer::unit_test::CartesianConnManager);
    connManager->initialize(comm,
        global_ordinal_t(nx*px),
        global_ordinal_t(ny*py),
        global_ordinal_t(nz*pz),
        px,py,pz,bx,by,bz);

    // *********************************** COMPUTE GLOBAL IDs OF VERTICES AND DOFs  ************************************

    // build the dof manager, and assocaite with the topology
    auto dofManager = Teuchos::rcp(new panzer::DOFManager);
    dofManager->setConnManager(connManager,*comm.getRawMpiComm());

    // add solution field to the element block
    Teuchos::RCP< Intrepid2::Basis<DeviceSpaceType, scalar_t,scalar_t> > basis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DeviceSpaceType,scalar_t,scalar_t>(degree));
    auto basisCardinality = basis->getCardinality();
    Teuchos::RCP<panzer::Intrepid2FieldPattern> fePattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
    dofManager->addField("block-0_0_0",fePattern);

    // try to get them all synced up
    comm.barrier();

    dofManager->buildGlobalUnknowns();

    local_ordinal_t numOwnedElems = nx*ny*nz;

    *outStream << "Number of elements on each processor nx X ny X nz: " << numOwnedElems << "\n";
    *outStream << "    nx" << "   ny" << "   nz\n";
    *outStream << std::setw(5) << nx <<
        std::setw(5) << ny <<
        std::setw(5) << nz << "\n\n";
    *outStream << "Number of processors px X py X pz: " << px*py*pz << "\n";
    *outStream << "    px" << "   py" << "   pz\n";
    *outStream << std::setw(5) << px <<
        std::setw(5) << py <<
        std::setw(5) << pz << "\n\n";

    global_ordinal_t totalNumElements = numOwnedElems*px*py*pz;
     *outStream << "Total number of elements: " << totalNumElements << ", number of DoFs per element: " << basisCardinality << "\n";

    // Print mesh information

    // Cube
    scalar_t leftX = 0.0, rightX = 1.0;
    scalar_t leftY = 0.0, rightY = 1.0;
    scalar_t leftZ = 0.0, rightZ = 1.0;

    // Mesh spacing
    scalar_t hx = (rightX-leftX)/((scalar_t)(nx*px*bx));
    scalar_t hy = (rightY-leftY)/((scalar_t)(ny*py*by));
    scalar_t hz = (rightZ-leftZ)/((scalar_t)(nz*pz*bz));


    // *********************************** COMPUTE COORDINATES OF PHYSICAL VERTICES  ************************************

    // Get coordinates of physical vertices
    DynRankView ConstructWithLabel(physVertexes, numOwnedElems, numNodesPerElem, dim);
    {
      auto physVertexesHost = Kokkos::create_mirror_view(physVertexes);
      DynRankView ConstructWithLabel(refVertices, numNodesPerElem, dim);
      Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,scalar_t,scalar_t> hexaLinearBasis;
      hexaLinearBasis.getDofCoords(refVertices);
      auto refVerticesHost = Kokkos::create_mirror_view(refVertices);
      Kokkos::deep_copy(refVerticesHost, refVertices);

      auto elemTriplet = connManager->getMyBrickElementsTriplet();
      double h[3] = {hx, hy, hz};

      for(int i=0; i<numOwnedElems; ++i) {
        elemTriplet =  connManager->computeLocalBrickElementGlobalTriplet(i,connManager->getMyBrickElementsTriplet(),connManager->getMyBrickOffsetTriplet());
        double offset[3] = {leftX + elemTriplet.x*hx+hx/2, leftY +elemTriplet.y*hy+hy/2, leftX +elemTriplet.z*hz+hz/2};
        for(int j=0; j<numNodesPerElem; ++j) {
          for(int k=0; k<dim; ++k)
            physVertexesHost(i,j,k) = offset[k]+h[k]/2.0*refVerticesHost(j,k);
        }
      }
      Kokkos::deep_copy(physVertexes, physVertexesHost);
    }

    meshTimer = Teuchos::null;


    // *********************************** COMPUTE ELEMENTS' ORIENTATION BASED ON GLOBAL IDs  ************************************


    auto localFeAssemblyTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Local Finite Element Assembly")));

    //compute global ids of element vertices
    DynRankViewGId ConstructWithLabel(elemNodesGID, numOwnedElems, numNodesPerElem);
    {
      auto elemNodesGID_host = Kokkos::create_mirror_view(elemNodesGID);

      for(int i=0; i<numOwnedElems; ++i) {
        const auto GIDs = connManager->getConnectivity(i);
        for(int j=0; j<numNodesPerElem; ++j) {
          elemNodesGID_host(i,j) = GIDs[j];
        }
      }
      Kokkos::deep_copy(elemNodesGID,elemNodesGID_host);
    }

    // compute orientations for cells (one time computation)
    Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpaceType> elemOrts("elemOrts", numOwnedElems);
    ots::getOrientation(elemOrts, elemNodesGID, hexa);


    // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************
    // Compute quadrature (cubature) points
    Intrepid2::DefaultCubatureFactory cubFactory;
    auto cellCub = cubFactory.create<DeviceSpaceType, scalar_t, scalar_t>(hexa.getBaseKey(), cubDegree);
    auto numQPoints = cellCub->getNumPoints();
    DynRankView ConstructWithLabel(quadPoints, numQPoints, dim);
    DynRankView ConstructWithLabel(weights, numQPoints);
    cellCub->getCubature(quadPoints, weights);


    //Compute physical Dof Coordinates and Reference coordinates
    DynRankView ConstructWithLabel(funAtQPoints, numOwnedElems, numQPoints);
    {
      DynRankView ConstructWithLabel(physQPoints, numOwnedElems, numQPoints, dim);
      ct::mapToPhysicalFrame(physQPoints,quadPoints,physVertexes,basis->getBaseCellTopology());
      EvalRhsFunctor<DynRankView> functor;
      functor.funAtPoints = funAtQPoints;
      functor.points = physQPoints;
      Kokkos::parallel_for("loop for evaluating the rhs at quadrature points", numOwnedElems,functor);
    }

    // compute oriented basis functions at quadrature points
    DynRankView ConstructWithLabel(basisValuesAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints);
    DynRankView ConstructWithLabel(transformedBasisValuesAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints);
    DynRankView basisValuesAtQPointsCells("inValues", numOwnedElems, basisCardinality, numQPoints);
    DynRankView ConstructWithLabel(basisValuesAtQPoints, basisCardinality, numQPoints);
    basis->getValues(basisValuesAtQPoints, quadPoints);
    rst::clone(basisValuesAtQPointsCells,basisValuesAtQPoints);

    // modify basis values to account for orientations
    ots::modifyBasisByOrientation(basisValuesAtQPointsOriented,
        basisValuesAtQPointsCells,
        elemOrts,
        basis.getRawPtr());

    // transform basis values
    fst::HGRADtransformVALUE(transformedBasisValuesAtQPointsOriented, basisValuesAtQPointsOriented);

    DynRankView ConstructWithLabel(basisGradsAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints, dim);
    DynRankView ConstructWithLabel(transformedBasisGradsAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints, dim);
    DynRankView basisGradsAtQPointsCells("inValues", numOwnedElems, basisCardinality, numQPoints, dim);
    DynRankView ConstructWithLabel(basisGradsAtQPoints, basisCardinality, numQPoints, dim);
    basis->getValues(basisGradsAtQPoints, quadPoints, Intrepid2::OPERATOR_GRAD);
    rst::clone(basisGradsAtQPointsCells,basisGradsAtQPoints);

    // modify basis values to account for orientations
    ots::modifyBasisByOrientation(basisGradsAtQPointsOriented,
        basisGradsAtQPointsCells,
        elemOrts,
        basis.getRawPtr());

    // map basis functions to reference (oriented) element
    DynRankView ConstructWithLabel(jacobianAtQPoints, numOwnedElems, numQPoints, dim, dim);
    DynRankView ConstructWithLabel(jacobianAtQPoints_inv, numOwnedElems, numQPoints, dim, dim);
    DynRankView ConstructWithLabel(jacobianAtQPoints_det, numOwnedElems, numQPoints);
    ct::setJacobian(jacobianAtQPoints, quadPoints, physVertexes, hexa);
    ct::setJacobianInv (jacobianAtQPoints_inv, jacobianAtQPoints);

    fst::HGRADtransformGRAD(transformedBasisGradsAtQPointsOriented, jacobianAtQPoints_inv, basisGradsAtQPointsOriented);

    // compute integrals to assembly local matrices
    DynRankView elemsMat("elemsMat", numOwnedElems, basisCardinality, basisCardinality),
        elemsRHS("elemsRHS", numOwnedElems, basisCardinality);

    DynRankView ConstructWithLabel(weightedTransformedBasisValuesAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints);
    DynRankView ConstructWithLabel(weightedTransformedBasisGradsAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints, dim);
    DynRankView ConstructWithLabel(cellWeights, numOwnedElems, numQPoints);
    rst::clone(cellWeights, weights);

    fst::multiplyMeasure(weightedTransformedBasisGradsAtQPointsOriented, cellWeights, transformedBasisGradsAtQPointsOriented);
    fst::multiplyMeasure(weightedTransformedBasisValuesAtQPointsOriented, cellWeights, transformedBasisValuesAtQPointsOriented);

    fst::integrate(elemsMat, transformedBasisGradsAtQPointsOriented, weightedTransformedBasisGradsAtQPointsOriented);
    fst::integrate(elemsMat, transformedBasisValuesAtQPointsOriented, weightedTransformedBasisValuesAtQPointsOriented, true);
    Kokkos::fence(); //make sure that funAtQPoints has been evaluated
    fst::integrate(elemsRHS, funAtQPoints, weightedTransformedBasisValuesAtQPointsOriented);

    localFeAssemblyTimer =  Teuchos::null;


    // ************************************ GENERATE GRAPH **************************************

    auto mapsTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Creating Maps")));

    auto globalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(dofManager);
    std::vector<global_ordinal_t> ownedIndices, ownedAndGhostedIndices;
    globalIndexer->getOwnedIndices(ownedIndices);
    Teuchos::RCP<const map_t> ownedMap = Teuchos::rcp(new map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedIndices,0,Teuchos::rcpFromRef(comm)));
    globalIndexer->getOwnedAndGhostedIndices(ownedAndGhostedIndices);
    Teuchos::RCP<const map_t> ownedAndGhostedMap = Teuchos::rcp(new const map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedAndGhostedIndices,0,Teuchos::rcpFromRef(comm)));

     *outStream << "Total number of DoFs: " << ownedMap->getGlobalNumElements() << ", number of owned DoFs: " << ownedMap->getLocalNumElements() << "\n";

    mapsTimer = Teuchos::null;
    auto graphGenerationTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Graph Generation")));

    //to compute the max number of nonzero in a row, we consider a patch of 8 hexas sharing a vertex
    int numVertices(27), numEdges(54), numFaces(36), numCells(8);
    auto numDofsPerVertex = basis->getDofCount(0,0);
    auto numDofsPerEdge = basis->getDofCount(1,0);
    auto numDofsPerFace = basis->getDofCount(2,0);
    auto numDofsPerCell = basis->getDofCount(3,0);
    auto maxNumRowEntries = numVertices*numDofsPerVertex+numEdges*numDofsPerEdge+
        numFaces*numDofsPerFace + numCells*numDofsPerCell;

    // this constructor ensures that the local ids in the owned+ghosted map and in the graph col map corresponds to the same global ids
    // in our case the owned row map is the same as the (owned) domain map
    auto feGraph = Teuchos::rcp(new fe_graph_t(ownedMap, ownedAndGhostedMap, maxNumRowEntries, ownedAndGhostedMap, Teuchos::null, ownedMap));

    Teuchos::Array<global_ordinal_t> globalIdsInRow(basisCardinality);
    const std::string blockId = "eblock-0_0_0";
    auto elmtOffsetKokkos = dofManager->getGIDFieldOffsetsKokkos(blockId,0);
    auto elmtOffsets_host = dofManager->getGIDFieldOffsets(blockId,0);

    // fill graph
    // for each element in the mesh...
    Tpetra::beginAssembly(*feGraph);
    for(int elemId=0; elemId<numOwnedElems; elemId++)
    {
      // Populate globalIdsInRow:
      // - Copy the global node ids for current element into an array.
      // - Since each element's contribution is a clique, we can re-use this for
      //   each row associated with this element's contribution.
      std::vector<global_ordinal_t> elementGIDs;
      dofManager->getElementGIDs(elemId, elementGIDs);
      for(int nodeId=0; nodeId<basisCardinality; nodeId++) {
        globalIdsInRow[nodeId] = elementGIDs[elmtOffsets_host[nodeId]];
      }

      // Add the contributions from the current row into the graph.
      // - For example, if Element 0 contains nodes [0,1,4,5, ...] then we insert the nodes:
      //   - node 0 inserts [0, 1, 4, 5, ...]
      //   - node 1 inserts [0, 1, 4, 5, ...]
      //   - node 4 inserts [0, 1, 4, 5, ...]
      //   - node 5 inserts [0, 1, 4, 5, ...]
      for(int nodeId=0; nodeId<basisCardinality; nodeId++)
      {
        feGraph->insertGlobalIndices(globalIdsInRow[nodeId], globalIdsInRow());
      }
    }
    Tpetra::endAssembly(*feGraph);

    graphGenerationTimer = Teuchos::null;

    // ************************************ MATRIX ASSEMBLY **************************************

    auto matrixAndRhsAllocationTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Allocation of Matrix and Rhs")));

    auto A = Teuchos::rcp(new fe_matrix_t(feGraph));
    auto b = Teuchos::rcp (new fe_multivector_t(ownedMap, feGraph->getImporter(), 1));

    matrixAndRhsAllocationTimer =  Teuchos::null;

    auto matrixAndRhsFillTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Fill of Matrix and Rhs")));
    Teuchos::Array<local_ordinal_t> columnLocalIds(basisCardinality);
    Teuchos::Array<global_ordinal_t> bLocalIds(basisCardinality);
    Teuchos::Array<scalar_t> columnScalarValues(basisCardinality);         // scalar values for each column

    auto localColMap  = A->getColMap()->getLocalMap();
    auto localMap  = ownedAndGhostedMap->getLocalMap();
    auto localMatrix  = A->getLocalMatrixDevice();
    auto localRHS     = b->getLocalViewDevice(Tpetra::Access::ReadWrite);

    //fill matrix
    // Loop over elements
    Tpetra::beginAssembly(*A, *b);

    std::vector<global_ordinal_t> elementGIDs(basisCardinality);
    auto elementLIDs = globalIndexer->getLIDs();

    Kokkos::parallel_for
      ("Assemble FE matrix and right-hand side",
       Kokkos::RangePolicy<DeviceSpaceType, int> (0, numOwnedElems),
       KOKKOS_LAMBDA (const size_t elemId) {
        // Get subviews
        auto elemRHS    = Kokkos::subview(elemsRHS,elemId, Kokkos::ALL());
        auto elemMat = Kokkos::subview(elemsMat,elemId, Kokkos::ALL(), Kokkos::ALL());
        auto elemLIds  = Kokkos::subview(elementLIDs,elemId, Kokkos::ALL());

        // For each node (row) on the current element
        for (local_ordinal_t nodeId = 0; nodeId < basisCardinality; ++nodeId) {
          const local_ordinal_t localRowId = elemLIds(elmtOffsetKokkos(nodeId));

          // Force atomics on sums
          for (local_ordinal_t colId = 0; colId < basisCardinality; ++colId)
            localMatrix.sumIntoValues (localRowId, &elemLIds(elmtOffsetKokkos(colId)), 1, &(elemMat(nodeId,colId)), true, true);

          Kokkos::atomic_add (&(localRHS(localRowId,0)), elemRHS(nodeId));
        }
      });

    // Release the device view
    localRHS = decltype(localRHS)("empty",0,0);

    Tpetra::endAssembly(*A, *b);

    matrixAndRhsFillTimer =  Teuchos::null;


    // ************************************ CODE VERIFICATION **************************************

    //interpolate analytic solution into finite element space
    DynRankView ConstructWithLabel(basisCoeffsLI, numOwnedElems, basisCardinality);
    {
      Teuchos::TimeMonitor liTimer =  *Teuchos::TimeMonitor::getNewTimer("Verification, locally interpolate analytic solution");
      DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);

      basis->getDofCoords(dofCoords);

      DynRankView ConstructWithLabel(funAtDofPoints, numOwnedElems, basisCardinality);
      {
        DynRankView ConstructWithLabel(physDofPoints, numOwnedElems, basisCardinality, dim);
        ct::mapToPhysicalFrame(physDofPoints,dofCoords,physVertexes,basis->getBaseCellTopology());
        EvalSolFunctor<DynRankView> functor;
        functor.funAtPoints = funAtDofPoints;
        functor.points = physDofPoints;
        Kokkos::parallel_for("loop for evaluating the function at DoF points", numOwnedElems,functor);
        Kokkos::fence(); //make sure that funAtDofPoints has been evaluated
      }

      li::getBasisCoeffs(basisCoeffsLI, funAtDofPoints, basis.getRawPtr(), elemOrts);
    }

    {
      Teuchos::TimeMonitor vTimer1 =  *Teuchos::TimeMonitor::getNewTimer("Verification, assemble solution");
      vector_t x(ownedMap); //solution
      auto basisCoeffsLIHost = Kokkos::create_mirror_view(basisCoeffsLI);
      Kokkos::deep_copy(basisCoeffsLIHost,basisCoeffsLI);
      for(int elemId=0; elemId<numOwnedElems; elemId++)
      {
        dofManager->getElementGIDs(elemId, elementGIDs);


        for(int nodeId=0; nodeId<basisCardinality; nodeId++)
        {
          global_ordinal_t gid = elementGIDs[elmtOffsets_host[nodeId]];
          if(ownedMap->isNodeGlobalElement(gid))
            x.replaceGlobalValue(gid, basisCoeffsLIHost(elemId, nodeId));
        }
      }

      {
        Teuchos::TimeMonitor vTimer2 =  *Teuchos::TimeMonitor::getNewTimer("Verification,compute rhs (matrix-vector product)");
        A->apply(x, *b, Teuchos::NO_TRANS, -1.0, 1.0);   // b - A x
      }
    }

    double res_l2_norm = b->getVector(0)->norm2();
    if((degree >= 4) && (res_l2_norm > 1e-10)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "Residual norm should be close to machine eps, but it is instead: " << res_l2_norm <<  "\n";
    }
    else {
      *outStream << "Residual l2 norm : " << res_l2_norm << "\n";
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
