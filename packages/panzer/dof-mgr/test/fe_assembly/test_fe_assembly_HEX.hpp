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
  \brief  Test for checking accuracy of interpolation-based projections for Hexaedral elements

  The test considers a uniform and structured hexahedral mesh of the cube [-1,1]^3, formed by
  N^3 hexas, and checks the accuracy of the HGRAD, HCURL, HDIV, HVOL projections of analytic
  target functions for increasing N.
  The accuracy is computed in the H^1, H^{curl}, H^{div} and L^2 norms respectively. The optimal
  order of convergence equates the basis degree.

  \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

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


#define Intrepid2_Experimental


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <array>
#include <set>
#include <random>
#include <algorithm>

namespace Intrepid2 {

namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }


  struct Fun {
    double
    KOKKOS_INLINE_FUNCTION
    operator()(const double& x, const double& y, const double& z) const {
      //const double pi = 3.14159265358979323846;
      return x*x*(x-1)*(x-1) + y*y*(y-1)*(y-1) + z*z*(z-1)*(z-1);
      //return cos(pi*x)*cos(pi*y)*cos(pi*z);
    }
  };

  struct FunLapl {
    double
    KOKKOS_INLINE_FUNCTION
    operator()(const double& x, const double& y, const double& z) const {
      //const double pi = 3.14159265358979323846;
      return  2*(x-1)*(x-1)+8*x*(x-1)+2*x*x + 2*(y-1)*(y-1)+8*y*(y-1)+2*y*y +2*(z-1)*(z-1)+8*z*(z-1)+2*z*z;
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

  typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
  typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
  typedef Kokkos::DynRankView<ValueType,HostSpaceType> DynRankViewHost;

  typedef Tpetra::Map<>::local_ordinal_type  local_ordinal_t;
  typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
  typedef Tpetra::Map<>::node_type           node_t;
  typedef Tpetra::Map<>             map_t;
  typedef Tpetra::FECrsGraph<local_ordinal_t,global_ordinal_t,node_t>      fe_graph_t;
  typedef ValueType scalar_t;
  typedef Tpetra::FECrsMatrix<scalar_t> fe_matrix_t;
  typedef Tpetra::FEMultiVector<scalar_t> fe_multivector_t;

  //typedef Kokkos::DynRankView<ordinal_type,DeviceSpaceType> DynRankViewInt;
  typedef Kokkos::DynRankView<global_ordinal_t,HostSpaceType> DynRankViewIntHost;


  int errorFlag = 0;


  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;
  typedef Experimental::LagrangianInterpolation<DeviceSpaceType> li;




#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  bool verbose = true;
  Teuchos::MpiComm<ordinal_type> comm(MPI_COMM_WORLD);
  Teuchos::RCP<Teuchos::FancyOStream> outStream = ((comm.getRank () == 0) && verbose) ?
      getFancyOStream(Teuchos::rcpFromRef (std::cout)) :
      getFancyOStream(Teuchos::rcp (new Teuchos::oblackholestream ()));

  *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
  *outStream << "\n";

  try {

    // ************************************ GET INPUTS **************************************


    constexpr ordinal_type dim = 3;
    ordinal_type degree = 4;
    local_ordinal_t nx = 2;
    local_ordinal_t ny            = nx;
    local_ordinal_t nz            = nx;
    ordinal_type np   = comm.getSize(); // number of processors
    ordinal_type px = std::cbrt(np), py = px, pz = px;
    constexpr ordinal_type bx=1, by=1, bz=1;  //blocks on each process. Here we assume there is only one block per process

    // timings output
    std::string timingsFile = "timings.yaml";
    // parse command line arguments
    Teuchos::CommandLineProcessor clp;
    clp.setOption("nx",&nx);
    clp.setOption("ny",&ny);
    clp.setOption("nz",&nz);
    clp.setOption("px",&px);
    clp.setOption("py",&py);
    clp.setOption("pz",&pz);
    clp.setOption("basis_degree",&degree);
    ordinal_type cub_degree = 2*degree;
    clp.setOption("cub_degree",&cub_degree);
    clp.setOption("timings-file",&timingsFile);
    auto cmdResult = clp.parse(argc,argv);
    if(cmdResult!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      clp.printHelpMessage(argv[0],std::cout);
      return -1;
    }

    // *********************************** MESH TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >() );

    // Get dimensions
    ordinal_type numNodesPerElem = hexa.getNodeCount();

    // build the topology
    auto connManager = Teuchos::rcp(new panzer::unit_test::CartesianConnManager);
    connManager->initialize(comm,
        global_ordinal_t(nx*px),
        global_ordinal_t(ny*py),
        global_ordinal_t(nz*pz),
        px,py,pz,bx,by,bz);


    // *********************************** COMPUTE GLOBAL IDs OF VERTICES AND DOFs  ************************************


    // build the dof manager, and assocaite with the topology
    auto dofManager = Teuchos::rcp(new panzer::DOFManager<local_ordinal_t, global_ordinal_t>);
    dofManager->setConnManager(connManager,*comm.getRawMpiComm());

    // add solution field to the element block
    Teuchos::RCP< Intrepid2::Basis<DeviceSpaceType, scalar_t,scalar_t> > basis = Teuchos::rcp(new Basis_HGRAD_HEX_Cn_FEM<DeviceSpaceType,scalar_t,scalar_t>(degree));
    Teuchos::RCP<panzer::Intrepid2FieldPattern> fe_pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
    dofManager->addField("block-0_0_0",fe_pattern);

    // try to get them all synced up
    comm.barrier();

    dofManager->buildGlobalUnknowns();


    //*outStream << "Generating mesh ... \n\n";

    //*outStream << "    nx" << "   ny" << "   nz\n";
    //*outStream << std::setw(5) << nx <<
    //    std::setw(5) << ny <<
    //    std::setw(5) << nz << "\n\n";

    // Print mesh information
    local_ordinal_t numOwnedElems = nx*ny*nz;

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
      Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,scalar_t,scalar_t> hexaLinearBasis;
      hexaLinearBasis.getDofCoords(refVertices);
      auto refVerticesHost = Kokkos::create_mirror_view(refVertices);   
      Kokkos::deep_copy(refVerticesHost, refVertices);
   
      auto elemTriplet = connManager->getMyElementsTriplet();
      double h[3] = {hx, hy, hz};

      for(ordinal_type i=0; i<numOwnedElems; ++i) {
        elemTriplet =  connManager->computeLocalElementGlobalTriplet(i,connManager->getMyElementsTriplet(),connManager->getMyOffsetTriplet());
        double offset[3] = {leftX + elemTriplet.x*hx+hx/2, leftY +elemTriplet.y*hy+hy/2, leftX +elemTriplet.z*hz+hz/2};
        for(ordinal_type j=0; j<numNodesPerElem; ++j) {
          for(ordinal_type k=0; k<dim; ++k)
            physVertexesHost(i,j,k) = offset[k]+h[k]/2.0*refVerticesHost(j,k);
        }
      }
      Kokkos::deep_copy(physVertexes, physVertexesHost);
    }


    // *********************************** COMPUTE ELEMENTS' ORIENTATION BASED ON GLOBAL IDs  ************************************

    //compute global ids of element vertices
    DynRankViewIntHost ConstructWithLabel(elemNodesGID, numOwnedElems, numNodesPerElem);
    {
      for(ordinal_type i=0; i<numOwnedElems; ++i) {
        const auto GIDs = connManager->getConnectivity(i);
        for(ordinal_type j=0; j<numNodesPerElem; ++j) {
          elemNodesGID(i,j) = GIDs[j];
        }
      }
    }

    // compute orientations for cells (one time computation)
    Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numOwnedElems);
    ots::getOrientation(elemOrts, elemNodesGID, hexa);


    // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************
    DefaultCubatureFactory cub_factory;
    auto cell_cub = cub_factory.create<DeviceSpaceType, scalar_t, scalar_t>(hexa.getBaseKey(), cub_degree);
    ordinal_type numQPoints = cell_cub->getNumPoints();
    DynRankView ConstructWithLabel(quadPoints, numQPoints, dim);
    DynRankView ConstructWithLabel(weights, numQPoints);
    cell_cub->getCubature(quadPoints, weights);
    ordinal_type basisCardinality = basis->getCardinality();



    //Compute physical Dof Coordinates and Reference coordinates
    DynRankView ConstructWithLabel(funAtQPoints, numOwnedElems, numQPoints);
    {
      DynRankView ConstructWithLabel(physQPoints, numOwnedElems, numQPoints, dim);
      ct::mapToPhysicalFrame(physQPoints,quadPoints,physVertexes,basis->getBaseCellTopology());
      EvalRhsFunctor<DynRankView> functor;
      functor.funAtPoints = funAtQPoints;
      functor.points = physQPoints;
      Kokkos::parallel_for(numOwnedElems,functor);
    }

    //check that fun values at reference points coincide with those computed using basis functions
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
    basis->getValues(basisGradsAtQPoints, quadPoints,OPERATOR_GRAD);
    rst::clone(basisGradsAtQPointsCells,basisGradsAtQPoints);

    // modify basis values to account for orientations
    ots::modifyBasisByOrientation(basisGradsAtQPointsOriented,
        basisGradsAtQPointsCells,
        elemOrts,
        basis.getRawPtr());

    DynRankView ConstructWithLabel(jacobianAtQPoints, numOwnedElems, numQPoints, dim, dim);
    DynRankView ConstructWithLabel(jacobianAtQPoints_inv, numOwnedElems, numQPoints, dim, dim);
    DynRankView ConstructWithLabel(jacobianAtQPoints_det, numOwnedElems, numQPoints);
    ct::setJacobian(jacobianAtQPoints, quadPoints, physVertexes, hexa);
    ct::setJacobianInv (jacobianAtQPoints_inv, jacobianAtQPoints);

    fst::HGRADtransformGRAD(transformedBasisGradsAtQPointsOriented, jacobianAtQPoints_inv, basisGradsAtQPointsOriented);

    DynRankView cellMat("cellMassMat", numOwnedElems, basisCardinality, basisCardinality),
        cellRhs("cellRhs", numOwnedElems, basisCardinality);

    DynRankView ConstructWithLabel(weightedTransformedBasisValuesAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints);
    DynRankView ConstructWithLabel(weightedTransformedBasisGradsAtQPointsOriented, numOwnedElems, basisCardinality, numQPoints, dim);
    DynRankView ConstructWithLabel(cellWeights, numOwnedElems, numQPoints);
    rst::clone(cellWeights, weights);

    fst::multiplyMeasure(weightedTransformedBasisGradsAtQPointsOriented, cellWeights, transformedBasisGradsAtQPointsOriented);
    fst::multiplyMeasure(weightedTransformedBasisValuesAtQPointsOriented, cellWeights, transformedBasisValuesAtQPointsOriented);

    fst::integrate(cellMat, transformedBasisGradsAtQPointsOriented, weightedTransformedBasisGradsAtQPointsOriented);
    fst::integrate(cellMat, transformedBasisValuesAtQPointsOriented, weightedTransformedBasisValuesAtQPointsOriented, true);
    fst::integrate(cellRhs, funAtQPoints, weightedTransformedBasisValuesAtQPointsOriented);


    // ************************************ GENERATE GRAPH **************************************

    auto ugi = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<local_ordinal_t, global_ordinal_t> >(dofManager);
    std::vector<global_ordinal_t> ownedIndices, ownedAndGhostedIndices;
    ugi->getOwnedIndices(ownedIndices);
    Teuchos::RCP<const map_t> owned_map = Teuchos::rcp(new map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedIndices,0,Teuchos::rcpFromRef(comm)));
    ugi->getOwnedAndGhostedIndices(ownedAndGhostedIndices);
    Teuchos::RCP<const map_t> ownedAndGhosted_map = Teuchos::rcp(new const map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedAndGhostedIndices,0,Teuchos::rcpFromRef(comm)));

    auto row_map = owned_map;
    auto domain_map = owned_map;
    auto range_map  = owned_map;

    auto fe_graph = Teuchos::rcp(new fe_graph_t(row_map, ownedAndGhosted_map, 8*basisCardinality));

    Teuchos::Array<global_ordinal_t> global_ids_in_row(basisCardinality);

    //auto elmtOffset = dofManager->getGIDFieldOffsetsKokkos("eblock-0_0_0",0);
    auto elmtOffset = dofManager->getGIDFieldOffsets("eblock-0_0_0",0);


    // for each element in the mesh...
    Tpetra::beginFill(*fe_graph);
    for(ordinal_type elem_id=0; elem_id<numOwnedElems; elem_id++)
    {
      // Populate global_ids_in_row:
      // - Copy the global node ids for current element into an array.
      // - Since each element's contribution is a clique, we can re-use this for
      //   each row associated with this element's contribution.
      std::vector<global_ordinal_t> elementGIDs;
      dofManager->getElementGIDs(elem_id, elementGIDs);
      for(ordinal_type node_id=0; node_id<basisCardinality; node_id++)
        global_ids_in_row[node_id] = elementGIDs[elmtOffset[node_id]];


      // Add the contributions from the current row into the graph.
      // - For example, if Element 0 contains nodes [0,1,4,5, ...] then we insert the nodes:
      //   - node 0 inserts [0, 1, 4, 5, ...]
      //   - node 1 inserts [0, 1, 4, 5, ...]
      //   - node 4 inserts [0, 1, 4, 5, ...]
      //   - node 5 inserts [0, 1, 4, 5, ...]
      for(ordinal_type node_id=0; node_id<basisCardinality; node_id++)
      {
        fe_graph->insertGlobalIndices(global_ids_in_row[node_id], global_ids_in_row());
      }
    }
    Tpetra::endFill(*fe_graph);


    // ************************************ MATRIX ASSEMBLY **************************************

    auto fe_matrix = Teuchos::rcp(new fe_matrix_t(fe_graph));
    auto rhs = Teuchos::rcp (new fe_multivector_t(domain_map, fe_graph->getImporter(), 1));

    Teuchos::Array<global_ordinal_t> column_global_ids(basisCardinality), column_local_ids(basisCardinality);  // global column ids list
    Teuchos::Array<scalar_t> column_scalar_values(basisCardinality);         // scalar values for each column

    // Loop over elements
    Tpetra::beginFill(*fe_matrix,*rhs);
    std::vector<global_ordinal_t> elementGIDs(basisCardinality);
    auto cellRhsHost = Kokkos::create_mirror_view(cellRhs);
    Kokkos::deep_copy(cellRhsHost,cellRhs);
    auto cellMatHost = Kokkos::create_mirror_view(cellMat);
    Kokkos::deep_copy(cellMatHost,cellMat);
    auto elementLIDs = ugi->getLIDs();
    auto elementLIDsHost = Kokkos::create_mirror_view(elementLIDs);
    Kokkos::deep_copy(elementLIDsHost,elementLIDs);
    for(ordinal_type elem_id=0; elem_id<numOwnedElems; elem_id++)
    {
      // Fill the global column ids array for this element
      dofManager->getElementGIDs(elem_id, elementGIDs);


      for(ordinal_type node_id=0; node_id<basisCardinality; node_id++)
        column_local_ids[node_id] = fe_matrix->getColMap()->getLocalElement(elementGIDs[elmtOffset[node_id]]);

      // For each node (row) on the current element:
      // - populate the values array
      // - add the values to the fe_matrix.

      for(ordinal_type node_id=0; node_id<basisCardinality; node_id++)
      {
        global_ordinal_t local_row_id = elementLIDsHost(elem_id, elmtOffset[node_id]);

        for(ordinal_type col_idx=0; col_idx<basisCardinality; col_idx++)
          column_scalar_values[col_idx] = cellMatHost(elem_id, node_id, col_idx);

        fe_matrix->sumIntoLocalValues(local_row_id, column_local_ids, column_scalar_values);
        rhs->sumIntoLocalValue(column_local_ids[node_id], 0, cellRhsHost(elem_id, node_id));
      }
    }

    Tpetra::endFill(*fe_matrix);
    Tpetra::endFill(*rhs);




    // ************************************ CODE VERIFICATION **************************************


    DynRankView ConstructWithLabel(basisCoeffsLI, numOwnedElems, basisCardinality);
    {
      DynRankView ConstructWithLabel(dofCoordsOriented, numOwnedElems, basisCardinality, dim);
      DynRankView ConstructWithLabel(dofCoeffsPhys, numOwnedElems, basisCardinality);

      li::getDofCoordsAndCoeffs(dofCoordsOriented,  dofCoeffsPhys, basis.getRawPtr(), POINTTYPE_EQUISPACED, elemOrts);

      //Compute physical Dof Coordinates

      DynRankView ConstructWithLabel(funAtDofPoints, numOwnedElems, basisCardinality);
      {
        DynRankView ConstructWithLabel(physDofPoints, numOwnedElems, numQPoints, dim);
        ct::mapToPhysicalFrame(physDofPoints,dofCoordsOriented,physVertexes,basis->getBaseCellTopology());
        EvalSolFunctor<DynRankView> functor;
        functor.funAtPoints = funAtDofPoints;
        functor.points = physDofPoints;
        Kokkos::parallel_for(numOwnedElems,functor);
      }

      li::getBasisCoeffs(basisCoeffsLI, funAtDofPoints, dofCoeffsPhys);
    }

    {
      Tpetra::Vector<scalar_t> sol(domain_map);
      auto basisCoeffsLIHost = Kokkos::create_mirror_view(basisCoeffsLI);
      Kokkos::deep_copy(basisCoeffsLIHost,basisCoeffsLI);
      for(ordinal_type elem_id=0; elem_id<numOwnedElems; elem_id++)
      {
        //auto elementLIDs = ugi->getElementLIDs(elem_id);
        dofManager->getElementGIDs(elem_id, elementGIDs);


        for(ordinal_type node_id=0; node_id<basisCardinality; node_id++)
        {
          global_ordinal_t gid = elementGIDs[elmtOffset[node_id]];
          if(domain_map->isNodeGlobalElement(gid))
            sol.replaceGlobalValue(gid, basisCoeffsLIHost(elem_id, node_id));
        }
      }

      fe_matrix->apply(sol, *rhs, Teuchos::NO_TRANS, -1.0, 1.0);   // b - A x
    }

    *outStream << "Residual norm : " << rhs->getVector(0)->norm2() <<std::endl;


    //row_map->describe(*outStream);


  } catch (std::exception err) {
    *outStream << " Exception\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    *outStream << "End Result: TEST FAILED = " << errorFlag << "\n";
  else
    *outStream << "End Result: TEST PASSED\n";

  return errorFlag;
}
}
}
