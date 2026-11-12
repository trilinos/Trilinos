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
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_IntegrationTools.hpp"

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_Interpolation.hpp"
#include "../../../dof-mgr/test/cartesian_topology/CartesianConnManager.cpp"

#include <Tpetra_Export.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_FECrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_Assembly_Helpers.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StackedTimer.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <array>
#include <set>
#include <random>
#include <algorithm>

#include "Thyra_VectorBase.hpp"
#include "Tpetra_MatrixFreeRowMatrix_decl.hpp"

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>

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

  struct SolveMetrics
  {
    double matrixFreeTotalTime;      // time in seconds for all matrix-free-specific operations (setup + solve)
    double assembledTotalTime;       // time in seconds for all assembled-matrix-specific operations (setup + solve)
    
    double assembledSolveTime;       // time in seconds for assembled-matrix solve (neglecting setup time)
    double matrixFreeSolveTime;      // time in seconds for matrix-free solve (neglecting setup time)
    
    double matrixFreeGemmFlops;      // estimated total operation count (multiplies and adds) involved in gemm() calls from PAMatrix
    double matrixFreeGemmTime;       // time in seconds spent in gemm() calls from PAMatrix
    double matrixFreeGemmThroughput; // estimated throughput in Gigaflops involved in gemm() calls from PAMatrix
    
    double matrixFreeExtractDiagonalTime; // time in seconds spent in PAMatrix::extractDiagonal()
    
    int assembledIterationCount;
    int matrixFreeIterationCount;
    
    bool assembledSolveSuccess;  // true indicates that the iterative solver converged
    bool matrixFreeSolveSuccess; // true indicates that the iterative solver converged
  };


using map_t = Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal>;

using local_ordinal_t = typename map_t::local_ordinal_type;
using global_ordinal_t = typename map_t::global_ordinal_type;
using node_t = typename map_t::node_type;
using fe_graph_t = Tpetra::FECrsGraph<local_ordinal_t,global_ordinal_t,node_t>;

std::pair<Teuchos::RCP<const map_t>, Teuchos::RCP<const map_t>>
buildMaps(Teuchos::RCP<panzer::DOFManager> dofManager, std::string& fieldName) {
  auto comm = dofManager->getComm();
  auto fieldPattern = dofManager->getFieldPattern(fieldName);
  auto basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(fieldPattern, true)->getIntrepidBasis();
  auto globalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(dofManager);
  std::vector<global_ordinal_t> ownedIndices, ownedAndGhostedIndices;
  globalIndexer->getOwnedIndices(ownedIndices);
  auto ownedMap = Teuchos::rcp(new map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedIndices,0,comm));
  globalIndexer->getOwnedAndGhostedIndices(ownedAndGhostedIndices);
  auto ownedAndGhostedMap = Teuchos::rcp(new const map_t(Teuchos::OrdinalTraits<global_ordinal_t>::invalid(),ownedAndGhostedIndices,0,comm));

  return std::make_pair(ownedMap, ownedAndGhostedMap);
}

Teuchos::RCP<fe_graph_t>
buildFEGraph(Teuchos::RCP<panzer::DOFManager> dofManager, std::string& fieldName, RCP<const map_t> ownedMap, RCP<const map_t> ownedAndGhostedMap) {

  auto fieldPattern = dofManager->getFieldPattern(fieldName);
  auto basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(fieldPattern, true)->getIntrepidBasis();

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
  Teuchos::RCP<fe_graph_t> feGraph = Teuchos::rcp(new fe_graph_t(ownedMap, ownedAndGhostedMap, maxNumRowEntries, ownedAndGhostedMap, Teuchos::null, ownedMap));

  auto basisCardinality = basis->getCardinality();
  Teuchos::Array<global_ordinal_t> globalIdsInRow(basisCardinality);
  // auto elmtOffsetKokkos = dofManager->getGIDFieldOffsetsKokkos(blockId,0);
  const std::string blockId = "eblock-0_0_0";
  auto elmtOffsets_host = dofManager->getGIDFieldOffsets(blockId,0);

  auto elementIds = dofManager->getElementBlock(blockId);
  int numOwnedElems = elementIds.size();

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
  return feGraph;
}

template<typename ValueType, typename DeviceType>
void fillMatrix(Teuchos::RCP<panzer::DOFManager> dofManager,
                std::string& fieldName,
                Teuchos::RCP<Intrepid2::CellGeometry<ValueType, 3, DeviceType> > & geometry,
                Kokkos::DynRankView<Intrepid2::Orientation,DeviceType> elemOrts,
                Teuchos::RCP<Tpetra::FECrsMatrix<ValueType, local_ordinal_t, global_ordinal_t,node_t>>& A_crs) {
  using scalar_t = ValueType;
  using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

  using ct = Intrepid2::CellTools<DeviceType>;
  using ots = Intrepid2::OrientationTools<DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using its = Intrepid2::IntegrationTools<DeviceType>;

  auto fieldPattern = dofManager->getFieldPattern(fieldName);
  auto basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(fieldPattern, true)->getIntrepidBasis();
  auto basisCardinality = basis->getCardinality();

  int cubDegree = 2*basis->getDegree();

  auto globalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(dofManager);

  const std::string blockId = "eblock-0_0_0";

  auto elementIds = dofManager->getElementBlock(blockId);
  int numOwnedElems = elementIds.size();

  // Get cell topology for base hexahedron
  typedef shards::CellTopology    CellTopology;
  // CellTopology topology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
  auto topology = basis->getBaseCellTopology();
  auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType, scalar_t, scalar_t>(topology.getBaseKey(),cubDegree);

  DynRankView elemsMat("elemsMat", numOwnedElems, basisCardinality, basisCardinality);
  {
    auto localFeAssemblyTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Local Finite Element Assembly")));

    // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************
    // Compute quadrature (cubature) points

    auto tensorQuadWeights = cubature->allocateCubatureWeights();
    Intrepid2::TensorPoints<scalar_t,DeviceType> tensorQuadPoints  = cubature->allocateCubaturePoints();
    cubature->getCubature(tensorQuadPoints, tensorQuadWeights);

    // compute oriented basis functions at quadrature points
    auto basisValuesAtQPoints = basis->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
    basis->getValues(basisValuesAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
    auto basisGradsAtQPoints = basis->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_GRAD);
    basis->getValues(basisGradsAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_GRAD);

    auto jacobian = geometry->allocateJacobianData(tensorQuadPoints);
    auto jacobianDet = ct::allocateJacobianDet(jacobian);
    auto jacobianInv = ct::allocateJacobianInv(jacobian);
    auto cellMeasures = geometry->allocateCellMeasure(jacobianDet, tensorQuadWeights);
    auto refData = geometry->getJacobianRefData(tensorQuadPoints);

    // compute jacobian and cell measures
    geometry->setJacobian(jacobian, tensorQuadPoints, refData);
    ct::setJacobianDet(jacobianDet, jacobian);
    ct::setJacobianInv(jacobianInv, jacobian);
    geometry->computeCellMeasure(cellMeasures, jacobianDet, tensorQuadWeights);

    // lazily-evaluated transformed values and gradients:
    auto transformedBasisValues = fst::getHGRADtransformVALUE(numOwnedElems, basisValuesAtQPoints);
    auto transformedBasisGradients = fst::getHGRADtransformGRAD(jacobianInv, basisGradsAtQPoints);

    // assemble the matrix: integrate and apply orientation
    auto integralData = its::allocateIntegralData(transformedBasisGradients, cellMeasures, transformedBasisGradients);

    bool sumInto = false;
    its::integrate(integralData, transformedBasisValues, cellMeasures, transformedBasisValues, sumInto);
    sumInto = true;
    its::integrate(integralData, transformedBasisGradients, cellMeasures, transformedBasisGradients, sumInto);

    //      {
    //        // DEBUGGING
    //        std::cout << "elemOrts(0): " << elemOrts(0) << std::endl;
    //        std::cout << "elemOrts(1): " << elemOrts(1) << std::endl;
    //      }

    ots::modifyMatrixByOrientation(elemsMat, integralData.getUnderlyingView(), elemOrts, basis.get(), basis.get());
  }

  // Loop over elements
  A_crs->beginAssembly();
  {
    auto localMatrix  = A_crs->getLocalMatrixDevice();
    auto elementLIDs = globalIndexer->getLIDs();
    auto elmtOffsetKokkos = dofManager->getGIDFieldOffsetsKokkos(blockId,0);

    Kokkos::parallel_for
      ("Assemble FE matrix",
       Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
       KOKKOS_LAMBDA (const size_t elemId) {
        // Get subviews
        auto elemMat = Kokkos::subview(elemsMat,elemId, Kokkos::ALL(), Kokkos::ALL());
        auto elemLIds  = Kokkos::subview(elementLIDs,elemId, Kokkos::ALL());

        // For each node (row) on the current element
        for (local_ordinal_t nodeId = 0; nodeId < basisCardinality; ++nodeId) {
          const local_ordinal_t localRowId = elemLIds(elmtOffsetKokkos(nodeId));

          // Force atomics on sums
          for (local_ordinal_t colId = 0; colId < basisCardinality; ++colId)
            localMatrix.sumIntoValues (localRowId, &elemLIds(elmtOffsetKokkos(colId)), 1, &(elemMat(nodeId,colId)), true, true);
        }
      });
  }
  A_crs->endAssembly();
}

template<typename ValueType, typename DeviceType>
SolveMetrics feAssemblyHex(const int &degree,
                  const local_ordinal_t &nx, const local_ordinal_t &ny, const local_ordinal_t &nz,
                  const int &px, const int &py, const int &pz,
                  const int &cubDegree,
                  const std::string &stratFileName,
                  const std::string &timingsFile,
                  const std::string &test_name,
                  const int &verbose)
{
  SolveMetrics solveMetrics;
  
  double gemmBaseFlops = Intrepid2::PAMatrix<DeviceType,double>::gemmFlopCount();
  double gemmBaseTime  = Intrepid2::PAMatrix<DeviceType,double>::gemmTimeSeconds();
  
  // ************************************ GET INPUTS **************************************
  constexpr local_ordinal_t dim = 3;
  constexpr int bx=1, by=1, bz=1;  //blocks on each process. Here we assume there is only one block per process

  // host_memory/execution/mirror_space deprecated for kokkos@3.7.00, removed after release
  // see https://github.com/kokkos/kokkos/pull/3973
  using exec_space = typename DeviceType::execution_space;
  using mem_space = typename DeviceType::memory_space;
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
      DeviceType,
      Kokkos::Device<do_not_use_host_execution_space,
                     do_not_use_host_memory_space>>;

  using HostSpaceType = typename host_mirror_space::execution_space;

  using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

  // using map_t = Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal>;

  // using local_ordinal_t = typename map_t::local_ordinal_type;
  // using global_ordinal_t = typename map_t::global_ordinal_type;
  // using node_t = typename map_t::node_type;
  // using fe_graph_t = Tpetra::FECrsGraph<local_ordinal_t,global_ordinal_t,node_t>;
  using scalar_t = ValueType;
  using fe_matrix_t = Tpetra::FECrsMatrix<scalar_t, local_ordinal_t, global_ordinal_t,node_t>;
  using fe_multivector_t = Tpetra::FEMultiVector<scalar_t, local_ordinal_t, global_ordinal_t,node_t>;
  using vector_t = Tpetra::Vector<scalar_t, local_ordinal_t, global_ordinal_t,node_t>;
  using multivector_t = Tpetra::MultiVector<scalar_t, local_ordinal_t, global_ordinal_t,node_t>;
  using rowmatrix_t = Tpetra::RowMatrix<scalar_t, local_ordinal_t, global_ordinal_t,node_t>;

  using DynRankViewGId = Kokkos::DynRankView<global_ordinal_t,DeviceType>;

  using ct = Intrepid2::CellTools<DeviceType>;
  using ots = Intrepid2::OrientationTools<DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using li = Intrepid2::LagrangianInterpolation<DeviceType>;
  using its = Intrepid2::IntegrationTools<DeviceType>;

  int errorFlag = 0;

  double standardAssemblyTotalTime = 0;
  double matrixFreeTotalTime = 0;
  double standardAssemblySolveTime = 0;
  double matrixFreeSolveTime = 0;
  double gemmThroughput = 0;


#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::RCP<const Teuchos::MpiComm<int> > comm
    = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

  //output stream/file
  Teuchos::RCP<Teuchos::FancyOStream> outStream = getFancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;

  stacked_timer = Teuchos::rcp(new Teuchos::StackedTimer("Matrix-free driver"));;
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  try {
    outStream = ((comm->getRank () == 0) && verbose) ?
      getFancyOStream(Teuchos::rcpFromRef (std::cout)) :
      getFancyOStream(Teuchos::rcp (new Teuchos::oblackholestream ()));

    *outStream << "DeviceSpace::  "; DeviceType().print_configuration(*outStream, false);
    *outStream << "HostSpace::    "; HostSpaceType().print_configuration(*outStream, false);
    *outStream << "\n";
    stacked_timer->setVerboseOstream(outStream);

    Teuchos::RCP<Teuchos::ParameterList> strat_params = Teuchos::rcp(new Teuchos::ParameterList("Stratimikos parameters"));
    Teuchos::RCP<Teuchos::ParameterList> strat_params_mf = Teuchos::rcp(new Teuchos::ParameterList("Stratimikos parameters"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(stratFileName, strat_params.ptr(), *comm);
    Teuchos::updateParametersFromXmlFileAndBroadcast(stratFileName, strat_params_mf.ptr(), *comm);

    local_ordinal_t numOwnedElems = nx*ny*nz;

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology topology(shards::getCellTopologyData<shards::Hexahedron<8> >() );

    // Get dimensions
    int numNodesPerElem = topology.getNodeCount();

    Teuchos::RCP<panzer::unit_test::CartesianConnManager> connManager;
    Teuchos::RCP<panzer::DOFManager> dofManager;
    Teuchos::RCP< Intrepid2::Basis<DeviceType, scalar_t,scalar_t> > basis;

    const std::string blockId = "eblock-0_0_0";
    DynRankView physVertices;
    Kokkos::DynRankView<Intrepid2::Orientation,DeviceType> elemOrts("elemOrts", numOwnedElems);

    Teuchos::RCP<const map_t> ownedMap, ownedAndGhostedMap;
    Teuchos::RCP<fe_graph_t> feGraph;
    Teuchos::RCP<fe_matrix_t> A_crs;
    Teuchos::RCP<fe_multivector_t> b;

    bool setupMultigrid = true;
    std::vector<int> pCoarsenSchedule;
    // Use all degrees between 1 and degree for now.
//    for (int deg = degree-1; deg>0; --deg) {
//      pCoarsenSchedule.push_back(deg);
//    }
    int pCoarse = std::max(degree/2, 1);
    pCoarsenSchedule.push_back(pCoarse);
    while (pCoarse > 1) {
      pCoarse = std::max(pCoarse/2, 1);
      pCoarsenSchedule.push_back(pCoarse);
    }
    std::vector<Teuchos::RCP<panzer::DOFManager> > dofManagers;

    {
      auto meshTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Mesh Generation")));

      // *********************************** MESH TOPOLOGY **********************************

      // build the topology
      connManager = Teuchos::rcp(new panzer::unit_test::CartesianConnManager);
      connManager->initialize(*comm,
                              global_ordinal_t(nx*px),
                              global_ordinal_t(ny*py),
                              global_ordinal_t(nz*pz),
                              px,py,pz,bx,by,bz);

      // *********************************** COMPUTE GLOBAL IDs OF VERTICES AND DOFs  ************************************

      // build the dof manager, and assocaite with the topology
      dofManager = Teuchos::rcp(new panzer::DOFManager);
      dofManager->setConnManager(connManager, *comm->getRawMpiComm());

      // add solution field to the element block
      using CG_DNBasis = Intrepid2::DerivedNodalBasisFamily<DeviceType,scalar_t,scalar_t>;
      basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_HEX(degree));
      auto basisCardinality = basis->getCardinality();
      Teuchos::RCP<panzer::Intrepid2FieldPattern> fePattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
      std::string fieldName = "order"+std::to_string(degree);
      dofManager->addField(fieldName,fePattern);

      if (setupMultigrid) {
        for (auto it = pCoarsenSchedule.begin(); it != pCoarsenSchedule.end(); ++it) {
          auto dofManager = Teuchos::rcp(new panzer::DOFManager);
          dofManager->setConnManager(connManager, *comm->getRawMpiComm());
          auto basis = Teuchos::rcp(new typename CG_DNBasis::HGRAD_HEX(*it));
          Teuchos::RCP<panzer::Intrepid2FieldPattern> fePattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
          std::string fieldName = "order"+std::to_string(*it);
          dofManager->addField(fieldName, fePattern);
          dofManager->buildGlobalUnknowns();
          dofManagers.push_back(dofManager);
        }
      }

      // try to get them all synced up
      comm->barrier();

      dofManager->buildGlobalUnknowns();

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
      physVertices = DynRankView("physVertices", numOwnedElems, numNodesPerElem, dim);
      {
        auto physVerticesHost = Kokkos::create_mirror_view(physVertices);
        DynRankView ConstructWithLabel(refVertices, numNodesPerElem, dim);
        Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType,scalar_t,scalar_t> hexaLinearBasis;
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
              physVerticesHost(i,j,k) = offset[k]+h[k]/2.0*refVerticesHost(j,k);
          }
        }
        Kokkos::deep_copy(physVertices, physVerticesHost);
      }
    }

    // *********************************** COMPUTE ELEMENTS' ORIENTATION BASED ON GLOBAL IDs  ************************************

    DynRankViewGId ConstructWithLabel(elemNodesGID, numOwnedElems, numNodesPerElem);
    {
      auto orientationsTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Orientations")));

      //compute global ids of element vertices
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
      ots::getOrientation(elemOrts, elemNodesGID, topology);
    }

    Kokkos::DynRankView<int,DeviceType> emptyView;
    // In principle we could pass elemNodesGID instead of emptyView, but then physVertices should be a 2d array
    // indexed by global ids instead of local ids, which is not efficient
    auto geometry = Teuchos::rcp(new Intrepid2::CellGeometry<scalar_t, dim, DeviceType>(topology, emptyView, physVertices));
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType, scalar_t, scalar_t>(topology.getBaseKey(),cubDegree);

    auto crsTotalTimer = Teuchos::TimeMonitor::getNewTimer("Crs-Specific Total Time");
    crsTotalTimer->start();
    auto basisCardinality = basis->getCardinality();
    DynRankView elemsMat("elemsMat", numOwnedElems, basisCardinality, basisCardinality);
    DynRankView elemsRHS("elemsRHS", numOwnedElems, basisCardinality);
    DynRankView elemsRHSTmp("elemsRHS", numOwnedElems, basisCardinality);
    {
      auto localFeAssemblyTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Local Finite Element Assembly")));

      // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************
      // Compute quadrature (cubature) points


      auto numQPoints = cubature->getNumPoints();
      auto tensorQuadWeights = cubature->allocateCubatureWeights();
      Intrepid2::TensorPoints<scalar_t,DeviceType> tensorQuadPoints  = cubature->allocateCubaturePoints();
      cubature->getCubature(tensorQuadPoints, tensorQuadWeights);

      // compute oriented basis functions at quadrature points
      auto basisCardinality = basis->getCardinality();
      auto basisValuesAtQPoints = basis->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
      basis->getValues(basisValuesAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
      auto basisGradsAtQPoints = basis->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_GRAD);
      basis->getValues(basisGradsAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_GRAD);

      auto jacobian = geometry->allocateJacobianData(tensorQuadPoints);
      auto jacobianDet = ct::allocateJacobianDet(jacobian);
      auto jacobianInv = ct::allocateJacobianInv(jacobian);
      auto cellMeasures = geometry->allocateCellMeasure(jacobianDet, tensorQuadWeights);
      auto refData = geometry->getJacobianRefData(tensorQuadPoints);

      // compute jacobian and cell measures
      geometry->setJacobian(jacobian, tensorQuadPoints, refData);
      ct::setJacobianDet(jacobianDet, jacobian);
      ct::setJacobianInv(jacobianInv, jacobian);
      geometry->computeCellMeasure(cellMeasures, jacobianDet, tensorQuadWeights);

      // lazily-evaluated transformed values and gradients:
      auto transformedBasisValues = fst::getHGRADtransformVALUE(numOwnedElems, basisValuesAtQPoints);
      auto transformedBasisGradients = fst::getHGRADtransformGRAD(jacobianInv, basisGradsAtQPoints);

      // assemble the matrix: integrate and apply orientation
      auto integralData = its::allocateIntegralData(transformedBasisGradients, cellMeasures, transformedBasisGradients);

      bool sumInto = false;
      its::integrate(integralData, transformedBasisValues, cellMeasures, transformedBasisValues, sumInto);
      sumInto = true;
      its::integrate(integralData, transformedBasisGradients, cellMeasures, transformedBasisGradients, sumInto);

//      {
//        // DEBUGGING
//        std::cout << "elemOrts(0): " << elemOrts(0) << std::endl;
//        std::cout << "elemOrts(1): " << elemOrts(1) << std::endl;
//      }

      ots::modifyMatrixByOrientation(elemsMat, integralData.getUnderlyingView(), elemOrts, basis.get(), basis.get());
      crsTotalTimer->stop();

    // ************************************ ASSEMBLY OF LOCAL ELEMENT RHS VECTORS **************************************

    //Compute physical points where to evaluate the function
    DynRankView ConstructWithLabel(funAtQPoints, numOwnedElems, numQPoints);
    {
        DynRankView quadPoints("quadPoints", tensorQuadPoints.extent_int(0), tensorQuadPoints.extent_int(1));
        tensorQuadPoints.copyPointsContainer(quadPoints, tensorQuadPoints);
        DynRankView ConstructWithLabel(physQPoints, numOwnedElems, numQPoints, dim);
        ct::mapToPhysicalFrame(physQPoints,quadPoints,physVertices,basis->getBaseCellTopology());
        EvalRhsFunctor<DynRankView> functor;
        functor.funAtPoints = funAtQPoints;
        functor.points = physQPoints;
        Kokkos::parallel_for("loop for evaluating the rhs at quadrature points",
          Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
          functor);
    }

    //compute the weighted basis functions
    DynRankView ConstructWithLabel(weightedTransformedBasisValuesAtQPoints, numOwnedElems, basisCardinality, numQPoints);
    auto policy = Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<3>>({0,0,0},{numOwnedElems, basisCardinality, numQPoints});
    Kokkos::parallel_for("compute weighted basis", policy,
    KOKKOS_LAMBDA (const int &cell, const int &field, const int &point) {
      weightedTransformedBasisValuesAtQPoints(cell,field,point) = transformedBasisValues(cell,field,point)*cellMeasures(cell,point);
    });


    // assemble the rhs: integrate and apply orientation
    fst::integrate(elemsRHSTmp, funAtQPoints, weightedTransformedBasisValuesAtQPoints);

    ots::modifyBasisByOrientation(elemsRHS, elemsRHSTmp, elemOrts, basis.get());
    }

    // ************************************ GENERATE GRAPH **************************************
    auto globalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(dofManager);
    {
      auto graphGenerationTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Graph Generation")));
      std::string fieldName = "order"+std::to_string(degree);
      auto maps = buildMaps(dofManager, fieldName);
      ownedMap = std::get<0>(maps);
      ownedAndGhostedMap = std::get<0>(maps);
      feGraph = buildFEGraph(dofManager, fieldName, ownedMap, ownedAndGhostedMap);
    }

    // ************************************ MATRIX ASSEMBLY **************************************

    crsTotalTimer->start();
    {
      auto matrixAndRhsAllocationTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Allocation of Matrix and Rhs")));

      A_crs = Teuchos::rcp(new fe_matrix_t(feGraph));
      b = Teuchos::rcp (new fe_multivector_t(ownedMap, feGraph->getImporter(), 1));
    }

    {
      auto matrixAndRhsFillTimer =  Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Fill of Matrix and Rhs")));

      //fill matrix
      // Loop over elements
      Tpetra::beginAssembly(*A_crs, *b);
      {
        auto localMatrix  = A_crs->getLocalMatrixDevice();
        auto localRHS     = b->getLocalViewDevice(Tpetra::Access::ReadWrite);
        auto elementLIDs = globalIndexer->getLIDs();
        auto elmtOffsetKokkos = dofManager->getGIDFieldOffsetsKokkos(blockId,0);

        Kokkos::parallel_for
          ("Assemble FE matrix and right-hand side",
           Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
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
      }
      Tpetra::endAssembly(*A_crs, *b);
    }
    crsTotalTimer->stop();

    auto mfTotalTimer = Teuchos::TimeMonitor::getNewTimer("MF-Specific Total Time");
    mfTotalTimer->start();
    Teuchos::RCP<rowmatrix_t> A_mf;
    {
      A_mf = Teuchos::rcp(new Tpetra::MatrixFreeRowMatrix<scalar_t, local_ordinal_t, global_ordinal_t>(ownedMap, ownedAndGhostedMap,
                                                                                                       basis, geometry,
                                                                                                       cubature, elemOrts,
                                                                                                       dofManager, globalIndexer));
    }
    mfTotalTimer->stop();

    // ************************************ MULTIGRID SETUP   **************************************

    Teuchos::RCP<Teuchos::ParameterList> muelu_user_data = Teuchos::rcp(new Teuchos::ParameterList("user data"));
    if (setupMultigrid) {

      Teuchos::RCP<panzer::DOFManager> fineDofManager = dofManager;
      Teuchos::RCP<panzer::DOFManager> coarseDofManager;
      int fine_degree = degree;
      std::string fineFieldName = "order"+std::to_string(degree);
      std::string coarseFieldName;
      Teuchos::RCP<const map_t> fineOwnedMap = ownedMap;
      Teuchos::RCP<const map_t> coarseOwnedMap;
      int mueluLevel = 0;
      int maxCoarseGridSize = Teuchos::as<int>(fineOwnedMap->getGlobalNumElements());

      for (auto it = pCoarsenSchedule.begin(); it != pCoarsenSchedule.end(); ++it) {

        int coarse_degree = *it;

        coarseFieldName = "order"+std::to_string(coarse_degree);
        coarseDofManager = dofManagers[mueluLevel];

        auto coarseMaps = buildMaps(coarseDofManager, coarseFieldName);
        coarseOwnedMap = std::get<0>(coarseMaps);
        auto coarseOwnedAndGhostedMap = std::get<1>(coarseMaps);
        maxCoarseGridSize = std::min(Teuchos::as<int>(fineOwnedMap->getGlobalNumElements()), maxCoarseGridSize);

        // grid transfer
        {
          Teuchos::TimeMonitor pTimer =  *Teuchos::TimeMonitor::getNewTimer("Assemble matrix prolongator between orders " +std::to_string(coarse_degree) +" and "+std::to_string(fine_degree));
          RCP<const Thyra::LinearOpBase<scalar_t> > P = panzer::buildInterpolation(connManager,
                                                                                   coarseDofManager, fineDofManager,
                                                                                   coarseFieldName, fineFieldName,
                                                                                   Intrepid2::OPERATOR_VALUE,
                                                                                   /*worksetSize=*/1000, /*forceVectorial=*/false,
                                                                                   /*useTpetra=*/true, /*matrixFree=*/false);
          strat_params->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("P", P);
        }
        {
          Teuchos::TimeMonitor pTimer =  *Teuchos::TimeMonitor::getNewTimer("Assemble matrix-free prolongator between orders " +std::to_string(coarse_degree) +" and "+std::to_string(fine_degree));
          // I set the worksetSize to 1. On serial backend this should not matter in terms of computations, but it saves a lot of memory.
          // Need to revisit this for parallel backends.
          RCP<const Thyra::LinearOpBase<scalar_t> > P = panzer::buildInterpolation(connManager,
                                                                                   coarseDofManager, fineDofManager,
                                                                                   coarseFieldName, fineFieldName,
                                                                                   Intrepid2::OPERATOR_VALUE,
                                                                                   /*worksetSize=*/1, /*forceVectorial=*/false,
                                                                                   /*useTpetra=*/true, /*matrixFree=*/true);
          strat_params_mf->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("P", P);
        }

        // coarse matrix
        auto thyra_vs_fine = Thyra::tpetraVectorSpace<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(fineOwnedMap);
        auto thyra_vs_coarse = Thyra::tpetraVectorSpace<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(coarseOwnedMap);
        Teuchos::RCP<const Thyra::LinearOpBase<scalar_t> > thyra_Ac, thyra_Ac_mf;
        {
          // matrix assembly
          Teuchos::TimeMonitor acTimer =  *Teuchos::TimeMonitor::getNewTimer("Assemble matrix for order " +std::to_string(coarse_degree));
          auto feGraph = buildFEGraph(coarseDofManager, coarseFieldName, coarseOwnedMap, coarseOwnedAndGhostedMap);
          auto Ac = Teuchos::rcp(new fe_matrix_t(feGraph));
          fillMatrix<ValueType, DeviceType>(coarseDofManager, coarseFieldName, geometry, elemOrts, Ac);
          thyra_Ac = Thyra::tpetraLinearOp<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(thyra_vs_coarse, thyra_vs_coarse, Ac);
        }

        if (coarse_degree > 1) {
          Teuchos::TimeMonitor acTimer =  *Teuchos::TimeMonitor::getNewTimer("Assemble matrix-free operator for order " +std::to_string(coarse_degree));
          // matrix-free assembly
          auto fieldPattern = coarseDofManager->getFieldPattern(coarseFieldName);
          auto coarseBasis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(fieldPattern, true)->getIntrepidBasis();
          auto coarseCubDegree = 2*coarseBasis->getDegree();
          auto coarseCubature = Intrepid2::DefaultCubatureFactory::create<DeviceType, scalar_t, scalar_t>(topology.getBaseKey(),coarseCubDegree);
          auto coarseGlobalIndexer = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer >(coarseDofManager);
          auto Ac = Teuchos::rcp(new Tpetra::MatrixFreeRowMatrix<scalar_t, local_ordinal_t, global_ordinal_t>(coarseOwnedMap, coarseOwnedAndGhostedMap,
                                                                                                              coarseBasis, geometry,
                                                                                                              coarseCubature, elemOrts,
                                                                                                              coarseDofManager, coarseGlobalIndexer));
          thyra_Ac_mf = Thyra::tpetraLinearOp<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(thyra_vs_coarse, thyra_vs_coarse, Ac);
        }
        strat_params->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("A", thyra_Ac);
        if (coarse_degree == 1) {
          strat_params_mf->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("A", thyra_Ac);
        } else {
          strat_params_mf->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("A", thyra_Ac_mf);
        }

        if (coarse_degree == 1) {
          { // Set up nullspace for problem
            auto coarseNullspace = Teuchos::rcp(new multivector_t(coarseOwnedMap, 1));
            coarseNullspace->putScalar(1.);
            auto thyra_coarseNullspace = Thyra::createMultiVector(coarseNullspace);
            strat_params->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("Nullspace", thyra_coarseNullspace);
          }
          { // Set up coordinates for problem

            // auto coarseCoords = Teuchos::rcp(new multivector_t(coarseOwnedMap, 1));
            auto fieldPattern = coarseDofManager->getFieldPattern(coarseFieldName);
            auto intrepid2FieldPattern = Teuchos::rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(fieldPattern, true);

            TEUCHOS_ASSERT(intrepid2FieldPattern->supportsInterpolatoryCoordinates());
            std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > data;
            Kokkos::DynRankView<double,PHX::Device> & fieldData = data[blockId];
            intrepid2FieldPattern->getInterpolatoryCoordinates(physVertices, fieldData);

            auto atv = Teuchos::rcp(new panzer::ArrayToFieldVector(coarseDofManager));
            auto coarseCoords = atv->template getDataVector<double>(std::string(coarseFieldName), data);

            auto thyra_coarseCoords = Thyra::createMultiVector(coarseCoords);
            strat_params->sublist("Preconditioner Types").sublist("MueLu").sublist("level " + std::to_string(mueluLevel+1)+" user data").set("Coordinates", thyra_coarseCoords);
          }
        }

        ++mueluLevel;
        fineDofManager = coarseDofManager;
        fine_degree = coarse_degree;
        fineFieldName = coarseFieldName;
        fineOwnedMap = coarseOwnedMap;
      }
      if (strat_params->sublist("Preconditioner Types").sublist("MueLu").isParameter("coarse: max size"))
        maxCoarseGridSize = std::min(strat_params->sublist("Preconditioner Types").sublist("MueLu").get<int>("coarse: max size"), maxCoarseGridSize);
      strat_params->sublist("Preconditioner Types").sublist("MueLu").set("coarse: max size", maxCoarseGridSize);
    }

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
        ct::mapToPhysicalFrame(physDofPoints,dofCoords,physVertices,basis->getBaseCellTopology());
        EvalSolFunctor<DynRankView> functor;
        functor.funAtPoints = funAtDofPoints;
        functor.points = physDofPoints;
        Kokkos::parallel_for("loop for evaluating the function at DoF points",
          Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
          functor);
        Kokkos::fence(); //make sure that funAtDofPoints has been evaluated
      }

      li::getBasisCoeffs(basisCoeffsLI, funAtDofPoints, basis.getRawPtr(), elemOrts);
    }

    Teuchos::RCP<vector_t> x;
    {
      Teuchos::TimeMonitor vTimer1 =  *Teuchos::TimeMonitor::getNewTimer("Verification, assemble solution");

      std::vector<global_ordinal_t> elementGIDs(basisCardinality);
      auto elmtOffsets_host = dofManager->getGIDFieldOffsets(blockId,0);

      x = Teuchos::rcp(new vector_t(ownedMap)); //solution
      auto basisCoeffsLIHost = Kokkos::create_mirror_view(basisCoeffsLI);
      Kokkos::deep_copy(basisCoeffsLIHost,basisCoeffsLI);
      for(int elemId=0; elemId<numOwnedElems; elemId++)
      {
        dofManager->getElementGIDs(elemId, elementGIDs);


        for(int nodeId=0; nodeId<basisCardinality; nodeId++)
        {
          global_ordinal_t gid = elementGIDs[elmtOffsets_host[nodeId]];
          if(ownedMap->isNodeGlobalElement(gid))
            x->replaceGlobalValue(gid, basisCoeffsLIHost(elemId, nodeId));
        }
      }
    }

    // Tpetra::MatrixMarket::Writer<crs_t>::writeDenseFile("b", *b);

    Teuchos::RCP<multivector_t> residual_crs = Teuchos::rcp(new multivector_t(x->getMap(), x->getNumVectors()));
    {
      Teuchos::TimeMonitor applyCrsTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix apply");
      A_crs->apply(*x, *residual_crs, Teuchos::NO_TRANS);
    }
    // Tpetra::MatrixMarket::Writer<crs_t>::writeDenseFile("b_crs", *residual_crs);
    {
      Teuchos::TimeMonitor applyCrsTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix apply");
      Tpetra::deep_copy(*residual_crs, *b);
      A_crs->apply(*x, *residual_crs, Teuchos::NO_TRANS, -1.0, 1.0);   // b - A x
    }
    // Tpetra::MatrixMarket::Writer<crs_t>::writeDenseFile("residual_crs", *residual_crs);

    double res_crs_l2_norm = residual_crs->getVector(0)->norm2();
    if((degree >= 4) && (res_crs_l2_norm > 1e-10)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "Residual norm should be close to machine eps, but it is instead: " << res_crs_l2_norm <<  "\n";
    }
    else {
      *outStream << "Residual matrix l2 norm : " << res_crs_l2_norm << "\n";
    }

    Teuchos::RCP<multivector_t> residual_mf = Teuchos::rcp(new multivector_t(x->getMap(), x->getNumVectors()));
    {
      Teuchos::TimeMonitor applyMFTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix-free apply");
      A_mf->apply(*x, *residual_mf, Teuchos::NO_TRANS);
    }
    // Tpetra::MatrixMarket::Writer<crs_t>::writeDenseFile("b_mf", *residual_mf);
    {
      Teuchos::TimeMonitor applyMFTimer = *Teuchos::TimeMonitor::getNewTimer("Matrix-free apply");
      Tpetra::deep_copy(*residual_mf, *b);
      A_mf->apply(*x, *residual_mf, Teuchos::NO_TRANS, -1.0, 1.0);
    }
    // Tpetra::MatrixMarket::Writer<crs_t>::writeDenseFile("residual_mf", *residual_mf);

    double res_mf_l2_norm = residual_mf->getVector(0)->norm2();
    if((degree >= 4) && (res_mf_l2_norm > 1e-10)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "Residual norm should be close to machine eps, but it is instead: " << res_mf_l2_norm <<  "\n";
    }
    else {
      *outStream << "Residual matrix-free l2 norm : " << res_mf_l2_norm << "\n";
    }

    {
      residual_crs->update(-1.0, *residual_mf, 1.0);
      double diff_l2_norm = residual_crs->getVector(0)->norm2();
      *outStream << "Difference matrix vs matrix-free l2 norm : " << diff_l2_norm << "\n";
    }


    // ************************************ SOLVE LINEAR SYSTEMS **************************************
    {
      *outStream << *strat_params << std::endl;
      Stratimikos::LinearSolverBuilder<scalar_t> linearSolverBuilder;
      Stratimikos::enableMueLu<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(linearSolverBuilder);

      auto thyra_vs = Thyra::tpetraVectorSpace<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(ownedMap);
      auto thyra_b = Thyra::constTpetraVector<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(thyra_vs, b->getVector(0));

      {
        crsTotalTimer->start();
        linearSolverBuilder.setParameterList(strat_params);
        Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<scalar_t> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);
        {
          Teuchos::RCP<Thyra::VectorBase<scalar_t> > thyra_x_crs;
          Teuchos::RCP<const Thyra::LinearOpBase<scalar_t> > thyra_A_crs;
          Teuchos::RCP<Thyra::LinearOpWithSolveBase<scalar_t> > thyra_inverse_A_crs;

          {
            Teuchos::TimeMonitor mfSolveTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix solver setup");
            thyra_x_crs = Thyra::createMember(*thyra_vs);
            thyra_x_crs->assign(0.);
            thyra_A_crs = Thyra::tpetraLinearOp<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(thyra_vs, thyra_vs, A_crs);
            thyra_inverse_A_crs = Thyra::linearOpWithSolve(*lowsFactory, thyra_A_crs);
          }

          crsTotalTimer->stop();
          double totalTimeBeforeSolve = crsTotalTimer->totalElapsedTime();
          crsTotalTimer->start();
          {
            Teuchos::TimeMonitor mfSolveTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix solve");
            Thyra::SolveStatus<scalar_t> status = Thyra::solve<scalar_t>(*thyra_inverse_A_crs, Thyra::NOTRANS, *thyra_b, thyra_x_crs.ptr());
            bool solveConverged = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
            solveMetrics.assembledIterationCount = status.extraParameters->template get<int>("Iteration Count", -1);
            solveMetrics.assembledSolveSuccess = solveConverged;
            errorFlag += !solveConverged;
          }
          crsTotalTimer->stop();
          standardAssemblySolveTime = crsTotalTimer->totalElapsedTime() - totalTimeBeforeSolve;
        }
      }

      {
        linearSolverBuilder.setParameterList(strat_params_mf);
        Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<scalar_t> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);
        mfTotalTimer->start();
        {
          Teuchos::RCP<Thyra::VectorBase<scalar_t> > thyra_x_mf;
          Teuchos::RCP<const Thyra::LinearOpBase<scalar_t> > thyra_A_mf;
          Teuchos::RCP<Thyra::LinearOpWithSolveBase<scalar_t> > thyra_inverse_A_mf;

          {
            Teuchos::TimeMonitor mfSolveTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix-free solver setup");
            thyra_x_mf = Thyra::createMember<scalar_t>(*thyra_vs);
            thyra_x_mf->assign(0.);
            thyra_A_mf = Thyra::tpetraLinearOp<scalar_t, local_ordinal_t, global_ordinal_t, node_t>(thyra_vs, thyra_vs, A_mf);
            thyra_inverse_A_mf = Thyra::linearOpWithSolve(*lowsFactory, thyra_A_mf);
          }

          mfTotalTimer->stop();
          double totalTimeBeforeSolve = mfTotalTimer->totalElapsedTime();
          mfTotalTimer->start();
          {
            Teuchos::TimeMonitor mfSolveTimer =  *Teuchos::TimeMonitor::getNewTimer("Matrix-free solve");
            Thyra::SolveStatus<scalar_t> status = Thyra::solve<scalar_t>(*thyra_inverse_A_mf, Thyra::NOTRANS, *thyra_b, thyra_x_mf.ptr());
            bool solveConverged = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
            solveMetrics.matrixFreeSolveSuccess = solveConverged;
            
            solveMetrics.matrixFreeIterationCount = status.extraParameters->template get<int>("Iteration Count", -1);
            errorFlag += !solveConverged;
          }
          mfTotalTimer->stop();
          matrixFreeSolveTime = mfTotalTimer->totalElapsedTime() - totalTimeBeforeSolve;
        }
      }
    }
    standardAssemblyTotalTime = crsTotalTimer->totalElapsedTime();
    matrixFreeTotalTime = mfTotalTimer->totalElapsedTime();
  } catch (const std::exception & err) {
    *outStream << " Exception\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  gemmThroughput = Intrepid2::PAMatrix<DeviceType,double>::gemmThroughputGFlops(gemmBaseFlops, gemmBaseTime);

  stacked_timer->stop("Matrix-free driver");
  
  {
    std::string computeDiagonalTimerLabel = "matrix-free computeDiagonal (not really matrix-free yet)";
    Teuchos::stat_map_type statData;
    std::vector<std::string> statNames;
    Teuchos::TimeMonitor::computeGlobalTimerStatistics(statData, statNames, Teuchos::Intersection, computeDiagonalTimerLabel);
    
    if (statData[computeDiagonalTimerLabel].size() > 0)
    {
      const double timeInSeconds = statData[computeDiagonalTimerLabel][0].first;
      solveMetrics.matrixFreeExtractDiagonalTime = timeInSeconds;
    }
    else
    {
      solveMetrics.matrixFreeExtractDiagonalTime = -1;
    }
  }
  
  solveMetrics.assembledTotalTime = standardAssemblyTotalTime;
  solveMetrics.assembledSolveTime = standardAssemblySolveTime;
  
  solveMetrics.matrixFreeTotalTime = matrixFreeTotalTime;
  solveMetrics.matrixFreeSolveTime = matrixFreeSolveTime;
  
  solveMetrics.matrixFreeGemmFlops      = Intrepid2::PAMatrix<DeviceType,double>::gemmFlopCount() - gemmBaseFlops;
  solveMetrics.matrixFreeGemmTime       = Intrepid2::PAMatrix<DeviceType,double>::gemmTimeSeconds();
  solveMetrics.matrixFreeGemmThroughput = gemmThroughput;
  
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = options.output_histogram = options.output_minmax = true;
  stacked_timer->report(*outStream, comm, options);
  auto xmlOut = stacked_timer->reportWatchrXML(test_name + ' ' + std::to_string(comm->getSize()) + " ranks", comm);
  if(xmlOut.length())
    *outStream << "\nAlso created Watchr performance report " << xmlOut << '\n';

  *outStream << "PAMatrix GEMM Throughput: " << solveMetrics.matrixFreeGemmThroughput << " GFlops.\n";
  
  *outStream << "Standard Assembly solve: " << solveMetrics.assembledSolveTime  << " seconds." << std::endl;
  *outStream << "Matrix-Free solve:       " << solveMetrics.matrixFreeSolveTime << " seconds." << std::endl;
  *outStream << "Solve speedup:           " << solveMetrics.assembledSolveTime / solveMetrics.matrixFreeSolveTime << std::endl << std::endl;
  
  *outStream << "Standard Assembly total: " << solveMetrics.assembledTotalTime  << " seconds." << std::endl;
  *outStream << "Matrix-Free total:       " << solveMetrics.matrixFreeTotalTime << " seconds." << std::endl;
  *outStream << "Overall speedup:         " << solveMetrics.assembledTotalTime / solveMetrics.matrixFreeTotalTime << std::endl << std::endl;
  
  *outStream << "Time to compute the diagonal for matrix-free (part of setup) " << solveMetrics.matrixFreeExtractDiagonalTime << " seconds." << std::endl << std::endl;
  
  if (! solveMetrics.assembledSolveSuccess)
  {
    *outStream << "Assembled-matrix solve did NOT converge after " << solveMetrics.assembledIterationCount << " iterations.\n";
  }
  else
  {
    *outStream << "Assembled-matrix solve converged after " << solveMetrics.assembledIterationCount << " iterations.\n";
  }
  if (! solveMetrics.matrixFreeSolveSuccess)
  {
    *outStream << "Matrix-free solve did NOT converge after " << solveMetrics.matrixFreeIterationCount << " iterations.\n";
  }
  else
  {
    *outStream << "Matrix-free solve converged after " << solveMetrics.matrixFreeIterationCount << " iterations.\n";
  }
  
  return solveMetrics;
}

template<typename ValueType, typename DeviceType>
int feAssemblyHex(int argc, char *argv[]) {
  Teuchos::RCP<const Teuchos::MpiComm<int> > comm
    = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

  // ************************************ GET INPUTS **************************************
  constexpr local_ordinal_t dim = 3;
  int degree = 4;
  local_ordinal_t nx = 2;
  local_ordinal_t ny            = nx;
  local_ordinal_t nz            = nx;
  int np   = comm->getSize(); // number of processors
  int px = std::cbrt(np); while(np%px!=0) --px;
  int py = std::sqrt(np/px); while(np%py!=0) --py;
  int pz = np/(px*py);
  constexpr int bx=1, by=1, bz=1;  //blocks on each process. Here we assume there is only one block per process

  int errorFlag = 0;
  int verbose = 1;

  std::string timingsFile = "";
  std::string test_name = "Matrix-free driver";

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
  clp.setOption("test-name", &test_name, "Name of test (for Watchr output)");
  clp.setOption("verbose", &verbose);
  std::string stratFileName = "stratimikos.xml";
  clp.setOption("stratimikosParams", &stratFileName);
  auto cmdResult = clp.parse(argc,argv);
  cubDegree = (cubDegree == -1) ? 2*degree : cubDegree;

  SolveMetrics solveMetrics;
  
  if(cmdResult!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
  {
    clp.printHelpMessage(argv[0], std::cout);
    errorFlag++;
  }
  else
  {
    solveMetrics = feAssemblyHex<ValueType, DeviceType>(degree, nx, ny, nz,
                                                        px, py, pz,
                                                        cubDegree,
                                                        stratFileName, timingsFile, test_name,
                                                        verbose);
    errorFlag += !(solveMetrics.assembledSolveSuccess && solveMetrics.matrixFreeSolveSuccess);
  }

  Teuchos::RCP<Teuchos::FancyOStream> outStream = getFancyOStream(Teuchos::rcpFromRef(std::cout));

  if (errorFlag != 0)
    *outStream << "End Result: TEST FAILED = " << errorFlag << "\n";
  else
    *outStream << "End Result: TEST PASSED\n";

  return errorFlag;
}
}
}
