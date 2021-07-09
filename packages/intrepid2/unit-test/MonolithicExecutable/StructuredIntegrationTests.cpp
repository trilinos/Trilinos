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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities.
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_Kernels.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_TensorArgumentIterator.hpp>
#include <Intrepid2_TestUtils.hpp>

#include <Intrepid2_CellGeometry.hpp>
#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedVectorData.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace
{
  using namespace Intrepid2;

//! version of integrate that performs a standard integration for affine meshes; does not take advantage of the tensor product structure at all
//! this version can be used to verify correctness of other versions
template<class Scalar, typename DeviceType>
void integrate_baseline(Data<Scalar,DeviceType> integrals, const TransformedVectorData<Scalar,DeviceType> vectorDataLeft,
                        const TensorData<Scalar,DeviceType> cellMeasures, const TransformedVectorData<Scalar,DeviceType> vectorDataRight)
{
  const int spaceDim       = vectorDataLeft.spaceDim();
  
  // use the CFPD operator() provided by the vector data objects; don't take advantage of tensor product structure at all
  const int numPoints = vectorDataLeft.numPoints();
  
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataLeft.rank()  != 4, std::invalid_argument, "vectorDataLeft must be of shape (C,F,P,D)");
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataRight.rank() != 4, std::invalid_argument, "vectorDataRight must be of shape (C,F,P,D)");
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(cellMeasures.rank()    != 2, std::invalid_argument, "cellMeasures must be of shape (C,P)");
  
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataLeft.spaceDim() != vectorDataRight.spaceDim(), std::invalid_argument, "vectorDataLeft and vectorDataRight must agree on the spatial dimension");
  
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataRight.extent_int(2) != vectorDataLeft.extent_int(2), std::invalid_argument, "vectorData point dimensions must match");
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(cellMeasures.extent_int(1) != vectorDataLeft.extent_int(2), std::invalid_argument,
                                           "cellMeasures point dimension (" + std::to_string(cellMeasures.extent_int(1)) +
                                           ") must match vectorData point dimension (" + std::to_string(vectorDataLeft.extent_int(2)) + ")");
  
//  printFunctor4(vectorDataLeft, std::cout, "vectorDataLeft");
//  printFunctor2(cellMeasures, std::cout, "cellMeasures");
  
  // integral data may have shape (C,F1,F2) or (if the variation type is CONSTANT in the cell dimension) shape (F1,F2)
  const int integralViewRank = integrals.getUnderlyingViewRank();
  
  using ExecutionSpace = typename DeviceType::execution_space;
  auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{integrals.getDataExtent(0),vectorDataLeft.numFields(),vectorDataRight.numFields()});
  Kokkos::parallel_for("fill expanded cell nodes", policy,
  KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &fieldOrdinalLeft, const int &fieldOrdinalRight)
  {
    Scalar integral=0;
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      Scalar pointSum = 0.0;
      for (int d=0; d<spaceDim; d++)
      {
        pointSum += vectorDataLeft(cellDataOrdinal,fieldOrdinalLeft,pointOrdinal,d) * vectorDataRight(cellDataOrdinal,fieldOrdinalRight,pointOrdinal,d);
      }
      integral += pointSum * cellMeasures(cellDataOrdinal,pointOrdinal);
    }
    if (integralViewRank == 3)
    {
      // shape (C,F1,F2)
      auto integralView = integrals.getUnderlyingView3();
      integralView(cellDataOrdinal,fieldOrdinalLeft,fieldOrdinalRight) = integral;
    }
    else
    {
      // shape (F1,F2)
      auto integralView = integrals.getUnderlyingView2();
      integralView(fieldOrdinalLeft,fieldOrdinalRight) = integral;
    }
  });
//  const int numFieldsLeft  = vectorDataLeft.numFields();
//  const int numFieldsRight = vectorDataRight.numFields();
//  int approximateFlopCount = (spaceDim*2 + 2) * numPoints * numFieldsLeft * numFieldsRight * cellDataExtent;
//  printView(integralView, std::cout, "stiffness in " + std::to_string(spaceDim) + "D");
//  std::cout << "\n\nApproximate flop count (baseline): " << approximateFlopCount << std::endl;
}

//! version that uses the classic Intrepid2 paths
template<class Scalar, class PointScalar, int spaceDim, typename DeviceType>
ScalarView<Scalar,DeviceType> performStandardQuadratureHypercube(int meshWidth, int polyOrder, int worksetSize)
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using CellTools = Intrepid2::CellTools<DeviceType>;
  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<ExecutionSpace>;
  
  using namespace std;
  // dimensions of the returned view are (C,F,F)
  auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;

  shards::CellTopology lineTopo = shards::getCellTopologyData< shards::Line<> >();
  shards::CellTopology cellTopo;
  if      (spaceDim == 1) cellTopo = shards::getCellTopologyData< shards::Line<>          >();
  else if (spaceDim == 2) cellTopo = shards::getCellTopologyData< shards::Quadrilateral<> >();
  else if (spaceDim == 3) cellTopo = shards::getCellTopologyData< shards::Hexahedron<>    >();
  
  auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<DeviceType> >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numHypercubes = 1;
  for (int d=0; d<spaceDim; d++)
  {
    numHypercubes *= meshWidth;
  }
  int numCells = numHypercubes;
  
  if (worksetSize > numCells) worksetSize = numCells;
  
  // local stiffness matrix:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  
  auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
  int numPoints = cubature->getNumPoints();
  ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
  ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
  
  cubature->getCubature(cubaturePoints, cubatureWeights);
  
  // Allocate some intermediate containers
  ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
  ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

  ScalarView<Scalar,DeviceType> transformedGradValues("transformed grad values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> transformedWeightedGradValues("transformed weighted grad values", worksetSize, numFields, numPoints, spaceDim);
  
  basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
  basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
  
  CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
  
  const int numNodesPerCell = cellNodes.numNodesPerCell();
  ScalarView<PointScalar,DeviceType> expandedCellNodes("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
  using ExecutionSpace = typename DeviceType::execution_space;
  auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numNodesPerCell});
  Kokkos::parallel_for("fill expanded cell nodes", policy,
  KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinal)
  {
    for (int d=0; d<spaceDim; d++)
    {
      expandedCellNodes(cellOrdinal,nodeOrdinal,d) = cellNodes(cellOrdinal,nodeOrdinal,d);
    }
  });
  
  // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell

  ScalarView<Scalar,DeviceType> cellMeasures("cell measures", worksetSize, numPoints);
  ScalarView<Scalar,DeviceType> jacobianDeterminant("jacobian determinant", worksetSize, numPoints);
  ScalarView<Scalar,DeviceType> jacobian("jacobian", worksetSize, numPoints, spaceDim, spaceDim);
  ScalarView<Scalar,DeviceType> jacobianInverse("jacobian inverse", worksetSize, numPoints, spaceDim, spaceDim);

  int cellOffset = 0;
  while (cellOffset < numCells)
  {
    int startCell         = cellOffset;
    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
    
    std::pair<int,int> cellRange = {startCell, startCell+numCellsInWorkset};
    auto cellWorkset = Kokkos::subview(expandedCellNodes, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    // note that the following will not work if numCellsInWorkset != worksetSize
    // (we would need to take an appropriate subview of jacobian, etc. containers)
    INTREPID2_TEST_FOR_EXCEPTION(numCellsInWorkset != worksetSize, std::invalid_argument, "workset size must evenly divide the number of cells!");
    CellTools::setJacobian(jacobian, cubaturePoints, cellWorkset, cellTopo);
    CellTools::setJacobianInv(jacobianInverse, jacobian);
    CellTools::setJacobianDet(jacobianDeterminant, jacobian);
    
    FunctionSpaceTools::computeCellMeasure(cellMeasures, jacobianDeterminant, cubatureWeights);
    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, jacobianInverse, basisGradValues);
    FunctionSpaceTools::multiplyMeasure(transformedWeightedGradValues, cellMeasures, transformedGradValues);
    
//    printView(transformedGradValues, std::cout, "transformedGradValues");
//    printView(cellMeasures, std::cout, "cellMeasures");
    
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    FunctionSpaceTools::integrate(cellStiffnessSubview, transformedGradValues, transformedWeightedGradValues);
    
    cellOffset += worksetSize;
  }
  return cellStiffness;
}

  template<class Scalar, typename DeviceType>
  void testIntegrateMatchesBaseline(const TransformedVectorData<Scalar,DeviceType> vectorDataLeft,
                                    const TensorData<Scalar,DeviceType> cellMeasures, const TransformedVectorData<Scalar,DeviceType> vectorDataRight,
                                    Teuchos::FancyOStream &out, bool &success)
  {
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
    using IntegrationTools = Intrepid2::IntegrationTools<DeviceType>;
    
    auto integralsBaseline  = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
    auto integralsIntegrate = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
    
    integrate_baseline(integralsBaseline, vectorDataLeft, cellMeasures, vectorDataRight);
    IntegrationTools::integrate(integralsIntegrate, vectorDataLeft, cellMeasures, vectorDataRight);
    
    const int integralsBaselineViewRank = integralsBaseline.getUnderlyingViewRank();
    
//    std::cout << "integralsBaselineView.extent_int(0) = " << integralsBaselineView.extent_int(0) << std::endl;
//    std::cout << "integralsBaselineView.extent_int(1) = " << integralsBaselineView.extent_int(1) << std::endl;
//    std::cout << "integralsBaselineView.extent_int(2) = " << integralsBaselineView.extent_int(2) << std::endl;
//
//    std::cout << "integralsBaselineView.rank() = " << integralsBaselineView.rank() << std::endl;
//
//    std::cout << "integralsBaselineView(0,0)  = " << integralsBaselineView(0,0) << std::endl;
//    std::cout << "integralsIntegrateView(0,0) = " << integralsIntegrateView(0,0) << std::endl;
    
    printFunctor3( integralsBaseline, out, "integralsBaseline");
    printFunctor3(integralsIntegrate, out, "integralsIntegrate");
//    printView(integralsBaselineView, out);
    
    if (integralsBaselineViewRank == 3)
    {
      auto integralsBaselineView  = integralsBaseline.getUnderlyingView3();
      auto integralsIntegrateView = integralsIntegrate.getUnderlyingView3();
      
      testFloatingEquality3(integralsBaselineView, integralsIntegrateView, relTol, absTol, out, success, "baseline integral", "sum factorized integral");
    }
    else
    {
      auto integralsBaselineView  = integralsBaseline.getUnderlyingView2();
      auto integralsIntegrateView = integralsIntegrate.getUnderlyingView2();
      
      testFloatingEquality2(integralsBaselineView, integralsIntegrateView, relTol, absTol, out, success, "baseline integral", "sum factorized integral");
    }
  }

  template<class Scalar, class PointScalar, int spaceDim>
  void testQuadratureHypercube(bool useAffinePath, int meshWidth, int polyOrder, int worksetSize, Teuchos::FancyOStream &out, bool &success)
  {
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
    using namespace std;
    using DeviceType = DefaultTestDeviceType;
    
    using IntegrationTools   = Intrepid2::IntegrationTools<DeviceType>;
    // dimensions of the returned view are (C,F,F)
    auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
    
    auto lineBasis = Intrepid2::getLineBasis< Intrepid2::NodalBasisFamily<DeviceType> >(fs, polyOrder);
    
    int numFields_1D = lineBasis->getCardinality();
    
    int numFields = 1;
    int numHypercubes = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numHypercubes *= meshWidth;
      numFields     *= numFields_1D;
    }
    int numCells = numHypercubes;
    
    if (worksetSize > numCells) worksetSize = numCells;
    
    // local stiffness matrix:
    ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
    
    shards::CellTopology lineTopo = shards::getCellTopologyData< shards::Line<> >();
    shards::CellTopology cellTopo;
    if      (spaceDim == 1) cellTopo = shards::getCellTopologyData< shards::Line<>          >();
    else if (spaceDim == 2) cellTopo = shards::getCellTopologyData< shards::Quadrilateral<> >();
    else if (spaceDim == 3) cellTopo = shards::getCellTopologyData< shards::Hexahedron<>    >();
    
    auto lineCubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(lineTopo,polyOrder*2);
    int numPoints_1D = lineCubature->getNumPoints();
    ScalarView<PointScalar,DeviceType> lineCubaturePoints("line cubature points",numPoints_1D,1);
    ScalarView<double,DeviceType> lineCubatureWeights("line cubature weights", numPoints_1D);
    
    lineCubature->getCubature(lineCubaturePoints, lineCubatureWeights);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,DeviceType> lineBasisValues    ("line basis values",      numFields_1D, numPoints_1D   );
    ScalarView<Scalar,DeviceType> lineBasisGradValues("line basis grad values", numFields_1D, numPoints_1D, 1);
    
    // for now, we use 1D values to build up the 2D or 3D gradients
    // eventually, TensorBasis should offer a getValues() variant that returns tensor basis data
    lineBasis->getValues(lineBasisValues,     lineCubaturePoints, Intrepid2::OPERATOR_VALUE );
    lineBasis->getValues(lineBasisGradValues, lineCubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    // drop the trivial space dimension in line gradient values:
    Kokkos::resize(lineBasisGradValues, numFields_1D, numPoints_1D);
      
    Kokkos::Array<TensorData<Scalar,DeviceType>, spaceDim> vectorComponents;
    
    for (int d=0; d<spaceDim; d++)
    {
      Kokkos::Array<Data<Scalar,DeviceType>, spaceDim> gradComponent_d;
      for (int d2=0; d2<spaceDim; d2++)
      {
        if (d2 == d) gradComponent_d[d2] = Data<Scalar,DeviceType>(lineBasisGradValues);
        else         gradComponent_d[d2] = Data<Scalar,DeviceType>(lineBasisValues);
      }
      vectorComponents[d] = TensorData<Scalar,DeviceType>(gradComponent_d);
    }
    VectorData<Scalar,DeviceType> gradientValues(vectorComponents, false); // false: not axis-aligned
    
    CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
    
    if (useAffinePath)
    {
      out << "Testing path with affine, grid-aligned CellGeometry.\n";
    }
    else
    {
      // make a "generic" copy of cellNodes, one that uses the (C,N), (N,D) node specification.  This will not know that the geometry is affine, grid-aligned, or uniform.
      const bool copyAffineness = false; // want to go through the non-affine geometry path
      CellGeometry<PointScalar,spaceDim,DeviceType> nodalCellGeometry = getNodalCellGeometry(cellNodes, copyAffineness);
      
      cellNodes = nodalCellGeometry;
      out << "Testing non-affine path.\n";
    }
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
    auto tensorCubatureWeights = cubature->allocateCubatureWeights();
    auto tensorCubaturePoints  = cubature->allocateCubaturePoints();
    
    cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
    
    // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
    
    int pointsPerCell = 1;
    for (int d=0; d<spaceDim; d++)
    {
      pointsPerCell *= numPoints_1D;
    }
    
    Data<PointScalar,DeviceType> jacobian = cellNodes.allocateJacobianData(tensorCubaturePoints);
    Data<PointScalar,DeviceType> jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
    Data<PointScalar,DeviceType> jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
    
    auto refData = cellNodes.getJacobianRefData(tensorCubaturePoints);
    cellNodes.setJacobian(jacobian, tensorCubaturePoints, refData);
    CellTools<DeviceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<DeviceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    using FunctionSpaceToolsDT = ::Intrepid2::FunctionSpaceTools<DeviceType>; // TODO: once FunctionSpaceTools has proper DeviceType support, change the earlier FunctionSpaceTools typedef, and use it here…
    auto transformedGradientValues = FunctionSpaceToolsDT::getHGRADtransformGRAD(jacobianInv, gradientValues);
    auto standardIntegrals = performStandardQuadratureHypercube<Scalar,PointScalar,spaceDim,DeviceType>(meshWidth, polyOrder, worksetSize);
    
    TensorData<PointScalar,DeviceType> cellMeasures = cellNodes.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
    if (!cellNodes.affine())
    {
      // if cellNodes is not (known to be) affine, then cellMeasures should not have a separate first component (indicating the cell dimension is separated, thanks to point-invariant cell Jacobian determinant)
      TEST_EQUALITY(cellMeasures.separateFirstComponent(), false);
    }
    cellNodes.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    
    auto integralData = IntegrationTools::allocateIntegralData(transformedGradientValues, cellMeasures, transformedGradientValues);
    
    IntegrationTools::integrate(integralData, transformedGradientValues, cellMeasures, transformedGradientValues);
    
    auto integralDataBaseline = IntegrationTools::allocateIntegralData(transformedGradientValues, cellMeasures, transformedGradientValues);
    
    integrate_baseline(integralDataBaseline, transformedGradientValues, cellMeasures, transformedGradientValues);
    
    out << "Comparing new integration path path with baseline integration…\n";
    testIntegrateMatchesBaseline(transformedGradientValues, cellMeasures, transformedGradientValues, out, success);
    
    out << "Comparing baseline to standard Intrepid2 integration…\n";
    testFloatingEquality3(standardIntegrals, integralDataBaseline, relTol, absTol, out, success, "standard Intrepid2 integral", "reduced data integral - baseline");
    
    out << "Comparing new integration path with standard Intrepid2 integration…\n";
    testFloatingEquality3(standardIntegrals, integralData, relTol, absTol, out, success, "standard Intrepid2 integral", "reduced data integral");
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_1D_p1_AffinePath
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_1D_p1_AffinePath )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 1;
    const int polyOrder = 1;
    const int spaceDim = 1;
    const int worksetSize = meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_1D_p1_GeneralPath
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_1D_p1_GeneralPath )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 4;
    const int polyOrder = 1;
    const int spaceDim = 1;
    const int worksetSize = meshWidth;
    const bool affinePath = false;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_1D_p4_GeneralPath
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_1D_p4_GeneralPath )
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 4;
  const int polyOrder = 4;
  const int spaceDim = 1;
  const int worksetSize = meshWidth;
  const bool affinePath = false;
  
  testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
}

// #pragma mark StructuredIntegration: QuadratureUniformMesh_1D_p2
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_1D_p2 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 1;
    const int polyOrder = 2;
    const int spaceDim = 1;
    const int worksetSize = meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_2D_p1
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_2D_p1 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 3;
    const int polyOrder = 1;
    const int spaceDim = 2;
    const int worksetSize = meshWidth * meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_2D_p1_GeneralPath
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_2D_p1_GeneralPath )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 2;
    const int polyOrder = 1;
    const int spaceDim = 2;
    const int worksetSize = meshWidth * meshWidth;
    const bool affinePath = false;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_2D_p2_GeneralPath
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_2D_p2_GeneralPath )
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int polyOrder = 2;
  const int spaceDim = 2;
  const int worksetSize = meshWidth * meshWidth;
  const bool affinePath = false;
  
  testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
}

// #pragma mark StructuredIntegration: QuadratureUniformMesh_2D_p2
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_2D_p2 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 3;
    const int polyOrder = 2;
    const int spaceDim = 2;
    const int worksetSize = meshWidth * meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p1
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p1 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 2;
    const int polyOrder = 1;
    const int spaceDim = 3;
    const int worksetSize = meshWidth * meshWidth * meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p1_GeneralPath
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p1_GeneralPath )
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int polyOrder = 1;
  const int spaceDim = 3;
  const int worksetSize = meshWidth * meshWidth * meshWidth;
  const bool affinePath = false;
  
  testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
}

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p2
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p2 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 2;
    const int polyOrder = 2;
    const int spaceDim = 3;
    const int worksetSize = meshWidth * meshWidth * meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p2_GeneralPath
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p2_GeneralPath )
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int polyOrder = 2;
  const int spaceDim = 3;
  const int worksetSize = meshWidth * meshWidth * meshWidth;
  const bool affinePath = false;
  
  testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
}

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p3
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p3 )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 2;
    const int polyOrder = 3;
    const int spaceDim = 3;
    const int worksetSize = meshWidth * meshWidth * meshWidth;
    const bool affinePath = true;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureUniformMesh_3D_p3_GeneralPath
  TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureUniformMesh_3D_p3_GeneralPath )
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 1;
    const int polyOrder = 3;
    const int spaceDim = 3;
    const int worksetSize = meshWidth * meshWidth * meshWidth;
    const bool affinePath = false;
    
    testQuadratureHypercube<DataScalar,PointScalar,spaceDim>(affinePath,meshWidth, polyOrder, worksetSize, out, success);
  }

// #pragma mark StructuredIntegration: QuadratureSynthetic_AxisAlignedPath_Case1
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_AxisAlignedPath_Case1 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numPoints = numComponentPoints * numComponentPoints;
  const int numFields = 1;
  
  Data<DataScalar,DeviceType> unitData(1.0, Kokkos::Array<int,2>{numFields,numComponentPoints});
  TensorData<DataScalar,DeviceType> unitTensorData(std::vector< Data<DataScalar,DeviceType> >{unitData,unitData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
  
  const int blockPlusDiagonalLastNonDiagonal = -1; // no non-diagonals; -1 is the default value, but I prefer to be explicit
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType, blockPlusDiagonalLastNonDiagonal);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {unitTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),unitTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> unitVectorData(vectorComponents);
  
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresCellComponent(1.0, Kokkos::Array<int,1>{numCells});
  Data<DataScalar,DeviceType> constantCellMeasuresPointComponent(1.0, Kokkos::Array<int,1>{numComponentPoints});
  std::vector< Data<DataScalar,DeviceType> > cellMeasuresComponents { constantCellMeasuresCellComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent};
  
  const bool separateFirstComponent = true; // so that constantCellMeasures advertises shape (C,P), instead of folding cell dimension into the point tensor product...
  TensorData<DataScalar,DeviceType> constantCellMeasures(cellMeasuresComponents, separateFirstComponent);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case1
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case1 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numPoints = 1;
  const int numFields = 1;
  
  Data<DataScalar,DeviceType> unitData(1.0, Kokkos::Array<int,2>{numFields,numPoints});
  TensorData<DataScalar,DeviceType> unitTensorData(std::vector< Data<DataScalar,DeviceType> >{unitData,unitData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {unitTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),unitTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> unitVectorData(vectorComponents);
  
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case2
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case2 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  const int numComponentFields = 1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numComponentFields);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numComponentFields);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 2.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numComponentFields,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,CONSTANT};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType> tensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim> family  {TensorData<DataScalar,DeviceType>(), tensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {family};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TransformedVectorData<DataScalar,DeviceType> transformedVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
//  printFunctor2(  fieldComponentData1, std::cout, "fieldComponentData1"); // (F,P)
//  printFunctor2(  fieldComponentData2, std::cout, "fieldComponentData2"); // (F,P)
//  printFunctor2(    nonzeroTensorData, std::cout, "nonzeroTensorData");   // (F,P)
//  printFunctor2(       zeroTensorData, std::cout, "zeroTensorData");      // (F,P)
//
//  printFunctor3(          vectorData, std::cout, "vectorData");             // (F,P,D)
//  printFunctor2( constantCellMeasures, std::cout, "constantCellMeasures");  // (C,P)
//  printFunctor4(transformedVectorData, std::cout, "transformedVectorData"); // (C,F,P,D)
  
  const int numFields = numComponentFields * numComponentFields;
  TEST_EQUALITY(numFields, vectorData.extent_int(0)); // (F,P,D)
  TEST_EQUALITY(numFields, transformedVectorData.extent_int(1)); // (C,F,P,D)
  
  testIntegrateMatchesBaseline(transformedVectorData, constantCellMeasures, transformedVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case3
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case3 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 2;
  const int numFields = 2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields,numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1.0;
  fieldComponentDataViewHost1(0,1) = 2.0;
  fieldComponentDataViewHost1(1,0) = 3.0;
  fieldComponentDataViewHost1(1,1) = 4.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields,numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1.0/1.0;
  fieldComponentDataViewHost2(0,1) = 1.0/2.0;
  fieldComponentDataViewHost2(1,0) = 1.0/3.0;
  fieldComponentDataViewHost2(1,1) = 1.0/4.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType> tensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {tensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),tensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case4
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case4 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 1;
  const int numFields1 = 2;
  const int numFields2 = 1;
  const int numFieldsPerFamily = numFields1 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType> nonzeroTensorData(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamily  {nonzeroTensorData,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamily {TensorData<DataScalar,DeviceType>(),nonzeroTensorData};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponents {firstFamily, secondFamily};
  
  VectorData<DataScalar,DeviceType> vectorData(vectorComponents);
  
  TEST_EQUALITY(numFieldsPerFamily, vectorData.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamily, vectorData.numFieldsInFamily(1));
  
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorData, constantCellMeasures, transformedUnitVectorData, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case5
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case5 )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 2;
  
  const int numCells  = 2;
  const int numComponentPoints = 1;
  const int numFields1 = 2;
  const int numFields2 = 3;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2;
  const int numFieldsPerFamilyRight = numFields2 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
    
  TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData2});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 2;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyLeft {TensorData<DataScalar,DeviceType>(),tensorDataLeft};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>()};
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyRight {TensorData<DataScalar,DeviceType>(),tensorDataRight};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight};
  
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
  TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_AxisAlignedPath_Case6_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_AxisAlignedPath_Case6_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 3;
  const int numFields3 = 2;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
//  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 1.0;
  fieldComponentDataViewHost3(1) = 0.25;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
   const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
   
   auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim);
   Kokkos::deep_copy(identityMatrixView, 1.0);
   
   Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
   Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, BLOCK_PLUS_DIAGONAL, BLOCK_PLUS_DIAGONAL};
   
   const int blockPlusDiagonalLastNonDiagonal = -1; // no non-diagonals; -1 is the default value, but I prefer to be explicit
   Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType, blockPlusDiagonalLastNonDiagonal);
  
  // confirm that the matrix is diagonal (this is required to follow the axis-aligned path):
  TEST_EQUALITY(true, explicitIdentityMatrix.isDiagonal());

  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,tensorDataLeft,tensorDataLeft};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,tensorDataRight,tensorDataRight};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType>  constantCellMeasuresCellComponent(1.0, Kokkos::Array<int,1>{numCells});
   Data<DataScalar,DeviceType> constantCellMeasuresPointComponent(1.0, Kokkos::Array<int,1>{numComponentPoints});
   std::vector< Data<DataScalar,DeviceType> > cellMeasuresComponents { constantCellMeasuresCellComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent, constantCellMeasuresPointComponent};
   
   const bool separateFirstComponent = true; // so that constantCellMeasures advertises shape (C,P), instead of folding cell dimension into the point tensor product...
   TensorData<DataScalar,DeviceType> constantCellMeasures(cellMeasuresComponents, separateFirstComponent);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case6_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case6_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 3;
  const int numFields3 = 2;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
//  fieldComponentDataViewHost1(1) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  fieldComponentDataViewHost2(1) = 1.0/3.0;
  fieldComponentDataViewHost2(2) = 1.0/4.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 0.5;
  fieldComponentDataViewHost3(1) = 0.25;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
   auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
   auto identityMatrixViewHost = getHostCopy(identityMatrixView);
   
   for (int d1=0; d1<spaceDim; d1++)
   {
     for (int d2=0; d2<spaceDim; d2++)
     {
       identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
     }
   }
   Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
   
   const int numPoints = numComponentPoints * numComponentPoints;
   Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
   Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
   
   Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
   
//   const int numFamilies = 3;
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyLeft {TensorData<DataScalar,DeviceType>(),tensorDataLeft,TensorData<DataScalar,DeviceType>()};
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyLeft  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataLeft};
//   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft, thirdFamilyLeft};
//
//   VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
//   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
//   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));
//
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyRight {TensorData<DataScalar,DeviceType>(),tensorDataRight,TensorData<DataScalar,DeviceType>()};
//   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyRight  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataRight};
//   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight, thirdFamilyRight};
//  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
//  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
//  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case7_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case7_3D )
{
  // super-simple case for 3D: symmetric data, 1 point, 1 field.
  
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFields3 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 1.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
   auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
   auto identityMatrixViewHost = getHostCopy(identityMatrixView);
   
   for (int d1=0; d1<spaceDim; d1++)
   {
     for (int d2=0; d2<spaceDim; d2++)
     {
       identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
     }
   }
   Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
   
   const int numPoints = numComponentPoints * numComponentPoints;
   Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
   Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
   
   Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
     
   const int numFamilies = 1;
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
   VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case8_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case8_3D )
{
  // like case 7, but multi-family
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 1;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFields3 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields3;
  const int numFieldsPerFamilyRight = numFields1 * numFields3 * numFields2;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("field component data 1", numFields1);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0) = 1.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("field component data 2", numFields2);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0) = 1.0/2.0;
  
  auto fieldComponentDataView3 = getFixedRankView<DataScalar>("field component data 3", numFields3);
  auto fieldComponentDataViewHost3 = Kokkos::create_mirror_view(fieldComponentDataView3);
  fieldComponentDataViewHost3(0) = 0.5;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  Kokkos::deep_copy(fieldComponentDataView3, fieldComponentDataViewHost3);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields3;
  Data<DataScalar,DeviceType> fieldComponentData3(fieldComponentDataView3,fieldComponentExtents,fieldComponentVariationTypes);
  
   TensorData<DataScalar,DeviceType>  tensorDataLeft(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData3});
   TensorData<DataScalar,DeviceType> tensorDataRight(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData3,fieldComponentData2});
   
   auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
   auto identityMatrixViewHost = getHostCopy(identityMatrixView);
   
   for (int d1=0; d1<spaceDim; d1++)
   {
     for (int d2=0; d2<spaceDim; d2++)
     {
       identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
     }
   }
   Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
   
   const int numPoints = numComponentPoints * numComponentPoints;
   Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
   Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
   
   Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
   
   const int numFamilies = 3;
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyLeft {TensorData<DataScalar,DeviceType>(),tensorDataLeft,TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyLeft  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataLeft};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, secondFamilyLeft, thirdFamilyLeft};

   VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
   TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(1));

   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > secondFamilyRight {TensorData<DataScalar,DeviceType>(),tensorDataRight,TensorData<DataScalar,DeviceType>()};
   Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > thirdFamilyRight  {TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>(),tensorDataRight};
   Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, secondFamilyRight, thirdFamilyRight};
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(1));
  
   TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
   Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
   TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
   
   testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case9_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case9_3D )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 3;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields1;
  const int numFieldsPerFamilyRight = numFields1 * numFields2 * numFields1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("imitation VALUE data", numFields1, numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1.0;
  fieldComponentDataViewHost1(0,1) = 2.0;
  fieldComponentDataViewHost1(0,2) = 3.0;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("imitation GRAD data", numFields2, numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1./1.0;
  fieldComponentDataViewHost2(0,1) = 1./2.0;
  fieldComponentDataViewHost2(0,2) = 1./3.0;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType>  tensorDataLeft1(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData1,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft2(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft3(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight1 = tensorDataLeft1;
  TensorData<DataScalar,DeviceType> tensorDataRight2 = tensorDataLeft2;
  TensorData<DataScalar,DeviceType> tensorDataRight3 = tensorDataLeft3;
   
  const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", numPoints, spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    for (int d1=0; d1<spaceDim; d1++)
    {
      for (int d2=0; d2<spaceDim; d2++)
      {
        identityMatrixViewHost(pointOrdinal,d1,d2) = (d1 == d2) ? 1.0 : 0.0;
      }
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, GENERAL, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft1,tensorDataLeft2,tensorDataLeft3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};

  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));

  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight1,tensorDataRight2,tensorDataRight3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  
  TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  Data<DataScalar,DeviceType> constantCellMeasuresData(1.0, Kokkos::Array<int,2>{numCells,numPoints});
  TensorData<DataScalar,DeviceType> constantCellMeasures(constantCellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, constantCellMeasures, transformedUnitVectorDataRight, out, success);
}

// #pragma mark StructuredIntegration: QuadratureSynthetic_GeneralPath_Case10_3D
TEUCHOS_UNIT_TEST( StructuredIntegration, QuadratureSynthetic_GeneralPath_Case10_3D )
{
  // test with variable quadrature weights
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  
  const int spaceDim = 3;
  
  const int numCells  = 1;
  const int numComponentPoints = 2;
  
  const int numFields1 = 1;
  const int numFields2 = 1;
  const int numFieldsPerFamilyLeft  = numFields1 * numFields2 * numFields1;
  const int numFieldsPerFamilyRight = numFields1 * numFields2 * numFields1;
  
  auto fieldComponentDataView1 = getFixedRankView<DataScalar>("imitation VALUE data", numFields1, numComponentPoints);
  auto fieldComponentDataViewHost1 = Kokkos::create_mirror_view(fieldComponentDataView1);
  fieldComponentDataViewHost1(0,0) = 1;
  fieldComponentDataViewHost1(0,1) = 1;
  
  auto fieldComponentDataView2 = getFixedRankView<DataScalar>("imitation GRAD data", numFields2, numComponentPoints);
  auto fieldComponentDataViewHost2 = Kokkos::create_mirror_view(fieldComponentDataView2);
  fieldComponentDataViewHost2(0,0) = 1;
  fieldComponentDataViewHost2(0,1) = 1;
  
  Kokkos::deep_copy(fieldComponentDataView1, fieldComponentDataViewHost1);
  Kokkos::deep_copy(fieldComponentDataView2, fieldComponentDataViewHost2);
  
  const int fieldComponentDataRank = 2;
  Kokkos::Array<int,fieldComponentDataRank> fieldComponentExtents {numFields1,numComponentPoints};
  Kokkos::Array<DataVariationType,fieldComponentDataRank> fieldComponentVariationTypes {GENERAL,GENERAL};
  Data<DataScalar,DeviceType> fieldComponentData1(fieldComponentDataView1,fieldComponentExtents,fieldComponentVariationTypes);
  fieldComponentExtents[0] = numFields2;
  Data<DataScalar,DeviceType> fieldComponentData2(fieldComponentDataView2,fieldComponentExtents,fieldComponentVariationTypes);
  
  TensorData<DataScalar,DeviceType>  tensorDataLeft1(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData2,fieldComponentData1,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft2(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData2,fieldComponentData1});
  TensorData<DataScalar,DeviceType>  tensorDataLeft3(std::vector< Data<DataScalar,DeviceType> >{fieldComponentData1,fieldComponentData1,fieldComponentData2});
  TensorData<DataScalar,DeviceType> tensorDataRight1 = tensorDataLeft1;
  TensorData<DataScalar,DeviceType> tensorDataRight2 = tensorDataLeft2;
  TensorData<DataScalar,DeviceType> tensorDataRight3 = tensorDataLeft3;
   
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  auto identityMatrixViewHost = getHostCopy(identityMatrixView);
  
  for (int d1=0; d1<spaceDim; d1++)
  {
    for (int d2=0; d2<spaceDim; d2++)
    {
      identityMatrixViewHost(d1,d2) = (d1 == d2) ? 1.0 : 0.0;
    }
  }
  Kokkos::deep_copy(identityMatrixView, identityMatrixViewHost);
  
  const int numPoints = numComponentPoints * numComponentPoints * numComponentPoints;
  Kokkos::Array<int,4> transformationExtents {numCells, numPoints, spaceDim, spaceDim};
  Kokkos::Array<DataVariationType,4> transformationVariationType {CONSTANT, CONSTANT, GENERAL, GENERAL};
  
  Data<DataScalar,DeviceType> explicitIdentityMatrix(identityMatrixView, transformationExtents, transformationVariationType);
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft1,tensorDataLeft2,tensorDataLeft3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};

  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));

  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight1,tensorDataRight2,tensorDataRight3};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
  VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
  TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
  
  TransformedVectorData<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedVectorData<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
  auto cellMeasures = getFixedRankView<DataScalar>("cellMeasures", numCells, numPoints);
  
  auto cellMeasuresHost = getHostCopy(cellMeasures);
  
  for (int cellOrdinal=0; cellOrdinal<numCells; cellOrdinal++)
  {
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      cellMeasuresHost(cellOrdinal,pointOrdinal) = (pointOrdinal == 0) ? 1 : 0; //1.0 / (pointOrdinal + 1.0);
    }
  }
  Kokkos::deep_copy(cellMeasures, cellMeasuresHost);
  
  Data<DataScalar,DeviceType> cellMeasuresData(cellMeasures, Kokkos::Array<int,2>{numCells,numPoints}, Kokkos::Array<DataVariationType,2>{GENERAL,GENERAL});
  TensorData<DataScalar,DeviceType> cellMeasuresTensorData(cellMeasuresData);
  
  testIntegrateMatchesBaseline(transformedUnitVectorDataLeft, cellMeasuresTensorData, transformedUnitVectorDataRight, out, success);
}

} // anonymous namespace
