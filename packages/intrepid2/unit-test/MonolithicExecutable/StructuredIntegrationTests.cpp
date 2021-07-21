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

#include "GRADGRADStandardAssembly.hpp"
#include "GRADGRADStructuredAssembly.hpp"
#include "H1StandardAssembly.hpp"
#include "H1StructuredAssembly.hpp"
#include "HDIVStandardAssembly.hpp"
#include "HDIVStructuredAssembly.hpp"
#include "HCURLStandardAssembly.hpp"
#include "HCURLStructuredAssembly.hpp"
#include "HVOLStandardAssembly.hpp"
#include "HVOLStructuredAssembly.hpp"

namespace
{
  using namespace Intrepid2;

  enum FormulationChoice
  {
    Poisson, // (grad, grad)
    Hgrad,   // (grad, grad) + (value, value)
    Hdiv,    // (div, div)   + (value, value)
    Hcurl,   // (curl, curl) + (value, value)
    L2       // (value, value)
  };

  enum AlgorithmChoice
  {
    Standard,
    AffineNonTensor,
    NonAffineTensor,
    AffineTensor,
    DiagonalJacobian,
    Uniform
  };

  enum BasisFamilyChoice
  {
    Nodal,
    Hierarchical,
    Serendipity
  };

  std::string to_string(AlgorithmChoice choice)
  {
    switch (choice) {
      case Standard:         return "Standard";
      case AffineNonTensor:  return "AffineNonTensor";
      case NonAffineTensor:  return "NonAffineTensor";
      case AffineTensor:     return "AffineTensor";
      case DiagonalJacobian: return "DiagonalJacobian";
      case Uniform:          return "Uniform";
      
      default:               return "Unknown AlgorithmChoice";
    }
  }

  std::string to_string(FormulationChoice choice)
  {
    switch (choice) {
      case Poisson: return "Poisson";
      case Hgrad:   return "Hgrad";
      case Hdiv:    return "Hdiv";
      case Hcurl:   return "Hcurl";
      case L2:      return "L2";
      
      default:      return "Unknown FormulationChoice";
    }
  }

  using namespace Intrepid2;

  template< typename PointScalar, int spaceDim, typename DeviceType >
  inline
  CellGeometry<PointScalar, spaceDim, DeviceType> getMesh(AlgorithmChoice algorithmChoice, const Kokkos::Array<int,spaceDim> &gridCellCounts)
  {
    Kokkos::Array<PointScalar,spaceDim> domainExtents;
    for (int d=0; d<spaceDim; d++)
    {
      domainExtents[d]  = 1.0;
    }
    auto uniformTensorGeometry = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(domainExtents, gridCellCounts);
    
    switch (algorithmChoice)
    {
      case Standard:
      case NonAffineTensor:
      {
        // Standard and non-affine tensor use the same geometry; the difference is how this is used in assembly
        const bool copyAffineness = false;
        auto genericGeometry = getNodalCellGeometry(uniformTensorGeometry, copyAffineness);
        return genericGeometry;
      }
      case Uniform:
        return uniformTensorGeometry;
      case AffineNonTensor:
      case AffineTensor:
      {
        const bool copyAffineness = true;
        auto affineNonTensorGeometry = getNodalCellGeometry(uniformTensorGeometry, copyAffineness);
        return affineNonTensorGeometry;
      }
      case DiagonalJacobian:
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "DiagonalJacobian case not yet implemented");
      }
    }
    return uniformTensorGeometry; // this line should be unreachable; included to avoid compiler warnings from nvcc
  }

  template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
  Intrepid2::ScalarView<Scalar,DeviceType> performStandardQuadrature(FormulationChoice formulation,
                                          Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                          double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
  {
    switch (formulation)
    {
      case Poisson:
        return performStandardQuadratureGRADGRAD<Scalar,BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hgrad:
        return performStandardQuadratureH1<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hdiv:
        return performStandardQuadratureHDIV<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hcurl:
        return performStandardQuadratureHCURL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case L2:
        return performStandardQuadratureHVOL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported formulation");
    }
  }

  template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
  Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadrature(FormulationChoice formulation,
                                            Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                            double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
  {
    switch (formulation)
    {
      case Poisson:
        return performStructuredQuadratureGRADGRAD<Scalar,BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hgrad:
        return performStructuredQuadratureH1<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hdiv:
        return performStructuredQuadratureHDIV<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case Hcurl:
        return performStructuredQuadratureHCURL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      case L2:
        return performStructuredQuadratureHVOL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported formulation");
    }
  }

//! version of integrate that performs a standard integration for affine meshes; does not take advantage of the tensor product structure at all
//! this version can be used to verify correctness of other versions
template<class Scalar, typename DeviceType>
void integrate_baseline(Data<Scalar,DeviceType> integrals, const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                        const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight)
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
                                           std::string("cellMeasures point dimension (" + std::to_string(cellMeasures.extent_int(1)) +
                                                       ") must match vectorData point dimension (" + std::to_string(vectorDataLeft.extent_int(2)) + ")").c_str());
  
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

////! version that uses the classic Intrepid2 paths
//template<class Scalar, class PointScalar, int spaceDim, typename DeviceType>
//ScalarView<Scalar,DeviceType> performStandardQuadratureHypercube(int meshWidth, int polyOrder, int worksetSize)
//{
//  using ExecutionSpace = typename DeviceType::execution_space;
//  using CellTools = Intrepid2::CellTools<DeviceType>;
//  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<DeviceType>;
//
//  using namespace std;
//  // dimensions of the returned view are (C,F,F)
//  auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
//
//  shards::CellTopology lineTopo = shards::getCellTopologyData< shards::Line<> >();
//  shards::CellTopology cellTopo;
//  if      (spaceDim == 1) cellTopo = shards::getCellTopologyData< shards::Line<>          >();
//  else if (spaceDim == 2) cellTopo = shards::getCellTopologyData< shards::Quadrilateral<> >();
//  else if (spaceDim == 3) cellTopo = shards::getCellTopologyData< shards::Hexahedron<>    >();
//
//  auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<DeviceType> >(cellTopo, fs, polyOrder);
//
//  int numFields = basis->getCardinality();
//  int numHypercubes = 1;
//  for (int d=0; d<spaceDim; d++)
//  {
//    numHypercubes *= meshWidth;
//  }
//  int numCells = numHypercubes;
//
//  if (worksetSize > numCells) worksetSize = numCells;
//
//  // local stiffness matrix:
//  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
//
//  auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
//  int numPoints = cubature->getNumPoints();
//  ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
//  ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
//
//  cubature->getCubature(cubaturePoints, cubatureWeights);
//
//  // Allocate some intermediate containers
//  ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
//  ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);
//
//  ScalarView<Scalar,DeviceType> transformedGradValues("transformed grad values", worksetSize, numFields, numPoints, spaceDim);
//  ScalarView<Scalar,DeviceType> transformedWeightedGradValues("transformed weighted grad values", worksetSize, numFields, numPoints, spaceDim);
//
//  basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
//  basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
//
//  CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
//
//  const int numNodesPerCell = cellNodes.numNodesPerCell();
//  ScalarView<PointScalar,DeviceType> expandedCellNodes("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
//  using ExecutionSpace = typename DeviceType::execution_space;
//  auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numNodesPerCell});
//  Kokkos::parallel_for("fill expanded cell nodes", policy,
//  KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinal)
//  {
//    for (int d=0; d<spaceDim; d++)
//    {
//      expandedCellNodes(cellOrdinal,nodeOrdinal,d) = cellNodes(cellOrdinal,nodeOrdinal,d);
//    }
//  });
//
//  // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
//
//  ScalarView<Scalar,DeviceType> cellMeasures("cell measures", worksetSize, numPoints);
//  ScalarView<Scalar,DeviceType> jacobianDeterminant("jacobian determinant", worksetSize, numPoints);
//  ScalarView<Scalar,DeviceType> jacobian("jacobian", worksetSize, numPoints, spaceDim, spaceDim);
//  ScalarView<Scalar,DeviceType> jacobianInverse("jacobian inverse", worksetSize, numPoints, spaceDim, spaceDim);
//
//  int cellOffset = 0;
//  while (cellOffset < numCells)
//  {
//    int startCell         = cellOffset;
//    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
//
//    std::pair<int,int> cellRange = {startCell, startCell+numCellsInWorkset};
//    auto cellWorkset = Kokkos::subview(expandedCellNodes, cellRange, Kokkos::ALL(), Kokkos::ALL());
//
//    // note that the following will not work if numCellsInWorkset != worksetSize
//    // (we would need to take an appropriate subview of jacobian, etc. containers)
//    INTREPID2_TEST_FOR_EXCEPTION(numCellsInWorkset != worksetSize, std::invalid_argument, "workset size must evenly divide the number of cells!");
//    CellTools::setJacobian(jacobian, cubaturePoints, cellWorkset, cellTopo);
//    CellTools::setJacobianInv(jacobianInverse, jacobian);
//    CellTools::setJacobianDet(jacobianDeterminant, jacobian);
//
//    FunctionSpaceTools::computeCellMeasure(cellMeasures, jacobianDeterminant, cubatureWeights);
//    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, jacobianInverse, basisGradValues);
//    FunctionSpaceTools::multiplyMeasure(transformedWeightedGradValues, cellMeasures, transformedGradValues);
//
////    printView(transformedGradValues, std::cout, "transformedGradValues");
////    printView(cellMeasures, std::cout, "cellMeasures");
//
//    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
//
//    FunctionSpaceTools::integrate(cellStiffnessSubview, transformedGradValues, transformedWeightedGradValues);
//
//    cellOffset += worksetSize;
//  }
//  return cellStiffness;
//}

  template<class Scalar, typename DeviceType>
  void testIntegrateMatchesBaseline(const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                                    const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight,
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

  template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
  void testQuadratureHypercube(int meshWidth, int polyOrder, int worksetSize,
                               const FormulationChoice &formulation, const AlgorithmChoice &algorithm,
                               const double &relTol, const double &absTol,
                               Teuchos::FancyOStream &out, bool &success)
  {
    using namespace std;
    
    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++)
    {
      gridCellCounts[d] = meshWidth;
    }
    
    auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(algorithm, gridCellCounts);
    double flopCountIntegration = 0, flopCountJacobian = 0;
    auto standardIntegrals = performStandardQuadrature<Scalar, BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian);
    
    auto structuredIntegrals = performStructuredQuadrature<Scalar, BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian);
    
    out << "Comparing standard Intrepid2 integration to new integration path…\n";
    testFloatingEquality3(standardIntegrals, structuredIntegrals, relTol, absTol, out, success, "standard Intrepid2 integral", "reduced data integral - baseline");
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

    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = NonAffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const AlgorithmChoice algorithm = NonAffineTensor;
  const FormulationChoice formulation = Poisson;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                      relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = NonAffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const AlgorithmChoice algorithm = NonAffineTensor;
  const FormulationChoice formulation = Poisson;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                      relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const AlgorithmChoice algorithm = NonAffineTensor;
  const FormulationChoice formulation = Poisson;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                      relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const AlgorithmChoice algorithm = NonAffineTensor;
  const FormulationChoice formulation = Poisson;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                      relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
    
    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = NonAffineTensor;
    const FormulationChoice formulation = Poisson;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
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
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,unitVectorData);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType> transformedVectorData(explicitIdentityMatrix,vectorData);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorData(explicitIdentityMatrix,vectorData);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
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
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
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
  
  const int numFamilies = 1;
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyLeft  {tensorDataLeft,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};
  
  VectorData<DataScalar,DeviceType> vectorDataLeft(vectorComponentsLeft);
  TEST_EQUALITY(numFieldsPerFamilyLeft, vectorDataLeft.numFieldsInFamily(0));
  
  Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim > firstFamilyRight  {tensorDataRight,TensorData<DataScalar,DeviceType>(),TensorData<DataScalar,DeviceType>()};
  Kokkos::Array< Kokkos::Array<TensorData<DataScalar,DeviceType>, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight};
   
   VectorData<DataScalar,DeviceType> vectorDataRight(vectorComponentsRight);
   TEST_EQUALITY(numFieldsPerFamilyRight, vectorDataRight.numFieldsInFamily(0));
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
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
   
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
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
  
   TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
   TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
   
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
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
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
  
  TransformedBasisValues<DataScalar,DeviceType>  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TransformedBasisValues<DataScalar,DeviceType> transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataRight);
  
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
