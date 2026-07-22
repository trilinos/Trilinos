// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   TransformedBasisValuesTests.cpp
    \brief  Tests against TransformedBasisValues.
    \author Created by Nate Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_DefaultCubatureFactory.hpp>

#include <Intrepid2_Kernels.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellGeometryTestUtils.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace
{
  using namespace Intrepid2;

  //! vectorTransform [in]: whether the basis values being transformed are vector valued.
  //! identityTransform [in]: whether the tested transformation should be the identity transform, indicated by an empty (invalid) Data object.
  //! matrixTransform [in]: if true, indicates that the transform has shape (C,P,D,D); otherwise, it has shape (C,P)
  void testMultiplyByPointwiseWeights(bool vectorValues, bool identityTransform, bool matrixTransform, Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    const int spaceDim     = 3;
    const int numCells     = 3;
    const int numFields_1D = 1; // since transform is independent of the fields, does not matter for the test
    const int numPoints_1D = 10;
    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    // Allocate 1D values container; the actual values will not matter for the test, just the shape of the containers.
    ScalarView<Scalar,DeviceType> values1D("dummy 1D values", numFields_1D, numPoints_1D );

    BasisValues<Scalar, DeviceType> basisValues;
    if (vectorValues)
    {
      Kokkos::Array<TensorData<Scalar,DeviceType>, spaceDim> vectorComponents;
      
      for (int d=0; d<spaceDim; d++)
      {
        Kokkos::Array<Data<Scalar,DeviceType>, spaceDim> gradComponent_d;
        for (int d2=0; d2<spaceDim; d2++)
        {
          gradComponent_d[d2] = Data<Scalar,DeviceType>(values1D);
        }
        vectorComponents[d] = TensorData<Scalar,DeviceType>(gradComponent_d);
      }
      VectorData<Scalar,DeviceType> gradientVectorData(vectorComponents, false); // false: not axis-aligned
      basisValues = BasisValues<Scalar, DeviceType>(gradientVectorData);
    }
    else
    {
      Kokkos::Array<Data<Scalar,DeviceType>, spaceDim> tensorComponents;
      for (int d2=0; d2<spaceDim; d2++)
      {
        tensorComponents[d2] = Data<Scalar,DeviceType>(values1D);
      }
      TensorData<Scalar,DeviceType> tensorData(tensorComponents);
      basisValues = BasisValues<Scalar, DeviceType>(tensorData);
    }
    TransformedBasisValues<Scalar, DeviceType> transformedValues;
    Data<Scalar,DeviceType> transform;
    Kokkos::View<Scalar**, DeviceType> weightsView("weights", numCells, numPoints);
    
    ViewType<Scalar,DeviceType> expectedResultsView;
    {
      auto transformMatrixView = getFixedRankView<Scalar>("transform matrix", numCells, numPoints, spaceDim, spaceDim);
      auto transformMatrixViewHost = getHostCopy(transformMatrixView);
      
      auto transformView = getFixedRankView<Scalar>("(C,P) transform", numCells, numPoints);
      auto transformViewHost = getHostCopy(transformView);
      
      auto weightsViewHost = getHostCopy(weightsView);
      
      // we compute the expected results for both the matrix case and the other (pointwise weight) case, and assign the appropriate one to expectedResultsView at bottom of this block.
      ViewType<Scalar,DeviceType> matrixResultsView("matrix results",numCells,numPoints,spaceDim,spaceDim);
      auto matrixResultsViewHost = getHostCopy(matrixResultsView);
      
      ViewType<Scalar,DeviceType> pointResultsView("results",numCells,numPoints);
      auto pointResultsViewHost = getHostCopy(pointResultsView);
      
      for (int cellOrdinal=0; cellOrdinal<numCells; cellOrdinal++)
      {
        for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
        {
          const Scalar cpWeight = (0.5 + cellOrdinal) * (pointOrdinal + 0.1);
          weightsViewHost(cellOrdinal,pointOrdinal) = cpWeight;
          
          Scalar transformValue = identityTransform ? 1.0 : cellOrdinal * (pointOrdinal + 1);
          transformViewHost(cellOrdinal, pointOrdinal) = transformValue;
          
          pointResultsViewHost(cellOrdinal, pointOrdinal) = transformValue * cpWeight;
          
          for (int d1=0; d1<spaceDim; d1++)
          {
            for (int d2=0; d2<spaceDim; d2++)
            {
              Scalar transformMatrixValue;
              if (identityTransform)
              {
                transformMatrixValue = (d1 == d2) ? 1.0 : 0.0;
              }
              else
              {
                transformMatrixValue = (0.8 * cellOrdinal) + d1 * .1 + d2 * .3 - 0.4 * pointOrdinal;
              }
              transformMatrixViewHost(cellOrdinal,pointOrdinal,d1,d2) = transformMatrixValue;
              matrixResultsViewHost(cellOrdinal,pointOrdinal,d1,d2) = transformMatrixValue * cpWeight;
            }
          }
        }
      }
      // copy all our host views to device:
      Kokkos::deep_copy(transformMatrixView, transformMatrixViewHost);
      Kokkos::deep_copy(transformView,       transformViewHost);
      Kokkos::deep_copy(weightsView,         weightsViewHost);
      Kokkos::deep_copy(matrixResultsView,   matrixResultsViewHost);
      Kokkos::deep_copy(pointResultsView,    pointResultsViewHost);
      
      if (matrixTransform && !identityTransform)
      {
        expectedResultsView = matrixResultsView;
        transform = Data<Scalar, DeviceType>(transformMatrixView);
      }
      else
      {
        // identity matrix transform results in a pointwise transform
        expectedResultsView = pointResultsView;
        transform = Data<Scalar, DeviceType>(transformView);
      }
    }
    
    if (identityTransform)
    {
      transformedValues = TransformedBasisValues<Scalar, DeviceType>(numCells, basisValues);
    }
    else
    {
      transformedValues = TransformedBasisValues<Scalar, DeviceType>(transform, basisValues);
    }
    
    transformedValues.multiplyByPointwiseWeights(weightsView);
    
    auto actualResults = transformedValues.transform();
    
    TEST_EQUALITY(actualResults.rank(), expectedResultsView.rank());
    
    const double relTol = 1e-15;
    const double absTol = 1e-13;
    if (actualResults.rank() == expectedResultsView.rank())
    {
      if (actualResults.rank() == 4)
      {
        testFloatingEquality4(actualResults, expectedResultsView, relTol, absTol, out, success);
      }
      else if (actualResults.rank() == 2)
      {
        testFloatingEquality2(actualResults, expectedResultsView, relTol, absTol, out, success);
      }
      else
      {
        success = false;
        out << "ERROR: unexpected rank for actualResults.\n";
      }
    }
    else
    {
      out << "TEST FAILURE: rank of actual transform does not match expected.\n";
    }
  }

  template<int spaceDim>
  void testVectorTransformation(const int &polyOrder, const int &meshWidth, Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    using PointScalar = double;
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
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
    VectorData<Scalar,DeviceType> gradientVectorData(vectorComponents, false); // false: not axis-aligned
    BasisValues<Scalar, DeviceType> gradientValues(gradientVectorData);
    
    CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
    
    // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
    
    int pointsPerCell = 1;
    for (int d=0; d<spaceDim; d++)
    {
      pointsPerCell *= numPoints_1D;
    }
    
    auto jacobian = cellNodes.allocateJacobianData(pointsPerCell);
    auto jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
    auto jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
    cellNodes.setJacobian(                   jacobian, pointsPerCell);
    CellTools<DeviceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<DeviceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    auto transformedGradientData = FunctionSpaceTools<DeviceType>::getHGRADtransformGRAD(jacobianInv, gradientValues);

    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    // now, compute transformed values in the classic, expanded way
    ScalarView<Scalar,DeviceType> expandedTransformedGradValues("transformed grad values", numCells, numFields, numPoints, spaceDim);
    
    auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<DeviceType> >(cellTopo, fs, polyOrder);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
    ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

    ScalarView<Scalar,DeviceType> transformedGradValues("transformed grad values", numCells, numFields, numPoints, spaceDim);
    ScalarView<Scalar,DeviceType> transformedWeightedGradValues("transformed weighted grad values", numCells, numFields, numPoints, spaceDim);
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
    TEST_EQUALITY( numPoints, cubature->getNumPoints());
    ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
    ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
    
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
    basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
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
    
    ScalarView<Scalar,DeviceType> expandedJacobian("jacobian", numCells, numPoints, spaceDim, spaceDim);
    ScalarView<Scalar,DeviceType> expandedJacobianInverse("jacobian inverse", numCells, numPoints, spaceDim, spaceDim);
    
    using CellTools = Intrepid2::CellTools<DeviceType>;
    using ExecutionSpace = typename DeviceType::execution_space;
    using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<DeviceType>;
    
    CellTools::setJacobian(expandedJacobian, cubaturePoints, expandedCellNodes, cellTopo);
    CellTools::setJacobianInv(expandedJacobianInverse, expandedJacobian);
    
    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, expandedJacobianInverse, basisGradValues);
    
    testFloatingEquality4(transformedGradValues, transformedGradientData, relTol, absTol, out, success);
  }

  // testVectorWeightedTransformation tests against a (C,P,D) transformation of a gradient field.
  template<int spaceDim>
  void testWeightedVectorTransformation(const int &polyOrder, const int &meshWidth, Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;
    using PointScalar = double;
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
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
    VectorData<Scalar,DeviceType> gradientVectorData(vectorComponents, false); // false: not axis-aligned
    BasisValues<Scalar, DeviceType> gradientValues(gradientVectorData);
    
    CellGeometry<PointScalar,spaceDim,DeviceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(1.0, meshWidth);
    
    // goal here is to do a vector-weighted Poisson; i.e. (f a_u \cdot grad u, a_v \cdot grad v) on each cell
    
    int pointsPerCell = 1;
    for (int d=0; d<spaceDim; d++)
    {
      pointsPerCell *= numPoints_1D;
    }
    
    auto jacobian = cellNodes.allocateJacobianData(pointsPerCell);
    auto jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
    auto jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
    cellNodes.setJacobian(                   jacobian, pointsPerCell);
    CellTools<DeviceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<DeviceType>::setJacobianInv(jacobianInv, jacobian);
    
    auto auView = getView<Scalar,DeviceType>("a_u", spaceDim);
    auto auViewHost = Kokkos::create_mirror(auView);
    double weight = 1.0;
    for (int d=0; d<spaceDim; d++)
    {
      auViewHost(d) = weight;
      weight /= 2.0;
    }
    Kokkos::deep_copy(auView, auViewHost);
    
    auto avView = getView<Scalar,DeviceType>("a_v", spaceDim);
    auto avViewHost = Kokkos::create_mirror(avView);
    weight = 0.5;
    for (int d=0; d<spaceDim; d++)
    {
      avViewHost(d) = weight;
      weight *= 2.0;
    }
    Kokkos::deep_copy(avView, avViewHost);
    
    Data<Scalar,DeviceType> au_data(auView, Kokkos::Array<int,3>{numCells,pointsPerCell,spaceDim}, Kokkos::Array<DataVariationType,3>{CONSTANT,CONSTANT,GENERAL});
    Data<Scalar,DeviceType> av_data(avView, Kokkos::Array<int,3>{numCells,pointsPerCell,spaceDim}, Kokkos::Array<DataVariationType,3>{CONSTANT,CONSTANT,GENERAL});
    
    auto uTransform = Data<Scalar,DeviceType>::allocateMatVecResult(jacobianInv, au_data, true);
    auto vTransform = Data<Scalar,DeviceType>::allocateMatVecResult(jacobianInv, av_data, true);
    
    uTransform.storeMatVec(jacobianInv, au_data, true); // true: transpose jacobianInv when multiplying
    vTransform.storeMatVec(jacobianInv, av_data, true); // true: transpose jacobianInv when multiplying
    
    Intrepid2::TransformedBasisValues<double, DeviceType> utransformedBasisGradients(uTransform, gradientValues);
    Intrepid2::TransformedBasisValues<double, DeviceType> vtransformedBasisGradients(vTransform, gradientValues);
    
    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    // now, compute transformed values in the classic, expanded way
    ScalarView<Scalar,DeviceType> expanded_uTransformedGradValues("transformed a_u dot grad values", numCells, numFields, numPoints);
    ScalarView<Scalar,DeviceType> expanded_vTransformedGradValues("transformed a_v dot grad values", numCells, numFields, numPoints);
    
    auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<DeviceType> >(cellTopo, fs, polyOrder);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
    ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

    ScalarView<Scalar,DeviceType> transformedGradValues("transformed grad values", numCells, numFields, numPoints, spaceDim);
    ScalarView<Scalar,DeviceType> transformedWeightedGradValues("transformed weighted grad values", numCells, numFields, numPoints, spaceDim);
    
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
    TEST_EQUALITY( numPoints, cubature->getNumPoints());
    ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
    ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
    
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
    basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
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
    
    ScalarView<Scalar,DeviceType> expandedJacobian("jacobian", numCells, numPoints, spaceDim, spaceDim);
    ScalarView<Scalar,DeviceType> expandedJacobianInverse("jacobian inverse", numCells, numPoints, spaceDim, spaceDim);
    
    using CellTools = Intrepid2::CellTools<DeviceType>;
    using ExecutionSpace = typename DeviceType::execution_space;
    using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<DeviceType>;
    
    CellTools::setJacobian(expandedJacobian, cubaturePoints, expandedCellNodes, cellTopo);
    CellTools::setJacobianInv(expandedJacobianInverse, expandedJacobian);
    
    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, expandedJacobianInverse, basisGradValues);
    
    auto policy3 = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numFields,numPoints});
    Kokkos::parallel_for("compute expanded_{u,v}TransformedGradValues", policy3,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
    {
      Scalar u_result = 0;
      Scalar v_result = 0;
      for (int d=0; d<spaceDim; d++)
      {
        u_result += auView(d) * transformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal,d);
        v_result += avView(d) * transformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal,d);
      }
      expanded_uTransformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal) = u_result;
      expanded_vTransformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal) = v_result;
    });
    
    testFloatingEquality3(expanded_uTransformedGradValues, utransformedBasisGradients, relTol, absTol, out, success);
    testFloatingEquality3(expanded_vTransformedGradValues, vtransformedBasisGradients, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Vector_Identity_Matrix )
  {
    const bool vectorValues      = true;
    const bool identityTransform = true;
    const bool matrixTransform   = true;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }
  
  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Vector_Identity_Pointwise )
  {
    const bool vectorValues      = true;
    const bool identityTransform = true;
    const bool matrixTransform   = false;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Vector_NonIdentity_Matrix )
  {
    const bool vectorValues      = true;
    const bool identityTransform = false;
    const bool matrixTransform   = true;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }
  
  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Vector_NonIdentity_Pointwise )
  {
    const bool vectorValues      = true;
    const bool identityTransform = false;
    const bool matrixTransform   = false;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Scalar_Identity_Matrix )
  {
    const bool vectorValues      = false;
    const bool identityTransform = true;
    const bool matrixTransform   = true;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }
  
  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Scalar_Identity_Pointwise )
  {
    const bool vectorValues      = false;
    const bool identityTransform = true;
    const bool matrixTransform   = false;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Scalar_NonIdentity_Matrix )
  {
    const bool vectorValues      = false;
    const bool identityTransform = false;
    const bool matrixTransform   = true;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }
  
  TEUCHOS_UNIT_TEST( TransformedBasisValues, MultiplyByPointWeights_Scalar_NonIdentity_Pointwise )
  {
    const bool vectorValues      = false;
    const bool identityTransform = false;
    const bool matrixTransform   = false;
    testMultiplyByPointwiseWeights(vectorValues, identityTransform, matrixTransform, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedVector_1D_p1 )
  {
    const int spaceDim = 1;
    const int polyOrder = 1;
    const int meshWidth = 10;
    testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }
   
 TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedVector_1D_p2 )
 {
   const int spaceDim = 1;
   const int polyOrder = 2;
   const int meshWidth = 10;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }
  
 TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedVector_2D_p1 )
 {
   const int spaceDim = 2;
   const int polyOrder = 1;
   const int meshWidth = 3;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }
  
 TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedVector_2D_p2 )
 {
   const int spaceDim = 2;
   const int polyOrder = 2;
   const int meshWidth = 3;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedWeightedVector_1D_p1 )
  {
    const int spaceDim = 1;
    const int polyOrder = 1;
    const int meshWidth = 10;
    testWeightedVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }
   
  TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedWeightedVector_1D_p2 )
  {
   const int spaceDim = 1;
   const int polyOrder = 2;
   const int meshWidth = 10;
   testWeightedVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedWeightedVector_2D_p1 )
  {
   const int spaceDim = 2;
   const int polyOrder = 1;
   const int meshWidth = 3;
   testWeightedVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedBasisValues, TransformedWeightedVector_2D_p2 )
  {
   const int spaceDim = 2;
   const int polyOrder = 2;
   const int meshWidth = 3;
   testWeightedVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }
} // anonymous namespace
