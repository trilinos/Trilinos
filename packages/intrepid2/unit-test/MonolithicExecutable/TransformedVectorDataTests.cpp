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


/** \file   TransformedVectorDataTests.cpp
    \brief  Tests against TransformedVectorData.
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
#include "Intrepid2_TransformedVectorData.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace
{
  using namespace Intrepid2;

  template<int spaceDim>
  void testVectorTransformation(const int &polyOrder, const int &meshWidth, Teuchos::FancyOStream &out, bool &success)
  {
    using Scalar = double;
    using PointScalar = double;
    
    const double relTol = 1e-12;
    const double absTol = 1e-12;
    
    auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
    
    auto lineBasis = Intrepid2::getLineBasis< Intrepid2::NodalBasisFamily<> >(fs, polyOrder);
    
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
    
    using ExecSpaceType = Kokkos::DefaultExecutionSpace;
    auto lineCubature = Intrepid2::DefaultCubatureFactory::create<ExecSpaceType>(lineTopo,polyOrder*2);
    int numPoints_1D = lineCubature->getNumPoints();
    ScalarView<PointScalar,ExecSpaceType> lineCubaturePoints("line cubature points",numPoints_1D,1);
    ScalarView<double,ExecSpaceType> lineCubatureWeights("line cubature weights", numPoints_1D);
    
    lineCubature->getCubature(lineCubaturePoints, lineCubatureWeights);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,ExecSpaceType> lineBasisValues    ("line basis values",      numFields_1D, numPoints_1D   );
    ScalarView<Scalar,ExecSpaceType> lineBasisGradValues("line basis grad values", numFields_1D, numPoints_1D, 1);
    
    // for now, we use 1D values to build up the 2D or 3D gradients
    // eventually, TensorBasis should offer a getValues() variant that returns tensor basis data
    lineBasis->getValues(lineBasisValues,     lineCubaturePoints, Intrepid2::OPERATOR_VALUE );
    lineBasis->getValues(lineBasisGradValues, lineCubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    // drop the trivial space dimension in line gradient values:
    Kokkos::resize(lineBasisGradValues, numFields_1D, numPoints_1D);
      
    Kokkos::Array<TensorData<Scalar,ExecSpaceType>, spaceDim> vectorComponents;
    
    for (int d=0; d<spaceDim; d++)
    {
      Kokkos::Array<Data<Scalar,ExecSpaceType>, spaceDim> gradComponent_d;
      for (int d2=0; d2<spaceDim; d2++)
      {
        if (d2 == d) gradComponent_d[d2] = Data<Scalar,ExecSpaceType>(lineBasisGradValues);
        else         gradComponent_d[d2] = Data<Scalar,ExecSpaceType>(lineBasisValues);
      }
      vectorComponents[d] = TensorData<Scalar,ExecSpaceType>(gradComponent_d);
    }
    VectorData<Scalar,ExecSpaceType> gradientValues(vectorComponents, false); // false: not axis-aligned
    
    CellGeometry<PointScalar,spaceDim,ExecSpaceType> cellNodes = uniformCartesianMesh<PointScalar,spaceDim,ExecSpaceType>(1.0, meshWidth);
    
    // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
    
    int pointsPerCell = 1;
    for (int d=0; d<spaceDim; d++)
    {
      pointsPerCell *= numPoints_1D;
    }
    
    auto jacobian = cellNodes.allocateJacobianData(pointsPerCell);
    auto jacobianDet = CellTools<ExecSpaceType>::allocateJacobianDet(jacobian);
    auto jacobianInv = CellTools<ExecSpaceType>::allocateJacobianInv(jacobian);
    cellNodes.setJacobian(                   jacobian, pointsPerCell);
    CellTools<ExecSpaceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<ExecSpaceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    auto transformedGradientData = FunctionSpaceTools<ExecSpaceType>::getHGRADtransformGRAD(jacobianInv, gradientValues);

    int numPoints = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numPoints *= numPoints_1D;
    }
    
    // now, compute transformed values in the classic, expanded way
    ScalarView<Scalar,ExecSpaceType> expandedTransformedGradValues("transformed grad values", numCells, numFields, numPoints, spaceDim);
    
    auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<> >(cellTopo, fs, polyOrder);
    
    // Allocate some intermediate containers
    ScalarView<Scalar,ExecSpaceType> basisValues    ("basis values", numFields, numPoints );
    ScalarView<Scalar,ExecSpaceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

    ScalarView<Scalar,ExecSpaceType> transformedGradValues("transformed grad values", numCells, numFields, numPoints, spaceDim);
    ScalarView<Scalar,ExecSpaceType> transformedWeightedGradValues("transformed weighted grad values", numCells, numFields, numPoints, spaceDim);
    
    using Kokkos::DefaultExecutionSpace;
    auto cubature = Intrepid2::DefaultCubatureFactory::create<DefaultExecutionSpace>(cellTopo,polyOrder*2);
    TEST_EQUALITY( numPoints, cubature->getNumPoints());
    ScalarView<PointScalar,ExecSpaceType> cubaturePoints("cubature points",numPoints,spaceDim);
    ScalarView<double,ExecSpaceType> cubatureWeights("cubature weights", numPoints);
    
    cubature->getCubature(cubaturePoints, cubatureWeights);
    
    basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
    basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
    
    const int numNodesPerCell = cellNodes.numNodesPerCell();
    ScalarView<PointScalar,ExecSpaceType> expandedCellNodes("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
    for (int cellOrdinal=0; cellOrdinal<numCells; cellOrdinal++)
    {
      for (int nodeOrdinal=0; nodeOrdinal<numNodesPerCell; nodeOrdinal++)
      {
        for (int d=0; d<spaceDim; d++)
        {
          expandedCellNodes(cellOrdinal,nodeOrdinal,d) = cellNodes(cellOrdinal,nodeOrdinal,d);
        }
      }
    }
    
    ScalarView<Scalar,ExecSpaceType> expandedJacobian("jacobian", numCells, numPoints, spaceDim, spaceDim);
    ScalarView<Scalar,ExecSpaceType> expandedJacobianInverse("jacobian inverse", numCells, numPoints, spaceDim, spaceDim);
    
    using CellTools = Intrepid2::CellTools<Kokkos::DefaultExecutionSpace>;
    using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<Kokkos::DefaultExecutionSpace>;
    
    CellTools::setJacobian(expandedJacobian, cubaturePoints, expandedCellNodes, cellTopo);
    CellTools::setJacobianInv(expandedJacobianInverse, expandedJacobian);
    
    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, expandedJacobianInverse, basisGradValues);
    
    testFloatingEquality4(transformedGradValues, transformedGradientData, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( TransformedVectorData, TransformedVector_1D_p1 )
  {
    const int spaceDim = 1;
    const int polyOrder = 1;
    const int meshWidth = 10;
    testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
  }
   
 TEUCHOS_UNIT_TEST( TransformedVectorData, TransformedVector_1D_p2 )
 {
   const int spaceDim = 1;
   const int polyOrder = 2;
   const int meshWidth = 10;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }
  
 TEUCHOS_UNIT_TEST( TransformedVectorData, TransformedVector_2D_p1 )
 {
   const int spaceDim = 2;
   const int polyOrder = 1;
   const int meshWidth = 3;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }
  
 TEUCHOS_UNIT_TEST( TransformedVectorData, TransformedVector_2D_p2 )
 {
   const int spaceDim = 2;
   const int polyOrder = 2;
   const int meshWidth = 3;
   testVectorTransformation<spaceDim>(polyOrder, meshWidth, out, success);
 }
} // anonymous namespace
