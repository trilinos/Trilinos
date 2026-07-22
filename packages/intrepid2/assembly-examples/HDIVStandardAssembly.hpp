// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  HDIVStandardAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 3/7/22.
//

#ifndef Intrepid2_HDIVStandardAssembly_hpp
#define Intrepid2_HDIVStandardAssembly_hpp

#include "JacobianFlopEstimate.hpp"

/** \file   HDIVStandardAssembly.hpp
    \brief  Locally assembles a matrix with the H(div) natural norm -- an array of shape (C,F,F), with formulation (e_i, e_j) + (div e_i, div e_j), using standard Intrepid2 methods; these do not algorithmically exploit geometric structure.
 */

//! Version that uses the classic, generic Intrepid2 paths.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStandardQuadratureHDIV(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry,
                                                                       const int &polyOrder, int worksetSize,
                                                                       double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  using namespace Intrepid2;
  using ExecutionSpace = typename DeviceType::execution_space;
  
  int numVertices = 1;
  for (int d=0; d<spaceDim; d++)
  {
    numVertices *= 2;
  }
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");
  initialSetupTimer->start();
  
  using CellTools = CellTools<DeviceType>;
  using FunctionSpaceTools = FunctionSpaceTools<DeviceType>;
  
  using namespace std;
  // dimensions of the returned view are (C,F,F)
  auto fs = FUNCTION_SPACE_HDIV;
  
  Intrepid2::ScalarView<Intrepid2::Orientation,DeviceType> orientations("orientations", geometry.numCells() );
  geometry.orientations(orientations, 0, -1);

  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis = getBasis< BasisFamily >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numCells = geometry.numCells();
  
  if (worksetSize > numCells) worksetSize = numCells;
  
  // local stiffness matrices:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  
  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
  int numPoints = cubature->getNumPoints();
  ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
  ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
  
  cubature->getCubature(cubaturePoints, cubatureWeights);
  
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numPoints);
  
  // Allocate some intermediate containers
  ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints, spaceDim );
  ScalarView<Scalar,DeviceType> basisDivValues("basis div values", numFields, numPoints);

  ScalarView<Scalar,DeviceType> unorientedTransformedDivValues("unoriented transformed div values", worksetSize, numFields, numPoints);
  ScalarView<Scalar,DeviceType> transformedDivValues("transformed div values", worksetSize, numFields, numPoints);
  ScalarView<Scalar,DeviceType> transformedWeightedDivValues("transformed weighted div values", worksetSize, numFields, numPoints);
  
  ScalarView<Scalar,DeviceType> unorientedTransformedBasisValues("unoriented transformed basis values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> transformedBasisValues("transformed basis values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> transformedWeightedBasisValues("transformed weighted basis values", worksetSize, numFields, numPoints, spaceDim);
  
  basis->getValues(basisValues,    cubaturePoints, OPERATOR_VALUE );
  basis->getValues(basisDivValues, cubaturePoints, OPERATOR_DIV  );
  
  const int numNodesPerCell = geometry.numNodesPerCell();
  ScalarView<PointScalar,DeviceType> expandedCellNodes("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0,numCells),
  KOKKOS_LAMBDA (const int &cellOrdinal) {
    for (int nodeOrdinal=0; nodeOrdinal<numNodesPerCell; nodeOrdinal++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        expandedCellNodes(cellOrdinal,nodeOrdinal,d) = geometry(cellOrdinal,nodeOrdinal,d);
      }
    }
  });
  
  ScalarView<Scalar,DeviceType> cellMeasures("cell measures", worksetSize, numPoints);
  ScalarView<Scalar,DeviceType> jacobianDeterminant("jacobian determinant", worksetSize, numPoints);
  ScalarView<Scalar,DeviceType> jacobian("jacobian", worksetSize, numPoints, spaceDim, spaceDim);
  ScalarView<Scalar,DeviceType> jacobianInverse("jacobian inverse", worksetSize, numPoints, spaceDim, spaceDim);

  initialSetupTimer->stop();
  
  transformIntegrateFlopCount  = 0;
  jacobianCellMeasureFlopCount  = numCells * flopsPerJacobianPerCell;    // jacobian itself
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianInvPerCell; // inverse
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianDetPerCell; // determinant
  jacobianCellMeasureFlopCount += numCells * numPoints; // cell measure: (C,P) gets weighted with cubature weights of shape (P)
  
  int cellOffset = 0;
  while (cellOffset < numCells)
  {
    int startCell         = cellOffset;
    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
    
    std::pair<int,int> cellRange = {startCell, startCell+numCellsInWorkset};
    auto cellWorkset = Kokkos::subview(expandedCellNodes, cellRange, Kokkos::ALL(), Kokkos::ALL());
    auto orientationsWorkset = Kokkos::subview(orientations, cellRange);
    
    if (numCellsInWorkset != worksetSize)
    {
      Kokkos::resize(jacobian,                         numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianInverse,                  numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianDeterminant,              numCellsInWorkset, numPoints);
      Kokkos::resize(cellMeasures,                     numCellsInWorkset, numPoints);
      Kokkos::resize(unorientedTransformedDivValues,   numCellsInWorkset, numFields, numPoints);
      Kokkos::resize(unorientedTransformedBasisValues, numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedBasisValues,           numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedDivValues,             numCellsInWorkset, numFields, numPoints);
      Kokkos::resize(transformedWeightedBasisValues,   numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedWeightedDivValues,     numCellsInWorkset, numFields, numPoints);
    }
    jacobianAndCellMeasureTimer->start();
    CellTools::setJacobian(jacobian, cubaturePoints, cellWorkset, cellTopo); // accounted for outside loop, as numCells * flopsPerJacobianPerCell.
    CellTools::setJacobianInv(jacobianInverse, jacobian);
    CellTools::setJacobianDet(jacobianDeterminant, jacobian);
    
    FunctionSpaceTools::computeCellMeasure(cellMeasures, jacobianDeterminant, cubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
        
    // because structured integration performs transformations within integrate(), to get a fairer comparison here we include the transformation calls.
    fstIntegrateCall->start();
    FunctionSpaceTools::HDIVtransformDIV(unorientedTransformedDivValues, jacobianDeterminant, basisDivValues);
    // we want to exclude orientation application in the core integration timing -- this time gets reported as "Other"
    fstIntegrateCall->stop();
    OrientationTools<DeviceType>::modifyBasisByOrientation(transformedDivValues, unorientedTransformedDivValues,
                                                           orientationsWorkset, basis.get());
    fstIntegrateCall->start();
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim) * (spaceDim - 1) * 2.0; // 2: one multiply, one add per (P,D) entry in the contraction.
    FunctionSpaceTools::multiplyMeasure(transformedWeightedDivValues, cellMeasures, transformedDivValues);
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim); // multiply each entry of transformedDivValues: one flop for each.
        
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    FunctionSpaceTools::integrate(cellStiffnessSubview, transformedDivValues, transformedWeightedDivValues);
    ExecutionSpace().fence();
    FunctionSpaceTools::HDIVtransformVALUE(unorientedTransformedBasisValues, jacobian, jacobianDeterminant, basisValues);
    // we want to exclude orientation application in the core integration timing -- this time gets reported as "Other"
    fstIntegrateCall->stop();
    OrientationTools<DeviceType>::modifyBasisByOrientation(transformedBasisValues, unorientedTransformedBasisValues,
                                                           orientationsWorkset, basis.get());
    fstIntegrateCall->start();
    FunctionSpaceTools::multiplyMeasure(transformedWeightedBasisValues, cellMeasures, transformedBasisValues);
    bool sumInto = true; // add the (value,value) integral to the (div,div) that we've already integrated
    FunctionSpaceTools::integrate(cellStiffnessSubview, transformedBasisValues, transformedWeightedBasisValues, sumInto);
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    // record flops for (div,div)
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numFields) * double(numPoints * 2 - 1); // one multiply, add for each (P) entry in the contraction.
    // for (value,value):
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numFields) * double(numPoints * spaceDim * 2 - 1); // 2: one multiply, one add per (P,D) entry in the contraction; but you do one less add than the number of multiplies, hence the minus 1.
    
    cellOffset += worksetSize;
  }
//  std::cout << "standard integration, approximateFlopCount: " << approximateFlopCount << std::endl;
  return cellStiffness;
}

#endif /* HDIVStandardAssembly_h */
