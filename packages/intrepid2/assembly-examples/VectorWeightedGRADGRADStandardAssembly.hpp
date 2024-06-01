//
//  VectorWeightedGRADGRADStandardAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/13/24.
//

#ifndef Intrepid2_VectorWeightedGRADGRADStandardAssembly_hpp
#define Intrepid2_VectorWeightedGRADGRADStandardAssembly_hpp

#include "JacobianFlopEstimate.hpp"
#include "Intrepid2_OrientationTools.hpp"

/** \file   VectorWeightedGRADGRADStandardAssembly.hpp
    \brief  Locally assembles a vector-weighted Poisson matrix -- an array of shape (C,F,F), with formulation (a dot grad e_i, b dot grad e_j), using standard Intrepid2 methods; these do not algorithmically exploit geometric structure.
 */

//! Version that uses the classic, generic Intrepid2 paths.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2 = spaceDim>
Intrepid2::ScalarView<Scalar,DeviceType> performStandardQuadratureVectorWeightedGRADGRAD(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry,
                                                                                         const int &polyOrder, int worksetSize,
                                                                                         Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight1,
                                                                                         Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight2,
                                                                                         double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  INTREPID2_TEST_FOR_EXCEPTION(vectorWeight1 == Teuchos::null, std::invalid_argument, "vectorWeight1 cannot be null");
  INTREPID2_TEST_FOR_EXCEPTION(vectorWeight2 == Teuchos::null, std::invalid_argument, "vectorWeight2 cannot be null");
  
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
  
  using CellTools = Intrepid2::CellTools<DeviceType>;
  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<DeviceType>;
  
  using namespace Intrepid2;
  
  using namespace std;
  // dimensions of the returned view are (C,F,F)
  auto fs = FUNCTION_SPACE_HGRAD;

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
  ScalarView<Scalar,DeviceType> basisValues    ("basis values", numFields, numPoints );
  ScalarView<Scalar,DeviceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

  ScalarView<Scalar,DeviceType> unorientedTransformedGradValues("unoriented transformed grad values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> transformedGradValues("transformed grad values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> transformedWeightedGradValues("transformed weighted grad values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,DeviceType> vectorWeightedTransformedGradValues("vector-weighted transformed grad values", worksetSize, numFields, numPoints);
  ScalarView<Scalar,DeviceType> vectorWeightedTransformedWeightedGradValues("vector-weighted transformed weighted grad values", worksetSize, numFields, numPoints);
  
  basis->getValues(basisValues,     cubaturePoints, OPERATOR_VALUE );
  basis->getValues(basisGradValues, cubaturePoints, OPERATOR_GRAD  );
  
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

  auto auView = getView<Scalar,DeviceType>("a_u", spaceDim);
  auto auViewHost = Kokkos::create_mirror(auView);

  for (int d=0; d<spaceDim; d++)
  {
    auViewHost(d) = (*vectorWeight1)[d];
  }
  Kokkos::deep_copy(auView, auViewHost);
  
  auto avView = getView<Scalar,DeviceType>("a_v", spaceDim);
  auto avViewHost = Kokkos::create_mirror(avView);
  for (int d=0; d<spaceDim; d++)
  {
    avViewHost(d) = (*vectorWeight2)[d];
  }
  Kokkos::deep_copy(avView, avViewHost);
  
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
    auto cellWorkset         = Kokkos::subview(expandedCellNodes, cellRange, Kokkos::ALL(), Kokkos::ALL());
    auto orientationsWorkset = Kokkos::subview(orientations, cellRange);
    
    if (numCellsInWorkset != worksetSize)
    {
      Kokkos::resize(jacobian,                        numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianInverse,                 numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianDeterminant,             numCellsInWorkset, numPoints);
      Kokkos::resize(cellMeasures,                    numCellsInWorkset, numPoints);
      Kokkos::resize(unorientedTransformedGradValues, numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedGradValues,           numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedWeightedGradValues,   numCellsInWorkset, numFields, numPoints, spaceDim);
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
    FunctionSpaceTools::HGRADtransformGRAD(unorientedTransformedGradValues, jacobianInverse, basisGradValues);
        // we want to exclude orientation application in the core integration timing -- this time gets reported as "Other"
    fstIntegrateCall->stop();
    OrientationTools<DeviceType>::modifyBasisByOrientation(transformedGradValues, unorientedTransformedGradValues,
                                                           orientationsWorkset, basis.get());
    fstIntegrateCall->start();
    
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim) * (spaceDim - 1) * 2.0; // 2: one multiply, one add per (P,D) entry in the contraction.
    FunctionSpaceTools::multiplyMeasure(transformedWeightedGradValues, cellMeasures, transformedGradValues);
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim); // multiply each entry of transformedGradValues: one flop for each.
        
    auto policy3 = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCellsInWorkset,numFields,numPoints});
    Kokkos::parallel_for("compute expanded_{u,v}TransformedGradValues", policy3,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
    {
      Scalar u_result = 0;
      Scalar v_result_weighted = 0;
      for (int d=0; d<spaceDim; d++)
      {
        u_result          += auView(d) *         transformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal,d);
        v_result_weighted += avView(d) * transformedWeightedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal,d);
      }
      vectorWeightedTransformedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal) = u_result;
      vectorWeightedTransformedWeightedGradValues(cellOrdinal,fieldOrdinal,pointOrdinal) = v_result_weighted;
    });
    
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim * 2 * 2); // 2 * 2: one multiply, one add per (D) entry, times 2 containers u and v
    
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    FunctionSpaceTools::integrate(cellStiffnessSubview, vectorWeightedTransformedGradValues, vectorWeightedTransformedWeightedGradValues);
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numFields) * double(numPoints * 2); // 2: one multiply, one add per P entry in the contraction.
    
    cellOffset += worksetSize;
  }
//  std::cout << "standard integration, approximateFlopCount: " << approximateFlopCount << std::endl;
  return cellStiffness;
}

#endif /* VectorWeightedGRADGRADStandardAssembly_h */
