//
//  HDIVStructuredAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 3/7/22.
//

#ifndef HDIVStructuredAssembly_h
#define HDIVStructuredAssembly_h

#include "JacobianFlopEstimate.hpp"

/** \file   HDIVStructuredAssembly.hpp
    \brief  Locally assembles a matrix with the H(div) natural norm -- an array of shape (C,F,F), with formulation (e_i, e_j) + (div e_i, div e_j), using "structured" Intrepid2 methods; these algorithmically exploit geometric structure as expressed in the provided CellGeometry.
 */

//! Version that takes advantage of new structured integration support, including sum factorization.  Computes H(div) norm: (div, div) + (value,value).
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadratureHDIV(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                                                         double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  using namespace Intrepid2;
  
  using ExecutionSpace = typename DeviceType::execution_space;
  
  int numVertices = 1;
  for (int d=0; d<spaceDim; d++)
  {
    numVertices *= 2;
  }
  
  auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");
  initialSetupTimer->start();
  using namespace std;
  using FunctionSpaceTools = FunctionSpaceTools<DeviceType>;
  using IntegrationTools   = IntegrationTools<DeviceType>;
  // dimensions of the returned view are (C,F,F)
  auto fs = FUNCTION_SPACE_HDIV;
  
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis = getBasis< BasisFamily >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numCells = geometry.numCells();
    
  // local stiffness matrix:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  
  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
  auto tensorCubatureWeights = cubature->allocateCubatureWeights();
  TensorPoints<PointScalar,DeviceType> tensorCubaturePoints  = cubature->allocateCubaturePoints();
  
  cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
  
  EOperator opDiv = OPERATOR_DIV;
  BasisValues<Scalar,DeviceType> divValues = basis->allocateBasisValues(tensorCubaturePoints, opDiv);
  basis->getValues(divValues, tensorCubaturePoints, opDiv);
  
  EOperator opValue = OPERATOR_VALUE;
  BasisValues<Scalar,DeviceType> basisValues = basis->allocateBasisValues(tensorCubaturePoints, opValue);
  basis->getValues(basisValues, tensorCubaturePoints, opValue);
  
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  
  Data<PointScalar,DeviceType> jacobian = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize);
  Data<PointScalar,DeviceType> jacobianDet        = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDetInverse = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDividedByJacobianDet = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize); // jacobianTimesJacobianDet has same underlying structure as jacobian
  Data<PointScalar,DeviceType> jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,DeviceType> cellMeasures = geometry.allocateCellMeasure(jacobianDetInverse, tensorCubatureWeights);
  
  Data<PointScalar,DeviceType> jacobianDetInverseExtended; // container with same underlying data as jacobianDet, but extended with CONSTANT type to have same logical shape as Jacobian
  // setup jacobianDetInverseExtended
  {
    auto variationTypes = jacobianDetInverse.getVariationTypes(); // defaults to CONSTANT in ranks beyond the rank of the container; this is what we want for our new extents
    auto extents        = jacobian.getExtents();
    
    jacobianDetInverseExtended = jacobianDetInverse.shallowCopy(jacobian.rank(), extents, variationTypes);
  }
  
  // lazily-evaluated transformed div values (temporary to allow integralData allocation)
  auto transformedDivValuesTemp = FunctionSpaceTools::getHDIVtransformDIV(jacobianDetInverse, divValues);
  auto integralData = IntegrationTools::allocateIntegralData(transformedDivValuesTemp, cellMeasures, transformedDivValuesTemp);
  
  const int numPoints = jacobian.getDataExtent(1); // data extent will be 1 for affine, numPoints for other cases
  
  // TODO: make the below determination accurate for diagonal/block-diagonal cases… (right now, will overcount)
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numPoints);
  
  transformIntegrateFlopCount = 0;
  jacobianCellMeasureFlopCount  = numCells * flopsPerJacobianPerCell;    // jacobian itself
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianInvPerCell; // inverse
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianDetPerCell; // determinant
  jacobianCellMeasureFlopCount += numCells * numPoints; // cell measure: (C,P) gets weighted with cubature weights of shape (P)
  
  auto refData = geometry.getJacobianRefData(tensorCubaturePoints);
  
  initialSetupTimer->stop();
  while (cellOffset < numCells)
  {
    int startCell         = cellOffset;
    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
    int endCell           = numCellsInWorkset + startCell;
    
    jacobianAndCellMeasureTimer->start();
    if (numCellsInWorkset != worksetSize)
    {
      const int CELL_DIM = 0; // first dimension corresponds to cell
      jacobian.setExtent(                    CELL_DIM, numCellsInWorkset);
      jacobianDet.setExtent(                 CELL_DIM, numCellsInWorkset);
      jacobianDetInverse.setExtent(          CELL_DIM, numCellsInWorkset);
      jacobianInv.setExtent(                 CELL_DIM, numCellsInWorkset);
      integralData.setExtent(                CELL_DIM, numCellsInWorkset);
      jacobianDividedByJacobianDet.setExtent(CELL_DIM, numCellsInWorkset);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<DeviceType>::setJacobianDet(   jacobianDet,        jacobian);
    CellTools<DeviceType>::setJacobianDetInv(jacobianDetInverse, jacobian);
    CellTools<DeviceType>::setJacobianInv(   jacobianInv,        jacobian);
    
    // compute the jacobian divided by its determinant
    jacobianDividedByJacobianDet.storeInPlaceProduct(jacobian,jacobianDetInverseExtended);
    
    // lazily-evaluated transformed div, values:
    auto transformedDivValues   = FunctionSpaceTools::getHDIVtransformDIV(jacobianDetInverse, divValues);
    auto transformedBasisValues = FunctionSpaceTools::getHDIVtransformVALUE(jacobianDividedByJacobianDet, basisValues);
    
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorksetDIVDIV = 0;
    double approximateFlopCountIntegrateWorksetVALUEVALUE = 0;
    fstIntegrateCall->start();
    // (div,div)
    IntegrationTools::integrate(integralData, transformedDivValues, cellMeasures, transformedDivValues, sumInto, &approximateFlopCountIntegrateWorksetDIVDIV);
    
    // (value,value)
    sumInto = true; // add to what we already accumulated for (div,div)
    IntegrationTools::integrate(integralData, transformedBasisValues, cellMeasures, transformedBasisValues, sumInto, &approximateFlopCountIntegrateWorksetVALUEVALUE);
    
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    // copy into cellStiffness container.  (Alternately, do something like allocateIntegralData, but outside this loop, and take a subview to construct the workset integralData.)
    if (integralData.getUnderlyingViewRank() == 3)
    {
      std::pair<int,int> cellRange = {startCell, endCell};
      auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
      Kokkos::deep_copy(cellStiffnessSubview, integralData.getUnderlyingView3());
    }
    else // underlying view rank is 2; copy to each cell in destination stiffness matrix
    {
      auto integralView2 = integralData.getUnderlyingView2();
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCellsInWorkset,numFields,numFields});
      Kokkos::parallel_for("copy uniform data to expanded container", policy,
                       KOKKOS_LAMBDA (const int &cellOrdinal, const int &leftFieldOrdinal, const int &rightFieldOrdinal) {
        cellStiffness(startCell + cellOrdinal, leftFieldOrdinal, rightFieldOrdinal) = integralView2(leftFieldOrdinal,rightFieldOrdinal);
      });
    }
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorksetDIVDIV + approximateFlopCountIntegrateWorksetVALUEVALUE;
    
    cellOffset += worksetSize;
  }
  return cellStiffness;
}

#endif /* HDIVStructuredAssembly_h */
