// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  H1StructuredAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 3/1/22.
//

#ifndef H1StructuredAssembly_h
#define H1StructuredAssembly_h

#include "JacobianFlopEstimate.hpp"

/** \file   H1StructuredAssembly.hpp
    \brief  Locally assembles a matrix with the H^1 natural norm -- an array of shape (C,F,F), with formulation (e_i, e_j) + (grad e_i, grad e_j), using "structured" Intrepid2 methods; these algorithmically exploit geometric structure as expressed in the provided CellGeometry.
 */

//! Version that takes advantage of new structured integration support, including sum factorization.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadratureH1(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
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
  auto fs = FUNCTION_SPACE_HGRAD;
  
  Intrepid2::ScalarView<Intrepid2::Orientation,DeviceType> orientations("orientations", geometry.numCells() );
  geometry.orientations(orientations, 0, -1);
  
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis = getBasis< BasisFamily >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numCells = geometry.numCells();
    
  // local stiffness matrix:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  ScalarView<Scalar,DeviceType> worksetCellStiffness("cell stiffness workset matrices",worksetSize,numFields,numFields);
  
  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder*2);
  auto tensorCubatureWeights = cubature->allocateCubatureWeights();
  TensorPoints<PointScalar,DeviceType> tensorCubaturePoints  = cubature->allocateCubaturePoints();
  
  cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
  
  EOperator opGrad = OPERATOR_GRAD;
  BasisValues<Scalar,DeviceType> gradientValues = basis->allocateBasisValues(tensorCubaturePoints, opGrad);
  basis->getValues(gradientValues, tensorCubaturePoints, opGrad);
  
  EOperator opValue = OPERATOR_VALUE;
  BasisValues<Scalar,DeviceType> basisValues = basis->allocateBasisValues(tensorCubaturePoints, opValue);
  basis->getValues(basisValues, tensorCubaturePoints, opValue);
  
  // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
    
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  
  Data<PointScalar,DeviceType> jacobian = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize);
  Data<PointScalar,DeviceType> jacobianDet = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,DeviceType> cellMeasures = geometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
  
  // lazily-evaluated transformed gradient values (temporary to allow integralData allocation)
  auto transformedGradientValuesTemp = FunctionSpaceTools::getHGRADtransformGRAD(jacobianInv, gradientValues);
  auto integralData = IntegrationTools::allocateIntegralData(transformedGradientValuesTemp, cellMeasures, transformedGradientValuesTemp);
  
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
      jacobian.setExtent(    CELL_DIM, numCellsInWorkset);
      jacobianDet.setExtent( CELL_DIM, numCellsInWorkset);
      jacobianInv.setExtent( CELL_DIM, numCellsInWorkset);
      integralData.setExtent(CELL_DIM, numCellsInWorkset);
      Kokkos::resize(worksetCellStiffness, numCellsInWorkset, numFields, numFields);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<DeviceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<DeviceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    auto transformedGradientValues = FunctionSpaceTools::getHGRADtransformGRAD(jacobianInv, gradientValues);
    
    auto transformedBasisValues = FunctionSpaceTools::getHGRADtransformVALUE(numCellsInWorkset, basisValues);
    
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorksetGRADGRAD = 0;
    double approximateFlopCountIntegrateWorksetVALUEVALUE = 0;
    fstIntegrateCall->start();
    // (grad,grad)
    IntegrationTools::integrate(integralData, transformedGradientValues, cellMeasures, transformedGradientValues, sumInto, &approximateFlopCountIntegrateWorksetGRADGRAD);
    
    // (value,value)
    sumInto = true; // add to what we already accumulated for (grad,grad)
    IntegrationTools::integrate(integralData, transformedBasisValues, cellMeasures, transformedBasisValues, sumInto, &approximateFlopCountIntegrateWorksetVALUEVALUE);
    
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    // modify integrals by orientations
    std::pair<int,int> cellRange = {startCell, endCell};
    auto orientationsWorkset = Kokkos::subview(orientations, cellRange);
    OrientationTools<DeviceType>::modifyMatrixByOrientation(worksetCellStiffness, integralData.getUnderlyingView(),
                                                            orientationsWorkset, basis.get(), basis.get());
    
    // copy into cellStiffness container.
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    Kokkos::deep_copy(cellStiffnessSubview, worksetCellStiffness);
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorksetGRADGRAD + approximateFlopCountIntegrateWorksetVALUEVALUE;
    
    cellOffset += worksetSize;
  }
  return cellStiffness;
}

#endif /* GRADGRADStructuredAssembly_h */
