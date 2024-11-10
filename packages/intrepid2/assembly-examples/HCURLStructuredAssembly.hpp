// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  HCURLStructuredAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 3/24/22.
//

#ifndef HCURLStructuredAssembly_h
#define HCURLStructuredAssembly_h

#include "JacobianFlopEstimate.hpp"

/** \file   HCURLStructuredAssembly.hpp
    \brief  Locally assembles a matrix with the H(curl) natural norm -- an array of shape (C,F,F), with formulation (e_i, e_j) + (curl e_i, curl e_j), using "structured" Intrepid2 methods; these algorithmically exploit geometric structure as expressed in the provided CellGeometry.
 */

//! Version that takes advantage of new structured integration support, including sum factorization.  Computes H(curl) norm: (curl, curl) + (value,value).
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadratureHCURL(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
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
  auto fs = FUNCTION_SPACE_HCURL;
  
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
  
  EOperator opCurl = OPERATOR_CURL;
  BasisValues<Scalar,DeviceType> curlValues = basis->allocateBasisValues(tensorCubaturePoints, opCurl);
  basis->getValues(curlValues, tensorCubaturePoints, opCurl);
  
  EOperator opValue = OPERATOR_VALUE;
  BasisValues<Scalar,DeviceType> basisValues = basis->allocateBasisValues(tensorCubaturePoints, opValue);
  basis->getValues(basisValues, tensorCubaturePoints, opValue);
  
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  
  Data<PointScalar,DeviceType> jacobian = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize);
  Data<PointScalar,DeviceType> jacobianDet        = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDetInverse = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDividedByJacobianDet = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize); // jacobianDividedByJacobianDet has same underlying structure as jacobian
  Data<PointScalar,DeviceType> jacobianInv = CellTools<DeviceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,DeviceType> cellMeasures = geometry.allocateCellMeasure(jacobianDetInverse, tensorCubatureWeights);
  
  // lazily-evaluated transformed curl values (temporary to allow integralData allocation)
  auto transformedCurlValuesTemp = FunctionSpaceTools::getHCURLtransformCURL(jacobianDividedByJacobianDet, curlValues);
  auto integralData = IntegrationTools::allocateIntegralData(transformedCurlValuesTemp, cellMeasures, transformedCurlValuesTemp);
  
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
      Kokkos::resize(worksetCellStiffness, numCellsInWorkset, numFields, numFields);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<DeviceType>::setJacobianDet(         jacobianDet,                  jacobian);
    CellTools<DeviceType>::setJacobianDetInv(      jacobianDetInverse,           jacobian);
    CellTools<DeviceType>::setJacobianInv(         jacobianInv,                  jacobian);
    CellTools<DeviceType>::setJacobianDividedByDet(jacobianDividedByJacobianDet, jacobian, jacobianDetInverse);
    
    // lazily-evaluated transformed curl, values:
    TransformedBasisValues<Scalar,DeviceType> transformedCurlValues;
    if (spaceDim == 2)
    {
      transformedCurlValues = FunctionSpaceTools::getHCURLtransformCURL(jacobianDetInverse, curlValues);
    }
    else
    {
      transformedCurlValues = FunctionSpaceTools::getHCURLtransformCURL(jacobianDividedByJacobianDet, curlValues);
    }
    auto transformedBasisValues = FunctionSpaceTools::getHCURLtransformVALUE(jacobianInv, basisValues);
    
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorksetCURLCURL = 0;
    double approximateFlopCountIntegrateWorksetVALUEVALUE = 0;
    fstIntegrateCall->start();
    // (curl,curl)
    IntegrationTools::integrate(integralData, transformedCurlValues, cellMeasures, transformedCurlValues, sumInto, &approximateFlopCountIntegrateWorksetCURLCURL);
    
    // (value,value)
    sumInto = true; // add to what we already accumulated for (curl,curl)
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
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorksetCURLCURL + approximateFlopCountIntegrateWorksetVALUEVALUE;
    
    cellOffset += worksetSize;
  }
  return cellStiffness;
}

#endif /* HCURLStructuredAssembly_h */
