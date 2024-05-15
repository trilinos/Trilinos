//
//  VectorWeightedGRADGRADStructuredAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/13/24.
//

#ifndef VectorWeightedGRADGRADStructuredAssembly_h
#define VectorWeightedGRADGRADStructuredAssembly_h

#include "JacobianFlopEstimate.hpp"
#include "Intrepid2_OrientationTools.hpp"

/** \file   VectorWeightedGRADGRADStructuredAssembly.hpp
    \brief  Locally assembles a vector-weighted Poisson matrix -- an array of shape (C,F,F), with formulation (a dot grad e_i, b dot grad e_j), using "structured" Intrepid2 methods; these algorithmically exploit geometric structure as expressed in the provided CellGeometry.
 */

//! Version that takes advantage of new structured integration support, including sum factorization.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadratureVectorWeightedGRADGRAD(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                                                                           Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight1,
                                                                                           Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight2,
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
  
  EOperator op = OPERATOR_GRAD;
  BasisValues<Scalar,DeviceType> gradientValues = basis->allocateBasisValues(tensorCubaturePoints, op);
  basis->getValues(gradientValues, tensorCubaturePoints, op);
  
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
  
  const int numJacobianDataPoints = jacobian.getDataExtent(1); // data extent will be 1 for affine, numPoints for other cases
  const int numPoints             = jacobian.extent_int(1); // logical point count
  
  // TODO: make the below determination accurate for diagonal/block-diagonal cases… (right now, will overcount)
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numJacobianDataPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numJacobianDataPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numJacobianDataPoints);
  
  transformIntegrateFlopCount = 0;
  jacobianCellMeasureFlopCount  = numCells * flopsPerJacobianPerCell;    // jacobian itself
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianInvPerCell; // inverse
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianDetPerCell; // determinant
  jacobianCellMeasureFlopCount += numCells * numJacobianDataPoints; // cell measure: (C,P) gets weighted with cubature weights of shape (P)
  
  auto refData = geometry.getJacobianRefData(tensorCubaturePoints);
  
  ScalarView<Scalar,DeviceType> auView("a_u", spaceDim);
  auto auViewHost = Kokkos::create_mirror(auView);
  
  for (int d=0; d<spaceDim; d++)
  {
    auViewHost(d) = (*vectorWeight1)[d];
  }
  Kokkos::deep_copy(auView, auViewHost);
  
  ScalarView<Scalar,DeviceType> avView("a_v", spaceDim);
  auto avViewHost = Kokkos::create_mirror(avView);
  
  for (int d=0; d<spaceDim; d++)
  {
    avViewHost(d) = (*vectorWeight2)[d];
  }
  Kokkos::deep_copy(avView, avViewHost);
  Data<Scalar,DeviceType> au_data(auView, Kokkos::Array<int,3>{worksetSize,numPoints,spaceDim}, Kokkos::Array<DataVariationType,3>{CONSTANT,CONSTANT,GENERAL});
  Data<Scalar,DeviceType> av_data(avView, Kokkos::Array<int,3>{worksetSize,numPoints,spaceDim}, Kokkos::Array<DataVariationType,3>{CONSTANT,CONSTANT,GENERAL});
  
  auto uTransform = Data<Scalar,DeviceType>::allocateMatVecResult(jacobianInv, au_data, true);
  auto vTransform = Data<Scalar,DeviceType>::allocateMatVecResult(jacobianInv, av_data, true);
  
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
      jacobian.setExtent    (CELL_DIM, numCellsInWorkset);
      jacobianDet.setExtent (CELL_DIM, numCellsInWorkset);
      jacobianInv.setExtent (CELL_DIM, numCellsInWorkset);
      integralData.setExtent(CELL_DIM, numCellsInWorkset);
      au_data.setExtent     (CELL_DIM, numCellsInWorkset);
      av_data.setExtent     (CELL_DIM, numCellsInWorkset);
      uTransform.setExtent  (CELL_DIM, numCellsInWorkset);
      vTransform.setExtent  (CELL_DIM, numCellsInWorkset);
      
      Kokkos::resize(worksetCellStiffness, numCellsInWorkset, numFields, numFields);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<DeviceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<DeviceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    uTransform.storeMatVec(jacobianInv, au_data, true); // true: transpose jacobianInv when multiplying
    vTransform.storeMatVec(jacobianInv, av_data, true); // true: transpose jacobianInv when multiplying
    
    Intrepid2::TransformedBasisValues<double, DeviceType> uTransformedGradientValues(uTransform, gradientValues);
    Intrepid2::TransformedBasisValues<double, DeviceType> vTransformedGradientValues(vTransform, gradientValues);
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorkset = 0;
    fstIntegrateCall->start();
    IntegrationTools::integrate(integralData, uTransformedGradientValues, cellMeasures, vTransformedGradientValues, sumInto, &approximateFlopCountIntegrateWorkset);
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
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorkset;
    
    cellOffset += worksetSize;
  }
  return cellStiffness;
}

#endif /* VectorWeightedGRADGRADStructuredAssembly_h */
