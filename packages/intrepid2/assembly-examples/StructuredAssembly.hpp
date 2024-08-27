// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  StructuredAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 2/28/23.
//

#ifndef StructuredAssembly_h
#define StructuredAssembly_h

#include "JacobianFlopEstimate.hpp"
#include "Intrepid2_OrientationTools.hpp"

/** \file   StructuredAssembly.hpp
    \brief  Locally assembles a matrix of one basis integrated cell-wise against another -- an array of shape (C,F1,F2), with formulation (op1(e1_i), op2(e2_j)), using "structured" Intrepid2 methods; these algorithmically exploit geometric structure as expressed in the provided CellGeometry.
 */

namespace {
  template<class Scalar, typename DeviceType>
  Intrepid2::TransformedBasisValues<Scalar,DeviceType>
  transform(const Intrepid2::BasisValues<Scalar,DeviceType> &refValues,
            const Intrepid2::EFunctionSpace &fs, const Intrepid2::EOperator &op,
            const Intrepid2::Data<Scalar,DeviceType> &jacobian,    const Intrepid2::Data<Scalar,DeviceType> &jacobianDet,
            const Intrepid2::Data<Scalar,DeviceType> &jacobianInv, const Intrepid2::Data<Scalar,DeviceType> &jacobianDetInv,
            const Intrepid2::Data<Scalar,DeviceType> &jacobianDividedByJacobianDet)
  {
    using namespace Intrepid2;
    using FST = Intrepid2::FunctionSpaceTools<DeviceType>;
    switch (fs)
    {
      case EFunctionSpace::FUNCTION_SPACE_HGRAD:
      {
        switch (op)
        {
          case OPERATOR_VALUE:
            return FST::getHGRADtransformVALUE(jacobian.extent_int(0), refValues); // jacobian.extent_int(0): numCells
          case OPERATOR_GRAD:
            return FST::getHGRADtransformGRAD(jacobianInv, refValues);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      case EFunctionSpace::FUNCTION_SPACE_HCURL:
      {
        switch(op)
        {
          case OPERATOR_CURL:
          {
            const int spaceDim = jacobian.extent_int(2); // jacobian has shape (C,P,D,D)
            if (spaceDim == 2)
            {
              return FST::getHCURLtransformCURL2D(jacobianDetInv, refValues);
            }
            else
            {
              return FST::getHCURLtransformCURL(jacobianDividedByJacobianDet, refValues);
            }
          }
          case OPERATOR_VALUE:
            return FST::getHCURLtransformVALUE(jacobianInv, refValues);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      case EFunctionSpace::FUNCTION_SPACE_HDIV:
      {
        switch(op)
        {
          case OPERATOR_DIV:
            return FST::getHDIVtransformDIV(jacobianDetInv, refValues);
          case OPERATOR_VALUE:
            return FST::getHDIVtransformVALUE(jacobianDividedByJacobianDet, refValues);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
      case EFunctionSpace::FUNCTION_SPACE_HVOL:
      {
        switch (op)
        {
          case OPERATOR_VALUE:
            return FST::getHVOLtransformVALUE(jacobianDetInv, refValues );
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
    }
  }
}

//! General assembly for two arbitrary bases and ops that takes advantage of the new structured integration support, including support for sum factorization.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2>  // spaceDim and spaceDim2 should agree in value (differ in type)
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredAssembly(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &worksetSize,
                                                                   const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight1,
                                                                   const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight2,
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
  using IntegrationTools   = IntegrationTools<DeviceType>;
  // dimensions of the returned view are (C,F,F)
  
  Intrepid2::ScalarView<Intrepid2::Orientation,DeviceType> orientations("orientations", geometry.numCells() );
  geometry.orientations(orientations, 0, -1);
  
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis1 = getBasis< BasisFamily >(cellTopo, fs1, polyOrder1);
  auto basis2 = getBasis< BasisFamily >(cellTopo, fs2, polyOrder2);
  
  int numFields1 = basis1->getCardinality();
  int numFields2 = basis2->getCardinality();
  int numCells = geometry.numCells();
    
  // local stiffness matrix:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields1,numFields2);
  ScalarView<Scalar,DeviceType> worksetCellStiffness("cell stiffness workset matrices",worksetSize,numFields1,numFields2);

  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder1 + polyOrder2);
  auto tensorCubatureWeights = cubature->allocateCubatureWeights();
  TensorPoints<PointScalar,DeviceType> tensorCubaturePoints  = cubature->allocateCubaturePoints();
  
  cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
  
  BasisValues<Scalar,DeviceType> basis1Values = basis1->allocateBasisValues(tensorCubaturePoints, op1);
  basis1->getValues(basis1Values, tensorCubaturePoints, op1);
  
  BasisValues<Scalar,DeviceType> basis2Values = basis2->allocateBasisValues(tensorCubaturePoints, op2);
  basis2->getValues(basis2Values, tensorCubaturePoints, op2);
  
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  
  Data<PointScalar,DeviceType> jacobian       = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize);
  Data<PointScalar,DeviceType> jacobianDividedByJacobianDet = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize); // jacobianDividedByJacobianDet has same underlying structure as jacobian
  Data<PointScalar,DeviceType> jacobianDet    = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDetInv = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianInv    = CellTools<DeviceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,DeviceType> cellMeasures = geometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
  
  // lazily-evaluated transformed basis values (temporary to allow integralData allocation)
  auto transformedBasis1ValuesTemp = transform(basis1Values, fs1, op1, jacobian, jacobianDet, jacobianInv, jacobianDetInv, jacobianDividedByJacobianDet);
  auto transformedBasis2ValuesTemp = transform(basis2Values, fs2, op2, jacobian, jacobianDet, jacobianInv, jacobianDetInv, jacobianDividedByJacobianDet);
  auto integralData = IntegrationTools::allocateIntegralData(transformedBasis1ValuesTemp, cellMeasures, transformedBasis2ValuesTemp);
  
  const int numJacobianDataPoints = jacobian.getDataExtent(1); // data extent will be 1 for affine, numPoints for other cases
  const int numPoints             = jacobian.extent_int(1); // number of logical points
  
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
      jacobian.setExtent                    (CELL_DIM, numCellsInWorkset);
      jacobianDividedByJacobianDet.setExtent(CELL_DIM, numCellsInWorkset);
      jacobianDet.setExtent                 (CELL_DIM, numCellsInWorkset);
      jacobianDetInv.setExtent              (CELL_DIM, numCellsInWorkset);
      jacobianInv.setExtent                 (CELL_DIM, numCellsInWorkset);
      integralData.setExtent                (CELL_DIM, numCellsInWorkset);
      Kokkos::resize(worksetCellStiffness, numCellsInWorkset, numFields1, numFields2);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<DeviceType>::setJacobianDet         (jacobianDet,    jacobian);
    CellTools<DeviceType>::setJacobianDetInv      (jacobianDetInv, jacobian);
    CellTools<DeviceType>::setJacobianInv         (jacobianInv,    jacobian);
    CellTools<DeviceType>::setJacobianDividedByDet(jacobianDividedByJacobianDet, jacobian, jacobianDetInv);
    
    // lazily-evaluated transformed gradient values:
    auto transformedBasis1Values = transform(basis1Values, fs1, op1, jacobian, jacobianDet, jacobianInv, jacobianDetInv, jacobianDividedByJacobianDet);
    auto transformedBasis2Values = transform(basis2Values, fs2, op2, jacobian, jacobianDet, jacobianInv, jacobianDetInv, jacobianDividedByJacobianDet);
    
    if (vectorWeight1 != Teuchos::null)
    {
      ScalarView<Scalar,DeviceType> auView("a_u", spaceDim);
      auto auViewHost = Kokkos::create_mirror(auView);
      for (int d=0; d<spaceDim; d++)
      {
        auViewHost(d) = (*vectorWeight1)[d];
      }
      Kokkos::deep_copy(auView, auViewHost);
      
      Kokkos::Array<int,3> extents {numCellsInWorkset,numPoints,spaceDim};
      Kokkos::Array<DataVariationType,3> variationTypes {CONSTANT, CONSTANT, GENERAL};
      
      Data<Scalar,DeviceType> au_data(auView, extents, variationTypes);
      auto uTransform = Data<Scalar,DeviceType>::allocateMatVecResult(transformedBasis1Values.transform(), au_data, true);
      uTransform.storeMatVec(transformedBasis1Values.transform(), au_data, true); // true: transpose basis transform when multiplying
      transformedBasis1Values = Intrepid2::TransformedBasisValues<double, DeviceType>(uTransform, basis1Values);
      
      // TODO: modify transformIntegrateFlopCount to include an estimate for above mat-vecs (but note that these will not be a dominant cost, especially at high order).
    }
    
    if (vectorWeight2 != Teuchos::null)
    {
      ScalarView<Scalar,DeviceType> avView("a_v", spaceDim);
      auto avViewHost = Kokkos::create_mirror(avView);
      
      for (int d=0; d<spaceDim; d++)
      {
        avViewHost(d) = (*vectorWeight2)[d];
      }
      Kokkos::deep_copy(avView, avViewHost);
      
      Kokkos::Array<int,3> extents {numCellsInWorkset,numPoints,spaceDim};
      Kokkos::Array<DataVariationType,3> variationTypes {CONSTANT, CONSTANT, GENERAL};
      
      Data<Scalar,DeviceType> av_data(avView, extents, variationTypes);
      auto vTransform = Data<Scalar,DeviceType>::allocateMatVecResult(transformedBasis2Values.transform(), av_data, true);
      vTransform.storeMatVec(transformedBasis2Values.transform(), av_data, true); // true: transpose basis transform when multiplying
      transformedBasis2Values = Intrepid2::TransformedBasisValues<double, DeviceType>(vTransform, basis2Values);
      
      // TODO: modify transformIntegrateFlopCount to include an estimate for above mat-vecs (but note that these will not be a dominant cost, especially at high order).
    }
    
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorkset = 0;
    fstIntegrateCall->start();
    IntegrationTools::integrate(integralData, transformedBasis1Values, cellMeasures, transformedBasis2Values, sumInto, &approximateFlopCountIntegrateWorkset);
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    // modify integrals by orientations
    std::pair<int,int> cellRange = {startCell, endCell};
    auto orientationsWorkset = Kokkos::subview(orientations, cellRange);
    OrientationTools<DeviceType>::modifyMatrixByOrientation(worksetCellStiffness, integralData.getUnderlyingView(),
                                                            orientationsWorkset, basis1.get(), basis2.get());
    
    // copy into cellStiffness container.
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    Kokkos::deep_copy(cellStiffnessSubview, worksetCellStiffness);
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorkset;
    
    cellOffset += worksetSize;
  }
  return cellStiffness;

}

//! General assembly for two arbitrary bases and ops that takes advantage of the new structured integration support, including support for sum factorization.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredAssembly(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &worksetSize,
                                                                   const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1,
                                                                   const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2,
                                                                   double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > nullVectorWeight = Teuchos::null;
  
  return performStructuredAssembly<Scalar,BasisFamily,PointScalar,spaceDim,DeviceType>(geometry, worksetSize,
                                                                                       polyOrder1, fs1, op1, nullVectorWeight,
                                                                                       polyOrder2, fs2, op2, nullVectorWeight,
                                                                                       transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
}

#endif /* StructuredAssembly_h */
