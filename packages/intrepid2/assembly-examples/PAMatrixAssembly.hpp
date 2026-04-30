// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  PAMatrixAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 10/17/24.
//

#ifndef PAMatrixAssembly_h
#define PAMatrixAssembly_h

#include "JacobianFlopEstimate.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_PAMatrix.hpp"

/** \file   PAMatrixAssembly.hpp
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
Intrepid2::PAMatrix<DeviceType,Scalar> constructPAMatrix(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry,
                                                         const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight1,
                                                         const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight2)
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
    
  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder1 + polyOrder2);
  auto tensorCubatureWeights = cubature->allocateCubatureWeights();
  TensorPoints<PointScalar,DeviceType> tensorCubaturePoints  = cubature->allocateCubaturePoints();
  
  cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
  
  BasisValues<Scalar,DeviceType> basis1Values = basis1->allocateBasisValues(tensorCubaturePoints, op1);
  basis1->getValues(basis1Values, tensorCubaturePoints, op1);
  basis1Values.setBasis(basis1);
  
  BasisValues<Scalar,DeviceType> basis2Values = basis2->allocateBasisValues(tensorCubaturePoints, op2);
  basis2->getValues(basis2Values, tensorCubaturePoints, op2);
  basis2Values.setBasis(basis2);
  
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  
  Data<PointScalar,DeviceType> jacobian       = geometry.allocateJacobianData(tensorCubaturePoints, 0, numCells);
  Data<PointScalar,DeviceType> jacobianDividedByJacobianDet = geometry.allocateJacobianData(tensorCubaturePoints, 0, numCells); // jacobianDividedByJacobianDet has same underlying structure as jacobian
  Data<PointScalar,DeviceType> jacobianDet    = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianDetInv = CellTools<DeviceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,DeviceType> jacobianInv    = CellTools<DeviceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,DeviceType> cellMeasures = geometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
  
  const int numJacobianDataPoints = jacobian.getDataExtent(1); // data extent will be 1 for affine, numPoints for other cases
  const int numPoints             = jacobian.extent_int(1); // number of logical points
  
  auto refData = geometry.getJacobianRefData(tensorCubaturePoints);
  
  jacobianAndCellMeasureTimer->start();
  
  geometry.setJacobian(jacobian, tensorCubaturePoints, refData);
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
    
    Kokkos::Array<int,3> extents {numCells,numPoints,spaceDim};
    Kokkos::Array<DataVariationType,3> variationTypes {CONSTANT, CONSTANT, GENERAL};
    
    Data<Scalar,DeviceType> au_data(auView, extents, variationTypes);
    auto uTransform = Data<Scalar,DeviceType>::allocateMatVecResult(transformedBasis1Values.transform(), au_data, true);
    uTransform.storeMatVec(transformedBasis1Values.transform(), au_data, true); // true: transpose basis transform when multiplying
    transformedBasis1Values = Intrepid2::TransformedBasisValues<double, DeviceType>(uTransform, basis1Values);
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
    
    Kokkos::Array<int,3> extents {numCells,numPoints,spaceDim};
    Kokkos::Array<DataVariationType,3> variationTypes {CONSTANT, CONSTANT, GENERAL};
    
    Data<Scalar,DeviceType> av_data(avView, extents, variationTypes);
    auto vTransform = Data<Scalar,DeviceType>::allocateMatVecResult(transformedBasis2Values.transform(), av_data, true);
    vTransform.storeMatVec(transformedBasis2Values.transform(), av_data, true); // true: transpose basis transform when multiplying
    transformedBasis2Values = Intrepid2::TransformedBasisValues<double, DeviceType>(vTransform, basis2Values);
  }
    
  geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
  ExecutionSpace().fence();
  jacobianAndCellMeasureTimer->stop();
    
  return PAMatrix<DeviceType,Scalar>(transformedBasis1Values, cellMeasures, transformedBasis2Values, orientations);
}

//! General assembly for two arbitrary bases and ops that takes advantage of the new structured integration support, including support for sum factorization.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::PAMatrix<DeviceType,Scalar> constructPAMatrix(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry,
                                                         const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1,
                                                         const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2)
{
  Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > nullVectorWeight = Teuchos::null;
  
  return constructPAMatrix<Scalar,BasisFamily,PointScalar,spaceDim,DeviceType>(geometry,
                                                                               polyOrder1, fs1, op1, nullVectorWeight,
                                                                               polyOrder2, fs2, op2, nullVectorWeight);
}

//! constructs a synthetic PAMatrix with custom (F,P) values.  Corresponds to a case with scalar-valued bases in each component dimension with with identity transforms and simple physical-space integration weights.
template<typename DeviceType, class Scalar, class ...ViewProperties>
Intrepid2::PAMatrix<DeviceType,Scalar> syntheticPAMatrix(std::vector<Kokkos::View<Scalar**,ViewProperties...> > basis1RefOps,
                                                         std::vector<Kokkos::View<Scalar**,ViewProperties...> > basis2RefOps,
                                                         Kokkos::View<Scalar**,ViewProperties...> cellMeasures)
{
  using Data                   = Intrepid2::Data                  <Scalar, DeviceType>;
  using DynRankView            = Intrepid2::ScalarView            <Scalar, DeviceType>;
  using TensorData             = Intrepid2::TensorData            <Scalar, DeviceType>;
  using BasisValues            = Intrepid2::BasisValues           <Scalar, DeviceType>;
  using TransformedBasisValues = Intrepid2::TransformedBasisValues<Scalar, DeviceType>;
  
  std::vector<Data> basis1DataComps, basis2DataComps;
  for (const auto &basis1RefOp : basis1RefOps)
  {
    basis1DataComps.push_back(Data(DynRankView(basis1RefOp)));
  }
  for (const auto &basis2RefOp : basis2RefOps)
  {
    basis2DataComps.push_back(Data(DynRankView(basis2RefOp)));
  }
  TensorData basis1TensorData(basis1DataComps);
  TensorData basis2TensorData(basis2DataComps);
  
  BasisValues basisValues1(basis1TensorData);
  BasisValues basisValues2(basis2TensorData);
  
//  template<size_t rank, class ...ViewProperties>
//  Data(Kokkos::View<DataScalar**,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
  
  const int numCells  = cellMeasures.extent_int(0);
  const int numPoints = cellMeasures.extent_int(1);
  
  Kokkos::Array<int,2> extents {numCells,numPoints};
  Kokkos::Array<Intrepid2::DataVariationType,2> variationType {Intrepid2::GENERAL,Intrepid2::GENERAL};
  Data cellMeasuresData(cellMeasures, extents, variationType);
  TensorData cellMeasuresTensorData(cellMeasuresData);
  
  TransformedBasisValues transformedBasisValues1(numCells, basisValues1); // will use identity transform
  TransformedBasisValues transformedBasisValues2(numCells, basisValues2); // will use identity transform
  
  return Intrepid2::PAMatrix<DeviceType,Scalar>(transformedBasisValues1, cellMeasuresTensorData, transformedBasisValues2);
}

#endif /* PAMatrixAssembly_h */
