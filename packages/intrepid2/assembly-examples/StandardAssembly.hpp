// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  StandardAssembly.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 2/23/23.
//

#ifndef Intrepid2_StandardAssembly_hpp
#define Intrepid2_StandardAssembly_hpp

#include "JacobianFlopEstimate.hpp"
#include "Intrepid2_OrientationTools.hpp"

/** \file   StandardAssembly.hpp
    \brief  Locally assembles a matrix of one basis integrated cell-wise against another -- an array of shape (C,F1,F2), with formulation (op1(e1_i), op2(e2_j)), using standard Intrepid2 methods; these do not algorithmically exploit geometric structure.
 */

namespace {
  template<class Scalar, typename DeviceType>
  void transform(Intrepid2::ScalarView<Scalar,DeviceType> &transformedValues, Intrepid2::ScalarView<Scalar,DeviceType> &refValues,
                 const Intrepid2::EFunctionSpace &fs, const Intrepid2::EOperator &op,
                 Intrepid2::ScalarView<Scalar,DeviceType> &jacobian, Intrepid2::ScalarView<Scalar,DeviceType> &jacobianDet,
                 Intrepid2::ScalarView<Scalar,DeviceType> &jacobianInv)
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
            FST::HGRADtransformVALUE(transformedValues, refValues);
            break;
          case OPERATOR_GRAD:
            FST::HGRADtransformGRAD(transformedValues, jacobianInv, refValues);
            break;
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      case EFunctionSpace::FUNCTION_SPACE_HCURL:
      {
        switch (op)
        {
          case OPERATOR_VALUE:
            FST::HCURLtransformVALUE(transformedValues, jacobianInv, refValues);
            break;
          case OPERATOR_CURL:
          {
            const int spaceDim = jacobian.extent_int(2); // jacobian has shape (C,P,D,D)
            if (spaceDim == 2)
            {
              // 2D curl
              FST::HCURLtransformCURL(transformedValues, jacobianDet, refValues);
            }
            else
            {
              // 3D curl
              FST::HCURLtransformCURL(transformedValues, jacobian, jacobianDet, refValues);
            }
          }
            break;
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      case EFunctionSpace::FUNCTION_SPACE_HDIV:
      {
        switch (op)
        {
          case OPERATOR_VALUE:
            FST::HDIVtransformVALUE(transformedValues, jacobian, jacobianDet, refValues);
            break;
          case OPERATOR_DIV:
            FST::HDIVtransformDIV(transformedValues, jacobianDet, refValues);
            break;
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported fs/op combination");
        }
      }
        break;
      case EFunctionSpace::FUNCTION_SPACE_HVOL:
      {
        switch (op)
        {
          case OPERATOR_VALUE:
            FST::HVOLtransformVALUE(transformedValues, jacobianDet, refValues);
            break;
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

//! General assembly for two arbitrary bases and ops that uses the classic, generic Intrepid2 paths.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2>  // spaceDim and spaceDim2 should agree on value (differ on type)
Intrepid2::ScalarView<Scalar,DeviceType> performStandardAssembly(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, int worksetSize,
                                                                 const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight1,
                                                                 const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim2> > vectorWeight2,
                                                                 double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
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
  
  using ViewType = Intrepid2::ScalarView<Scalar,DeviceType>;
  
  using namespace Intrepid2;
  
  using namespace std;
  // dimensions of the returned view are (C,F1,F2)
  
  Intrepid2::ScalarView<Intrepid2::Orientation,DeviceType> orientations("orientations", geometry.numCells() );
  geometry.orientations(orientations, 0, -1);
  
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis1 = getBasis< BasisFamily >(cellTopo, fs1, polyOrder1);
  auto basis2 = getBasis< BasisFamily >(cellTopo, fs2, polyOrder2);
  
  int numFields1 = basis1->getCardinality();
  int numFields2 = basis2->getCardinality();
  int numCells = geometry.numCells();
  
  if (worksetSize > numCells) worksetSize = numCells;
  
  // local stiffness matrices:
  ScalarView<Scalar,DeviceType> cellStiffness("cell stiffness matrices",numCells,numFields1,numFields2);
  
  auto cubature = DefaultCubatureFactory::create<DeviceType>(cellTopo,polyOrder1 + polyOrder2);
  int numPoints = cubature->getNumPoints();
  ScalarView<PointScalar,DeviceType> cubaturePoints("cubature points",numPoints,spaceDim);
  ScalarView<double,DeviceType> cubatureWeights("cubature weights", numPoints);
  
  cubature->getCubature(cubaturePoints, cubatureWeights);
  
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numPoints);
    
  // Allocate some intermediate containers
  ViewType basis1Values = basis1->allocateOutputView(numPoints, op1); // (F1,P[,D])
  ViewType basis2Values = basis2->allocateOutputView(numPoints, op2); // (F2,P[,D])
  
  ViewType orientedValues1, transformedValues1, ultimateValues1;
  ViewType orientedValues2, transformedValues2, ultimateValues2, ultimateWeightedValues2;
  
  int ultimateBasis1Rank, ultimateBasis2Rank;
  if (basis1Values.rank() == 2)
  {
    // the un-transformed values have shape (F,P): scalar values
    // if vector weights supplied, these will increase the rank
    ultimateBasis1Rank = (vectorWeight1 == Teuchos::null) ? 3 : 4; // (C,F,P) or (C,F,P,D)
  }
  else if (basis1Values.rank() == 3)
  {
    // the un-transformed values have shape (F,P,D): vector values
    // if vector weights supplied, these will decrease the rank (we interpret as a dot product)
    ultimateBasis1Rank = (vectorWeight1 == Teuchos::null) ? 4 : 3; // (C,F,P,D) or (C,F,P)
  }
  if (basis2Values.rank() == 2)
  {
    // the un-transformed values have shape (F,P): scalar values
    // if vector weights supplied, these will increase the rank
    ultimateBasis2Rank = (vectorWeight2 == Teuchos::null) ? 3 : 4; // (C,F,P) or (C,F,P,D)
  }
  else if (basis2Values.rank() == 3)
  {
    // the un-transformed values have shape (F,P,D): vector values
    // if vector weights supplied, these will decrease the rank (we interpret as a dot product)
    ultimateBasis2Rank = (vectorWeight2 == Teuchos::null) ? 4 : 3; // (C,F,P,D) or (C,F,P)
  }
  
  INTREPID2_TEST_FOR_EXCEPTION(ultimateBasis1Rank != ultimateBasis2Rank, std::invalid_argument, "basis1 and basis2 must agree on their rank under the respective operators");
  
  if (basis1Values.rank() == 2)
  {
    orientedValues1 = ViewType("oriented values 1", worksetSize, numFields1, numPoints);
    transformedValues1 = ViewType("transformed values 1", worksetSize, numFields1, numPoints);
  }
  else
  {
    orientedValues1 = ViewType("oriented values 1", worksetSize, numFields1, numPoints, spaceDim);
    transformedValues1 = ViewType("transformed values 1", worksetSize, numFields1, numPoints, spaceDim);
  }
  if (basis2Values.rank() == 2)
  {
    orientedValues2 = ViewType("oriented values 2", worksetSize, numFields2, numPoints);
    transformedValues2 = ViewType("transformed values 2", worksetSize, numFields2, numPoints);
  }
  else
  {
    orientedValues2 = ViewType("oriented values 2", worksetSize, numFields2, numPoints, spaceDim);
    transformedValues2 = ViewType("transformed values 2", worksetSize, numFields2, numPoints, spaceDim);
  }
  
  const bool scalarValued = (ultimateBasis1Rank == 3); // (C,F1,P): scalar-valued
  if (scalarValued)
  {
    ultimateValues1 = ViewType("ultimate values 1", worksetSize, numFields1, numPoints);
    ultimateValues2 = ViewType("ultimate values 2", worksetSize, numFields2, numPoints);
    
    ultimateWeightedValues2 = ViewType("ultimate weighted values 2", worksetSize, numFields2, numPoints);
  }
  else // (F1, P, D)
  {
    ultimateValues1 = ViewType("ultimate values 1", worksetSize, numFields1, numPoints, spaceDim);
    ultimateValues2 = ViewType("ultimate values 2", worksetSize, numFields2, numPoints, spaceDim);
    
    ultimateWeightedValues2 = ViewType("ultimate weighted values 2", worksetSize, numFields2, numPoints, spaceDim);
  }
    
  basis1->getValues(basis1Values, cubaturePoints, op1 );
  basis2->getValues(basis2Values, cubaturePoints, op2 );
  
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
  
  ViewType cellMeasures("cell measures", worksetSize, numPoints);
  ViewType jacobianDeterminant("jacobian determinant", worksetSize, numPoints);
  ViewType jacobian("jacobian", worksetSize, numPoints, spaceDim, spaceDim);
  ViewType jacobianInverse("jacobian inverse", worksetSize, numPoints, spaceDim, spaceDim);
  
  // Views used for vector-weighted case:
  ViewType scalarTransformedValues1        ("scalar transformed values 1", worksetSize, numFields1, numPoints);
  ViewType scalarTransformedWeightedValues2("scalar transformed weighted values 2", worksetSize, numFields2, numPoints);

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
      Kokkos::resize(jacobian,            numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianInverse,     numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianDeterminant, numCellsInWorkset, numPoints);
      Kokkos::resize(cellMeasures,        numCellsInWorkset, numPoints);
      Kokkos::resize(jacobianDeterminant, numCellsInWorkset, numPoints);
      
      Kokkos::resize(scalarTransformedValues1,         numCellsInWorkset, numFields1, numPoints);
      Kokkos::resize(scalarTransformedWeightedValues2, numCellsInWorkset, numFields2, numPoints);
      
      if (basis1Values.rank() == 2)
      {
        Kokkos::resize(orientedValues1,    numCellsInWorkset, numFields1, numPoints);
        Kokkos::resize(transformedValues1, numCellsInWorkset, numFields1, numPoints);
      }
      else
      {
        Kokkos::resize(orientedValues1,    numCellsInWorkset, numFields1, numPoints, spaceDim);
        Kokkos::resize(transformedValues1, numCellsInWorkset, numFields1, numPoints, spaceDim);
      }
      if (basis2Values.rank() == 2)
      {
        Kokkos::resize(orientedValues2,    numCellsInWorkset, numFields2, numPoints);
        Kokkos::resize(transformedValues2, numCellsInWorkset, numFields2, numPoints);
      }
      else
      {
        Kokkos::resize(orientedValues2,    numCellsInWorkset, numFields2, numPoints, spaceDim);
        Kokkos::resize(transformedValues2, numCellsInWorkset, numFields2, numPoints, spaceDim);
      }
      
      if (scalarValued)
      {
        Kokkos::resize(ultimateValues1, numCellsInWorkset, numFields1, numPoints);
        Kokkos::resize(ultimateValues2, numCellsInWorkset, numFields2, numPoints);
        
        Kokkos::resize(ultimateWeightedValues2, numCellsInWorkset, numFields2, numPoints);
      }
      else // (F1, P, D)
      {
        ultimateValues1 = ViewType("ultimate values 1", worksetSize, numFields1, numPoints, spaceDim);
        ultimateValues2 = ViewType("ultimate values 2", worksetSize, numFields2, numPoints, spaceDim);
        
        ultimateWeightedValues2 = ViewType("ultimate weighted values 2", worksetSize, numFields2, numPoints, spaceDim);
      }
    }
    jacobianAndCellMeasureTimer->start();
    CellTools::setJacobian(jacobian, cubaturePoints, cellWorkset, cellTopo); // accounted for outside loop, as numCells * flopsPerJacobianPerCell.
    CellTools::setJacobianInv(jacobianInverse, jacobian);
    CellTools::setJacobianDet(jacobianDeterminant, jacobian);
    
    FunctionSpaceTools::computeCellMeasure(cellMeasures, jacobianDeterminant, cubatureWeights);
    ExecutionSpace().fence();
    jacobianAndCellMeasureTimer->stop();
    
    OrientationTools<DeviceType>::modifyBasisByOrientation(orientedValues1, basis1Values, orientationsWorkset, basis1.get());
    OrientationTools<DeviceType>::modifyBasisByOrientation(orientedValues2, basis2Values, orientationsWorkset, basis2.get());
    
    // because structured integration performs transformations within integrate(), to get a fairer comparison here we include the transformation calls.
    fstIntegrateCall->start();
    transform(transformedValues1, orientedValues1, fs1, op1, jacobian, jacobianDeterminant, jacobianInverse);
    transform(transformedValues2, orientedValues2, fs2, op2, jacobian, jacobianDeterminant, jacobianInverse);
        
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    if (vectorWeight1 != Teuchos::null)
    {
      auto uWeight = *vectorWeight1;
      
      auto policy3 = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCellsInWorkset,numFields1,numPoints});
      if (transformedValues1.rank() == 4)
      {
        Kokkos::parallel_for("compute ultimateValues1", policy3,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
        {
          Scalar u_result = 0;
          for (int d=0; d<spaceDim; d++)
          {
            u_result += uWeight[d] * transformedValues1(cellOrdinal,fieldOrdinal,pointOrdinal,d);
          }
          ultimateValues1(cellOrdinal,fieldOrdinal,pointOrdinal) = u_result;
        });
      }
      else // transformedValues1.rank() == 3
      {
        Kokkos::parallel_for("compute ultimateValues1", policy3,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
        {
          const Scalar & value1 = transformedValues1(cellOrdinal,fieldOrdinal,pointOrdinal);
          for (int d=0; d<spaceDim; d++)
          {
            ultimateValues1(cellOrdinal,fieldOrdinal,pointOrdinal, d) = value1 * uWeight[d];
          }
        });
      }
    }
    else
    {
      ultimateValues1 = transformedValues1;
    }
    
    if (vectorWeight2 != Teuchos::null)
    {
      auto vWeight = *vectorWeight2;
      
      auto policy3 = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCellsInWorkset,numFields2,numPoints});
      if (transformedValues2.rank() == 4)
      {
        Kokkos::parallel_for("compute ultimateValues2", policy3,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
        {
          Scalar v_result = 0;
          for (int d=0; d<spaceDim; d++)
          {
            v_result += vWeight[d] * transformedValues2(cellOrdinal,fieldOrdinal,pointOrdinal,d);
          }
          ultimateValues2(cellOrdinal,fieldOrdinal,pointOrdinal) = v_result;
        });
      }
      else // transformedValues2.rank() == 3
      {
        Kokkos::parallel_for("compute ultimateValues2", policy3,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal)
        {
          const Scalar & value2 = transformedValues2(cellOrdinal,fieldOrdinal,pointOrdinal);
          for (int d=0; d<spaceDim; d++)
          {
            ultimateValues2(cellOrdinal,fieldOrdinal,pointOrdinal, d) = value2 * vWeight[d];
          }
        });
      }
    }
    else
    {
      ultimateValues2 = transformedValues2;
    }
    
    FunctionSpaceTools::multiplyMeasure(ultimateWeightedValues2, cellMeasures, ultimateValues2);
    transformIntegrateFlopCount += ultimateValues2.size(); // multiply each entry of ultimateValues2: one flop for each.
    
    FunctionSpaceTools::integrate(cellStiffnessSubview, ultimateValues1, ultimateWeightedValues2);
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields1+numFields2) * double(numPoints) * double(spaceDim) * (spaceDim - 1) * 2.0; // 2: one multiply, one add per (P,D) entry in the contraction.
    
    ExecutionSpace().fence();
    fstIntegrateCall->stop();
    
    // TODO: correct these flop count sums based on whether they are scalar or vector-valued (and -- above -- whether we do HGRAD_transform_VALUE, which is just a copy operation, no flops).
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields1) * double(numFields2) * double(numPoints * spaceDim * 2 - 1); // 2: one multiply, one add per (P,D) entry in the contraction.
    
    cellOffset += worksetSize;
  }
//  std::cout << "standard integration, approximateFlopCount: " << approximateFlopCount << std::endl;
  return cellStiffness;
}

//! General assembly for two arbitrary bases and ops that uses the classic, generic Intrepid2 paths.
template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStandardAssembly(Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, int worksetSize,
                                                                 const int &polyOrder1, const Intrepid2::EFunctionSpace &fs1, const Intrepid2::EOperator &op1,
                                                                 const int &polyOrder2, const Intrepid2::EFunctionSpace &fs2, const Intrepid2::EOperator &op2,
                                                                 double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > nullVectorWeight = Teuchos::null;
  
  return performStandardAssembly<Scalar,BasisFamily,PointScalar,spaceDim,DeviceType>(geometry, worksetSize,
                                                                                     polyOrder1, fs1, op1, nullVectorWeight,
                                                                                     polyOrder2, fs2, op2, nullVectorWeight,
                                                                                     transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
}
#endif /* StandardAssembly_hpp */
