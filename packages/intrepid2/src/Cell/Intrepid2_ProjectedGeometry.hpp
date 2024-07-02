// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectedGeometry.hpp
    \brief  Allows generation of geometry degrees of freedom based on a provided map from straight-edged mesh domain to curvilinear mesh domain.

    \author Nathan V. Roberts
*/
#ifndef Intrepid2_ProjectedGeometry_h
#define Intrepid2_ProjectedGeometry_h

#include "Intrepid2_ScalarView.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"

#include "Intrepid2_TestUtils.hpp"

namespace Intrepid2
{
/** \class Intrepid2::ProjectedGeometry
    \brief Allows generation of geometry degrees of freedom based on a provided map from straight-edged mesh domain to curvilinear mesh domain.
*/
  template<int spaceDim, typename PointScalar, typename DeviceType>
  class ProjectedGeometry
  {
  public:
    using ViewType      = ScalarView<      PointScalar, DeviceType>;
    using ConstViewType = ScalarView<const PointScalar, DeviceType>;
    using BasisPtr      = Teuchos::RCP< Basis<DeviceType,PointScalar,PointScalar> >;
    
    /** \brief Generate geometry degrees of freedom based on a provided map from straight-edged mesh domain to curvilinear mesh domain.
       \param [out] projectedBasisNodes - the projected geometry degrees of freedom
       \param [in] targetHGradBasis - the (potentially higher-order) H^1 basis to use to define the basis nodes
       \param [in] flatCellGeometry - the low-order geometry on which the exact geometry functions are defined
       \param [in] exactGeometry - a function defined on the domain described by flatCellGeometry, with values corresponding to the (higher-order) physical geometry.
       \param [in] exactGeometryGradient - the gradient of exactGeometry.
        For space dimension D, exactGeometry is a functor, accepting a Kokkos::Array<PointScalar,D> as first argument, and a component ordinal representing the coordinate dimension d (0 <= d < D).  The return is a PointScalar.  exactGeometryGradient is similar, except that it also accepts a second dimension component d2, indicating the direction of the derivative to be evaluated.
     \see Intrepid2_ProjectedGeometryExamples.hpp for sample implementations of exactGeometry and exactGeometryGradient.
    */
    template<class ExactGeometry, class ExactGeometryGradient>
    static void projectOntoHGRADBasis(ViewType projectedBasisNodes, BasisPtr targetHGradBasis, CellGeometry<PointScalar,spaceDim,DeviceType> flatCellGeometry,
                                      const ExactGeometry &exactGeometry, const ExactGeometryGradient &exactGeometryGradient)
    {
      const ordinal_type numCells = flatCellGeometry.extent_int(0); // (C,N,D)
      
      INTREPID2_TEST_FOR_EXCEPTION(spaceDim != targetHGradBasis->getBaseCellTopology().getDimension(), std::invalid_argument, "spaceDim must match the cell topology on which target basis is defined");
      INTREPID2_TEST_FOR_EXCEPTION(projectedBasisNodes.rank() != 3, std::invalid_argument, "projectedBasisNodes must have shape (C,F,D)");
      INTREPID2_TEST_FOR_EXCEPTION(projectedBasisNodes.extent_int(0) != numCells, std::invalid_argument, "cell counts must match in projectedBasisNodes and cellNodesToMap");
      INTREPID2_TEST_FOR_EXCEPTION(projectedBasisNodes.extent_int(1) != targetHGradBasis->getCardinality(), std::invalid_argument, "projectedBasisNodes must have shape (C,F,D)");
      INTREPID2_TEST_FOR_EXCEPTION(projectedBasisNodes.extent_int(2) != spaceDim, std::invalid_argument, "projectedBasisNodes must have shape (C,F,D)");
      
      using ExecutionSpace = typename DeviceType::execution_space;
      using ProjectionTools  = ProjectionTools<DeviceType>; 
      using ProjectionStruct = ProjectionStruct<DeviceType,PointScalar>;
      
      ProjectionStruct projectionStruct;
      ordinal_type targetQuadratureDegree(targetHGradBasis->getDegree()), targetDerivativeQuadratureDegree(targetHGradBasis->getDegree());
      projectionStruct.createHGradProjectionStruct(targetHGradBasis, targetQuadratureDegree, targetDerivativeQuadratureDegree);
    
      auto evaluationPointsRefSpace = projectionStruct.getAllEvalPoints();
      auto evaluationGradPointsRefSpace = projectionStruct.getAllDerivEvalPoints();
      const ordinal_type numPoints     = evaluationPointsRefSpace.extent(0);
      const ordinal_type numGradPoints = evaluationGradPointsRefSpace.extent(0);

      auto elementOrientations = flatCellGeometry.getOrientations();
      
//      printFunctor1(elementOrientations, std::cout);
      
      // the evaluation points are all still in reference space; map to physical space:
      ViewType evaluationPoints    ("ProjectedGeometry evaluation points (value)",    numCells, numPoints,     spaceDim);
      ViewType evaluationGradPoints("ProjectedGeometry evaluation points (gradient)", numCells, numGradPoints, spaceDim);
  
      using CellTools = CellTools<DeviceType>;
      BasisPtr hgradLinearBasisForFlatGeometry = flatCellGeometry.basisForNodes();
      if (numPoints > 0)
      {
        CellTools::mapToPhysicalFrame(evaluationPoints, evaluationPointsRefSpace, flatCellGeometry, hgradLinearBasisForFlatGeometry);
      }
      if (numGradPoints > 0)
      {
        CellTools::mapToPhysicalFrame(evaluationGradPoints, evaluationGradPointsRefSpace, flatCellGeometry, hgradLinearBasisForFlatGeometry);
      }
      
      auto refData = flatCellGeometry.getJacobianRefData(evaluationGradPoints);
      
      // evaluate, transform, and project in each component
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},  {numCells,numPoints});
      auto gradPolicy  = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numGradPoints,spaceDim});
      
      ViewType evaluationValues    ("exact geometry values",    numCells, numPoints);
      ViewType evaluationGradients ("exact geometry gradients", numCells, numGradPoints, spaceDim);
      
//      printView(evaluationPoints, std::cout, "evaluationPoints");
      
      for (int comp=0; comp<spaceDim; comp++)
      {
        Kokkos::parallel_for("evaluate geometry function for projection", policy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal) {
          Kokkos::Array<PointScalar,spaceDim> point;
          for (int d=0; d<spaceDim; d++)
          {
            point[d] = evaluationPoints(cellOrdinal,pointOrdinal,d);
          }
          evaluationValues(cellOrdinal,pointOrdinal) = exactGeometry(point,comp);
        });
        
//        printView(evaluationValues, std::cout, "evaluationValues");
        
        // projection occurs in ref space, so we need to apply inverse of the pullback
        // HGRADtransformVALUE is identity, so evaluationValues above is correct
        // HGRADtransformGRAD  is multiplication by inverse of Jacobian, so here we want to multiply by Jacobian
        
        auto gradPointsJacobians = flatCellGeometry.allocateJacobianData(evaluationGradPoints);
        flatCellGeometry.setJacobian(gradPointsJacobians,evaluationGradPoints,refData);
        
        Kokkos::parallel_for("evaluate geometry gradients for projection", gradPolicy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal, const int &d2) {
          Kokkos::Array<PointScalar,spaceDim> point;
          for (int d=0; d<spaceDim; d++)
          {
            point[d] = evaluationGradPoints(cellOrdinal,pointOrdinal,d);
          }
          evaluationGradients(cellOrdinal,pointOrdinal,d2) = exactGeometryGradient(point,comp,d2);
        });
        
        // apply Jacobian
        Data<PointScalar,DeviceType> gradientData(evaluationGradients);
        auto transformedGradientData = Data<PointScalar,DeviceType>::allocateMatVecResult(gradPointsJacobians,gradientData);
        
        transformedGradientData.storeMatVec(gradPointsJacobians,gradientData);
        
        auto projectedBasisNodesForComp = Kokkos::subview(projectedBasisNodes,Kokkos::ALL(),Kokkos::ALL(),comp);
        
        ProjectionTools::getHGradBasisCoeffs(projectedBasisNodesForComp,
                                             evaluationValues,
                                             transformedGradientData.getUnderlyingView(),
                                             elementOrientations,
                                             targetHGradBasis.get(),
                                             &projectionStruct);
      }
    }
  };
}

#endif /* Intrepid2_ProjectedGeometry_h */
