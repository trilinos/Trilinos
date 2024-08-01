// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking accuracy of interpolation-based projections for quadrilateral elements

    The test considers a uniform and structured quadrilateral mesh of the square [-1,1]^2, formed by
    N^2 quads, and checks the accuracy of the HGRAD, HCURL, HDIV, HVOL projections of analytic
    target functions for for Hierarchical and Nodal basis functions as N increases.
    The accuracy is computed in the H^1, H^{curl}, H^{div} and L^2 norms respectively. The optimal
    order of convergence equates the basis degree.

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "struct_mesh_utils.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <array>
#include <set>
#include <random>
#include <algorithm>

namespace Intrepid2 {

namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceType>
int ConvergenceQuad(const bool verbose) {

  using  DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;
  using  DynRankViewInt = Kokkos::DynRankView<ordinal_type,DeviceType>;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  using ExecSpaceType = typename DeviceType::execution_space;
  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

  *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
  *outStream << "\n";

  int errorFlag = 0;
  const ValueType relTol = 1e-4;

  struct Fun {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y) {
      return sin(x*2)*sin(y*2)+sin(x*y*8);
    }
  };

  struct GradFun {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const int comp) {
      switch (comp) {
      case 0:
        return cos(x*2)*sin(y*2)*2+cos(x*y*8)*y*8;
      case 1:
        return sin(x*2)*cos(y*2)*2+cos(x*y*8)*x*8;
      default:
        return 0;
      }
    }
  };

  struct FunCurl {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const int comp) {
      ValueType f0 = sin(x*2)*sin(y*2)+sin(x*y*8);
      ValueType f1 = cos(x*2)*cos(y*2);
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0;
      case 1:
        return f1;
      default:
        return 0;
      }
    }
  };

  struct CurlFunCurl {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y) {
      ValueType gf0[2] = {cos(x*2)*sin(y*2)*2+cos(x*y*8)*y*8, sin(x*2)*cos(y*2)*2+cos(x*y*8)*x*8};
      ValueType gf1[2] = {-sin(x*2)*cos(y*2)*2, -cos(x*2)*sin(y*2)*2};
      return gf1[0] - gf0[1];
    }
  };

  struct FunDiv {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const int comp) {
      ValueType f0 = sin(x*2)*sin(y*2)+sin(x*y*8);
      ValueType f1 = cos(x*2)*cos(y*2);
      //fun = f + a x
      switch (comp) {
      case 0:
        return f0;
      case 1:
        return f1;
      default:
        return 0;
      }
    }
  };

  struct DivFunDiv {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y) {
      ValueType gxf0 = cos(x*2)*sin(y*2)*2+cos(x*y*8)*y*8;
      ValueType gyf1 = -cos(x*2)*sin(y*2)*2;
      return gxf0+gyf1;
    }
  };

  using ct = CellTools<DeviceType>;
  using ots = OrientationTools<DeviceType>;
  using pts = ProjectionTools<DeviceType>;
  using ProjStruct = ProjectionStruct<DeviceType,ValueType>;
  using rst = RealSpaceTools<DeviceType>;
  using fst = FunctionSpaceTools<DeviceType>;

  constexpr ordinal_type dim = 2;
  const ordinal_type basisDegree = 3;
  ordinal_type cub_degree = 7;



  // ************************************ GET INPUTS **************************************


  int NX = 2;
  constexpr int numRefinements = 2;

  // Expected values of the projection errors in H1, HCurl, HDiv and L2 norms for HGRAD, HDIV, HCURL and HVOL elements respectively.
  // These values have been computed running the code with numRefinements=4 and the convergence rates are close to the optimal ones.
  // Note that these values are independent of the basis choice (Hierarchical or Nodal) as long as they generate the same functional space.
  // We currently only test two mesh refinements to make the test run faster, so this is used as a regression test rather than
  // a convergence test, but the test can be use for verifying optimal accuracy as well.
  ValueType hgradNorm[numRefinements], hcurlNorm[numRefinements], hdivNorm[numRefinements], hvolNorm[numRefinements];

  ValueType hgrad_errors[4] = {4.06298, 1.12792, 0.127858, 0.0164394};
  ValueType hcurl_errors[4] = {3.06114, 0.976496, 0.114555, 0.0150085};
  ValueType hdiv_errors[4] = {4.11097, 0.997027, 0.115616, 0.0150505};
  ValueType hvol_errors[4] = {0.896361, 0.104481, 0.0189055, 0.00239322};
  
  ValueType hgrad_errors_L2[4] = {4.04405, 1.13086, 0.128483, 0.0165344};
  ValueType hcurl_errors_L2[4] = {3.13114, 0.980243, 0.114986, 0.0150622};
  ValueType hdiv_errors_L2[4] = {3.83504, 1.00134, 0.117383, 0.0152593};
  ValueType hvol_errors_L2[4] = {0.896361, 0.104481, 0.0189055, 0.00239322};

  for(int iter= 0; iter<numRefinements; iter++, NX *= 2) {
    int NY            = NX;
    bool randomMesh    = 0;  // 1 if mesh randomizer is to be used 0 if not

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base quadrilateral
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Get dimensions
    int numNodesPerElem = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    DynRankView nodeCoord;
    DynRankViewInt elemNodes;
    createStructMesh(nodeCoord, elemNodes, cellTopo, NX, NY, -1, randomMesh, *outStream);

    //computing vertices coords
    ordinal_type numElems = elemNodes.extent(0);
    DynRankView ConstructWithLabel(physVertexes, numElems, numNodesPerElem, dim);
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
    KOKKOS_LAMBDA (const int &i) {
      for(ordinal_type j=0; j<numNodesPerElem; ++j)
        for(ordinal_type k=0; k<dim; ++k)
          physVertexes(i,j,k) = nodeCoord(elemNodes(i,j),k);
    });
    ExecSpaceType().fence();

    DefaultCubatureFactory cub_factory;
    auto cell_cub = cub_factory.create<DeviceType, ValueType, ValueType>(cellTopo.getBaseKey(), cub_degree);
    ordinal_type numRefCoords = cell_cub->getNumPoints();
    DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
    DynRankView ConstructWithLabel(weights, numRefCoords);
    cell_cub->getCubature(refPoints, weights);

    using basisType = Basis<DeviceType,ValueType,ValueType>;
    using CG_NBasis = NodalBasisFamily<DeviceType,ValueType,ValueType>;
    using CG_HBasis = HierarchicalBasisFamily<DeviceType,ValueType,ValueType>;
    //using CG_DNBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;


    std::vector<bool> useL2Proj;
    useL2Proj.push_back(false);
    useL2Proj.push_back(true); // use L2 projection for all the bases

    *outStream
    << "===============================================================================\n"
    << "|                                                                             |\n"
    << "|                 Test 1 (Convergence - HGRAD)                                |\n"
    << "|                                                                             |\n"
    << "===============================================================================\n";

    try {
      //compute reference points
      Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(basisDegree,POINTTYPE_WARPBLEND); //used only for computing reference points

      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, cellTopo);

      for (auto useL2Projection:useL2Proj) { //

      std::vector<basisType*> basis_set;
      basis_set.push_back(new typename  CG_NBasis::HGRAD_QUAD(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HGRAD_QUAD(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(linearBasisValuesAtRefCoords, numNodesPerElem, numRefCoords);
          linearBasis.getValues(linearBasisValuesAtRefCoords, refPoints);
          ExecSpaceType().fence();
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type d=0; d<dim; ++d)
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*linearBasisValuesAtRefCoords(k,j);
          });
          ExecSpaceType().fence();
        }

        DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
        DynRankView ConstructWithLabel(funGradAtPhysRefCoords, numElems, numRefCoords, dim);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i) {
          Fun fun;
          GradFun gradFun;
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
            for(ordinal_type d=0; d<dim; ++d)
              funGradAtPhysRefCoords(i,j,d) = gradFun(physRefCoords(i,j,0), physRefCoords(i,j,1),d);
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HGRAD
        DynRankView ConstructWithLabel(basisCoeffsHGrad, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree());

          ProjStruct projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHGradProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }
          
          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationGradPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numGradPoints = evaluationGradPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints);
          DynRankView ConstructWithLabel(targetGradAtEvalPoints, numElems, numGradPoints, dim);

          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalGradPoints, numElems, numGradPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalGradPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalGradPoint = Kokkos::subview(linearBasisValuesAtEvalGradPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numGradPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationGradPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalGradPoint, evalGradPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalGradPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalGradPoint(k);
              }
            });
            ExecSpaceType().fence();
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numElems, numGradPoints, dim, dim);
          if(numGradPoints>0)
            ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, cellTopo);

          Kokkos::deep_copy(targetGradAtEvalPoints,0.);
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &ic) {
            Fun fun;
            GradFun gradFun;
            for(int i=0;i<numPoints;i++) {
              targetAtEvalPoints(ic,i) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
            }
            for(int i=0;i<numGradPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetGradAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*gradFun(physEvalGradPoints(ic,i,0), physEvalGradPoints(ic,i,1), d);
            }
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHGrad,
                targetAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHGradBasisCoeffs(basisCoeffsHGrad,
                targetAtEvalPoints,
                targetGradAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          }
        }

        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView basisValuesAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords);

        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values
        deep_copy(transformedBasisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsOriented);

        DynRankView ConstructWithLabel(basisGradsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisGradsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView basisGradsAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords, dim);

        DynRankView ConstructWithLabel(basisGradsAtRefCoords, basisCardinality, numRefCoords, dim);
        basis.getValues(basisGradsAtRefCoords, refPoints,OPERATOR_GRAD);
        rst::clone(basisGradsAtRefCoordsCells,basisGradsAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisGradsAtRefCoordsOriented,
            basisGradsAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values to the reference element (pullback)
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_inv, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numElems, numRefCoords);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, cellTopo);
        ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
        ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
        fst::HGRADtransformGRAD(transformedBasisGradsAtRefCoordsOriented,
            jacobianAtRefCoords_inv,
            basisGradsAtRefCoordsOriented);

        DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords);
        DynRankView ConstructWithLabel(funGradAtRefCoordsOriented, numElems, numRefCoords,dim);

        //compute error of projection in H1 norm
        ValueType norm2(0);
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i, double &norm2Update) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k) {
              projectedFunAtRefCoords(i,j) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
              for (ordinal_type d=0; d<dim; ++d)
                funGradAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisGradsAtRefCoordsOriented(i,k,j,d);
            }
            norm2Update += (funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j))*
                (funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j))*
                weights(j)*jacobianAtRefCoords_det(i,j);
            for (ordinal_type d=0; d<dim; ++d)
              norm2Update += (funGradAtPhysRefCoords(i,j,d) - funGradAtRefCoordsOriented(i,j,d))*
                (funGradAtPhysRefCoords(i,j,d) - funGradAtRefCoordsOriented(i,j,d))*
                weights(j)*jacobianAtRefCoords_det(i,j);
          }
        }, norm2);

        ExecSpaceType().fence();

        hgradNorm[iter] =  std::sqrt(norm2);
        auto expected_error = useL2Projection ? hgrad_errors_L2[iter] : hgrad_errors[iter];
        if(std::abs(hgradNorm[iter]-expected_error)/expected_error > relTol){
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "For N = " << NX << ", computed error (" << hgradNorm[iter] << ") is different than expected one (" << expected_error << ")";
          *outStream << std::endl;
        }
        delete basisPtr;
      }
      if(useL2Projection)
        *outStream << "HGRAD Error (L2 Projection): " << hgradNorm[iter] <<std::endl;
      else
        *outStream << "HGRAD Error (HGrad Projection): " << hgradNorm[iter] <<std::endl;
      }
    } catch (std::exception &err) {
      std::cout << " Exeption\n";
      *outStream << err.what() << "\n\n";
      errorFlag = -1000;
    }

    *outStream
    << "===============================================================================\n"
    << "|                                                                             |\n"
    << "|                 Test 2 (Convergence - HCURL)                                |\n"
    << "|                                                                             |\n"
    << "===============================================================================\n";


    try {
      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, cellTopo);

      for (auto useL2Projection:useL2Proj) { 

      std::vector<basisType*> basis_set;
      basis_set.push_back(new typename  CG_NBasis::HCURL_QUAD(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HCURL_QUAD(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(linearBasisValuesAtRefCoords, cellTopo.getNodeCount(), numRefCoords);
          linearBasis.getValues(linearBasisValuesAtRefCoords, refPoints);
          ExecSpaceType().fence();
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type d=0; d<dim; ++d)
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*linearBasisValuesAtRefCoords(k,j);
          });
          ExecSpaceType().fence();
        }
        //check function reproducibility

        DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(funCurlAtPhysRefCoords, numElems, numRefCoords);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i) {
          FunCurl fun;
          CurlFunCurl curlFun;
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<dim; ++k)
              funAtRefCoords(i,j,k) =         fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
            funCurlAtPhysRefCoords(i,j) = curlFun(physRefCoords(i,j,0), physRefCoords(i,j,1));
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HCURL
        DynRankView ConstructWithLabel(basisCoeffsHCurl, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(cub_degree),targetDerivCubDegree(cub_degree-1);

          ProjStruct projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHCurlProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationCurlPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numCurlPoints = evaluationCurlPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numElems, numCurlPoints);


          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalCurlPoints, numElems, numCurlPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalCurlPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalCurlPoint = Kokkos::subview(linearBasisValuesAtEvalCurlPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numCurlPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationCurlPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalCurlPoint, evalGradPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalCurlPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalCurlPoint(k);
              }
            });
            ExecSpaceType().fence();
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, cellTopo);

          DynRankView ConstructWithLabel(jacobianCurl_det, numElems, numCurlPoints);
          if(numCurlPoints>0){
            DynRankView ConstructWithLabel(jacobianCurl, numElems, numCurlPoints, dim, dim);
            ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, cellTopo);
            ct::setJacobianDet (jacobianCurl_det, jacobianCurl);
          }

          Kokkos::deep_copy(targetCurlAtEvalPoints,0.);
          Kokkos::deep_copy(targetAtEvalPoints,0.);
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &ic) {
            FunCurl fun;
            CurlFunCurl curlFun;
            for(int i=0;i<numPoints;i++) {
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
            }
            for(int i=0;i<numCurlPoints;i++)
              targetCurlAtEvalPoints(ic,i) += jacobianCurl_det(ic,i)*curlFun(physEvalCurlPoints(ic,i,0), physEvalCurlPoints(ic,i,1));
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHCurl,
                targetAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
                targetAtEvalPoints,
                targetCurlAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          }
        }

        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView basisValuesAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords, dim);


        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values to the reference element (pullback)
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_inv, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numElems, numRefCoords);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, cellTopo);
        ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
        ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
        fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
            jacobianAtRefCoords_inv,
            basisValuesAtRefCoordsOriented);


        DynRankView ConstructWithLabel(basisCurlsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView ConstructWithLabel(transformedBasisCurlsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView basisCurlsAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords);

        DynRankView ConstructWithLabel(basisCurlsAtRefCoords, basisCardinality, numRefCoords);
        basis.getValues(basisCurlsAtRefCoords, refPoints,OPERATOR_CURL);
        rst::clone(basisCurlsAtRefCoordsCells,basisCurlsAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisCurlsAtRefCoordsOriented,
            basisCurlsAtRefCoordsCells,
            elemOrts,
            &basis);


        fst::HVOLtransformVALUE(transformedBasisCurlsAtRefCoordsOriented,
            jacobianAtRefCoords_det,
            basisCurlsAtRefCoordsOriented);

        DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(funCurlAtRefCoordsOriented, numElems, numRefCoords);

        //compute error of projection in HCURL norm
        ValueType norm2(0);
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i, double &norm2Update) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k) {
              for(ordinal_type d=0; d<dim; ++d)
                projectedFunAtRefCoords(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
              funCurlAtRefCoordsOriented(i,j) += basisCoeffsHCurl(i,k)*transformedBasisCurlsAtRefCoordsOriented(i,k,j);
            }
            
            for(ordinal_type d=0; d<dim; ++d)
              norm2Update += (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                weights(j)*jacobianAtRefCoords_det(i,j);

            norm2Update += (funCurlAtPhysRefCoords(i,j) - funCurlAtRefCoordsOriented(i,j))*
              (funCurlAtPhysRefCoords(i,j) - funCurlAtRefCoordsOriented(i,j))*
              weights(j)*jacobianAtRefCoords_det(i,j);
          }
        },norm2);
        ExecSpaceType().fence();

        hcurlNorm[iter] =  std::sqrt(norm2);
        auto expected_error = useL2Projection ? hcurl_errors_L2[iter] : hcurl_errors[iter];
        if(std::abs(hcurlNorm[iter]-expected_error)/expected_error > relTol){
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "For N = " << NX << ", computed error (" << hcurlNorm[iter] << ") is different than expected one (" << expected_error << ")";
          *outStream << std::endl;
        }
        delete basisPtr;
      }
      if(useL2Projection)
        *outStream << "HCURL Error (L2 Projection): " << hcurlNorm[iter] <<std::endl;
      else
        *outStream << "HCURL Error (HCurl Projection): " << hcurlNorm[iter] <<std::endl;
      }
    } catch (std::exception &err) {
      std::cout << " Exeption\n";
      *outStream << err.what() << "\n\n";
      errorFlag = -1000;
    }


    *outStream
    << "===============================================================================\n"
    << "|                                                                             |\n"
    << "|                 Test 3 (Convergence - HDIV)                                 |\n"
    << "|                                                                             |\n"
    << "===============================================================================\n";


    try {
      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, cellTopo);

      for (auto useL2Projection:useL2Proj) { //

      std::vector<basisType*> basis_set;
      basis_set.push_back(new typename  CG_NBasis::HDIV_QUAD(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HDIV_QUAD(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numElems, basisCardinality, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(linearBasisValuesAtRefCoords, numNodesPerElem, numRefCoords);
          linearBasis.getValues(linearBasisValuesAtRefCoords, refPoints);
          ExecSpaceType().fence();
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type d=0; d<dim; ++d)
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*linearBasisValuesAtRefCoords(k,j);
          });
          ExecSpaceType().fence();
        }

        DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(funDivAtPhysRefCoords, numElems, numRefCoords);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i) {
          FunDiv fun;
          DivFunDiv funDiv;
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            funDivAtPhysRefCoords(i,j) = funDiv(physRefCoords(i,j,0), physRefCoords(i,j,1));
            for(ordinal_type k=0; k<dim; ++k)
              funAtRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree()-1);

          ProjStruct projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHDivProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationDivPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numDivPoints = evaluationDivPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numElems, numDivPoints);


          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalDivPoints, numElems, numDivPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalDivPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalDivPoint = Kokkos::subview(linearBasisValuesAtEvalDivPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numDivPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationDivPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalDivPoint, evalGradPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalDivPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalDivPoint(k);
              }
            });
            ExecSpaceType().fence();
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
          DynRankView ConstructWithLabel(jacobian_inv, numElems, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, cellTopo);
          ct::setJacobianDet (jacobian_det, jacobian);
          ct::setJacobianInv (jacobian_inv, jacobian);

          DynRankView ConstructWithLabel(jacobianDiv_det, numElems, numDivPoints);
          if(numDivPoints>0){
            DynRankView ConstructWithLabel(jacobianDiv, numElems, numDivPoints, dim, dim);
            ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, cellTopo);
            ct::setJacobianDet (jacobianDiv_det, jacobianDiv);
          }

          Kokkos::deep_copy(targetDivAtEvalPoints,0.);
          Kokkos::deep_copy(targetAtEvalPoints,0.);
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &ic) {
            FunDiv fun;
            DivFunDiv divFun;
            for(int i=0;i<numPoints;i++) {
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
            }
            for(int i=0;i<numDivPoints;i++) {
              targetDivAtEvalPoints(ic,i) += jacobianDiv_det(ic,i)*divFun(physEvalDivPoints(ic,i,0), physEvalDivPoints(ic,i,1));
            }
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHDiv,
                targetAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHDivBasisCoeffs(basisCoeffsHDiv,
                targetAtEvalPoints,
                targetDivAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          }
        }

        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
        DynRankView basisValuesAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords, dim);

        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values to the reference element (pullback)
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numElems, numRefCoords);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, cellTopo);
        ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
        fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
            jacobianAtRefCoords,
            jacobianAtRefCoords_det,
            basisValuesAtRefCoordsOriented);

        DynRankView ConstructWithLabel(basisDivsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView ConstructWithLabel(transformedBasisDivsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView basisDivsAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords);

        DynRankView ConstructWithLabel(basisDivsAtRefCoords, basisCardinality, numRefCoords);
        basis.getValues(basisDivsAtRefCoords, refPoints,OPERATOR_DIV);
        rst::clone(basisDivsAtRefCoordsCells,basisDivsAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisDivsAtRefCoordsOriented,
            basisDivsAtRefCoordsCells,
            elemOrts,
            &basis);


        fst::HDIVtransformDIV(transformedBasisDivsAtRefCoordsOriented,
            jacobianAtRefCoords_det,
            basisDivsAtRefCoordsOriented);


        DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(funDivAtRefCoordsOriented, numElems, numRefCoords);

        //compute error of projection in HDIV norm
        ValueType norm2(0);
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i, double &norm2Update) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k) {
              for(ordinal_type d=0; d<dim; ++d)
                projectedFunAtRefCoords(i,j,d) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
              funDivAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisDivsAtRefCoordsOriented(i,k,j);
            }

            for(ordinal_type d=0; d<dim; ++d) {
              norm2Update += (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                  (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                  weights(j)*jacobianAtRefCoords_det(i,j);
            }
            norm2Update += (funDivAtPhysRefCoords(i,j) - funDivAtRefCoordsOriented(i,j))*
                (funDivAtPhysRefCoords(i,j) - funDivAtRefCoordsOriented(i,j))*
                weights(j)*jacobianAtRefCoords_det(i,j);
          }
        },norm2);
        ExecSpaceType().fence();
        hdivNorm[iter] = std::sqrt(norm2);
        auto expected_error = useL2Projection ? hdiv_errors_L2[iter] : hdiv_errors[iter];
        if(std::abs(hdivNorm[iter]-expected_error)/expected_error > relTol){
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "For N = " << NX << ", computed error (" << hdivNorm[iter] << ") is different than expected one (" << expected_error << ")";
          *outStream << std::endl;
        }
        delete basisPtr;
      }
      if(useL2Projection)
        *outStream << "HDIV Error (L2 Projection): " << hdivNorm[iter] <<std::endl;
      else
        *outStream << "HDIV Error (HDiv Projection): " << hdivNorm[iter] <<std::endl;
      }
    } catch (std::exception &err) {
      std::cout << " Exeption\n";
      *outStream << err.what() << "\n\n";
      errorFlag = -1000;
    }



    *outStream
    << "===============================================================================\n"
    << "|                                                                             |\n"
    << "|                 Test 4 (Convergence - HVOL)                                 |\n"
    << "|                                                                             |\n"
    << "===============================================================================\n";


    try {
      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, cellTopo);

      for (auto useL2Projection:useL2Proj) { //
      std::vector<basisType*> basis_set;
      basis_set.push_back(new typename  CG_NBasis::HVOL_QUAD(basisDegree-1));
      basis_set.push_back(new typename  CG_HBasis::HVOL_QUAD(basisDegree-1));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(linearBasisValuesAtRefCoords, numNodesPerElem, numRefCoords);
          linearBasis.getValues(linearBasisValuesAtRefCoords, refPoints);
          ExecSpaceType().fence();
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type d=0; d<dim; ++d)
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*linearBasisValuesAtRefCoords(k,j);
          });
          ExecSpaceType().fence();
        }

        //check function reproducibility
        DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i) {
          Fun fun;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HVOL
        DynRankView ConstructWithLabel(basisCoeffsHVol, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree());

          ProjStruct projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHVolProjectionStruct(&basis, targetCubDegree);
          }

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints);


          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
            ExecSpaceType().fence();
          }

          //transform the target function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, cellTopo);
          ct::setJacobianDet (jacobian_det, jacobian);

          Kokkos::deep_copy(targetAtEvalPoints,0.);
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
          KOKKOS_LAMBDA (const int &ic) {
            Fun fun;
            for(int i=0;i<numPoints;i++)
              targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
          });
          ExecSpaceType().fence();
          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHVol,
                targetAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHVolBasisCoeffs(basisCoeffsHVol,
                targetAtEvalPoints,
                elemOrts,
                &basis,
                &projStruct);
          }
        }


        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numElems, basisCardinality, numRefCoords);
        DynRankView basisValuesAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords);

        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values to the reference element (pullback)
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numElems, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numElems, numRefCoords);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, cellTopo);
        ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
        fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
            jacobianAtRefCoords_det,
            basisValuesAtRefCoordsOriented);

        DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords);

        //compute error of projection in L2 norm
        ValueType norm2(0);
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i, double &norm2Update) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k)
              projectedFunAtRefCoords(i,j) += basisCoeffsHVol(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            norm2Update += (funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j))*
                (funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j))*
                weights(j)*jacobianAtRefCoords_det(i,j);
          }
        },norm2);
        ExecSpaceType().fence();

        hvolNorm[iter] =  std::sqrt(norm2);
        auto expected_error = useL2Projection ? hvol_errors_L2[iter] : hvol_errors[iter];
        if(std::abs(hvolNorm[iter]-expected_error)/expected_error > relTol){
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "For N = " << NX << ", computed error (" << hvolNorm[iter] << ") is different than expected one (" << expected_error << ")";
          *outStream << std::endl;
        }
        delete basisPtr;
      }
      if(useL2Projection)
        *outStream << "HVOL Error (L2 Projection): " << hvolNorm[iter] <<std::endl;
      else
        *outStream << "HVOL Error (HVol Projection): " << hvolNorm[iter] <<std::endl;
      }
    } catch (std::exception &err) {
      std::cout << " Exception\n";
      *outStream << err.what() << "\n\n";
      errorFlag = -1000;
    }
  }
/*
  *outStream << "\nHGRAD ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    *outStream << " " << hgradNorm[iter];
  *outStream << "\nHCURL ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    *outStream << " " << hcurlNorm[iter];
  *outStream << "\nHDIV ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    *outStream << " " << hdivNorm[iter];
  *outStream << "\nHVOL ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    *outStream << " " << hvolNorm[iter];
  *outStream << std::endl;
 */

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED = " << errorFlag << "\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  return errorFlag;
}
}
}

