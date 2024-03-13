// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief Test the DeRahm commutativity property of projection-based interpolation for Tri elements

    The test considers two triangles in the physical space sharing a common face.
    In order to test significant configurations, we consider 3 mappings of the reference triangle
    to the first (physical) triangular, so that the common face is mapped from all the 3 sides
    of the reference triangular.
    Then, for each of the mappings, the global ids of the vertices of the common side are permuted.

    The test checks that, for elements of different degree,
    1. projecting a target function f into the HGRAD space, and taking the gradient of the projection
       is equivalent to taking the gradient of f and then project it to the HCURL space.

       projecting a target function f into the HGRAD space, and taking the curl of the projection
       is equivalent to taking the curl of f and then project it to the HDIV space.

    2. projecting a target vector function f into the HCURL space, and taking the curl of the projection
       is equivalent to taking the curl of f and then project it to the HVOL space.

    3. projecting a target vector function f into the HDIV space, and taking the divergence of the projection
       is equivalent to taking the divergence of f and then project it to the HVOL space.

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


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
int DeRhamCommutativityTri(const bool verbose) {

  using ExecSpaceType = typename DeviceType::execution_space;
  using MemSpaceType = typename DeviceType::memory_space;

  using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;
  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;
  using DynRankViewIntHost = Kokkos::DynRankView<ordinal_type,HostSpaceType>;
  
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
  *outStream << "\n";

  int errorFlag = 0;
  const ValueType tol = tolerence();

  //target functions and their derivatives

  struct Fun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      return (x*x*x*x+1)*(y-2);
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 5;}
  };

  struct GradFun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp) {
      switch (comp) {
      case 0:
        return 4*x*x*x*(y - 2);
      case 1:
        return x*x*x*x + 1;
      default:
        return 0;
      }
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 4;}
  };

  struct CurlFun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp) {
      GradFun df;
      return std::pow(-1,comp)*df(x,y,1-comp);
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return GradFun::degree();}
  };


  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a = 2*x*y+x*x*x;
      ValueType f0 = 5+y+x*x;
      ValueType f1 = -7-2+x+y*y;

      //fun = f + a x
      switch (comp) {
      case 0:
        return f0 + a*x;
      case 1:
        return f1 + a*y;
      default:
        return 0;
      }
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 4;}
  };

  struct DivFunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      ValueType a = 2*x*y+x*x*x;
      ValueType ga[2] = {2*y+3*x*x,2*x};
      ValueType divf = 2*x+2*y;

      //fun = div f + x \cdot \nabla a + 2 a
      return divf + (x*ga[0]+y*ga[1]) + 2*a;
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 3;}
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a = (y-7)*x*x*x*x;
      ValueType f0 = 2+x+x*y;
      ValueType f1 = 3;

      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 -a*y;//2*x+y-z + (x+2*(y+z);
      case 1:
        return f1 + a*x;//y+2*(z+x);
      default:
        return 0;
      }
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 5;}
  };

  struct CurlFunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      ValueType a = (y-7)*x*x*x*x;
      ValueType ga[2] = {4*x*x*x*(y-7),x*x*x*x};
      ValueType gf0[2] = {1.0+y,x};
      ValueType gf1[2] = {0, 0};

      //fun = curl f +  2 a - (x \cdot \nabla) a
      return gf1[0] - gf0[1] + 2*a + (x*ga[0]+y*ga[1]);
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 5;}
  };

  using ct = CellTools<DeviceType>;
  using ots = OrientationTools<DeviceType>;
  using pts = ProjectionTools<DeviceType>;
  using ProjStruct = ProjectionStruct<DeviceType,ValueType>;
  using rst = RealSpaceTools<DeviceType>;
  using fst = FunctionSpaceTools<DeviceType>;

  constexpr ordinal_type dim = 2;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 3;
  constexpr ordinal_type numTotalVertexes = 4;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{0,0},{0,1},{1,0},{1,1}};
  ordinal_type tris_orig[numCells][numElemVertexes] = {{0,1,2},{1,2,3}};  
  ordinal_type tris_rotated[numCells][numElemVertexes];

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|           Test 1 (DeRham Commutativity - HGRAD-HCURL, HGRAD-HDIV)           |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    //reordering of nodes to explore different orientations

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3};

    do {

      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tris[numCells][numElemVertexes];
      std::copy(&tris_orig[0][0], &tris_orig[0][0]+numCells*numElemVertexes, &tris_rotated[0][0]);
      for (ordinal_type shift=0; shift<3; ++shift) {
        std::rotate_copy(&tris_orig[0][0], &tris_orig[0][0]+shift, &tris_orig[0][0]+3, &tris_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tris[i][j] = reorder[tris_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tri 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[0][j] << " ";
        *outStream << "] and Tri 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        const ordinal_type numNodesPerElem = tri.getNodeCount();

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hPhysVertexes(i,j,k) = vertices[tris[i][j]][k];
        Kokkos::deep_copy(physVertexes,hPhysVertexes);

        //compute reference points
        Basis_HGRAD_TRI_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewIntHost elemNodesHost(&tris[0][0], 2, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tri);

        Basis_HGRAD_TRI_Cn_FEM<DeviceType,ValueType,ValueType> basis(order);
        Basis_HCURL_TRI_In_FEM<DeviceType,ValueType,ValueType> basisHCurl(order);
        Basis_HDIV_TRI_In_FEM<DeviceType,ValueType,ValueType> basisHDiv(order);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHCurlCardinality = basisHCurl.getCardinality();
        ordinal_type basisHDivCardinality = basisHDiv.getCardinality();

        // compute projection-based interpolation of fun into HGRAD
        DynRankView ConstructWithLabel(basisCoeffsHGrad, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(Fun::degree()),targetDerivCubDegree(GradFun::degree());

          ProjStruct projStruct;
          projStruct.createHGradProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationGradPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numGradPoints = evaluationGradPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
          DynRankView ConstructWithLabel(targetGradAtEvalPoints, numCells, numGradPoints, dim);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalGradPoints, numCells, numGradPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalGradPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              Fun fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                targetAtEvalPoints(i,j) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1));
              }

              GradFun gradFun;
              auto basisValuesAtEvalCurlPoints = Kokkos::subview(linearBasisValuesAtEvalGradPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numGradPoints; ++j) {
                auto evalCurlPoint = Kokkos::subview(evaluationGradPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalCurlPoints, evalCurlPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalGradPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalCurlPoints(k);
                for(int d=0;d<dim;d++)
                  targetGradAtEvalPoints(i,j,d) = gradFun(physEvalGradPoints(i,j,0), physEvalGradPoints(i,j,1), d);

              }
            });
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints);
          fst::mapHGradDataFromPhysToRef(refTargetAtEvalPoints,targetAtEvalPoints);

          DynRankView ConstructWithLabel(refTargetGradAtEvalPoints, numCells, numGradPoints, dim);
          if(numGradPoints>0) {
            DynRankView ConstructWithLabel(jacobian, numCells, numGradPoints, dim, dim);
            ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, tri);
            fst::mapHCurlDataFromPhysToRef(refTargetGradAtEvalPoints,jacobian,targetGradAtEvalPoints);
          }

          pts::getHGradBasisCoeffs(basisCoeffsHGrad,
              refTargetAtEvalPoints,
              refTargetGradAtEvalPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of gradient of fun into HCURL
        DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisHCurl.getCardinality());
        {
          ordinal_type targetCubDegree(GradFun::degree()),targetDerivCubDegree(0);

          ProjStruct projStruct;
          projStruct.createHCurlProjectionStruct(&basisHCurl, targetCubDegree, targetDerivCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationDivPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numDivPoints = evaluationDivPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              GradFun fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                for(ordinal_type d=0; d<dim; ++d)
                  targetAtEvalPoints(i,j,d) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1),d);
              }
            });
          }

          //transform the target function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            fst::mapHCurlDataFromPhysToRef(refTargetAtEvalPoints,jacobian,targetAtEvalPoints);
          }

          DynRankView ConstructWithLabel(refTargetDivAtEvalPoints, numCells, numDivPoints); //zero, curl of grad
          pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
              refTargetAtEvalPoints,
              refTargetDivAtEvalPoints,
              elemOrts,
              &basisHCurl,
              &projStruct);
        }


        // compute gradient of the fun HGRAD projection at reference points
        DynRankView ConstructWithLabel(gradProjFunAtRefCoordsOriented, numCells, numRefCoords, dim);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
          basis.getValues(basisValuesAtRefCoords, refPoints, OPERATOR_GRAD);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basis);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianInv (jacobian_inv, jacobian);
          fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_inv,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  gradProjFunAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          });
        }

        // compute HCURL projection of fun gradient at reference points
        DynRankView ConstructWithLabel(projGradFunAtRefCoordsOriented, numCells, numRefCoords, dim);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisHCurlCardinality, numRefCoords, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisHCurlCardinality, numRefCoords, dim);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisHCurlCardinality, numRefCoords, dim);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisHCurlCardinality, numRefCoords, dim);
          basisHCurl.getValues(basisValuesAtRefCoords, refPoints);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basisHCurl);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numRefCoords, dim, dim);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianInv (jacobian_inv, jacobian);
          fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_inv,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHCurlCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  projGradFunAtRefCoordsOriented(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          });
        }

        // compare the gradient of the target HGRAD projection and the HCURL projection of the gradient of the target at reference points
        auto hostGradProjFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gradProjFunAtRefCoordsOriented);
        auto hostProjGradFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), projGradFunAtRefCoordsOriented);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              error = std::max(std::abs( hostGradProjFun(i,j,d) - hostProjGradFun(i,j,d)), error);
            }
          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Gradient of projection is different than projection of gradient at reference points " << i << "\n";
            *outStream << "Gradient of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostGradProjFun(i,j,0) << "," << hostGradProjFun(i,j,1) << "," << hostGradProjFun(i,j,2)   << ")";
            *outStream << "\nProjection of gradient at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostProjGradFun(i,j,0) << "," << hostProjGradFun(i,j,1) << "," << hostProjGradFun(i,j,2)  << ")";
            *outStream << std::endl;
          }
        }

        // compute projection-based interpolation of curl of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisHDiv.getCardinality());
        {
          ordinal_type targetCubDegree(GradFun::degree()),targetDerivCubDegree(0);

          ProjStruct projStruct;
          projStruct.createHDivProjectionStruct(&basisHDiv, targetCubDegree, targetDerivCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationDivPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numDivPoints = evaluationDivPoints.extent(0);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              CurlFun fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                for(ordinal_type d=0; d<dim; ++d)
                  targetAtEvalPoints(i,j,d) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1),d);
              }
            });
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);

            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            ct::setJacobianInv(jacobian_inv, jacobian);
            ct::setJacobianDet(jacobian_det, jacobian);
            fst::mapHDivDataFromPhysToRef(refTargetAtEvalPoints,jacobian_inv,jacobian_det,targetAtEvalPoints);
          }

          DynRankView ConstructWithLabel(refTargetDivAtEvalPoints, numCells, numDivPoints); //zero, div of curl
          pts::getHDivBasisCoeffs(basisCoeffsHDiv,
              refTargetAtEvalPoints,
              refTargetDivAtEvalPoints,
              elemOrts,
              &basisHDiv,
              &projStruct);
        }

        // compute curl of the fun HCURL projection at reference points
        DynRankView ConstructWithLabel(curlProjFunAtRefCoordsOriented, numCells, numRefCoords, dim);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
          basis.getValues(basisValuesAtRefCoords, refPoints, OPERATOR_CURL);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basis);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HGRADtransformCURL(transformedBasisValuesAtRefCoordsOriented,
              jacobian, jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  curlProjFunAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          });
        }

        // compute HDIV projection of fun curl at reference points
        DynRankView ConstructWithLabel(projCurlFunAtRefCoordsOriented, numCells, numRefCoords, dim);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisHDivCardinality, numRefCoords, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisHDivCardinality, numRefCoords, dim);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisHDivCardinality, numRefCoords, dim);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisHDivCardinality, numRefCoords, dim);
          basisHDiv.getValues(basisValuesAtRefCoords, refPoints);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basisHDiv);

          // transform basis to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian, jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHDivCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  projCurlFunAtRefCoordsOriented(i,j,d) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          });
        }

        // compare the curl of the target HGRAD projection and the HDIV projection of the curl of the target at reference points
        auto hostCurlProjFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), curlProjFunAtRefCoordsOriented);
        auto hostProjCurlFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), projCurlFunAtRefCoordsOriented);

        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              error = std::max(std::abs( hostCurlProjFun(i,j,d) - hostProjCurlFun(i,j,d)), error);
            }
          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Curl of projection is different than projection of curl at reference points " << i << "\n";
            *outStream << "curl of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostCurlProjFun(i,j,0) << "," << hostCurlProjFun(i,j,1)   << ")";
            *outStream << "\nProjection of curl at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostProjCurlFun(i,j,0) << "," << hostProjCurlFun(i,j,1)  << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception &err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }



  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|              Test 2 (DeRham Commutativity - HCURL-HVOL)                     |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tris[numCells][numElemVertexes];
      std::copy(&tris_orig[0][0], &tris_orig[0][0]+numCells*numElemVertexes, &tris_rotated[0][0]);
      for (ordinal_type shift=0; shift<3; ++shift) {
        std::rotate_copy(&tris_orig[0][0], &tris_orig[0][0]+shift, &tris_orig[0][0]+3, &tris_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tris[i][j] = reorder[tris_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tri 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[0][j] << " ";
        *outStream << "] and Tri 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        const ordinal_type numNodesPerElem = tri.getNodeCount();

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hostPhysVertexes(i,j,k) = vertices[tris[i][j]][k];
        deep_copy(physVertexes, hostPhysVertexes);


        //compute reference points
        Basis_HGRAD_TRI_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewIntHost elemNodesHost(&tris[0][0], 2, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tri);


        Basis_HCURL_TRI_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        Basis_HVOL_TRI_Cn_FEM<DeviceType,ValueType,ValueType> basisHVol(order-1);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHVolCardinality = basisHVol.getCardinality();

        // compute projection-based interpolation of fun into HVOL
        DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(FunCurl::degree()),targetDerivCubDegree(CurlFunCurl::degree());


          ProjStruct projStruct;
          projStruct.createHCurlProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationCurlPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numCurlPoints = evaluationCurlPoints.extent(0);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numCurlPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalCurlPoints, numCells, numCurlPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalCurlPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              FunCurl fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                for(ordinal_type d=0; d<dim; ++d)
                  targetAtEvalPoints(i,j,d) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1),d);
              }

              CurlFunCurl curlFun;
              auto basisValuesAtEvalCurlPoints = Kokkos::subview(linearBasisValuesAtEvalCurlPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numCurlPoints; ++j) {
                auto evalCurlPoint = Kokkos::subview(evaluationCurlPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalCurlPoints, evalCurlPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalCurlPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalCurlPoints(k);
                targetCurlAtEvalPoints(i,j) = curlFun(physEvalCurlPoints(i,j,0), physEvalCurlPoints(i,j,1));
              }
            });
          }



          //transform the function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            fst::mapHCurlDataFromPhysToRef(refTargetAtEvalPoints,jacobian,targetAtEvalPoints);
          }

          DynRankView ConstructWithLabel(refTargetCurlAtEvalPoints, numCells, numCurlPoints);
          if(numCurlPoints>0) {

            DynRankView ConstructWithLabel(jacobianCurl, numCells, numCurlPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobianCurl_det, numCells, numCurlPoints);
            ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, tri);
            ct::setJacobianDet (jacobianCurl_det, jacobianCurl);

            fst::mapHVolDataFromPhysToRef(refTargetCurlAtEvalPoints,jacobianCurl_det,targetCurlAtEvalPoints);
          }

          pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
              refTargetAtEvalPoints,
              refTargetCurlAtEvalPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of curl of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisHVol.getCardinality());
        {
          ordinal_type targetCubDegree(CurlFunCurl::degree());

          ProjStruct projStruct;
          projStruct.createHVolProjectionStruct(&basisHVol, targetCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);



          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              CurlFunCurl fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                targetAtEvalPoints(i,j) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1));
              }
            });
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            ct::setJacobianInv(jacobian_inv,jacobian);
            ct::setJacobianDet (jacobian_det, jacobian);
            fst::mapHVolDataFromPhysToRef(refTargetAtEvalPoints,jacobian_det,targetAtEvalPoints);
          }

          pts::getHVolBasisCoeffs(basisCoeffsHDiv,
              refTargetAtEvalPoints,
              elemOrts,
              &basisHVol,
              &projStruct);
        }


        DynRankView ConstructWithLabel(curlProjFunAtRefCoordsOriented, numCells, numRefCoords);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
          basis.getValues(basisValuesAtRefCoords, refPoints, OPERATOR_CURL);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basis);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HCURLtransformCURL(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                curlProjFunAtRefCoordsOriented(i,j) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          });
        }

        // compute HVOL projection of fun curl at reference points
        DynRankView ConstructWithLabel(projCurlFunAtRefCoordsOriented, numCells, numRefCoords);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisHVolCardinality, numRefCoords);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisHVolCardinality, numRefCoords);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisHVolCardinality, numRefCoords);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisHVolCardinality, numRefCoords);
          basisHVol.getValues(basisValuesAtRefCoords, refPoints);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basisHVol);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHVolCardinality; ++k)
                projCurlFunAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          });
        }

        // compare the curl of the target HCURL projection and the HVOL projection of the curl of the target at reference points
        auto hostCurlProjFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), curlProjFunAtRefCoordsOriented);
        auto hostProjCurlFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), projCurlFunAtRefCoordsOriented);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            error = std::max(std::abs( hostCurlProjFun(i,j) - hostProjCurlFun(i,j)), error);

          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Curl of projection is different than projection of curl at reference points " << i << "\n";
            *outStream << "Curl of rojection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostCurlProjFun(i,j)  << ")";
            *outStream << "\nProjection of curl at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostProjCurlFun(i,j)  << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception &err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }


  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                Test 3 (DeRham Commutativity - HDIV-HVOL)                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3};

    do {

      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tris[numCells][numElemVertexes];
      std::copy(&tris_orig[0][0], &tris_orig[0][0]+numCells*numElemVertexes, &tris_rotated[0][0]);
      for (ordinal_type shift=0; shift<3; ++shift) {
        std::rotate_copy(&tris_orig[0][0], &tris_orig[0][0]+shift, &tris_orig[0][0]+3, &tris_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tris[i][j] = reorder[tris_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tri 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[0][j] << " ";
        *outStream << "] and Tri 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tris[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        const ordinal_type numNodesPerElem = tri.getNodeCount();

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hostPhysVertexes(i,j,k) = vertices[tris[i][j]][k];
        deep_copy(physVertexes, hostPhysVertexes);

        // compute orientations for cells (one time computation)
        DynRankViewIntHost elemNodesHost(&tris[0][0], numCells, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tri);

        //compute reference points
        Basis_HGRAD_TRI_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_TRI_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        Basis_HVOL_TRI_Cn_FEM<DeviceType,ValueType,ValueType> basisHVol(order-1);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHVolCardinality = basisHVol.getCardinality();

        // compute projection-based interpolation of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(FunDiv::degree()),targetDerivCubDegree(DivFunDiv::degree());


          ProjStruct projStruct;
          projStruct.createHDivProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          auto evaluationDivPoints = projStruct.getAllDerivEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0), numDivPoints = evaluationDivPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalDivPoints, numCells, numDivPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalDivPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              FunDiv fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                for(ordinal_type d=0; d<dim; ++d)
                  targetAtEvalPoints(i,j,d) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1),d);
              }

              DivFunDiv divFun;
              auto basisValuesAtEvalDivPoints = Kokkos::subview(linearBasisValuesAtEvalDivPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numDivPoints; ++j) {
                auto evalDivPoint = Kokkos::subview(evaluationDivPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalDivPoints, evalDivPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalDivPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalDivPoints(k);
                targetDivAtEvalPoints(i,j) = divFun(physEvalDivPoints(i,j,0), physEvalDivPoints(i,j,1));
              }
            });
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
            DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            ct::setJacobianDet (jacobian_det, jacobian);
            ct::setJacobianInv (jacobian_inv, jacobian);
            fst::mapHDivDataFromPhysToRef(refTargetAtEvalPoints,jacobian_inv,jacobian_det,targetAtEvalPoints);
          }

          DynRankView ConstructWithLabel(refTargetDivAtEvalPoints, numCells, numDivPoints);
          if(numDivPoints) {
            DynRankView ConstructWithLabel(jacobianDiv, numCells, numDivPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobianDiv_det, numCells, numDivPoints);
            ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, tri);
            ct::setJacobianDet (jacobianDiv_det, jacobianDiv);
            fst::mapHVolDataFromPhysToRef(refTargetDivAtEvalPoints,jacobianDiv_det,targetDivAtEvalPoints);
          }

          pts::getHDivBasisCoeffs(basisCoeffsHDiv,
              refTargetAtEvalPoints,
              refTargetDivAtEvalPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of divergence of fun into HVOL
        DynRankView ConstructWithLabel(basisCoeffsHVol, numCells, basisHVol.getCardinality());
        {
          ordinal_type targetCubDegree(DivFunDiv::degree());

          ProjStruct projStruct;
          projStruct.createHVolProjectionStruct(&basisHVol, targetCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoints, numCells, numNodesPerElem);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              DivFunDiv fun;
              auto basisValuesAtEvalPoints = Kokkos::subview(linearBasisValuesAtEvalPoints,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j) {
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoints, evalPoint);
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    physEvalPoints(i,j,d) += physVertexes(i,k,d)*basisValuesAtEvalPoints(k);
                targetAtEvalPoints(i,j) = fun(physEvalPoints(i,j,0), physEvalPoints(i,j,1));
              }
            });
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(refTargetAtEvalPoints, numCells, numPoints);
          {
            DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
            DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
            ct::setJacobian(jacobian, evaluationPoints, physVertexes, tri);
            ct::setJacobianDet (jacobian_det, jacobian);
            fst::mapHVolDataFromPhysToRef(refTargetAtEvalPoints,jacobian_det,targetAtEvalPoints);
          }

          pts::getHVolBasisCoeffs(basisCoeffsHVol,
              refTargetAtEvalPoints,
              elemOrts,
              &basisHVol,
              &projStruct);
        }

        // compute divergence of the fun HDIV projection at reference points
        DynRankView ConstructWithLabel(divProjFunAtRefCoordsOriented, numCells, numRefCoords);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
          basis.getValues(basisValuesAtRefCoords, refPoints, OPERATOR_DIV);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basis);

          // transform basis values to reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HDIVtransformDIV(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                divProjFunAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          });
        }

        // compute HVOL projection of fun divergence at reference points
        DynRankView ConstructWithLabel(projDivFunAtRefCoordsOriented, numCells, numRefCoords);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisHVolCardinality, numRefCoords);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisHVolCardinality, numRefCoords);
          DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisHVolCardinality, numRefCoords);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisHVolCardinality, numRefCoords);
          basisHVol.getValues(basisValuesAtRefCoords, refPoints);
          rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisValuesAtRefCoordsCells,
              elemOrts,
              &basisHVol);

          // transform basis values
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tri);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHVolCardinality; ++k)
                projDivFunAtRefCoordsOriented(i,j) += basisCoeffsHVol(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          });
        }

        // compare the divergence of the target HDIV projection and the HVOL projection of the divergence of the target at reference points
        auto hostDivProjFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), divProjFunAtRefCoordsOriented);
        auto hostProjDivFun = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), projDivFunAtRefCoordsOriented);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            error = std::max(std::abs( hostDivProjFun(i,j) - hostProjDivFun(i,j)), error);

          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Divergence of projection is different than projection of divergence at reference points " << i << "\n";
            *outStream << "Divergence of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostDivProjFun(i,j)    << ")";
            *outStream << "\nProjection of divergence at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << hostProjDivFun(i,j)   << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception &err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }


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

