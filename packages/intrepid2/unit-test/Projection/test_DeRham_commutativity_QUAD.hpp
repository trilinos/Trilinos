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
    \brief Test the DeRahm commutativity property of projection-based interpolation for Quad elements

    The test considers two quadrilaterals in the physical space sharing a common face.
    In order to test significant configurations, we consider 4 mappings of the reference quadrilateral
    to the first (physical) quadrilateral, so that the common face is mapped from all the 4 sides
    of the reference quadrilateral.
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
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#define Intrepid2_Experimental


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
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceSpaceType>
int DeRhamCommutativityQuad(const bool verbose) {

  typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
  typedef Kokkos::DynRankView<ordinal_type,DeviceSpaceType> DynRankViewInt;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

  *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
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
    degree() {return 4;}
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
    degree() {return 4;}
  };

  typedef std::array<ordinal_type,2> edgeType;
  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef Experimental::ProjectionTools<DeviceSpaceType> pts;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;

  constexpr ordinal_type dim = 2;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 4;
  constexpr ordinal_type numTotalVertexes = 6;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1},{1,-1},{1,1},{-1,1}, {-1,2},{1,2}};
  ordinal_type quads_orig[numCells][numElemVertexes] = {{0,1,2,3},{2,3,4,5}};  
  edgeType common_edge = {{2,3}};
  ordinal_type quads_rotated[numCells][numElemVertexes];

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|           Test 1 (DeRham Commutativity - HGRAD-HCURL, HGRAD-HDIV)           |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    //reordering of nodes to explore different orientations

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5};

    do {

      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type quads[numCells][numElemVertexes];
      std::copy(&quads_orig[0][0], &quads_orig[0][0]+numCells*numElemVertexes, &quads_rotated[0][0]);
      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&quads_orig[0][0], &quads_orig[0][0]+shift, &quads_orig[0][0]+4, &quads_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            quads[i][j] = reorder[quads_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Quad 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[0][j] << " ";
        *outStream << "] and Quad 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, quad.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<quad.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[quads[i][j]][k];

        //compute reference points
        Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], 2, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);

        Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HCURL_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basisHCurl(order);
        Basis_HDIV_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basisHDiv(order);

        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHCurlCardinality = basisHCurl.getCardinality();
        ordinal_type basisHDivCardinality = basisHDiv.getCardinality();


        // compute projection-based interpolation of fun into HGRAD
        DynRankView ConstructWithLabel(basisCoeffsHGrad, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(Fun::degree()),targetDerivCubDegree(GradFun::degree());

          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHGradProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numGradPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationGradPoints, numCells, numGradPoints, dim);

          pts::getHGradEvaluationPoints(evaluationPoints,
              evaluationGradPoints,
              elemOrts,
              &basis,
              &projStruct);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
          DynRankView ConstructWithLabel(targetGradAtEvalPoints, numCells, numGradPoints, dim);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalGradPoints, numCells, numGradPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalGradPoints, quad.getNodeCount(), numGradPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numGradPoints>0)
                quadLinearBasis.getValues(quadLinearBasisValuesAtEvalGradPoints, Kokkos::subview(evaluationGradPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numGradPoints; ++j)
                    physEvalGradPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalGradPoints(k,j);
                }
              }
            }
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numGradPoints, dim, dim);
          if(numGradPoints>0)
            ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, quad);

          Fun fun;
          GradFun gradFun;
          Kokkos::deep_copy(targetGradAtEvalPoints,0.);
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++) {
              targetAtEvalPoints(ic,i) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
            }
            for(int i=0;i<numGradPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetGradAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*gradFun(physEvalGradPoints(ic,i,0), physEvalGradPoints(ic,i,1), d);//funHGradCoeffs(k)
            }

          }

          pts::getHGradBasisCoeffs(basisCoeffsHGrad,
              targetAtEvalPoints,
              targetGradAtEvalPoints,
              evaluationPoints,
              evaluationGradPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of gradient of fun into HCURL
        DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisHCurl.getCardinality());
        {
          ordinal_type targetCubDegree(GradFun::degree()),targetDerivCubDegree(0);

          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHCurlProjectionStruct(&basisHCurl, targetCubDegree, targetDerivCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationDivPoints, numCells, numDivPoints, dim);

          pts::getHCurlEvaluationPoints(evaluationPoints,
              evaluationDivPoints,
              elemOrts,
              &basisHCurl,
              &projStruct);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the target function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);

          GradFun fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
          }

          pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
              targetAtEvalPoints,
              targetDivAtEvalPoints,
              evaluationPoints,
              evaluationDivPoints,
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianInv (jacobian_inv, jacobian);
          fst::HGRADtransformGRAD(transformedBasisValuesAtRefCoordsOriented,
              jacobian_inv,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  gradProjFunAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          }
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianInv (jacobian_inv, jacobian);
          fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_inv,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHCurlCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  projGradFunAtRefCoordsOriented(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          }
        }

        // compare the gradient of the target HGRAD projection and the HCURL projection of the gradient of the target at reference points
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              error = std::max(std::abs( gradProjFunAtRefCoordsOriented(i,j,d) - projGradFunAtRefCoordsOriented(i,j,d)), error);
            }
          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Gradient of projection is different than projection of gradient at reference points " << i << "\n";
            *outStream << "Gradient of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << gradProjFunAtRefCoordsOriented(i,j,0) << "," << gradProjFunAtRefCoordsOriented(i,j,1) << "," << gradProjFunAtRefCoordsOriented(i,j,2)   << ")";
            *outStream << "\nProjection of gradient at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << projGradFunAtRefCoordsOriented(i,j,0) << "," << projGradFunAtRefCoordsOriented(i,j,1) << "," << projGradFunAtRefCoordsOriented(i,j,2)  << ")";
            *outStream << std::endl;
          }
        }

        // compute projection-based interpolation of curl of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisHDiv.getCardinality());
        {
          ordinal_type targetCubDegree(GradFun::degree()),targetDerivCubDegree(0);

          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHDivProjectionStruct(&basisHDiv, targetCubDegree, targetDerivCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationDivPoints, numCells, numDivPoints, dim);

          pts::getHDivEvaluationPoints(evaluationPoints,
              evaluationDivPoints,
              elemOrts,
              &basisHDiv,
              &projStruct);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numDivPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);

          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
          ct::setJacobianInv(jacobian_inv, jacobian);
          ct::setJacobianDet(jacobian_det, jacobian);

          CurlFun fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
          }

          pts::getHDivBasisCoeffs(basisCoeffsHDiv,
              targetAtEvalPoints,
              targetCurlAtEvalPoints,
              evaluationPoints,
              evaluationDivPoints,
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HGRADtransformCURL(transformedBasisValuesAtRefCoordsOriented,
              jacobian, jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  curlProjFunAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          }
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian, jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHDivCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  projCurlFunAtRefCoordsOriented(i,j,d) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          }
        }

        // compare the curl of the target HGRAD projection and the HDIV projection of the curl of the target at reference points
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              error = std::max(std::abs( curlProjFunAtRefCoordsOriented(i,j,d) - projCurlFunAtRefCoordsOriented(i,j,d)), error);
            }
          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Curl of projection is different than projection of curl at reference points " << i << "\n";
            *outStream << "curl of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << curlProjFunAtRefCoordsOriented(i,j,0) << "," << curlProjFunAtRefCoordsOriented(i,j,1)   << ")";
            *outStream << "\nProjection of curl at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << projCurlFunAtRefCoordsOriented(i,j,0) << "," << projCurlFunAtRefCoordsOriented(i,j,1)  << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception err) {
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
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5};

    do {

      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type quads[numCells][numElemVertexes];
      std::copy(&quads_orig[0][0], &quads_orig[0][0]+numCells*numElemVertexes, &quads_rotated[0][0]);
      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&quads_orig[0][0], &quads_orig[0][0]+shift, &quads_orig[0][0]+4, &quads_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            quads[i][j] = reorder[quads_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Quad 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[0][j] << " ";
        *outStream << "] and Quad 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, quad.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<quad.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[quads[i][j]][k];

        //compute reference points
        Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], 2, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);


        Basis_HCURL_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HVOL_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basisHVol(order-1);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHVolCardinality = basisHVol.getCardinality();

        DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(FunCurl::degree()),targetDerivCubDegree(CurlFunCurl::degree());


          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHCurlProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numCurlPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationCurlPoints, numCells, numCurlPoints, dim);

          pts::getHCurlEvaluationPoints(evaluationPoints,
              evaluationCurlPoints,
              elemOrts,
              &basis,
              &projStruct);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numCurlPoints);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalCurlPoints, numCells, numCurlPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalCurlPoints, quad.getNodeCount(), numCurlPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numCurlPoints>0)
                quadLinearBasis.getValues(quadLinearBasisValuesAtEvalCurlPoints, Kokkos::subview(evaluationCurlPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numCurlPoints; ++j)
                    physEvalCurlPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalCurlPoints(k,j);
                }
              }
            }
          }

          //transform the function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);

          DynRankView ConstructWithLabel(jacobianCurl, numCells, numCurlPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianCurl_inv, numCells, numCurlPoints, dim, dim);
          //  DynRankView ConstructWithLabel(jacobian_invT, numCells, numCurlPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianCurl_det, numCells, numCurlPoints);
          if(numCurlPoints>0) {
            ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, quad);
            ct::setJacobianInv (jacobianCurl_inv, jacobianCurl);
            ct::setJacobianDet (jacobianCurl_det, jacobianCurl);
          }

          FunCurl fun;
          CurlFunCurl curlFun;
          Kokkos::deep_copy(targetCurlAtEvalPoints,0.);
          Kokkos::deep_copy(targetAtEvalPoints,0.);
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++) {
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
            }
            for(int i=0;i<numCurlPoints;i++) {
              targetCurlAtEvalPoints(ic,i) += jacobianCurl_det(ic,i)*curlFun(physEvalCurlPoints(ic,i,0), physEvalCurlPoints(ic,i,1));
            }

          }

          pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
              targetAtEvalPoints,
              targetCurlAtEvalPoints,
              evaluationPoints,
              evaluationCurlPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of curl of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisHVol.getCardinality());
        {
          ordinal_type targetCubDegree(CurlFunCurl::degree());

          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHVolProjectionStruct(&basisHVol, targetCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);

          pts::getHVolEvaluationPoints(evaluationPoints,
              elemOrts,
              &basisHVol,
              &projStruct);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
          ct::setJacobianInv(jacobian_inv,jacobian);
          ct::setJacobianDet (jacobian_det, jacobian);

          CurlFunCurl fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
          }

          pts::getHVolBasisCoeffs(basisCoeffsHDiv,
              targetAtEvalPoints,
              evaluationPoints,
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HCURLtransformCURL(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                curlProjFunAtRefCoordsOriented(i,j) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          }
        }

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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHVolCardinality; ++k)
                projCurlFunAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          }
        }

        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            error = std::max(std::abs( curlProjFunAtRefCoordsOriented(i,j) - projCurlFunAtRefCoordsOriented(i,j)), error);

          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Curl of projection is different than projection of curl at reference points " << i << "\n";
            *outStream << "Curl of rojection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << curlProjFunAtRefCoordsOriented(i,j)  << ")";
            *outStream << "\nProjection of curl at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << projCurlFunAtRefCoordsOriented(i,j)  << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception err) {
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
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5};

    do {

      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type quads[numCells][numElemVertexes];
      std::copy(&quads_orig[0][0], &quads_orig[0][0]+numCells*numElemVertexes, &quads_rotated[0][0]);
      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&quads_orig[0][0], &quads_orig[0][0]+shift, &quads_orig[0][0]+4, &quads_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            quads[i][j] = reorder[quads_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Quad 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[0][j] << " ";
        *outStream << "] and Quad 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, quad.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<quad.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[quads[i][j]][k];


        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);

        //compute reference points
        Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HVOL_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basisHVol(order-1);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHVolCardinality = basisHVol.getCardinality();

        // compute projection-based interpolation of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(FunDiv::degree()),targetDerivCubDegree(DivFunDiv::degree());


          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHDivProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();

          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationDivPoints, numCells, numDivPoints, dim);


          pts::getHDivEvaluationPoints(evaluationPoints,
              evaluationDivPoints,
              elemOrts,
              &basis,
              &projStruct);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);
          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalDivPoints, numCells, numDivPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalDivPoints, quad.getNodeCount(), numDivPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numDivPoints)
                quadLinearBasis.getValues(quadLinearBasisValuesAtEvalDivPoints, Kokkos::subview(evaluationDivPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numDivPoints; ++j)
                    physEvalDivPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalDivPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          ct::setJacobianInv (jacobian_inv, jacobian);

          DynRankView ConstructWithLabel(jacobianDiv, numCells, numDivPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianDiv_det, numCells, numDivPoints);
          if(numDivPoints) {
            ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, quad);
            ct::setJacobianDet (jacobianDiv_det, jacobianDiv);
          }

          FunDiv fun;
          DivFunDiv divFun;
          Kokkos::deep_copy(targetDivAtEvalPoints,0.);
          Kokkos::deep_copy(targetAtEvalPoints,0.);
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++) {
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
            }
            for(int i=0;i<numDivPoints;i++)
              targetDivAtEvalPoints(ic,i) += jacobianDiv_det(ic,i)*divFun(physEvalDivPoints(ic,i,0), physEvalDivPoints(ic,i,1));//funHGradCoeffs(k)
          }

          pts::getHDivBasisCoeffs(basisCoeffsHDiv,
              targetAtEvalPoints,
              targetDivAtEvalPoints,
              evaluationPoints,
              evaluationDivPoints,
              elemOrts,
              &basis,
              &projStruct);
        }

        // compute projection-based interpolation of divergence of fun into HVOL
        DynRankView ConstructWithLabel(basisCoeffsHVol, numCells, basisHVol.getCardinality());
        {
          ordinal_type targetCubDegree(DivFunDiv::degree());

          Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
          projStruct.createHVolProjectionStruct(&basisHVol, targetCubDegree);

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationDivPoints, numCells, numDivPoints, dim);

          pts::getHVolEvaluationPoints(evaluationPoints,
              elemOrts,
              &basisHVol,
              &projStruct);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);

          DivFunDiv fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
          }

          pts::getHVolBasisCoeffs(basisCoeffsHVol,
              targetAtEvalPoints,
              evaluationPoints,
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HDIVtransformDIV(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                divProjFunAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          }
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
          ct::setJacobian(jacobian, refPoints, physVertexes, quad);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
              jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisHVolCardinality; ++k)
                projDivFunAtRefCoordsOriented(i,j) += basisCoeffsHVol(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            }
          }
        }

        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            error = std::max(std::abs( divProjFunAtRefCoordsOriented(i,j) - projDivFunAtRefCoordsOriented(i,j)), error);

          if(error>10000*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE! " << error << "\n";
            *outStream << "Divergence of projection is different than projection of divergence at reference points " << i << "\n";
            *outStream << "Divergence of projection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << divProjFunAtRefCoordsOriented(i,j)    << ")";
            *outStream << "\nProjection of divergence at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << projDivFunAtRefCoordsOriented(i,j)   << ")";
            *outStream << std::endl;
          }
        }
      }
    } while(0);//std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

  } catch (std::exception err) {
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

