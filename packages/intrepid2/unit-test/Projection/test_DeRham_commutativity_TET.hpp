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
    \brief Test the DeRahm commutativity property of projection-based interpolation for tet elements

    The test considers two tetrahedra in the physical space sharing a common face.
    In order to test significant configurations, we consider 4 mappings of the reference tetrahedron
    to the first (physical) tetrahedron, so that the common face is mapped from all the 4 faces
    of the reference tetrahedron.
    Then, for each of the mappings, the global ids of the vertices of the common face are permuted.

    The test checks that, for elements of different degree,
    1. projecting a target function f into the HGRAD space, and taking the gradient of the projection
       is equivalent to taking the gradient of f and then project it to the HCURL space.

    2. projecting a target vector function f into the HCURL space, and taking the curl of the projection
       is equivalent to taking the curl of f and then project it to the HDIV space.

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
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_TET_Cn_FEM.hpp"
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
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceSpaceType>
int DeRhamCommutativityTet(const bool verbose) {

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
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return (x*x*x*x+1)*(y-2);//*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 5;}
  };

  struct GradFun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      switch (comp) {
      case 0:
        return 4*x*x*x*(y - 2);
      case 1:
        return x*x*x*x + 1;
      case 2:
        return 0;
      default:
        return 0;
      }
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 4;}
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a = 2*x*y+x*x*x;
      ValueType f0 = 5+y+x*x+z*z;
      ValueType f1 = -7-2*z+x+y*y+z*z;
      ValueType f2 = 0.5+z*z+x*x;
      //fun = f + a x
      switch (comp) {
      case 0:
        return f0 + a*x;
      case 1:
        return f1 + a*y;
      case 2:
        return f2 + a*z;
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
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      ValueType a = 2*x*y+x*x*x;
      ValueType ga[3] = {2*y+3*x*x,2*x,0};
      ValueType divf = 2*x+2*y+2*z;

      //fun = div f + x \cdot \nabla a + 3 a
      return divf + (x*ga[0]+y*ga[1]+z*ga[2]) + 3*a;
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 3;}
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a0 = y-7+z*z*z*z;
      ValueType a1 = 2*z-1+z*x;
      ValueType a2 = z-2+x*x;
      ValueType f0 = 2+x+z+x*y;
      ValueType f1 = 3-3*z;
      ValueType f2 = -5+x;

      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 + (a1*z-a2*y);//2*x+y-z + (x+2*(y+z);
      case 1:
        return f1 + (a2*x-a0*z);//y+2*(z+x);
      case 2:
        return f2 + (a0*y-a1*x);//z+2*(x+y);
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
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a0 = y-7+z*z*z*z;
      ValueType a1 = 2*z-1+z*x;
      ValueType a2 = z-2+x*x;
      ValueType ga0[3] = {0,1,4*z*z*z};
      ValueType ga1[3] = {z, 0, 2+x};
      ValueType ga2[3] = {2*x, 0, 1.0};
      ValueType gf0[3] = {1.0+y,x,1};
      ValueType gf1[3] = {0, 0, -3};
      ValueType gf2[3] = {1, 0, 0};
      ValueType diva = ga0[0]+ga1[1]+ga2[2];

      //fun = curl f +  3 a - (x \cdot \nabla) a -  (div a) x - a
      switch (comp) {
      case 0:
        return gf2[1] - gf1[2] + 2*a0 + (x*ga0[0]+y*ga0[1]+z*ga0[2]) - x*diva;//2*x+y-z + (x+2*(y+z);
      case 1:
        return gf0[2] - gf2[0] + 2*a1 + (x*ga1[0]+y*ga1[1]+z*ga1[2]) - y*diva;//y+2*(z+x);
      case 2:
        return gf1[0] - gf0[1] + 2*a2 + (x*ga2[0]+y*ga2[1]+z*ga2[2]) - z*diva;//z+2*(x+y);
      default:
        return 0;
      }
    }
    static ordinal_type
    KOKKOS_INLINE_FUNCTION
    degree() {return 4;}
  };

  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef Experimental::ProjectionTools<DeviceSpaceType> pts;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 4;
  constexpr ordinal_type numTotalVertexes = 5;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
  ordinal_type tets_orig[numCells][numElemVertexes] = {{0,1,2,3},{1,2,3,4}};  //common face is {1,2,3}
  ordinal_type tets_rotated[numCells][numElemVertexes];

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|              Test 1 (DeRham Commutativity - HGRAD-HCURL)                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    //reordering of nodes to explore different orientations

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tets[numCells][numElemVertexes];
      std::copy(&tets_orig[0][0], &tets_orig[0][0]+numCells*numElemVertexes, &tets_rotated[0][0]);

      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&tets_orig[0][0], &tets_orig[0][0]+shift, &tets_orig[0][0]+4, &tets_rotated[0][0]);

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tets[i][j] = reorder[tets_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tet 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[0][j] << " ";
        *outStream << "] and Tet 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, tet.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<tet.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[tets[i][j]][k];

        //compute reference points
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&tets[0][0], 2, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);

        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HCURL_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basisHCurl(order);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHCurlCardinality = basisHCurl.getCardinality();


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
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalGradPoints, tet.getNodeCount(), numGradPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numGradPoints>0)
                tetLinearBasis.getValues(tetLinearBasisValuesAtEvalGradPoints, Kokkos::subview(evaluationGradPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numGradPoints; ++j)
                    physEvalGradPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalGradPoints(k,j);
                }
              }
            }
          }

          //transform the target function and its derivative to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numGradPoints, dim, dim);
          if(numGradPoints>0)
            ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, tet);

          Fun fun;
          GradFun gradFun;
          Kokkos::deep_copy(targetGradAtEvalPoints,0.);
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++) {
              targetAtEvalPoints(ic,i) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2));
            }
            for(int i=0;i<numGradPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetGradAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*gradFun(physEvalGradPoints(ic,i,0), physEvalGradPoints(ic,i,1), physEvalGradPoints(ic,i,2), d);//funHGradCoeffs(k)
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
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints, dim);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the target function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);

          GradFun fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
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
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
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
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
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
          if(error>1000*tol) {
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
  << "|              Test 2 (DeRham Commutativity - HCURL-HDIV)                     |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tets[numCells][numElemVertexes];
      std::copy(&tets_orig[0][0], &tets_orig[0][0]+numCells*numElemVertexes, &tets_rotated[0][0]);

      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&tets_orig[0][0], &tets_orig[0][0]+shift, &tets_orig[0][0]+4, &tets_rotated[0][0]);

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tets[i][j] = reorder[tets_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d){
            vertices[i][d] = vertices_orig[orderback[i]][d];
          }


        *outStream <<  "Considering Tet 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[0][j] << " ";
        *outStream << "] and Tet 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, tet.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<tet.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[tets[i][j]][k];

        //compute reference points
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&tets[0][0], 2, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);


        Basis_HCURL_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HDIV_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basisHDiv(order);
        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type basisHDivCardinality = basisHDiv.getCardinality();

        // compute projection-based interpolation of fun into HCURL
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
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numCurlPoints, dim);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalCurlPoints, numCells, numCurlPoints, dim);
          {
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalCurlPoints, tet.getNodeCount(), numCurlPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numCurlPoints>0)
                tetLinearBasis.getValues(tetLinearBasisValuesAtEvalCurlPoints, Kokkos::subview(evaluationCurlPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numCurlPoints; ++j)
                    physEvalCurlPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalCurlPoints(k,j);
                }
              }
            }
          }

          //need to transform dofCoeff to physical space (they transform as normals)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);

          DynRankView ConstructWithLabel(jacobianCurl, numCells, numCurlPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianCurl_inv, numCells, numCurlPoints, dim, dim);
          //  DynRankView ConstructWithLabel(jacobian_invT, numCells, numCurlPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianCurl_det, numCells, numCurlPoints);
          if(numCurlPoints>0) {
            ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, tet);
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
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
            }
            for(int i=0;i<numCurlPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetCurlAtEvalPoints(ic,i,j) += jacobianCurl_det(ic,i)*jacobianCurl_inv(ic,i,j,d)*curlFun(physEvalCurlPoints(ic,i,0), physEvalCurlPoints(ic,i,1), physEvalCurlPoints(ic,i,2), d);//funHGradCoeffs(k)
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
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisHDiv.getCardinality());
        {
          ordinal_type targetCubDegree(CurlFunCurl::degree()),targetDerivCubDegree(0);

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
          DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);

          DynRankView ConstructWithLabel(physEvalPoints, numCells, numPoints, dim);
          {
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);
          ct::setJacobianInv(jacobian_inv,jacobian);
          ct::setJacobianDet (jacobian_det, jacobian);

          CurlFunCurl fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              for(int j=0;j<dim;j++)
                for(int d=0;d<dim;d++)
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
          }

          pts::getHDivBasisCoeffs(basisCoeffsHDiv,
              targetAtEvalPoints,
              targetDivAtEvalPoints,
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
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
          ct::setJacobianDet (jacobian_det, jacobian);
          fst::HCURLtransformCURL(transformedBasisValuesAtRefCoordsOriented,
              jacobian, jacobian_det,
              basisValuesAtRefCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  curlProjFunAtRefCoordsOriented(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            }
          }
        }

        // compute HDIV projection of fun curl at reference points
        DynRankView ConstructWithLabel(projCurlFunAtRefCoordsOriented, numCells, numRefCoords, dim);
        {
          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisHDivCardinality, numRefCoords, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisHDivCardinality, numRefCoords, dim);
          DynRankView basisHCurlValuesAtRefCoordsCells("inValues", numCells, basisHDivCardinality, numRefCoords, dim);

          DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisHDivCardinality, numRefCoords, dim);
          basisHDiv.getValues(basisValuesAtRefCoords, refPoints);
          rst::clone(basisHCurlValuesAtRefCoordsCells,basisValuesAtRefCoords);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
              basisHCurlValuesAtRefCoordsCells,
              elemOrts,
              &basisHDiv);

          // transform basis values to the reference element (pullback)
          DynRankView ConstructWithLabel(jacobian, numCells, numRefCoords, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numRefCoords);
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
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
            *outStream << "Curl of rojection at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << curlProjFunAtRefCoordsOriented(i,j,0) << "," << curlProjFunAtRefCoordsOriented(i,j,1) << "," << curlProjFunAtRefCoordsOriented(i,j,2)   << ")";
            *outStream << "\nProjection of curl at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << projCurlFunAtRefCoordsOriented(i,j,0) << "," << projCurlFunAtRefCoordsOriented(i,j,1) << "," << projCurlFunAtRefCoordsOriented(i,j,2)  << ")";
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
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type tets[numCells][numElemVertexes];
      std::copy(&tets_orig[0][0], &tets_orig[0][0]+numCells*numElemVertexes, &tets_rotated[0][0]);

      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&tets_orig[0][0], &tets_orig[0][0]+shift, &tets_orig[0][0]+4, &tets_rotated[0][0]);

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            tets[i][j] = reorder[tets_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tet 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[0][j] << " ";
        *outStream << "] and Tet 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << tets[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
        shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, tet.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<tet.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[tets[i][j]][k];

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&tets[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);

        //compute reference points
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HVOL_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basisHVol(order-1);
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
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalDivPoints, tet.getNodeCount(), numDivPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              if(numDivPoints)
                tetLinearBasis.getValues(tetLinearBasisValuesAtEvalDivPoints, Kokkos::subview(evaluationDivPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                  for(ordinal_type j=0; j<numDivPoints; ++j)
                    physEvalDivPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalDivPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);
          ct::setJacobianDet (jacobian_det, jacobian);
          ct::setJacobianInv (jacobian_inv, jacobian);

          DynRankView ConstructWithLabel(jacobianDiv, numCells, numDivPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianDiv_det, numCells, numDivPoints);
          if(numDivPoints) {
            ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, tet);
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
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
            }
            for(int i=0;i<numDivPoints;i++) {
              targetDivAtEvalPoints(ic,i) += jacobianDiv_det(ic,i)*divFun(physEvalDivPoints(ic,i,0), physEvalDivPoints(ic,i,1), physEvalDivPoints(ic,i,2));//funHGradCoeffs(k)
            }

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
            Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(tetLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);

            for(ordinal_type i=0; i<numCells; ++i) {
              tetLinearBasis.getValues(tetLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
              for(ordinal_type d=0; d<dim; ++d) {
                for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                  for(ordinal_type j=0; j<numPoints; ++j)
                    physEvalPoints(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtEvalPoints(k,j);
                }
              }
            }
          }

          //transform the function to the reference element (inverse of pullback operator)
          DynRankView ConstructWithLabel(jacobian, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_inv, numCells, numPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, numPoints);
          ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);
          ct::setJacobianDet (jacobian_det, jacobian);

          DivFunDiv fun;
          for(int ic=0; ic<numCells; ic++) {
            for(int i=0;i<numPoints;i++)
              targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2));
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
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
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
          ct::setJacobian(jacobian, refPoints, physVertexes, tet);
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

