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
    \brief  Test for checking accuracy of interpolation-based projections for tetrahedral elements

    The test considers a structured tetrahedral mesh of the cube [-1,1]^3, formed by first
    building an hexahedral mesh with N^3 hexas and then splitting each hexa in 5 tetrahedra.
    The test checks the accuracy of the HGRAD, HCURL, HDIV, HVOL projections of analytic
    target functions for increasing N.
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
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceSpaceType>
int ConvergenceTet(const bool verbose) {

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
  const ValueType relTol = 1e-4;

  struct Fun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
    }
  };

  struct GradFun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      switch (comp) {
      case 0:
        return  cos(x*2)*sin(y*2)*sin(z*2)*2+cos(x*y*z*8)*y*z*8;
      case 1:
        return sin(x*2)*cos(y*2)*sin(z*2)*2+cos(x*y*z*8)*x*z*8;
      case 2:
        return sin(x*2)*sin(y*2)*cos(z*2)*2+cos(x*y*z*8)*x*y*8;
      default:
        return 0;
      }
    }
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType f0 = sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
      ValueType f1 = cos(x*2)*cos(y*2)*cos(z*2);
      ValueType f2 = cos(x*y*z*8);
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0;
      case 1:
        return f1;
      case 2:
        return f2;
      default:
        return 0;
      }
    }
  };

  struct CurlFunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType gf0[3] = {cos(x*2)*sin(y*2)*sin(z*2)*2+cos(x*y*z*8)*y*z*8, sin(x*2)*cos(y*2)*sin(z*2)*2+cos(x*y*z*8)*x*z*8, sin(x*2)*sin(y*2)*cos(z*2)*2+cos(x*y*z*8)*x*y*8};
      ValueType gf1[3] = {-sin(x*2)*cos(y*2)*cos(z*2)*2, -cos(x*2)*sin(y*2)*cos(z*2)*2, -cos(x*2)*cos(y*2)*sin(z*2)*2};
      ValueType gf2[3] = {-sin(x*y*z*8)*y*z*8, -sin(x*y*z*8)*x*z*8, -sin(x*y*z*8)*x*y*8};
      switch (comp) {
      case 0:
        return gf2[1] - gf1[2];
        break;
      case 1:
        return gf0[2] - gf2[0];
        break;
      case 2:
        return gf1[0] - gf0[1];
        break;
      default:
        return 0;
      }
    }
  };


  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType f0 = sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
      ValueType f1 = cos(x*2)*cos(y*2)*cos(z*2);
      ValueType f2 = cos(x*y*z*8);
      //fun = f + a x
      switch (comp) {
      case 0:
        return f0;
      case 1:
        return f1;
      case 2:
        return f2;
      default:
        return 0;
      }
    }
  };

  struct DivFunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      ValueType gxf0 = cos(x*2)*sin(y*2)*sin(z*2)*2+cos(x*y*z*8)*y*z*8;
      ValueType gyf1 = -cos(x*2)*sin(y*2)*cos(z*2)*2;
      ValueType gzf2 = -sin(x*y*z*8)*x*y*8;
      return gxf0+gyf1+gzf2;
    }
  };

  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef Experimental::ProjectionTools<DeviceSpaceType> pts;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;

  constexpr ordinal_type dim = 3;
  const ordinal_type order = 3;
  ordinal_type cub_degree = 7;





  // ************************************ GET INPUTS **************************************


  int NX = 2;
  constexpr int numRefinements = 2;

  // Expected values of the projection errors in H1, Hcurl, Hdiv and L2 norms for HGRAD, HDIV, HCURL and HVOL elements respectively.
  // These values have been computed running the code with numRefinements=4 and the convergence rates are close to the optimal ones.
  // We currently only test two mesh refinements to make the test run faster, so this is used as a regression test rather than
  // a convergence test, but the test can be use for verifying optimal accuracy as well.
  ValueType hgradNorm[numRefinements], hcurlNorm[numRefinements], hdivNorm[numRefinements], hvolNorm[numRefinements];

  ValueType hgrad_errors[4] = {7.3278, 2.34301, 0.352119, 0.0488729};
  ValueType hcurl_errors[4] = {5.85967, 1.76969, 0.298946, 0.0404068};
  ValueType hdiv_errors[4] = {4.68569, 1.22549, 0.196687, 0.0262624};
  ValueType hvol_errors[4] = {0.424385, 0.150827, 0.0231164, 0.00301102};

  for(int iter= 0; iter<numRefinements; iter++, NX *= 2) {
    int NY            = NX;
    int NZ            = NX;
    bool randomMesh    = 0;  // 1 if mesh randomizer is to be used 0 if not

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base tetrahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >() );

    // Get dimensions
    int numNodesPerElem = tet.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    *outStream << "Generating mesh ... \n\n";

    *outStream << "    NX" << "   NY" << "   NZ\n";
    *outStream << std::setw(5) << NX <<
        std::setw(5) << NY <<
        std::setw(5) << NZ << "\n\n";

    // Print mesh information
    int numElems = NX*NY*NZ*5;
    int numNodes = (NX+1)*(NY+1)*(NZ+1);
    int numEdges = (NX)*(NY + 1)*(NZ + 1) + (NX + 1)*(NY)*(NZ + 1) + (NX + 1)*(NY + 1)*(NZ);
    int numFaces = (NX)*(NY)*(NZ + 1) + (NX)*(NY + 1)*(NZ) + (NX + 1)*(NY)*(NZ);
    *outStream << " Number of Elements: " << numElems << " \n";
    *outStream << "    Number of Nodes: " << numNodes << " \n";
    *outStream << "    Number of Edges: " << numEdges << " \n";
    *outStream << "    Number of Faces: " << numFaces << " \n\n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

    // Get nodal coordinates
    DynRankView ConstructWithLabel(nodeCoord, numNodes, dim);
    int inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
          nodeCoord(inode,0) = leftX + (double)i*hx;
          nodeCoord(inode,1) = leftY + (double)j*hy;
          nodeCoord(inode,2) = leftZ + (double)k*hz;
          inode++;
        }
      }
    }


    // Perturb mesh coordinates (only interior nodes)
    if (randomMesh){
      for (int k=1; k<NZ; k++) {
        for (int j=1; j<NY; j++) {
          for (int i=1; i<NX; i++) {
            int inode = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
            // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0;
            // limit variation to 1/4 edge length
            nodeCoord(inode,0) = nodeCoord(inode,0) + 0.125 * hx * rx;
            nodeCoord(inode,1) = nodeCoord(inode,1) + 0.125 * hy * ry;
            nodeCoord(inode,2) = nodeCoord(inode,2) + 0.125 * hz * rz;
          }
        }
      }
    }

    // Element to Node map
    DynRankViewInt ConstructWithLabel(elemNodes, numElems, numNodesPerElem);
    int ielem = 0;
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          auto v0 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          auto v1 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          auto v2 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          auto v3 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          auto v4 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          auto v5 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          auto v6 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          auto v7 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;

          elemNodes(ielem,0) = v0;
          elemNodes(ielem,1) = v1;
          elemNodes(ielem,2) = v2;
          elemNodes(ielem,3) = v5;
          ielem++;

          elemNodes(ielem,0) = v0;
          elemNodes(ielem,1) = v2;
          elemNodes(ielem,2) = v7;
          elemNodes(ielem,3) = v5;
          ielem++;

          elemNodes(ielem,0) = v0;
          elemNodes(ielem,1) = v2;
          elemNodes(ielem,2) = v3;
          elemNodes(ielem,3) = v7;
          ielem++;

          elemNodes(ielem,0) = v0;
          elemNodes(ielem,1) = v5;
          elemNodes(ielem,2) = v7;
          elemNodes(ielem,3) = v4;
          ielem++;

          elemNodes(ielem,0) = v2;
          elemNodes(ielem,1) = v7;
          elemNodes(ielem,2) = v5;
          elemNodes(ielem,3) = v6;
          ielem++;
        }
      }
    }


    //computing vertices coords
    DynRankView ConstructWithLabel(physVertexes, numElems, tet.getNodeCount(), dim);
    for(ordinal_type i=0; i<numElems; ++i)
      for(std::size_t j=0; j<tet.getNodeCount(); ++j)
        for(ordinal_type k=0; k<dim; ++k)
          physVertexes(i,j,k) = nodeCoord(elemNodes(i,j),k);


    DefaultCubatureFactory cub_factory;
    auto cell_cub = cub_factory.create<DeviceSpaceType, ValueType, ValueType>(tet.getBaseKey(), cub_degree);
    ordinal_type numRefCoords = cell_cub->getNumPoints();
    DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
    DynRankView ConstructWithLabel(weights, numRefCoords);
    cell_cub->getCubature(refPoints, weights);




    *outStream
    << "===============================================================================\n"
    << "|                                                                             |\n"
    << "|                 Test 1 (Convergence - HGRAD)                                |\n"
    << "|                                                                             |\n"
    << "===============================================================================\n";


    try {
      //compute reference points
      Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points

      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, tet);

      Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*tetLinearBasisValuesAtRefCoords(k,j);
      }

      Fun fun;
      GradFun gradFun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
      DynRankView ConstructWithLabel(funGradAtPhysRefCoords, numElems, numRefCoords, dim);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
          for(ordinal_type d=0; d<dim; ++d)
            funGradAtPhysRefCoords(i,j,d) = gradFun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2),d);
        }
      }

      // compute projection-based interpolation of fun into HGRAD
      DynRankView ConstructWithLabel(basisCoeffsHGrad, numElems, basisCardinality);
      {
        ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree());

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createHGradProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numGradPoints = projStruct.getNumTargetDerivEvalPoints();

        DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(evaluationGradPoints, numElems, numGradPoints, dim);


        pts::getHGradEvaluationPoints(evaluationPoints,
            evaluationGradPoints,
            elemOrts,
            &basis,
            &projStruct);


        DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints);
        DynRankView ConstructWithLabel(targetGradAtEvalPoints, numElems, numGradPoints, dim);

        DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(physEvalGradPoints, numElems, numGradPoints, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalGradPoints, tet.getNodeCount(), numGradPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            if(numGradPoints>0)
              hexLinearBasis.getValues(hexLinearBasisValuesAtEvalGradPoints, Kokkos::subview(evaluationGradPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numGradPoints; ++j)
                  physEvalGradPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalGradPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numGradPoints, dim, dim);
        if(numGradPoints>0)
          ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, tet);

        GradFun gradFun;
        Kokkos::deep_copy(targetGradAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
      ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
      ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
      fst::HCURLtransformVALUE(transformedBasisGradsAtRefCoordsOriented,
          jacobianAtRefCoords_inv,
          basisGradsAtRefCoordsOriented);

      DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords);
      DynRankView ConstructWithLabel(funGradAtRefCoordsOriented, numElems, numRefCoords,dim);

      //compute error of projection in H1 norm
      ValueType norm2(0);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<basisCardinality; ++k) {
            projectedFunAtRefCoords(i,j) += basisCoeffsHGrad(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
            for (ordinal_type d=0; d<dim; ++d)
              funGradAtRefCoordsOriented(i,j,d) += basisCoeffsHGrad(i,k)*transformedBasisGradsAtRefCoordsOriented(i,k,j,d);
          }
          norm2 += std::pow(funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j),2)*weights(j)*jacobianAtRefCoords_det(i,j);
          for (ordinal_type d=0; d<dim; ++d)
            norm2 += std::pow(funGradAtPhysRefCoords(i,j,d) - funGradAtRefCoordsOriented(i,j,d),2)*weights(j)*jacobianAtRefCoords_det(i,j);
        }
      }
      hgradNorm[iter] =  std::sqrt(norm2);
      auto expected_error = hgrad_errors[iter];
      if(std::abs(hgradNorm[iter]-expected_error)/expected_error > relTol){
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "For N = " << NX << ", computed error (" << hgradNorm[iter] << ") is different than expected one (" << expected_error << ")";
        *outStream << std::endl;
      }
      *outStream << "HGRAD Error: " << hgradNorm[iter] <<std::endl;

    } catch (std::exception err) {
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
      Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, tet);

      Basis_HCURL_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*tetLinearBasisValuesAtRefCoords(k,j);
      }

      //check function reproducibility
      FunCurl fun;
      CurlFunCurl curlFun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(funCurlAtPhysRefCoords, numElems, numRefCoords, dim);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<dim; ++k) {
            funAtRefCoords(i,j,k) =         fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
            funCurlAtPhysRefCoords(i,j,k) = curlFun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
          }
        }
      }

      // compute projection-based interpolation of fun into HCURL
      DynRankView ConstructWithLabel(basisCoeffsHCurl, numElems, basisCardinality);
      {
        ordinal_type targetCubDegree(cub_degree),targetDerivCubDegree(cub_degree-1);

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createHCurlProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numCurlPoints = projStruct.getNumTargetDerivEvalPoints();
        DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(evaluationCurlPoints, numElems, numCurlPoints, dim);
        pts::getHCurlEvaluationPoints(evaluationPoints,
            evaluationCurlPoints,
            elemOrts,
            &basis,
            &projStruct);


        DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numElems, numCurlPoints, dim);


        DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(physEvalCurlPoints, numElems, numCurlPoints, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalCurlPoints, tet.getNodeCount(), numCurlPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalCurlPoints, Kokkos::subview(evaluationCurlPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numCurlPoints; ++j)
                  physEvalCurlPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalCurlPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);

        DynRankView ConstructWithLabel(jacobianCurl, numElems, numCurlPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobianCurl_inv, numElems, numCurlPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobianCurl_det, numElems, numCurlPoints);
        ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, tet);
        ct::setJacobianInv (jacobianCurl_inv, jacobianCurl);
        ct::setJacobianDet (jacobianCurl_det, jacobianCurl);

        CurlFunCurl curlFun;
        Kokkos::deep_copy(targetCurlAtEvalPoints,0.);
        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
      ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
      ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
      fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
          jacobianAtRefCoords_inv,
          basisValuesAtRefCoordsOriented);


      DynRankView ConstructWithLabel(basisCurlsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
      DynRankView ConstructWithLabel(transformedBasisCurlsAtRefCoordsOriented, numElems, basisCardinality, numRefCoords, dim);
      DynRankView basisCurlsAtRefCoordsCells("inValues", numElems, basisCardinality, numRefCoords, dim);

      DynRankView ConstructWithLabel(basisCurlsAtRefCoords, basisCardinality, numRefCoords, dim);
      basis.getValues(basisCurlsAtRefCoords, refPoints,OPERATOR_CURL);
      rst::clone(basisCurlsAtRefCoordsCells,basisCurlsAtRefCoords);

      // modify basis values to account for orientations
      ots::modifyBasisByOrientation(basisCurlsAtRefCoordsOriented,
          basisCurlsAtRefCoordsCells,
          elemOrts,
          &basis);


      fst::HCURLtransformCURL(transformedBasisCurlsAtRefCoordsOriented,
          jacobianAtRefCoords,
          jacobianAtRefCoords_det,
          basisCurlsAtRefCoordsOriented);


      DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(funCurlAtRefCoordsOriented, numElems, numRefCoords,dim);

      //compute error of projection in HCURL norm
      ValueType norm2(0);
      for(ordinal_type i=0; i<numElems; ++i) {
        ValueType error=0;
        for(ordinal_type j=0; j<numRefCoords; ++j)
          for(ordinal_type d=0; d<dim; ++d) {
            for(ordinal_type k=0; k<basisCardinality; ++k) {
              projectedFunAtRefCoords(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
              funCurlAtRefCoordsOriented(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisCurlsAtRefCoordsOriented(i,k,j,d);
            }

            norm2 += std::pow(funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d),2)*weights(j)*jacobianAtRefCoords_det(i,j);
            norm2 += std::pow(funCurlAtPhysRefCoords(i,j,d) - funCurlAtRefCoordsOriented(i,j,d),2)*weights(j)*jacobianAtRefCoords_det(i,j);
            error = std::max(std::abs( funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d)), error);
          }
      }
      hcurlNorm[iter] =  std::sqrt(norm2);
      auto expected_error = hcurl_errors[iter];
      if(std::abs(hcurlNorm[iter]-expected_error)/expected_error > relTol){
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "For N = " << NX << ", computed error (" << hcurlNorm[iter] << ") is different than expected one (" << expected_error << ")";
        *outStream << std::endl;
      }
      *outStream << "HCURL Error: " << hcurlNorm[iter] <<std::endl;
    } catch (std::exception err) {
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
      Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, tet);

      Basis_HDIV_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(physDofCoords, numElems, basisCardinality, dim);
      {
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*tetLinearBasisValuesAtRefCoords(k,j);
      }

      FunDiv fun;
      DivFunDiv funDiv;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(funDivAtPhysRefCoords, numElems, numRefCoords);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          funDivAtPhysRefCoords(i,j) = funDiv(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
          for(ordinal_type k=0; k<dim; ++k)
            funAtRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
        }
      }

      // compute projection-based interpolation of fun into HDIV
      DynRankView ConstructWithLabel(basisCoeffsHDiv, numElems, basisCardinality);
      {
        ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree()-1);

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createHDivProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();

        DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(evaluationDivPoints, numElems, numDivPoints, dim);

        pts::getHDivEvaluationPoints(evaluationPoints,
            evaluationDivPoints,
            elemOrts,
            &basis,
            &projStruct);

        DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(targetDivAtEvalPoints, numElems, numDivPoints);


        DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(physEvalDivPoints, numElems, numDivPoints, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalDivPoints, tet.getNodeCount(), numDivPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalDivPoints, Kokkos::subview(evaluationDivPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numDivPoints; ++j)
                  physEvalDivPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalDivPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
        DynRankView ConstructWithLabel(jacobian_inv, numElems, numPoints, dim, dim);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);
        ct::setJacobianDet (jacobian_det, jacobian);
        ct::setJacobianInv (jacobian_inv, jacobian);

        DynRankView ConstructWithLabel(jacobianDiv, numElems, numDivPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobianDiv_det, numElems, numDivPoints);
        ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, tet);
        ct::setJacobianDet (jacobianDiv_det, jacobianDiv);


        DivFunDiv divFun;
        Kokkos::deep_copy(targetDivAtEvalPoints,0.);
        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
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
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<basisCardinality; ++k) {
            for(ordinal_type d=0; d<dim; ++d)
              projectedFunAtRefCoords(i,j,d) += basisCoeffsHDiv(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            funDivAtRefCoordsOriented(i,j) += basisCoeffsHDiv(i,k)*transformedBasisDivsAtRefCoordsOriented(i,k,j);
          }

          for(ordinal_type d=0; d<dim; ++d) {
            norm2 += std::pow(funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d),2)*weights(j)*jacobianAtRefCoords_det(i,j);
          }
          norm2 += std::pow(funDivAtPhysRefCoords(i,j) - funDivAtRefCoordsOriented(i,j),2)*weights(j)*jacobianAtRefCoords_det(i,j);
        }
      }
      hdivNorm[iter] = std::sqrt(norm2);
      auto expected_error = hdiv_errors[iter];
      if(std::abs(hdivNorm[iter]-expected_error)/expected_error > relTol){
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "For N = " << NX << ", computed error (" << hdivNorm[iter] << ") is different than expected one (" << expected_error << ")";
        *outStream << std::endl;
      }
      *outStream << "HDIV Error: " << hdivNorm[iter] <<std::endl;

    } catch (std::exception err) {
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
      Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, tet);

      Basis_HVOL_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order-1);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*tetLinearBasisValuesAtRefCoords(k,j);
      }

      //check function reproducibility
      Fun fun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j)
          funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
      }

      // compute projection-based interpolation of fun into HVOL
      DynRankView ConstructWithLabel(basisCoeffsHVol, numElems, basisCardinality);
      {
        ordinal_type targetCubDegree(basis.getDegree());

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createHVolProjectionStruct(&basis, targetCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints();

        DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);

        pts::getHVolEvaluationPoints(evaluationPoints,
            elemOrts,
            &basis,
            &projStruct);

        DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints);


        DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtEvalPoints, tet.getNodeCount(), numPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            hexLinearBasis.getValues(hexLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<tet.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*hexLinearBasisValuesAtEvalPoints(k,j);
              }
            }
          }
        }

        //transform the target function to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, tet);
        ct::setJacobianDet (jacobian_det, jacobian);

        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
          for(int i=0;i<numPoints;i++)
            targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2));
        }

        pts::getHVolBasisCoeffs(basisCoeffsHVol,
            targetAtEvalPoints,
            evaluationPoints,
            elemOrts,
            &basis,
            &projStruct);
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
      ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
      fst::HVOLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
          jacobianAtRefCoords_det,
          basisValuesAtRefCoordsOriented);

      DynRankView ConstructWithLabel(projectedFunAtRefCoords, numElems, numRefCoords);

      //compute error of projection in L2 norm
      ValueType norm2(0);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<basisCardinality; ++k)
            projectedFunAtRefCoords(i,j) += basisCoeffsHVol(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);
          norm2 += std::pow(funAtRefCoords(i,j) - projectedFunAtRefCoords(i,j),2)*weights(j)*jacobianAtRefCoords_det(i,j);
        }
      }
      hvolNorm[iter] =  std::sqrt(norm2);
      auto expected_error = hvol_errors[iter];
      if(std::abs(hvolNorm[iter]-expected_error)/expected_error > relTol){
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "For N = " << NX << ", computed error (" << hvolNorm[iter] << ") is different than expected one (" << expected_error << ")";
        *outStream << std::endl;
      }
      *outStream << "HVOL Error: " << hvolNorm[iter] <<std::endl;

    } catch (std::exception err) {
      std::cout << " Exeption\n";
      *outStream << err.what() << "\n\n";
      errorFlag = -1000;
    }
  }

  std::cout<< "\nHGRAD ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    std::cout << " " << hgradNorm[iter];
  std::cout<< "\nHCURL ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    std::cout << " " << hcurlNorm[iter];
  std::cout<< "\nHDIV ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    std::cout << " " << hdivNorm[iter];
  std::cout<< "\nHVOL ERROR:";
  for(int iter = 0; iter<numRefinements; iter++)
    std::cout << " " << hvolNorm[iter];
  std::cout<<std::endl;

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

