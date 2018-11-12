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
    \brief  Test for checking accuracy of interpolation-based projections for quadrilateral elements

    The test considers a uniform and structured quadrilateral mesh of the square [-1,1]^2, formed by
    N^2 quads, and checks the accuracy of the HGRAD, HCURL, HDIV, HVOL projections of analytic
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
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
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
int ConvergenceQuad(const bool verbose) {

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
    operator()(const ValueType& x, const ValueType& y) {
      return sin(x*2)*sin(y*2)+sin(x*y*8);
    }
  };

  struct GradFun {
    ValueType
    KOKKOS_INLINE_FUNCTION
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
    ValueType
    KOKKOS_INLINE_FUNCTION
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
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      ValueType gf0[2] = {cos(x*2)*sin(y*2)*2+cos(x*y*8)*y*8, sin(x*2)*cos(y*2)*2+cos(x*y*8)*x*8};
      ValueType gf1[2] = {-sin(x*2)*cos(y*2)*2, -cos(x*2)*sin(y*2)*2};
      return gf1[0] - gf0[1];
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
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
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      ValueType gxf0 = cos(x*2)*sin(y*2)*2+cos(x*y*8)*y*8;
      ValueType gyf1 = -cos(x*2)*sin(y*2)*2;
      return gxf0+gyf1;
    }
  };

  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef Experimental::ProjectionTools<DeviceSpaceType> pts;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;

  constexpr ordinal_type dim = 2;
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

  ValueType hgrad_errors[4] = {4.06298, 1.12792, 0.127858, 0.0164394};
  ValueType hcurl_errors[4] = {3.06114, 0.976496, 0.114555, 0.0150085};
  ValueType hdiv_errors[4] = {4.11097, 0.997027, 0.115616, 0.0150505};
  ValueType hvol_errors[4] = {0.896361, 0.104481, 0.0189055, 0.00239322};

  for(int iter= 0; iter<numRefinements; iter++, NX *= 2) {
    int NY            = NX;
    bool randomMesh    = 0;  // 1 if mesh randomizer is to be used 0 if not

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base quadrilateral
    typedef shards::CellTopology    CellTopology;
    CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Get dimensions
    int numNodesPerElem = quad.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    *outStream << "Generating mesh ... \n\n";

    *outStream << "    NX" << "   NY\n";
    *outStream << std::setw(5) << NX <<
        std::setw(5) << NY << "\n\n";

    // Print mesh information
    int numElems = NX*NY;
    int numNodes = (NX+1)*(NY+1);
    *outStream << " Number of Elements: " << numElems << " \n";
    *outStream << "    Number of Nodes: " << numNodes << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);

    // Get nodal coordinates
    DynRankView ConstructWithLabel(nodeCoord, numNodes, dim);
    int inode = 0;
    for (int j=0; j<NY+1; j++) {
      for (int i=0; i<NX+1; i++) {
        nodeCoord(inode,0) = leftX + (double)i*hx;
        nodeCoord(inode,1) = leftY + (double)j*hy;
        inode++;
      }
    }


    // Perturb mesh coordinates (only interior nodes)
    if (randomMesh){
      for (int j=1; j<NY; j++) {
        for (int i=1; i<NX; i++) {
          int inode = i + j * (NX + 1);
          // random numbers between -1.0 and 1.0
          double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
          double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
          // limit variation to 1/4 edge length
          nodeCoord(inode,0) = nodeCoord(inode,0) + 0.125 * hx * rx;
          nodeCoord(inode,1) = nodeCoord(inode,1) + 0.125 * hy * ry;
        }
      }
    }

    // Element to Node map
    DynRankViewInt ConstructWithLabel(elemNodes, numElems, numNodesPerElem);
    int ielem = 0;

    for (int j=0; j<NY; j++) {
      for (int i=0; i<NX; i++) {
        elemNodes(ielem,0) = (NX + 1)*j + i;
        elemNodes(ielem,1) = (NX + 1)*j + i + 1;
        elemNodes(ielem,2) = (NX + 1)*(j + 1) + i + 1;
        elemNodes(ielem,3) = (NX + 1)*(j + 1) + i;
        ielem++;
      }
    }


    //computing vertices coords
    DynRankView ConstructWithLabel(physVertexes, numElems, quad.getNodeCount(), dim);
    for(ordinal_type i=0; i<numElems; ++i) {
      for(std::size_t j=0; j<quad.getNodeCount(); ++j)
        for(ordinal_type k=0; k<dim; ++k)
          physVertexes(i,j,k) = nodeCoord(elemNodes(i,j),k);
    }



    DefaultCubatureFactory cub_factory;
    auto cell_cub = cub_factory.create<DeviceSpaceType, ValueType, ValueType>(quad.getBaseKey(), cub_degree);
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
      Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points

      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, quad);

      Basis_HGRAD_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
        quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtRefCoords(k,j);
      }

      Fun fun;
      GradFun gradFun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
      DynRankView ConstructWithLabel(funGradAtPhysRefCoords, numElems, numRefCoords, dim);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
          for(ordinal_type d=0; d<dim; ++d)
            funGradAtPhysRefCoords(i,j,d) = gradFun(physRefCoords(i,j,0), physRefCoords(i,j,1),d);
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
          Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalGradPoints, quad.getNodeCount(), numGradPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            if(numGradPoints>0)
              quadLinearBasis.getValues(quadLinearBasisValuesAtEvalGradPoints, Kokkos::subview(evaluationGradPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numGradPoints; ++j)
                  physEvalGradPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalGradPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numGradPoints, dim, dim);
        if(numGradPoints>0)
          ct::setJacobian(jacobian, evaluationGradPoints, physVertexes, quad);

        GradFun gradFun;
        Kokkos::deep_copy(targetGradAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
          for(int i=0;i<numPoints;i++) {
            targetAtEvalPoints(ic,i) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
          }
          for(int i=0;i<numGradPoints;i++) {
            for(int d=0;d<dim;d++)
              for(int j=0;j<dim;j++)
                targetGradAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*gradFun(physEvalGradPoints(ic,i,0), physEvalGradPoints(ic,i,1), d);
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
      ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
      ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
      fst::HGRADtransformGRAD(transformedBasisGradsAtRefCoordsOriented,
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
      ots::getOrientation(elemOrts, elemNodes, quad);

      Basis_HCURL_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
        quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtRefCoords(k,j);
      }

      //check function reproducibility
      FunCurl fun;
      CurlFunCurl curlFun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(funCurlAtPhysRefCoords, numElems, numRefCoords);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<dim; ++k)
            funAtRefCoords(i,j,k) =         fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
          funCurlAtPhysRefCoords(i,j) = curlFun(physRefCoords(i,j,0), physRefCoords(i,j,1));
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
        DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numElems, numCurlPoints);


        DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
        DynRankView ConstructWithLabel(physEvalCurlPoints, numElems, numCurlPoints, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalCurlPoints, quad.getNodeCount(), numCurlPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalCurlPoints, Kokkos::subview(evaluationCurlPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numCurlPoints; ++j)
                  physEvalCurlPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalCurlPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);

        DynRankView ConstructWithLabel(jacobianCurl, numElems, numCurlPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobianCurl_det, numElems, numCurlPoints);
        ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, quad);
        ct::setJacobianDet (jacobianCurl_det, jacobianCurl);

        CurlFunCurl curlFun;
        Kokkos::deep_copy(targetCurlAtEvalPoints,0.);
        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
          for(int i=0;i<numPoints;i++) {
            for(int j=0;j<dim;j++)
              for(int d=0;d<dim;d++)
                targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
          }
          for(int i=0;i<numCurlPoints;i++)
            targetCurlAtEvalPoints(ic,i) += jacobianCurl_det(ic,i)*curlFun(physEvalCurlPoints(ic,i,0), physEvalCurlPoints(ic,i,1));
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
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
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for(ordinal_type k=0; k<basisCardinality; ++k) {
            for(ordinal_type d=0; d<dim; ++d)
              projectedFunAtRefCoords(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
            funCurlAtRefCoordsOriented(i,j) += basisCoeffsHCurl(i,k)*transformedBasisCurlsAtRefCoordsOriented(i,k,j);
          }
          for(ordinal_type d=0; d<dim; ++d)
            norm2 += std::pow(funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d),2)*weights(j)*jacobianAtRefCoords_det(i,j);

          norm2 += std::pow(funCurlAtPhysRefCoords(i,j) - funCurlAtRefCoordsOriented(i,j),2)*weights(j)*jacobianAtRefCoords_det(i,j);
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
      ots::getOrientation(elemOrts, elemNodes, quad);

      Basis_HDIV_QUAD_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(physDofCoords, numElems, basisCardinality, dim);
      {
        Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
        quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtRefCoords(k,j);
      }

      FunDiv fun;
      DivFunDiv funDiv;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
      DynRankView ConstructWithLabel(funDivAtPhysRefCoords, numElems, numRefCoords);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          funDivAtPhysRefCoords(i,j) = funDiv(physRefCoords(i,j,0), physRefCoords(i,j,1));
          for(ordinal_type k=0; k<dim; ++k)
            funAtRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
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
          Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalDivPoints, quad.getNodeCount(), numDivPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalDivPoints, Kokkos::subview(evaluationDivPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalPoints(k,j);
                for(ordinal_type j=0; j<numDivPoints; ++j)
                  physEvalDivPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalDivPoints(k,j);
              }
            }
          }
        }

        //transform the target function and its derivative to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
        DynRankView ConstructWithLabel(jacobian_inv, numElems, numPoints, dim, dim);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
        ct::setJacobianDet (jacobian_det, jacobian);
        ct::setJacobianInv (jacobian_inv, jacobian);

        DynRankView ConstructWithLabel(jacobianDiv, numElems, numDivPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobianDiv_det, numElems, numDivPoints);
        ct::setJacobian(jacobianDiv, evaluationDivPoints, physVertexes, quad);
        ct::setJacobianDet (jacobianDiv_det, jacobianDiv);


        DivFunDiv divFun;
        Kokkos::deep_copy(targetDivAtEvalPoints,0.);
        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
          for(int i=0;i<numPoints;i++) {
            for(int j=0;j<dim;j++)
              for(int d=0;d<dim;d++)
                targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1),d);
          }
          for(int i=0;i<numDivPoints;i++) {
            targetDivAtEvalPoints(ic,i) += jacobianDiv_det(ic,i)*divFun(physEvalDivPoints(ic,i,0), physEvalDivPoints(ic,i,1));
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
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
      ots::getOrientation(elemOrts, elemNodes, quad);

      Basis_HVOL_QUAD_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order-1);
      ordinal_type basisCardinality = basis.getCardinality();

      //Compute physical Dof Coordinates and Reference coordinates
      DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
      {
        Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
        DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
        quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
        for(ordinal_type i=0; i<numElems; ++i)
          for(ordinal_type d=0; d<dim; ++d)
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                physRefCoords(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtRefCoords(k,j);
      }

      //check function reproducibility
      Fun fun;
      DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords);
      for(ordinal_type i=0; i<numElems; ++i) {
        for(ordinal_type j=0; j<numRefCoords; ++j)
          funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
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
          Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtEvalPoints, quad.getNodeCount(), numPoints);

          for(ordinal_type i=0; i<numElems; ++i) {
            quadLinearBasis.getValues(quadLinearBasisValuesAtEvalPoints, Kokkos::subview(evaluationPoints,i,Kokkos::ALL(),Kokkos::ALL()));
            for(ordinal_type d=0; d<dim; ++d) {
              for(std::size_t k=0; k<quad.getNodeCount(); ++k) {
                for(ordinal_type j=0; j<numPoints; ++j)
                  physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*quadLinearBasisValuesAtEvalPoints(k,j);
              }
            }
          }
        }

        //transform the target function to the reference element (inverse of pullback operator)
        DynRankView ConstructWithLabel(jacobian, numElems, numPoints, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numElems, numPoints);
        ct::setJacobian(jacobian, evaluationPoints, physVertexes, quad);
        ct::setJacobianDet (jacobian_det, jacobian);

        Kokkos::deep_copy(targetAtEvalPoints,0.);
        for(int ic=0; ic<numElems; ic++) {
          for(int i=0;i<numPoints;i++)
            targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1));
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
      ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
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

