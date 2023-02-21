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
    target functions for Hierarchical and Nodal basis functions as N increases.
    The accuracy is computed in the H^1, H^{curl}, H^{div} and L^2 norms respectively.
    The optimal order of convergence equates the basis degree.

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

template<typename ValueType, typename DeviceType>
int ConvergenceTet(const bool verbose) {

  using ExecSpaceType = typename DeviceType::execution_space;

  typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
  typedef Kokkos::DynRankView<ordinal_type,DeviceType> DynRankViewInt;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

  *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
  *outStream << "\n";

  int errorFlag = 0;
  const ValueType relTol = 1e-4;

  struct Fun {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
    }
  };

  struct GradFun {
    KOKKOS_INLINE_FUNCTION
    ValueType
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
    KOKKOS_INLINE_FUNCTION
    ValueType
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
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType gf0[3] = {cos(x*2)*sin(y*2)*sin(z*2)*2+cos(x*y*z*8)*y*z*8, sin(x*2)*cos(y*2)*sin(z*2)*2+cos(x*y*z*8)*x*z*8, sin(x*2)*sin(y*2)*cos(z*2)*2+cos(x*y*z*8)*x*y*8};
      ValueType gf1[3] = {-sin(x*2)*cos(y*2)*cos(z*2)*2, -cos(x*2)*sin(y*2)*cos(z*2)*2, -cos(x*2)*cos(y*2)*sin(z*2)*2};
      ValueType gf2[3] = {-sin(x*y*z*8)*y*z*8, -sin(x*y*z*8)*x*z*8, -sin(x*y*z*8)*x*y*8};
      switch (comp) {
      case 0:
        return gf2[1] - gf1[2];
      case 1:
        return gf0[2] - gf2[0];
      case 2:
        return gf1[0] - gf0[1];
      default:
        return 0;
      }
    }
  };


  struct FunDiv {
    KOKKOS_INLINE_FUNCTION
    ValueType
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
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      ValueType gxf0 = cos(x*2)*sin(y*2)*sin(z*2)*2+cos(x*y*z*8)*y*z*8;
      ValueType gyf1 = -cos(x*2)*sin(y*2)*cos(z*2)*2;
      ValueType gzf2 = -sin(x*y*z*8)*x*y*8;
      return gxf0+gyf1+gzf2;
    }
  };

  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef Experimental::ProjectionTools<DeviceType> pts;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;

  constexpr ordinal_type dim = 3;
  const ordinal_type basisDegree = 3;
  ordinal_type cub_degree = 7;



  // ************************************ GET INPUTS **************************************


  int NX = 2;
  constexpr int numRefinements = 2;

  // Expected values of the projection errors in H1, Hcurl, Hdiv and L2 norms for HGRAD, HDIV, HCURL and HVOL elements respectively.
  // These values have been computed running the code with numRefinements=4 and the convergence rates are close to the optimal ones.
  // Note that these values are independent of the basis choice (Hierarchical or Nodal) as long as they generate the same functional space.
  // We currently only test two mesh refinements to make the test run faster, so this is used as a regression test rather than
  // a convergence test, but the test can be use for verifying optimal accuracy as well.
  ValueType hgradNorm[numRefinements], hcurlNorm[numRefinements], hdivNorm[numRefinements], hvolNorm[numRefinements];

  ValueType hgrad_errors[4] = {7.3278, 2.34301, 0.352119, 0.0488729};
  ValueType hcurl_errors[4] = {5.85967, 1.76969, 0.298946, 0.0404068};
  ValueType hdiv_errors[4] = {4.68569, 1.22549, 0.196687, 0.0262624};
  ValueType hvol_errors[4] = {0.424385, 0.150827, 0.0231164, 0.00301102};

  ValueType hgrad_errors_L2[4] = {7.23725, 2.20317, 0.362514, 0.049889};
  ValueType hcurl_errors_L2[4] = {6.24416, 2.00036, 0.428591, 0.090222};
  ValueType hdiv_errors_L2[4] = {4.60503, 1.2758, 0.260201, 0.0561971};
  ValueType hvol_errors_L2[4] = {0.424385, 0.150827, 0.0231164, 0.00301102};

  for(int iter= 0; iter<numRefinements; iter++, NX *= 2) {
    int NY            = NX;
    int NZ            = NX;
    bool randomMesh    = 0;  // 1 if mesh randomizer is to be used 0 if not

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base tetrahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Tetrahedron<4> >() );

    // Get dimensions
    ordinal_type numNodesPerElem = cellTopo.getNodeCount();

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
    auto hNodeCoord = Kokkos::create_mirror_view(nodeCoord);
    for (int k=0, inode = 0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
          hNodeCoord(inode,0) = leftX + (double)i*hx;
          hNodeCoord(inode,1) = leftY + (double)j*hy;
          hNodeCoord(inode,2) = leftZ + (double)k*hz;
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
            hNodeCoord(inode,0) = hNodeCoord(inode,0) + 0.125 * hx * rx;
            hNodeCoord(inode,1) = hNodeCoord(inode,1) + 0.125 * hy * ry;
            hNodeCoord(inode,2) = hNodeCoord(inode,2) + 0.125 * hz * rz;
          }
        }
      }
    }
    deep_copy(nodeCoord,hNodeCoord);

    // Element to Node map
    DynRankViewInt ConstructWithLabel(elemNodes, numElems, numNodesPerElem);
    auto hElemNodes = Kokkos::create_mirror_view(elemNodes);
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

          hElemNodes(ielem,0) = v0;
          hElemNodes(ielem,1) = v1;
          hElemNodes(ielem,2) = v2;
          hElemNodes(ielem,3) = v5;
          ielem++;

          hElemNodes(ielem,0) = v0;
          hElemNodes(ielem,1) = v2;
          hElemNodes(ielem,2) = v7;
          hElemNodes(ielem,3) = v5;
          ielem++;

          hElemNodes(ielem,0) = v0;
          hElemNodes(ielem,1) = v2;
          hElemNodes(ielem,2) = v3;
          hElemNodes(ielem,3) = v7;
          ielem++;

          hElemNodes(ielem,0) = v0;
          hElemNodes(ielem,1) = v5;
          hElemNodes(ielem,2) = v7;
          hElemNodes(ielem,3) = v4;
          ielem++;

          hElemNodes(ielem,0) = v2;
          hElemNodes(ielem,1) = v7;
          hElemNodes(ielem,2) = v5;
          hElemNodes(ielem,3) = v6;
          ielem++;
        }
      }
    }
    deep_copy(elemNodes,hElemNodes);


    //computing vertices coords
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
    //using CG_NBasis = NodalBasisFamily<DeviceType,ValueType,ValueType>;
    using CG_HBasis = HierarchicalBasisFamily<DeviceType,ValueType,ValueType>;
    using CG_DNBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;


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
      Basis_HGRAD_TET_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(basisDegree,POINTTYPE_WARPBLEND); //used only for computing reference points

      // compute orientations for cells (one time computation)
      Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numElems);
      ots::getOrientation(elemOrts, elemNodes, cellTopo);

      for (auto useL2Projection:useL2Proj) { //

      std::vector<basisType*> basis_set;
      //basis_set.push_back(new typename  CG_NBasis::HGRAD_TET(basisDegree));
      basis_set.push_back(new typename  CG_DNBasis::HGRAD_TET(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HGRAD_TET(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
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
            funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
            for(ordinal_type d=0; d<dim; ++d)
              funGradAtPhysRefCoords(i,j,d) = gradFun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2),d);
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HGRAD
        DynRankView ConstructWithLabel(basisCoeffsHGrad, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree());

          Experimental::ProjectionStruct<DeviceType,ValueType> projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHGradProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }
          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numGradPoints = projStruct.getNumTargetDerivEvalPoints();

          DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationGradPoints, numElems, numGradPoints, dim);

          if(useL2Projection) {
            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
          pts::getHGradEvaluationPoints(evaluationPoints,
              evaluationGradPoints,
              elemOrts,
              &basis,
              &projStruct);
          }


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
                auto evalPoint = Kokkos::subview(evaluationPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalGradPoint = Kokkos::subview(linearBasisValuesAtEvalGradPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numGradPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationGradPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalGradPoint, evalGradPoint);
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
              targetAtEvalPoints(ic,i) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2));
            }
            for(int i=0;i<numGradPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetGradAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*gradFun(physEvalGradPoints(ic,i,0), physEvalGradPoints(ic,i,1), physEvalGradPoints(ic,i,2), d);//funHGradCoeffs(k)
            }
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHGrad,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHGradBasisCoeffs(basisCoeffsHGrad,
                targetAtEvalPoints,
                targetGradAtEvalPoints,
                evaluationPoints,
                evaluationGradPoints,
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
        fst::HCURLtransformVALUE(transformedBasisGradsAtRefCoordsOriented,
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
      //basis_set.push_back(new typename  CG_NBasis::HCURL_TET(basisDegree));
      basis_set.push_back(new typename  CG_DNBasis::HCURL_TET(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HCURL_TET(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
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

        DynRankView ConstructWithLabel(funAtRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(funCurlAtPhysRefCoords, numElems, numRefCoords, dim);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i) {
          FunCurl fun;
          CurlFunCurl curlFun;
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<dim; ++k) {
              funAtRefCoords(i,j,k) =         fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
              funCurlAtPhysRefCoords(i,j,k) = curlFun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
            }
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HCURL
        DynRankView ConstructWithLabel(basisCoeffsHCurl, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(cub_degree),targetDerivCubDegree(cub_degree-1);

          Experimental::ProjectionStruct<DeviceType,ValueType> projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHCurlProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }
          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numCurlPoints = projStruct.getNumTargetDerivEvalPoints();
          DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationCurlPoints, numElems, numCurlPoints, dim);
          if(useL2Projection) {
            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHCurlEvaluationPoints(evaluationPoints,
                evaluationCurlPoints,
                elemOrts,
                &basis,
                &projStruct);
          }

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numElems, numCurlPoints, dim);


          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(physEvalCurlPoints, numElems, numCurlPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalCurlPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalCurlPoint = Kokkos::subview(linearBasisValuesAtEvalCurlPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numCurlPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationCurlPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalCurlPoint, evalGradPoint);
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


          DynRankView ConstructWithLabel(jacobianCurl_inv, numElems, numCurlPoints, dim, dim);
          DynRankView ConstructWithLabel(jacobianCurl_det, numElems, numCurlPoints);
          if(numCurlPoints>0){
            DynRankView ConstructWithLabel(jacobianCurl, numElems, numCurlPoints, dim, dim);
            ct::setJacobian(jacobianCurl, evaluationCurlPoints, physVertexes, cellTopo);
            ct::setJacobianInv (jacobianCurl_inv, jacobianCurl);
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
                  targetAtEvalPoints(ic,i,j) += jacobian(ic,i,d,j)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
            }
            for(int i=0;i<numCurlPoints;i++) {
              for(int d=0;d<dim;d++)
                for(int j=0;j<dim;j++)
                  targetCurlAtEvalPoints(ic,i,j) += jacobianCurl_det(ic,i)*jacobianCurl_inv(ic,i,j,d)*curlFun(physEvalCurlPoints(ic,i,0), physEvalCurlPoints(ic,i,1), physEvalCurlPoints(ic,i,2), d);//funHGradCoeffs(k)
            }
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHCurl,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
                targetAtEvalPoints,
                targetCurlAtEvalPoints,
                evaluationPoints,
                evaluationCurlPoints,
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
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const int &i, double &norm2Update) {
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type k=0; k<basisCardinality; ++k) {
                projectedFunAtRefCoords(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);
                funCurlAtRefCoordsOriented(i,j,d) += basisCoeffsHCurl(i,k)*transformedBasisCurlsAtRefCoordsOriented(i,k,j,d);
              }

              norm2Update += (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                  (funAtRefCoords(i,j,d) - projectedFunAtRefCoords(i,j,d))*
                  weights(j)*jacobianAtRefCoords_det(i,j);
              norm2Update += (funCurlAtPhysRefCoords(i,j,d) - funCurlAtRefCoordsOriented(i,j,d))*
                  (funCurlAtPhysRefCoords(i,j,d) - funCurlAtRefCoordsOriented(i,j,d))*
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
      //basis_set.push_back(new typename  CG_NBasis::HDIV_TET(basisDegree));
      basis_set.push_back(new typename  CG_DNBasis::HDIV_TET(basisDegree));
      basis_set.push_back(new typename  CG_HBasis::HDIV_TET(basisDegree));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numElems, basisCardinality, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
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
            funDivAtPhysRefCoords(i,j) = funDiv(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
            for(ordinal_type k=0; k<dim; ++k)
              funAtRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
          }
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HDIV
        DynRankView ConstructWithLabel(basisCoeffsHDiv, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree()),targetDerivCubDegree(basis.getDegree()-1);

          Experimental::ProjectionStruct<DeviceType,ValueType> projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHDivProjectionStruct(&basis, targetCubDegree, targetDerivCubDegree);
          }

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();

          DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
          DynRankView ConstructWithLabel(evaluationDivPoints, numElems, numDivPoints, dim);
          if(useL2Projection) {
            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHDivEvaluationPoints(evaluationPoints,
                evaluationDivPoints,
                elemOrts,
                &basis,
                &projStruct);
          }

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
                auto evalPoint = Kokkos::subview(evaluationPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoord(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }

              auto basisValuesAtEvalDivPoint = Kokkos::subview(linearBasisValuesAtEvalDivPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numDivPoints; ++j) {
                auto evalGradPoint = Kokkos::subview(evaluationDivPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalDivPoint, evalGradPoint);
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
                  targetAtEvalPoints(ic,i,j) += jacobian_det(ic,i)*jacobian_inv(ic,i,j,d)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2),d);
            }
            for(int i=0;i<numDivPoints;i++) {
              targetDivAtEvalPoints(ic,i) += jacobianDiv_det(ic,i)*divFun(physEvalDivPoints(ic,i,0), physEvalDivPoints(ic,i,1), physEvalDivPoints(ic,i,2));//funHGradCoeffs(k)
            }
          });
          ExecSpaceType().fence();

          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHDiv,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHDivBasisCoeffs(basisCoeffsHDiv,
                targetAtEvalPoints,
                targetDivAtEvalPoints,
                evaluationPoints,
                evaluationDivPoints,
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
      //basis_set.push_back(new typename  CG_NBasis::HVOL_TET(basisDegree-1));
      basis_set.push_back(new typename  CG_DNBasis::HVOL_TET(basisDegree-1));
      basis_set.push_back(new typename  CG_HBasis::HVOL_TET(basisDegree-1));

      for (auto basisPtr:basis_set) {
        auto& basis = *basisPtr;
        *outStream << " " << basis.getName() << std::endl;
        ordinal_type basisCardinality = basis.getCardinality();

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numElems, numRefCoords, dim);
        {
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> linearBasis; //used for computing physical coordinates
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
            funAtRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
        });
        ExecSpaceType().fence();

        // compute projection-based interpolation of fun into HVOL
        DynRankView ConstructWithLabel(basisCoeffsHVol, numElems, basisCardinality);
        {
          ordinal_type targetCubDegree(basis.getDegree());

          Experimental::ProjectionStruct<DeviceType,ValueType> projStruct;
          if(useL2Projection) {
            projStruct.createL2ProjectionStruct(&basis, targetCubDegree);
          } else {
            projStruct.createHVolProjectionStruct(&basis, targetCubDegree);
          }

          ordinal_type numPoints = projStruct.getNumTargetEvalPoints();

          DynRankView ConstructWithLabel(evaluationPoints, numElems, numPoints, dim);
          if(useL2Projection) {
            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHVolEvaluationPoints(evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          }

          DynRankView ConstructWithLabel(targetAtEvalPoints, numElems, numPoints);


          DynRankView ConstructWithLabel(physEvalPoints, numElems, numPoints, dim);
          {
            DynRankView ConstructWithLabel(linearBasisValuesAtEvalPoint, numElems, numNodesPerElem);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const int &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,i,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
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
              targetAtEvalPoints(ic,i) += jacobian_det(ic,i)*fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2));
          });
          ExecSpaceType().fence();
          if(useL2Projection) {
            pts::getL2BasisCoeffs(basisCoeffsHVol,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                &basis,
                &projStruct);
          } else {
            pts::getHVolBasisCoeffs(basisCoeffsHVol,
                targetAtEvalPoints,
                evaluationPoints,
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

