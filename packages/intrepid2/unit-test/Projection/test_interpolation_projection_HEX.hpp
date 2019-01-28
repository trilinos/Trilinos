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
    \brief Test interpolation and projection capabilities for Hexaedral elements

    The test considers two hexahedra in the physical space sharing a common face. 
    In order to test significant configurations, we consider 6 mappings of the reference hexahedron
    to the first (physical) hexahedron, so that the common face is mapped from all the 6 faces
    of the reference hexahedron.
    Then, for each of the mappings, the global ids of the vertices of the common face are permuted.

    The test considers HGRAD, HCURL, HDIV and HVOL, of different degree, and for each of them checks that
    the Lagrangian interpolation, the interpolation-based projection, and the L2 projection, reproduce the
    target function to round-off errors when the target function is in the corresponding finite element space.

    Also, for the Lagrangian Interpolations, it checks that:
    1. that the Kronecker property holds for the oriented basis evaluated at the oriented DOF coordinates.
    2. that the basis coefficients located at the common faces/edges are the same when computed on the
       first and second hexahedron.

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"


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
int InterpolationProjectionHex(const bool verbose) {

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
    operator()(const ordinal_type& degree, const ValueType& x, const ValueType& y, const ValueType& z) {
      return std::pow(x+y-z,degree-1)*(x+y-2);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ordinal_type& degree, const ValueType& x, const ValueType& y, const ValueType& z, const int comp) {
      ValueType a = std::pow(x-y+z, degree-1);
      ValueType f0 = 3;
      ValueType f1 = std::pow(y, degree-1);
      ValueType f2 = std::pow(x+z, degree-1);
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
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ordinal_type& degree, const ValueType& x, const ValueType& y, const ValueType& z, const int comp) {
      ValueType a0 = std::pow(x-y+z, degree-1);
      ValueType a1 = std::pow(2-y+z, degree-1);
      ValueType a2 = std::pow(x-1, degree-1);
      ValueType f0 = 3;
      ValueType f1 = std::pow(y, degree-1);
      ValueType f2 = std::pow(x+z, degree-1);
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 + (a1*z-a2*y);
      case 1:
        return f1 + (a2*x-a0*z);
      case 2:
        return f2 + (a0*y-a1*x);
      default:
        return 0;
      }
    }
  };

  typedef std::array<ordinal_type,2> edgeType;
  typedef std::array<ordinal_type,4> faceType;
  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef Experimental::ProjectionTools<DeviceSpaceType> pts;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;
  typedef Experimental::LagrangianInterpolation<DeviceSpaceType> li;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 8;
  constexpr ordinal_type numTotalVertexes = 12;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}, {-1,-1,2},{1,-1,2},{1,1,2},{-1,1,2}};
  ordinal_type hexas_orig[numCells][numElemVertexes] = {{0,1,2,3,4,5,6,7},{4,5,6,7,8,9,10,11}};  //common face is {4,5,6,7}
  faceType common_face = {{4,5,6,7}};
  faceType faceLeft = {{0, 3, 7, 4}};
  faceType faceRight = {{1, 2, 6, 5}};
  faceType faceFront = {{0, 4, 5, 1}};
  faceType faceBack = {{2, 3, 7, 6}};
  ordinal_type hexas_rotated[numCells][numElemVertexes];
  faceType faceLeftOriented, faceRightOriented, faceBackOriented, faceFrontOriented;

  std::set<edgeType> common_edges;
  common_edges.insert(edgeType({{4,5}})); common_edges.insert(edgeType({{5,6}})); common_edges.insert(edgeType({{6,7}})); common_edges.insert(edgeType({{4,7}}));
  const ordinal_type max_degree = 4;

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 1 (Orientation - HGRAD)                                |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    //reordering of nodes to explore different orientations

    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);
      for (ordinal_type shift=0; shift<1; ++shift) {
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceLeft[ii]] = hexas_orig[0][faceLeftOriented[ii]];
            hexas_rotated[0][faceRight[ii]] = hexas_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceFront[ii]] = hexas_orig[0][faceFrontOriented[ii]];
            hexas_rotated[0][faceBack[ii]] = hexas_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            hexas[i][j] = reorder[hexas_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hex.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hex.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];

        //computing common face and edges
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][4];

        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<hex.getSideCount(); ++is) {
              for (std::size_t k=0; k<hex.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hex.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<hex.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<hex.getNodeCount(1,ie); ++k)
                edge[k]= hexas_rotated[i][hex.getNodeMap(1,ie,k)];
              std::sort(edge.begin(),edge.end());
              auto it=common_edges.find(edge);
              if(it !=common_edges.end()){
                auto edge_lid = std::distance(common_edges.begin(),it);
                edgeIndexes[i][edge_lid]=ie;
              }
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hex);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {

          Teuchos::RCP<Basis<DeviceSpaceType,ValueType,ValueType>> basisPtr;

          if(degree==1)
            basisPtr = Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,ValueType,ValueType>());
          else
            basisPtr = Teuchos::rcp(new Basis_HGRAD_HEX_Cn_FEM<DeviceSpaceType,ValueType,ValueType>(degree,POINTTYPE_WARPBLEND));

          ordinal_type basisCardinality = basisPtr->getCardinality();

          //compute DofCoords Oriented
          DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality);
          DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality);
          DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

          //compute Lagrangian Interpolation of fun
          {
            li::getDofCoordsAndCoeffs(dofCoordsOriented,  dofCoeffsPhys, basisPtr.get(), POINTTYPE_WARPBLEND, elemOrts);

            //Compute physical Dof Coordinates
            {
              Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
              DynRankView ConstructWithLabel(hexLinearBasisValuesAtDofCoords, numCells, hex.getNodeCount(), basisCardinality);
              for(ordinal_type i=0; i<numCells; ++i)
                for(ordinal_type d=0; d<dim; ++d) {
                  auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                  auto outView =Kokkos::subview( hexLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
                  hexLinearBasis.getValues(outView, inView);

                  for(ordinal_type j=0; j<basisCardinality; ++j)
                    for(std::size_t k=0; k<hex.getNodeCount(); ++k)
                      physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexLinearBasisValuesAtDofCoords(i,k,j);
                }
            }

            Fun fun;

            for(ordinal_type i=0; i<numCells; ++i) {
              for(ordinal_type j=0; j<basisCardinality; ++j)
                funAtDofCoords(i,j) += fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2));
            }

            li::getBasisCoeffs(basisCoeffsLI, funAtDofCoords, dofCoeffsPhys);
          }

          //Testing Kronecker property of basis functions
          {
            for(ordinal_type i=0; i<numCells; ++i) {
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
              DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL());
              basisPtr->getValues(outView, inView);

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr.get());

              // transform basis values
              fst::HGRADtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoordsOriented);


              for(ordinal_type k=0; k<basisCardinality; ++k) {
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  ValueType dofValue = transformedBasisValuesAtDofCoordsOriented(i,k,j) * dofCoeffsPhys(i,j);
                  if ( k==j && std::abs( dofValue - 1.0 ) > 100*tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not have unit value at its node (" << dofValue <<")\n";
                  }
                  if ( k!=j && std::abs( dofValue ) > tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
                  }
                }
              }
            }
          }

          //check that fun values are consistent on common face dofs
          {
            bool areDifferent(false);
            auto numFaceDOFs = basisPtr->getDofCount(2,0);
            for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
              areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                  - basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
            }

            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
              *outStream << "Function DOFs for Hex 0 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
              *outStream << "\nFunction DOFs for Hex 1 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
              *outStream << std::endl;
            }
          }

          //check that fun values are consistent on common edges dofs
          {
            bool areDifferent(false);
            auto numEdgeDOFs = basisPtr->getDofCount(1,0);
            for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
              for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
                areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                    - basisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
              }
              if(areDifferent) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                *outStream << "Function DOFs for Hex 0 are:";
                for(ordinal_type j=0;j<numEdgeDOFs;j++)
                  *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j));
                *outStream << "\nFunction DOFs for Hex 1 are:";
                for(ordinal_type j=0;j<numEdgeDOFs;j++)
                  *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j));
                *outStream << std::endl;
              }
            }
          }

          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
          DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
          DynRankView basisValuesAtDofCoordsCells("inValues", numCells, basisCardinality, basisCardinality);

          for (ordinal_type ic = 0; ic < numCells; ++ic)
            basisPtr->getValues(Kokkos::subview(basisValuesAtDofCoordsCells, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(dofCoordsOriented, ic, Kokkos::ALL(), Kokkos::ALL()));

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoordsCells,
              elemOrts,
              basisPtr.get());

          // transform basis values
          // transform basis values
          fst::HGRADtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
              basisValuesAtDofCoordsOriented);

          DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i) {
            ValueType error=0;
            for(ordinal_type j=0; j<basisCardinality; ++j) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                funAtDofCoordsOriented(i,j) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j);

              error = std::max(std::abs( funAtDofCoords(i,j) - funAtDofCoordsOriented(i,j)), error);
            }

            if(error>100*tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
              *outStream << "Function values at reference points are:\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoords(i,j)  << ")";
              *outStream << "\nFunction values at reference points computed using basis functions are\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoordsOriented(i,j)  << ")";
              *outStream << std::endl;
            }
          }


          //compute projection-based interpolation of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsHGrad, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree());

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createHGradProjectionStruct(basisPtr.get(), targetCubDegree, targetDerivCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numGradPoints = projStruct.getNumTargetDerivEvalPoints();
            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
            DynRankView ConstructWithLabel(evaluationGradPoints, numCells, numGradPoints, dim);


            pts::getHGradEvaluationPoints(evaluationPoints,
                evaluationGradPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);


            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
            DynRankView ConstructWithLabel(targetGradAtEvalPoints, numCells, numGradPoints, dim);

            DynRankView ConstructWithLabel(hgradBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
            DynRankView ConstructWithLabel(hgradBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints);
            for(int ic=0; ic<numCells; ic++)
              basisPtr->getValues(Kokkos::subview(hgradBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hgradBasisAtEvaluationPoints,
                hgradBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());

            DynRankView ConstructWithLabel(gradOfHGradBasisAtEvaluationPoints, numCells, basisCardinality , numGradPoints, dim);
            if(numGradPoints>0) {
              DynRankView ConstructWithLabel(gradOfHGradBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numGradPoints, dim);
              for(int ic=0; ic<numCells; ic++)
                basisPtr->getValues(Kokkos::subview(gradOfHGradBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationGradPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_GRAD);
              ots::modifyBasisByOrientation(gradOfHGradBasisAtEvaluationPoints,
                  gradOfHGradBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr.get());
            }


            for(int ic=0; ic<numCells; ic++) {
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hgradBasisAtEvaluationPoints(ic,k,i);
              }
              for(int i=0;i<numGradPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetGradAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*gradOfHGradBasisAtEvaluationPoints(ic,k,i,d);//funHGradCoeffs(k)
              }

            }

            pts::getHGradBasisCoeffs(basisCoeffsHGrad,
                targetAtEvalPoints,
                targetGradAtEvalPoints,
                evaluationPoints,
                evaluationGradPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian nterpolation are the same as those of the projection-based interpolation
          {
            ValueType diffErr(0);

            for(int k=0;k<basisCardinality;k++) {
              for(int ic=0; ic<numCells; ic++)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsHGrad(ic,k)));
            }

            if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HGRAD_C" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
            }
          }

          //compute L2 projection of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree());

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createL2ProjectionStruct(basisPtr.get(), targetCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints();
            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);

            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);


            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
            DynRankView ConstructWithLabel(hgradBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
            DynRankView ConstructWithLabel(hgradBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints);
            for(int ic=0; ic<numCells; ic++)
              basisPtr->getValues(Kokkos::subview(hgradBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hgradBasisAtEvaluationPoints,
                hgradBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());


            for(int ic=0; ic<numCells; ic++) {
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hgradBasisAtEvaluationPoints(ic,k,i);
              }
            }

            pts::getL2BasisCoeffs(basisCoeffsL2,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
          {
            ValueType diffErr =0;
            for(int k=0;k<basisCardinality;k++) {
              for(int ic=0; ic<numCells; ic++)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsL2(ic,k)));
            }

            if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HGRAD_C" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
            }
          }
        }
      }
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

  } catch (std::exception err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }




  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 2 (Orientation - HCURL)                                |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);

      for (ordinal_type shift=0; shift<1; ++shift) {
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceLeft[ii]] = hexas_orig[0][faceLeftOriented[ii]];
            hexas_rotated[0][faceRight[ii]] = hexas_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceFront[ii]] = hexas_orig[0][faceFrontOriented[ii]];
            hexas_rotated[0][faceBack[ii]] = hexas_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            hexas[i][j] = reorder[hexas_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d){
            vertices[i][d] = vertices_orig[orderback[i]][d];
          }


        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hex.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hex.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];



        //computing edges and tangents
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][4];

        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<hex.getSideCount(); ++is) {
              for (std::size_t k=0; k<hex.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hex.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<hex.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<hex.getNodeCount(1,ie); ++k)
                edge[k]= hexas_rotated[i][hex.getNodeMap(1,ie,k)];
              std::sort(edge.begin(),edge.end());
              auto it=common_edges.find(edge);
              if(it !=common_edges.end()){
                auto edge_lid = std::distance(common_edges.begin(),it);
                edgeIndexes[i][edge_lid]=ie;
              }
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hex);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {

          Teuchos::RCP<Basis<DeviceSpaceType,ValueType,ValueType> > basisPtr;
          if(degree==1)
            basisPtr = Teuchos::rcp(new Basis_HCURL_HEX_I1_FEM<DeviceSpaceType,ValueType,ValueType>());
          else
            basisPtr = Teuchos::rcp(new Basis_HCURL_HEX_In_FEM<DeviceSpaceType,ValueType,ValueType>(degree));

          ordinal_type basisCardinality = basisPtr->getCardinality();

          //compute DofCoords Oriented
          DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(dofCoeffs, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

          //compute Lagrangian Interpolation of fun
          {
            li::getDofCoordsAndCoeffs(dofCoordsOriented, dofCoeffs, basisPtr.get(),POINTTYPE_EQUISPACED, elemOrts);

            //Compute physical Dof Coordinates

            {
              Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
              DynRankView ConstructWithLabel(hexLinearBasisValuesAtDofCoords, numCells, hex.getNodeCount(), basisCardinality);
              for(ordinal_type i=0; i<numCells; ++i)
                for(ordinal_type d=0; d<dim; ++d) {
                  auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                  auto outView =Kokkos::subview( hexLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
                  hexLinearBasis.getValues(outView, inView);

                  for(ordinal_type j=0; j<basisCardinality; ++j)
                    for(std::size_t k=0; k<hex.getNodeCount(); ++k)
                      physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexLinearBasisValuesAtDofCoords(i,k,j);
                }
            }

            DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
            ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hex);

            FunCurl fun;
            DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality, dim);
            for(ordinal_type i=0; i<numCells; ++i) {
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                for(ordinal_type k=0; k<dim; ++k)
                  funAtDofCoords(i,j,k) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k);
                for(ordinal_type k=0; k<dim; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    fwdFunAtDofCoords(i,j,k) += jacobian(i,j,d,k)*funAtDofCoords(i,j,d);
              }
            }

            li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, dofCoeffs);
          }


          //Testing Kronecker property of basis functions
          {
            for(ordinal_type i=0; i<numCells; ++i) {
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              basisPtr->getValues(outView, inView);

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr.get());

              for(ordinal_type k=0; k<basisCardinality; ++k) {
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  ValueType dofValue=0;
                  for(ordinal_type d=0; d<dim; ++d)
                    dofValue += basisValuesAtDofCoordsOriented(i,k,j,d) * dofCoeffs(i,j,d);
                  if ( k==j && std::abs( dofValue - 1.0 ) > 100*tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not have unit value at its node (" << dofValue <<")\n";
                  }
                  if ( k!=j && std::abs( dofValue ) > tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
                  }
                }
              }
            }
          }

          //check that fun values are consistent on common edges dofs
          {
            bool areDifferent(false);
            auto numEdgeDOFs = basisPtr->getDofCount(1,0);
            for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
              for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
                areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                    - basisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
              }
              if(areDifferent) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                *outStream << "Function DOFs for Hex 0 are:";
                for(ordinal_type j=0;j<numEdgeDOFs;j++)
                  *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j));
                *outStream << "\nFunction DOFs for Hex 1 are:";
                for(ordinal_type j=0;j<numEdgeDOFs;j++)
                  *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j));
                *outStream << std::endl;
              }
            }
          }

          //check that fun values are consistent on common face dofs
          {
            bool areDifferent(false);
            auto numFaceDOFs = basisPtr->getDofCount(2,0);
            for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
              areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                  - basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
            }

            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
              *outStream << "Function DOFs for Hex 0 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
              *outStream << "\nFunction DOFs for Hex 1 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
              *outStream << std::endl;
            }
          }

          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          DynRankView basisValuesAtDofCoordsCells("inValues", numCells, basisCardinality, basisCardinality, dim);

          for (ordinal_type ic = 0; ic < numCells; ++ic)
            basisPtr->getValues(Kokkos::subview(basisValuesAtDofCoordsCells, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(dofCoordsOriented, ic, Kokkos::ALL(), Kokkos::ALL()));

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoordsCells,
              elemOrts,
              basisPtr.get());

          // transform basis values
          DynRankView ConstructWithLabel(jacobianAtDofCoords, numCells, basisCardinality, dim, dim);
          DynRankView ConstructWithLabel(jacobianAtDofCoords_inv, numCells, basisCardinality, dim, dim);
          ct::setJacobian(jacobianAtDofCoords, dofCoordsOriented, physVertexes, hex);
          ct::setJacobianInv (jacobianAtDofCoords_inv, jacobianAtDofCoords);
          fst::HCURLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
              jacobianAtDofCoords_inv,
              basisValuesAtDofCoordsOriented);
          DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality, dim);
          for(ordinal_type i=0; i<numCells; ++i) {
            ValueType error=0;
            for(ordinal_type j=0; j<basisCardinality; ++j)
              for(ordinal_type d=0; d<dim; ++d) {
                for(ordinal_type k=0; k<basisCardinality; ++k)
                  funAtDofCoordsOriented(i,j,d) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j,d);

                error = std::max(std::abs( funAtDofCoords(i,j,d) - funAtDofCoordsOriented(i,j,d)), error);
              }

            if(error>100*tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
              *outStream << "Function values at reference points are:\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoords(i,j,0) << "," << funAtDofCoords(i,j,1) << ", " << funAtDofCoords(i,j,2) << ")";
              *outStream << "\nFunction values at reference points computed using basis functions are\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoordsOriented(i,j,0) << "," << funAtDofCoordsOriented(i,j,1) << ", " << funAtDofCoordsOriented(i,j,2) << ")";
              *outStream << std::endl;
            }
          }

          //compute projection-based interpolation of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree()-1);

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createHCurlProjectionStruct(basisPtr.get(), targetCubDegree, targetDerivCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numCurlPoints = projStruct.getNumTargetDerivEvalPoints();
            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
            DynRankView ConstructWithLabel(evaluationCurlPoints, numCells, numCurlPoints, dim);

            pts::getHCurlEvaluationPoints(evaluationPoints,
                evaluationCurlPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);


            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
            DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numCurlPoints, dim);

            DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
            DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints, dim);
            for(int ic=0; ic<numCells; ic++)
              basisPtr->getValues(Kokkos::subview(hcurlBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hcurlBasisAtEvaluationPoints,
                hcurlBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());

            DynRankView ConstructWithLabel(curlOfHCurlBasisAtEvaluationPoints, numCells, basisCardinality , numCurlPoints, dim);
            if(numCurlPoints>0) {
              DynRankView ConstructWithLabel(curlOfHCurlBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numCurlPoints, dim);
              for(int ic=0; ic<numCells; ic++)
                basisPtr->getValues(Kokkos::subview(curlOfHCurlBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationCurlPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_CURL);
              ots::modifyBasisByOrientation(curlOfHCurlBasisAtEvaluationPoints,
                  curlOfHCurlBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr.get());
            }


            for(int ic=0; ic<numCells; ic++){
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hcurlBasisAtEvaluationPoints(ic,k,i,d);
              }
              for(int i=0;i<numCurlPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetCurlAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*curlOfHCurlBasisAtEvaluationPoints(ic,k,i,d);//funHCurlCoeffs(k)
              }
            }

            pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
                targetAtEvalPoints,
                targetCurlAtEvalPoints,
                evaluationPoints,
                evaluationCurlPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
          {
            ValueType diffErr(0);

            for(int k=0;k<basisCardinality;k++) {
              for(int ic=0; ic<numCells; ic++)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsHCurl(ic,k)));
            }


            if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HCURL_I" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
            }
          }

          //compute L2 projection of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree());

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createL2ProjectionStruct(basisPtr.get(), targetCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints();
            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);

            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);


            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

            DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
            DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints, dim);
            for(int ic=0; ic<numCells; ic++)
              basisPtr->getValues(Kokkos::subview(hcurlBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hcurlBasisAtEvaluationPoints,
                hcurlBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());

            for(int ic=0; ic<numCells; ic++){
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hcurlBasisAtEvaluationPoints(ic,k,i,d);
              }
            }

            pts::getL2BasisCoeffs(basisCoeffsL2,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
          {
            ValueType diffErr = 0;
            for(int k=0;k<basisCardinality;k++) {
              for(int ic=0; ic<numCells; ic++)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsL2(ic,k)));
            }

            if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HCURL_I" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the functon."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
            }
          }
        }
      }
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

  } catch (std::exception err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }



  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 3 (Orientation - HDIV)                                 |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);

      for (ordinal_type shift=0; shift<1; ++shift) {
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceLeft[ii]] = hexas_orig[0][faceLeftOriented[ii]];
            hexas_rotated[0][faceRight[ii]] = hexas_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceFront[ii]] = hexas_orig[0][faceFrontOriented[ii]];
            hexas_rotated[0][faceBack[ii]] = hexas_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            hexas[i][j] = reorder[hexas_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hex.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hex.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];


        //computing edges and tangents

        ordinal_type faceIndex[numCells];
        {
          faceType face={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<hex.getSideCount(); ++is) {
              for (std::size_t k=0; k<hex.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hex.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hex);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {

          Teuchos::RCP<Basis<DeviceSpaceType,ValueType,ValueType> > basisPtr;
          if(degree==1)
            basisPtr = Teuchos::rcp(new Basis_HDIV_HEX_I1_FEM<DeviceSpaceType,ValueType,ValueType>());
          else
            basisPtr = Teuchos::rcp(new Basis_HDIV_HEX_In_FEM<DeviceSpaceType,ValueType,ValueType>(degree));

          ordinal_type basisCardinality = basisPtr->getCardinality();

          //compute DofCoords Oriented

          DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(dofCoeffs, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

          //compute Lagrangian Interpolation of fun
          {

            li::getDofCoordsAndCoeffs(dofCoordsOriented,  dofCoeffs, basisPtr.get(), POINTTYPE_WARPBLEND, elemOrts);

            //Compute physical Dof Coordinates
            Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
            DynRankView ConstructWithLabel(hexLinearBasisValuesAtDofCoords, numCells, hex.getNodeCount(), basisCardinality);
            for(ordinal_type i=0; i<numCells; ++i)
              for(ordinal_type d=0; d<dim; ++d) {
                auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                auto outView =Kokkos::subview( hexLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
                hexLinearBasis.getValues(outView, inView);

                for(ordinal_type j=0; j<basisCardinality; ++j)
                  for(std::size_t k=0; k<hex.getNodeCount(); ++k)
                    physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexLinearBasisValuesAtDofCoords(i,k,j);
              }

            //need to transform dofCoeff to physical space (they transform as normals)
            DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
            DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
            DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
            ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hex);
            ct::setJacobianInv (jacobian_inv, jacobian);
            ct::setJacobianDet (jacobian_det, jacobian);


            FunDiv fun;
            DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality, dim);
            for(ordinal_type i=0; i<numCells; ++i) {
              for(ordinal_type j=0; j<basisCardinality; ++j){
                for(ordinal_type k=0; k<dim; ++k)
                  funAtDofCoords(i,j,k) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k);
                for(ordinal_type k=0; k<dim; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    fwdFunAtDofCoords(i,j,k) += jacobian_det(i,j)*jacobian_inv(i,j,k,d)*funAtDofCoords(i,j,d);

              }
            }
            li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, dofCoeffs);
          }

          //Testing Kronecker property of basis functions
          {
            for(ordinal_type i=0; i<numCells; ++i) {
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              basisPtr->getValues(outView, inView);

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr.get());

              DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
              DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
              ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hex);
              ct::setJacobianDet (jacobian_det, jacobian);

              for(ordinal_type k=0; k<basisCardinality; ++k) {
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  ValueType dofValue=0;
                  for(ordinal_type d=0; d<dim; ++d)
                    dofValue += basisValuesAtDofCoordsOriented(i,k,j,d) * dofCoeffs(i,j,d);
                  if ( k==j && std::abs( dofValue - 1.0 ) > 100*tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not have unit value at its node (" << dofValue <<")\n";
                  }
                  if ( k!=j && std::abs( dofValue ) > tol ) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    *outStream << " Basis function " << k << " of cell " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
                  }
                }
              }
            }
          }

          //check that fun values are consistent on common face dofs
          {
            bool areDifferent(false);
            auto numFaceDOFs = basisPtr->getDofCount(2,0);
            for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
              areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                  - basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
            }

            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
              *outStream << "Function DOFs for Hex 0 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
              *outStream << "\nFunction DOFs for Hex 1 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
              *outStream << std::endl;
            }
          }

          //check that fun values are consistent on common face dofs
          {
            bool areDifferent(false);
            auto numFaceDOFs = basisPtr->getDofCount(2,0);
            for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
              areDifferent = std::abs(basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                  - basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
            }

            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
              *outStream << "Function DOFs for Hex 0 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
              *outStream << "\nFunction DOFs for Hex 1 are:";
              for(ordinal_type j=0;j<numFaceDOFs;j++)
                *outStream << " " << basisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
              *outStream << std::endl;
            }
          }

          //check that fun values at reference points coincide with those computed using basis functions
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          DynRankView basisValuesAtDofCoordsCells("inValues", numCells, basisCardinality, basisCardinality, dim);

          for (ordinal_type ic = 0; ic < numCells; ++ic)
            basisPtr->getValues(Kokkos::subview(basisValuesAtDofCoordsCells, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(dofCoordsOriented, ic, Kokkos::ALL(), Kokkos::ALL()));

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoordsCells,
              elemOrts,
              basisPtr.get());

          // transform basis values
          DynRankView ConstructWithLabel(jacobianAtDofCoords, numCells, basisCardinality, dim, dim);
          DynRankView ConstructWithLabel(jacobianAtDofCoords_det, numCells, basisCardinality);
          ct::setJacobian(jacobianAtDofCoords, dofCoordsOriented, physVertexes, hex);
          ct::setJacobianDet (jacobianAtDofCoords_det, jacobianAtDofCoords);
          fst::HDIVtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
              jacobianAtDofCoords,
              jacobianAtDofCoords_det,
              basisValuesAtDofCoordsOriented);
          DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality, dim);
          for(ordinal_type i=0; i<numCells; ++i) {
            ValueType error=0;
            for(ordinal_type j=0; j<basisCardinality; ++j)
              for(ordinal_type d=0; d<dim; ++d) {
                for(ordinal_type k=0; k<basisCardinality; ++k)
                  funAtDofCoordsOriented(i,j,d) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j,d);

                error = std::max(std::abs( funAtDofCoords(i,j,d) - funAtDofCoordsOriented(i,j,d)), error);
              }

            if(error>100*tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
              *outStream << "Function values at reference points are:\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoords(i,j,0) << "," << funAtDofCoords(i,j,1) << ", " << funAtDofCoords(i,j,2) << ")";
              *outStream << "\nFunction values at reference points computed using basis functions are\n";
              for(ordinal_type j=0; j<basisCardinality; ++j)
                *outStream << " (" << funAtDofCoordsOriented(i,j,0) << "," << funAtDofCoordsOriented(i,j,1) << ", " << funAtDofCoordsOriented(i,j,2) << ")";
              *outStream << std::endl;
            }
          }

          //compute projection-based interpolation of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree()-1);

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createHDivProjectionStruct(basisPtr.get(), targetCubDegree, targetDerivCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints(), numDivPoints = projStruct.getNumTargetDerivEvalPoints();

            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);
            DynRankView ConstructWithLabel(evaluationDivPoints, numCells, numDivPoints, dim);

            pts::getHDivEvaluationPoints(evaluationPoints,
                evaluationDivPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);

            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
            DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);

            DynRankView ConstructWithLabel(hdivBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
            DynRankView ConstructWithLabel(hdivBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints, dim);
            for(ordinal_type ic=0; ic<numCells; ++ic)
              basisPtr->getValues(Kokkos::subview(hdivBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hdivBasisAtEvaluationPoints,
                hdivBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());

            DynRankView ConstructWithLabel(divOfHDivBasisAtEvaluationPoints, numCells, basisCardinality , numDivPoints);
            if(numDivPoints>0) {
              DynRankView ConstructWithLabel(divOfHDivBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numDivPoints);
              for(ordinal_type ic=0; ic<numCells; ++ic)
                basisPtr->getValues(Kokkos::subview(divOfHDivBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationDivPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_DIV);
              ots::modifyBasisByOrientation(divOfHDivBasisAtEvaluationPoints,
                  divOfHDivBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr.get());
            }



            for(ordinal_type ic=0; ic<numCells; ++ic) {
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hdivBasisAtEvaluationPoints(ic,k,i,d);
              }
              for(int i=0;i<numDivPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  targetDivAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*divOfHDivBasisAtEvaluationPoints(ic,k,i);//basisCoeffsLI(k)
              }
            }

            pts::getHDivBasisCoeffs(basisCoeffsHDiv,
                targetAtEvalPoints,
                targetDivAtEvalPoints,
                evaluationPoints,
                evaluationDivPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
          {
            ValueType diffErr(0);
            for(int k=0;k<basisCardinality;k++) {
              //std::cout << "["<< basisCoeffsLI(0,k) << " " <<  basisCoeffsHDiv(0,k) << "] [" << basisCoeffsLI(1,k) << " " <<  basisCoeffsHDiv(1,k) << "]" <<std::endl;
              for(ordinal_type ic=0; ic<numCells; ++ic)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsHDiv(ic,k)));
            }

            if(diffErr > pow(8, degree)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HDIV_I" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << " (tol: " << pow(8, degree)*tol << ")" << std::endl;
            }
          }


          //compute L2 projection interpolation of the Lagrangian interpolation
          DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
          {
            ordinal_type targetCubDegree(basisPtr->getDegree());

            Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
            projStruct.createL2ProjectionStruct(basisPtr.get(), targetCubDegree);

            ordinal_type numPoints = projStruct.getNumTargetEvalPoints();

            DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);

            pts::getL2EvaluationPoints(evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);

            DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

            DynRankView ConstructWithLabel(hdivBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
            DynRankView ConstructWithLabel(hdivBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints, dim);
            for(ordinal_type ic=0; ic<numCells; ++ic)
              basisPtr->getValues(Kokkos::subview(hdivBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
            ots::modifyBasisByOrientation(hdivBasisAtEvaluationPoints,
                hdivBasisAtEvaluationPointsNonOriented,
                elemOrts,
                basisPtr.get());

            for(ordinal_type ic=0; ic<numCells; ++ic) {
              for(int i=0;i<numPoints;i++) {
                for(int k=0;k<basisCardinality;k++)
                  for(int d=0;d<dim;d++)
                    targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hdivBasisAtEvaluationPoints(ic,k,i,d);
              }
            }

            pts::getL2BasisCoeffs(basisCoeffsL2,
                targetAtEvalPoints,
                evaluationPoints,
                elemOrts,
                basisPtr.get(),
                &projStruct);
          }

          //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
          {
            ValueType diffErr = 0;
            for(int k=0;k<basisCardinality;k++) {
              //std::cout << "["<< basisCoeffsLI(0,k) << " " <<  basisCoeffsHDiv(0,k) << "] [" << basisCoeffsLI(1,k) << " " <<  basisCoeffsHDiv(1,k) << "]" <<std::endl;
              for(ordinal_type ic=0; ic<numCells; ++ic)
                diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsL2(ic,k)));
            }

            if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "HDIV_I" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the function."<<
                  "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
            }
          }
        }
      }
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

  } catch (std::exception err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }


  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 4 (Orientation - HVOL)                                 |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {


    ValueType vertices[numTotalVertexes][dim];
    ordinal_type hexas[numCells][numElemVertexes];

    for(ordinal_type i=0; i<numCells;++i)
      for(ordinal_type j=0; j<numElemVertexes;++j)
        hexas[i][j] = hexas_orig[i][j];

    for(ordinal_type i=0; i<numTotalVertexes;++i)
      for(ordinal_type d=0; d<dim;++d)
        vertices[i][d] = vertices_orig[i][d];

    *outStream <<  "Considering Hex 0: [ ";
    for(ordinal_type j=0; j<numElemVertexes;++j)
      *outStream << hexas[0][j] << " ";
    *outStream << "] and Hex 1: [ ";
    for(ordinal_type j=0; j<numElemVertexes;++j)
      *outStream << hexas[1][j] << " ";
    *outStream << "]\n";

    shards::CellTopology hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
    shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
    shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

    //computing vertices coords
    DynRankView ConstructWithLabel(physVertexes, numCells, hex.getNodeCount(), dim);
    for(ordinal_type i=0; i<numCells; ++i)
      for(std::size_t j=0; j<hex.getNodeCount(); ++j)
        for(ordinal_type k=0; k<dim; ++k)
          physVertexes(i,j,k) = vertices[hexas[i][j]][k];

    // compute orientations for cells (one time computation)
    DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
    Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
    ots::getOrientation(elemOrts, elemNodes, hex);

    for (ordinal_type degree=1; degree <= max_degree; degree++) {

      Teuchos::RCP<Basis<DeviceSpaceType,ValueType,ValueType> > basisPtr;
      if(degree==1)
        basisPtr = Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceSpaceType,ValueType,ValueType>(hex));
      else
        basisPtr = Teuchos::rcp(new Basis_HVOL_HEX_Cn_FEM<DeviceSpaceType,ValueType,ValueType>(degree));
      ordinal_type basisCardinality = basisPtr->getCardinality();

      //compute DofCoords Oriented
      DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
      DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality);
      DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
      DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality);
      DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

      //compute Lagrangian Interpolation of fun
      {
        li::getDofCoordsAndCoeffs(dofCoordsOriented,  dofCoeffsPhys, basisPtr.get(), POINTTYPE_EQUISPACED, elemOrts);

        //Compute physical Dof Coordinates
        {
          Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,ValueType,ValueType> hexLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexLinearBasisValuesAtDofCoords, numCells, hex.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( hexLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              hexLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<hex.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexLinearBasisValuesAtDofCoords(i,k,j);
            }
        }

        //need to transform dofCoeff to physical space (they transform as normals)
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hex);
        ct::setJacobianDet (jacobian_det, jacobian);

        Fun fun;

        DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<basisCardinality; ++j) {
            funAtDofCoords(i,j) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2));
            fwdFunAtDofCoords(i,j) = jacobian_det(i,j)*funAtDofCoords(i,j);
          }

        li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, dofCoeffsPhys);
      }

      //Testing Kronecker property of basis functions
      {
        for(ordinal_type i=0; i<numCells; ++i) {
          DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality);
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
          auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
          auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL());
          basisPtr->getValues(outView, inView);

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoords,
              elemOrts,
              basisPtr.get());

          for(ordinal_type k=0; k<basisCardinality; ++k) {
            for(ordinal_type j=0; j<basisCardinality; ++j){
              ValueType dofValue = basisValuesAtDofCoordsOriented(i,k,j) * dofCoeffsPhys(i,j);
              if ( k==j && std::abs( dofValue - 1.0 ) > 100*tol ) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Basis function " << k << " of cell " << i << " does not have unit value at its node (" << dofValue <<")\n";
              }
              if ( k!=j && std::abs( dofValue ) > tol ) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Basis function " << k << " of cell " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
              }
            }
          }
        }
      }

      //check that fun values at reference points coincide with those computed using basis functions
      DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
      DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
      DynRankView basisValuesAtDofCoordsCells("inValues", numCells, basisCardinality, basisCardinality);

      for (ordinal_type ic = 0; ic < numCells; ++ic)
        basisPtr->getValues(Kokkos::subview(basisValuesAtDofCoordsCells, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(dofCoordsOriented, ic, Kokkos::ALL(), Kokkos::ALL()));

      // modify basis values to account for orientations
      ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
          basisValuesAtDofCoordsCells,
          elemOrts,
          basisPtr.get());

      // transform basis values
      DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
      DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
      ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hex);
      ct::setJacobianDet (jacobian_det, jacobian);
      fst::HVOLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
          jacobian_det,
          basisValuesAtDofCoordsOriented);

      DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality);
      for(ordinal_type i=0; i<numCells; ++i) {
        ValueType error=0;
        for(ordinal_type j=0; j<basisCardinality; ++j) {
          for(ordinal_type k=0; k<basisCardinality; ++k)
            funAtDofCoordsOriented(i,j) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j);

          error = std::max(std::abs( funAtDofCoords(i,j) - funAtDofCoordsOriented(i,j)), error);
        }

        if(error>100*tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
          *outStream << "Function values at reference points are:\n";
          for(ordinal_type j=0; j<basisCardinality; ++j)
            *outStream << " (" << funAtDofCoords(i,j)  << ")";
          *outStream << "\nFunction values at reference points computed using basis functions are\n";
          for(ordinal_type j=0; j<basisCardinality; ++j)
            *outStream << " (" << funAtDofCoordsOriented(i,j)  << ")";
          *outStream << std::endl;
        }
      }

      //compute projection-based interpolation of the Lagrangian interpolation
      DynRankView ConstructWithLabel(basisCoeffsHVol, numCells, basisCardinality);
      {
        ordinal_type targetCubDegree(basisPtr->getDegree());

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createHVolProjectionStruct(basisPtr.get(), targetCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints();
        DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);


        pts::getHVolEvaluationPoints(evaluationPoints,
            elemOrts,
            basisPtr.get(),
            &projStruct);


        DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

        DynRankView ConstructWithLabel(hvolBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
        DynRankView ConstructWithLabel(hvolBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints);
        for(int ic=0; ic<numCells; ic++)
          basisPtr->getValues(Kokkos::subview(hvolBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
        ots::modifyBasisByOrientation(hvolBasisAtEvaluationPoints,
            hvolBasisAtEvaluationPointsNonOriented,
            elemOrts,
            basisPtr.get());

        for(int ic=0; ic<numCells; ic++) {
          for(int i=0;i<numPoints;i++) {
            for(int k=0;k<basisCardinality;k++)
              targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hvolBasisAtEvaluationPoints(ic,k,i);
          }
        }

        pts::getHVolBasisCoeffs(basisCoeffsHVol,
            targetAtEvalPoints,
            evaluationPoints,
            elemOrts,
            basisPtr.get(),
            &projStruct);
      }

      //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
      {
        ValueType diffErr(0);
        for(int k=0;k<basisCardinality;k++) {
          for(int ic=0; ic<numCells; ic++)
            diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsHVol(ic,k)));
        }

        //Check that the two representations of the gradient of ifun are consistent
        if(diffErr > pow(16, degree)*tol) { //heuristic relation on how round-off error depends on degree
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "HGRAD_C" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
              "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
        }
      }

      //compute L2 projection of the Lagrangian interpolation
      DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
      {
        ordinal_type targetCubDegree(basisPtr->getDegree());

        Experimental::ProjectionStruct<DeviceSpaceType,ValueType> projStruct;
        projStruct.createL2ProjectionStruct(basisPtr.get(), targetCubDegree);

        ordinal_type numPoints = projStruct.getNumTargetEvalPoints();
        DynRankView ConstructWithLabel(evaluationPoints, numCells, numPoints, dim);


        pts::getL2EvaluationPoints(evaluationPoints,
            elemOrts,
            basisPtr.get(),
            &projStruct);


        DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

        DynRankView ConstructWithLabel(hvolBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
        DynRankView ConstructWithLabel(hvolBasisAtEvaluationPointsNonOriented, numCells, basisCardinality , numPoints);
        for(int ic=0; ic<numCells; ic++)
          basisPtr->getValues(Kokkos::subview(hvolBasisAtEvaluationPointsNonOriented, ic, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()), OPERATOR_VALUE);
        ots::modifyBasisByOrientation(hvolBasisAtEvaluationPoints,
            hvolBasisAtEvaluationPointsNonOriented,
            elemOrts,
            basisPtr.get());

        for(int ic=0; ic<numCells; ic++) {
          for(int i=0;i<numPoints;i++) {
            for(int k=0;k<basisCardinality;k++)
              targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hvolBasisAtEvaluationPoints(ic,k,i);
          }
        }

        pts::getL2BasisCoeffs(basisCoeffsL2,
            targetAtEvalPoints,
            evaluationPoints,
            elemOrts,
            basisPtr.get(),
            &projStruct);
      }

      //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
      {
        ValueType diffErr = 0;
        for(int k=0;k<basisCardinality;k++) {
          for(int ic=0; ic<numCells; ic++)
            diffErr = std::max(diffErr, std::abs(basisCoeffsLI(ic,k) - basisCoeffsL2(ic,k)));
        }

        if(diffErr > pow(16, degree)*tol) { //heuristic relation on how round-off error depends on degree
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "HGRAD_C" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the functon."<<
              "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
        }
      }
    }
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

