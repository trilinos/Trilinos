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
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
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
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceSpaceType>
int OrientationTet(const bool verbose) {

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

  struct Fun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return (x+1)*(y-2)*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      switch (comp) {
      case 0:
        return 5+2*x;
        break;
      case 1:
        return -7+2*y;
        break;
      case 2:
        return 0.5+2*z;
        break;
      default:
        return 0;
      }
    }
  };

  typedef std::array<ordinal_type,2> edgeType;
  typedef std::array<ordinal_type,3> faceType;
  typedef CellTools<DeviceSpaceType> ct;
  typedef OrientationTools<DeviceSpaceType> ots;
  typedef RealSpaceTools<DeviceSpaceType> rst;
  typedef FunctionSpaceTools<DeviceSpaceType> fst;
  typedef ArrayTools<DeviceSpaceType> at;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numTetVertexes = 4;

  ValueType  vertices_orig[5][dim] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
  ordinal_type tets_orig[numCells][numTetVertexes] = {{0,1,2,3},{1,2,3,4}};  //common face is {1,2,3}
  faceType common_face = {1,2,3};
  std::set<edgeType> common_edges;
  common_edges.insert(edgeType({1,2})); common_edges.insert(edgeType({1,3})); common_edges.insert(edgeType({2,3}));


  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 1 (Orientation - HGRAD)                                |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    const ordinal_type order = 4;
    ordinal_type reorder[5] = {0,1,2,3,4};

    do {
      ordinal_type orderback[5];
      for(ordinal_type i=0;i<5;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[5][dim];
      ordinal_type tets[numCells][numTetVertexes];
      for (ordinal_type shift=0; shift<4; ++shift) {
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numTetVertexes;++j)
            tets[i][j] = reorder[tets_orig[i][j]];

        for(ordinal_type i=0; i<5;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tet 0: [ ";
        for(ordinal_type j=0; j<numTetVertexes;++j)
          *outStream << tets[0][j] << " ";
        *outStream << "] and Tet 1: [ ";
        for(ordinal_type j=0; j<numTetVertexes;++j)
          *outStream << tets[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND);
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        Basis_HGRAD_LINE_C1_FEM<DeviceSpaceType,ValueType,ValueType> lineLinearBasis;
        Basis_HGRAD_TRI_C1_FEM<DeviceSpaceType,ValueType,ValueType> faceLinearBasis;


        faceType face={};
        edgeType edge={};
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][3];
        for(ordinal_type i=0; i<numCells; ++i) {
          for (std::size_t is=0; is<tet.getSideCount(); ++is) {
            for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
              face[k]= tets_orig[i][tet.getNodeMap(2,is,k)];
            std::sort(face.begin(),face.end());
            if(face == common_face)
              faceIndex[i]=is;
          }
          //ordinal_type edge_count=0;
          for (std::size_t ie=0; ie<tet.getEdgeCount(); ++ie) {
            for (std::size_t k=0; k<tet.getNodeCount(1,ie); ++k)
              edge[k]= tets_orig[i][tet.getNodeMap(1,ie,k)];
            std::sort(edge.begin(),edge.end());
            auto it=common_edges.find(edge);
            if(it !=common_edges.end())
              edgeIndexes[i][std::distance(common_edges.begin(),it)]=ie;
          }
        }

        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(physVertexes, numCells, tet.getNodeCount(), dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtDofCoords, tet.getNodeCount(), basisCardinality);
        DynRankView ConstructWithLabel(funAtPhysDofCoords, numCells, basisCardinality);
        DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords);
        DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords);
        DynRankView ConstructWithLabel(funAtPhysDofCoordsOriented, numCells, basisCardinality);


        const auto allTags = basis.getAllDofTags();

        basis.getDofCoords(dofCoords);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<tet.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[tets[i][j]][k];

        warpBasis.getDofCoords(refPoints);

        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        tetLinearBasis.getValues(tetLinearBasisValuesAtDofCoords, dofCoords);

        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type d=0; d<dim; ++d) {
            for(ordinal_type j=0; j<basisCardinality; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtRefCoords(k,j);
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physDofCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtDofCoords(k,j);
          }



//        auto edgeTopo = shards::getCellTopologyData<shards::Line<2> >();
//        auto faceTopo = shards::getCellTopologyData<shards::Triangle<3> >();
//        for(int i=0; i<numCells; ++i)
//          for(int j=0; j<basisCardinality; ++j) {
//            if(allTags(j,0) == 0) { //point
//              ordinal_type point_id = allTags(j,1);
//              for(int d=0; d<dim; ++d)
//                physDofCoords(i,j, d) = vertices[tets[i][point_id]][d];
//            }
//            if(allTags(j,0) == 1) { //edge
//              ordinal_type edge_id = allTags(j,1);
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numEdgeDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numEdgeDofs, 1);
//              PointTools::getLattice( refPoints,
//                  edgeTopo,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValuesAtDofCoords, lineLinearBasis.getCardinality(), numEdgeDofs);
//              lineLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//              for (ordinal_type k=0; k<tet.getNodeCount(1,edge_id); ++k)
//                edge[k]= tets[i][tet.getNodeMap(1,edge_id,k)];
//              std::sort(edge.begin(),edge.end());
//              //   std::cout << "[" << edge[0] << ", " << edge[1] <<"]";// <<std::endl;
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<lineLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[edge[k]][d]*basisValuesAtDofCoords(k,ldof_id);
//              // std::cout << " (" << physDofCoords(i,j,0) << ", " << physDofCoords(i,j,1) << ", " << physDofCoords(i,j,2) <<")" <<std::endl;
//            }
//            else if(allTags(j,0) == 2) { //face
//              ordinal_type face_id = allTags(j,1);
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numFaceDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numFaceDofs, 2);
//              PointTools::getLattice( refPoints,
//                  faceTopo,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValuesAtDofCoords, faceLinearBasis.getCardinality(), numFaceDofs);
//              faceLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//              for (ordinal_type k=0; k<tet.getNodeCount(2,face_id); ++k)
//                face[k]= tets[i][tet.getNodeMap(2,face_id,k)];
//              std::sort(face.begin(),face.end());
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<faceLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[face[k]][d]*basisValuesAtDofCoords(k,ldof_id);
//            }
//            else if (allTags(j,0) == 3) { //tet
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numElemDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numElemDofs, dim);
//              PointTools::getLattice( refPoints,
//                  tet,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValuesAtDofCoords, tetLinearBasis.getCardinality(), numElemDofs);
//              tetLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<faceLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[tets[i][k]][d]*basisValuesAtDofCoords(k,ldof_id);
//            }
//          }



        Fun fun;
        for(ordinal_type i=0; i<numCells; ++i) {
          for(ordinal_type j=0; j<numRefCoords; ++j)
            funAtPhysRefCoords(i,j) = fun(physRefCoords(i,j,0) ,physRefCoords(i,j,1) ,physRefCoords(i,j,2));
          for(ordinal_type j=0; j<basisCardinality; ++j)
            funAtPhysDofCoords(i,j) =fun(physDofCoords(i,j,0),physDofCoords(i,j,1),physDofCoords(i,j,2));
          //   std::cout << funAtDofCoords(i,j) << " ";
          //   std::cout << "(" << physDofCoords(i,j,0) << ", " << physDofCoords(i,j,1) << ", " << physDofCoords(i,j,2) << ")" << std::endl;;
          //     std::cout << std::endl;
        }

        DynRankViewInt elemNodes(&tets[0][0], 2, 4);

        // compute orientations for cells (one time computation)
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);

        // cell specific modified basis
        DynRankView basisValuesAtRefCoordsOriented("outValues", numCells, basisCardinality, numRefCoords);

        //Kokkos::deep_copy(funAtDofCoords,0);
        ots::modifyBasisByOrientation(funAtPhysDofCoordsOriented, funAtPhysDofCoords,elemOrts, &basis);

        DynRankView basisValuesAtRefCoordsCelles("inValues", numCells, basisCardinality, numRefCoords);
        rst::clone(basisValuesAtRefCoordsCelles,basisValuesAtRefCoords);

        // modify refValues accounting for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCelles,
            elemOrts,
            &basis);

        auto ordinalToTag = basis.getAllDofTags();
        auto tagToOrdinal = basis.getAllDofOrdinal();
        auto numEdgeDOFs = ordinalToTag(tagToOrdinal(1, 0, 0), 3);
        auto numFaceDOFs = ordinalToTag(tagToOrdinal(2, 0, 0), 3);

        //check that fun values at reference points
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k){
              //auto a = (std::abs(basisValuesAtRefCoordsCelles(i,j,k)) > 1e-10) ? basisValuesAtRefCoordsCelles(i,j,k) : 0.0;
              //auto b = (std::abs(basisValuesAtRefDofCoordsOriented(i,j,k)) > 1e-10) ? basisValuesAtRefDofCoordsOriented(i,j,k) : 0.0;
              //std::cout <<  "(" << a << ", " << b <<") ";
              funAtRefCoordsOriented(i,j) += funAtPhysDofCoordsOriented(i,k)*basisValuesAtRefCoordsOriented(i,k,j);
            }
            error = std::max(std::abs( funAtPhysRefCoords(i,j) - funAtRefCoordsOriented(i,j)), error);
          }
          if(error>100*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function values at reference points differ from those computed using basis functions of Tet " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " " << funAtPhysRefCoords(i,j);
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " " << funAtRefCoordsOriented(i,j);
            *outStream << std::endl;
          }
        }

        //check that fun values are consistent on common face dofs
        {
          bool areDifferent(false);
          for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funAtPhysDofCoordsOriented(0,basis.getDofOrdinal(2,faceIndex[0],j))
                - funAtPhysDofCoordsOriented(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
          }
          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
            *outStream << "Function DOFs for Tet 0 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funAtPhysDofCoordsOriented(0,basis.getDofOrdinal(2,faceIndex[0],j));
            *outStream << "\nFunction DOFs for Tet 1 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funAtPhysDofCoordsOriented(1,basis.getDofOrdinal(2,faceIndex[1],j));
            *outStream << std::endl;
          }
        }

        //check that fun values are consistent on common edges dofs
        {
          bool areDifferent(false);
          for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
            for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
              areDifferent = std::abs(funAtPhysDofCoordsOriented(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                  - funAtPhysDofCoordsOriented(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
            }
            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common edge " << iEdge << " computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
              *outStream << "Function DOFs for Tet 0 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funAtPhysDofCoordsOriented(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
              *outStream << "\nFunction DOFs for Tet 1 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funAtPhysDofCoordsOriented(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j));
              *outStream << std::endl;
            }
          }
        }
        std::rotate(&tets_orig[0][0], &tets_orig[0][0]+1, &tets_orig[0][0]+4);
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

    const ordinal_type order = 1;
    ordinal_type reorder[5] = {0,1,2,3,4};

    do {
      ordinal_type orderback[5];
      for(ordinal_type i=0;i<5;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[5][dim];
      ordinal_type tets[numCells][numTetVertexes];
      for (ordinal_type shift=0; shift<4; ++shift) {
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numTetVertexes;++j)
            tets[i][j] = reorder[tets_orig[i][j]];

        for(ordinal_type i=0; i<5;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Tet 0: [ ";
        for(ordinal_type j=0; j<numTetVertexes;++j)
          *outStream << tets[0][j] << " ";
        *outStream << "] and Tet 1: [ ";
        for(ordinal_type j=0; j<numTetVertexes;++j)
          *outStream << tets[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
        Basis_HDIV_TET_In_FEM<DeviceSpaceType,ValueType,ValueType> basis(order);
        Basis_HGRAD_TET_Cn_FEM<DeviceSpaceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND);
        Basis_HGRAD_TET_C1_FEM<DeviceSpaceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
        Basis_HGRAD_LINE_C1_FEM<DeviceSpaceType,ValueType,ValueType> lineLinearBasis;
        Basis_HGRAD_TRI_C1_FEM<DeviceSpaceType,ValueType,ValueType> faceLinearBasis;


        faceType face={};
        edgeType edge={};
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][3];
        //bool faceOrientation[numCells][4];
        for(ordinal_type i=0; i<numCells; ++i) {
          for (std::size_t is=0; is<tet.getSideCount(); ++is) {
            for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
              face[k]= tets_orig[i][tet.getNodeMap(2,is,k)];
            std::sort(face.begin(),face.end());
            if(face == common_face)
              faceIndex[i]=is;

//            for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
//              face[k]= tets[i][tet.getNodeMap(2,is,k)];
//            std::rotate(face.begin(), std::min_element(face.begin(),face.end()),face.end());
//            auto sorted_face=face;
//            std::sort(sorted_face.begin(),sorted_face.end());
//            faceOrientation[i][is] = (face == sorted_face);
          }
          //ordinal_type edge_count=0;
          for (std::size_t ie=0; ie<tet.getEdgeCount(); ++ie) {
            for (std::size_t k=0; k<tet.getNodeCount(1,ie); ++k)
              edge[k]= tets_orig[i][tet.getNodeMap(1,ie,k)];
            std::sort(edge.begin(),edge.end());
            auto it=common_edges.find(edge);
            if(it !=common_edges.end())
              edgeIndexes[i][std::distance(common_edges.begin(),it)]=ie;
          }
        }

        ordinal_type basisCardinality = basis.getCardinality();
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(physVertexes, numCells, tet.getNodeCount(), dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(tetLinearBasisValuesAtDofCoords, tet.getNodeCount(), basisCardinality);
        DynRankView ConstructWithLabel(funDofs, numCells, basisCardinality);
        DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(funDofsOriented, numCells, basisCardinality);


        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_invT, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);

        DynRankView ConstructWithLabel(jacobianAtRefCoord, numCells, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoord_det, numCells, numRefCoords);




        const auto allTags = basis.getAllDofTags();

        basis.getDofCoords(dofCoords);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<tet.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[tets[i][j]][k];

        ct::setJacobian(jacobian, dofCoords, physVertexes, tet);
        ct::setJacobianDet (jacobian_det, jacobian);
        ct::setJacobianInv (jacobian_inv, jacobian);
        rst::transpose(jacobian_invT,jacobian_inv);

        warpBasis.getDofCoords(refPoints);
        ct::setJacobian(jacobianAtRefCoord, refPoints, physVertexes, tet);
        ct::setJacobianDet (jacobianAtRefCoord_det, jacobianAtRefCoord);

        tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        tetLinearBasis.getValues(tetLinearBasisValuesAtDofCoords, dofCoords);

        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type d=0; d<dim; ++d) {
            for(ordinal_type j=0; j<numRefCoords; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physRefCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtRefCoords(k,j);
            for(ordinal_type j=0; j<basisCardinality; ++j)
              for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                physDofCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtDofCoords(k,j);
          }


//
//        auto edgeTopo = shards::getCellTopologyData<shards::Line<2> >();
//        auto faceTopo = shards::getCellTopologyData<shards::Triangle<3> >();
//        for(int i=0; i<numCells; ++i)
//          for(int j=0; j<basisCardinality; ++j) {
//            if(allTags(j,0) == 0) { //point
//              ordinal_type point_id = allTags(j,1);
//              for(int d=0; d<dim; ++d)
//                physDofCoords(i,j, d) = vertices[tets[i][point_id]][d];
//            }
//            if(allTags(j,0) == 1) { //edge
//              ordinal_type edge_id = allTags(j,1);
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numEdgeDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numEdgeDofs, 1);
//              PointTools::getLattice( refPoints,
//                  edgeTopo,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValutetLinearBasisValuesAtDofCoordsearBasis.getCardinality(), numEdgeDofs);
//              lineLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//              for (ordinal_type k=0; k<tet.getNodeCount(1,edge_id); ++k)
//                edge[k]= tets[i][tet.getNodeMap(1,edge_id,k)];
//              std::sort(edge.begin(),edge.end());
//              //   std::cout << "[" << edge[0] << ", " << edge[1] <<"]";// <<std::endl;
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<lineLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[edge[k]][d]*basisValuesAtDofCoords(k,ldof_id);
//              // std::cout << " (" << physDofCoords(i,j,0) << ", " << physDofCoords(i,j,1) << ", " << physDofCoords(i,j,2) <<")" <<std::endl;
//            }
//            else if(allTags(j,0) == 2) { //face
//              ordinal_type face_id = allTags(j,1);
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numFaceDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numFaceDofs, 2);
//              PointTools::getLattice( refPoints,
//                  faceTopo,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValuesAtDofCoords, faceLinearBasis.getCardinality(), numFaceDofs);
//              faceLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//              for (ordinal_type k=0; k<tet.getNodeCount(2,face_id); ++k)
//                face[k]= tets[i][tet.getNodeMap(2,face_id,k)];
//              std::sort(face.begin(),face.end());
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<faceLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[face[k]][d]*basisValuesAtDofCoords(k,ldof_id);
//            }
//            else if (allTags(j,0) == 3) { //tet
//              ordinal_type ldof_id = allTags(j,2);
//              ordinal_type numElemDofs = allTags(j,3);
//              DynRankView ConstructWithLabel(refPoints, numElemDofs, dim);
//              PointTools::getLattice( refPoints,
//                  tet,
//                  order, 1);
//              DynRankView ConstructWithLabel(basisValuesAtDofCoords, tetLinearBasis.getCardinality(), numElemDofs);
//              tetLinearBasis.getValues(basisValuesAtDofCoords, refPoints, OPERATOR_VALUE);
//
//              for(int d=0; d<dim; ++d)
//                for(int k=0; k<faceLinearBasis.getCardinality(); ++k)
//                  physDofCoords(i,j,d) += vertices[tets[i][k]][d]*basisValuesAtDofCoords(k,ldof_id);
//            }
//          }



        FunDiv fun;
        DynRankView ConstructWithLabel(dofCoeffs, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);

        //need to tranform dofCoeff to physical space (they transform as normals)
        basis.getDofCoeffs(dofCoeffs);
        rst::clone(dofCoeffsPhys,dofCoeffs);
        rst::matvec(dofCoeffsTmp, jacobian_invT, dofCoeffsPhys);
        at::scalarMultiplyDataData(dofCoeffsPhys, jacobian_det, dofCoeffsTmp);

        for(ordinal_type i=0; i<numCells; ++i) {
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
          for(ordinal_type j=0; j<basisCardinality; ++j)
            for(ordinal_type k=0; k<dim; ++k) {
              funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k) * dofCoeffsPhys(i,j,k);
            }
        }

        DynRankViewInt elemNodes(&tets[0][0], 2, 4);

        // compute orientations for cells (one time computation)
        Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);

        // cell specific modified basis
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);

        ots::modifyBasisByOrientation(funDofsOriented, funDofs,elemOrts, &basis);

        DynRankView basisValuesAtRefCoordsCelles("inValues", numCells, basisCardinality, numRefCoords, dim);
        rst::clone(basisValuesAtRefCoordsCelles,basisValuesAtRefCoords);

        // modify refValues accounting for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCelles,
            elemOrts,
            &basis);

        fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
                            jacobianAtRefCoord,
                            jacobianAtRefCoord_det,
                            basisValuesAtRefCoordsOriented);

        auto numEdgeDOFs = basis.getDofCount(1,0);
        auto numFaceDOFs = basis.getDofCount(2,0);

        //check that fun values at reference points
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                funAtRefCoordsOriented(i,j,d) += funDofsOriented(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);

            error = std::max(std::abs( funAtPhysRefCoords(i,j,d) - funAtRefCoordsOriented(i,j,d)), error);
          }
          if(error>100*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function values at reference points differ from those computed using basis functions of Tet " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtPhysRefCoords(i,j,0) << "," << funAtPhysRefCoords(i,j,1) << ", " << funAtPhysRefCoords(i,j,2) << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtRefCoordsOriented(i,j,0) << "," << funAtRefCoordsOriented(i,j,1) << ", " << funAtRefCoordsOriented(i,j,2) << ")";
            *outStream << std::endl;
          }
        }

        //check that fun values are consistent on common face dofs
        {
          bool areDifferent(false);
          for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofsOriented(0,basis.getDofOrdinal(2,faceIndex[0],j))
                - funDofsOriented(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
          }
          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
            *outStream << "Function DOFs for Tet 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofsOriented(0,basis.getDofOrdinal(2,faceIndex[0],j));
            *outStream << "\nFunction DOFs for Tet 1 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofsOriented(1,basis.getDofOrdinal(2,faceIndex[1],j));
            *outStream << std::endl;
          }
        }

        //check that fun values are consistent on common edges dofs
        {
          bool areDifferent(false);
          for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
            for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
              areDifferent = std::abs(funDofsOriented(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                  - funDofsOriented(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
            }
            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common edge " << iEdge << " computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
              *outStream << "Function DOFs for Tet 0 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funDofsOriented(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
              *outStream << "\nFunction DOFs for Tet 1 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funDofsOriented(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j));
              *outStream << std::endl;
            }
          }
        }
        std::rotate(&tets_orig[0][0], &tets_orig[0][0]+1, &tets_orig[0][0]+4);
      }
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

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

