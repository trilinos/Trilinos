// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief Test interpolation and projection capabilities for Hexaedral elements

    The test considers two hexahedra in the physical space sharing a common face. 
    In order to test significant configurations, we consider 6 rotations of the reference hexahedron
    to the first (physical) hexahedron, so that the common face is mapped from all the 6 faces
    of the reference hexahedron.
    Then, for each of the mappings, the global ids of the vertices of the common face are permuted (8 permutations).

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
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"


//this allows to reduce cost of the tests.
//undefine when debugging/developing
#define RANDOMLY_PICK_ELEM_PERMUTATION

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
int InterpolationProjectionHex(const bool verbose) {

  using ExecSpaceType = typename DeviceType::execution_space;
  using MemSpaceType = typename DeviceType::memory_space;

  using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

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

  using DynRankViewIntHost = Kokkos::DynRankView<ordinal_type,HostSpaceType>;

  *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
  *outStream << "\n";

  int errorFlag = 0;
  const ValueType tol = tolerence();

  bool pickTest = false;
  int elemPermutation=0, sharedSidePermutation=0;

#ifdef  RANDOMLY_PICK_ELEM_PERMUTATION
  pickTest = true;
  /* initialize random seed: */
  std::srand (std::time(NULL));
  int configuration = std::rand() % 48;
  elemPermutation = configuration % 6;
  sharedSidePermutation = configuration/6;
  *outStream << "Randomly picked configuration (cellTopo premutation, shared face permutation): (" << elemPermutation << ", " <<sharedSidePermutation << ")" << std::endl;
#endif

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

  using edgeType = std::array<ordinal_type,2>;
  using faceType = std::array<ordinal_type,4>;
  using ct = CellTools<DeviceType>;
  using ots = OrientationTools<DeviceType>;
  using fst = FunctionSpaceTools<DeviceType>;
  using pts = ProjectionTools<DeviceType>;
  using li = LagrangianInterpolation<DeviceType>;
  using ProjStruct = ProjectionStruct<DeviceType,ValueType>;
  using lt = LagrangianTools<DeviceType>;
  using  basisType = Basis<DeviceType,ValueType,ValueType>;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 8;
  constexpr ordinal_type numTotalVertexes = 12;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}, {-1,-1,2},{1,-1,2},{1,1,2},{-1,1,2}};
  ordinal_type cells_orig[numCells][numElemVertexes] = {{0,1,2,3,4,5,6,7},{4,5,6,7,8,9,10,11}};  //common face is {4,5,6,7}
  faceType common_face = {{4,5,6,7}};
  faceType faceLeft = {{0, 3, 7, 4}};
  faceType faceRight = {{1, 2, 6, 5}};
  faceType faceFront = {{0, 4, 5, 1}};
  faceType faceBack = {{2, 3, 7, 6}};
  ordinal_type cells_rotated[numCells][numElemVertexes];
  faceType faceLeftOriented, faceRightOriented, faceBackOriented, faceFrontOriented;

  std::set<edgeType> common_edges;
  common_edges.insert(edgeType({{4,5}})); common_edges.insert(edgeType({{5,6}})); common_edges.insert(edgeType({{6,7}})); common_edges.insert(edgeType({{4,7}}));
  const ordinal_type max_degree = 4;

  using CG_NBasis = NodalBasisFamily<DeviceType,ValueType,ValueType>;
  using CG_DNBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;
  std::vector<basisType*> basis_set;

  shards::CellTopology cellTopo(shards::getCellTopologyData<shards::Hexahedron<8> >());
  ordinal_type numNodesPerElem = cellTopo.getNodeCount();

  using faceShape = shards::Quadrilateral<4>;
  const CellTopologyData * const faceTopoData = shards::getCellTopologyData<faceShape >();

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Test 1 (Orientation - HGRAD)                                |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";


  try {

    //reordering of nodes to explore different orientations

    for(int sharedSideCount = 0; sharedSideCount<faceShape::permutation_count; sharedSideCount++) {
      if((sharedSideCount != sharedSidePermutation) && pickTest)
        continue;

      ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};
      ordinal_type orderback[numTotalVertexes];
      
      for ( unsigned i = 0; i < faceShape::node_count; ++i ) {
        reorder[common_face[i]] = common_face[faceTopoData->permutation[sharedSideCount].node[i]];
      }

      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type cells[numCells][numElemVertexes];
      std::copy(&cells_orig[0][0], &cells_orig[0][0]+numCells*numElemVertexes, &cells_rotated[0][0]);
      for (ordinal_type shift=0; shift<6; ++shift) {
        if(pickTest && (shift != elemPermutation))
          continue;

        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceLeft[ii]] = cells_orig[0][faceLeftOriented[ii]];
            cells_rotated[0][faceRight[ii]] = cells_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceFront[ii]] = cells_orig[0][faceFrontOriented[ii]];
            cells_rotated[0][faceBack[ii]] = cells_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            cells[i][j] = reorder[cells_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[1][j] << " ";
        *outStream << "]\n";

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hostPhysVertexes(i,j,k) = vertices[cells[i][j]][k];
        deep_copy(physVertexes, hostPhysVertexes);

        //computing common face and edges
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][4];

        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<cellTopo.getSideCount(); ++is) {
              for (std::size_t k=0; k<cellTopo.getNodeCount(2,is); ++k)
                face[k]= cells_rotated[i][cellTopo.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<cellTopo.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<cellTopo.getNodeCount(1,ie); ++k)
                edge[k]= cells_rotated[i][cellTopo.getNodeMap(1,ie,k)];
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
        DynRankViewIntHost elemNodesHost(&cells[0][0], numCells, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, cellTopo);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {
          basis_set.clear();
          if(degree==1)
            basis_set.push_back(new Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType>());
          if(degree==2)
            basis_set.push_back(new Basis_HGRAD_HEX_C2_FEM<DeviceType,ValueType,ValueType>());
          basis_set.push_back(new typename  CG_NBasis::HGRAD_HEX(degree,POINTTYPE_EQUISPACED));
          basis_set.push_back(new typename  CG_DNBasis::HGRAD_HEX(degree,POINTTYPE_WARPBLEND));

          for (auto basisPtr:basis_set) {

            auto name = basisPtr->getName();
            *outStream << " " << name <<  ": " << degree << std::endl;
            ordinal_type basisCardinality = basisPtr->getCardinality();

            //compute DofCoords Oriented
            DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
            DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
            DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality);
            DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

            //compute Lagrangian Interpolation of fun
            {
              basisPtr->getDofCoords(dofCoords);
              //Compute physical Dof Coordinates
              DynRankView ConstructWithLabel(linearBasisValuesAtDofCoords, numNodesPerElem, basisCardinality);
              Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> linearBasis;
              linearBasis.getValues(linearBasisValuesAtDofCoords,dofCoords);

              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &i) {
                Fun fun;
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    for(ordinal_type d=0; d<dim; ++d)
                      physDofCoords(i,j,d) += physVertexes(i,k,d)*linearBasisValuesAtDofCoords(k,j);

                  funAtDofCoords(i,j) += fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2));
                }
              });

              li::getBasisCoeffs(basisCoeffsLI, funAtDofCoords, basisPtr, elemOrts);
              Kokkos::fence();
            }

            //Testing Kronecker property of basis functions
            {
              DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
              DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality);
              lt::getOrientedDofCoords(dofCoordsOriented, basisPtr, elemOrts);
              lt::getOrientedDofCoeffs(dofCoeffsPhys, basisPtr, elemOrts);
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
              for(ordinal_type i=0; i<numCells; ++i) {
                auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL());
                basisPtr->getValues(outView, inView);
              }

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr);

              auto hostBasisValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisValuesAtDofCoordsOriented);
              auto hostDofCoeffsPhys = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dofCoeffsPhys);
              ExecSpaceType().fence();
              for(ordinal_type i=0; i<numCells; ++i) {
                for(ordinal_type k=0; k<basisCardinality; ++k) {
                  for(ordinal_type j=0; j<basisCardinality; ++j){
                    ValueType dofValue = hostBasisValues(i,k,j) * hostDofCoeffsPhys(i,j);
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
            auto hostBasisCoeffsLI = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsLI);
            {
              bool areDifferent(false);
              auto numFaceDOFs = basisPtr->getDofCount(2,0);
              for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
                areDifferent = std::abs(hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                    - hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
              }

              if(areDifferent) {
                auto hostPhysDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), physDofCoords);
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                *outStream << "Function DOFs for Hex 0 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
                *outStream << "\nFunction DOFs for Hex 1 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
                *outStream << std::endl;
              }
            }

            //check that fun values are consistent on common edges dofs
            {
              bool areDifferent(false);
              auto numEdgeDOFs = basisPtr->getDofCount(1,0);
              for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
                for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
                  areDifferent = std::abs(hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                      - hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
                }
                if(areDifferent)
                {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                  *outStream << "Function DOFs for Hex 0 are:";
                  for(ordinal_type j=0;j<numEdgeDOFs;j++)
                    *outStream << " " << hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j));
                  *outStream << "\nFunction DOFs for Hex 1 are:";
                  for(ordinal_type j=0;j<numEdgeDOFs;j++)
                    *outStream << " " << hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j));
                  *outStream << std::endl;
                }
              }
            }

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
            DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
            DynRankView basisValuesAtDofCoords("inValues", basisCardinality, basisCardinality);

            basisPtr->getValues(basisValuesAtDofCoords, dofCoords);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                basisValuesAtDofCoords,
                elemOrts,
                basisPtr);

            // transform basis (pullback)
            fst::HGRADtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
                basisValuesAtDofCoordsOriented);

            DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &i) {
              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(ordinal_type k=0; k<basisCardinality; ++k)
                  funAtDofCoordsOriented(i,j) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j);
            });

            auto hostFunAtDofCoordsOriented = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoordsOriented);
            auto hostFunAtDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoords);
            for(ordinal_type i=0; i<numCells; ++i) {
              ValueType error=0;
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                error = std::max(std::abs( hostFunAtDofCoords(i,j) - hostFunAtDofCoordsOriented(i,j)), error);
              }

              if(error>100*tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
                *outStream << "Function values at reference points are:\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoords(i,j)  << ")";
                *outStream << "\nFunction values at reference points computed using basis functions are\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoordsOriented(i,j)  << ")";
                *outStream << std::endl;
              }
            }

            //compute projection-based interpolation of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsHGrad, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree());

              ProjStruct projStruct;
              projStruct.createHGradProjectionStruct(basisPtr, targetCubDegree, targetDerivCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              auto evaluationGradPoints = projStruct.getAllDerivEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0), numGradPoints = evaluationGradPoints.extent(0);

              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
              DynRankView ConstructWithLabel(targetGradAtEvalPoints, numCells, numGradPoints, dim);

              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
              basisPtr->getValues(hgradBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hgradBasisAtEvaluationPoints,
                  hgradBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              DynRankView ConstructWithLabel(gradOfHGradBasisAtEvaluationPoints, numCells, basisCardinality , numGradPoints, dim);
              if(numGradPoints>0) {
                DynRankView ConstructWithLabel(gradOfHGradBasisAtEvaluationPointsNonOriented, basisCardinality , numGradPoints, dim);
                basisPtr->getValues(gradOfHGradBasisAtEvaluationPointsNonOriented,evaluationGradPoints, OPERATOR_GRAD);
                ots::modifyBasisByOrientation(gradOfHGradBasisAtEvaluationPoints,
                    gradOfHGradBasisAtEvaluationPointsNonOriented,
                    elemOrts,
                    basisPtr);
              }


              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hgradBasisAtEvaluationPoints(ic,k,i);
                }
                for(int i=0;i<numGradPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetGradAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*gradOfHGradBasisAtEvaluationPoints(ic,k,i,d);//funHGradCoeffs(k)
                }
              });

              pts::getHGradBasisCoeffs(basisCoeffsHGrad,
                  targetAtEvalPoints,
                  targetGradAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian nterpolation are the same as those of the projection-based interpolation
            {
              ValueType diffErr(0);
              auto hostBasisCoeffsHGrad = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsHGrad);

              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsHGrad(ic,k)));
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

              ProjStruct projStruct;
              projStruct.createL2ProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);


              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
              basisPtr->getValues(hgradBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hgradBasisAtEvaluationPoints,
                  hgradBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);


              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hgradBasisAtEvaluationPoints(ic,k,i);
                }
              });

              pts::getL2BasisCoeffs(basisCoeffsL2,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }
            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr =0;
              auto hostBasisCoeffsL2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2);

              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2(ic,k)));
              }

              if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HGRAD_C" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the functon."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }

            //compute DG L2 projection of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsL2DG, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree());

              ProjStruct projStruct;
              projStruct.createL2DGProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);


              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);
              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
              DynRankView ConstructWithLabel(hgradBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
              basisPtr->getValues(hgradBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hgradBasisAtEvaluationPoints,
                  hgradBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);


              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {

                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hgradBasisAtEvaluationPoints(ic,k,i);
                }
              });

              pts::getL2DGBasisCoeffs(basisCoeffsL2DG,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr =0;
              auto hostBasisCoeffsL2DG = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2DG);
              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2DG(ic,k)));
              }

              if(diffErr > 1e4*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HGRAD_C" << degree << ": The weights recovered with the L2DG optimization are different than the one used for generating the functon."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }
            delete basisPtr;
          }
        }
      }
    } //reorder vertices of common face

  } catch (std::exception &err) {
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

    for(int sharedSideCount = 0; sharedSideCount<faceShape::permutation_count; sharedSideCount++) {
      if((sharedSideCount != sharedSidePermutation) && pickTest)
        continue;

      ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};
      ordinal_type orderback[numTotalVertexes];
      
      for ( unsigned i = 0; i < faceShape::node_count; ++i ) {
        reorder[common_face[i]] = common_face[faceTopoData->permutation[sharedSideCount].node[i]];
      }

      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type cells[numCells][numElemVertexes];
      std::copy(&cells_orig[0][0], &cells_orig[0][0]+numCells*numElemVertexes, &cells_rotated[0][0]);

      for (ordinal_type shift=0; shift<6; ++shift) {
        if(pickTest && (shift != elemPermutation))
          continue;
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceLeft[ii]] = cells_orig[0][faceLeftOriented[ii]];
            cells_rotated[0][faceRight[ii]] = cells_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceFront[ii]] = cells_orig[0][faceFrontOriented[ii]];
            cells_rotated[0][faceBack[ii]] = cells_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            cells[i][j] = reorder[cells_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d){
            vertices[i][d] = vertices_orig[orderback[i]][d];
          }


        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[1][j] << " ";
        *outStream << "]\n";

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hostPhysVertexes(i,j,k) = vertices[cells[i][j]][k];
        deep_copy(physVertexes, hostPhysVertexes);

        //computing edges and tangents
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][4];

        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<cellTopo.getSideCount(); ++is) {
              for (std::size_t k=0; k<cellTopo.getNodeCount(2,is); ++k)
                face[k]= cells_rotated[i][cellTopo.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<cellTopo.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<cellTopo.getNodeCount(1,ie); ++k)
                edge[k]= cells_rotated[i][cellTopo.getNodeMap(1,ie,k)];
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
        DynRankViewIntHost elemNodesHost(&cells[0][0], numCells, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, cellTopo);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {

          basis_set.clear();
          if(degree==1)
            basis_set.push_back(new Basis_HCURL_HEX_I1_FEM<DeviceType,ValueType,ValueType>());
          basis_set.push_back(new typename  CG_NBasis::HCURL_HEX(degree,POINTTYPE_WARPBLEND));
          basis_set.push_back(new typename  CG_DNBasis::HCURL_HEX(degree,POINTTYPE_EQUISPACED));

          for (auto basisPtr:basis_set) {

            auto name = basisPtr->getName();
            *outStream << " " << name <<  ": " << degree << std::endl;

            ordinal_type basisCardinality = basisPtr->getCardinality();

            //compute DofCoords Oriented
            DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
            DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
            DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality, dim);
            DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

            //compute Lagrangian Interpolation of fun
            {
              basisPtr->getDofCoords(dofCoords);

              //Compute physical Dof Coordinates

              DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
              ct::setJacobian(jacobian, dofCoords, physVertexes, cellTopo);

              Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> linearBasis;
              DynRankView ConstructWithLabel(linearBasisValuesAtDofCoords, numNodesPerElem, basisCardinality);
              linearBasis.getValues(linearBasisValuesAtDofCoords, dofCoords);

              DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality, dim);
              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &i) {
                FunCurl fun;
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    for(ordinal_type d=0; d<dim; ++d)
                      physDofCoords(i,j,d) += physVertexes(i,k,d)*linearBasisValuesAtDofCoords(k,j);

                  for(ordinal_type k=0; k<dim; ++k)
                    funAtDofCoords(i,j,k) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k);
                  for(ordinal_type k=0; k<dim; ++k)
                    for(ordinal_type d=0; d<dim; ++d)
                      fwdFunAtDofCoords(i,j,k) += jacobian(i,j,d,k)*funAtDofCoords(i,j,d);
                }
              });

              li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, basisPtr, elemOrts);
            }


            //Testing Kronecker property of basis functions
            {
              DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
              DynRankView ConstructWithLabel(dofCoeffs, numCells, basisCardinality, dim);
              lt::getOrientedDofCoords(dofCoordsOriented, basisPtr, elemOrts);
              lt::getOrientedDofCoeffs(dofCoeffs, basisPtr, elemOrts);
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
              for(ordinal_type i=0; i<numCells; ++i) {
                auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
                basisPtr->getValues(outView, inView);
              }

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr);

              auto hostBasisValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisValuesAtDofCoordsOriented);
              auto hostDofCoeffs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dofCoeffs);
              for(ordinal_type i=0; i<numCells; ++i) {
                for(ordinal_type k=0; k<basisCardinality; ++k) {
                  for(ordinal_type j=0; j<basisCardinality; ++j){
                    ValueType dofValue=0;
                    for(ordinal_type d=0; d<dim; ++d)
                      dofValue += hostBasisValues(i,k,j,d) * hostDofCoeffs(i,j,d);
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
            auto hostBasisCoeffsLI = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsLI);
            {
              bool areDifferent(false);
              auto numEdgeDOFs = basisPtr->getDofCount(1,0);
              for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
                for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
                  areDifferent = std::abs(hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                      - hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
                }
                if(areDifferent) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                  *outStream << "Function DOFs for Hex 0 are:";
                  for(ordinal_type j=0;j<numEdgeDOFs;j++)
                    *outStream << " " << hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(1,edgeIndexes[0][iEdge],j));
                  *outStream << "\nFunction DOFs for Hex 1 are:";
                  for(ordinal_type j=0;j<numEdgeDOFs;j++)
                    *outStream << " " << hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(1,edgeIndexes[1][iEdge],j));
                  *outStream << std::endl;
                }
              }
            }

            //check that fun values are consistent on common face dofs
            {
              bool areDifferent(false);
              auto numFaceDOFs = basisPtr->getDofCount(2,0);
              for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
                areDifferent = std::abs(hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                    - hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
              }

              if(areDifferent) {
                auto hostPhysDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), physDofCoords);
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                *outStream << "Function DOFs for Hex 0 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
                *outStream << "\nFunction DOFs for Hex 1 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
                *outStream << std::endl;
              }
            }

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            DynRankView basisValuesAtDofCoords("inValues", basisCardinality, basisCardinality, dim);

            basisPtr->getValues(basisValuesAtDofCoords, dofCoords);            

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                basisValuesAtDofCoords,
                elemOrts,
                basisPtr);

            // transform basis values
            DynRankView ConstructWithLabel(jacobianAtDofCoords, numCells, basisCardinality, dim, dim);
            DynRankView ConstructWithLabel(jacobianAtDofCoords_inv, numCells, basisCardinality, dim, dim);
            ct::setJacobian(jacobianAtDofCoords, dofCoords, physVertexes, cellTopo);
            ct::setJacobianInv (jacobianAtDofCoords_inv, jacobianAtDofCoords);
            fst::HCURLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
                jacobianAtDofCoords_inv,
                basisValuesAtDofCoordsOriented);
            DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality, dim);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(ordinal_type d=0; d<dim; ++d) {
                  for(ordinal_type k=0; k<basisCardinality; ++k)
                    funAtDofCoordsOriented(i,j,d) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j,d);
                }
            });

            auto hostFunAtDofCoordsOriented = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoordsOriented);
            auto hostFunAtDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoords);
            for(ordinal_type i=0; i<numCells; ++i) {
              ValueType error=0;
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                for(ordinal_type d=0; d<dim; ++d)
                  error = std::max(std::abs( hostFunAtDofCoords(i,j,d) - hostFunAtDofCoordsOriented(i,j,d)), error);
              }

              if(error>100*tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
                *outStream << "Function values at reference points are:\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoords(i,j,0) << "," << hostFunAtDofCoords(i,j,1) << ", " << hostFunAtDofCoords(i,j,2) << ")";
                *outStream << "\nFunction values at reference points computed using basis functions are\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoordsOriented(i,j,0) << "," << hostFunAtDofCoordsOriented(i,j,1) << ", " << hostFunAtDofCoordsOriented(i,j,2) << ")";
                *outStream << std::endl;
              }
            }


            //compute projection-based interpolation of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsHCurl, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree()-1);

              ProjStruct projStruct;
              projStruct.createHCurlProjectionStruct(basisPtr, targetCubDegree, targetDerivCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              auto evaluationCurlPoints = projStruct.getAllDerivEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0), numCurlPoints = evaluationCurlPoints.extent(0);


              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
              DynRankView ConstructWithLabel(targetCurlAtEvalPoints, numCells, numCurlPoints, dim);

              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hcurlBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hcurlBasisAtEvaluationPoints,
                  hcurlBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              DynRankView ConstructWithLabel(curlOfHCurlBasisAtEvaluationPoints, numCells, basisCardinality , numCurlPoints, dim);
              if(numCurlPoints>0) {
                DynRankView ConstructWithLabel(curlOfHCurlBasisAtEvaluationPointsNonOriented, basisCardinality , numCurlPoints, dim);
                basisPtr->getValues(curlOfHCurlBasisAtEvaluationPointsNonOriented, evaluationCurlPoints, OPERATOR_CURL);
                ots::modifyBasisByOrientation(curlOfHCurlBasisAtEvaluationPoints,
                    curlOfHCurlBasisAtEvaluationPointsNonOriented,
                    elemOrts,
                    basisPtr);
              }


              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
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
              });

              pts::getHCurlBasisCoeffs(basisCoeffsHCurl,
                  targetAtEvalPoints,
                  targetCurlAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
            {
              ValueType diffErr(0);
              auto hostBasisCoeffsHCurl = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsHCurl);
              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsHCurl(ic,k)));
              }


              if(diffErr > pow(10, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HCURL_I" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << " " <<pow(7, degree-1)*tol <<  std::endl;
              }
            }

            //compute L2 projection of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree());

              ProjStruct projStruct;
              projStruct.createL2ProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);


              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hcurlBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hcurlBasisAtEvaluationPoints,
                  hcurlBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hcurlBasisAtEvaluationPoints(ic,k,i,d);
                }
              });

              pts::getL2BasisCoeffs(basisCoeffsL2,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr = 0;
              auto hostBasisCoeffsL2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2);
              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2(ic,k)));
              }

              if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HCURL_I" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the functon."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }

            //compute L2DG projection of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsL2DG, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree());

              ProjStruct projStruct;
              projStruct.createL2DGProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);


              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hcurlBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hcurlBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);

              ots::modifyBasisByOrientation(hcurlBasisAtEvaluationPoints,
                  hcurlBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hcurlBasisAtEvaluationPoints(ic,k,i,d);
                }
              });

              pts::getL2DGBasisCoeffs(basisCoeffsL2DG,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr = 0;
              auto hostBasisCoeffsL2DG = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2DG);
              for(int k=0;k<basisCardinality;k++) {
                for(int ic=0; ic<numCells; ic++)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2DG(ic,k)));
              }

              if(diffErr > 1e4*tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HCURL_I" << degree << ": The weights recovered with the L2DG optimization are different than the one used for generating the functon."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }

            delete basisPtr;
          }
        }
      }
    }  //reorder vertices of common face

  } catch (std::exception &err) {
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

    for(int sharedSideCount = 0; sharedSideCount<faceShape::permutation_count; sharedSideCount++) {
      if((sharedSideCount != sharedSidePermutation) && pickTest)
        continue;

      ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};
      ordinal_type orderback[numTotalVertexes];
      
      for ( unsigned i = 0; i < faceShape::node_count; ++i ) {
        reorder[common_face[i]] = common_face[faceTopoData->permutation[sharedSideCount].node[i]];
      }

      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type cells[numCells][numElemVertexes];
      std::copy(&cells_orig[0][0], &cells_orig[0][0]+numCells*numElemVertexes, &cells_rotated[0][0]);

      for (ordinal_type shift=0; shift<6; ++shift) {
        if(pickTest && (shift != elemPermutation))
          continue;
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceLeft[ii]] = cells_orig[0][faceLeftOriented[ii]];
            cells_rotated[0][faceRight[ii]] = cells_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            cells_rotated[0][faceFront[ii]] = cells_orig[0][faceFrontOriented[ii]];
            cells_rotated[0][faceBack[ii]] = cells_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            cells[i][j] = reorder[cells_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << cells[1][j] << " ";
        *outStream << "]\n";

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
        auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
        for(ordinal_type i=0; i<numCells; ++i)
          for(ordinal_type j=0; j<numNodesPerElem; ++j)
            for(ordinal_type k=0; k<dim; ++k)
              hostPhysVertexes(i,j,k) = vertices[cells[i][j]][k];
        deep_copy(physVertexes, hostPhysVertexes);


        //computing edges and tangents

        ordinal_type faceIndex[numCells];
        {
          faceType face={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<cellTopo.getSideCount(); ++is) {
              for (std::size_t k=0; k<cellTopo.getNodeCount(2,is); ++k)
                face[k]= cells_rotated[i][cellTopo.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewIntHost elemNodesHost(&cells[0][0], numCells, numElemVertexes);
        auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, cellTopo);

        for (ordinal_type degree=1; degree <= max_degree; degree++) {

          basis_set.clear();
          if(degree==1)
            basis_set.push_back(new Basis_HDIV_HEX_I1_FEM<DeviceType,ValueType,ValueType>());
          basis_set.push_back(new typename  CG_NBasis::HDIV_HEX(degree,POINTTYPE_EQUISPACED));
          basis_set.push_back(new typename  CG_DNBasis::HDIV_HEX(degree,POINTTYPE_WARPBLEND));

          for (auto basisPtr:basis_set) {

            auto name = basisPtr->getName();
            *outStream << " " << name << ": "<< degree << std::endl;
            ordinal_type basisCardinality = basisPtr->getCardinality();

            //compute DofCoords Oriented

            DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
            DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
            DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality, dim);
            DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

            //compute Lagrangian Interpolation of fun
            {
              basisPtr->getDofCoords(dofCoords);
              

              //need to transform dofCoeff to physical space (they transform as normals)
              DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
              DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
              DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
              ct::setJacobian(jacobian, dofCoords, physVertexes, cellTopo);
              ct::setJacobianInv (jacobian_inv, jacobian);
              ct::setJacobianDet (jacobian_det, jacobian);
              
              //Compute physical Dof Coordinates
              DynRankView ConstructWithLabel(linearBasisValuesAtDofCoord, numNodesPerElem, basisCardinality);
              Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> linearBasis;
              linearBasis.getValues(linearBasisValuesAtDofCoord,dofCoords);

              DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality, dim);
              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &i) {
                FunDiv fun;
                for(ordinal_type j=0; j<basisCardinality; ++j){
                  for(ordinal_type k=0; k<numNodesPerElem; ++k)
                    for(ordinal_type d=0; d<dim; ++d)
                      physDofCoords(i,j,d) += physVertexes(i,k,d)*linearBasisValuesAtDofCoord(k,j);

                  for(ordinal_type k=0; k<dim; ++k)
                    funAtDofCoords(i,j,k) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k);
                  for(ordinal_type k=0; k<dim; ++k)
                    for(ordinal_type d=0; d<dim; ++d)
                      fwdFunAtDofCoords(i,j,k) += jacobian_det(i,j)*jacobian_inv(i,j,k,d)*funAtDofCoords(i,j,d);
                }
              });

              li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, basisPtr, elemOrts);
            }

            //Testing Kronecker property of basis functions
            {
              DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
              DynRankView ConstructWithLabel(dofCoeffs, numCells, basisCardinality, dim);
              lt::getOrientedDofCoords(dofCoordsOriented, basisPtr, elemOrts);
              lt::getOrientedDofCoeffs(dofCoeffs, basisPtr, elemOrts);
              DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
              DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
              for(ordinal_type i=0; i<numCells; ++i) {
                auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
                auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
                basisPtr->getValues(outView, inView);
              }

              // modify basis values to account for orientations
              ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                  basisValuesAtDofCoords,
                  elemOrts,
                  basisPtr);

              DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
              DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
              ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, cellTopo);
              ct::setJacobianDet (jacobian_det, jacobian);

              auto hostBasisValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisValuesAtDofCoordsOriented);
              auto hostDofCoeffs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dofCoeffs);
              for(ordinal_type i=0; i<numCells; ++i) {
                for(ordinal_type k=0; k<basisCardinality; ++k) {
                  for(ordinal_type j=0; j<basisCardinality; ++j){
                    ValueType dofValue=0;
                    for(ordinal_type d=0; d<dim; ++d)
                      dofValue += hostBasisValues(i,k,j,d) * hostDofCoeffs(i,j,d);
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
            auto hostBasisCoeffsLI = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsLI);
            {
              bool areDifferent(false);
              auto numFaceDOFs = basisPtr->getDofCount(2,0);
              for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
                areDifferent = std::abs(hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j))
                    - hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
              }

              if(areDifferent) {
                errorFlag++;
                auto hostPhysDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), physDofCoords);
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
                *outStream << "Function DOFs for Hex 0 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(0,basisPtr->getDofOrdinal(2,faceIndex[0],j)) << " | (" << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),0) << "," << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),1) << ", " << hostPhysDofCoords(0,basisPtr->getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
                *outStream << "\nFunction DOFs for Hex 1 are:";
                for(ordinal_type j=0;j<numFaceDOFs;j++)
                  *outStream << " " << hostBasisCoeffsLI(1,basisPtr->getDofOrdinal(2,faceIndex[1],j))<< " | (" << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),0) << "," << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),1) << ", " << hostPhysDofCoords(1,basisPtr->getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
                *outStream << std::endl;
              }
            }

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            DynRankView basisValuesAtDofCoords("inValues", basisCardinality, basisCardinality, dim);

            basisPtr->getValues(basisValuesAtDofCoords, dofCoords);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                basisValuesAtDofCoords,
                elemOrts,
                basisPtr);

            // transform basis values
            DynRankView ConstructWithLabel(jacobianAtDofCoords, numCells, basisCardinality, dim, dim);
            DynRankView ConstructWithLabel(jacobianAtDofCoords_det, numCells, basisCardinality);
            ct::setJacobian(jacobianAtDofCoords, dofCoords, physVertexes, cellTopo);
            ct::setJacobianDet (jacobianAtDofCoords_det, jacobianAtDofCoords);
            fst::HDIVtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
                jacobianAtDofCoords,
                jacobianAtDofCoords_det,
                basisValuesAtDofCoordsOriented);
            DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality, dim);
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
            KOKKOS_LAMBDA (const int &i) {
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                for(ordinal_type d=0; d<dim; ++d)
                  for(ordinal_type k=0; k<basisCardinality; ++k)
                    funAtDofCoordsOriented(i,j,d) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j,d);
              }
            });

            auto hostFunAtDofCoordsOriented = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoordsOriented);
            auto hostFunAtDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoords);
            for(ordinal_type i=0; i<numCells; ++i) {
              ValueType error=0;
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                for(ordinal_type d=0; d<dim; ++d)
                  error = std::max(std::abs( hostFunAtDofCoords(i,j,d) - hostFunAtDofCoordsOriented(i,j,d)), error);
              }

              if(error>100*tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
                *outStream << "Function values at reference points are:\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoords(i,j,0) << "," << hostFunAtDofCoords(i,j,1) << ", " << hostFunAtDofCoords(i,j,2) << ")";
                *outStream << "\nFunction values at reference points computed using basis functions are\n";
                for(ordinal_type j=0; j<basisCardinality; ++j)
                  *outStream << " (" << hostFunAtDofCoordsOriented(i,j,0) << "," << hostFunAtDofCoordsOriented(i,j,1) << ", " << hostFunAtDofCoordsOriented(i,j,2) << ")";
                *outStream << std::endl;
              }
            }

            //compute projection-based interpolation of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsHDiv, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree()),targetDerivCubDegree(basisPtr->getDegree()-1);

              ProjStruct projStruct;
              projStruct.createHDivProjectionStruct(basisPtr, targetCubDegree, targetDerivCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              auto evaluationDivPoints = projStruct.getAllDerivEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0), numDivPoints = evaluationDivPoints.extent(0);

              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);
              DynRankView ConstructWithLabel(targetDivAtEvalPoints, numCells, numDivPoints);

              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hdivBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hdivBasisAtEvaluationPoints,
                  hdivBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              DynRankView ConstructWithLabel(divOfHDivBasisAtEvaluationPoints, numCells, basisCardinality , numDivPoints);
              if(numDivPoints>0) {
                DynRankView ConstructWithLabel(divOfHDivBasisAtEvaluationPointsNonOriented, basisCardinality , numDivPoints);
                basisPtr->getValues(divOfHDivBasisAtEvaluationPointsNonOriented, evaluationDivPoints, OPERATOR_DIV);
                ots::modifyBasisByOrientation(divOfHDivBasisAtEvaluationPoints,
                    divOfHDivBasisAtEvaluationPointsNonOriented,
                    elemOrts,
                    basisPtr);
              }

              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hdivBasisAtEvaluationPoints(ic,k,i,d);
                }
                for(int i=0;i<numDivPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    targetDivAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*divOfHDivBasisAtEvaluationPoints(ic,k,i);//basisCoeffsLI(k)
                }
              });

              pts::getHDivBasisCoeffs(basisCoeffsHDiv,
                  targetAtEvalPoints,
                  targetDivAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
            {
              ValueType diffErr(0);
              auto hostBasisCoeffsHDiv = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsHDiv);
              for(int k=0;k<basisCardinality;k++) {
                //std::cout << "["<< basisCoeffsLI(0,k) << " " <<  basisCoeffsHDiv(0,k) << "] [" << basisCoeffsLI(1,k) << " " <<  basisCoeffsHDiv(1,k) << "]" <<std::endl;
                for(ordinal_type ic=0; ic<numCells; ++ic)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsHDiv(ic,k)));
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

              ProjStruct projStruct;
              projStruct.createL2ProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);

              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hdivBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hdivBasisAtEvaluationPoints,
                  hdivBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);

              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hdivBasisAtEvaluationPoints(ic,k,i,d);
                }
              });

              pts::getL2BasisCoeffs(basisCoeffsL2,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr = 0;
              auto hostBasisCoeffsL2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2);
              for(int k=0;k<basisCardinality;k++) {
                //std::cout << "["<< basisCoeffsLI(0,k) << " " <<  basisCoeffsHDiv(0,k) << "] [" << basisCoeffsLI(1,k) << " " <<  basisCoeffsHDiv(1,k) << "]" <<std::endl;
                for(ordinal_type ic=0; ic<numCells; ++ic)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2(ic,k)));
              }

              if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HDIV_I" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the function."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }

            //compute DG L2 projection interpolation of the Lagrangian interpolation
            DynRankView ConstructWithLabel(basisCoeffsL2DG, numCells, basisCardinality);
            {
              ordinal_type targetCubDegree(basisPtr->getDegree());

              ProjStruct projStruct;
              projStruct.createL2DGProjectionStruct(basisPtr, targetCubDegree);

              auto evaluationPoints = projStruct.getAllEvalPoints();
              ordinal_type numPoints = evaluationPoints.extent(0);

              DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints, dim);

              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPoints, numCells, basisCardinality , numPoints, dim);
              DynRankView ConstructWithLabel(hdivBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints, dim);
              basisPtr->getValues(hdivBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
              ots::modifyBasisByOrientation(hdivBasisAtEvaluationPoints,
                  hdivBasisAtEvaluationPointsNonOriented,
                  elemOrts,
                  basisPtr);
                  
              Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
              KOKKOS_LAMBDA (const int &ic) {
                for(int i=0;i<numPoints;i++) {
                  for(int k=0;k<basisCardinality;k++)
                    for(int d=0;d<dim;d++)
                      targetAtEvalPoints(ic,i,d) += basisCoeffsLI(ic,k)*hdivBasisAtEvaluationPoints(ic,k,i,d);
                }
              });

              pts::getL2DGBasisCoeffs(basisCoeffsL2DG,
                  targetAtEvalPoints,
                  elemOrts,
                  basisPtr,
                  &projStruct);
            }

            //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
            {
              ValueType diffErr = 0;
              auto hostBasisCoeffsL2DG = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2DG);
              for(int k=0;k<basisCardinality;k++) {
                //std::cout << "["<< basisCoeffsLI(0,k) << " " <<  basisCoeffsHDiv(0,k) << "] [" << basisCoeffsLI(1,k) << " " <<  basisCoeffsHDiv(1,k) << "]" <<std::endl;
                for(ordinal_type ic=0; ic<numCells; ++ic)
                  diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2DG(ic,k)));
              }

              if(diffErr > pow(7, degree-1)*tol) { //heuristic relation on how round-off error depends on degree
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "HDIV_I" << degree << ": The weights recovered with the L2DG optimization are different than the one used for generating the function."<<
                    "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
              }
            }

            delete basisPtr;
          }
        }
      }
    }  //reorder vertices of common face

  } catch (std::exception &err) {
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
    ordinal_type cells[numCells][numElemVertexes];

    for(ordinal_type i=0; i<numCells;++i)
      for(ordinal_type j=0; j<numElemVertexes;++j)
        cells[i][j] = cells_orig[i][j];

    for(ordinal_type i=0; i<numTotalVertexes;++i)
      for(ordinal_type d=0; d<dim;++d)
        vertices[i][d] = vertices_orig[i][d];

    *outStream <<  "Considering Hex 0: [ ";
    for(ordinal_type j=0; j<numElemVertexes;++j)
      *outStream << cells[0][j] << " ";
    *outStream << "] and Hex 1: [ ";
    for(ordinal_type j=0; j<numElemVertexes;++j)
      *outStream << cells[1][j] << " ";
    *outStream << "]\n";

    //computing vertices coords
    DynRankView ConstructWithLabel(physVertexes, numCells, numNodesPerElem, dim);
    auto hostPhysVertexes = Kokkos::create_mirror_view(physVertexes);
    for(ordinal_type i=0; i<numCells; ++i)
      for(ordinal_type j=0; j<numNodesPerElem; ++j)
        for(ordinal_type k=0; k<dim; ++k)
          hostPhysVertexes(i,j,k) = vertices[cells[i][j]][k];
    deep_copy(physVertexes, hostPhysVertexes);

    // compute orientations for cells (one time computation)
    DynRankViewIntHost elemNodesHost(&cells[0][0], numCells, numElemVertexes);
    auto elemNodes = Kokkos::create_mirror_view_and_copy(MemSpaceType(),elemNodesHost);
    Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
    ots::getOrientation(elemOrts, elemNodes, cellTopo);

    for (ordinal_type degree=1; degree <= max_degree; degree++) {

      basis_set.clear();
      if(degree==1)
        basis_set.push_back(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(cellTopo));
      basis_set.push_back(new typename  CG_NBasis::HVOL_HEX(degree,POINTTYPE_EQUISPACED));
      basis_set.push_back(new typename  CG_DNBasis::HVOL_HEX(degree,POINTTYPE_WARPBLEND));

      for (auto basisPtr:basis_set) {

        auto name = basisPtr->getName();
        *outStream << " " << name << ": "<< degree << std::endl;

        ordinal_type basisCardinality = basisPtr->getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(funAtDofCoords, numCells, basisCardinality);
        DynRankView ConstructWithLabel(basisCoeffsLI, numCells, basisCardinality);

        //compute Lagrangian Interpolation of fun
        {
          basisPtr->getDofCoords(dofCoords);

          //need to transform dofCoeff to physical space (they transform as normals)
          DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
          DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
          ct::setJacobian(jacobian, dofCoords, physVertexes, cellTopo);
          ct::setJacobianDet (jacobian_det, jacobian);

          //Compute physical Dof Coordinates
          DynRankView ConstructWithLabel(linearBasisValuesAtDofCoord, numNodesPerElem, basisCardinality);
          Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> linearBasis;
          linearBasis.getValues(linearBasisValuesAtDofCoord, dofCoords);

          DynRankView ConstructWithLabel(fwdFunAtDofCoords, numCells, basisCardinality);
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &i) {
            Fun fun;
            for(ordinal_type j=0; j<basisCardinality; ++j){
              for(ordinal_type k=0; k<numNodesPerElem; ++k)
                for(ordinal_type d=0; d<dim; ++d)
                  physDofCoords(i,j,d) += physVertexes(i,k,d)*linearBasisValuesAtDofCoord(k,j);

              funAtDofCoords(i,j) = fun(degree, physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2));
              fwdFunAtDofCoords(i,j) = jacobian_det(i,j)*funAtDofCoords(i,j);
            }
          });

          li::getBasisCoeffs(basisCoeffsLI, fwdFunAtDofCoords, basisPtr, elemOrts);
        }

        //Testing Kronecker property of basis functions
        {
          DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
          DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality);
          lt::getOrientedDofCoords(dofCoordsOriented, basisPtr, elemOrts);
          lt::getOrientedDofCoeffs(dofCoeffsPhys, basisPtr, elemOrts);
          DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality);
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i) {
            auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
            auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL());
            basisPtr->getValues(outView, inView);
          }

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoords,
              elemOrts,
              basisPtr);

          auto hostBasisValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisValuesAtDofCoordsOriented);
          auto hostDofCoeffsPhys = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dofCoeffsPhys);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type k=0; k<basisCardinality; ++k) {
              for(ordinal_type j=0; j<basisCardinality; ++j){
                ValueType dofValue = hostBasisValues(i,k,j) * hostDofCoeffsPhys(i,j);
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
        DynRankView basisValuesAtDofCoords("inValues", basisCardinality, basisCardinality);

        basisPtr->getValues(basisValuesAtDofCoords, dofCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
            basisValuesAtDofCoords,
            elemOrts,
            basisPtr);

        // transform basis values
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
        ct::setJacobian(jacobian, dofCoords, physVertexes, cellTopo);
        ct::setJacobianDet (jacobian_det, jacobian);
        fst::HVOLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
            jacobian_det,
            basisValuesAtDofCoordsOriented);

        DynRankView ConstructWithLabel(funAtDofCoordsOriented, numCells, basisCardinality);
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
        KOKKOS_LAMBDA (const int &i) {
          for(ordinal_type j=0; j<basisCardinality; ++j) {
            for(ordinal_type k=0; k<basisCardinality; ++k)
              funAtDofCoordsOriented(i,j) += basisCoeffsLI(i,k)*transformedBasisValuesAtDofCoordsOriented(i,k,j);
          }
        });

        auto hostFunAtDofCoordsOriented = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoordsOriented);
        auto hostFunAtDofCoords = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), funAtDofCoords);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<basisCardinality; ++j) {
            error = std::max(std::abs( hostFunAtDofCoords(i,j) - hostFunAtDofCoordsOriented(i,j)), error);
          }

          if(error>100*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<basisCardinality; ++j)
              *outStream << " (" << hostFunAtDofCoords(i,j)  << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<basisCardinality; ++j)
              *outStream << " (" << hostFunAtDofCoordsOriented(i,j)  << ")";
            *outStream << std::endl;
          }
        }

        //compute projection-based interpolation of the Lagrangian interpolation
        DynRankView ConstructWithLabel(basisCoeffsHVol, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(basisPtr->getDegree());

          ProjStruct projStruct;
          projStruct.createHVolProjectionStruct(basisPtr, targetCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
          basisPtr->getValues(hvolBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
          ots::modifyBasisByOrientation(hvolBasisAtEvaluationPoints,
              hvolBasisAtEvaluationPointsNonOriented,
              elemOrts,
              basisPtr);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &ic) {
            for(int i=0;i<numPoints;i++) {
              for(int k=0;k<basisCardinality;k++)
                targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hvolBasisAtEvaluationPoints(ic,k,i);
            }
          });

          pts::getHVolBasisCoeffs(basisCoeffsHVol,
              targetAtEvalPoints,
              elemOrts,
              basisPtr,
              &projStruct);
        }

        //check that the basis coefficients of the Lagrangian interpolation are the same as those of the projection-based interpolation
        auto hostBasisCoeffsLI = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsLI);
        {
          ValueType diffErr(0);
          auto hostBasisCoeffsHVol = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsHVol);
          for(int k=0;k<basisCardinality;k++) {
            for(int ic=0; ic<numCells; ic++)
              diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsHVol(ic,k)));
          }

          //Check that the two representations of the gradient of ifun are consistent
          if(diffErr > pow(20, degree)*tol) { //heuristic relation on how round-off error depends on degree
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "HVOL_C" << degree << ": The weights recovered with the optimization are different than the one used for generating the functon."<<
                "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
          }
        }

        //compute L2 projection of the Lagrangian interpolation
        DynRankView ConstructWithLabel(basisCoeffsL2, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(basisPtr->getDegree());

          ProjStruct projStruct;
          projStruct.createL2ProjectionStruct(basisPtr, targetCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);


          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
          basisPtr->getValues(hvolBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
          ots::modifyBasisByOrientation(hvolBasisAtEvaluationPoints,
              hvolBasisAtEvaluationPointsNonOriented,
              elemOrts,
              basisPtr);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &ic) {
            for(int i=0;i<numPoints;i++) {
              for(int k=0;k<basisCardinality;k++)
                targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hvolBasisAtEvaluationPoints(ic,k,i);
            }
          });

          pts::getL2BasisCoeffs(basisCoeffsL2,
              targetAtEvalPoints,
              elemOrts,
              basisPtr,
              &projStruct);
        }

        //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
        {
          ValueType diffErr = 0;
          auto hostBasisCoeffsL2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2);
          for(int k=0;k<basisCardinality;k++) {
            for(int ic=0; ic<numCells; ic++)
              diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2(ic,k)));
          }

          if(diffErr > pow(20, degree)*tol) { //heuristic relation on how round-off error depends on degree
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "HVOL_C" << degree << ": The weights recovered with the L2 optimization are different than the one used for generating the functon."<<
                "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
          }
        }

        //compute DG L2 projection of the Lagrangian interpolation
        DynRankView ConstructWithLabel(basisCoeffsL2DG, numCells, basisCardinality);
        {
          ordinal_type targetCubDegree(basisPtr->getDegree());

          ProjStruct projStruct;
          projStruct.createL2DGProjectionStruct(basisPtr, targetCubDegree);

          auto evaluationPoints = projStruct.getAllEvalPoints();
          ordinal_type numPoints = evaluationPoints.extent(0);

          DynRankView ConstructWithLabel(targetAtEvalPoints, numCells, numPoints);

          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPoints, numCells, basisCardinality , numPoints);
          DynRankView ConstructWithLabel(hvolBasisAtEvaluationPointsNonOriented, basisCardinality , numPoints);
          basisPtr->getValues(hvolBasisAtEvaluationPointsNonOriented, evaluationPoints, OPERATOR_VALUE);
          ots::modifyBasisByOrientation(hvolBasisAtEvaluationPoints,
              hvolBasisAtEvaluationPointsNonOriented,
              elemOrts,
              basisPtr);

          Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
          KOKKOS_LAMBDA (const int &ic) {
            for(int i=0;i<numPoints;i++) {
              for(int k=0;k<basisCardinality;k++)
                targetAtEvalPoints(ic,i) += basisCoeffsLI(ic,k)*hvolBasisAtEvaluationPoints(ic,k,i);
            }
          });

          pts::getL2DGBasisCoeffs(basisCoeffsL2DG,
              targetAtEvalPoints,
              basisPtr,
              &projStruct);
        }

        //check that the basis coefficients of the Lagrangian interpolation are the same as those of the L2 projection
        {
          ValueType diffErr = 0;
          auto hostBasisCoeffsL2DG = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisCoeffsL2DG);
          for(int k=0;k<basisCardinality;k++) {
            for(int ic=0; ic<numCells; ic++)
              diffErr = std::max(diffErr, std::abs(hostBasisCoeffsLI(ic,k) - hostBasisCoeffsL2DG(ic,k)));
          }

          if(diffErr > pow(20, degree)*tol) { //heuristic relation on how round-off error depends on degree
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "HVOL_C" << degree << ": The weights recovered with the L2DG optimization are different than the one used for generating the functon."<<
                "\nThe max The infinite norm of the difference between the weights is: " <<  diffErr << std::endl;
          }
        }

        delete basisPtr;
      }
    }
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

