// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking orientation tools for Hexahedral elements
 
    The test considers two hexahedra in the physical space sharing a common face. 
    In order to test significant configurations, we consider 6 mappings of the reference hexahedron
    to the first (physical) hexahedron, so that the common face is mapped from all the 6 faces
    of the reference hexahedron.
    Then, for each of the mappings, the global ids of the vertices of the common face are permuted.
    This gives a total of 144 combinations
    
    The test considers HGRAD, HCURL and HDIV elements and for each of them checks the following:
    1. that the physical (oriented) dof coefficients related to faces and edges are proportional 
       to the tangents (for HCURL) of the faces/edges or to the normals (for HDIV) of the faces.
    2. that the Kronecker property of the oriented basis evaluated at the physiscal dof coordinates holds.
    3. that given a function, the dofs of the function located at the common faces/edges are the same
       when computed on the first and second hexahedron.
    4. that a function, belonging to the considered H-space, is exactly reproduced at given points
       using the oriented basis functions and the degrees of freedom of the function.
 
    

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"
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
int OrientationHex(const bool verbose) {

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

  int errorFlag = 0;
  const ValueType tol = tolerence();

  struct Fun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return (x+1)*(y-2);//*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a = 2*x*y+x*x;
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
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a0 = y-7+z*z;
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
  };

  typedef std::array<ordinal_type,2> edgeType;
  typedef std::array<ordinal_type,4> faceType;
  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef Impl::OrientationTools iots;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;
  typedef ArrayTools<DeviceType> at;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 8;
  constexpr ordinal_type numElemEdges = 12;
  constexpr ordinal_type numElemFaces = 6;
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


  *outStream
 << "===============================================================================\n"
 << "|                                                                             |\n"
 << "|                 Test 1 (Orientation - HGRAD)                                |\n"
 << "|                                                                             |\n"
 << "===============================================================================\n";


 try {



   const ordinal_type order = 3;
   ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

   do {
     ordinal_type orderback[numTotalVertexes];
     for(ordinal_type i=0;i<numTotalVertexes;++i) {
       orderback[reorder[i]]=i;
     }
     ValueType vertices[numTotalVertexes][dim];
     ordinal_type hexas[numCells][numElemVertexes];
     std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);
     for (ordinal_type shift=0; shift<6; ++shift) {
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

       shards::CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >());
       shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
       shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

       //computing vertices coords
       DynRankView ConstructWithLabel(physVertexes, numCells, hexa.getNodeCount(), dim);
       for(ordinal_type i=0; i<numCells; ++i)
         for(std::size_t j=0; j<hexa.getNodeCount(); ++j)
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
           for (std::size_t is=0; is<hexa.getSideCount(); ++is) {
             for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
               face[k]= hexas_rotated[i][hexa.getNodeMap(2,is,k)];
             
             //rotate and flip
             auto minElPtr= std::min_element(face.begin(), face.end());
             std::rotate(face.begin(),minElPtr,face.end());
             if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}
             
             if(face == common_face) faceIndex[i]=is;
           }
           //compute edges' tangents
           for (std::size_t ie=0; ie<hexa.getEdgeCount(); ++ie) {
             for (std::size_t k=0; k<hexa.getNodeCount(1,ie); ++k)
               edge[k]= hexas_rotated[i][hexa.getNodeMap(1,ie,k)];
             std::sort(edge.begin(),edge.end());
             auto it=common_edges.find(edge);
             if(it !=common_edges.end()){
               auto edge_lid = std::distance(common_edges.begin(),it);
               edgeIndexes[i][edge_lid]=ie;
             }
           }
         }
       }

       //compute reference points
       Basis_HGRAD_HEX_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
       ordinal_type numRefCoords = warpBasis.getCardinality();
       DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
       warpBasis.getDofCoords(refPoints);

       // compute orientations for cells (one time computation)
       DynRankViewInt elemNodes(&hexas[0][0], 2, numElemVertexes);
       Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
       ots::getOrientation(elemOrts, elemNodes, hexa);

       Basis_HGRAD_HEX_Cn_FEM<DeviceType,ValueType,ValueType> basis(order);
       ordinal_type basisCardinality = basis.getCardinality();

       //compute DofCoords Oriented
       DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
       DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
       basis.getDofCoords(dofCoords);
       {
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type ivertex=0; ivertex<numElemVertexes; ++ivertex) {
             auto idof = basis.getDofOrdinal(0, ivertex, 0);
             for(ordinal_type d=0; d <dim; ++d)
               dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
           }


         Basis_HGRAD_LINE_Cn_FEM<DeviceType,ValueType,ValueType> lineBasis(order);
         ordinal_type lineBasisCardinality = lineBasis.getCardinality();
         DynRankView ConstructWithLabel(lineDofCoords, lineBasisCardinality, dim-2);
         lineBasis.getDofCoords(lineDofCoords);
         DynRankView ConstructWithLabel(lineDofCoordsOriented,  lineBasisCardinality, dim-2);
         DynRankView ConstructWithLabel(lineDofCoordsOriented3d,  lineBasisCardinality, dim);
         auto numEdgeDOFs = basis.getDofCount(1,0);
         for(ordinal_type i=0; i<numCells; ++i) {
           ordinal_type eOrt[numElemEdges];
           elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);
           for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
             iots::mapToModifiedReference(lineDofCoordsOriented,lineDofCoords,line,eOrt[iedge]);
             ct::mapToReferenceSubcell(lineDofCoordsOriented3d, lineDofCoordsOriented, dim-2, iedge, hexa);

             for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
               auto idof = basis.getDofOrdinal(1, iedge, j);
               auto linedof = lineBasis.getDofOrdinal(1, 0, j);

               for(ordinal_type d=0; d <dim; ++d)
                 dofCoordsOriented(i,idof,d) = lineDofCoordsOriented3d(linedof,d);
             }
           }
         }
        
         Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> quadBasis(order);
         ordinal_type quadBasisCardinality = quadBasis.getCardinality();
         ordinal_type numInternalDofs = quadBasis.getDofCount(dim-1,0);
         DynRankView ConstructWithLabel(quadDofCoords, quadBasisCardinality, dim-1);
         DynRankView ConstructWithLabel(quadInternalDofCoords, numInternalDofs, dim-1);
         quadBasis.getDofCoords(quadDofCoords);
         for(ordinal_type i=0; i<numInternalDofs; ++i)
           for(ordinal_type d=0; d <dim-1; ++d)
             quadInternalDofCoords(i,d) = quadDofCoords(quadBasis.getDofOrdinal(dim-1, 0, i),d);

         DynRankView ConstructWithLabel(quadInternalDofCoordsOriented,  numInternalDofs, dim-1);
         DynRankView ConstructWithLabel(quadDofCoordsOriented3d, numInternalDofs, dim);
         ordinal_type fOrt[numElemFaces];
         for(ordinal_type i=0; i<numCells; ++i) {
           elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
           for(ordinal_type iface=0; iface<numElemFaces; iface++) {
             ordinal_type ort = fOrt[iface];
             iots::mapToModifiedReference(quadInternalDofCoordsOriented,quadInternalDofCoords,quad,ort);
             ct::mapToReferenceSubcell(quadDofCoordsOriented3d, quadInternalDofCoordsOriented, dim-1, iface, hexa);
             DynRankView ConstructWithLabel(quadDofs,  quadBasis.getCardinality(), quadBasisCardinality);

             for(ordinal_type j=0; j<numInternalDofs; ++j) {
               auto idof = basis.getDofOrdinal(2, iface, j);
               for(ordinal_type d=0; d <dim; ++d)
                 dofCoordsOriented(i,idof,d) = quadDofCoordsOriented3d(j,d);
             }
           }
         }

         auto numElemDOFs = basis.getDofCount(3,0);
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type j=0; j<numElemDOFs; ++j) {
             auto idof = basis.getDofOrdinal(3, 0, j);
             for(ordinal_type d=0; d <dim; ++d)
               dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
           }
       }

       //Compute physical Dof Coordinates and Reference coordinates
       DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
       DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
       {
         Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> hexaLinearBasis; //used for computing physical coordinates
         DynRankView ConstructWithLabel(hexaLinearBasisValuesAtRefCoords, hexa.getNodeCount(), numRefCoords);
         hexaLinearBasis.getValues(hexaLinearBasisValuesAtRefCoords, refPoints);
         DynRankView ConstructWithLabel(hexaLinearBasisValuesAtDofCoords, numCells, hexa.getNodeCount(), basisCardinality);
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type d=0; d<dim; ++d) {
             for(ordinal_type j=0; j<numRefCoords; ++j)
               for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                 physRefCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtRefCoords(k,j);

             auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
             auto outView =Kokkos::subview( hexaLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
             hexaLinearBasis.getValues(outView, inView);

             for(ordinal_type j=0; j<basisCardinality; ++j)
               for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                 physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtDofCoords(i,k,j);
           }
       }

       DynRankView ConstructWithLabel(dofCoeffs, basisCardinality);
       DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality);

       basis.getDofCoeffs(dofCoeffs);
       rst::clone(dofCoeffsPhys,dofCoeffs);

       //Testing Kronecker property of basis functions
       {
         DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality);
         DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
         DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality);
         for(ordinal_type i=0; i<numCells; ++i) {
           auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
           auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL());
           basis.getValues(outView, inView);
         }

         // modify basis values to account for orientations
         ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                                       basisValuesAtDofCoords,
                                       elemOrts,
                                       &basis);

         // transform basis values
         deep_copy(transformedBasisValuesAtDofCoordsOriented,
                   basisValuesAtDofCoordsOriented);

         for(ordinal_type i=0; i<numCells; ++i) {
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

       //check function reproducibility
       Fun fun;
       DynRankView ConstructWithLabel(funDofs, numCells, basisCardinality);
       DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords);
       for(ordinal_type i=0; i<numCells; ++i) {
         for(ordinal_type j=0; j<numRefCoords; ++j)
             funAtPhysRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
         for(ordinal_type j=0; j<basisCardinality; ++j)
             funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2)) * dofCoeffsPhys(i,j);
       }


       //check that fun values are consistent on common face dofs
       {
         bool areDifferent(false);
         auto numFaceDOFs = basis.getDofCount(2,0);
         for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
           areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j))
               - funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
         }

         if(areDifferent) {
           errorFlag++;
           *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
           *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
           *outStream << "Function DOFs for Hex 0 are:";
           for(ordinal_type j=0;j<numFaceDOFs;j++)
             *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
           *outStream << "\nFunction DOFs for Hex 1 are:";
           for(ordinal_type j=0;j<numFaceDOFs;j++)
             *outStream << " " << funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
           *outStream << std::endl;
         }
       }

       //check that fun values are consistent on common edges dofs
       {
         bool areDifferent(false);
         auto numEdgeDOFs = basis.getDofCount(1,0);
         for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
           for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
             areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                 - funDofs(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
           }
           if(areDifferent) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
             *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
             *outStream << "Function DOFs for Hex 0 are:";
             for(ordinal_type j=0;j<numEdgeDOFs;j++)
               *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
             *outStream << "\nFunction DOFs for Hex 1 are:";
             for(ordinal_type j=0;j<numEdgeDOFs;j++)
               *outStream << " " << funDofs(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j));
             *outStream << std::endl;
           }
         }
       }

       //check that fun values at reference points coincide with those computed using basis functions
       DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
       DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);

       DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
       basis.getValues(basisValuesAtRefCoords, refPoints);

       // modify basis values to account for orientations
       ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
                                     basisValuesAtRefCoords,
                                     elemOrts,
                                     &basis);

       // transform basis values
       deep_copy(transformedBasisValuesAtRefCoordsOriented,
           basisValuesAtRefCoordsOriented);

       DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords);
       for(ordinal_type i=0; i<numCells; ++i) {
         ValueType error=0;
         for(ordinal_type j=0; j<numRefCoords; ++j) {
           for(ordinal_type k=0; k<basisCardinality; ++k)
             funAtRefCoordsOriented(i,j) += funDofs(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j);

           error = std::max(std::abs( funAtPhysRefCoords(i,j) - funAtRefCoordsOriented(i,j)), error);
         }

         if(error>100*tol) {
           errorFlag++;
           *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
           *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
           *outStream << "Function values at reference points are:\n";
           for(ordinal_type j=0; j<numRefCoords; ++j)
             *outStream << " (" << funAtPhysRefCoords(i,j)  << ")";
           *outStream << "\nFunction values at reference points computed using basis functions are\n";
           for(ordinal_type j=0; j<numRefCoords; ++j)
             *outStream << " (" << funAtRefCoordsOriented(i,j)  << ")";
           *outStream << std::endl;
         }
       }

     }
   } while(std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

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

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);
      for (ordinal_type shift=0; shift<6; ++shift) {
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

        shards::CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hexa.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hexa.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];



        //computing edges and tangents
        ValueType edgeTan[numCells][numElemEdges][3];
        ValueType faceTanU[numCells][numElemFaces][3];
        ValueType faceTanV[numCells][numElemFaces][3];
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][4];

        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<hexa.getSideCount(); ++is) {
              for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hexa.getNodeMap(2,is,k)];
             
              //rotate and flip to compute common face
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}
              
              if(face == common_face) faceIndex[i]=is;

              //compute tangents using global ordering
              for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
                face[k]= reorder[face[k]];

              //rotate and flip
              minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}
              
              for(int d=0; d<dim; ++d) {
                faceTanU[i][is][d] = vertices[face[1]][d]-vertices[face[0]][d];
                faceTanV[i][is][d] = vertices[face[3]][d]-vertices[face[0]][d];
              }
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<hexa.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<hexa.getNodeCount(1,ie); ++k)
                edge[k]= hexas_rotated[i][hexa.getNodeMap(1,ie,k)];

              //compute common edge        
              std::sort(edge.begin(),edge.end());
              auto it=common_edges.find(edge);
              if(it !=common_edges.end()){
                auto edge_lid = std::distance(common_edges.begin(),it);
                edgeIndexes[i][edge_lid]=ie;
              }
           
              //compute edge tangent using global numbering
              for (std::size_t k=0; k<hexa.getNodeCount(1,ie); ++k)
                edge[k] = reorder[edge[k]];
              std::sort(edge.begin(),edge.end());
              for(int d=0; d<dim; ++d)
                edgeTan[i][ie][d] = vertices[edge[1]][d]-vertices[edge[0]][d];
            }
          }
        }

        //compute reference points
        Basis_HGRAD_HEX_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], 2, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hexa);


        Basis_HCURL_HEX_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          Basis_HVOL_LINE_Cn_FEM<DeviceType,ValueType,ValueType> lineBasis(order-1,POINTTYPE_GAUSS);
          ordinal_type lineBasisCardinality = lineBasis.getCardinality();
          DynRankView ConstructWithLabel(lineDofCoords, lineBasisCardinality, dim-2);
          lineBasis.getDofCoords(lineDofCoords);
          DynRankView ConstructWithLabel(lineDofCoordsOriented,  lineBasisCardinality, dim-2);
          DynRankView ConstructWithLabel(lineDofCoordsOriented3d,  lineBasisCardinality, dim);
          auto numEdgeDOFs = basis.getDofCount(1,0);
          for(ordinal_type i=0; i<numCells; ++i) {
            ordinal_type eOrt[numElemEdges];
            elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);
            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
              iots::mapToModifiedReference(lineDofCoordsOriented,lineDofCoords,line,eOrt[iedge]);
              ct::mapToReferenceSubcell(lineDofCoordsOriented3d, lineDofCoordsOriented, dim-2, iedge, hexa);

              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                auto linedof = lineBasis.getDofOrdinal(1, 0, j);

                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = lineDofCoordsOriented3d(linedof,d);
              }
            }
          }

          Basis_HCURL_QUAD_In_FEM<DeviceType,ValueType,ValueType> quadBasis(order);
          ordinal_type quadBasisCardinality = quadBasis.getCardinality();
          ordinal_type numInternalDofs = quadBasis.getDofCount(dim-1,0);
          DynRankView ConstructWithLabel(quadDofCoords, quadBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(quadInternalDofCoords, numInternalDofs, dim-1);
          quadBasis.getDofCoords(quadDofCoords);
          for(ordinal_type i=0; i<numInternalDofs; ++i)
            for(ordinal_type d=0; d <dim-1; ++d)
              quadInternalDofCoords(i,d) = quadDofCoords(quadBasis.getDofOrdinal(dim-1, 0, i),d);

          DynRankView ConstructWithLabel(quadInternalDofCoordsOriented,  numInternalDofs, dim-1);
          DynRankView ConstructWithLabel(quadDofCoordsOriented3d, numInternalDofs, dim);
          ordinal_type fOrt[numElemFaces];
          for(ordinal_type i=0; i<numCells; ++i) {
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              ordinal_type ort = fOrt[iface];
              iots::mapToModifiedReference(quadInternalDofCoordsOriented,quadInternalDofCoords,quad,ort);
              ct::mapToReferenceSubcell(quadDofCoordsOriented3d, quadInternalDofCoordsOriented, dim-1, iface, hexa);

              for(ordinal_type j=0; j<numInternalDofs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iface, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = quadDofCoordsOriented3d(j,d);
              }
            }
          }

          auto numElemDOFs = basis.getDofCount(3,0);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type j=0; j<numElemDOFs; ++j) {
              auto idof = basis.getDofOrdinal(3, 0, j);
              for(ordinal_type d=0; d <dim; ++d)
                dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
            }
        }

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        {
          Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> hexaLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexaLinearBasisValuesAtRefCoords, hexa.getNodeCount(), numRefCoords);
          hexaLinearBasis.getValues(hexaLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(hexaLinearBasisValuesAtDofCoords, numCells, hexa.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( hexaLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              hexaLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtDofCoords(i,k,j);
            }
        }

        DynRankView ConstructWithLabel(dofCoeffs, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);


        basis.getDofCoeffs(dofCoeffs);
        rst::clone(dofCoeffsTmp,dofCoeffs);

        { //orient DofCoeffs
          auto numEdgeDOFs = basis.getDofCount(1,0);
          auto numFaceDOFs = basis.getDofCount(2,0);

          DynRankView ConstructWithLabel(refEdgeTan,  dim);
          DynRankView ConstructWithLabel(refFaceTangents, dim, 2);


          DynRankView ortJacobian("ortJacobian", dim-1, dim-1);
          Basis_HCURL_QUAD_In_FEM<DeviceType,ValueType,ValueType> quadBasis(order);
          DynRankView ConstructWithLabel(quadDofCoeffs, quadBasis.getCardinality(), dim-1);
          quadBasis.getDofCoeffs(quadDofCoeffs);
          for(ordinal_type i=0; i<numCells; ++i) {

            ordinal_type fOrt[numElemFaces];
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              auto refFaceTanU = Kokkos::subview(refFaceTangents, Kokkos::ALL, 0);
              auto refFaceTanV = Kokkos::subview(refFaceTangents, Kokkos::ALL, 1);
              ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, iface, hexa);
              Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, quad, fOrt[iface]);
              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
                auto jdof = quadBasis.getDofOrdinal(dim-1, 0, j);
                for(ordinal_type d=0; d <dim; ++d) {
                  dofCoeffsTmp(i,idof,d) = 0;
                  for(ordinal_type k=0; k <dim-1; ++k)
                    for(ordinal_type l=0; l <dim-1; ++l)
                      dofCoeffsTmp(i,idof,d) += refFaceTangents(d,l)*ortJacobian(l,k)*quadDofCoeffs(jdof,k);
                }
              }
            }

            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {

              ordinal_type eOrt[numElemEdges];
              elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);

              ct::getReferenceEdgeTangent(refEdgeTan, iedge, hexa);
              ValueType edgeTan3d[3] = {};
              for(ordinal_type d=0; d <dim; ++d)
                edgeTan3d[d] = refEdgeTan(d)*((eOrt[iedge] == 0) ? 1 : -1);

              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                for(ordinal_type d=0; d <dim; ++d){
                  dofCoeffsTmp(i,idof,d) = edgeTan3d[d];
                }
              }
            }
          }
        }

        //need to transform dofCoeff to physical space (they transform as tangents)
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hexa);
        ct::setJacobianInv (jacobian_inv, jacobian);
        rst::matvec(dofCoeffsPhys, jacobian, dofCoeffsTmp);

        //check whether dofCoeffs related to faces and edges are proportional to faces' and edges' tangents
        {
          auto numEdgeDOFs = basis.getDofCount(1,0);
          auto numFaceDOFs = basis.getDofCount(2,0);
          const ValueType refEdgeLength = 2.0;
          for(ordinal_type i=0; i<numCells; ++i) {

            for(ordinal_type iface=0; iface<numElemFaces; iface++)
              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                ordinal_type itan = j/(numFaceDOFs/2);     //Convention here is different than for tets!
                auto idof = basis.getDofOrdinal(2, iface, j);
                ValueType tangent[dim];
                bool areDifferent = false;
                for(ordinal_type d=0; d <dim; ++d) {
                  tangent[d] = (itan == 0) ? faceTanU[i][iface][d] : faceTanV[i][iface][d];
                  if(std::abs(dofCoeffsPhys(i,idof,d)*refEdgeLength - tangent[d])>tol)
                    areDifferent = true;
                }
                if(areDifferent) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << "Coefficients of Dof " << idof << " at cell "  << i << ", scaled by the reference edge length, are NOT equivalent to tangent " << itan << " of face " << iface << "\n";
                  *outStream << "Dof Coefficients are: (" << dofCoeffsPhys(i,idof,0) << ", " << dofCoeffsPhys(i,idof,1) << ", " << dofCoeffsPhys(i,idof,2) << ")\n";
                  *outStream << "Face tangent is: (" << tangent[0] << ", " << tangent[1] << ", " << tangent[2] << ")\n";
                }
              }

            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++)
              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                ValueType tangent[dim];
                bool areDifferent = false;
                for(ordinal_type d=0; d <dim; ++d) {
                  tangent[d] = edgeTan[i][iedge][d];
                  if(std::abs(dofCoeffsPhys(i,idof,d)*refEdgeLength - tangent[d])>tol)
                    areDifferent = true;
                }
                if(areDifferent) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << "Coefficients of Dof " << idof << " at cell "  << i << ", scaled by the reference edge length, are NOT equivalent to the tangent of edge " << iedge << "\n";
                  *outStream << "Dof Coefficients are: (" << dofCoeffsPhys(i,idof,0) << ", " << dofCoeffsPhys(i,idof,1) << ", " << dofCoeffsPhys(i,idof,2) << ")\n";
                  *outStream << "Face tangent is: (" << tangent[0] << ", " << tangent[1] << ", " << tangent[2] << ")\n";
                }
              }
          }
        }

        //Testing Kronecker property of basis functions
        {

          DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
          DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
          for(ordinal_type i=0; i<numCells; ++i) {
            auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
            auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
            basis.getValues(outView, inView);
          }

          // modify basis values to account for orientations
          ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
              basisValuesAtDofCoords,
              elemOrts,
              &basis);

          // transform basis values
          fst::HCURLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
              jacobian_inv,
              basisValuesAtDofCoordsOriented);

          for(ordinal_type i=0; i<numCells; ++i) {

            for(ordinal_type k=0; k<basisCardinality; ++k) {
              for(ordinal_type j=0; j<basisCardinality; ++j){
                ValueType dofValue=0;
                for(ordinal_type d=0; d<dim; ++d)
                  dofValue += transformedBasisValuesAtDofCoordsOriented(i,k,j,d) * dofCoeffsPhys(i,j,d);
                if ( k==j && std::abs( dofValue - 1.0 ) > 100*tol ) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function " << k << " of cell " << i << " does not have unit value at its node (" << dofValue <<")\n";
                  *outStream << " Fist tag: " << basis.getDofTag(k)[0] << " " << basis.getDofTag(k)[1]<< " " <<  basis.getDofTag(k)[2] << " Second tag: " <<
                      basis.getDofTag(i)[0] << " " << basis.getDofTag(i)[1]<< " " <<  basis.getDofTag(i)[2] <<")\n";
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


        //check function reproducibility
        FunCurl fun;
        DynRankView ConstructWithLabel(funDofs, numCells, basisCardinality);
        DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords, dim);
        for(ordinal_type i=0; i<numCells; ++i) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<dim; ++k)
              funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
          }
          for(ordinal_type j=0; j<basisCardinality; ++j)
            for(ordinal_type k=0; k<dim; ++k) {
              funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k) * dofCoeffsPhys(i,j,k);
            }
        }

        //check that fun values are consistent on common edges dofs
        {
          bool areDifferent(false);
          auto numEdgeDOFs = basis.getDofCount(1,0);
          for(std::size_t iEdge=0;iEdge<common_edges.size();iEdge++) {
            for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
              areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j))
                  - funDofs(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
            }
            if(areDifferent) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << "Function DOFs on common edge " << iEdge << " computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
              *outStream << "Function DOFs for Hex 0 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
              *outStream << "\nFunction DOFs for Hex 1 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funDofs(1,basis.getDofOrdinal(1,edgeIndexes[1][iEdge],j));
              *outStream << std::endl;
            }
          }
        }

        //check that fun values are consistent on common face dofs
        {
          bool areDifferent(false);
          auto numFaceDOFs = basis.getDofCount(2,0);
          for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j))
                - funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
          }

          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
            *outStream << "Function DOFs for Hex 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Hex 1 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
            *outStream << std::endl;
          }
        }

        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
        DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);


        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numCells, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_inv, numCells, numRefCoords, dim, dim);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, hexa);
        ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
        fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
            jacobianAtRefCoords_inv,
            basisValuesAtRefCoordsOriented);
        DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords, dim);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                funAtRefCoordsOriented(i,j,d) += funDofs(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);

              error = std::max(std::abs( funAtPhysRefCoords(i,j,d) - funAtRefCoordsOriented(i,j,d)), error);
            }

          if(error>100*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtPhysRefCoords(i,j,0) << "," << funAtPhysRefCoords(i,j,1) << ", " << funAtPhysRefCoords(i,j,2) << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtRefCoordsOriented(i,j,0) << "," << funAtRefCoordsOriented(i,j,1) << ", " << funAtRefCoordsOriented(i,j,2) << ")";
            *outStream << std::endl;
          }
        }

      }
    } while(std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

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

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);
      for (ordinal_type shift=0; shift<6; ++shift) {
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

        shards::CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hexa.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hexa.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];


        //computing edges and tangents

        ValueType faceNormal[numCells][numElemFaces][3];
        ordinal_type faceIndex[numCells];
        {
          faceType face={};
          ValueType faceTanU[3];
          ValueType faceTanV[3];
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' normals
            for (std::size_t is=0; is<hexa.getSideCount(); ++is) {
              for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hexa.getNodeMap(2,is,k)];
             
              //rotate and flip and compute common face id
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;

              //Compute face tangents using global numbering
              for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
                face[k]= reorder[face[k]];
              
              //rotate and flip
              minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}
              
              for(int d=0; d<dim; ++d) {
                faceTanU[d] = vertices[face[1]][d]-vertices[face[0]][d];
                faceTanV[d] = vertices[face[2]][d]-vertices[face[0]][d];
              }
              faceNormal[i][is][0] = faceTanU[1]*faceTanV[2] - faceTanU[2]*faceTanV[1];
              faceNormal[i][is][1] = faceTanU[2]*faceTanV[0] - faceTanU[0]*faceTanV[2];
              faceNormal[i][is][2] = faceTanU[0]*faceTanV[1] - faceTanU[1]*faceTanV[0];
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hexa);

        //compute reference points
        Basis_HGRAD_HEX_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_HEX_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          //Basis_HVOL_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> quadBasis(order-1);
          Basis_HVOL_QUAD_Cn_FEM<Kokkos::DefaultHostExecutionSpace> quadBasis(order-1, POINTTYPE_GAUSS);
          ordinal_type quadBasisCardinality = quadBasis.getCardinality();
          ordinal_type numInternalDofs = quadBasis.getDofCount(dim-1,0);
          DynRankView ConstructWithLabel(quadDofCoords, quadBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(quadInternalDofCoords, numInternalDofs, dim-1);
          quadBasis.getDofCoords(quadDofCoords);
          for(ordinal_type i=0; i<numInternalDofs; ++i)
            for(ordinal_type d=0; d <dim-1; ++d)
              quadInternalDofCoords(i,d) = quadDofCoords(quadBasis.getDofOrdinal(dim-1, 0, i),d);

          DynRankView ConstructWithLabel(quadInternalDofCoordsOriented,  numInternalDofs, dim-1);
          DynRankView ConstructWithLabel(quadDofCoordsOriented3d, numInternalDofs, dim);
          ordinal_type fOrt[numElemFaces];
          for(ordinal_type i=0; i<numCells; ++i) {
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              ordinal_type ort = fOrt[iface];
              iots::mapToModifiedReference(quadInternalDofCoordsOriented,quadInternalDofCoords,quad,ort);
              ct::mapToReferenceSubcell(quadDofCoordsOriented3d, quadInternalDofCoordsOriented, dim-1, iface, hexa);
              for(ordinal_type j=0; j<numInternalDofs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iface, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = quadDofCoordsOriented3d(j,d);
              }
            }
          }

          auto numElemDOFs = basis.getDofCount(3,0);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type j=0; j<numElemDOFs; ++j) {
              auto idof = basis.getDofOrdinal(3, 0, j);
              for(ordinal_type d=0; d <dim; ++d)
                dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
            }
        }

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        {
          Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> hexaLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexaLinearBasisValuesAtRefCoords, hexa.getNodeCount(), numRefCoords);
          hexaLinearBasis.getValues(hexaLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(hexaLinearBasisValuesAtDofCoords, numCells, hexa.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( hexaLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              hexaLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtDofCoords(i,k,j);
            }
        }

        DynRankView ConstructWithLabel(dofCoeffs, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsOriented, numCells, basisCardinality, dim);

        basis.getDofCoeffs(dofCoeffs);
        rst::clone(dofCoeffsOriented,dofCoeffs);

        { //orient DofCoeffs
          auto numFaceDOFs = basis.getDofCount(2,0);

          DynRankView ConstructWithLabel(refFaceNormal,  dim);
          DynRankView ortJacobian("ortJacobian", dim-1, dim-1);

          for(ordinal_type i=0; i<numCells; ++i) {
            ordinal_type fOrt[numElemFaces];
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, quad, fOrt[iface]);
              auto ortJacobianDet = ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);

              ct::getReferenceFaceNormal(refFaceNormal, iface, hexa);

              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoeffsOriented(i,idof,d) = refFaceNormal(d)*ortJacobianDet;
              }
            }
          }
        }

        //need to transform dofCoeff to physical space (they transform as normals)
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_invT, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, hexa);
        ct::setJacobianInv (jacobian_inv, jacobian);
        ct::setJacobianDet (jacobian_det, jacobian);
        rst::transpose(jacobian_invT,jacobian_inv);

        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);
        rst::matvec(dofCoeffsTmp, jacobian_invT, dofCoeffsOriented);
        at::scalarMultiplyDataData(dofCoeffsPhys, jacobian_det, dofCoeffsTmp);

        //dofCoeffsPhys do not account for orientation. We overwrite dofCoeffs of tangents and edges.
        {

          const ValueType refQuadArea = 4.0;
          auto numFaceDOFs = basis.getDofCount(2,0);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type iface=0; iface<numElemFaces; iface++)
              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
              ValueType normal[dim];
              bool areDifferent = false;
              for(ordinal_type d=0; d <dim; ++d) {
                normal[d] = faceNormal[i][iface][d];
                if(std::abs(dofCoeffsPhys(i,idof,d)*refQuadArea- normal[d])>tol)
                  areDifferent = true;
              }
              if(areDifferent) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Coefficients of Dof " << idof << " at cell "  << i << ", scaled by the area of the reference Quad, are NOT equivalent to the normal of face " << iface << "\n";
                *outStream << "Dof Coefficients are: (" << dofCoeffsPhys(i,idof,0) << ", " << dofCoeffsPhys(i,idof,1) << ", " << dofCoeffsPhys(i,idof,2) << ")\n";
                *outStream << "Face normal is: (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")\n";
              }
            }
          }
        }

        //Testing Kronecker property of basis functions
        {
          for(ordinal_type i=0; i<numCells; ++i) {
            DynRankView ConstructWithLabel(basisValuesAtDofCoords, numCells, basisCardinality, basisCardinality, dim);
            DynRankView ConstructWithLabel(basisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtDofCoordsOriented, numCells, basisCardinality, basisCardinality, dim);
            auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
            auto outView =Kokkos::subview( basisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
            basis.getValues(outView, inView);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtDofCoordsOriented,
                basisValuesAtDofCoords,
                elemOrts,
                &basis);

            // transform basis values
            fst::HDIVtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
              jacobian,
              jacobian_det,
              basisValuesAtDofCoordsOriented);

            for(ordinal_type k=0; k<basisCardinality; ++k) {
              for(ordinal_type j=0; j<basisCardinality; ++j){
                ValueType dofValue=0;
                for(ordinal_type d=0; d<dim; ++d)
                  dofValue += transformedBasisValuesAtDofCoordsOriented(i,k,j,d) * dofCoeffsPhys(i,j,d);
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

        //check function reproducibility
        FunDiv fun;
        DynRankView ConstructWithLabel(funDofs, numCells, basisCardinality);
        DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords, dim);
        for(ordinal_type i=0; i<numCells; ++i) {
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            for(ordinal_type k=0; k<dim; ++k)
              funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
          }
          for(ordinal_type j=0; j<basisCardinality; ++j)
            for(ordinal_type k=0; k<dim; ++k) {
              funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), physDofCoords(i,j,2), k) * dofCoeffsPhys(i,j,k);
            }
        }

        //check that fun values are consistent on common face dofs
        {
          bool areDifferent(false);
          auto numFaceDOFs = basis.getDofCount(2,0);
          for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j))
                - funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
          }

          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
            *outStream << "Function DOFs for Hex 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Hex 1 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
            *outStream << std::endl;
          }
        }

        //check that fun values are consistent on common face dofs
        {
          bool areDifferent(false);
          auto numFaceDOFs = basis.getDofCount(2,0);
          for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j))
                - funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
          }

          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
            *outStream << "Function DOFs for Hex 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Hex 1 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(1,basis.getDofOrdinal(2,faceIndex[1],j))<< " | (" << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),0) << "," << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),1) << ", " << physDofCoords(1,basis.getDofOrdinal(2,faceIndex[1],j),2) << ") ||";
            *outStream << std::endl;
          }
        }

        //check that fun values at reference points coincide with those computed using basis functions
        DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
        DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
        DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);

        DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
        basis.getValues(basisValuesAtRefCoords, refPoints);
        rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

        // modify basis values to account for orientations
        ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
            basisValuesAtRefCoordsCells,
            elemOrts,
            &basis);

        // transform basis values
        DynRankView ConstructWithLabel(jacobianAtRefCoords, numCells, numRefCoords, dim, dim);
        DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numCells, numRefCoords);
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, hexa);
        ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
        fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
            jacobianAtRefCoords,
            jacobianAtRefCoords_det,
            basisValuesAtRefCoordsOriented);
        DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords, dim);
        for(ordinal_type i=0; i<numCells; ++i) {
          ValueType error=0;
          for(ordinal_type j=0; j<numRefCoords; ++j)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                funAtRefCoordsOriented(i,j,d) += funDofs(i,k)*transformedBasisValuesAtRefCoordsOriented(i,k,j,d);

              error = std::max(std::abs( funAtPhysRefCoords(i,j,d) - funAtRefCoordsOriented(i,j,d)), error);
            }

          if(error>100*tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function values at reference points differ from those computed using basis functions of Hex " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtPhysRefCoords(i,j,0) << "," << funAtPhysRefCoords(i,j,1) << ", " << funAtPhysRefCoords(i,j,2) << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtRefCoordsOriented(i,j,0) << "," << funAtRefCoordsOriented(i,j,1) << ", " << funAtRefCoordsOriented(i,j,2) << ")";
            *outStream << std::endl;
          }
        }

      }
    } while(std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

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

