// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking orientation tools for Tetrahedral elements
 
    The test considers two tetrahedra in the physical space sharing a common face. 
    In order to test significant configurations, we consider 4 mappings of the reference tetrahedron 
    to the first (physical) tetrahedron, so that the common face is mapped from all the 4 faces 
    of the reference tetrahedron.
    Then, for each of the mappings, the global ids of the vertices of the common faces are permuted.
    This gives a total of 24 combinations
    
    The test considers HGRAD, HCURL and HDIV elements and for each of them checks the following:
    1. that the physical (oriented) dof coefficients related to faces and edges are proportional
       to the tangents (for HCURL) of the physical faces/edges or to the normals (for HDIV) of physical faces.
    2. that the Kronecker property of the oriented basis evaluated at the physiscal dof coordinates holds.
    3. that given a function, the dofs of the function located at the common faces/edges
       are the same when computed one the first and second tetrahedron.
    4. that a function, belonging to the considered H-space, is exactly at given points using the
       oriented basis functions and the degrees of freedom of the function.
 
    

    \author Created by Mauro Perego
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
#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"
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
int OrientationTet(const bool verbose) {

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
  typedef std::array<ordinal_type,3> faceType;
  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef Impl::OrientationTools iots;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;
  typedef ArrayTools<DeviceType> at;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 4;
  constexpr ordinal_type numElemEdges = 6;
  constexpr ordinal_type numElemFaces = 4;
  constexpr ordinal_type numTotalVertexes = 5;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
  ordinal_type tets_orig[numCells][numElemVertexes] = {{0,1,2,3},{1,2,3,4}};  //common face is {1,2,3}
  ordinal_type tets_rotated[numCells][numElemVertexes];
  faceType common_face = {{1,2,3}};
  std::set<edgeType> common_edges;
  common_edges.insert(edgeType({{1,2}})); common_edges.insert(edgeType({{1,3}})); common_edges.insert(edgeType({{2,3}}));


  *outStream
 << "===============================================================================\n"
 << "|                                                                             |\n"
 << "|                 Test 1 (Orientation - HGRAD)                                |\n"
 << "|                                                                             |\n"
 << "===============================================================================\n";


 try {

   const ordinal_type order = 4;
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

       //computing common face and edges
       ordinal_type faceIndex[numCells];
       ordinal_type edgeIndexes[numCells][3];

       {
         faceType face={};
         edgeType edge={};
         //bool faceOrientation[numCells][4];
         for(ordinal_type i=0; i<numCells; ++i) {
           //compute faces' tangents
           for (std::size_t is=0; is<tet.getSideCount(); ++is) {
             for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
               face[k]= tets_rotated[i][tet.getNodeMap(2,is,k)];
             std::sort(face.begin(),face.end());
             if(face == common_face) faceIndex[i]=is;
           }
           //compute edges' tangents
           for (std::size_t ie=0; ie<tet.getEdgeCount(); ++ie) {
             for (std::size_t k=0; k<tet.getNodeCount(1,ie); ++k)
               edge[k]= tets_rotated[i][tet.getNodeMap(1,ie,k)];
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
       Basis_HGRAD_TET_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
       ordinal_type numRefCoords = warpBasis.getCardinality();
       DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
       warpBasis.getDofCoords(refPoints);

       // compute orientations for cells (one time computation)
       DynRankViewInt elemNodes(&tets[0][0], numCells, numElemVertexes);
       Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
       ots::getOrientation(elemOrts, elemNodes, tet);

       Basis_HGRAD_TET_Cn_FEM<DeviceType,ValueType,ValueType> basis(order);
       ordinal_type basisCardinality = basis.getCardinality();

       //compute DofCoords Oriented
       DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
       DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
       basis.getDofCoords(dofCoords);
       {
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type ivertex=0; ivertex<4; ++ivertex) {
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
             ct::mapToReferenceSubcell(lineDofCoordsOriented3d, lineDofCoordsOriented, dim-2, iedge, tet);

             for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
               auto idof = basis.getDofOrdinal(1, iedge, j);
               auto linedof = lineBasis.getDofOrdinal(1, 0, j);

               for(ordinal_type d=0; d <dim; ++d)
                 dofCoordsOriented(i,idof,d) = lineDofCoordsOriented3d(linedof,d);
             }
           }
         }
        
         Basis_HGRAD_TRI_Cn_FEM<DeviceType,ValueType,ValueType> triBasis(order);
         ordinal_type triBasisCardinality = triBasis.getCardinality();
         ordinal_type numInternalDofs = triBasis.getDofCount(dim-1,0);
         DynRankView ConstructWithLabel(triDofCoords, triBasisCardinality, dim-1);
         DynRankView ConstructWithLabel(triInternalDofCoords, numInternalDofs, dim-1);
         triBasis.getDofCoords(triDofCoords);
         for(ordinal_type i=0; i<numInternalDofs; ++i)
           for(ordinal_type d=0; d <dim-1; ++d)
             triInternalDofCoords(i,d) = triDofCoords(triBasis.getDofOrdinal(dim-1, 0, i),d);

         DynRankView ConstructWithLabel(triInternalDofCoordsOriented,  numInternalDofs, dim-1);
         DynRankView ConstructWithLabel(triDofCoordsOriented3d, numInternalDofs, dim);
         ordinal_type fOrt[numElemFaces];
         for(ordinal_type i=0; i<numCells; ++i) {
           elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
           for(ordinal_type iface=0; iface<numElemFaces; iface++) {
             ordinal_type ort = fOrt[iface];
             iots::mapToModifiedReference(triInternalDofCoordsOriented,triInternalDofCoords,tri,ort);
             ct::mapToReferenceSubcell(triDofCoordsOriented3d, triInternalDofCoordsOriented, dim-1, iface, tet);

             for(ordinal_type j=0; j<numInternalDofs; ++j) {
               auto idof = basis.getDofOrdinal(2, iface, j);
               for(ordinal_type d=0; d <dim; ++d)
                 dofCoordsOriented(i,idof,d) = triDofCoordsOriented3d(j,d);
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
         Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
         DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
         tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
         DynRankView ConstructWithLabel(tetLinearBasisValuesAtDofCoords, numCells, tet.getNodeCount(), basisCardinality);
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type d=0; d<dim; ++d) {
             for(ordinal_type j=0; j<numRefCoords; ++j)
               for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                 physRefCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtRefCoords(k,j);

             auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
             auto outView =Kokkos::subview( tetLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
             tetLinearBasis.getValues(outView, inView);

             for(ordinal_type j=0; j<basisCardinality; ++j)
               for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                 physDofCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtDofCoords(i,k,j);
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

       //check function reproducbility
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
           *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
           *outStream << "Function DOFs for Tet 0 are:";
           for(ordinal_type j=0;j<numFaceDOFs;j++)
             *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
           *outStream << "\nFunction DOFs for Tet 1 are:";
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
             *outStream << "Function DOFs on common edge " << iEdge << " computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
             *outStream << "Function DOFs for Tet 0 are:";
             for(ordinal_type j=0;j<numEdgeDOFs;j++)
               *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
             *outStream << "\nFunction DOFs for Tet 1 are:";
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
           *outStream << "Function values at reference points differ from those computed using basis functions of Tet " << i << "\n";
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
   } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

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



        //computing edges and faces
        ordinal_type faceIndex[numCells];
        ordinal_type edgeIndexes[numCells][3];
        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute faces' tangents
            for (std::size_t is=0; is<tet.getSideCount(); ++is) {
              for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
                face[k]= tets_rotated[i][tet.getNodeMap(2,is,k)];
              
              //sort face and compute common face
              std::sort(face.begin(),face.end());

              if(face == common_face) faceIndex[i]=is;

              //compute tangents using global ordering
              for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
                face[k]= reorder[face[k]];

              std::sort(face.begin(),face.end());
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<tet.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<tet.getNodeCount(1,ie); ++k)
                edge[k]= tets_rotated[i][tet.getNodeMap(1,ie,k)];

              //compute common edge        
              std::sort(edge.begin(),edge.end());
              auto it=common_edges.find(edge);
              if(it !=common_edges.end()){
                auto edge_lid = std::distance(common_edges.begin(),it);
                edgeIndexes[i][edge_lid]=ie;
              }
           
              //compute edge tangent using global numbering
              for (std::size_t k=0; k<tet.getNodeCount(1,ie); ++k)
                edge[k] = reorder[edge[k]];
              std::sort(edge.begin(),edge.end());
            }
          }
        }

        //compute reference points
        Basis_HGRAD_TET_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&tets[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);


        Basis_HCURL_TET_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          Basis_HVOL_LINE_Cn_FEM<DeviceType,ValueType,ValueType> lineBasis(order-1);
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
              ct::mapToReferenceSubcell(lineDofCoordsOriented3d, lineDofCoordsOriented, dim-2, iedge, tet);

              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                auto linedof = lineBasis.getDofOrdinal(1, 0, j);

                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = lineDofCoordsOriented3d(linedof,d);
              }
            }
          }

          Basis_HCURL_TRI_In_FEM<DeviceType,ValueType,ValueType> triBasis(order);
          ordinal_type triBasisCardinality = triBasis.getCardinality();
          ordinal_type numInternalDofs = triBasis.getDofCount(dim-1,0);
          DynRankView ConstructWithLabel(triDofCoords, triBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(triInternalDofCoords, numInternalDofs, dim-1);
          triBasis.getDofCoords(triDofCoords);
          for(ordinal_type i=0; i<numInternalDofs; ++i)
            for(ordinal_type d=0; d <dim-1; ++d)
              triInternalDofCoords(i,d) = triDofCoords(triBasis.getDofOrdinal(dim-1, 0, i),d);

          DynRankView ConstructWithLabel(triInternalDofCoordsOriented,  numInternalDofs, dim-1);
          DynRankView ConstructWithLabel(triDofCoordsOriented3d, numInternalDofs, dim);
          ordinal_type fOrt[numElemFaces];
          for(ordinal_type i=0; i<numCells; ++i) {
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              ordinal_type ort = fOrt[iface];
              iots::mapToModifiedReference(triInternalDofCoordsOriented,triInternalDofCoords,tri,ort);
              ct::mapToReferenceSubcell(triDofCoordsOriented3d, triInternalDofCoordsOriented, dim-1, iface, tet);

              for(ordinal_type j=0; j<numInternalDofs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = triDofCoordsOriented3d(j,d);
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
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
          tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(tetLinearBasisValuesAtDofCoords, numCells, tet.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( tetLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              tetLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtDofCoords(i,k,j);
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
          Basis_HCURL_TRI_In_FEM<DeviceType,ValueType,ValueType> triBasis(order);
          DynRankView ConstructWithLabel(triDofCoeffs, triBasis.getCardinality(), dim-1);
          triBasis.getDofCoeffs(triDofCoeffs);
          for(ordinal_type i=0; i<numCells; ++i) {
            ordinal_type fOrt[numElemFaces];
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              auto refFaceTanU = Kokkos::subview(refFaceTangents, Kokkos::ALL, 0);
              auto refFaceTanV = Kokkos::subview(refFaceTangents, Kokkos::ALL, 1);
              ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, iface, tet);
              Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, tri, fOrt[iface]);
              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
                auto jdof = triBasis.getDofOrdinal(dim-1, 0, j);
                for(ordinal_type d=0; d <dim; ++d) {
                  dofCoeffsTmp(i,idof,d) = 0;
                  for(ordinal_type k=0; k <dim-1; ++k)
                    for(ordinal_type l=0; l <dim-1; ++l)
                      dofCoeffsTmp(i,idof,d) += refFaceTangents(d,l)*ortJacobian(l,k)*triDofCoeffs(jdof,k);
                }
              }
            }

            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {

              ordinal_type eOrt[numElemEdges];
              elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);

              ct::getReferenceEdgeTangent(refEdgeTan, iedge, tet);
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
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, tet);
        ct::setJacobianInv (jacobian_inv, jacobian);
        rst::matvec(dofCoeffsPhys, jacobian, dofCoeffsTmp);

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
            fst::HCURLtransformVALUE(transformedBasisValuesAtDofCoordsOriented,
                jacobian_inv,
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
              *outStream << "Function DOFs on common edge " << iEdge << " computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
              *outStream << "Function DOFs for Tet 0 are:";
              for(ordinal_type j=0;j<numEdgeDOFs;j++)
                *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndexes[0][iEdge],j));
              *outStream << "\nFunction DOFs for Tet 1 are:";
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
            *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
            *outStream << "Function DOFs for Tet 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Tet 1 are:";
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
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
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
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

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
            for (std::size_t is=0; is<tet.getSideCount(); ++is) {
              for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
                face[k]= tets_rotated[i][tet.getNodeMap(2,is,k)];
              std::sort(face.begin(),face.end());
              if(face == common_face) faceIndex[i]=is;

              //Compute face tangents using global numbering
              for (std::size_t k=0; k<tet.getNodeCount(2,is); ++k)
                face[k]= reorder[face[k]];
              std::sort(face.begin(),face.end());
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
        DynRankViewInt elemNodes(&tets[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, tet);

        //compute reference points
        Basis_HGRAD_TET_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_TET_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          Basis_HVOL_TRI_Cn_FEM<DeviceType,ValueType,ValueType> triBasis(order-1);
          ordinal_type triBasisCardinality = triBasis.getCardinality();
          ordinal_type numInternalDofs = triBasis.getDofCount(dim-1,0);
          DynRankView ConstructWithLabel(triDofCoords, triBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(triInternalDofCoords, numInternalDofs, dim-1);
          triBasis.getDofCoords(triDofCoords);
       
          for(ordinal_type i=0; i<numInternalDofs; ++i)
            for(ordinal_type d=0; d <dim-1; ++d)
              triInternalDofCoords(i,d) = triDofCoords(triBasis.getDofOrdinal(dim-1, 0, i),d);

          DynRankView ConstructWithLabel(triInternalDofCoordsOriented,  numInternalDofs, dim-1);
          DynRankView ConstructWithLabel(triDofCoordsOriented3d, numInternalDofs, dim);
          ordinal_type fOrt[numElemFaces];
          for(ordinal_type i=0; i<numCells; ++i) {
            elemOrts(i).getFaceOrientation(fOrt, numElemFaces);
            for(ordinal_type iface=0; iface<numElemFaces; iface++) {
              ordinal_type ort = fOrt[iface];
              iots::mapToModifiedReference(triInternalDofCoordsOriented,triInternalDofCoords,tri,ort);
              ct::mapToReferenceSubcell(triDofCoordsOriented3d, triInternalDofCoordsOriented, dim-1, iface, tet);
              for(ordinal_type j=0; j<numInternalDofs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iface, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = triDofCoordsOriented3d(j,d);
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
          Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType> tetLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(tetLinearBasisValuesAtRefCoords, tet.getNodeCount(), numRefCoords);
          tetLinearBasis.getValues(tetLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(tetLinearBasisValuesAtDofCoords, numCells, tet.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( tetLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              tetLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<tet.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[tets[i][k]][d]*tetLinearBasisValuesAtDofCoords(i,k,j);
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
              Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, tri, fOrt[iface]);
              auto ortJacobianDet = ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);

              ct::getReferenceFaceNormal(refFaceNormal, iface, tet);

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
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, tet);
        ct::setJacobianInv (jacobian_inv, jacobian);
        ct::setJacobianDet (jacobian_det, jacobian);
        rst::transpose(jacobian_invT,jacobian_inv);

        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);
        rst::matvec(dofCoeffsTmp, jacobian_invT, dofCoeffsOriented);
        at::scalarMultiplyDataData(dofCoeffsPhys, jacobian_det, dofCoeffsTmp);

        //dofCoeffsPhys do not account for orientation. We overwrite dofCoeffs of tangents and edges.
        {
          auto numFaceDOFs = basis.getDofCount(2,0);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type iface=0; iface<numElemFaces; iface++)
              for(ordinal_type j=0; j<numFaceDOFs; ++j) {
                auto idof = basis.getDofOrdinal(2, iface, j);
              ValueType normal[dim];
              bool areDifferent = false;
              for(ordinal_type d=0; d <dim; ++d) {
                normal[d] = faceNormal[i][iface][d];
                if(std::abs(dofCoeffsPhys(i,idof,d)- normal[d])>tol)
                  areDifferent = true;
              }
              if(areDifferent) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Coefficients of Dof " << idof << " at cell "  << i << " are NOT equivalent to the normal of face " << iface << "\n";
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
            *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
            *outStream << "Function DOFs for Tet 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Tet 1 are:";
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
            *outStream << "Function DOFs on common face computed using Tet 0 basis functions are not consistent with those computed using Tet 1\n";
            *outStream << "Function DOFs for Tet 0 are:";
            for(ordinal_type j=0;j<numFaceDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(2,faceIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),1) << ", " << physDofCoords(0,basis.getDofOrdinal(2,faceIndex[0],j),2) << ") ||";
            *outStream << "\nFunction DOFs for Tet 1 are:";
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
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, tet);
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
    } while(std::next_permutation(&reorder[0]+1, &reorder[0]+4)); //reorder vertices of common face

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

