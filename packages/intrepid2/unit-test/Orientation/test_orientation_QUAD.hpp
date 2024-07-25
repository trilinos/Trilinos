// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking orientation tools for quadrialateral elements
 
    The test considers two quadrilaterals in the physical space sharing a common edge. 
    In order to test significant configurations, we consider 4 mappings of the reference quadrilateral
    to the first (physical) quadrilateral, so that the common edge is mapped from all the 4 edges
    of the reference quadrilateral.
    Then, for each of the mappings, the global ids of the vertices of the common edge are permuted.
    This gives a total of 8 combinations
    
    The test considers HGRAD, HCURL and HDIV elements and for each of them checks the following:
    1. that the physical (oriented) dof coefficients related to edges are proportional to the
       tangents (for HCURL) of the edges or to the normals (for HDIV) of the edges.
    2. that the Kronecker property of the oriented basis evaluated at the physiscal dof coordinates holds.
    3. that given a function, the dofs of the function located at the common edge are the same
       when computed on the first and second quadrilateral.
    4. that a function, belinging to the considered H-space, is exactly reproduced at given points
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
int OrientationQuad(const bool verbose) {

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
    operator()(const ValueType& x, const ValueType& y) {
      return (x+1)*(y-2);//*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a = 2*x*y+x*x;
      ValueType f0 = 5+y+x*x;
      ValueType f1 = -7-2*x+y*y;
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
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a2 = y-7+y*y;
      ValueType f0 = 2+x+x*y;
      ValueType f1 = 3-3*x;
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 - a2*y;
      case 1:
        return f1 + a2*x;
      default:
        return 0;
      }
    }
  };

  typedef std::array<ordinal_type,2> edgeType;
  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef Impl::OrientationTools iots;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;
  typedef ArrayTools<DeviceType> at;

  constexpr ordinal_type dim = 2;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 4;
  constexpr ordinal_type numElemEdges = 4;
  constexpr ordinal_type numTotalVertexes = 6;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1},{1,-1},{1,1},{-1,1}, {-1,2},{1,2}};
  ordinal_type quads_orig[numCells][numElemVertexes] = {{0,1,2,3},{2,3,4,5}};  
  edgeType common_edge = {{2,3}};
  ordinal_type quads_rotated[numCells][numElemVertexes];

  *outStream
 << "===============================================================================\n"
 << "|                                                                             |\n"
 << "|                 Test 1 (Orientation - HGRAD)                                |\n"
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



       //computing common and edge
       ordinal_type edgeIndex[numCells];

       {
         edgeType edge={};
         for(ordinal_type i=0; i<numCells; ++i) {
           //compute common edge
           for (std::size_t ie=0; ie<quad.getSideCount(); ++ie) {
             for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
               edge[k]= quads_rotated[i][quad.getNodeMap(1,ie,k)];
             
             if(edge == common_edge) edgeIndex[i]=ie;
           }
         }
       }

       //compute reference points
       Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
       ordinal_type numRefCoords = warpBasis.getCardinality();
       DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
       warpBasis.getDofCoords(refPoints);

       // compute orientations for cells (one time computation)
       DynRankViewInt elemNodes(&quads[0][0], numCells, numElemVertexes);
       Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
       ots::getOrientation(elemOrts, elemNodes, quad);

       Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> basis(order);
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
         DynRankView ConstructWithLabel(lineDofCoords, lineBasisCardinality, dim-1);
         lineBasis.getDofCoords(lineDofCoords);
         DynRankView ConstructWithLabel(lineDofCoordsOriented,  lineBasisCardinality, dim-1);
         DynRankView ConstructWithLabel(lineDofCoordsOriented2d,  lineBasisCardinality, dim);
         auto numEdgeDOFs = basis.getDofCount(1,0);
         for(ordinal_type i=0; i<numCells; ++i) {
           ordinal_type eOrt[numElemEdges];
           elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);
           for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
             iots::mapToModifiedReference(lineDofCoordsOriented,lineDofCoords,line,eOrt[iedge]);
             ct::mapToReferenceSubcell(lineDofCoordsOriented2d, lineDofCoordsOriented, dim-1, iedge, quad);

             for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
               auto idof = basis.getDofOrdinal(1, iedge, j);
               auto linedof = lineBasis.getDofOrdinal(1, 0, j);

               for(ordinal_type d=0; d <dim; ++d)
                 dofCoordsOriented(i,idof,d) = lineDofCoordsOriented2d(linedof,d);
             }
           }
         }

         auto numElemDOFs = basis.getDofCount(dim,0);
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type j=0; j<numElemDOFs; ++j) {
             auto idof = basis.getDofOrdinal(dim, 0, j);
             for(ordinal_type d=0; d <dim; ++d)
               dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
           }
       }

       //Compute physical Dof Coordinates and Reference coordinates
       DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
       DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
       {
         Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
         DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
         quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
         DynRankView ConstructWithLabel(quadLinearBasisValuesAtDofCoords, numCells, quad.getNodeCount(), basisCardinality);
         for(ordinal_type i=0; i<numCells; ++i)
           for(ordinal_type d=0; d<dim; ++d) {
             for(ordinal_type j=0; j<numRefCoords; ++j)
               for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                 physRefCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtRefCoords(k,j);

             auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
             auto outView =Kokkos::subview( quadLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
             quadLinearBasis.getValues(outView, inView);

             for(ordinal_type j=0; j<basisCardinality; ++j)
               for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                 physDofCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtDofCoords(i,k,j);
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
             funAtPhysRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
         for(ordinal_type j=0; j<basisCardinality; ++j)
             funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1)) * dofCoeffsPhys(i,j);
       }


       //check that fun values are consistent on common edges dofs
       {
         bool areDifferent(false);
         auto numEdgeDOFs = basis.getDofCount(1,0);
         
         for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
           areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(1,edgeIndex[0],j))
               - funDofs(1,basis.getDofOrdinal(1,edgeIndex[1],j))) > 10*tol;
         }
         if(areDifferent) {
           errorFlag++;
           *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
           *outStream << "Function DOFs on common edge computed using Quad 0 basis functions are not consistent with those computed using Quad 1\n";
           *outStream << "Function DOFs for Quad 0 are:";
           for(ordinal_type j=0;j<numEdgeDOFs;j++)
             *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndex[0],j));
           *outStream << "\nFunction DOFs for Quad 1 are:";
           for(ordinal_type j=0;j<numEdgeDOFs;j++)
             *outStream << " " << funDofs(1,basis.getDofOrdinal(1,edgeIndex[1],j));
           *outStream << std::endl;
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
           *outStream << "Function values at reference points differ from those computed using basis functions of Quad " << i << "\n";
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
   } while(std::next_permutation(&reorder[0]+2, &reorder[0]+4)); //reorder vertices of common edge

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



        //computing edges and tangents
        ValueType edgeTan[numCells][numElemEdges][dim];
        ordinal_type edgeIndex[numCells];

        {
          edgeType edge={};
          //bool edgeOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {

            //compute edges' tangents
            for (std::size_t ie=0; ie<quad.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
                edge[k]= quads_rotated[i][quad.getNodeMap(1,ie,k)];

              //compute common edge
              if(edge == common_edge)
                edgeIndex[i]=ie;
           
              //compute edge tangent using global numbering
              for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
                edge[k] = reorder[edge[k]];
              std::sort(edge.begin(),edge.end());
              for(int d=0; d<dim; ++d)
                edgeTan[i][ie][d] = vertices[edge[1]][d]-vertices[edge[0]][d];
            }
          }
        }

        //compute reference points
        Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);


        Basis_HCURL_QUAD_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          Basis_HVOL_LINE_Cn_FEM<DeviceType,ValueType,ValueType> lineBasis(order-1,POINTTYPE_GAUSS);
          ordinal_type lineBasisCardinality = lineBasis.getCardinality();
          DynRankView ConstructWithLabel(lineDofCoords, lineBasisCardinality, dim-1);
          lineBasis.getDofCoords(lineDofCoords);
          DynRankView ConstructWithLabel(lineDofCoordsOriented,  lineBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(lineDofCoordsOriented2d,  lineBasisCardinality, dim);
          auto numEdgeDOFs = basis.getDofCount(1,0);
          for(ordinal_type i=0; i<numCells; ++i) {
            ordinal_type eOrt[numElemEdges];
            elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);
            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
              iots::mapToModifiedReference(lineDofCoordsOriented,lineDofCoords,line,eOrt[iedge]);
              ct::mapToReferenceSubcell(lineDofCoordsOriented2d, lineDofCoordsOriented, dim-1, iedge, quad);

              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                auto linedof = lineBasis.getDofOrdinal(1, 0, j);

                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = lineDofCoordsOriented2d(linedof,d);
              }
            }
          }

          auto numElemDOFs = basis.getDofCount(dim,0);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type j=0; j<numElemDOFs; ++j) {
              auto idof = basis.getDofOrdinal(dim, 0, j);
              for(ordinal_type d=0; d <dim; ++d)
                dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
            }
        }

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
          quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtDofCoords, numCells, quad.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( quadLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              quadLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtDofCoords(i,k,j);
            }
        }

        DynRankView ConstructWithLabel(dofCoeffs, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);


        basis.getDofCoeffs(dofCoeffs);
        rst::clone(dofCoeffsTmp,dofCoeffs);

        { //orient DofCoeffs
          auto numEdgeDOFs = basis.getDofCount(1,0);
          DynRankView ConstructWithLabel(refEdgeTan,  dim);

          for(ordinal_type i=0; i<numCells; ++i) {

            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {

              ordinal_type eOrt[numElemEdges];
              elemOrts(i).getEdgeOrientation(eOrt, numElemEdges);

              ct::getReferenceEdgeTangent(refEdgeTan, iedge, quad);
              ValueType edgeTan2d[dim] = {};
              for(ordinal_type d=0; d <dim; ++d)
                edgeTan2d[d] = refEdgeTan(d)*((eOrt[iedge] == 0) ? 1 : -1);

              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(1, iedge, j);
                for(ordinal_type d=0; d <dim; ++d){
                  dofCoeffsTmp(i,idof,d) = edgeTan2d[d];
                }
              }
            }
          }
        }

        //need to transform dofCoeff to physical space (they transform as tangents)
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, quad);
        ct::setJacobianInv (jacobian_inv, jacobian);
        rst::matvec(dofCoeffsPhys, jacobian, dofCoeffsTmp);

        //check whether dofCoeffs related to edges are proportional to edges' tangents
        {
          const ValueType refEdgeLength = 2.0;
          auto numEdgeDOFs = basis.getDofCount(1,0);
          for(ordinal_type i=0; i<numCells; ++i) {

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
                  *outStream << "Coefficients of Dof " << idof << " at cell "  << i << ", scaled by the length of the reference edge, are NOT equivalent to the tangent of edge " << iedge << "\n";
                  *outStream << "Dof Coefficients are: (" << dofCoeffsPhys(i,idof,0) << ", " << dofCoeffsPhys(i,idof,1)  << ")\n";
                  *outStream << "Edge tangent is: (" << tangent[0] << ", " << tangent[1] << ")\n";
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
              funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
          }
          for(ordinal_type j=0; j<basisCardinality; ++j)
            for(ordinal_type k=0; k<dim; ++k) {
              funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), k) * dofCoeffsPhys(i,j,k);
            }
        }

        //check that fun values are consistent on common edges dofs
        {
          bool areDifferent(false);
          auto numEdgeDOFs = basis.getDofCount(1,0);
         
          for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(1,edgeIndex[0],j))
                - funDofs(1,basis.getDofOrdinal(1,edgeIndex[1],j))) > 10*tol;
          }
          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on the common edge computed using Quad 0 basis functions are not consistent with those computed using Quad 1\n";
            *outStream << "Function DOFs for Quad 0 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(1,edgeIndex[0],j));
            *outStream << "\nFunction DOFs for Quad 1 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funDofs(1,basis.getDofOrdinal(1,edgeIndex[1],j));
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
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
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
            *outStream << "Function values at reference points differ from those computed using basis functions of Quad " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtPhysRefCoords(i,j,0) << "," << funAtPhysRefCoords(i,j,1) << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtRefCoordsOriented(i,j,0) << "," << funAtRefCoordsOriented(i,j,1) << ")";
            *outStream << std::endl;
          }
        }

      }
    } while(std::next_permutation(&reorder[0]+2, &reorder[0]+4)); //reorder vertices of common edge

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


        //computing edges and tangents

        ValueType edgeNormal[numCells][numElemEdges][dim];
        ordinal_type edgeIndex[numCells];
        {
          edgeType edge={};
          ValueType edgeTan[dim];
          for(ordinal_type i=0; i<numCells; ++i) {
            //compute edges' normals
            for (std::size_t ie=0; ie<quad.getSideCount(); ++ie) {
              for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
                edge[k]= quads_rotated[i][quad.getNodeMap(1,ie,k)];

              if(edge == common_edge) edgeIndex[i]=ie;

              //Compute edge tangents using global numbering
              for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
                edge[k]= reorder[edge[k]];

              std::sort(edge.begin(), edge.end());
              
              for(int d=0; d<dim; ++d)
                edgeTan[d] = vertices[edge[1]][d]-vertices[edge[0]][d];
              
              edgeNormal[i][ie][0] = edgeTan[1];
              edgeNormal[i][ie][1] = -edgeTan[0];
            }
          }
        }

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);

        //compute reference points
        Basis_HGRAD_QUAD_Cn_FEM<DeviceType,ValueType,ValueType> warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        Basis_HDIV_QUAD_In_FEM<DeviceType,ValueType,ValueType> basis(order);
        ordinal_type basisCardinality = basis.getCardinality();

        //compute DofCoords Oriented
        DynRankView ConstructWithLabel(dofCoords, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoordsOriented, numCells, basisCardinality, dim);
        basis.getDofCoords(dofCoords);
        {
          Basis_HVOL_LINE_Cn_FEM<Kokkos::DefaultHostExecutionSpace> sideBasis(order-1, POINTTYPE_GAUSS);
          ordinal_type sideBasisCardinality = sideBasis.getCardinality();
          ordinal_type numInternalDofs = sideBasis.getDofCount(dim-1,0);
          DynRankView ConstructWithLabel(sideDofCoords, sideBasisCardinality, dim-1);
          DynRankView ConstructWithLabel(sideInternalDofCoords, numInternalDofs, dim-1);
          sideBasis.getDofCoords(sideDofCoords);
          for(ordinal_type i=0; i<numInternalDofs; ++i)
            for(ordinal_type d=0; d <dim-1; ++d)
              sideInternalDofCoords(i,d) = sideDofCoords(sideBasis.getDofOrdinal(dim-1, 0, i),d);

          DynRankView ConstructWithLabel(sideInternalDofCoordsOriented,  numInternalDofs, dim-1);
          DynRankView ConstructWithLabel(sideDofCoordsOriented2d, numInternalDofs, dim);
          ordinal_type fOrt[numElemEdges];
          for(ordinal_type i=0; i<numCells; ++i) {
            elemOrts(i).getEdgeOrientation(fOrt, numElemEdges);
            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
              ordinal_type ort = fOrt[iedge];
              iots::mapToModifiedReference(sideInternalDofCoordsOriented,sideInternalDofCoords,line,ort);
              ct::mapToReferenceSubcell(sideDofCoordsOriented2d, sideInternalDofCoordsOriented, dim-1, iedge, quad);
              for(ordinal_type j=0; j<numInternalDofs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iedge, j);
                for(ordinal_type d=0; d <dim; ++d)
                  dofCoordsOriented(i,idof,d) = sideDofCoordsOriented2d(j,d);
              }
            }
          }

          auto numElemDOFs = basis.getDofCount(dim,0);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type j=0; j<numElemDOFs; ++j) {
              auto idof = basis.getDofOrdinal(dim, 0, j);
              for(ordinal_type d=0; d <dim; ++d)
                dofCoordsOriented(i,idof,d) = dofCoords(idof,d);
            }
        }

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        DynRankView ConstructWithLabel(physDofCoords, numCells, basisCardinality, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
          quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtDofCoords, numCells, quad.getNodeCount(), basisCardinality);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtRefCoords(k,j);

              auto inView = Kokkos::subview( dofCoordsOriented,i,Kokkos::ALL(),Kokkos::ALL());
              auto outView =Kokkos::subview( quadLinearBasisValuesAtDofCoords,i,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
              quadLinearBasis.getValues(outView, inView);

              for(ordinal_type j=0; j<basisCardinality; ++j)
                for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                  physDofCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtDofCoords(i,k,j);
            }
        }

        DynRankView ConstructWithLabel(dofCoeffs, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsPhys, numCells, basisCardinality, dim);
        DynRankView ConstructWithLabel(dofCoeffsOriented, numCells, basisCardinality, dim);

        basis.getDofCoeffs(dofCoeffs);
        rst::clone(dofCoeffsOriented,dofCoeffs);

        { //orient DofCoeffs
          auto numEdgeDOFs = basis.getDofCount(dim-1,0);

          const ordinal_type c[2] = {  1, -1 };

          DynRankView ConstructWithLabel(refEdgeNormal,  dim);

          for(ordinal_type i=0; i<numCells; ++i) {
            ordinal_type fOrt[numElemEdges];
            elemOrts(i).getEdgeOrientation(fOrt, numElemEdges);
            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++) {
              ct::getReferenceSideNormal(refEdgeNormal, iedge, quad);
              ValueType edgeNormal2d[dim];
              for(ordinal_type d=0; d <dim; ++d)
                edgeNormal2d[d] = refEdgeNormal(d)*c[fOrt[iedge]];
              
              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iedge, j);

                for(ordinal_type d=0; d <dim; ++d)
                  dofCoeffsOriented(i,idof,d) = edgeNormal2d[d];
              }
            }
          }
        }

        //need to transform dofCoeff to physical space (they transform as normals)
        DynRankView ConstructWithLabel(jacobian, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_inv, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_invT, numCells, basisCardinality, dim, dim);
        DynRankView ConstructWithLabel(jacobian_det, numCells, basisCardinality);
        ct::setJacobian(jacobian, dofCoordsOriented, physVertexes, quad);
        ct::setJacobianInv (jacobian_inv, jacobian);
        ct::setJacobianDet (jacobian_det, jacobian);
        rst::transpose(jacobian_invT,jacobian_inv);

        DynRankView ConstructWithLabel(dofCoeffsTmp, numCells, basisCardinality, dim);
        rst::matvec(dofCoeffsTmp, jacobian_invT, dofCoeffsOriented);
        at::scalarMultiplyDataData(dofCoeffsPhys, jacobian_det, dofCoeffsTmp);

        //dofCoeffsPhys do not account for orientation. We overwrite dofCoeffs of tangents and edges.
        {
          const ValueType refEdgeLength = 2.0;
          auto numEdgeDOFs = basis.getDofCount(dim-1,0);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type iedge=0; iedge<numElemEdges; iedge++)
              for(ordinal_type j=0; j<numEdgeDOFs; ++j) {
                auto idof = basis.getDofOrdinal(dim-1, iedge, j);
              ValueType normal[dim];
              bool areDifferent = false;
              for(ordinal_type d=0; d <dim; ++d) {
                normal[d] = edgeNormal[i][iedge][d];
                if(std::abs(dofCoeffsPhys(i,idof,d)*refEdgeLength - normal[d])>tol)
                  areDifferent = true;
              }
              if(areDifferent) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "Coefficients of Dof " << idof << " at cell "  << i << ", scaled by the length of the reference edge, are NOT equivalent to the normal of edge " << iedge << "\n";
                *outStream << "Dof Coefficients are: (" << dofCoeffsPhys(i,idof,0) << ", " << dofCoeffsPhys(i,idof,1) <<  ")\n";
                *outStream << "Edge normal is: (" << normal[0] << ", " << normal[1]  << ")\n";
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
              funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
          }
          for(ordinal_type j=0; j<basisCardinality; ++j)
            for(ordinal_type k=0; k<dim; ++k) {
              funDofs(i,j) += fun(physDofCoords(i,j,0), physDofCoords(i,j,1), k) * dofCoeffsPhys(i,j,k);
            }
        }

        //check that fun values are consistent on common edge dofs
        {
          bool areDifferent(false);
          auto numEdgeDOFs = basis.getDofCount(dim-1,0);
          for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
            areDifferent = std::abs(funDofs(0,basis.getDofOrdinal(dim-1,edgeIndex[0],j))
                - funDofs(1,basis.getDofOrdinal(dim-1,edgeIndex[1],j))) > 10*tol;
          }

          if(areDifferent) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Function DOFs on common edge computed using Quad 0 basis functions are not consistent with those computed using Quad 1\n";
            *outStream << "Function DOFs for Quad 0 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funDofs(0,basis.getDofOrdinal(dim-1,edgeIndex[0],j)) << " | (" << physDofCoords(0,basis.getDofOrdinal(dim-1,edgeIndex[0],j),0) << "," << physDofCoords(0,basis.getDofOrdinal(dim-1,edgeIndex[0],j),1) << ") ||";
            *outStream << "\nFunction DOFs for Quad 1 are:";
            for(ordinal_type j=0;j<numEdgeDOFs;j++)
              *outStream << " " << funDofs(1,basis.getDofOrdinal(dim-1,edgeIndex[1],j))<< " | (" << physDofCoords(1,basis.getDofOrdinal(dim-1,edgeIndex[1],j),0) << "," << physDofCoords(1,basis.getDofOrdinal(dim-1,edgeIndex[1],j),1)  << ") ||";
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
        ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
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
            *outStream << "Function values at reference points differ from those computed using basis functions of Quad " << i << "\n";
            *outStream << "Function values at reference points are:\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtPhysRefCoords(i,j,0) << "," << funAtPhysRefCoords(i,j,1) << ")";
            *outStream << "\nFunction values at reference points computed using basis functions are\n";
            for(ordinal_type j=0; j<numRefCoords; ++j)
              *outStream << " (" << funAtRefCoordsOriented(i,j,0) << "," << funAtRefCoordsOriented(i,j,1) << ")";
            *outStream << std::endl;
          }
        }

      }
    } while(std::next_permutation(&reorder[0]+2, &reorder[0]+4)); //reorder vertices of common edge

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

