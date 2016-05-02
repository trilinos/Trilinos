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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                             Control Volume Coordinates                                     //
  //                                                                                            //
  //============================================================================================//

  template<class Scalar>
  template<class ArrayCVCoord, class ArrayCellCoord>
  void CellTools<Scalar>::getSubCVCoords(ArrayCVCoord & subCVCoords,
                                         const ArrayCellCoord & cellCoords,
                                         const shards::CellTopology& primaryCell)
  {

    // get array dimensions
    index_type numCells        = static_cast<index_type>(cellCoords.dimension(0));
    index_type numNodesPerCell = static_cast<index_type>(cellCoords.dimension(1));
    index_type spaceDim        = static_cast<index_type>(cellCoords.dimension(2));

    // num edges per primary cell
    int numEdgesPerCell = primaryCell.getEdgeCount();

    // num faces per primary cell
    int numFacesPerCell = 0;
    if (spaceDim > 2){
      numFacesPerCell = primaryCell.getFaceCount();
    }

    // get cell centroids
    FieldContainer<Scalar> barycenter(numCells,spaceDim);
    getBarycenter(barycenter,cellCoords);

    // loop over cells
    for (index_type icell = 0; icell < numCells; icell++){

      // get primary edge midpoints
      FieldContainer<Scalar> edgeMidpts(numEdgesPerCell,spaceDim);
      for (int iedge = 0; iedge < numEdgesPerCell; iedge++){
        for (index_type idim = 0; idim < spaceDim; idim++){

          int node0 = primaryCell.getNodeMap(1,iedge,0);
          int node1 = primaryCell.getNodeMap(1,iedge,1);
          edgeMidpts(iedge,idim) = (cellCoords(icell,node0,idim) +
                                    cellCoords(icell,node1,idim))/2.0;

        } // end loop over dimensions
      } // end loop over cell edges

      // get primary face midpoints in 3-D
      int numNodesPerFace;
      FieldContainer<Scalar> faceMidpts(numFacesPerCell,spaceDim);
      if (spaceDim > 2) {
        for (int iface = 0; iface < numFacesPerCell; iface++){
          numNodesPerFace = primaryCell.getNodeCount(2,iface);

          for (int idim = 0; idim < spaceDim; idim++){

            for (int inode0 = 0; inode0 < numNodesPerFace; inode0++) {
              int node1 = primaryCell.getNodeMap(2,iface,inode0);
              faceMidpts(iface,idim) += cellCoords(icell,node1,idim)/numNodesPerFace;
            }

          } // end loop over dimensions
        } // end loop over cell faces
      }
      // define coordinates for subcontrol volumes
      switch(primaryCell.getKey() ) {

        // 2-d  parent cells
      case shards::Triangle<3>::key:
      case shards::Quadrilateral<4>::key:

        for (int inode = 0; inode < numNodesPerCell; inode++){
          for (index_type idim = 0; idim < spaceDim; idim++){

            // set first node to primary cell node
            subCVCoords(icell,inode,0,idim) = cellCoords(icell,inode,idim);

            // set second node to adjacent edge midpoint
            subCVCoords(icell,inode,1,idim) = edgeMidpts(inode,idim);

            // set third node to cell barycenter
            subCVCoords(icell,inode,2,idim) = barycenter(icell,idim);

            // set fourth node to other adjacent edge midpoint
            int jnode = numNodesPerCell-1;
            if (inode > 0) jnode = inode - 1;
            subCVCoords(icell,inode,3,idim) = edgeMidpts(jnode,idim);

          } // dim loop
        } // node loop

        break;

      case shards::Hexahedron<8>::key:

        for (index_type idim = 0; idim < spaceDim; idim++){

          // loop over the horizontal quads that define the subcontrol volume coords
          for (int icount = 0; icount < 4; icount++){

            // set first node of bottom hex to primary cell node
            // and fifth node of upper hex
            subCVCoords(icell,icount,0,idim) = cellCoords(icell,icount,idim);
            subCVCoords(icell,icount+4,4,idim) = cellCoords(icell,icount+4,idim);

            // set second node of bottom hex to adjacent edge midpoint
            // and sixth node of upper hex
            subCVCoords(icell,icount,1,idim) = edgeMidpts(icount,idim);
            subCVCoords(icell,icount+4,5,idim) = edgeMidpts(icount+4,idim);

            // set third node of bottom hex to bottom face midpoint (number 4)
            // and seventh node of upper hex to top face midpoint
            subCVCoords(icell,icount,2,idim) = faceMidpts(4,idim);
            subCVCoords(icell,icount+4,6,idim) = faceMidpts(5,idim);

            // set fourth node of bottom hex to other adjacent edge midpoint
            // and eight node of upper hex to other adjacent edge midpoint
            int jcount = 3;
            if (icount > 0) jcount = icount - 1;
            subCVCoords(icell,icount,3,idim) = edgeMidpts(jcount,idim);
            subCVCoords(icell,icount+4,7,idim) = edgeMidpts(jcount+4,idim);

            // set fifth node to vertical edge
            // same as first node of upper hex
            subCVCoords(icell,icount,4,idim) = edgeMidpts(icount+numNodesPerCell,idim);
            subCVCoords(icell,icount+4,0,idim) = edgeMidpts(icount+numNodesPerCell,idim);

            // set sixth node to adjacent face midpoint
            // same as second node of upper hex
            subCVCoords(icell,icount,5,idim) = faceMidpts(icount,idim);
            subCVCoords(icell,icount+4,1,idim) = faceMidpts(icount,idim);

            // set seventh node to barycenter
            // same as third node of upper hex
            subCVCoords(icell,icount,6,idim) = barycenter(icell,idim);
            subCVCoords(icell,icount+4,2,idim) = barycenter(icell,idim);

            // set eighth node to other adjacent face midpoint
            // same as fourth node of upper hex
            jcount = 3;
            if (icount > 0) jcount = icount - 1;
            subCVCoords(icell,icount,7,idim) = faceMidpts(jcount,idim);
            subCVCoords(icell,icount+4,3,idim) = faceMidpts(jcount,idim);

          } // count loop

        } // dim loop

        break;

      case shards::Tetrahedron<4>::key:

        for (index_type idim = 0; idim < spaceDim; idim++){

          // loop over the three bottom nodes
          for (int icount = 0; icount < 3; icount++){

            // set first node of bottom hex to primary cell node
            subCVCoords(icell,icount,0,idim) = cellCoords(icell,icount,idim);

            // set second node of bottom hex to adjacent edge midpoint
            subCVCoords(icell,icount,1,idim) = edgeMidpts(icount,idim);

            // set third node of bottom hex to bottom face midpoint (number 3)
            subCVCoords(icell,icount,2,idim) = faceMidpts(3,idim);

            // set fourth node of bottom hex to other adjacent edge midpoint
            int jcount = 2;
            if (icount > 0) jcount = icount - 1;
            subCVCoords(icell,icount,3,idim) = edgeMidpts(jcount,idim);

            // set fifth node to vertical edge
            subCVCoords(icell,icount,4,idim) = edgeMidpts(icount+3,idim);

            // set sixth node to adjacent face midpoint
            subCVCoords(icell,icount,5,idim) = faceMidpts(icount,idim);

            // set seventh node to barycenter
            subCVCoords(icell,icount,6,idim) = barycenter(icell,idim);

            // set eighth node to other adjacent face midpoint
            jcount = 2;
            if (icount > 0) jcount = icount - 1;
            subCVCoords(icell,icount,7,idim) = faceMidpts(jcount,idim);

          } //count loop

            // Control volume attached to fourth node
          // set first node of bottom hex to primary cell node
          subCVCoords(icell,3,0,idim) = cellCoords(icell,3,idim);

          // set second node of bottom hex to adjacent edge midpoint
          subCVCoords(icell,3,1,idim) = edgeMidpts(3,idim);

          // set third node of bottom hex to bottom face midpoint (number 3)
          subCVCoords(icell,3,2,idim) = faceMidpts(2,idim);

          // set fourth node of bottom hex to other adjacent edge midpoint
          subCVCoords(icell,3,3,idim) = edgeMidpts(5,idim);

          // set fifth node to vertical edge
          subCVCoords(icell,3,4,idim) = edgeMidpts(4,idim);

          // set sixth node to adjacent face midpoint
          subCVCoords(icell,3,5,idim) = faceMidpts(0,idim);

          // set seventh node to barycenter
          subCVCoords(icell,3,6,idim) = barycenter(icell,idim);

          // set eighth node to other adjacent face midpoint
          subCVCoords(icell,3,7,idim) = faceMidpts(1,idim);

        } // dim loop

        break;

      default:
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                      ">>> ERROR (getSubCVCoords: invalid cell topology.");
      } // cell key

    } // cell loop

  } // getSubCVCoords

  template<class Scalar>
  template<class ArrayCent, class ArrayCellCoord>
  void CellTools<Scalar>::getBarycenter(ArrayCent & barycenter, const ArrayCellCoord & cellCoords)
  {
    // get array dimensions
    index_type numCells        = static_cast<index_type>(cellCoords.dimension(0));
    index_type numVertsPerCell = static_cast<index_type>(cellCoords.dimension(1));
    index_type spaceDim        = static_cast<index_type>(cellCoords.dimension(2));

    if (spaceDim == 2)
      {
        // Method for general polygons
        for (index_type icell = 0; icell < numCells; icell++){

          FieldContainer<Scalar> cell_centroid(spaceDim);
          Scalar area = 0;

          for (index_type inode = 0; inode < numVertsPerCell; inode++){

            index_type jnode = inode + 1;
            if (jnode >= numVertsPerCell) {
              jnode = 0;
            }

            Scalar area_mult = cellCoords(icell,inode,0)*cellCoords(icell,jnode,1)
              - cellCoords(icell,jnode,0)*cellCoords(icell,inode,1);
            cell_centroid(0) += (cellCoords(icell,inode,0) + cellCoords(icell,jnode,0))*area_mult;
            cell_centroid(1) += (cellCoords(icell,inode,1) + cellCoords(icell,jnode,1))*area_mult;

            area += 0.5*area_mult;
          }

          barycenter(icell,0) = cell_centroid(0)/(6.0*area);
          barycenter(icell,1) = cell_centroid(1)/(6.0*area);
        }

      }
    else
      {
        // This method works fine for simplices, but for other 3-d shapes
        // is not precisely accurate. Could replace with approximate integration
        // perhaps.
        for (index_type icell = 0; icell < numCells; icell++){

          FieldContainer<Scalar> cell_centroid(spaceDim);

          for (index_type inode = 0; inode < numVertsPerCell; inode++){
            for (index_type idim = 0; idim < spaceDim; idim++){
              cell_centroid(idim) += cellCoords(icell,inode,idim)/numVertsPerCell;
            }
          }
          for (index_type idim = 0; idim < spaceDim; idim++){
            barycenter(icell,idim) = cell_centroid(idim);
          }
        }
      }

  } 
}

#endif
