// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Orientation.hpp
    \brief  Header file for the Intrepid2::Orientation class.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_ORIENTATION_HPP__
#define __INTREPID2_ORIENTATION_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"

namespace Intrepid2 {

  /**
    \brief Orientation encoding and decoding

     Use input and output as pointer arrays which assumes that input/output
     are located on the stack and contiguous.
  */
  class Orientation {
  private:
    template<typename cellVertViewType>
    static void getCellVertexMap(typename cellVertViewType::non_const_value_type *subCellVerts,
                                  ordinal_type &numVerts,
                                  const shards::CellTopology cellTopo,
                                  const cellVertViewType cellVerts,
                                  const ordinal_type subCellDim,
                                  const ordinal_type subCellOrd);
    
    // orientation is always computed on the right-handed coordinates
    template<typename subCellVertType>
    static ordinal_type getOrientation(const subCellVertType subCellVerts[],
                                       const ordinal_type numVerts);


  public:

    /**  Computes and returns the Orientation of a cell, based on vertices' Ids.
         For a 3D cell, the orientation is computed only for cell faces and edges, not for the cell itself.
         For a 2D cell, the orientation of edges is always computed, and the orientation of the cell itself 
         is computed only if the parameter <var><b>isSide</b></var> is true.
         For a 1D cell, the orientation of the cell is computed only when <var><b>isSide</b></var> is true.
    
        \param  cellTopo   [in] - Cell Topology
        \param  cellVerts  [in] - rank-1 view with vertices' Ids
        \param  isSide     [in] - boolean, whether the cell is a side (of a parent cell) for orientation purposes
    */
    template<typename cellVertViewType>
    static Orientation getOrientation(const shards::CellTopology cellTopo,
                                      const cellVertViewType cellVerts,
                                      bool isSide = false);
    
    /**  Returns the ordinal of an edge w.r.t. the parent face, given the ordinals of the
         edge and the face w.r.t. the parent 3D cell.

        \param  subsubcellOrd  [in] - ordinal of the edge w.r.t. the parent (3D) cell
        \param  subcellOrd     [in] - ordinal of the face w.r.t. the parent (3D) cell
        \param  cellTopo       [in] - Cell Topology
    */
    static ordinal_type getEdgeOrdinalOfFace(const ordinal_type subsubcellOrd,
                                             const ordinal_type subcellOrd,
                                             const shards::CellTopology cellTopo);


    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSubcellTangents
    template<typename refTanViewType>
    static void getReferenceEdgeTangent(const refTanViewType &tanE,
                                        const ordinal_type subcellOrd,
                                        const shards::CellTopology cellTopo,
                                        const ordinal_type ort,
                                        const bool is_normalize = true);
    */
    
    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSubcellTangents
    template<typename refTanViewType>
    static void getReferenceFaceTangents(const refTanViewType &tanU,
                                         const refTanViewType &tanV,
                                         const ordinal_type subcellOrd,
                                         const shards::CellTopology cellTopo,
                                         const ordinal_type ort,
                                         const bool is_normalize = true);
    */

    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSideTangentsAndNormal
    template<typename refNormalViewType>
    static void getReferenceFaceNormal(const refNormalViewType &normalV,
                                       const ordinal_type subcellOrd,
                                       const shards::CellTopology cellTopo,
                                       const ordinal_type ort,
                                       const bool is_normalize = true);
    */

      
  private:
    ordinal_type _edgeOrt, _faceOrt;
    
  public:
    /**  Constructor, initializing the orientation so that it's aligned to the reference cell
    */
    KOKKOS_INLINE_FUNCTION
    Orientation();

    /**  Copy constructor
    */
    KOKKOS_DEFAULTED_FUNCTION
    Orientation(const Orientation &b) = default;

    /**  Whether the cell is aligned to reference
    */
    KOKKOS_INLINE_FUNCTION
    bool isAlignedToReference() const;

    /**  Set the orientation of all edges of a cell
        \param  numEdge  [in] - number of edges of the cell
        \param  edgeOrt  [in] - integer array containing the orientations of each edge
    */
    KOKKOS_INLINE_FUNCTION
    void setEdgeOrientation(const ordinal_type numEdge, const ordinal_type edgeOrt[]);

    /**  Get the orientation of the first <var><b>numEdge</b></var> edges of a cell
         \param  edgeOrt [out] - integer array containing the orientations of each edge
         \param  numEdge  [in] - number of edges
    */
    KOKKOS_INLINE_FUNCTION
    void getEdgeOrientation(ordinal_type *edgeOrt, const ordinal_type numEdge) const;

    /**  Set the orientation of all faces of a cell
        \param  numFace  [in] - number of faces of the cell
        \param  faceOrt  [in] - integer array containing the orientations of each face
    */
    KOKKOS_INLINE_FUNCTION
    void setFaceOrientation(const ordinal_type numFace, const ordinal_type faceOrt[]);

    /**  Get the orientation of the first <var><b>numFace</b></var> faces of a cell
         \param  faceOrt [out] - integer array containing the orientations of each face
         \param  numFace  [in] - number of faces
    */
    KOKKOS_INLINE_FUNCTION
    void getFaceOrientation(ordinal_type *faceOrt, const ordinal_type numFace) const;
    
    /**  Returns the orientation converted to a string
    */
    inline std::string to_string() const;
  };
}

inline std::ostream& operator<<(std::ostream& os, const ::Intrepid2::Orientation& orientation)
{
  os << orientation.to_string();
  return os;
}

// include templated function definitions
#include "Intrepid2_OrientationDef.hpp"

#endif
