// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefControlVolume.hpp
    \brief  Definition file for the control volume functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_CONTROL_VOLUME_HPP__
#define __INTREPID2_CELLTOOLS_DEF_CONTROL_VOLUME_HPP__

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

  namespace FunctorCellTools {

    /** \brief  Computes barycenter of polygonal cells
         \param  center            [out] - cell barycenter
         \param  verts             [in]  - cell vertices
    */
    template<typename centerViewType, typename vertViewType>
    KOKKOS_INLINE_FUNCTION
    void getBarycenterPolygon2D(       centerViewType center,
                                 const vertViewType verts) {
      // the enumeration already assumes the ordering of vertices (circling around the polygon)
      const ordinal_type nvert = verts.extent(0);

      center(0) = 0;
      center(1) = 0;
      typename centerViewType::value_type area = 0;
      for (ordinal_type i=0;i<nvert;++i) {
        const ordinal_type j = (i + 1)%nvert;
        const auto scale = verts(i,0)*verts(j,1) - verts(j,0)*verts(i,1);
        center(0) += (verts(i,0) + verts(j,0))*scale;
        center(1) += (verts(i,1) + verts(j,1))*scale;
        area += 0.5*scale;
      }
      center(0) /= (6.0*area);
      center(1) /= (6.0*area);
    }
    
    template<typename midPointViewType, typename nodeMapViewType, typename vertViewType>
    KOKKOS_INLINE_FUNCTION      
    void getMidPoints(       midPointViewType midpts,
                       const nodeMapViewType map,
                       const vertViewType verts) {
      const ordinal_type npts = map.extent(0);
      const ordinal_type dim  = verts.extent(1);
      
      for (ordinal_type i=0;i<npts;++i) {
        // first entry is the number of subcell vertices
        const ordinal_type nvert_per_subcell = map(i, 0);
        for (ordinal_type j=0;j<dim;++j) {
          midpts(i,j) = 0;
          for (ordinal_type k=1;k<=nvert_per_subcell;++k)
            midpts(i,j) += verts(map(i,k),j);
          midpts(i,j) /= nvert_per_subcell;
        }
      }
    }

    template<typename subcvCoordViewType,
             typename cellCoordViewType,
             typename mapViewType>
    /**
     \brief Functor for calculation of sub-control volume coordinates on polygons see Intrepid2::CellTools for more
    */
    struct F_getSubcvCoords_Polygon2D {
            subcvCoordViewType _subcvCoords;
      const cellCoordViewType  _cellCoords;
      const mapViewType        _edgeMap;
      
      KOKKOS_INLINE_FUNCTION
      F_getSubcvCoords_Polygon2D( subcvCoordViewType subcvCoords_,
                                  cellCoordViewType  cellCoords_,
                                  mapViewType edgeMap_ )
        : _subcvCoords(subcvCoords_), _cellCoords(cellCoords_), _edgeMap(edgeMap_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {
        // vertices of cell (P,D)
        const auto verts = Kokkos::subdynrankview( _cellCoords, cell, 
                                                   Kokkos::ALL(), Kokkos::ALL() );
        const ordinal_type nvert = verts.extent(0);
        const ordinal_type dim   = verts.extent(1);

        // control volume coords (N,P,D), here N corresponds to cell vertices
        auto cvCoords = Kokkos::subdynrankview( _subcvCoords, cell, 
                                                Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        // work space for barycenter and midpoints on edges
        typedef typename subcvCoordViewType::value_type value_type;
        value_type buf_center[2], buf_midpts[4*2];
        Kokkos::View<value_type*,Kokkos::AnonymousSpace> center(buf_center, 2);
        Kokkos::View<value_type**,Kokkos::AnonymousSpace> midpts(buf_midpts, 4, 2);

        getBarycenterPolygon2D(center, verts);
        getMidPoints(midpts, _edgeMap, verts);

        for (ordinal_type i=0;i<nvert;++i) {
          for (ordinal_type j=0;j<dim;++j) {
            // control volume is always quad
            cvCoords(i, 0, j) = verts(i, j);
            cvCoords(i, 1, j) = midpts(i, j);
            cvCoords(i, 2, j) = center(j);
            cvCoords(i, 3, j) = midpts((i+nvert-1)%nvert, j);
          }
        }
      }
    };

    template<typename subcvCoordViewType,
             typename cellCoordViewType,
             typename mapViewType>
    /**
     \brief Functor for calculation of sub-control volume coordinates on hexahedra see Intrepid2::CellTools for more
    */
    struct F_getSubcvCoords_Hexahedron {
            subcvCoordViewType _subcvCoords;
      const cellCoordViewType  _cellCoords;
      const mapViewType        _edgeMap, _faceMap;
      
      KOKKOS_INLINE_FUNCTION
      F_getSubcvCoords_Hexahedron( subcvCoordViewType subcvCoords_,
                                   cellCoordViewType  cellCoords_,
                                   mapViewType edgeMap_,
                                   mapViewType faceMap_ )
        : _subcvCoords(subcvCoords_), _cellCoords(cellCoords_), _edgeMap(edgeMap_), _faceMap(faceMap_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {
        // vertices of cell (P,D)
        const auto verts = Kokkos::subdynrankview( _cellCoords, cell, 
                                                   Kokkos::ALL(), Kokkos::ALL() );
        const ordinal_type nvert = verts.extent(0);
        const ordinal_type dim   = verts.extent(1);

        // control volume coords (N,P,D), here N corresponds to cell vertices
        auto cvCoords = Kokkos::subdynrankview( _subcvCoords, cell, 
                                                Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        //  work space for barycenter and midpoints on edges
        typedef typename subcvCoordViewType::value_type value_type;
        value_type buf_center[3], buf_edge_midpts[12*3], buf_face_midpts[6*3];
        Kokkos::View<value_type*,Kokkos::AnonymousSpace> center(buf_center, 3);
        Kokkos::View<value_type**,Kokkos::AnonymousSpace> edge_midpts(buf_edge_midpts, 12, 3);
        Kokkos::View<value_type**,Kokkos::AnonymousSpace> face_midpts(buf_face_midpts,  6, 3);

        // find barycenter
        //Warning! I think this assumes the Hexa is affinely mapped from the reference Hexa
        for (ordinal_type j=0;j<3;++j) {
          center(j) = 0;
          for (ordinal_type i=0;i<nvert;++i) 
            center(j) += verts(i,j);
          center(j) /= nvert;
        }

        getMidPoints(edge_midpts, _edgeMap, verts);
        getMidPoints(face_midpts, _faceMap, verts);

        for (ordinal_type i=0;i<4;++i) {
          const ordinal_type ii = (i+4-1)%4;
          for (ordinal_type j=0;j<dim;++j) {
            
            // set first node of bottom hex to primary cell node
            // and fifth node of upper hex
            cvCoords(i,  0,j) = verts(i,  j);
            cvCoords(i+4,4,j) = verts(i+4,j);

            // set second node of bottom hex to adjacent edge midpoint
            // and sixth node of upper hex
            cvCoords(i,  1,j) = edge_midpts(i,  j);
            cvCoords(i+4,5,j) = edge_midpts(i+4,j);

            // set third node of bottom hex to bottom face midpoint (number 4)
            // and seventh node of upper hex to top face midpoint
            cvCoords(i,  2,j) = face_midpts(4,j);
            cvCoords(i+4,6,j) = face_midpts(5,j);

            // set fourth node of bottom hex to other adjacent edge midpoint
            // and eight node of upper hex to other adjacent edge midpoint
            cvCoords(i,  3,j) = edge_midpts(ii,  j);
            cvCoords(i+4,7,j) = edge_midpts(ii+4,j);

            // set fifth node to vertical edge
            // same as first node of upper hex
            cvCoords(i,  4,j) = edge_midpts(i+8,j);
            cvCoords(i+4,0,j) = edge_midpts(i+8,j);

            // set sixth node to adjacent face midpoint
            // same as second node of upper hex
            cvCoords(i,  5,j) = face_midpts(i,j);
            cvCoords(i+4,1,j) = face_midpts(i,j);

            // set seventh node to barycenter
            // same as third node of upper hex
            cvCoords(i,  6,j) = center(j);
            cvCoords(i+4,2,j) = center(j);

            // set eighth node to other adjacent face midpoint
            // same as fourth node of upper hex
            cvCoords(i,  7,j) = face_midpts(ii,j);
            cvCoords(i+4,3,j) = face_midpts(ii,j);
          } 
        } 
      }
    };


    template<typename subcvCoordViewType,
             typename cellCoordViewType,
             typename mapViewType>
    /**
     \brief Functor for calculation of sub-control volume coordinates on tetrahedra see Intrepid2::CellTools for more
    */
    struct F_getSubcvCoords_Tetrahedron {
            subcvCoordViewType _subcvCoords;
      const cellCoordViewType  _cellCoords;
      const mapViewType        _edgeMap, _faceMap;
      
      KOKKOS_INLINE_FUNCTION
      F_getSubcvCoords_Tetrahedron( subcvCoordViewType subcvCoords_,
                                    cellCoordViewType  cellCoords_,
                                    mapViewType edgeMap_,
                                    mapViewType faceMap_ )
        : _subcvCoords(subcvCoords_), _cellCoords(cellCoords_), _edgeMap(edgeMap_), _faceMap(faceMap_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {
        // ** vertices of cell (P,D)
        const auto verts = Kokkos::subdynrankview( _cellCoords, cell, 
                                                   Kokkos::ALL(), Kokkos::ALL() );
        const ordinal_type nvert = verts.extent(0);
        const ordinal_type dim   = verts.extent(1);

        //  control volume coords (N,P,D), here N corresponds to cell vertices
        auto cvCoords = Kokkos::subdynrankview( _subcvCoords, cell, 
                                                Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        //  work space for barycenter and midpoints on edges
        typedef typename subcvCoordViewType::value_type value_type;
        value_type buf_center[3], buf_edge_midpts[6*3], buf_face_midpts[4*3];
        Kokkos::View<value_type*,Kokkos::AnonymousSpace> center(buf_center, 3);
        Kokkos::View<value_type**,Kokkos::AnonymousSpace> edge_midpts(buf_edge_midpts,  6, 3);
        Kokkos::View<value_type**,Kokkos::AnonymousSpace> face_midpts(buf_face_midpts,  4, 3);

        // find barycenter
        for (ordinal_type j=0;j<3;++j) {
          center(j) = 0;
          for (ordinal_type i=0;i<nvert;++i) 
            center(j) += verts(i,j);
          center(j) /= nvert;
        }
        
        getMidPoints(edge_midpts, _edgeMap, verts);
        getMidPoints(face_midpts, _faceMap, verts);

        for (ordinal_type i=0;i<3;++i) {
          const ordinal_type ii = (i+3-1)%3;
          for (ordinal_type j=0;j<dim;++j) {
            // set first node of bottom hex to primary cell node
            cvCoords(i,0,j) = verts(i,j);

            // set second node of bottom hex to adjacent edge midpoint
            cvCoords(i,1,j) = edge_midpts(i,j);

            // set third node of bottom hex to bottom face midpoint (number 3)
            cvCoords(i,2,j) = face_midpts(3,j);

            // set fourth node of bottom hex to other adjacent edge midpoint
            cvCoords(i,3,j) = edge_midpts(ii,j);

            // set fifth node to vertical edge
            cvCoords(i,4,j) = edge_midpts(i+3,j);

            // set sixth node to adjacent face midpoint
            cvCoords(i,5,j) = face_midpts(i,j);

            // set seventh node to barycenter
            cvCoords(i,6,j) = center(j);

            // set eighth node to other adjacent face midpoint
            cvCoords(i,7,j) = face_midpts(ii,j);
          }
        }

        for (ordinal_type j=0;j<dim;++j) {
          // Control volume attached to fourth node
          // set first node of bottom hex to primary cell node
          cvCoords(3,0,j) = verts(3,j);

          // set second node of bottom hex to adjacent edge midpoint
          cvCoords(3,1,j) = edge_midpts(3,j);

          // set third node of bottom hex to bottom face midpoint (number 3)
          cvCoords(3,2,j) = face_midpts(2,j);

          // set fourth node of bottom hex to other adjacent edge midpoint
          cvCoords(3,3,j) = edge_midpts(5,j);

          // set fifth node to vertical edge
          cvCoords(3,4,j) = edge_midpts(4,j);

          // set sixth node to adjacent face midpoint
          cvCoords(3,5,j) = face_midpts(0,j);

          // set seventh node to barycenter
          cvCoords(3,6,j) = center(j);

          // set eighth node to other adjacent face midpoint
          cvCoords(3,7,j) = face_midpts(1,j);
        }
      }
    };

  }
  
  template<typename DeviceType>
  template<typename subcvCoordValueType, class ...subcvCoordProperties,
           typename cellCoordValueType,  class ...cellCoordProperties>
  void
  CellTools<DeviceType>::
  getSubcvCoords(       Kokkos::DynRankView<subcvCoordValueType,subcvCoordProperties...> subcvCoords,
                  const Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...>   cellCoords,
                  const shards::CellTopology primaryCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(primaryCell), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getSubcvCoords): the primary cell must have a reference cell." );

    INTREPID2_TEST_FOR_EXCEPTION( cellCoords.extent(1) != primaryCell.getVertexCount(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getSubcvCoords): cell coords dimension(1) does not match to # of vertices of the cell." );

    INTREPID2_TEST_FOR_EXCEPTION( cellCoords.extent(2) != primaryCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getSubcvCoords): cell coords dimension(2) does not match to the dimension of the cell." );
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(subcvCoords)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(cellCoords)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::getSubcvCoords(..): input/output views' memory spaces are not compatible with DeviceType");


    // get array dimensions
    const ordinal_type numCells = cellCoords.extent(0);
    //const ordinal_type numVerts = cellCoords.extent(1);
    const ordinal_type spaceDim = cellCoords.extent(2);

    // construct edge and face map for the cell type
    const ordinal_type numEdge = primaryCell.getSubcellCount(1);
    Kokkos::View<ordinal_type**,DeviceType> edgeMap("CellTools::getSubcvCoords::edgeMap", numEdge, 3);
    auto edgeMapHost = Kokkos::create_mirror_view(edgeMap);
    for (ordinal_type i=0;i<numEdge;++i) {
      edgeMapHost(i,0) = primaryCell.getNodeCount(1, i);
      for (ordinal_type j=0;j<edgeMapHost(i,0);++j)
        edgeMapHost(i,j+1) = primaryCell.getNodeMap(1, i, j);
    }

    const ordinal_type numFace = (spaceDim > 2 ? primaryCell.getSubcellCount(2) : 0);
    Kokkos::View<ordinal_type**,DeviceType> faceMap("CellTools::getSubcvCoords::faceMap", numFace, 5);
    auto faceMapHost = Kokkos::create_mirror_view(faceMap);
    for (ordinal_type i=0;i<numFace;++i) {
      faceMapHost(i,0) = primaryCell.getNodeCount(2, i);
      for (ordinal_type j=0;j<faceMapHost(i,0);++j)
        faceMapHost(i,j+1) = primaryCell.getNodeMap(2, i, j);
    }

    Kokkos::deep_copy(edgeMap, edgeMapHost);
    Kokkos::deep_copy(faceMap, faceMapHost);

    // parallel run
    using subcvCoordViewType = decltype(subcvCoords);
    using cellCoordViewType = decltype(cellCoords);
    using mapViewType = decltype(edgeMap);

    const auto loopSize = numCells;
    Kokkos::RangePolicy<ExecSpaceType, Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

    switch (primaryCell.getKey()) {
    case shards::Triangle<3>::key:
    case shards::Quadrilateral<4>::key: {
      // 2D polygon
      typedef FunctorCellTools::F_getSubcvCoords_Polygon2D<subcvCoordViewType,cellCoordViewType,mapViewType> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(subcvCoords, cellCoords, edgeMap) );
      break;
    }
    case shards::Hexahedron<8>::key: {
      typedef FunctorCellTools::F_getSubcvCoords_Hexahedron<subcvCoordViewType,cellCoordViewType,mapViewType> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(subcvCoords, cellCoords, edgeMap, faceMap) );
      break;
    }
    case shards::Tetrahedron<4>::key: {
      typedef FunctorCellTools::F_getSubcvCoords_Tetrahedron<subcvCoordViewType,cellCoordViewType,mapViewType> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(subcvCoords, cellCoords, edgeMap, faceMap) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::getSubcvCoords: the give cell topology is not supported.");
    }
    }
  }
}

#endif
