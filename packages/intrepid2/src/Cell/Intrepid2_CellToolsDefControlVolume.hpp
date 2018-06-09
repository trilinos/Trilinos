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

    /** \brief  Computes cell barycenter
         \param  center            [out] - cell barycenter
         \param  verts             [in]  - cell vertices
    */
    template<typename centerViewType, typename vertViewType>
    KOKKOS_INLINE_FUNCTION
    void getBaryCenter(       centerViewType center,
                        const vertViewType verts) {
      // the enumeration already assumes the ordering of vertices (circling around the polygon)
      const ordinal_type nvert = verts.extent(0);
      const ordinal_type dim   = verts.extent(1);
      
      switch (dim) {
      case 2: {
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
        break;
      }
      case 3: {
        // This method works fine for simplices, but for other 3-d shapes
        // is not precisely accurate. Could replace with approximate integration
        // perhaps.
        for (ordinal_type j=0;j<dim;++j) {
          center(j) = 0;
          for (ordinal_type i=0;i<nvert;++i) 
            center(j) += verts(i,j);
          center(j) /= nvert;
        }
        break;
      }
      }
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
        Kokkos::View<value_type*,Kokkos::Impl::ActiveExecutionMemorySpace> center(buf_center, 2);
        Kokkos::View<value_type**,Kokkos::Impl::ActiveExecutionMemorySpace> midpts(buf_midpts, 4, 2);

        getBaryCenter(center, verts);
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
        // const ordinal_type nvert = verts.extent(0);
        const ordinal_type dim   = verts.extent(1);

        // control volume coords (N,P,D), here N corresponds to cell vertices
        auto cvCoords = Kokkos::subdynrankview( _subcvCoords, cell, 
                                                Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        //  work space for barycenter and midpoints on edges
        typedef typename subcvCoordViewType::value_type value_type;
        value_type buf_center[3], buf_edge_midpts[12*3], buf_face_midpts[6*3];
        Kokkos::View<value_type*,Kokkos::Impl::ActiveExecutionMemorySpace> center(buf_center, 3);
        Kokkos::View<value_type**,Kokkos::Impl::ActiveExecutionMemorySpace> edge_midpts(buf_edge_midpts, 12, 3);
        Kokkos::View<value_type**,Kokkos::Impl::ActiveExecutionMemorySpace> face_midpts(buf_face_midpts,  6, 3);

        getBaryCenter(center, verts);
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
        //const ordinal_type nvert = verts.extent(0);
        const ordinal_type dim   = verts.extent(1);

        //  control volume coords (N,P,D), here N corresponds to cell vertices
        auto cvCoords = Kokkos::subdynrankview( _subcvCoords, cell, 
                                                Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        //  work space for barycenter and midpoints on edges
        typedef typename subcvCoordViewType::value_type value_type;
        value_type buf_center[3], buf_edge_midpts[6*3], buf_face_midpts[4*3];
        Kokkos::View<value_type*,Kokkos::Impl::ActiveExecutionMemorySpace> center(buf_center, 3);
        Kokkos::View<value_type**,Kokkos::Impl::ActiveExecutionMemorySpace> edge_midpts(buf_edge_midpts,  6, 3);
        Kokkos::View<value_type**,Kokkos::Impl::ActiveExecutionMemorySpace> face_midpts(buf_face_midpts,  4, 3);

        getBaryCenter(center, verts);
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
  
  template<typename SpT>
  template<typename subcvCoordValueType, class ...subcvCoordProperties,
           typename cellCoordValueType,  class ...cellCoordProperties>
  void
  CellTools<SpT>::
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

    // get array dimensions
    const ordinal_type numCells = cellCoords.extent(0);
    //const ordinal_type numVerts = cellCoords.extent(1);
    const ordinal_type spaceDim = cellCoords.extent(2);

    // construct edge and face map for the cell type
    const ordinal_type numEdge = primaryCell.getSubcellCount(1);
    Kokkos::View<ordinal_type**,Kokkos::LayoutRight,Kokkos::HostSpace> edgeMapHost("CellTools::getSubcvCoords::edgeMapHost", numEdge, 3);
    for (ordinal_type i=0;i<numEdge;++i) {
      edgeMapHost(i,0) = primaryCell.getNodeCount(1, i);
      for (ordinal_type j=0;j<edgeMapHost(i,0);++j)
        edgeMapHost(i,j+1) = primaryCell.getNodeMap(1, i, j);
    }

    const ordinal_type numFace = (spaceDim > 2 ? primaryCell.getSubcellCount(2) : 0);
    Kokkos::View<ordinal_type**,Kokkos::LayoutRight,Kokkos::HostSpace> faceMapHost("CellTools::getSubcvCoords::faceMapHost", numFace, 5);
    for (ordinal_type i=0;i<numFace;++i) {
      faceMapHost(i,0) = primaryCell.getNodeCount(2, i);
      for (ordinal_type j=0;j<faceMapHost(i,0);++j)
        faceMapHost(i,j+1) = primaryCell.getNodeMap(2, i, j);
    }

    // create mirror to device
    auto edgeMap = Kokkos::create_mirror_view(typename SpT::memory_space(), edgeMapHost);
    auto faceMap = Kokkos::create_mirror_view(typename SpT::memory_space(), faceMapHost);

    Kokkos::deep_copy(edgeMap, edgeMapHost);
    Kokkos::deep_copy(faceMap, faceMapHost);

    // parallel run
    typedef Kokkos::DynRankView<subcvCoordValueType,subcvCoordProperties...> subcvCoordViewType;
    typedef Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...>   cellCoordViewType;
    typedef Kokkos::View<ordinal_type**,Kokkos::LayoutRight,SpT>             mapViewType;

    typedef typename ExecSpace<typename subcvCoordViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    const auto loopSize = numCells;
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

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
