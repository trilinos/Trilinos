// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellTopologyTags.hpp
    \brief  Definition of cell topology information.

    \verbatim
    CellTopologyTags - compile time information
    template<int N = # of nodal points>
    struct CellTopo {
       base_cell_topology_type - base cell topology type
       dimension - dimension of the topology
       numVert - # of vertices
       numEdge - # of edges
       numFace - # of faces 
       numIntr - # of interior (volume or face according to dimension)
       coords - coordinate points
    
    bool checkPointInclusion(point, threshold)
       - point : input view (D)
       - threshold : extra margin that include cell area (minus threshold is also possible)
    \endverbatim

    \author Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOPOLOGYTAGS_HPP__
#define __INTREPID2_CELLTOPOLOGYTAGS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

//#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
//#include "Intrepid2_Kernels.hpp"

namespace Intrepid2 {

  namespace Impl {


    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Line;

    /** 
     \brief Line topology, 2 nodes
    */
    template<>
    struct Line<2> {
      typedef struct Line<2> base_cell_topology_type;
      enum : int { dimension = 1,
                   numNode = 2,
                   numVert = 2,
                   numEdge = 0,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[2][3]{ {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0} };

      // base topology has this check method
      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
        return (minus_one <= point(0) && point(0) <= plus_one);
      }
    };
    
    /** 
     \brief Line topology, 3 nodes
    */
    template<>
    struct Line<3> {
      typedef struct Line<2> base_cell_topology_type;
      enum : int { dimension = 1,
                   numNode = 3,
                   numVert = 2,
                   numEdge = 0,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[3][3]{ {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };
    
    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Triangle;
    
    /** 
     \brief Triangle topology, 3 nodes
    */
    template<>
    struct Triangle<3> {
      typedef struct Triangle<3> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 3,
                   numVert = 3,
                   numEdge = 3,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[3][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double distance = max( max( -point(0), -point(1) ), point(0) + point(1) - 1.0 );
        return distance < threshold;
      }        
    };
    
    /** 
     \brief Triangle topology, 4 nodes
    */
    template<>
    struct Triangle<4> {
      typedef struct Triangle<3> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 4,
                   numVert = 3,
                   numEdge = 3,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[4][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 1.0/3.0, 1.0/3.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };
    
    /** 
     \brief Triangle topology, 6 nodes
    */
    template<>
    struct Triangle<6> {
      typedef struct Triangle<3> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 6,
                   numVert = 3,
                   numEdge = 3,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[6][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
                                { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Quadrilateral;
    
    /** 
     \brief Quadrilateral topology, 4 nodes
    */
    template<>
    struct Quadrilateral<4> {
      typedef struct Quadrilateral<4> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 4,
                   numVert = 4,
                   numEdge = 4,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[4][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
        return ((minus_one <= point(0) && point(0) <= plus_one) &&
                (minus_one <= point(1) && point(1) <= plus_one));
      }
    };
    
    /** 
     \brief Quadrilateral topology, 8 nodes
    */
    template<>
    struct Quadrilateral<8> {
      typedef struct Quadrilateral<4> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 8,
                   numVert = 4,
                   numEdge = 4,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[8][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
                                { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    /** 
     \brief Quadrilateral topology, 9 nodes
    */
    template<>
    struct Quadrilateral<9> {
      typedef struct Quadrilateral<4> base_cell_topology_type;
      enum : int { dimension = 2,
                   numNode = 9,
                   numVert = 4,
                   numEdge = 4,
                   numFace = 0,
                   numIntr = 1 };
      static constexpr double coords[9][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
                                            { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Tetrahedron;
    
    /** 
     \brief Tetrahedron topology, 4 nodes
    */
    template<>
    struct Tetrahedron<4> {
      typedef struct Tetrahedron<4> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 4,
                   numVert = 4,
                   numEdge = 6,
                   numFace = 4,
                   numIntr = 1 };
      static constexpr double coords[4][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point,
                          const double threshold) {
        const double distance = max( max(-point(0),-point(1)),
                                     max(-point(2), point(0) + point(1) + point(2) - 1) );

        return distance < threshold;
      }
    };

    /** 
     \brief Tetrahedron topology, 8 nodes
    */
    template<>
    struct Tetrahedron<8> {
      typedef struct Tetrahedron<4> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 8,
                   numVert = 4,
                   numEdge = 6,
                   numFace = 4,
                   numIntr = 1 };
      static constexpr double coords[8][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
                                { 1/3, 0.0, 1/3}, { 1/3, 1/3, 1/3}, { 1/3, 1/3, 0.0}, { 0.0, 1/3, 1/3} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    /** 
     \brief Tetrahedron topology, 10 nodes
    */
    template<>
    struct Tetrahedron<10> {
      typedef struct Tetrahedron<4> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 10,
                   numVert = 4,
                   numEdge = 6,
                   numFace = 4,
                   numIntr = 1 };
      static constexpr double coords[10][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
                                 { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    /** 
     \brief Tetrahedron topology, 11 nodes
    */
    template<>
    struct Tetrahedron<11> {
      typedef struct Tetrahedron<4> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 11,
                   numVert = 4,
                   numEdge = 6,
                   numFace = 4,
                   numIntr = 1 };
      static constexpr double coords[11][3]{ { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
                                 { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5} };
      
      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Hexahedron;
    
    /** 
     \brief Hexahedron topology, 8 nodes
    */
    template<>
    struct Hexahedron<8> {
      typedef struct Hexahedron<8> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 8,
                   numVert = 8,
                   numEdge = 12,
                   numFace = 6,
                   numIntr = 1 };
      static constexpr double coords[8][3]{ {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
                                {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
        return ((minus_one <= point(0) && point(0) <= plus_one) &&
                (minus_one <= point(1) && point(1) <= plus_one) &&
                (minus_one <= point(2) && point(2) <= plus_one));
      }
    };

    /** 
     \brief Hexahedron topology, 20 nodes
    */
    template<>
    struct Hexahedron<20> {
      typedef struct Hexahedron<8> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 20,
                   numVert = 8,
                   numEdge = 12,
                   numFace = 6,
                   numIntr = 1 };
      static constexpr double coords[20][3]{ {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
                                 {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
                                 { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
                                 {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
                                 { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }

    };

    /** 
     \brief Hexahedron topology, 27 nodes
    */
    template<>
    struct Hexahedron<27> {
      typedef struct Hexahedron<8> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 27,
                   numVert = 8,
                   numEdge = 12,
                   numFace = 6,
                   numIntr = 1 };
      static constexpr double coords[27][3]{ {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
                                 {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
                                 { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
                                 {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
                                 { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
                                 { 0.0, 0.0, 0.0},
                                 { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Pyramid;

    /** 
     \brief Pyramid topology, 5 nodes
    */
    template<>
    struct Pyramid<5> {
      typedef struct Pyramid<5> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 5,
                   numVert = 5,
                   numEdge = 8,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[5][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double minus_one = -1.0 - threshold, plus_one = 1.0 + threshold, minus_zero = -threshold;
        const double left  = minus_one + point(2);
        const double right =  plus_one - point(2);
        return ((left       <= point(0) && point(0) <= right) &&
                (left       <= point(1) && point(1) <= right) &&
                (minus_zero <= point(2) && point(2) <= plus_one));
      }
    };

    /** 
     \brief Pyramid topology, 13 nodes
    */
    template<>
    struct Pyramid<13> {
      typedef struct Pyramid<5> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 13,
                   numVert = 5,
                   numEdge = 8,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[13][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
                                 { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
                                 {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };

    /** 
     \brief Pyramid topology, 14 nodes
    */
    template<>
    struct Pyramid<14> {
      typedef struct Pyramid<5> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 14,
                   numVert = 5,
                   numEdge = 8,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[14][3]{ {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
                                 { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
                                 {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };
    
    // ---------------------------------------------------------------------------------------

    template<int N>
    struct Wedge;

    /** 
     \brief Wedge topology, 6 nodes
    */
    template<>
    struct Wedge<6> {
      typedef struct Wedge<6> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 6,
                   numVert = 6,
                   numEdge = 9,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[6][3]{ { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        const double minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
        const double distance = max( max( -point(0), -point(1) ), point(0) + point(1) - 1 );
        return (distance < threshold && (minus_one <= point(2) && point(2) <= plus_one));
      }
    };

    /** 
     \brief Wedge topology, 15 nodes
    */
    template<>
    struct Wedge<15> {
      typedef struct Wedge<6> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 15,
                   numVert = 6,
                   numEdge = 9,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[15][3]{ { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
                                 { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
                                 { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0} };

      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };
    
    /** 
     \brief Wedge topology, 18 nodes
    */
    template<>
    struct Wedge<18> {
      typedef struct Wedge<6> base_cell_topology_type;
      enum : int { dimension = 3,
                   numNode = 18,
                   numVert = 6,
                   numEdge = 9,
                   numFace = 5,
                   numIntr = 1 };
      static constexpr double coords[18][3]{ { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
                                 { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
                                 { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
                                 { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0} };

      
      template<typename PointViewType>
      KOKKOS_INLINE_FUNCTION
      static bool
      checkPointInclusion(const PointViewType &point, 
                          const double threshold) {
        return base_cell_topology_type::checkPointInclusion(point, threshold);
      }
    };
    
  }
}

#endif

