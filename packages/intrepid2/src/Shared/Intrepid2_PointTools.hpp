// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PointTools.hpp
    \brief  Header file for Intrepid2::PointTools class to provide utilities for
            barycentric coordinates, equispaced lattices, and
            warp-blend point distrubtions.
    \author Created by R. Kirby
 */

#ifndef __INTREPID2_POINTTOOLS_HPP__
#define __INTREPID2_POINTTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_Polylib.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {

/** \class Intrepid2::PointTools
      \brief Utility class that provides methods for calculating
      distributions of points on different cells

      Simplicial lattices in PointTools are sets
      of points with certain ordering properties.
      They are used for defining degrees of freedom
      for higher order finite elements.

      Each lattice has an "order".  In
      general, this is the same as the cardinality
      of the polynomial space of degree "order".
      In terms of binomial coefficients, this is
      binomial(order+d,order) for the simplex in
      d dimensions.  On the line,
      the size is order+1.  On the triangle
      and tetrahedron, there are
      (order+1)(order+2)/2 and (order+1)(order+2)(order+3)/6,
      respectively.

      The points are ordered lexicographically from low to
      high in increasing spatial dimension.  For example,
      the line lattice of order 3 looks like:

      \verbatim
      x--x--x--x
      \endverbatim

      where "x" denotes a point location.
      These are ordered from left to right, so that
      the points are labeled:

      \verbatim
      0--1--2--3
      \endverbatim

      The triangular lattice of order 3 is 

      \verbatim
      x
      |\
      | \
      x  x
      |   \
      |    \
      x  x  x
      |      \
      |       \
      x--x--x--x
      \endverbatim

      The ordering starts in the bottom left and
      increases first from left to right.  The ordering
      is 

      \verbatim
      9
      |\
      | \
      8  7
      |   \
      |    \
      4  5  6
      |      \
      |       \
      0--1--2--3
      \endverbatim

      Tetrahedral lattices are similar but difficult to
      draw with ASCII art.

      Each lattice also has an "offset", which indicates
      a number of layers of points on the bounary taken away.
      All of the lattices above have a 0 offest.  In Intrepid,
      typically only offset = 0 or 1 will be used.  The offset=1
      case is used to generate sets of points properly inside
      a given simplex.  These are used, for example, to construct
      points internal to an edge or face for H(curl) and H(div)
      finite elements.

      For example, for a line lattice with order = 3 and
      offset = 1, the points will look like

      \verbatim
      ---x--x---
      \endverbatim

      and a triangle with order=3 and offset=1 will 
      contain a single point 

      \verbatim
      .
      |\
      | \
      |  \
      |   \
      |    \
      |  x  \
      |      \
      |       \
      |--------\
      \endverbatim

      When points on lattices with nonzero offset are numbered,
      the are numbered contiguously from 0, so that the line and
      triangle above are respectively

      \verbatim
      ---0--1---
      \endverbatim

      \verbatim
      .
      |\
      | \
      |  \
      |   \
      |    \
      |  0  \
      |      \
      |       \
      |--------\
      \endverbatim

      Additionally, two types of point distributions are currently support.
      The points may be on an equispaced lattice, which is easy to compute
      but can lead to numerical ill-conditioning in finite element bases
      and stiffness matrices.  Alternatively, the warp-blend points of
      Warburton are provided on each lattice (which are just the
      Gauss-Lobatto points on the line).

 */
class PointTools {
public:
  /** \brief Computes the number of points in a lattice of a given order
        on a simplex (currently disabled for
        other cell types).
        If offset == 0,
        the lattice will include only include the vertex points if order == 1,
        and will include edge midpoints if order == 2, and so on.  
        In particular, this is the dimension of polynomials of degree "order"
        on the given simplex.
        The offset argument is used to indicate that the layer of points on the
        boundary is omitted (if offset == 1).  For greater offsets, more layers
        are omitteed.

        \param  cellType     [in]  - type of reference cell (currently only supports the simplex)
        \param  order        [in]  - order of the lattice
        \param  offset       [in]  - the number of boundary layers to omit

   */
  inline
  static ordinal_type
  getLatticeSize( const shards::CellTopology cellType,
      const ordinal_type order,
      const ordinal_type offset = 0 );

  /** \brief Computes a lattice of points of a given
        order on a reference simplex, quadrilateral or hexahedron
        (currently disabled for other cell types).
        The output array is (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  cellType    [in]  - type of reference cell (currently only supports the simplex)
        \param  order       [in]  - number of points per edge, minus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced and
        warp/blend points are supported
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getLattice(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const shards::CellTopology cellType,
      const ordinal_type order,
      const ordinal_type offset = 0 ,
      const EPointType pointType = POINTTYPE_EQUISPACED );

  /** \brief Computes a lattice of points of a given
        order on a reference line.
        The output array is (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per edge, minus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced and
        warp/blend points are supported
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getLatticeLine(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 ,
      const EPointType pointType = POINTTYPE_EQUISPACED );

  /** \brief Computes a lattice of points of a given
        order on a reference triangle.
        The output array is (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per edge, minus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced and
        warp/blend points are supported
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getLatticeTriangle(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 ,
      const EPointType pointType = POINTTYPE_EQUISPACED );

  /** \brief Computes a lattice of points of a given
        order on a reference tetrahedron.
        The output array is (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per edge, minus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced and
        warp/blend points are supported
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getLatticeTetrahedron(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 ,
      const EPointType pointType = POINTTYPE_EQUISPACED );

  /** \brief Computes a lattice of points of a given order on a reference pyramid.
        The output array is (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per edge, minus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced points are supported.
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getLatticePyramid(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 ,
      const EPointType pointType = POINTTYPE_EQUISPACED );

  /** Retrieves the Gauss-Legendre points from PolyLib, but lets us
	do it in an arbitrary ArrayType.
	\param  points      [out] - Output array of point coords (P,)
	\param  order       [out] - number of Gauss points - 1
   */
  template<typename pointValueType, class ...pointProperties>
  static void getGaussPoints(      Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order );

private:

  /** \brief Computes an equispaced lattice of a given
        order on the reference line [-1,1].  The output array is
        (P,1), where
        \code
        P - number of points per cell
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
        minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
        in per boundary

   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getEquispacedLatticeLine(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 );


  /** \brief Returns the Gauss-Lobatto points of a given
        order on the reference line [-1,1].  The output array is
        (P,1), where
        \code
        P - number of points
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
        minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
        in per boundary

   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getWarpBlendLatticeLine(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order ,
      const ordinal_type offset = 0 );


  // /* \brief Converts Cartesian coordinates to barycentric coordinates
  //     on a batch of triangles.
  //     The input array cartValues is (C,P,2)
  //     The output array baryValues is (C,P,3).
  //     The input array vertices is (C,3,2), where
  //     \code
  //     C - num. integration domains
  //     P - number of points per cell
  //     \endcode

  //     \param  baryValues      [out] - Output array of barycentric coords
  //     \param  cartValues      [in]  - Input array of Cartesian coords
  //     \param  vertices        [out] - Vertices of each cell.

  // */
  // template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  // static void cartToBaryTriangle( ArrayTypeOut & baryValues ,
  //                                 const ArrayTypeIn1 & cartValues ,
  //                                 const ArrayTypeIn2 & vertices );


  // /* \brief Converts barycentric coordinates to Cartesian coordinates
  //     on a batch of triangles.
  //     The input array baryValues is (C,P,3)
  //     The output array cartValues is (C,P,2).
  //     The input array vertices is (C,3,2), where
  //     \code
  //     C - num. integration domains
  //     P - number of points per cell
  //     D - is the spatial dimension
  //     \endcode

  //     \param  baryValues      [out] - Output array of barycentric coords
  //     \param  cartValues      [in]  - Input array of Cartesian coords
  //     \param  vertices        [out] - Vertices of each cell.

  // */
  // template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  // static void baryToCartTriangle( ArrayTypeOut & cartValues ,
  //                                 const ArrayTypeIn1 & baryValues ,
  //                                 const ArrayTypeIn2 & vertices );


  // /* \brief Converts Cartesian coordinates to barycentric coordinates
  //     on a batch of tetrahedra.
  //     The input array cartValues is (C,P,3)
  //     The output array baryValues is (C,P,4).
  //     The input array vertices is (C,4,3), where
  //     \code
  //     C - num. integration domains
  //     P - number of points per cell
  //     D - is the spatial dimension
  //     \endcode

  //     \param  baryValues      [out] - Output array of barycentric coords
  //     \param  cartValues      [in]  - Input array of Cartesian coords
  //     \param  vertices        [out] - Vertices of each cell.

  // */
  // template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  // static void cartToBaryTetrahedron( ArrayTypeOut & baryValues ,
  //                                    const ArrayTypeIn1 & cartValues ,
  //                                    const ArrayTypeIn2 & vertices );

  // /* \brief Converts barycentric coordinates to Cartesian coordinates
  //     on a batch of tetrahedra.
  //     The input array baryValues is (C,P,4)
  //     The output array cartValues is (C,P,3).
  //     The input array vertices is (C,4,3), where
  //     \code
  //     C - num. integration domains
  //     P - number of points per cell
  //     D - is the spatial dimension
  //     \endcode

  //     \param  baryValues      [out] - Output array of barycentric coords
  //     \param  cartValues      [in]  - Input array of Cartesian coords
  //     \param  vertices        [out] - Vertices of each cell.

  // */
  // template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  // static void baryToCartTetrahedron( ArrayTypeOut & cartValues ,
  //                                    const ArrayTypeIn1 & baryValues ,
  //                                    const ArrayTypeIn2 & vertices );




   /** \brief Computes an equispaced lattice of a given
       order on the reference triangle.  The output array is
       (P,2), where
       \code
       P - number of points
       \endcode

       \param  points      [out] - Output Kokkos::DynRankView of point coords
       \param  order       [in]  - The lattice has order + 1 points,
       minus any skipped by offset
       \param  offset      [in]  - Number of points on boundary to skip coming
       in from boundary

   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getEquispacedLatticeTriangle(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
      const ordinal_type order,
      const ordinal_type offset = 0 );

   /** \brief Computes an equispaced lattice of a given
       order on the reference tetrahedron.  The output array is
       (P,3), where
       \code
       P - number of points
       \endcode

       \param  points      [out] - Output Kokkos::DynRankView of point coords
       \param  order       [in]  - The lattice has (order + 1) * (order + 2) / 2 points, minus any skipped by offset
       \param  offset      [in]  - Number of points on boundary to skip coming
       in from boundary
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getEquispacedLatticeTetrahedron( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                                   const ordinal_type order ,
                                   const ordinal_type offset = 0 );


  /** \brief Computes an equispaced lattice of a given
      order on the reference pyramid.  The output array is
      (P,3), where
      \code
      P - number of points
      \endcode

      \param  points      [out] - Output Kokkos::DynRankView of point coords
      \param  order       [in]  - The lattice has (order + 1) * (order + 2) * (2*order + 3) / 6 points, minus any skipped by offset
      \param  offset      [in]  - Number of points on boundary to skip coming in from boundary

  */
 template<typename pointValueType, class ...pointProperties>
 static void
 getEquispacedLatticePyramid(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
     const ordinal_type order,
     const ordinal_type offset = 0 );


   /** \brief interpolates Warburton's warp function on the line
       \param  warp        [out] - Kokkos::DynRankView of the amount to warp each point
       \param  order       [in]  - The polynomial order
       \param  xnodes      [in]  - Kokkos::DynRankView of node locations to interpolate
       \param  xout        [in]  - Kokkos::DynRankView of the warpfunction at xout, +/- 1 roots deflated
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  warpFactor( Kokkos::DynRankView<pointValueType,pointProperties...>  warp ,
              const ordinal_type order ,
              const Kokkos::DynRankView<pointValueType,pointProperties...> xnodes ,
              const Kokkos::DynRankView<pointValueType,pointProperties...> xout
              );

   /** \brief Returns Warburton's warp-blend points of a given
       order on the reference triangle.  The output array is
       (P,2), where
       \code
       P - number of points
       \endcode

       \param  points      [out] - Kokkos::DynRankView of point coords
       \param  order       [in]  - The lattice has order + 1 points,
       minus any skipped by offset
       \param  offset      [in]  - Number of points on boundary to skip coming
       in per boundary

   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getWarpBlendLatticeTriangle( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                               const ordinal_type order ,
                               const ordinal_type offset = 0 );

   /** \brief Returns Warburton's warp-blend points of a given
       order on the reference tetrahedron.  The output array is
       (P,3), where
       \code
       P - number of points
       \endcode

       \param  points      [out] - Kokkos::DynRankView of point coords
       \param  order       [in]  - The lattice has order + 1 points,
       minus any skipped by offset
       \param  offset      [in]  - Number of points on boundary to skip coming
       in per boundary
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  getWarpBlendLatticeTetrahedron( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                                  const ordinal_type order ,
                                  const ordinal_type offset = 0 );


   /** \brief This is used internally to compute the tetrahedral warp-blend points one each face
       \param dxy     [out] - Kokkos::DynRankView containing the amount to shift each point in the x and y direction
       \param order   [in] - the order of lattice
       \param pval    [in] - the "alpha" term in the warping function
       \param L1      [in] - Kokkos::DynRankView of the first barycentric coordinate of the input points
       \param L2      [in] - Kokkos::DynRankView of the second barycentric coordinate of the input points
       \param L3      [in] - Kokkos::DynRankView of the third barycentric coordinate of the input points
       \param L4      [in] - Kokkos::DynRankView of the fourth barycentric coordinate of the input points
   */

  template<typename pointValueType, class ...pointProperties>
  static void
  warpShiftFace3D( Kokkos::DynRankView<pointValueType,pointProperties...> dxy,
                   const ordinal_type order ,
                   const pointValueType pval ,
                   const Kokkos::DynRankView<pointValueType,pointProperties...> L1,
                   const Kokkos::DynRankView<pointValueType,pointProperties...> L2,
                   const Kokkos::DynRankView<pointValueType,pointProperties...> L3,
                   const Kokkos::DynRankView<pointValueType,pointProperties...> L4
                   );

   /** \brief Used internally to evaluate the point shift for warp-blend points on faces of tets
       \param dxy     [out] - Kokkos::DynRankView containing the amount to shift each point in the x and y direction
       \param order   [in] - the order of lattice
       \param pval    [in] - the "alpha" term in the warping function
       \param L1      [in] - Kokkos::DynRankView of the first barycentric coordinate of the input points
       \param L2      [in] - Kokkos::DynRankView of the second barycentric coordinate of the input points
       \param L3      [in] - Kokkos::DynRankView of the third barycentric coordinate of the input points
   */

  template<typename pointValueType, class ...pointProperties>
  static void
  evalshift( Kokkos::DynRankView<pointValueType,pointProperties...> dxy ,
             const ordinal_type order ,
             const pointValueType pval ,
             const Kokkos::DynRankView<pointValueType,pointProperties...> L1 ,
             const Kokkos::DynRankView<pointValueType,pointProperties...> L2 ,
             const Kokkos::DynRankView<pointValueType,pointProperties...> L3
             );

   /** \brief Used internally to compute the warp on edges of a triangle in warp-blend points
       \param warp      [out] - a 1d Kokkos::DynRankView containing the amount to move each point
       \param order     [in]  - the order of the lattice
       \param xnodes    [in]  - Kokkos::DynRankView of the points to warp to, typically the Gauss-Lobatto points
       \param xout      [in]  - Kokkos::DynRankView of the equispaced points on the edge
   */
  template<typename pointValueType, class ...pointProperties>
  static void
  evalwarp( Kokkos::DynRankView<pointValueType,pointProperties...> warp ,
            const ordinal_type order ,
            const Kokkos::DynRankView<pointValueType,pointProperties...> xnodes ,
            const Kokkos::DynRankView<pointValueType,pointProperties...> xout );
};

} 

#include "Intrepid2_PointToolsDef.hpp"

#endif
