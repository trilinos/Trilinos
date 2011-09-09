// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_PointTools.hpp
    \brief  Header file for utility class to provide point tools,
            such as barycentric coordinates, equispaced lattices, and
            warp-blend point distrubtions.
    \author Created by R. Kirby
*/

#ifndef INTREPID_POINTTOOLS_HPP
#define INTREPID_POINTTOOLS_HPP

#include "Shards_CellTopology.hpp"
#include "Teuchos_TestForException.hpp"
#include "Intrepid_Polylib.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include <stdexcept>

namespace Intrepid {
  
  /** \class Intrepid::PointTools
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
   static inline int getLatticeSize( const shards::CellTopology& cellType ,
                                     const int order ,
                                     const int offset = 0 )
    {
      switch( cellType.getKey() ) {
      case shards::Tetrahedron<4>::key:
      case shards::Tetrahedron<8>::key:
      case shards::Tetrahedron<10>::key:
        {
          const int effectiveOrder = order - 4 * offset;
          if (effectiveOrder < 0) return 0;
          else return (effectiveOrder+1)*(effectiveOrder+2)*(effectiveOrder+3)/6;
        }
        break;
      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key: 
        {
          const int effectiveOrder = order - 3 * offset;
          if (effectiveOrder < 0) return 0;
          else return (effectiveOrder+1)*(effectiveOrder+2)/2;
        }
        break;
      case shards::Line<2>::key:
      case shards::Line<3>::key:
        {
          const int effectiveOrder = order - 2 * offset;
          if (effectiveOrder < 0) return 0;
          else return (effectiveOrder+1);
        }
        break;
      default:
        TEST_FOR_EXCEPTION( true , std::invalid_argument ,
                            ">>> ERROR (Intrepid::PointTools::getLatticeSize): Illegal cell type" );
      }
    }


    /** \brief Computes a lattice of points of a given
              order on a reference simplex (currently disabled for
              other cell types).  The output array is
              (P,D), where
        \code
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  pts         [out] - Output array of point coords
        \param  cellType    [in]  - type of reference cell (currently only supports the simplex)
        \param  order       [in]  - number of points per side, plus 1
        \param  pointType   [in]  - flag for point distribution.  Currently equispaced and
                                    warp/blend points are supported
        \param  offset      [in]  - Number of points on boundary to skip
    */
    template<class Scalar, class ArrayType>
    static void getLattice( ArrayType &pts ,
                            const shards::CellTopology& cellType ,
                            const int order ,
                            const int offset = 0 ,
                            const EPointType pointType = POINTTYPE_EQUISPACED );

    /** Retrieves the Gauss-Legendre points from PolyLib, but lets us
	do it in an arbitrary ArrayType.
	\param  pts         [out] - Output array of point coords (P,)
	\param  order       [out] - number of Gauss points - 1
    */
 
    template<class Scalar, class ArrayType>
    static void getGaussPoints( ArrayType &pts ,
				const int order );


  private:
    /** \brief Converts Cartesian coordinates to barycentric coordinates
              on a batch of triangles.  
              The input array cartValues is (C,P,2)
              The output array baryValues is (C,P,3).
              The input array vertices is (C,3,2), where
        \code
          C - num. integration domains
          P - number of points per cell
        \endcode

        \param  baryValues      [out] - Output array of barycentric coords
        \param  cartValues      [in]  - Input array of Cartesian coords
        \param  vertices        [out] - Vertices of each cell.

    */
    template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
    static void cartToBaryTriangle( ArrayTypeOut & baryValues ,
                                    const ArrayTypeIn1 & cartValues ,
                                    const ArrayTypeIn2 & vertices );


    /** \brief Converts barycentric coordinates to Cartesian coordinates
              on a batch of triangles.  
              The input array baryValues is (C,P,3)
              The output array cartValues is (C,P,2).
              The input array vertices is (C,3,2), where
        \code
          C - num. integration domains
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  baryValues      [out] - Output array of barycentric coords
        \param  cartValues      [in]  - Input array of Cartesian coords
        \param  vertices        [out] - Vertices of each cell.

    */
    template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
    static void baryToCartTriangle( ArrayTypeOut & cartValues ,
                                    const ArrayTypeIn1 & baryValues ,
                                    const ArrayTypeIn2 & vertices );


    /** \brief Converts Cartesian coordinates to barycentric coordinates
              on a batch of tetrahedra.  
              The input array cartValues is (C,P,3)
              The output array baryValues is (C,P,4).
              The input array vertices is (C,4,3), where
        \code
          C - num. integration domains
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  baryValues      [out] - Output array of barycentric coords
        \param  cartValues      [in]  - Input array of Cartesian coords
        \param  vertices        [out] - Vertices of each cell.

    */
    template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
    static void cartToBaryTetrahedron( ArrayTypeOut & baryValues ,
                                      const ArrayTypeIn1 & cartValues ,
                                      const ArrayTypeIn2 & vertices );

    /** \brief Converts barycentric coordinates to Cartesian coordinates
              on a batch of tetrahedra.  
              The input array baryValues is (C,P,4)
              The output array cartValues is (C,P,3).
              The input array vertices is (C,4,3), where
        \code
          C - num. integration domains
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  baryValues      [out] - Output array of barycentric coords
        \param  cartValues      [in]  - Input array of Cartesian coords
        \param  vertices        [out] - Vertices of each cell.

    */
    template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
    static void baryToCartTetrahedron( ArrayTypeOut & cartValues ,
                                      const ArrayTypeIn1 & baryValues ,
                                      const ArrayTypeIn2 & vertices );


    /** \brief Computes an equispaced lattice of a given
              order on a reference simplex (currently disabled for
              other cell types).  The output array is
              (P,D), where
        \code
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per side, plus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  cellType    [in]  - type of reference cell (currently only supports the simplex)

    */
  template<class Scalar, class ArrayType>
  static void getEquispacedLattice( const shards::CellTopology& cellType ,
                                    ArrayType &points ,
                                    const int order ,
                                    const int offset = 0 );


    /** \brief Computes a warped lattice (ala Warburton's warp-blend points of a given 
        order on a reference simplex (currently disabled for
        other cell types).  The output array is
        (P,D), where
        \code
        P - number of points per cell
        D - is the spatial dimension
        \endcode
        
        \param  points      [out] - Output array of point coords
        \param  order       [in]  - number of points per side, plus 1
        \param  offset      [in]  - Number of points on boundary to skip
        \param  cellType    [in]  - type of reference cell (currently only supports the simplex)
        
    */
    template<class Scalar, class ArrayType>
    static void getWarpBlendLattice( const shards::CellTopology& cellType ,
                                    ArrayType &points ,
                                    const int order ,
                                    const int offset = 0);
    

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
    template<class Scalar, class ArrayType>
    static void getEquispacedLatticeLine( ArrayType &points ,
                                          const int order ,
                                          const int offset = 0 );

    /** \brief Computes an equispaced lattice of a given
              order on the reference triangle.  The output array is
              (P,2), where
        \code
          P - number of points, which is 
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
                                    minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
                                    in from boundary

    */
    template<class Scalar, class ArrayType>
    static void getEquispacedLatticeTriangle( ArrayType &points ,
                                              const int order ,
                                              const int offset = 0 );

    /** \brief Computes an equispaced lattice of a given
              order on the reference tetrahedron.  The output array is
              (P,3), where
        \code
          P - number of points, which is 
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
                                    minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
                                    in from boundary

    */
    template<class Scalar, class ArrayType>
    static void getEquispacedLatticeTetrahedron( ArrayType &points ,
                                                const int order ,
                                                const int offset = 0 );

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
    template<class Scalar, class ArrayType>
    static void getWarpBlendLatticeLine( ArrayType &points ,
                                        const int order ,
                                        const int offset = 0 );

    /** \brief interpolates Warburton's warp function on the line
        \param  order       [in]  - The polynomial order
        \param  xnodes      [in]  - vector of node locations to interpolate
        \param  xout        [in]  - warpfunction at xout, +/- 1 roots deflated
        \param  warp        [out] - the amount to warp each point
    */
    template<class Scalar, class ArrayType>
    static void warpFactor( const int order ,
                            const ArrayType &xnodes ,
                            const ArrayType &xout ,
                            ArrayType &warp );

    /** \brief Returns Warburton's warp-blend points of a given
              order on the reference triangle.  The output array is
              (P,2), where
        \code
          P - number of points
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
                                    minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
                                    in per boundary

    */
    template<class Scalar, class ArrayType>
    static void getWarpBlendLatticeTriangle(ArrayType &points ,
                                            const int order ,
                                            const int offset = 0 );

    /** \brief Returns Warburton's warp-blend points of a given
              order on the reference tetrahedron.  The output array is
              (P,3), where
        \code
          P - number of points
        \endcode

        \param  points      [out] - Output array of point coords
        \param  order       [in]  - The lattice has order + 1 points,
                                    minus any skipped by offset
        \param  offset      [in]  - Number of points on boundary to skip coming
                                    in per boundary
    */
    template<class Scalar, class ArrayType>
    static void getWarpBlendLatticeTetrahedron( ArrayType &points ,
                                                const int order ,
                                                const int offset = 0 );


    /** \brief This is used internally to compute the tetrahedral warp-blend points one each face
        \param order   [in] - the order of lattice
        \param pval    [in] - the "alpha" term in the warping function
        \param L1      [in] - the first barycentric coordinate of the input points
        \param L2      [in] - the second barycentric coordinate of the input points
        \param L3      [in] - the third barycentric coordinate of the input points
        \param L4      [in] - the fourth barycentric coordinate of the input points
        \param dxy     [out] - contains the amount to shift each point in the x and y direction 
    */

  template<class Scalar, class ArrayType>
  static void warpShiftFace3D( const int order ,
                        const Scalar pval ,
                        const ArrayType &L1,
                        const ArrayType &L2,
                        const ArrayType &L3,
                        const ArrayType &L4,
                        ArrayType &dxy);

    /** \brief Used internally to evaluate the point shift for warp-blend points on faces of tets 
        \param order   [in] - the order of lattice
        \param pval    [in] - the "alpha" term in the warping function
        \param L1      [in] - the first barycentric coordinate of the input points
        \param L2      [in] - the second barycentric coordinate of the input points
        \param L3      [in] - the third barycentric coordinate of the input points
        \param dxy     [out] - contains the amount to shift each point in the x and y direction 
    */

  template<class Scalar, class ArrayType>
  static void evalshift( const int order ,
                  const Scalar pval ,
                  const ArrayType &L1 ,
                  const ArrayType &L2 ,
                  const ArrayType &L3 ,
                  ArrayType &dxy );

    /** \brief Used internally to compute the warp on edges of a triangle in warp-blend points
        \param warp      [out] - a 1d array containing the amount to move each point
        \param order     [in]  - the order of the lattice 
        \param xnodes    [in]  - the points to warp to, typically the Gauss-Lobatto points
        \param xout      [in]  - the equispaced points on the edge
    */
  template<class Scalar, class ArrayType>
  static void evalwarp( ArrayType &warp ,
                const int order ,
                const ArrayType &xnodes ,
                const ArrayType &xout );


  }; // end class PointTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_PointToolsDef.hpp>

#endif
