// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu
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

  */
  class PointTools {
  public:

    /** \brief Converts Cartesian coordinates to barycentric coordinates
	       on a batch of simplices (triangle or tetrahedron).  
               The input array cartValues is (C,P,D)
               The output array baryValues is (C,P,D+1).
               The input array vertices is (C,D+1,D), where
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
    static void cartToBary( ArrayTypeOut & baryValues ,
	   	            const ArrayTypeIn1 & cartValues ,
			    const ArrayTypeIn2 & vertices ,
                            const shards::CellTopology& cellType );

    /** \brief Converts barycentric coordinates to Cartesian coordinates
	       on a batch of triangles.  
               The input array baryValues is (C,P,D+1)
               The output array cartValues is (C,P,D).
               The input array vertices is (C,D+1,D), where
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
    static void baryToCart( ArrayTypeOut & cartValues ,
	   	            const ArrayTypeIn1 & baryValues ,
			    const ArrayTypeIn2 & vertices ,
                            const shards::CellTopology& cellType );


    /** \brief Computes the number of pointsin a lattice of a given order
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
    static int getLatticeSize( const shards::CellTopology& cellType ,
			       const int order ,
			       const int offset = 0 );

    /** \brief Computes a lattice of points of a given
               order on a reference simplex (currently disabled for
               other cell types).  The output array is
               (P,D), where
        \code
          P - number of points per cell
          D - is the spatial dimension
        \endcode

        \param  points      [out] - Output array of point coords
        \param  cellTYpe    [in]  - type of reference cell (currently only supports the simplex)
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
        \param  cellTYpe    [in]  - type of reference cell (currently only supports the simplex)

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
        \param  cellTYpe    [in]  - type of reference cell (currently only supports the simplex)
	
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
	\param  xnodes      [in] - vector of node locations to interpolate
	\param  xout        [out]  - warpfunction at xout, \pm 1 roots deflated

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
    static void getWarpBlendLatticeTriangle( ArrayType &points ,
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

    

  template<class Scalar, class ArrayType>
  static void warpShiftFace3D( const int order ,
			const Scalar pval ,
			const ArrayType &L1,
			const ArrayType &L2,
			const ArrayType &L3,
			const ArrayType &L4,
			ArrayType &dxy);

  template<class Scalar, class ArrayType>
  static void evalshift( const int order ,
		  const Scalar pval ,
		  const ArrayType &L1 ,
		  const ArrayType &L2 ,
		  const ArrayType &l3 ,
		  ArrayType &dxy );

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
