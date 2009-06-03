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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_PointToolsDef.hpp
    \brief  Definition file for utilities for barycentric coordinates and lattices
    \author Created by R. Kirby
*/


namespace Intrepid {

  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::cartToBary( ArrayTypeOut & baryValues ,
				      const ArrayTypeIn1 & cartValues ,
				      const ArrayTypeIn2 & vertices ,
				      const shards::CellTopology& cellType )
  {
    switch (cellType.getKey()) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
      cartToBaryTetrahedron( baryValues ,
			     cartValues ,
			     vertices );
      break;
    case shards::Triangle<3>:
    case shards::Triangle<4>:
    case shards::Triangle<6>:
      cartToBaryTriangle( baryValues ,
			  cartValues ,
			  vertices );
    break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  ">>> ERROR (Intrepid::PointTools::cartToBary): Illegal cell type" );
    }
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::baryToCart( ArrayTypeOut & cartValues ,
				      const ArrayTypeIn1 & baryValues ,
				      const ArrayTypeIn2 & vertices ,
				      const shards::CellTopology& cellType )
  {
    switch (cellType.getKey()) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
      baryToCartTetrahedron( baryValues ,
			  cartValues ,
			  vertices );
      break;
    case shards::Triangle<3>:
    case shards::Triangle<4>:
    case shards::Triangle<6>:
      baryToCartTriangle( baryValues ,
			  cartValues ,
			  vertices );
    break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  ">>> ERROR (Intrepid::PointTools::baryToCart): Illegal cell type" );
    }
  }


  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::cartToBaryTriangle( ArrayTypeOut & baryValues ,
					      const ArrayTypeIn1 & cartValues ,
					      const ArrayTypeIn2 & vertices )
  {
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::cartToBaryTetrahedron( ArrayTypeOut & baryValues ,
						 const ArrayTypeIn1 & cartValues ,
						 const ArrayTypeIn2 & vertices )
  {
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::baryToCartTriangle( ArrayTypeOut & baryValues ,
					      const ArrayTypeIn1 & cartValues ,
					      const ArrayTypeIn2 & vertices )
  {
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
  static void PointTools::baryToCartTetrahedron( ArrayTypeOut & baryValues ,
						 const ArrayTypeIn1 & cartValues ,
						 const ArrayTypeIn2 & vertices )
  {
  }

  static int PointTools::getLatticeSize( const shards::CellTopology& cellType ,
					 const int order ,
					 const int offset = 0 )
  {
    switch( cellType.getKey() ) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
      const int effectiveOrder = order - 4 * offset;
      if (effectiveOrder < 0) return 0;
      else return (effectiveOrder+1)*(effectiveOrder+2)*(effectiveOrder+3)/6;
      break;
    case shards::Triangle<3>::key:
    case shards::Triangle<4>::key:
    case shards::Triangle<6>::key:
      const int effectiveOrder = order - 3 * offset;
      if (effectiveOrder < 0) return 0;
      else return (effectiveOrder+1)*(effectiveOrder+2)/2;
      break;
    case shards::Line<2>::key:
    case shards::Line<3>::key:
      const int effectiveOrder = order - 2 * offset;
      if (effectiveOrder < 0) return 0;
      else return (effectiveOrder+1);
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  ">>> ERROR (Intrepid::PointTools::getLatticeSize): Illegal cell type" );

    }
  }

  
  template<class Scalar, class ArrayType>
  static void PointTools::getEquispacedLattice( ArrayType &points ,
						const int order ,
						const int offset ,
						const shards::CellTopology& cellType )
  {
    switch (cellType.getKey()) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
      getEquispacedLatticeTetrahedron( points , order , offset );
      break;
    case shards::Triangle<3>::key:
    case shards::Triangle<4>::key:
    case shards::Triangle<6>::key:
      getEquispacedLatticeTriangle( points , order , offset );
      break;
    case shards::Line<2>::key:
    case shards::Line<3>::key:
      getEquispacedLatticeLine( points , order , offset );
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  ">>> ERROR (Intrepid::PointTools::getEquispacedLattice): Illegal cell type" );
    }
     
   }

   template<class Scalar, class ArrayType>
   static void PointTools::getWarpBlendLattice( ArrayType &points ,
						 const int order ,
						 const int offset ,
						 const shards::CellTopology& cellType )
   {
    switch (cellType.getKey()) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
      getWarpBlendLatticeTetrahedron( points , order , offset );
      break;
    case shards::Triangle<3>::key:
    case shards::Triangle<4>::key:
    case shards::Triangle<6>::key:
      getWarpBlendLatticeTriangle( points , order , offset );
      break;
    case shards::Line<2>::key:
    case shards::Line<3>::key:
      getWarpBlendLatticeLine( points , order , offset );
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  ">>> ERROR (Intrepid::PointTools::getWarpBlendLattice): Illegal cell type" );
    }
     
   }

   template<class Scalar, class ArrayType>
   static void PointTools::getEquispacedLatticeLine( ArrayType &points ,
						     const int order ,
						     const int offset )
   {
     TEST_FOR_EXCEPTION( order <= 0 ,
			 std::invalid_argument ,
			 ">>> ERROR (Intrepid::PointTools::getEquispacedLatticeLine): order must be positive" );
     const Scalar h = 2 / order;

     for (int i=offset;i<=order-offset;i++) {
       points(i-offset,0) = -1.0 + h * (Scalar) i;
     }

     return;
   }

   template<class Scalar, class ArrayType>
   static void PointTools::getEquispacedLatticeTriangle( ArrayType &points ,
							 const int order ,
							 const int offset )
   {
     TEST_FOR_EXCEPTION( order <= 0 ,
			 std::invalid_argument ,
			 ">>> ERROR (Intrepid::PointTools::getEquispacedLatticeLine): order must be positive" );

     const Scalar h = 1.0 / order;
     int cur = 0;

     for (int i=interior;i<=n-interior;i++) {
       for (int j=interior;j<=n-ii-interior;j++) {
	 points(cur,0) = (Scalar)0.0 + (Scalar) j * h ;
	 points(cur,1) = (Scalar)0.0 + (Scalar) i * h;
	 cur++;
       }
     }

     return;
   }

   template<class Scalar, class ArrayType>
   static void PointTools::getEquispacedLatticeTetrahedron( ArrayType &points ,
							    const int order ,
							    const int offset )
   {
     TEST_FOR_EXCEPTION( order <= 0 ,
			 std::invalid_argument ,
			 ">>> ERROR (Intrepid::PointTools::getEquispacedLatticeLine): order must be positive" );

     const Scalar h = 1.0 / order;
     int cur = 0;

     for (int i=interior;i<=n-interior;i++) {
       for (int j=interior;j<=n-ii-interior;j++) {
	 for (int k=interior;k<=n-ii-jj-interior;j++) {
	   points(cur,0) = (Scalar) k * h;
	   points(cur,1) = (Scalar) j * h;
	   points(cur,2) = (Scalar) i * h;
	   cur++;
	 }
       }
     }

     return;
   }





} // namespace Intrepid
