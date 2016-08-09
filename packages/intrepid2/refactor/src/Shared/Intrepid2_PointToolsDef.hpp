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

/** \file   Intrepid_PointToolsDef.hpp
    \brief  Definition file for utilities for barycentric coordinates and lattices
    \author Created by R. Kirby
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_POINTTOOLS_DEF_HPP__
#define __INTREPID2_POINTTOOLS_DEF_HPP__

#ifdef _MSC_VER
#include "winmath.h"
#endif


namespace Intrepid2 {

  // -----------------------------------------------------------------------------------------  
  // Front interface
  // -----------------------------------------------------------------------------------------
  
  inline
  ordinal_type
  PointTools::
  getLatticeSize( const shards::CellTopology cellType,
                  const ordinal_type order,
                  const ordinal_type offset ) {
#ifdef HAVE_INTREPID2_DEBUG    
    INTREPID2_TEST_FOR_EXCEPTION( order < 0 || offset < 0,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getLatticeSize): order and offset must be positive values." );
#endif
    ordinal_type r_val = 0;
    switch (cellType.getKey()) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key: {
      const auto effectiveOrder = order - 4 * offset;
      r_val = (effectiveOrder < 0 ? 0 :(effectiveOrder+1)*(effectiveOrder+2)*(effectiveOrder+3)/6);
      break;
    }
    case shards::Triangle<3>::key:
    case shards::Triangle<4>::key:
    case shards::Triangle<6>::key: {
      const auto effectiveOrder = order - 3 * offset;
      r_val = (effectiveOrder < 0 ? 0 : (effectiveOrder+1)*(effectiveOrder+2)/2);
      break;
    }
    case shards::Line<2>::key:
    case shards::Line<3>::key: {
      const auto effectiveOrder = order - 2 * offset;
      r_val = (effectiveOrder < 0 ? 0 : (effectiveOrder+1));
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
                                    ">>> ERROR (Intrepid2::PointTools::getLatticeSize): the specified cell type is not supported." );
    }
    }
    return r_val;
  }

  template<typename pointValueType, class ...pointProperties>
  void 
  PointTools::
  getLattice( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
              const shards::CellTopology cell,
              const ordinal_type         order,
              const ordinal_type         offset,
              const EPointType           pointType ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( points.rank() != 2,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getLattice): points rank must be 2." );
    INTREPID2_TEST_FOR_EXCEPTION( order < 0 || offset < 0,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getLattice): order and offset must be positive values." );

    const auto latticeSize = getLatticeSize( cell, order, offset );
    const auto spaceDim = cell.getDimension();
    
    INTREPID2_TEST_FOR_EXCEPTION( points.dimension(0) != latticeSize ||
                                  points.dimension(1) != spaceDim,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getLattice): dimension does not match to lattice size." );
#endif

    // const auto latticeSize = getLatticeSize( cell, order, offset );
    // const auto spaceDim = cell.getDimension();
    
    // // the interface assumes that the input array follows the cell definition
    // // so, let's match all dimensions according to the cell specification
    // typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
    // auto pts = Kokkos::subdynrankview( points, 
    //                                    range_type(0, latticeSize), 
    //                                    range_type(0, spaceDim) );   
    switch (pointType) {
    case POINTTYPE_EQUISPACED:  getEquispacedLattice( points, cell, order, offset ); break;
    case POINTTYPE_WARPBLEND:   getWarpBlendLattice ( points, cell, order, offset ); break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true ,
                                    std::invalid_argument ,
                                    ">>> ERROR (PointTools::getLattice): invalid EPointType." );
    }
    }
  }

  template<typename pointValueType, class ...pointProperties>
  void 
  PointTools::
  getGaussPoints( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
                  const ordinal_type order ) {
#ifdef HAVE_INTREPID2_DEBUG    
    INTREPID2_TEST_FOR_EXCEPTION( points.rank() != 2,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getGaussPoints): points rank must be 1." );
    INTREPID2_TEST_FOR_EXCEPTION( order < 0,
                                  std::invalid_argument ,
                                  ">>> ERROR (PointTools::getGaussPoints): order must be positive value." );
#endif
    const ordinal_type np = order + 1;
    const double alpha = 0.0, beta = 0.0;
    
    // until view and dynrankview inter-operatible, we use views in a consistent way
    Kokkos::View<pointValueType*,Kokkos::HostSpace> 
      zHost("PointTools::getGaussPoints::z", np), 
      wHost("PointTools::getGaussPoints::w", np);
    
    // sequential means that the code is decorated with KOKKOS_INLINE_FUNCTION 
    // and it does not invoke parallel for inside (cheap operation), which means 
    // that gpu memory is not accessible unless this is called inside of functor.
    Polylib::Serial::Cubature<POLYTYPE_GAUSS>::getValues(zHost, wHost, np, alpha, beta);

    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
    auto pts = Kokkos::subdynrankview( points, range_type(0,np), 0 );
    // should be fixed after view and dynrankview are inter-operatible
    auto z   = Kokkos::DynRankView<pointValueType,Kokkos::HostSpace>(zHost.data(), np);
    Kokkos::deep_copy(pts, z);
  }

  // -----------------------------------------------------------------------------------------  
  // Internal implementation
  // -----------------------------------------------------------------------------------------
  
  template<typename pointValueType, class ...pointProperties>
  void PointTools::
  getEquispacedLattice( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
                        const shards::CellTopology cell,
                        const ordinal_type order,
                        const ordinal_type offset ) {
    switch (cell.getKey()) {
    // case shards::Tetrahedron<4>::key:
    // case shards::Tetrahedron<8>::key:
    // case shards::Tetrahedron<10>::key: getEquispacedLatticeTetrahedron( points, order, offset );  break;
    // case shards::Triangle<3>::key:
    // case shards::Triangle<4>::key:
    // case shards::Triangle<6>::key:     getEquispacedLatticeTriangle   ( points, order, offset );  break;
    case shards::Line<2>::key:
    case shards::Line<3>::key:         getEquispacedLatticeLine       ( points, order, offset );  break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
                                    ">>> ERROR (Intrepid2::PointTools::getEquispacedLattice): the specified cell type is not supported." );
    }
    }
    
  }

  template<typename pointValueType, class ...pointProperties>
  void PointTools::
  getWarpBlendLattice( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
                       const shards::CellTopology cell,
                       const ordinal_type order,
                       const ordinal_type offset ) {
    switch (cell.getKey()) {
    // case shards::Tetrahedron<4>::key:
    // case shards::Tetrahedron<8>::key:
    // case shards::Tetrahedron<10>::key: getWarpBlendLatticeTetrahedron( points, order, offset );  break;
    // case shards::Triangle<3>::key:
    // case shards::Triangle<4>::key:
    // case shards::Triangle<6>::key:     getWarpBlendLatticeTriangle   ( points, order, offset );  break;
    case shards::Line<2>::key:
    case shards::Line<3>::key:         getWarpBlendLatticeLine       ( points, order, offset );  break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
                                    ">>> ERROR (Intrepid2::PointTools::getWarpBlendLattice): the specified cell type is not supported." );
    }
    }
  }

  // -----------------------------------------------------------------------------------------

  template<typename pointValueType, class ...pointProperties>
  void
  PointTools::
  getEquispacedLatticeLine( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
                            const ordinal_type order,
                            const ordinal_type offset ) {
    auto pointsHost = Kokkos::create_mirror_view(Kokkos::HostSpace::memory_space(), points);
    
    if (order == 0) 
      pointsHost(0,0) = 0.0;
    else {
      const pointValueType h = 2.0 / order;
      const ordinal_type ibeg = offset, iend = order-offset+1;
      for (auto i=ibeg;i<iend;++i) 
	pointsHost(i-ibeg, 0) = -1.0 + h * (pointValueType) i;
    }

    Kokkos::deep_copy(points, pointsHost);
  }

  template<typename pointValueType, class ...pointProperties>
  void 
  PointTools::
  getWarpBlendLatticeLine( /**/  Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const ordinal_type order,
                           const ordinal_type offset ) {
    // order is order of polynomial degree.  The Gauss-Lobatto points are accurate
    // up to degree 2*i-1
    const ordinal_type np = order + 1;
    const double alpha = 0.0, beta = 0.0;
    const ordinal_type s = np - 2*offset;
    
    if (s > 0) {
      // until view and dynrankview inter-operatible, we use views in a consistent way
      Kokkos::View<pointValueType*,Kokkos::HostSpace> 
        zHost("PointTools::getGaussPoints::z", np), 
        wHost("PointTools::getGaussPoints::w", np);
      
      // sequential means that the code is decorated with KOKKOS_INLINE_FUNCTION 
      // and it does not invoke parallel for inside (cheap operation), which means 
      // that gpu memory is not accessible unless this is called inside of functor.
      Polylib::Serial::Cubature<POLYTYPE_GAUSS_LOBATTO>::getValues(zHost, wHost, np, alpha, beta);
      
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      auto pts = Kokkos::subdynrankview( points, range_type(0, s), 0 );
      
      // this should be fixed after view and dynrankview is interoperatable
      auto z   = Kokkos::DynRankView<pointValueType,Kokkos::HostSpace>(zHost.data() + offset, np-offset);
      
      Kokkos::deep_copy(pts, z);
    } 
  }

  // template<class Scalar, class ArrayType>
  // void PointTools::getEquispacedLatticeTriangle( ArrayType &points ,
  //                                                const int order ,
  //                                                const int offset )
  // {
  //   TEUCHOS_TEST_FOR_EXCEPTION( order <= 0 ,
  //                       std::invalid_argument ,
  //                       ">>> ERROR (Intrepid2::PointTools::getEquispacedLatticeLine): order must be positive" );

  //   const Scalar h = 1.0 / order;
  //   int cur = 0;

  //   for (int i=offset;i<=order-offset;i++) {
  //     for (int j=offset;j<=order-i-offset;j++) {
  //       points(cur,0) = (Scalar)0.0 + (Scalar) j * h ;
  //       points(cur,1) = (Scalar)0.0 + (Scalar) i * h;
  //       cur++;
  //     }
  //   }

  //   return;
  // }

  // template<class Scalar, class ArrayType>
  // void PointTools::getEquispacedLatticeTetrahedron( ArrayType &points ,
  //                                                   const int order ,
  //                                                   const int offset )
  // {
  //   TEUCHOS_TEST_FOR_EXCEPTION( (order <= 0) ,
  //                       std::invalid_argument ,
  //                       ">>> ERROR (Intrepid2::PointTools::getEquispacedLatticeTetrahedron): order must be positive" );

  //   const Scalar h = 1.0 / order;
  //   int cur = 0;

  //   for (int i=offset;i<=order-offset;i++) {
  //     for (int j=offset;j<=order-i-offset;j++) {
  //       for (int k=offset;k<=order-i-j-offset;k++) {
  //         points(cur,0) = (Scalar) k * h;
  //         points(cur,1) = (Scalar) j * h;
  //         points(cur,2) = (Scalar) i * h;
  //         cur++;
  //       }
  //     }
  //   }

  //   return;
  // }


  // template<class Scalar, class ArrayType>
  // void PointTools::warpFactor( const int order , 
  //                             const ArrayType &xnodes ,
  //                             const ArrayType &xout ,
  //                             ArrayType &warp)
  // {
  //   TEUCHOS_TEST_FOR_EXCEPTION( ( warp.dimension(0) != xout.dimension(0) ) ,
  //                       std::invalid_argument ,
  //                       ">>> ERROR (PointTools::warpFactor): xout and warp must be same size." );

  //   warp.initialize();

  //   FieldContainer<Scalar> d( xout.dimension(0) );
  //   d.initialize();

  //   FieldContainer<Scalar> xeq( order + 1 ,1);
  //   PointTools::getEquispacedLatticeLine<Scalar,ArrayType>( xeq , order , 0 );
  //   xeq.resize( order + 1 );

  //   TEUCHOS_TEST_FOR_EXCEPTION( ( xeq.dimension(0) != xnodes.dimension(0) ) ,
  //                       std::invalid_argument ,
  //                       ">>> ERROR (PointTools::warpFactor): xeq and xnodes must be same size." );
    
  //   for (int i=0;i<=order;i++) {

  //     for (int k=0;k<d.dimension(0);k++) {
  //       d(k) = xnodes(i) - xeq(i);
  //     }

  //     for (int j=1;j<order;j++) {
  //       if (i != j) {
  //         for (int k=0;k<d.dimension(0);k++) {
  //           d(k) = d(k) * ( (xout(k)-xeq(j)) / (xeq(i)-xeq(j)) );
  //         }
  //       }
  //     }
      
  //     // deflate end roots
  //     if ( i != 0 ) {
  //       for (int k=0;k<d.dimension(0);k++) {
  //         d(k) = -d(k) / (xeq(i) - xeq(0));
  //       }
  //     }

  //     if (i != order ) {
  //       for (int k=0;k<d.dimension(0);k++) {
  //         d(k) = d(k) / (xeq(i) - xeq(order));
  //       }
  //     }

  //     for (int k=0;k<d.dimension(0);k++) {
  //       warp(k) += d(k);
  //     }

  //   }


  //   return;
  // }

  // template<class Scalar, class ArrayType>
  // void PointTools::getWarpBlendLatticeTriangle( ArrayType &points ,
  //                                               const int order ,
  //                                               const int offset  )
  // {
  //   /* get Gauss-Lobatto points */

  //   Intrepid2::FieldContainer<Scalar> gaussX( order + 1 , 1 );
    
  //   PointTools::getWarpBlendLatticeLine<Scalar,FieldContainer<Scalar> >( gaussX , order , 0 );
    
  //   gaussX.resize(gaussX.dimension(0));

  //   Scalar alpopt[] = {0.0000,0.0000,1.4152,0.1001,0.2751,0.9800,1.0999,
  //                       1.2832,1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};

  //   Scalar alpha;

  //   if (order >= 1 && order < 16) {
  //     alpha = alpopt[order-1];
  //   }
  //   else {
  //     alpha = 5.0 / 3.0;
  //   }

  //   const int p = order; /* switch to Warburton's notation */
  //   int N = (p+1)*(p+2)/2;
    
  //   /* equidistributed nodes on equilateral triangle */
  //   Intrepid2::FieldContainer<Scalar> L1( N );
  //   Intrepid2::FieldContainer<Scalar> L2( N );
  //   Intrepid2::FieldContainer<Scalar> L3( N );     
  //   Intrepid2::FieldContainer<Scalar> X(N);
  //   Intrepid2::FieldContainer<Scalar> Y(N);

  //   int sk = 0;
  //   for (int n=1;n<=p+1;n++) {
  //     for (int m=1;m<=p+2-n;m++) {
  //       L1(sk) = (n-1) / (Scalar)p;
  //       L3(sk) = (m-1) / (Scalar)p;
  //       L2(sk) = 1.0 - L1(sk) - L3(sk);
  //       sk++;
  //     }
  //   }
    
  //   for (int n=0;n<N;n++) {
  //     X(n) = -L2(n) + L3(n);
  //     Y(n) = (-L2(n) - L3(n) + 2*L1(n))/1.7320508075688772;
  //   }

  //   /* get blending function for each node at each edge */
  //   Intrepid2::FieldContainer<Scalar> blend1(N);
  //   Intrepid2::FieldContainer<Scalar> blend2(N);
  //   Intrepid2::FieldContainer<Scalar> blend3(N);
    
  //   for (int n=0;n<N;n++) {
  //     blend1(n) = 4.0 * L2(n) * L3(n);
  //     blend2(n) = 4.0 * L1(n) * L3(n);
  //     blend3(n) = 4.0 * L1(n) * L2(n);
  //   }
    
  //   /* get difference of each barycentric coordinate */
  //   Intrepid2::FieldContainer<Scalar> L3mL2(N);
  //   Intrepid2::FieldContainer<Scalar> L1mL3(N);
  //   Intrepid2::FieldContainer<Scalar> L2mL1(N);

  //   for (int k=0;k<N;k++) {
  //     L3mL2(k) = L3(k)-L2(k);
  //     L1mL3(k) = L1(k)-L3(k);
  //     L2mL1(k) = L2(k)-L1(k);
  //   }

  //   FieldContainer<Scalar> warpfactor1(N);
  //   FieldContainer<Scalar> warpfactor2(N);
  //   FieldContainer<Scalar> warpfactor3(N);
    
  //   warpFactor<Scalar,FieldContainer<Scalar> >( order , gaussX , L3mL2 , warpfactor1 );
  //   warpFactor<Scalar,FieldContainer<Scalar> >( order , gaussX , L1mL3 , warpfactor2 );
  //   warpFactor<Scalar,FieldContainer<Scalar> >( order , gaussX , L2mL1 , warpfactor3 );

  //   FieldContainer<Scalar> warp1(N);
  //   FieldContainer<Scalar> warp2(N);
  //   FieldContainer<Scalar> warp3(N);

  //   for (int k=0;k<N;k++) {
  //     warp1(k) = blend1(k) * warpfactor1(k) *
  //       ( 1.0 + alpha * alpha * L1(k) * L1(k) );
  //     warp2(k) = blend2(k) * warpfactor2(k) *
  //       ( 1.0 + alpha * alpha * L2(k) * L2(k) );
  //     warp3(k) = blend3(k) * warpfactor3(k) *
  //       ( 1.0 + alpha * alpha * L3(k) * L3(k) );
  //   }

  //   for (int k=0;k<N;k++) {
  //     X(k) += 1.0 * warp1(k) + cos( 2.0 * M_PI / 3.0 ) * warp2(k) + cos(4*M_PI/3.0) * warp3(k);
  //     Y(k) += 0.0 * warp1(k) + sin( 2.0 * M_PI / 3.0 ) * warp2(k) + sin( 4*M_PI/3.0) * warp3(k);
  //   }

  //   FieldContainer<Scalar> warXY(N,2);
    
  //   for (int k=0;k<N;k++) {
  //     warXY(k,0) = X(k);
  //     warXY(k,1) = Y(k);
  //   }


  //   // finally, convert the warp-blend points to the correct triangle
  //   FieldContainer<Scalar> warburtonVerts(1,3,2);
  //   warburtonVerts(0,0,0) = -1.0;
  //   warburtonVerts(0,0,1) = -1.0/sqrt(3.0);
  //   warburtonVerts(0,1,0) = 1.0;
  //   warburtonVerts(0,1,1) = -1.0/sqrt(3.0);
  //   warburtonVerts(0,2,0) = 0.0;
  //   warburtonVerts(0,2,1) = 2.0/sqrt(3.0);

  //   FieldContainer<Scalar> refPts(N,2);

  //   Intrepid2::CellTools<Scalar>::mapToReferenceFrame( refPts ,
  //                                                     warXY ,
  //                                                     warburtonVerts ,
  //                                                     shards::getCellTopologyData< shards::Triangle<3> >(),
  //                                                     0 );

  //   // now write from refPts into points, taking care of offset
  //   int noffcur = 0;  // index into refPts
  //   int offcur = 0;   // index int points
  //   for (int i=0;i<=order;i++) {
  //     for (int j=0;j<=order-i;j++) {
  //       if ( (i >= offset) && (i <= order-offset) &&
  //             (j >= offset) && (j <= order-i-offset) ) {
  //         points(offcur,0) = refPts(noffcur,0);
  //         points(offcur,1) = refPts(noffcur,1);
  //         offcur++;
  //       }
  //       noffcur++;
  //     }
  //   }

  //   return;
  // }
  

  // template<class Scalar, class ArrayType>
  // void PointTools::warpShiftFace3D( const int order ,
  //                                   const Scalar pval ,
  //                                   const ArrayType &L1,
  //                                   const ArrayType &L2,
  //                                   const ArrayType &L3,
  //                                   const ArrayType &L4,
  //                                   ArrayType &dxy)
  // {
  //   evalshift<Scalar,ArrayType>(order,pval,L2,L3,L4,dxy);
  //   return;
  // }

  // template<class Scalar, class ArrayType>
  // void PointTools::evalshift( const int order ,
  //                             const Scalar pval ,
  //                             const ArrayType &L1 ,
  //                             const ArrayType &L2 ,
  //                             const ArrayType &L3 ,
  //                             ArrayType &dxy )
  // {
  //   // get Gauss-Lobatto-nodes
  //   FieldContainer<Scalar> gaussX(order+1,1);
  //   PointTools::getWarpBlendLatticeLine<Scalar,FieldContainer<Scalar> >( gaussX , order , 0 );
  //   gaussX.resize(order+1);
  //   const int N = L1.dimension(0);
    
  //   // Warburton code reverses them
  //   for (int k=0;k<=order;k++) {
  //     gaussX(k) *= -1.0;
  //   }

  //   // blending function at each node for each edge
  //   FieldContainer<Scalar> blend1(N);
  //   FieldContainer<Scalar> blend2(N);
  //   FieldContainer<Scalar> blend3(N);

  //   for (int i=0;i<N;i++) {
  //     blend1(i) = L2(i) * L3(i);
  //     blend2(i) = L1(i) * L3(i);
  //     blend3(i) = L1(i) * L2(i);
  //   }

  //   // amount of warp for each node for each edge
  //   FieldContainer<Scalar> warpfactor1(N);
  //   FieldContainer<Scalar> warpfactor2(N);
  //   FieldContainer<Scalar> warpfactor3(N);

  //   // differences of barycentric coordinates 
  //   FieldContainer<Scalar> L3mL2(N);
  //   FieldContainer<Scalar> L1mL3(N);
  //   FieldContainer<Scalar> L2mL1(N);
    
  //   for (int k=0;k<N;k++) {
  //     L3mL2(k) = L3(k)-L2(k);
  //     L1mL3(k) = L1(k)-L3(k);
  //     L2mL1(k) = L2(k)-L1(k);
  //   }
    
  //   evalwarp<Scalar,FieldContainer<Scalar> >( warpfactor1 , order , gaussX , L3mL2 );
  //   evalwarp<Scalar,FieldContainer<Scalar> >( warpfactor2 , order , gaussX , L1mL3 );
  //   evalwarp<Scalar,FieldContainer<Scalar> >( warpfactor3 , order , gaussX , L2mL1 );
    
  //   for (int k=0;k<N;k++) {
  //     warpfactor1(k) *= 4.0;
  //     warpfactor2(k) *= 4.0;
  //     warpfactor3(k) *= 4.0;      
  //   }

  //   FieldContainer<Scalar> warp1(N);
  //   FieldContainer<Scalar> warp2(N);
  //   FieldContainer<Scalar> warp3(N);
    
  //   for (int k=0;k<N;k++) {
  //     warp1(k) = blend1(k) * warpfactor1(k) *
  //       ( 1.0 + pval * pval * L1(k) * L1(k) );
  //     warp2(k) = blend2(k) * warpfactor2(k) *
  //       ( 1.0 + pval * pval * L2(k) * L2(k) );
  //     warp3(k) = blend3(k) * warpfactor3(k) *
  //       ( 1.0 + pval * pval * L3(k) * L3(k) );
  //   }

  //   for (int k=0;k<N;k++) {
  //     dxy(k,0) = 1.0 * warp1(k) + cos( 2.0 * M_PI / 3.0 ) * warp2(k) + cos( 4.0*M_PI/3.0 ) * warp3(k);
  //     dxy(k,1) = 0.0 * warp1(k) + sin( 2.0 * M_PI / 3.0 ) * warp2(k) + sin( 4.0*M_PI/3.0 ) * warp3(k);
  //   }

  //   return;

  // }

  // /* one-d edge warping function */
  // template<class Scalar, class ArrayType>
  // void PointTools::evalwarp( ArrayType &warp ,
  //                           const int order ,
  //                           const ArrayType &xnodes ,
  //                           const ArrayType &xout )
  // {
  //   FieldContainer<Scalar> xeq(order+1);
  //   FieldContainer<Scalar> d(xout.dimension(0));

  //   d.initialize();

  //   for (int i=0;i<=order;i++) {
  //     xeq(i) = -1.0 + 2.0 * ( order - i ) / order;
  //   }



  //   for (int i=0;i<=order;i++) {
  //     d.initialize( xnodes(i) - xeq(i) );
  //     for (int j=1;j<order;j++) {
  //       if (i!=j) {
  //         for (int k=0;k<d.dimension(0);k++) {
  //           d(k) = d(k) * (xout(k)-xeq(j))/(xeq(i)-xeq(j));
  //         }
  //       }
  //     }
  //     if (i!=0) {
  //       for (int k=0;k<d.dimension(0);k++) {
  //         d(k) = -d(k)/(xeq(i)-xeq(0));
  //       }
  //     }
  //     if (i!=order) {
  //       for (int k=0;k<d.dimension(0);k++) {
  //         d(k) = d(k)/(xeq(i)-xeq(order));
  //       }
  //     }
      
  //     for (int k=0;k<d.dimension(0);k++) {
  //       warp(k) += d(k);
  //     } 
  //   }    

  //   return;
  // }


  // template<class Scalar, class ArrayType>
  // void PointTools::getWarpBlendLatticeTetrahedron(ArrayType &points ,
  //                                                 const int order ,
  //                                                 const int offset  )
  // {
  //   Scalar alphastore[] = { 0,0,0,0.1002, 1.1332,1.5608,1.3413,1.2577,1.1603,
  //                           1.10153,0.6080,0.4523,0.8856,0.8717,0.9655};
  //   Scalar alpha;

  //   if (order <= 15) {
  //     alpha = alphastore[order-1]; 
  //   }
  //   else {
  //     alpha = 1.0;
  //   }

  //   const int N = (order+1)*(order+2)*(order+3)/6;
  //   Scalar tol = 1.e-10;

  //   FieldContainer<Scalar> shift(N,3);
  //   shift.initialize();

  //   /* create 3d equidistributed nodes on Warburton tet */
  //   FieldContainer<Scalar> equipoints(N,3);
  //   int sk = 0;
  //   for (int n=0;n<=order;n++) {
  //     for (int m=0;m<=order-n;m++) {
  //       for (int q=0;q<=order-n-m;q++) {
  //         equipoints(sk,0) = -1.0 + (q * 2.0 ) / order;
  //         equipoints(sk,1) = -1.0 + (m * 2.0 ) / order;
  //         equipoints(sk,2) = -1.0 + (n * 2.0 ) / order;
  //         sk++;
  //       }
  //     }
  //   }
    

  //   /* construct barycentric coordinates for equispaced lattice */
  //   FieldContainer<Scalar> L1(N);
  //   FieldContainer<Scalar> L2(N);
  //   FieldContainer<Scalar> L3(N);
  //   FieldContainer<Scalar> L4(N);
  //   for (int i=0;i<N;i++) {
  //     L1(i) = (1.0 + equipoints(i,2)) / 2.0;
  //     L2(i) = (1.0 + equipoints(i,1)) / 2.0;
  //     L3(i) = -(1.0 + equipoints(i,0) + equipoints(i,1) + equipoints(i,2)) / 2.0;
  //     L4(i) = (1.0 + equipoints(i,0)) / 2.0;
  //   }
    
    
  //   /* vertices of Warburton tet */
  //   FieldContainer<Scalar> warVerts(4,3);
  //   warVerts(0,0) = -1.0;
  //   warVerts(0,1) = -1.0/sqrt(3.0);
  //   warVerts(0,2) = -1.0/sqrt(6.0);
  //   warVerts(1,0) = 1.0;
  //   warVerts(1,1) = -1.0/sqrt(3.0);
  //   warVerts(1,2) = -1.0/sqrt(6.0);
  //   warVerts(2,0) = 0.0;
  //   warVerts(2,1) = 2.0 / sqrt(3.0);
  //   warVerts(2,2) = -1.0/sqrt(6.0);
  //   warVerts(3,0) = 0.0;
  //   warVerts(3,1) = 0.0;
  //   warVerts(3,2) = 3.0 / sqrt(6.0);


  //   /* tangents to faces */
  //   FieldContainer<Scalar> t1(4,3);
  //   FieldContainer<Scalar> t2(4,3);
  //   for (int i=0;i<3;i++) {
  //     t1(0,i) = warVerts(1,i) - warVerts(0,i);
  //     t1(1,i) = warVerts(1,i) - warVerts(0,i);
  //     t1(2,i) = warVerts(2,i) - warVerts(1,i);
  //     t1(3,i) = warVerts(2,i) - warVerts(0,i);
  //     t2(0,i) = warVerts(2,i) - 0.5 * ( warVerts(0,i) + warVerts(1,i) );
  //     t2(1,i) = warVerts(3,i) - 0.5 * ( warVerts(0,i) + warVerts(1,i) );
  //     t2(2,i) = warVerts(3,i) - 0.5 * ( warVerts(1,i) + warVerts(2,i) );
  //     t2(3,i) = warVerts(3,i) - 0.5 * ( warVerts(0,i) + warVerts(2,i) );
  //   }

  //   /* normalize tangents */
  //   for (int n=0;n<4;n++) {
  //     /* Compute norm of t1(n) and t2(n) */
  //     Scalar normt1n = 0.0;
  //     Scalar normt2n = 0.0;
  //     for (int i=0;i<3;i++) {
  //       normt1n += (t1(n,i) * t1(n,i));
  //       normt2n += (t2(n,i) * t2(n,i));
  //     }
  //     normt1n = sqrt(normt1n);
  //     normt2n = sqrt(normt2n);
  //     /* normalize each tangent now */
  //     for (int i=0;i<3;i++) {
  //       t1(n,i) /= normt1n;
  //       t2(n,i) /= normt2n;
  //     }
  //   }

  //   /* undeformed coordinates */
  //   FieldContainer<Scalar> XYZ(N,3);
  //   for (int i=0;i<N;i++) {
  //     for (int j=0;j<3;j++) {
  //       XYZ(i,j) = L3(i)*warVerts(0,j) + L4(i)*warVerts(1,j) + L2(i)*warVerts(2,j) + L1(i)*warVerts(3,j);
  //     }
  //   }

  //   for (int face=1;face<=4;face++) {
  //     FieldContainer<Scalar> La, Lb, Lc, Ld;
  //     FieldContainer<Scalar> warp(N,2);
  //     FieldContainer<Scalar> blend(N);
  //     FieldContainer<Scalar> denom(N);
  //     switch (face) {
  //     case 1:
  //       La = L1; Lb = L2; Lc = L3; Ld = L4; break;
  //     case 2:
  //       La = L2; Lb = L1; Lc = L3; Ld = L4; break;
  //     case 3:
  //       La = L3; Lb = L1; Lc = L4; Ld = L2; break;
  //     case 4:
  //       La = L4; Lb = L1; Lc = L3; Ld = L2; break;
  //     }
      
  //     /* get warp tangential to face */
  //     warpShiftFace3D<Scalar,FieldContainer<Scalar> >(order,alpha,La,Lb,Lc,Ld,warp);
      
  //     for (int k=0;k<N;k++) {
  //       blend(k) = Lb(k) * Lc(k) * Ld(k);
  //     }

  //     for (int k=0;k<N;k++) {
  //       denom(k) = (Lb(k) + 0.5 * La(k)) * (Lc(k) + 0.5*La(k)) * (Ld(k) + 0.5 * La(k));
  //     }

  //     for (int k=0;k<N;k++) {
  //       if (denom(k) > tol) {
  //         blend(k) *= ( 1.0 + alpha * alpha * La(k) * La(k) ) / denom(k);
  //       }
  //     }  


  //     // compute warp and blend
  //     for (int k=0;k<N;k++) {
  //       for (int j=0;j<3;j++) {
  //         shift(k,j) = shift(k,j) + blend(k) * warp(k,0) * t1(face-1,j)
  //           + blend(k) * warp(k,1) * t2(face-1,j);
  //       }
  //     }

  //     for (int k=0;k<N;k++) {
  //       if (La(k) < tol && ( Lb(k) < tol || Lc(k) < tol || Ld(k) < tol )) {
  //         for (int j=0;j<3;j++) {
  //           shift(k,j) = warp(k,0) * t1(face-1,j) + warp(k,1) * t2(face-1,j);
  //         }
  //       }
  //     }
      
  //   }

  //   FieldContainer<Scalar> updatedPoints(N,3);
  //   for (int k=0;k<N;k++) {
  //     for (int j=0;j<3;j++) {
  //       updatedPoints(k,j) = XYZ(k,j) + shift(k,j);
  //     }
  //   }

  //   warVerts.resize( 1 , 4 , 3 );

  //   // now we convert to Pavel's reference triangle!
  //   FieldContainer<Scalar> refPts(N,3);
  //   CellTools<Scalar>::mapToReferenceFrame( refPts ,updatedPoints ,
  //                                           warVerts ,
  //                                           shards::getCellTopologyData<shards::Tetrahedron<4> >() ,
  //                                           0 );

  //   // now write from refPts into points, taking offset into account
  //   int noffcur = 0;
  //   int offcur = 0;
  //   for (int i=0;i<=order;i++) {
  //     for (int j=0;j<=order-i;j++) {
  //       for (int k=0;k<=order-i-j;k++) {
  //         if ( (i >= offset) && (i <= order-offset) &&
  //             (j >= offset) && (j <= order-i-offset) &&
  //             (k >= offset) && (k <= order-i-j-offset) ) {
  //           points(offcur,0) = refPts(noffcur,0);
  //           points(offcur,1) = refPts(noffcur,1);
  //           points(offcur,2) = refPts(noffcur,2);
  //           offcur++;
  //         }
  //         noffcur++;
  //       }
  //     }
  //   }
                                            


  // }
  

} // face Intrepid
#endif
