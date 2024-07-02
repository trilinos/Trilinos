// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PointToolsDef.hpp
    \brief  Definition file for point tool utilities for barycentric coordinates and lattices
    \author Created by R. Kirby
            Kokkorized by Kyungjoo Kim
 */
#ifndef __INTREPID2_POINTTOOLS_DEF_HPP__
#define __INTREPID2_POINTTOOLS_DEF_HPP__

#if defined(_MSC_VER) || defined(_WIN32) && defined(__ICL)
// M_PI, M_SQRT2, etc. are hidden in MSVC by #ifdef _USE_MATH_DEFINES
  #ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
  #endif
  #include <math.h>
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
  switch (cellType.getBaseKey()) {
  case shards::Tetrahedron<>::key: {
    const auto effectiveOrder = order - 4 * offset;
    r_val = (effectiveOrder < 0 ? 0 :(effectiveOrder+1)*(effectiveOrder+2)*(effectiveOrder+3)/6);
    break;
  }
  case shards::Triangle<>::key: {
    const auto effectiveOrder = order - 3 * offset;
    r_val = (effectiveOrder < 0 ? 0 : (effectiveOrder+1)*(effectiveOrder+2)/2);
    break;
  }
  case shards::Line<>::key: {
    const auto effectiveOrder = order - 2 * offset;
    r_val = (effectiveOrder < 0 ? 0 : (effectiveOrder+1));
    break;
  }
  case shards::Quadrilateral<>::key: {
    const auto effectiveOrder = order - 2 * offset;
    r_val = std::pow(effectiveOrder < 0 ? 0 : (effectiveOrder+1),2);
    break;
  }
  case shards::Hexahedron<>::key: {
    const auto effectiveOrder = order - 2 * offset;
    r_val = std::pow(effectiveOrder < 0 ? 0 : (effectiveOrder+1),3);
    break;
  }
  case shards::Pyramid<>::key: {
    const auto effectiveOrder = order - 2 * offset;
    r_val = (effectiveOrder < 0 ? 0 : (effectiveOrder+1)*(effectiveOrder+2)*(2*effectiveOrder+3)/6);
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
getLattice(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
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

  const size_type latticeSize = getLatticeSize( cell, order, offset );
  const size_type spaceDim = cell.getDimension();

  INTREPID2_TEST_FOR_EXCEPTION( points.extent(0) != latticeSize ||
      points.extent(1) != spaceDim,
      std::invalid_argument ,
      ">>> ERROR (PointTools::getLattice): dimension does not match to lattice size." );
#endif

  switch (cell.getBaseKey()) {
  case shards::Tetrahedron<>::key: getLatticeTetrahedron( points, order, offset, pointType );  break;
  case shards::Pyramid<>::key:     getLatticePyramid    ( points, order, offset, pointType );  break;
  case shards::Triangle<>::key:    getLatticeTriangle   ( points, order, offset, pointType );  break;
  case shards::Line<>::key:        getLatticeLine       ( points, order, offset, pointType );  break;
  case shards::Quadrilateral<>::key:   {
    auto hostPoints = Kokkos::create_mirror_view(points);
    shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());
    const ordinal_type numPoints = getLatticeSize( line, order, offset );
    auto linePoints = getMatchingViewWithLabel(hostPoints, "linePoints", numPoints, 1);
    getLatticeLine( linePoints, order, offset, pointType );
    ordinal_type idx=0;
    for (ordinal_type j=0; j<numPoints; ++j) {
      for (ordinal_type i=0; i<numPoints; ++i, ++idx) {
        hostPoints(idx,0) = linePoints(i,0);
        hostPoints(idx,1) = linePoints(j,0);
      }
    }
    Kokkos::deep_copy(points,hostPoints);
  }
  break;
  case shards::Hexahedron<>::key:   {
    auto hostPoints = Kokkos::create_mirror_view(points);
    shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());
    const ordinal_type numPoints = getLatticeSize( line, order, offset );
    auto linePoints = getMatchingViewWithLabel(hostPoints, "linePoints", numPoints, 1);
    getLatticeLine( linePoints, order, offset, pointType );
    ordinal_type idx=0;
    for (ordinal_type k=0; k<numPoints; ++k) {
      for (ordinal_type j=0; j<numPoints; ++j) {
        for (ordinal_type i=0; i<numPoints; ++i, ++idx) {
          hostPoints(idx,0) = linePoints(i,0);
          hostPoints(idx,1) = linePoints(j,0);
          hostPoints(idx,2) = linePoints(k,0);
        }
      }
    }
    Kokkos::deep_copy(points,hostPoints);
  }
  break;
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
        ">>> ERROR (Intrepid2::PointTools::getLattice): the specified cell type is not supported." );
  }
  }
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getGaussPoints(      Kokkos::DynRankView<pointValueType,pointProperties...> points,
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
  auto pts = Kokkos::subview( points, range_type(0,np), 0 );
  // should be fixed after view and dynrankview are inter-operatible
  auto z   = Kokkos::DynRankView<pointValueType,Kokkos::HostSpace>(zHost.data(), np);
  Kokkos::deep_copy(pts, z);
}

// -----------------------------------------------------------------------------------------
// Internal implementation
// -----------------------------------------------------------------------------------------

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getLatticeLine(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
    const ordinal_type         order,
    const ordinal_type         offset,
    const EPointType           pointType ) {
  switch (pointType) {
  case POINTTYPE_EQUISPACED:  getEquispacedLatticeLine( points, order, offset ); break;
  case POINTTYPE_WARPBLEND:   getWarpBlendLatticeLine( points, order, offset ); break;
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
getLatticeTriangle(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
    const ordinal_type         order,
    const ordinal_type         offset,
    const EPointType           pointType ) {
  switch (pointType) {
  case POINTTYPE_EQUISPACED:  getEquispacedLatticeTriangle( points, order, offset ); break;
  case POINTTYPE_WARPBLEND:   getWarpBlendLatticeTriangle ( points, order, offset ); break;
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
getLatticeTetrahedron(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
    const ordinal_type         order,
    const ordinal_type         offset,
    const EPointType           pointType ) {
  switch (pointType) {
    case POINTTYPE_EQUISPACED:  getEquispacedLatticeTetrahedron( points, order, offset ); break;
    case POINTTYPE_WARPBLEND:   getWarpBlendLatticeTetrahedron ( points, order, offset ); break;
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
getLatticePyramid(     Kokkos::DynRankView<pointValueType,pointProperties...> points,
                  const ordinal_type order,
                  const ordinal_type offset,
                  const EPointType pointType )
{
  switch (pointType) {
    case POINTTYPE_EQUISPACED:  getEquispacedLatticePyramid( points, order, offset ); break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true ,
          std::invalid_argument ,
          ">>> ERROR (PointTools::getLattice): invalid EPointType." );
    }
  }
}

// -----------------------------------------------------------------------------------------

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getEquispacedLatticeLine(      Kokkos::DynRankView<pointValueType,pointProperties...> points,
    const ordinal_type order,
    const ordinal_type offset ) {
  auto pointsHost = Kokkos::create_mirror_view(points);

  if (order == 0)
    pointsHost(0,0) = 0.0;
  else {
    const pointValueType h = 2.0 / order;
    const ordinal_type ibeg = offset, iend = order-offset+1;
    for (ordinal_type i=ibeg;i<iend;++i)
      pointsHost(i-ibeg, 0) = -1.0 + h * (pointValueType) i;
  }

  Kokkos::deep_copy(points, pointsHost);
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getWarpBlendLatticeLine(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
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
    auto pts = Kokkos::subview( points, range_type(0, s), 0 );

    // this should be fixed after view and dynrankview is interoperatable
    auto z   = Kokkos::DynRankView<pointValueType,Kokkos::HostSpace>(zHost.data() + offset, np-offset);

    const auto common_range = range_type(0, std::min(pts.extent(0), z.extent(0)));
    Kokkos::deep_copy(Kokkos::subview(pts, common_range), 
                      Kokkos::subview(z,   common_range));
  }
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getEquispacedLatticeTriangle(       Kokkos::DynRankView<pointValueType,pointProperties...> points,
    const ordinal_type order,
    const ordinal_type offset ) {
    TEUCHOS_TEST_FOR_EXCEPTION( order <= 0 ,
      std::invalid_argument ,
      ">>> ERROR (Intrepid2::PointTools::getEquispacedLatticeTriangle): order must be positive" );

  auto pointsHost = Kokkos::create_mirror_view(points);

  const pointValueType h = 1.0 / order;
  ordinal_type cur = 0;

  for (ordinal_type i=offset;i<=order-offset;i++) {
    for (ordinal_type j=offset;j<=order-i-offset;j++) {
      pointsHost(cur,0) =  j * h ;
      pointsHost(cur,1) =  i * h;
      cur++;
    }
  }
  Kokkos::deep_copy(points, pointsHost);
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getEquispacedLatticeTetrahedron( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                                 const ordinal_type order ,
                                 const ordinal_type offset ) {
  TEUCHOS_TEST_FOR_EXCEPTION( (order <= 0) ,
    std::invalid_argument ,
    ">>> ERROR (Intrepid2::PointTools::getEquispacedLatticeTetrahedron): order must be positive" );

  auto pointsHost = Kokkos::create_mirror_view(points);

  const pointValueType h = 1.0 / order;
  ordinal_type cur = 0;

  for (ordinal_type i=offset;i<=order-offset;i++) {
    for (ordinal_type j=offset;j<=order-i-offset;j++) {
      for (ordinal_type k=offset;k<=order-i-j-offset;k++) {
        pointsHost(cur,0) = k * h;
        pointsHost(cur,1) = j * h;
        pointsHost(cur,2) = i * h;
        cur++;
      }
    }
  }
  Kokkos::deep_copy(points, pointsHost);
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getEquispacedLatticePyramid( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                             const ordinal_type order ,
                             const ordinal_type offset ) {
  TEUCHOS_TEST_FOR_EXCEPTION( (order <= 0) ,
    std::invalid_argument ,
    ">>> ERROR (Intrepid2::PointTools::getEquispacedLatticePyramid): order must be positive" );

  auto pointsHost = Kokkos::create_mirror_view(points);

  const pointValueType h = 1.0 / order;
  ordinal_type cur = 0;

  for (ordinal_type i=offset;i<=order-offset;i++) { // z dimension (goes from 0 to 1)
    pointValueType z = i * h;
    for (ordinal_type j=offset;j<=order-i-offset;j++) { // y (goes from -(1-z) to 1-z)
      for (ordinal_type k=offset;k<=order-i-offset;k++) { // x (goes from -(1-z) to 1-z)
        pointsHost(cur,0) = (2 * k * h - 1.) * (1-z);
        pointsHost(cur,1) = (2 * j * h - 1.) * (1-z);
        pointsHost(cur,2) = z;
        cur++;
      }
    }
  }
  Kokkos::deep_copy(points, pointsHost);
}

template<typename pointValueType, class ...pointProperties>
void
PointTools::
warpFactor( Kokkos::DynRankView<pointValueType,pointProperties...> warp,
            const ordinal_type order,
            const Kokkos::DynRankView<pointValueType,pointProperties...> xnodes ,
            const Kokkos::DynRankView<pointValueType,pointProperties...> xout
            )
 {
   TEUCHOS_TEST_FOR_EXCEPTION( ( warp.extent(0) != xout.extent(0) ) ,
                       std::invalid_argument ,
                       ">>> ERROR (PointTools::warpFactor): xout and warp must be same size." );

   Kokkos::deep_copy(warp, pointValueType(0.0));

   ordinal_type xout_dim0 = xout.extent(0);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> d("d", xout_dim0 );

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> xeq_("xeq", order + 1 ,1);
   PointTools::getEquispacedLatticeLine( xeq_ , order , 0 );
   const auto xeq = Kokkos::subview(xeq_, Kokkos::ALL(),0);

   TEUCHOS_TEST_FOR_EXCEPTION( ( xeq.extent(0) != xnodes.extent(0) ) ,
                       std::invalid_argument ,
                       ">>> ERROR (PointTools::warpFactor): xeq and xnodes must be same size." );


   for (ordinal_type i=0;i<=order;i++) {

     Kokkos::deep_copy(d, xnodes(i,0) - xeq(i));

     for (ordinal_type j=1;j<order;j++) {
       if (i != j) {
         for (ordinal_type k=0;k<xout_dim0;k++) {
           d(k) = d(k) * ( (xout(k)-xeq(j)) / (xeq(i)-xeq(j)) );
         }
       }
     }

     // deflate end roots
     if ( i != 0 ) {
       for (ordinal_type k=0;k<xout_dim0;k++) {
         d(k) = -d(k) / (xeq(i) - xeq(0));
       }
     }

     if (i != order ) {
       for (ordinal_type k=0;k<xout_dim0;k++) {
         d(k) = d(k) / (xeq(i) - xeq(order));
       }
     }

     for (ordinal_type k=0;k<xout_dim0;k++) {
       warp(k) += d(k);
     }

   }
 }

template<typename pointValueType, class ...pointProperties>
void
PointTools::
getWarpBlendLatticeTriangle( Kokkos::DynRankView<pointValueType,pointProperties...> points,
                             const ordinal_type order ,
                             const ordinal_type offset)
 {
   /* get Gauss-Lobatto points */

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> gaussX("gaussX", order + 1 , 1 );

   PointTools::getWarpBlendLatticeLine( gaussX , order , 0 );
   //auto gaussX = Kokkos::subdynrankview(gaussX_, Kokkos::ALL(),0);

 //  gaussX.resize(gaussX.extent(0));

   pointValueType alpopt[] = {0.0000,0.0000,1.4152,0.1001,0.2751,0.9800,1.0999,
                       1.2832,1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};

   pointValueType alpha;

   if (order >= 1 && order < 16) {
     alpha = alpopt[order-1];
   }
   else {
     alpha = 5.0 / 3.0;
   }

   const ordinal_type p = order; /* switch to Warburton's notation */
   ordinal_type N = (p+1)*(p+2)/2;

   /* equidistributed nodes on equilateral triangle */
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> L1("L1", N );
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace>  L2("L2", N );
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace>  L3("L3", N );
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace>  X("X", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace>  Y("Y", N);

   ordinal_type sk = 0;
   for (ordinal_type n=1;n<=p+1;n++) {
     for (ordinal_type m=1;m<=p+2-n;m++) {
       L1(sk) = (n-1) / (pointValueType)p;
       L3(sk) = (m-1) / (pointValueType)p;
       L2(sk) = 1.0 - L1(sk) - L3(sk);
       sk++;
     }
   }

   for (ordinal_type n=0;n<N;n++) {
     X(n) = -L2(n) + L3(n);
     Y(n) = (-L2(n) - L3(n) + 2*L1(n))/1.7320508075688772;
   }

   /* get blending function for each node at each edge */
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> blend1("blend1", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> blend2("blend2", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> blend3("blend3", N);

   for (ordinal_type n=0;n<N;n++) {
     blend1(n) = 4.0 * L2(n) * L3(n);
     blend2(n) = 4.0 * L1(n) * L3(n);
     blend3(n) = 4.0 * L1(n) * L2(n);
   }

   /* get difference of each barycentric coordinate */
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> L3mL2("L3mL2", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> L1mL3("L1mL3", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> L2mL1("L2mL1", N);

   for (ordinal_type k=0;k<N;k++) {
     L3mL2(k) = L3(k)-L2(k);
     L1mL3(k) = L1(k)-L3(k);
     L2mL1(k) = L2(k)-L1(k);
   }

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warpfactor1("warpfactor1", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warpfactor2("warpfactor2", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warpfactor3("warpfactor3", N);

   warpFactor( warpfactor1, order , gaussX , L3mL2 );
   warpFactor( warpfactor2, order , gaussX , L1mL3 );
   warpFactor( warpfactor3, order , gaussX , L2mL1 );

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warp1("warp1", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warp2("warp2", N);
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warp3("warp3", N);

   for (ordinal_type k=0;k<N;k++) {
     warp1(k) = blend1(k) * warpfactor1(k) *
       ( 1.0 + alpha * alpha * L1(k) * L1(k) );
     warp2(k) = blend2(k) * warpfactor2(k) *
       ( 1.0 + alpha * alpha * L2(k) * L2(k) );
     warp3(k) = blend3(k) * warpfactor3(k) *
       ( 1.0 + alpha * alpha * L3(k) * L3(k) );
   }

   for (ordinal_type k=0;k<N;k++) {
     X(k) += 1.0 * warp1(k) + cos( 2.0 * M_PI / 3.0 ) * warp2(k) + cos(4*M_PI/3.0) * warp3(k);
     Y(k) += 0.0 * warp1(k) + sin( 2.0 * M_PI / 3.0 ) * warp2(k) + sin( 4*M_PI/3.0) * warp3(k);
   }

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warXY("warXY", 1, N,2);

   for (ordinal_type k=0;k<N;k++) {
     warXY(0, k,0) = X(k);
     warXY(0, k,1) = Y(k);
   }


   // finally, convert the warp-blend points to the correct triangle
   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> warburtonVerts("warburtonVerts", 1,3,2);
   warburtonVerts(0,0,0) = -1.0;
   warburtonVerts(0,0,1) = -1.0/std::sqrt(3.0);
   warburtonVerts(0,1,0) = 1.0;
   warburtonVerts(0,1,1) = -1.0/std::sqrt(3.0);
   warburtonVerts(0,2,0) = 0.0;
   warburtonVerts(0,2,1) = 2.0/std::sqrt(3.0);

   Kokkos::DynRankView<pointValueType, Kokkos::DefaultHostExecutionSpace> refPts("refPts", 1, N,2);

   Intrepid2::CellTools<Kokkos::HostSpace>::mapToReferenceFrame( refPts ,
                                              warXY ,
                                              warburtonVerts ,
                                              shards::getCellTopologyData< shards::Triangle<3> >()
                                              );


   auto pointsHost = Kokkos::create_mirror_view(points);
   // now write from refPts into points, taking care of offset
   ordinal_type noffcur = 0;  // index into refPts
   ordinal_type offcur = 0;   // index ordinal_type points
   for (ordinal_type i=0;i<=order;i++) {
     for (ordinal_type j=0;j<=order-i;j++) {
       if ( (i >= offset) && (i <= order-offset) &&
             (j >= offset) && (j <= order-i-offset) ) {
         pointsHost(offcur,0) = refPts(0, noffcur,0);
         pointsHost(offcur,1) = refPts(0, noffcur,1);
         offcur++;
       }
       noffcur++;
     }
   }

   Kokkos::deep_copy(points, pointsHost);

 }


template<typename pointValueType, class ...pointProperties>
void
PointTools::
warpShiftFace3D( Kokkos::DynRankView<pointValueType,pointProperties...>  dxy,
                 const ordinal_type order ,
                 const pointValueType pval ,
                 const Kokkos::DynRankView<pointValueType,pointProperties...>  /* L1 */,
                 const Kokkos::DynRankView<pointValueType,pointProperties...>  L2,
                 const Kokkos::DynRankView<pointValueType,pointProperties...>  L3,
                 const Kokkos::DynRankView<pointValueType,pointProperties...> L4
                )
 {
   evalshift(dxy,order,pval,L2,L3,L4);
   return;
 }

template<typename pointValueType, class ...pointProperties>
void
PointTools::
evalshift( Kokkos::DynRankView<pointValueType,pointProperties...>  dxy,
           const ordinal_type order ,
           const pointValueType pval ,
           const Kokkos::DynRankView<pointValueType,pointProperties...>  L1 ,
           const Kokkos::DynRankView<pointValueType,pointProperties...>  L2 ,
           const Kokkos::DynRankView<pointValueType,pointProperties...>  L3
           )
 {
   // get Gauss-Lobatto-nodes
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> gaussX("gaussX",order+1,1);
   PointTools::getWarpBlendLatticeLine( gaussX , order , 0 );
   //gaussX.resize(order+1);
   const ordinal_type N = L1.extent(0);

   // Warburton code reverses them
   for (ordinal_type k=0;k<=order;k++) {
     gaussX(k,0) *= -1.0;
   }

   // blending function at each node for each edge
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> blend1("blend1",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> blend2("blend2",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> blend3("blend3",N);

   for (ordinal_type i=0;i<N;i++) {
     blend1(i) = L2(i) * L3(i);
     blend2(i) = L1(i) * L3(i);
     blend3(i) = L1(i) * L2(i);
   }

   // amount of warp for each node for each edge
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warpfactor1("warpfactor1",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warpfactor2("warpfactor2",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warpfactor3("warpfactor3",N);

   // differences of barycentric coordinates
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L3mL2("L3mL2",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L1mL3("L1mL3",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L2mL1("L2mL1",N);

   for (ordinal_type k=0;k<N;k++) {
     L3mL2(k) = L3(k)-L2(k);
     L1mL3(k) = L1(k)-L3(k);
     L2mL1(k) = L2(k)-L1(k);
   }

   evalwarp( warpfactor1 , order , gaussX , L3mL2 );
   evalwarp( warpfactor2 , order , gaussX , L1mL3 );
   evalwarp( warpfactor3 , order , gaussX , L2mL1 );

   for (ordinal_type k=0;k<N;k++) {
     warpfactor1(k) *= 4.0;
     warpfactor2(k) *= 4.0;
     warpfactor3(k) *= 4.0;
   }

   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warp1("warp1",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warp2("warp2",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warp3("warp3",N);

   for (ordinal_type k=0;k<N;k++) {
     warp1(k) = blend1(k) * warpfactor1(k) *
       ( 1.0 + pval * pval * L1(k) * L1(k) );
     warp2(k) = blend2(k) * warpfactor2(k) *
       ( 1.0 + pval * pval * L2(k) * L2(k) );
     warp3(k) = blend3(k) * warpfactor3(k) *
       ( 1.0 + pval * pval * L3(k) * L3(k) );
   }

   for (ordinal_type k=0;k<N;k++) {
     dxy(k,0) = 1.0 * warp1(k) + cos( 2.0 * M_PI / 3.0 ) * warp2(k) + cos( 4.0*M_PI/3.0 ) * warp3(k);
     dxy(k,1) = 0.0 * warp1(k) + sin( 2.0 * M_PI / 3.0 ) * warp2(k) + sin( 4.0*M_PI/3.0 ) * warp3(k);
   }
 }

 /* one-d edge warping function */
template<typename pointValueType, class ...pointProperties>
void
PointTools::
evalwarp(Kokkos::DynRankView<pointValueType,pointProperties...>  warp ,
         const ordinal_type order ,
         const Kokkos::DynRankView<pointValueType,pointProperties...>  xnodes,
         const Kokkos::DynRankView<pointValueType,pointProperties...>  xout )
 {
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> xeq("xeq",order+1);

   ordinal_type xout_dim0 = xout.extent(0);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> d("d",xout_dim0);

   //Kokkos::deep_copy(d, 0.0);

   for (ordinal_type i=0;i<=order;i++) {
     xeq(i) = -1.0 + 2.0 * ( order - i ) / order;
   }

   for (ordinal_type i=0;i<=order;i++) {
     Kokkos::deep_copy(d, xnodes(i,0) - xeq(i));
     for (ordinal_type j=1;j<order;j++) {
       if (i!=j) {
         for (ordinal_type k=0;k<xout_dim0;k++) {
           d(k) = d(k) * (xout(k)-xeq(j))/(xeq(i)-xeq(j));
         }
       }
     }
     if (i!=0) {
       for (ordinal_type k=0;k<xout_dim0;k++) {
         d(k) = -d(k)/(xeq(i)-xeq(0));
       }
     }
     if (i!=order) {
       for (ordinal_type k=0;k<xout_dim0;k++) {
         d(k) = d(k)/(xeq(i)-xeq(order));
       }
     }

     for (ordinal_type k=0;k<xout_dim0;k++) {
       warp(k) += d(k);
     }
   }
 }


template<typename pointValueType, class ...pointProperties>
void
PointTools::
getWarpBlendLatticeTetrahedron(Kokkos::DynRankView<pointValueType,pointProperties...> points,
                               const ordinal_type order ,
                               const ordinal_type offset  )
 {
   pointValueType alphastore[] = { 0,0,0,0.1002, 1.1332,1.5608,1.3413,1.2577,1.1603,
                           1.10153,0.6080,0.4523,0.8856,0.8717,0.9655};
   pointValueType alpha;

   if (order <= 15) {
     alpha = alphastore[std::max(order-1,0)];
   }
   else {
     alpha = 1.0;
   }

   const ordinal_type N = (order+1)*(order+2)*(order+3)/6;
   pointValueType tol = 1.e-10;

   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> shift("shift",N,3);
   Kokkos::deep_copy(shift,0.0);

   /* create 3d equidistributed nodes on Warburton tet */
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> equipoints("equipoints", N,3);
   ordinal_type sk = 0;
   for (ordinal_type n=0;n<=order;n++) {
     for (ordinal_type m=0;m<=order-n;m++) {
       for (ordinal_type q=0;q<=order-n-m;q++) {
         equipoints(sk,0) = -1.0 + (q * 2.0 ) / order;
         equipoints(sk,1) = -1.0 + (m * 2.0 ) / order;
         equipoints(sk,2) = -1.0 + (n * 2.0 ) / order;
         sk++;
       }
     }
   }


   /* construct barycentric coordinates for equispaced lattice */
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L1("L1",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L2("L2",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L3("L3",N);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> L4("L4",N);
   for (ordinal_type i=0;i<N;i++) {
     L1(i) = (1.0 + equipoints(i,2)) / 2.0;
     L2(i) = (1.0 + equipoints(i,1)) / 2.0;
     L3(i) = -(1.0 + equipoints(i,0) + equipoints(i,1) + equipoints(i,2)) / 2.0;
     L4(i) = (1.0 + equipoints(i,0)) / 2.0;
   }


   /* vertices of Warburton tet */
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warVerts_("warVerts",1,4,3);
   auto warVerts = Kokkos::subview(warVerts_,0,Kokkos::ALL(),Kokkos::ALL());
   warVerts(0,0) = -1.0;
   warVerts(0,1) = -1.0/sqrt(3.0);
   warVerts(0,2) = -1.0/sqrt(6.0);
   warVerts(1,0) = 1.0;
   warVerts(1,1) = -1.0/sqrt(3.0);
   warVerts(1,2) = -1.0/sqrt(6.0);
   warVerts(2,0) = 0.0;
   warVerts(2,1) = 2.0 / sqrt(3.0);
   warVerts(2,2) = -1.0/sqrt(6.0);
   warVerts(3,0) = 0.0;
   warVerts(3,1) = 0.0;
   warVerts(3,2) = 3.0 / sqrt(6.0);


   /* tangents to faces */
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> t1("t1",4,3);
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> t2("t2",4,3);
   for (ordinal_type i=0;i<3;i++) {
     t1(0,i) = warVerts(1,i) - warVerts(0,i);
     t1(1,i) = warVerts(1,i) - warVerts(0,i);
     t1(2,i) = warVerts(2,i) - warVerts(1,i);
     t1(3,i) = warVerts(2,i) - warVerts(0,i);
     t2(0,i) = warVerts(2,i) - 0.5 * ( warVerts(0,i) + warVerts(1,i) );
     t2(1,i) = warVerts(3,i) - 0.5 * ( warVerts(0,i) + warVerts(1,i) );
     t2(2,i) = warVerts(3,i) - 0.5 * ( warVerts(1,i) + warVerts(2,i) );
     t2(3,i) = warVerts(3,i) - 0.5 * ( warVerts(0,i) + warVerts(2,i) );
   }

   /* normalize tangents */
   for (ordinal_type n=0;n<4;n++) {
     /* Compute norm of t1(n) and t2(n) */
     pointValueType normt1n = 0.0;
     pointValueType normt2n = 0.0;
     for (ordinal_type i=0;i<3;i++) {
       normt1n += (t1(n,i) * t1(n,i));
       normt2n += (t2(n,i) * t2(n,i));
     }
     normt1n = sqrt(normt1n);
     normt2n = sqrt(normt2n);
     /* normalize each tangent now */
     for (ordinal_type i=0;i<3;i++) {
       t1(n,i) /= normt1n;
       t2(n,i) /= normt2n;
     }
   }

   /* undeformed coordinates */
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> XYZ("XYZ",N,3);
   for (ordinal_type i=0;i<N;i++) {
     for (ordinal_type j=0;j<3;j++) {
       XYZ(i,j) = L3(i)*warVerts(0,j) + L4(i)*warVerts(1,j) + L2(i)*warVerts(2,j) + L1(i)*warVerts(3,j);
     }
   }

   for (ordinal_type face=1;face<=4;face++) {
     Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> La, Lb, Lc, Ld;
     Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> warp("warp",N,2);
     Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> blend("blend",N);
     Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> denom("denom",N);
     switch (face) {
     case 1:
       La = L1; Lb = L2; Lc = L3; Ld = L4; break;
     case 2:
       La = L2; Lb = L1; Lc = L3; Ld = L4; break;
     case 3:
       La = L3; Lb = L1; Lc = L4; Ld = L2; break;
     case 4:
       La = L4; Lb = L1; Lc = L3; Ld = L2; break;
     }

     /* get warp tangential to face */
     warpShiftFace3D(warp,order,alpha,La,Lb,Lc,Ld);

     for (ordinal_type k=0;k<N;k++) {
       blend(k) = Lb(k) * Lc(k) * Ld(k);
     }

     for (ordinal_type k=0;k<N;k++) {
       denom(k) = (Lb(k) + 0.5 * La(k)) * (Lc(k) + 0.5*La(k)) * (Ld(k) + 0.5 * La(k));
     }

     for (ordinal_type k=0;k<N;k++) {
       if (denom(k) > tol) {
         blend(k) *= ( 1.0 + alpha * alpha * La(k) * La(k) ) / denom(k);
       }
     }


     // compute warp and blend
     for (ordinal_type k=0;k<N;k++) {
       for (ordinal_type j=0;j<3;j++) {
         shift(k,j) = shift(k,j) + blend(k) * warp(k,0) * t1(face-1,j)
           + blend(k) * warp(k,1) * t2(face-1,j);
       }
     }

     for (ordinal_type k=0;k<N;k++) {
       if (La(k) < tol && ( Lb(k) < tol || Lc(k) < tol || Ld(k) < tol )) {
         for (ordinal_type j=0;j<3;j++) {
           shift(k,j) = warp(k,0) * t1(face-1,j) + warp(k,1) * t2(face-1,j);
         }
       }
     }

   }

   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> updatedPoints("updatedPoints",1,N,3);
   for (ordinal_type k=0;k<N;k++) {
     for (ordinal_type j=0;j<3;j++) {
       updatedPoints(0,k,j) = XYZ(k,j) + shift(k,j);
     }
   }

   //warVerts.resize( 1 , 4 , 3 );

   // now we convert to Pavel's reference triangle!
   Kokkos::DynRankView<pointValueType,Kokkos::DefaultHostExecutionSpace> refPts("refPts",1,N,3);
   CellTools<Kokkos::HostSpace>::mapToReferenceFrame( refPts ,updatedPoints ,
                                           warVerts_ ,
                                           shards::getCellTopologyData<shards::Tetrahedron<4> >()
                                           );

   auto pointsHost = Kokkos::create_mirror_view(points);
   // now write from refPts into points, taking offset into account
   ordinal_type noffcur = 0;
   ordinal_type offcur = 0;
   for (ordinal_type i=0;i<=order;i++) {
     for (ordinal_type j=0;j<=order-i;j++) {
       for (ordinal_type k=0;k<=order-i-j;k++) {
         if ( (i >= offset) && (i <= order-offset) &&
             (j >= offset) && (j <= order-i-offset) &&
             (k >= offset) && (k <= order-i-j-offset) ) {
           pointsHost(offcur,0) = refPts(0,noffcur,0);
           pointsHost(offcur,1) = refPts(0,noffcur,1);
           pointsHost(offcur,2) = refPts(0,noffcur,2);
           offcur++;
         }
         noffcur++;
       }
     }
   }

   Kokkos::deep_copy(points, pointsHost);
 }


} // face Intrepid
#endif
