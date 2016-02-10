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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef INTREPID2_ORIENTATIONTOOLSDEF_HPP
#define INTREPID2_ORIENTATIONTOOLSDEF_HPP

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

namespace Intrepid2 {

  template<class Scalar>
  template<class ArrayPoint>
  void OrientationTools<Scalar>::getTriangleLatticePointsByTopology(ArrayPoint &                     outPts,
                                                                    const ArrayPoint &               refPts,
                                                                    const Basis<Scalar,ArrayPoint> & basis) {
    const int nv = basis.getNumVertexDofs();
    const int ne = basis.getNumEdgeDofs();
    const int ni = basis.getNumInteriorDofs();
    const int nn = 3*nv + ne*3 + ni;

    // group (vertex, edges, interior)
    int eoff0 = 3*nv, eoff1 = eoff0+ne, eoff2 = eoff1+ne, ioff = eoff2+ne;

    // group dofs with bicentric coordiantes
    const double eps = INTREPID2_EPSILON;

    for (int k=0;k<nn;++k) {
      const Scalar pt[] = { refPts(k, 0),
                            refPts(k, 1) };

      const Scalar lambda[] = { (1.0 - pt[0] - pt[1])*(1.0 - pt[0] - pt[1]),
                                pt[0]*pt[0],
                                pt[1]*pt[1] };

      const bool vflag[] = { (lambda[0] - 1.0 < eps),
                             (lambda[1] - 1.0 < eps),
                             (lambda[2] - 1.0 < eps) };

      const bool eflag[] = { (lambda[2] < eps),
                             (lambda[0] < eps),
                             (lambda[1] < eps) };

      if      ( vflag[0] &&  eflag[1] &&  eflag[2]) { outPts(0,0) = refPts(k,0); outPts(0,1) = refPts(k,1); } //vert0
      else if ( eflag[0] &&  vflag[1] &&  eflag[2]) { outPts(1,0) = refPts(k,0); outPts(1,1) = refPts(k,1); } //vert1
      else if ( eflag[0] &&  eflag[1] &&  vflag[2]) { outPts(2,0) = refPts(k,0); outPts(2,1) = refPts(k,1); } //vert2
      else if ( eflag[0] && !eflag[1] && !eflag[2]) { outPts(eoff0, 0) = refPts(k,0); outPts(eoff0++, 1) = refPts(k,1); } //edge0 
      else if (!eflag[0] &&  eflag[1] && !eflag[2]) { outPts(eoff1, 0) = refPts(k,0); outPts(eoff1++, 1) = refPts(k,1); } //edge1 
      else if (!eflag[0] && !eflag[1] &&  eflag[2]) { outPts(eoff2, 0) = refPts(k,0); outPts(eoff2++, 1) = refPts(k,1); } //edge2
      else                                          { outPts(ioff,  0) = refPts(k,0); outPts(ioff++,  1) = refPts(k,1); } //interior
    }
  }

  
  template<class Scalar>
  void OrientationTools<Scalar>::getModifiedLinePoint(double &ot, 
                                                      const double pt,
                                                      const int ort) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( !( -1.0 <= pt && pt <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedLinePoint): pt is out of range [-1, 1].");  
#endif
  
    switch (ort) {
    case 0: ot =   pt; break;
    case 1: ot = - pt; break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::invalid_argument,
                                 ">>> ERROR (Intrepid2::OrientationTools::getModifiedLinePoint): Invalid orientation number (0--1)." );
    }
  }

  template<class Scalar>
  void OrientationTools<Scalar>::getModifiedTrianglePoint(double &ot0,     
                                                          double &ot1, 
                                                          const double pt0, 
                                                          const double pt1,
                                                          const int ort) {
    const double lambda[3] = { 1.0 - pt0 - pt1, 
                               pt0,
                               pt1 };
  
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( !( 0.0 <= lambda[0] && lambda[0] <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): Bicentric coordinate (lamba[0]) is out of range [0, 1].");  
  
    TEUCHOS_TEST_FOR_EXCEPTION( !( 0.0 <= lambda[1] && lambda[1] <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): Bicentric coordinate (lamba[1]) is out of range [0, 1].");  
  
    TEUCHOS_TEST_FOR_EXCEPTION( !( 0.0 <= lambda[2] && lambda[2] <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): Bicentric coordinate (lamba[2]) is out of range [0, 1].");  
#endif

    switch (ort) {
    case 0: ot0 = lambda[1]; ot1 = lambda[2]; break;
    case 1: ot0 = lambda[2]; ot1 = lambda[0]; break;
    case 2: ot0 = lambda[0]; ot1 = lambda[1]; break;
    case 3: ot0 = lambda[2]; ot1 = lambda[1]; break;
    case 4: ot0 = lambda[0]; ot1 = lambda[2]; break;
    case 5: ot0 = lambda[1]; ot1 = lambda[0]; break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::invalid_argument,
                                 ">>> ERROR (Intrepid2::OrientationTools::getModifiedTrianglePoint): Invalid orientation number (0--5)." );
    }
  }

  template<class Scalar>
  void OrientationTools<Scalar>::getModifiedQuadrilateralPoint(double &ot0,     
                                                               double &ot1, 
                                                               const double pt0, 
                                                               const double pt1,
                                                               const int ort) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( !( -1.0 <= pt0 && pt0 <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): pt0 is out of range [-1, 1].");  
  
    TEUCHOS_TEST_FOR_EXCEPTION( !( -1.0 <= pt1 && pt1 <= 1.0 ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): pt1 is out of range [-1, 1].");  
#endif

    const double lambda[2][2] = { { pt0, -pt0 }, 
                                  { pt1, -pt1 } };
  
    switch (ort) {
    case 0: ot0 = lambda[0][0]; ot1 = lambda[1][0]; break;
    case 1: ot0 = lambda[1][0]; ot1 = lambda[0][1]; break;
    case 2: ot0 = lambda[0][1]; ot1 = lambda[0][1]; break;
    case 3: ot0 = lambda[1][1]; ot1 = lambda[0][0]; break;
    case 4: ot0 = lambda[1][0]; ot1 = lambda[0][0]; break;
    case 5: ot0 = lambda[0][0]; ot1 = lambda[1][1]; break;
    case 6: ot0 = lambda[1][1]; ot1 = lambda[1][1]; break;
    case 7: ot0 = lambda[0][1]; ot1 = lambda[1][0]; break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::invalid_argument,
                                 ">>> ERROR (Intrepid2::OrientationTools::getModifiedQuadrilateralPoint): Invalid orientation number (0--7)." );
    }
  }

  // -- Public interface

  template<class Scalar>
  template<class ArrayPoint>
  void OrientationTools<Scalar>::getLatticePointsByTopology(ArrayPoint &                     outPoints,
                                                            const ArrayPoint &               refPoints,
                                                            const Basis<Scalar,ArrayPoint> & basis) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( !(hasReferenceCell(cellTopo) ), std::invalid_argument, 
                                ">>> ERROR (Intrepid::OrientationTools::getLatticePointsByTopology): populate point array by a topological order.");

    TEUCHOS_TEST_FOR_EXCEPTION( !( outPoints.dimension(0) == refPoints.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): size of input and output point arrays does not match each other.");
#endif
    
    const auto key = basis.getBaseCellTopology().getBaseCellTopologyData()->key;
    switch (key) {
    case shards::Line<>::key : {
      break;
    }
    case shards::Triangle<>::key : {
      getTriangleLatticePointsByTopology(outPoints, 
                                         refPoints,
                                         basis);
      break;
    }
    case shards::Quadrilateral<>::key : {
      break;
    }
    default: 
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::invalid_argument,
                                 ">>> ERROR (Intrepid2::OrientationTools::getLatticePointsByTopology): Invalid cell topology." );
    }
  }

  template<class Scalar>
  template<class ArrayPoint>
  void OrientationTools<Scalar>::mapToModifiedReference(ArrayPoint &                  ortPoints,
                                                        const ArrayPoint &            refPoints,
                                                        const shards::CellTopology &  cellTopo,
                                                        const int                     cellOrt) {
#ifdef HAVE_INTREPID_DEBUG
    {
      const int cellDim = cellTopo.getDimension();
      TEUCHOS_TEST_FOR_EXCEPTION( !(hasReferenceCell(cellTopo) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): the specified cell topology does not have a reference cell.");
    
      TEUCHOS_TEST_FOR_EXCEPTION( !( (1 <= cellDim) && (cellDim <= 2 ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): method defined only for 1 and 2-dimensional subcells.");  

      TEUCHOS_TEST_FOR_EXCEPTION( !( ortPoints.dimension(0) == refPoints.dimension(0) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): size of input and output point arrays does not match each other.");
    }
#endif
    
    // Apply the parametrization map to every point in parameter domain
    const size_t numPts  = static_cast<size_t>(ortPoints.dimension(0));
    const auto key = cellTopo.getBaseCellTopologyData()->key ;
    switch (key) {
    case shards::Line<>::key : {
      for (size_t pt=0;pt<numPts;++pt) 
        getModifiedLinePoint(ortPoints(pt, 0),
                             refPoints(pt, 0), 
                             cellOrt);
      break;
    }
    case shards::Triangle<>::key : {
      for (size_t pt=0;pt<numPts;++pt) 
        getModifiedTrianglePoint(ortPoints(pt, 0), ortPoints(pt, 1),
                                 refPoints(pt, 0), refPoints(pt, 1), 
                                 cellOrt);
      break;
    }
    case shards::Quadrilateral<>::key : {
      for (size_t pt=0;pt<numPts;++pt) 
        getModifiedQuadrilateralPoint(ortPoints(pt, 0), ortPoints(pt, 1),
                                      refPoints(pt, 0), refPoints(pt, 1), 
                                      cellOrt);
      break;
    }
    default: 
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::invalid_argument,
                                 ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): Invalid cell topology." );
    }
  }
}
#endif
  
#endif
