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


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

//#include "Teuchos_LAPACK.hpp"
namespace Intrepid2 {

  namespace Impl {

    // ------------------------------------------------------------------------------------
    // Modified points according to orientations
    //
    //
    template<typename VT>
    inline
    void
    OrientationTools::
    getModifiedLinePoint(VT &ot,
                         const VT pt,
                         const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt && pt <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedLinePoint): "
                                "Input point is out of range [-1, 1].");
#endif
      
      switch (ort) {
      case 0: ot =  pt; break;
      case 1: ot = -pt; break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedLinePoint): "
                                  "Orientation is invalid (0--1)." );
      }
    }
    
    template<typename VT>
    inline
    void
    OrientationTools::getModifiedTrianglePoint(VT &ot0,
                                               VT &ot1,
                                               const VT pt0,
                                               const VT pt1,
                                               const ordinal_type ort) {
      const VT lambda[3] = { 1.0 - pt0 - pt1,
                                 pt0,
                                 pt1 };
      
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[0] && lambda[0] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[0]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[1] && lambda[1] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[1]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[2] && lambda[2] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): "
                                "Computed bicentric coordinate (lamba[2]) is out of range [0, 1].");
#endif
      
      switch (ort) {
      case 0: ot0 = lambda[1]; ot1 = lambda[2]; break;
      case 1: ot0 = lambda[0]; ot1 = lambda[1]; break;
      case 2: ot0 = lambda[2]; ot1 = lambda[0]; break;
        
      case 3: ot0 = lambda[2]; ot1 = lambda[1]; break;
      case 4: ot0 = lambda[0]; ot1 = lambda[2]; break;
      case 5: ot0 = lambda[1]; ot1 = lambda[0]; break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedTrianglePoint): " \
                                  "Orientation is invalid (0--5)." );
      }
    }

    template<typename VT>
    inline
    void
    OrientationTools::getModifiedQuadrilateralPoint(VT &ot0,
                                                    VT &ot1,
                                                    const VT pt0,
                                                    const VT pt1,
                                                    const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt0 && pt0 <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): " \
                                "Input point(0) is out of range [-1, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt1 && pt1 <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): " \
                                "Input point(1) is out of range [-1, 1].");
#endif
      
      const VT lambda[2][2] = { { pt0, -pt0 },
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
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedQuadrilateralPoint): " \
                                  "Orientation is invalid (0--7)." );
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // CoeffMatrix
  //
  //
  template<typename SpT, typename VT>
  OrientationTools<SpT,VT>::CoeffMatrix::
  CoeffMatrix()
    : _m(0), _n(0), _ap(), _aj(), _ax() { }

  template<typename SpT, typename VT>
  inline
  void
  OrientationTools<SpT,VT>::CoeffMatrix::
  createInternalArrays(const ordinal_type m,
                       const ordinal_type n,
                       const size_type nnz) {
    _m = m;
    _n = n;

    _ap = Kokkos::View<size_type*,   SpT>("OrientationTools::CoeffMatrix::RowPtrArray", m+1);
    _aj = Kokkos::View<ordinal_type*,SpT>("OrientationTools::CoeffMatrix::ColsArray",   nnz);
    _ax = Kokkos::View<VT*,          SpT>("OrientationTools::CoeffMAtrix::ValuesArray", nnz);
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  ordinal_type
  OrientationTools<SpT,VT>::CoeffMatrix::
  NumRows() const {
    return _m;
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  ordinal_type
  OrientationTools<SpT,VT>::CoeffMatrix::
  NumCols() const {
    return _n;
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  size_type
  OrientationTools<SpT,VT>::CoeffMatrix::
  RowPtr(const ordinal_type i) const {
    return _ap(i);
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  Kokkos::View<ordinal_type*,SpT>
  OrientationTools<SpT,VT>::CoeffMatrix::
  ColsInRow(const ordinal_type i) const {
    return Kokkos::subview(_aj, Kokkos::pair<ordinal_type,ordinal_type>(_ap(i), _ap(i+1)));
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  Kokkos::View<VT*,SpT>
  OrientationTools<SpT,VT>::CoeffMatrix::
  ValuesInRow(const ordinal_type i) const {
    return Kokkos::subview(_ax, Kokkos::pair<ordinal_type,ordinal_type>(_ap(i), _ap(i+1)));
  }

  template<typename SpT, typename VT>
  KOKKOS_INLINE_FUNCTION
  ordinal_type
  OrientationTools<SpT,VT>::CoeffMatrix::
  NumNonZerosInRow(const ordinal_type i) const {
    return (_ap(i+1) - _ap(i));
  }














































  template<typename SpT,
           typename VT>
  inline 
  void 
  OrientationTools<SpT,VT>::
  getEdgeCoeffMatrix_HGRAD(OrientationTools<SpT,VT>::CoeffMatrix &C,
                           const Basis<SpT,VT,VT> lineBasis,
                           const Basis<SpT,VT,VT> cellBasis,
                           const ordinal_type edgeId,
                           const ordinal_type edgeOrt) {
    typedef typename OrientationTools<Scalar>::CoeffMatrix CoeffMatrixType;

    // check lookup table
    bool found = false;
    if (found) {
      // C = foundCoeffMatrix'
    } else {

      // populate points on a line and map to subcell
      shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      shards::CellTopology lineTopo = lineBasis.getBaseCellTopology();

      const unsigned ordinal_type cellDim  = cellTopo.getDimension();
      const unsigned ordinal_type lineDim  = lineTopo.getDimension();
      const unsigned ordinal_type degree = cellBasis.getDegree();

      const unsigned ordinal_type numCellBasis = cellBasis.getCardinality();
      const unsigned ordinal_type numLineBasis = lineBasis.getCardinality();

      const ordinal_type ordEdge = cellBasis.getDofOrdinal(lineDim, edgeId, 0);
      const unsigned ordinal_type ndofEdge = cellBasis.getDofTag(ordEdge)[3];

      // reference points between (-1 , 1)
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !(ndofEdge == PointTools::getLatticeSize<Scalar>(lineTopo, degree, 1)),
                                  std::logic_error,
                                  ">>> ERROR (Intrepid::OrientationTools::getEdgeCoeffMatrix_HGRAD): " \
                                  "The number of DOFs does not match to the number of collocation points.");
#endif
      ArrayType refPtsLine(ndofEdge, lineDim);
      PointTools::getLattice<Scalar>(refPtsLine,
                                     lineTopo,
                                     degree,
                                     1,
                                     POINTTYPE_EQUISPACED);

      // modified points with orientation
      ArrayType ortPtsLine(ndofEdge, lineDim);
      mapToModifiedReference(ortPtsLine,
                             refPtsLine,
                             lineTopo,
                             edgeOrt);

      // map to reference coordinates
      ArrayType refPtsCell(ndofEdge, cellDim);
      CellTools<Scalar>::mapToReferenceSubcell(refPtsCell,
                                               refPtsLine,
                                               lineDim,
                                               edgeId,
                                               cellTopo);

      // temporary storage to evaluate vanila basis on reference points
      // basis is not reordered yet
      ArrayType tmpValues(numCellBasis, ndofEdge);

      // reorder values by topology
      ArrayType refValues(numCellBasis, ndofEdge);
      cellBasis.getValues(tmpValues, refPtsCell, OPERATOR_VALUE);
      getBasisFunctionsByTopology(refValues,
                                  tmpValues,
                                  cellBasis);

      // reorder values by topology
      ArrayType outValues(numLineBasis, ndofEdge);
      lineBasis.getValues(tmpValues, ortPtsLine, OPERATOR_VALUE);
      getBasisFunctionsByTopology(outValues,
                                  tmpValues,
                                  lineBasis);

      // compute offset
      unsigned ordinal_type offCell = 0, offLine = 2;
      {
        const unsigned ordinal_type numVert = cellTopo.getVertexCount();
        for (ordinal_type i=0;i<numVert;++i) {
          const unsigned ordinal_type ord = cellBasis.getDofOrdinal(0, i, 0);
          offCell += cellBasis.getDofTag(ord)[3]; // ndof of this vertex
        }
        for (ordinal_type i=0;i<edgeId;++i) {
          const unsigned ordinal_type ord = cellBasis.getDofOrdinal(1, i, 0);
          offCell += cellBasis.getDofTag(ord)[3]; // ndof of this edge
        }
      }

      // construct collocation matrix
      DenseMatrixType
        refMat(ndofEdge, ndofEdge),
        ortMat(ndofEdge, ndofEdge),
        pivVec(ndofEdge, 1);

      // transpose the matrix
      for (ordinal_type i=0;i<ndofEdge;++i) {
        for (ordinal_type j=0;j<ndofEdge;++j) {
          refMat.Value(j,i) = refValues(i+offCell, j);
          ortMat.Value(j,i) = outValues(i+offLine, j);
        }
      }

      // solve the system
      Teuchos::LAPACK<ordinal_type,Scalar> lapack;
      ordinal_type info = 0;
      lapack.GESV(ndofEdge, ndofEdge,
                  refMat.ValuePtr(),
                  refMat.ColStride(),
                  (ordinal_type*)pivVec.ValuePtr(),
                  ortMat.ValuePtr(),
                  ortMat.ColStride(),
                  &info);

      if (info) {
        std::stringstream ss;
        ss << ">>> ERROR (Intrepid::OrientationTools::getEdgeCoeffMatrix_HGRAD): "
           << "LAPACK return with error code: "
           << info;
        INTREPID2_TEST_FOR_ABORT( true, std::runtime_error, ss.str() );
      }

      CoeffMatrixType R;
      const bool transpose = reverse;

      // sparcify
      R.import(ortMat, transpose);

      // done!!
      C = R;
    }
  }

  template<typename SpT>
  template<class ArrayType>
  void
  OrientationTools<Scalar>::getTriangleCoeffMatrix_HGRAD(OrientationTools<Scalar>::CoeffMatrix & C,
                                                         const Basis<Scalar,ArrayType> &         faceBasis,
                                                         const Basis<Scalar,ArrayType> &         cellBasis,
                                                         const ordinal_type                               faceId,
                                                         const ordinal_type                               faceOrt) {
    typedef typename OrientationTools<Scalar>::DenseMatrix DenseMatrixType;
    typedef typename OrientationTools<Scalar>::CoeffMatrix CoeffMatrixType;

    // check lookup table
    bool found = false;
    if (found) {
      // C = foundCoeffMatrix'
    } else {
      // populate points on a line and map to subcell
      shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      shards::CellTopology faceTopo = faceBasis.getBaseCellTopology();

      // if the face is left-handed system, the orientation should be re-enumerated
      const ordinal_type leftHanded = cellTopo.getNodeMap(2, faceId, 1) > cellTopo.getNodeMap(2, faceId, 2);
      const ordinal_type leftOrt[] = { 0, 2, 1, 3, 5, 4 };
      const ordinal_type ort = (leftHanded ? leftOrt[faceOrt] : faceOrt);

      const unsigned ordinal_type cellDim  = cellTopo.getDimension();
      const unsigned ordinal_type faceDim  = faceTopo.getDimension();
      const unsigned ordinal_type degree = cellBasis.getDegree();

      const unsigned ordinal_type numCellBasis = cellBasis.getCardinality();
      const unsigned ordinal_type numFaceBasis = faceBasis.getCardinality();

      const ordinal_type ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
      const unsigned ordinal_type ndofFace = cellBasis.getDofTag(ordFace)[3];

      // reference points in triangle
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !(ndofFace == PointTools::getLatticeSize<Scalar>(faceTopo, degree, 1)),
                                  std::logic_error,
                                  ">>> ERROR (Intrepid::OrientationTools::getTriangleCoeffMatrix_HGRAD): " \
                                  "The number of DOFs does not match to the number of collocation points.");
#endif
      ArrayType refPtsFace(ndofFace, faceDim);
      PointTools::getLattice<Scalar>(refPtsFace,
                                     faceTopo,
                                     degree,
                                     1,
                                     POINTTYPE_EQUISPACED);

      // modified points with orientation
      ArrayType ortPtsFace(ndofFace, faceDim);
      mapToModifiedReference(ortPtsFace,
                             refPtsFace,
                             faceTopo,
                             ort);

      // map to reference coordinates
      ArrayType refPtsCell(ndofFace, cellDim);
      CellTools<Scalar>::mapToReferenceSubcell(refPtsCell,
                                               refPtsFace,
                                               faceDim,
                                               faceId,
                                               cellTopo);

      // temporary storage to evaluate vanila basis on reference points
      // basis is not reordered yet
      ArrayType tmpValues(numCellBasis, ndofFace);

      // reorder values by topology
      ArrayType refValues(numCellBasis, ndofFace);
      cellBasis.getValues(tmpValues, refPtsCell, OPERATOR_VALUE);
      getBasisFunctionsByTopology(refValues,
                                  tmpValues,
                                  cellBasis);

      // reorder values by topology
      ArrayType outValues(numFaceBasis, ndofFace);
      faceBasis.getValues(tmpValues, ortPtsFace, OPERATOR_VALUE);
      getBasisFunctionsByTopology(outValues,
                                  tmpValues,
                                  faceBasis);

      // compute offset
      unsigned ordinal_type offCell = 0;
      {
        const unsigned ordinal_type numVert = cellTopo.getVertexCount();
        for (ordinal_type i=0;i<numVert;++i) {
          const unsigned ordinal_type ord = cellBasis.getDofOrdinal(0, i, 0);
          offCell += cellBasis.getDofTag(ord)[3]; // ndof of this vertex
        }
        const unsigned ordinal_type numEdge = cellTopo.getEdgeCount();
        for (ordinal_type i=0;i<numEdge;++i) {
          const unsigned ordinal_type ord = cellBasis.getDofOrdinal(1, i, 0);
          offCell += cellBasis.getDofTag(ord)[3]; // ndof of this edge
        }
        for (ordinal_type i=0;i<faceId;++i) {
          const unsigned ordinal_type ord = cellBasis.getDofOrdinal(2, i, 0);
          offCell += cellBasis.getDofTag(ord)[3]; // ndof of this face
        }
      }
      unsigned ordinal_type offFace = 0;
      {
        const unsigned ordinal_type numVert = faceTopo.getVertexCount();
        for (ordinal_type i=0;i<numVert;++i) {
          const unsigned ordinal_type ord = faceBasis.getDofOrdinal(0, i, 0);
          offFace += faceBasis.getDofTag(ord)[3]; // ndof of this vertex
        }
        const unsigned ordinal_type numEdge = faceTopo.getEdgeCount();
        for (ordinal_type i=0;i<numEdge;++i) {
          const unsigned ordinal_type ord = faceBasis.getDofOrdinal(1, i, 0);
          offFace += faceBasis.getDofTag(ord)[3]; // ndof of this edge
        }
      }

      // construct collocation matrix
      DenseMatrixType
        refMat(ndofFace, ndofFace),
        ortMat(ndofFace, ndofFace),
        pivVec(ndofFace, 1);

      // transpose the matrix
      for (ordinal_type i=0;i<ndofFace;++i) {
        for (ordinal_type j=0;j<ndofFace;++j) {
          refMat.Value(j,i) = refValues(i+offCell, j);
          ortMat.Value(j,i) = outValues(i+offFace, j);
        }
      }

      // solve the system
      Teuchos::LAPACK<ordinal_type,Scalar> lapack;
      ordinal_type info = 0;
      lapack.GESV(ndofFace, ndofFace,
                  refMat.ValuePtr(),
                  refMat.ColStride(),
                  (ordinal_type*)pivVec.ValuePtr(),
                  ortMat.ValuePtr(),
                  ortMat.ColStride(),
                  &info);

      if (info) {
        std::stringstream ss;
        ss << ">>> ERROR (Intrepid::OrientationTools::getTriangleCoeffMatrix_HGRAD): "
           << "LAPACK return with error code: "
           << info;
        INTREPID2_TEST_FOR_ABORT( true, std::runtime_error, ss.str() );
      }

      CoeffMatrixType R;
      const bool transpose = reverse;

      // sparcify
      R.import(ortMat, transpose);

      // done!!
      C = R;
    }
  }

  template<typename SpT>
  template<class ArrayType>
  void OrientationTools<Scalar>::applyCoeffMatrix(ArrayType &                                   outValues,
                                                  const ArrayType &                             refValues,
                                                  const OrientationTools<Scalar>::CoeffMatrix & C,
                                                  const unsigned ordinal_type                            offset,
                                                  const unsigned ordinal_type                            numDofs) {
    const unsigned ordinal_type numPts = refValues.dimension(1);
    if (C.NumRows() && C.NumCols()) {
      const ordinal_type rank = refValues.rank();
      switch (rank) {
      case 2: {
        for (auto i=0;i<numDofs;++i) {
          const auto nnz      = C.NumNonZerosInRow(i);
          const auto colsPtr  = C.ColsInRow(i);
          const auto valsPtr  = C.ValuesInRow(i);
          const auto ii = i + offset;

          // sparse mat-vec
          for (auto j=0;j<numPts;++j) {
            Scalar temp = 0.0;
            for (ordinal_type p=0;p<nnz;++p)
              temp += valsPtr[p]*refValues(colsPtr[p] + offset, j);
            outValues(ii, j) = temp;
          }
        }
        break;
      }
      case 3: {
        const auto dimVal = refValues.dimension(2);
        for (auto i=0;i<numDofs;++i) {
          const auto nnz     = C.NumNonZerosInRow(i);
          const auto colsPtr = C.ColsInRow(i);
          const auto valsPtr = C.ValuesInRow(i);
          const auto ii = i + offset;

          // sparse mat-vec
          for (auto j=0;j<numPts;++j)
            for (auto k=0;k<dimVal;++k) {
              Scalar temp = 0.0;
              for (ordinal_type p=0;p<nnz;++p)
                temp += valsPtr[p]*refValues(colsPtr[p] + offset, j, k);
              outValues(ii, j, k) = temp;
            }
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::applyCoeffMatrix): " \
                                    "The rank of refValues is not 2 or 3.");
      }
      }
    } else {
      copyBasisValues(outValues,
                      refValues,
                      offset, numDofs,
                      0,      numPts);
    }
  }

  template<typename SpT>
  template<class ArrayType>
  void OrientationTools<Scalar>::copyBasisValues(ArrayType &        outValues,
                                                 const ArrayType &  refValues,
                                                 const unsigned ordinal_type offRows, const unsigned ordinal_type numRows,
                                                 const unsigned ordinal_type offCols, const unsigned ordinal_type numCols) {
    if (numRows && numCols) {
      const ordinal_type rank = refValues.rank();
      switch (rank) {
      case 2: {
        for (auto i=0;i<numRows;++i)
          for (auto j=0;j<numCols;++j)
            outValues(i + offRows, j + offCols) = refValues(i + offRows, j + offCols);
        break;
      }
      case 3: {
        const ordinal_type dimVal = refValues.dimension(2);
        for (auto i=0;i<numRows;++i)
          for (auto j=0;j<numCols;++j) {
            const auto ii = i + offRows, jj = j + offCols;
            for (auto k=0;k<dimVal;++k)
              outValues(ii, jj, k) = refValues(ii, jj, k);
          }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::copyBasis): " \
                                    "The rank of refValues is not 2 or 3.");
      }
      }
    }
  }


  // ------------------------------------------------------------------------------------
  // Public interface
  //
  //
  template<typename SpT>
  template<class ArrayType>
  bool OrientationTools<Scalar>::isLeftHandedCell(const ArrayType & pts) {
    // From all the tests, nodes seems to be fed as 1 dimensional array
    // with 1 x npts x ndim
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( pts.dimension(0) != 1, std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::isLeftHandedCell): " \
                                "Node point array is supposed to have 1 dimensional array.");
#endif
    const ordinal_type dim = pts.dimension(2);
    Scalar det = 0.0;
    switch (dim) {
    case 2: {
      // need 3 points (origin, x end point, y end point)
      const double v[2][2] = { { pts(0,1,0) - pts(0,0,0), pts(0,1,1) - pts(0,0,1) },
                               { pts(0,2,0) - pts(0,0,0), pts(0,2,1) - pts(0,0,1) } };

      det = (v[0][0]*v[1][1] - v[1][0]*v[0][1]);
      break;
    }
    case 3: {
      // need 4 points (origin, x end point, y end point, z end point)
      const double v[3][3] = { { pts(0,1,0) - pts(0,0,0), pts(0,1,1) - pts(0,0,1), pts(0,1,2) - pts(0,0,2) },
                               { pts(0,2,0) - pts(0,0,0), pts(0,2,1) - pts(0,0,1), pts(0,2,2) - pts(0,0,2) },
                               { pts(0,3,0) - pts(0,0,0), pts(0,3,1) - pts(0,0,1), pts(0,3,2) - pts(0,0,2) } };

      det = (v[0][0] * v[1][1] * v[2][2] +
             v[0][1] * v[1][2] * v[2][0] +
             v[0][2] * v[1][0] * v[2][1] -
             v[0][2] * v[1][1] * v[2][0] -
             v[0][0] * v[1][2] * v[2][1] -
             v[0][1] * v[1][0] * v[2][2]);
      break;
    }
    default:{
      INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
                                  ">>> ERROR (Intrepid::Orientation::setLeftHandedFlag): " \
                                  "Dimension of points must be 2 or 3");
    }
    }
    return (det < 0.0);
  }

  template<typename SpT>
  template<class ArrayType>
  void OrientationTools<Scalar>::mapToModifiedReference(ArrayType &                   outPoints,
                                                        const ArrayType &             refPoints,
                                                        const shards::CellTopology &  cellTopo,
                                                        const ordinal_type                     cellOrt) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      const ordinal_type cellDim = cellTopo.getDimension();
      INTREPID2_TEST_FOR_ABORT( !(hasReferenceCell(cellTopo) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): " \
                                  "The specified cell topology does not have a reference cell.");

      INTREPID2_TEST_FOR_ABORT( !( (1 <= cellDim) && (cellDim <= 2 ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): " \
                                  "Method defined only for 1 and 2-dimensional subcells.");

      INTREPID2_TEST_FOR_ABORT( !( outPoints.dimension(0) == refPoints.dimension(0) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): " \
                                  "Size of input and output point arrays does not match each other.");
    }
#endif

    // Apply the parametrization map to every point in parameter domain
    const size_type numPts  = static_cast<size_t>(outPoints.dimension(0));
    const auto key = cellTopo.getBaseCellTopologyData()->key ;
    switch (key) {
    case shards::Line<>::key : {
      for (size_type pt=0;pt<numPts;++pt)
        getModifiedLinePoint(outPoints(pt, 0),
                             refPoints(pt, 0),
                             cellOrt);
      break;
    }
    case shards::Triangle<>::key : {
      for (size_type pt=0;pt<numPts;++pt)
        getModifiedTrianglePoint(outPoints(pt, 0), outPoints(pt, 1),
                                 refPoints(pt, 0), refPoints(pt, 1),
                                 cellOrt);
      break;
    }
    case shards::Quadrilateral<>::key : {
      for (size_type pt=0;pt<numPts;++pt)
        getModifiedQuadrilateralPoint(outPoints(pt, 0), outPoints(pt, 1),
                                      refPoints(pt, 0), refPoints(pt, 1),
                                      cellOrt);
      break;
    }
    default:
      INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): " \
                                  "Invalid cell topology." );
    }
  }

  template<typename SpT>
  template<class ArrayType>
  void OrientationTools<Scalar>::getBasisFunctionsByTopology(ArrayType &                     outValues,
                                                             const ArrayType &               refValues,
                                                             const Basis<Scalar,ArrayType> & basis) {
    // cardinality
    const unsigned ordinal_type numBasis = basis.getCardinality();
    const unsigned ordinal_type numPts   = refValues.dimension(1);

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !( numBasis <= outValues.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getBasisFunctionsByTopology): " \
                                "Basis cardinality is bigger than outValues dimension(0).");

    INTREPID2_TEST_FOR_ABORT( !( numBasis <= refValues.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getBasisFunctionsByTopology): " \
                                "Basis cardinality is bigger than refValues dimension(0).");
#endif

    // topology
    shards::CellTopology cellTopo = basis.getBaseCellTopology();

    const unsigned ordinal_type numVert = cellTopo.getVertexCount();
    const unsigned ordinal_type numEdge = cellTopo.getEdgeCount();
    const unsigned ordinal_type numFace = cellTopo.getFaceCount();

    // offset storage for topological nodes
    ordinal_type offVert[9] = {}, offEdge[13] = {}, offFace[7] = {}, offIntr = 0, numBasisCounted = 0;

    // loop over tags
    if (numBasisCounted < numBasis) {
      for (ordinal_type i=0;i<numVert;++i) {
        const unsigned ordinal_type ord = basis.getDofOrdinal(0, i, 0);
        const unsigned ordinal_type ndof = basis.getDofTag( ord )[3];
        numBasisCounted += ndof;
        offVert[i+1] = ndof;
      }
    }
    if (numBasisCounted < numBasis) {
      for (ordinal_type i=0;i<numEdge;++i) {
        const unsigned ordinal_type ord = basis.getDofOrdinal(1, i, 0);
        const unsigned ordinal_type ndof = basis.getDofTag( ord )[3];
        numBasisCounted += ndof;
        offEdge[i+1] = ndof;
      }
    }
    if (numBasisCounted < numBasis) {
      for (ordinal_type i=0;i<numFace;++i) {
        const unsigned ordinal_type ord = basis.getDofOrdinal(2, i, 0);
        const unsigned ordinal_type ndof = basis.getDofTag( ord )[3];
        numBasisCounted += ndof;
        offFace[i+1] = ndof;
      }
    }

    // prefixsum
    offVert[0] = 0;
    for (ordinal_type i=0;i<numVert;++i)
      offVert[i+1] += offVert[i];

    offEdge[0] = offVert[numVert];
    for (ordinal_type i=0;i<numEdge;++i)
      offEdge[i+1] += offEdge[i];

    offFace[0] = offEdge[numEdge];
    for (ordinal_type i=0;i<numFace;++i)
      offFace[i+1] += offFace[i];

    // for 2D, offIntr is same as offFace[0]
    offIntr = offFace[numFace];

    // group basis by topology : vertex, edge, face, interior
    for (ordinal_type i=0;i<numBasis;++i) {
      const auto &tag = basis.getDofTag(i);
      const unsigned ordinal_type dimSubCell = tag[0];
      const unsigned ordinal_type ordSubCell = tag[1];
      const unsigned ordinal_type dofSubCell = tag[2];

      unsigned ordinal_type offset = 0;
      switch (dimSubCell) {
      case 0: offset = offVert[ordSubCell]; break;
      case 1: offset = offEdge[ordSubCell]; break;
      case 2: offset = offFace[ordSubCell]; break;
      case 3: offset = offIntr;             break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getBasisFunctionsByTopology): " \
                                    "Tag has invalid information.");
      }

      const ordinal_type rank = refValues.rank();
      switch (rank) {
      case 2: {
        for (ordinal_type j=0;j<numPts;++j)
          outValues(offset + dofSubCell, j) = refValues(i, j);
        break;
      }
      case 3: {
        const auto dimVal = refValues.dimension(2);
        for (ordinal_type j=0;j<numPts;++j)
          for (ordinal_type k=0;k<dimVal;++k)
            outValues(offset + dofSubCell, j, k) = refValues(i, j, k);
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::getBasisFunctionsByTopology): " \
                                    "The rank of refValues is not 2 or 3.");
      }
      }
    }
  }

  template<typename SpT>
  template<class ArrayType>
  void OrientationTools<Scalar>::getModifiedBasisFunctions(ArrayType &                        outValues,
                                                           const ArrayType &                  refValues,
                                                           const BasisSet<Scalar,ArrayType> & basisSet,
                                                           const Orientation                  ort) {
    const Basis<Scalar,ArrayType>& cellBasis = basisSet.getCellBasis();
    const EFunctionSpace space = basisSet.getFunctionSpace();

    const ordinal_type numBasis = cellBasis.getCardinality();
    const ordinal_type numPts = refValues.dimension(1);


#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !( numBasis <= outValues.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
                                "Basis cardinality is bigger than outValues dimension(0).");

    INTREPID2_TEST_FOR_ABORT( !( numBasis <= refValues.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions: " \
                                "Basis cardinality is bigger than refValues dimension(0).");

    INTREPID2_TEST_FOR_ABORT( !( refValues.dimension(0) <= outValues.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions: " \
                                "Dimension(0) in outValues is less than the dimension(0) in refValues.");

    INTREPID2_TEST_FOR_ABORT( !( numPts <= outValues.dimension(1) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
                                "Dimension(1) in refValues is greater than the number of points on outValues.");
#endif
    shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();

    // early return after simply copied
    const auto key = cellTopo.getBaseCellTopologyData()->key;
    if (key   == shards::Line<>::key ||
        space == FUNCTION_SPACE_HVOL) {
      const unsigned ordinal_type numRows = refValues.dimension(0);
      const unsigned ordinal_type numCols = refValues.dimension(1);
      copyBasisValues(outValues,
                      refValues,
                      0, numRows,
                      0, numCols);
      return;
    }

    // topology structure,
    // note that topology substructure dimension may be different from basis subdimension
    const unsigned ordinal_type cellDim = cellTopo.getDimension();
    const unsigned ordinal_type numVert = cellTopo.getVertexCount();
    const unsigned ordinal_type numEdge = cellTopo.getEdgeCount();
    const unsigned ordinal_type numFace = cellTopo.getFaceCount();
    const unsigned ordinal_type numIntr = 1;

    // vertex copy
    unsigned ordinal_type offset = 0;
    if (offset < numBasis) {
      for (ordinal_type i=0;i<numVert;++i) {
        const ordinal_type ord = cellBasis.getDofOrdinal(0, i, 0);
        const unsigned ordinal_type ndof= cellBasis.getDofTag(ord)[3];
        if (ndof)
          copyBasisValues(outValues,
                          refValues,
                          offset, ndof,
                          0,      numPts);
        offset += ndof;
      }
    }

    // edge rotation
    ordinal_type ortEdge[12];
    ort.getEdgeOrientation(ortEdge, numEdge);

    if (offset < numBasis) {
      for (ordinal_type i=0;i<numEdge;++i) {
        const ordinal_type ord = cellBasis.getDofOrdinal(1, i, 0);
        const unsigned ordinal_type ndof = cellBasis.getDofTag(ord)[3];
        const Basis<Scalar,ArrayType>& lineBasis = basisSet.getLineBasis();
        if (ndof) {
          typename OrientationTools<Scalar>::CoeffMatrix C;
          switch (space) {
          case FUNCTION_SPACE_HGRAD:
            OrientationTools<Scalar>::getEdgeCoeffMatrix_HGRAD(C,
                                                               lineBasis,
                                                               cellBasis,
                                                               i,
                                                               ortEdge[i]);
            break;
          case FUNCTION_SPACE_HCURL:
            break;
          case FUNCTION_SPACE_HDIV:
            break;
          default:
            INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
                                        ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
                                        "Functions space is invalid.");
          }

          if (verbose) {
            *verboseStreamPtr << " -- Edge ID : " << i << ", ndof = " << ndof << ", " << EFunctionSpaceToString(space) << "\n";
            C.showMe(*verboseStreamPtr);
          }
          applyCoeffMatrix(outValues,
                           refValues,
                           C, offset, ndof);

        }
        offset += ndof;
      }
    }

    // face rotation
    ordinal_type ortFace[6];
    ort.getFaceOrientation(ortFace, numFace);

    if (offset < numBasis) {
      for (ordinal_type i=0;i<numFace;++i) {
        const ordinal_type ord = cellBasis.getDofOrdinal(2, i, 0);
        const unsigned ordinal_type ndof = cellBasis.getDofTag(ord)[3];
        if (ndof) {
          typename OrientationTools<Scalar>::CoeffMatrix C;
          switch (space) {
          case FUNCTION_SPACE_HGRAD: {
            shards::CellTopology faceTopo(cellTopo.getCellTopologyData(2, i));
            const unsigned ordinal_type key = faceTopo.getBaseCellTopologyData()->key;
            if        (key == shards::Triangle<>::key) {
              const Basis<Scalar,ArrayType>& trigBasis = basisSet.getTriangleBasis();
              OrientationTools<Scalar>::getTriangleCoeffMatrix_HGRAD(C,
                                                                     trigBasis,
                                                                     cellBasis,
                                                                     i,
                                                                     ortFace[i]);
            } else if (key == shards::Quadrilateral<>::key) {
              const Basis<Scalar,ArrayType>& quadBasis = basisSet.getQuadrilateralBasis();
              OrientationTools<Scalar>::getQuadrilateralCoeffMatrix_HGRAD(C,
                                                                          quadBasis,
                                                                          cellBasis,
                                                                          i,
                                                                          ortFace[i]);
            } else {
              INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
                                          "Face topology is invalid.");
            }
            break;
          }
          case FUNCTION_SPACE_HCURL:
            break;
          case FUNCTION_SPACE_HDIV:
            break;
          default:
            INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
                                        ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
                                        "Functions space is invalid.");
          }
          if (verbose) {
            *verboseStreamPtr << " -- Face ID : " << i << ", ndof = " << ndof << ", " << EFunctionSpaceToString(space) << "\n";
            C.showMe(*verboseStreamPtr);
          }
          applyCoeffMatrix(outValues,
                           refValues,
                           C, offset, ndof);
        }
        offset += ndof;
      }
    }

    // interior copy
    if (offset < numBasis) {
      for (ordinal_type intrId=0;intrId<numIntr;++intrId) {
        const ordinal_type ord = cellBasis.getDofOrdinal(cellDim, intrId, 0);
        const unsigned ordinal_type ndof = cellBasis.getDofTag(ord)[3];
        if (ndof)
          copyBasisValues(outValues,
                          refValues,
                          offset, ndof,
                          0,      numPts);
        offset += ndof;
      }
    }
  }

  template<typename SpT>
  bool OrientationTools<Scalar>::verbose = false;

  // apply transposed map (reverse mapping for test only)
  template<typename SpT>
  bool OrientationTools<Scalar>::reverse = false;

  template<typename SpT>
  std::ostream* OrientationTools<Scalar>::verboseStreamPtr = &std::cout;

}

#endif




//   template<typename SpT>
//   void
//   OrientationTools<SpT>::CoeffMatrix::import(const OrientationTools<SpT>::DenseMatrix &b,
//                                                 const bool transpose) {
// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_ABORT( !( NumRows() == b.NumRows() && NumCols() == b.NumCols() ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::Orientation::CoeffMatrix::import): "
//                                 "Matrix dimensions are not matched");
// #endif
//     // count size
//     const SpT eps = 1.0e-8;
//     const ordinal_type nrows = (transpose ? b.NumCols() : b.NumRows());
//     const ordinal_type ncols = (transpose ? b.NumRows() : b.NumCols());
//     size_type nnz = b.countNumNonZeros(eps);
//     createInternalArrays(nrows, ncols, nnz);

//     // construct sparse array
//     nnz = 0;
//     for (ordinal_type i=0;i<nrows;++i) {
//       _ap(i) = nnz;
//       for (ordinal_type j=0;j<ncols;++j) {
//         const SpT val  = (transpose ? b.Value(j,i) : b.Value(i,j));
//         const SpT val2 = val*val;

//         // consider it as a nonzero entry
//         if (val2 > eps) {
//           _aj(nnz) = j;
//           _ax(nnz) = val;
//           ++nnz;
//         }
//       }
//     }
//     _ap(nrows) = nnz;
//   }

//   template<typename SpT>
//   std::ostream&
//   OrientationTools<SpT>::CoeffMatrix::showMe(std::ostream &os) const {
//     std::ofstream prec;
//     prec.copyfmt(os);

//     os.precision(3);

//     os << " -- OrientationTools::CoeffMatrix -- " << std::endl
//        << "    # of Rows          = " << _m << std::endl
//        << "    # of Cols          = " << _n << std::endl
//        << std::endl
//        << "    RowPtrArray length = " << _ap.dimension_0() << std::endl
//        << "    ColArray    length = " << _aj.dimension_0() << std::endl
//        << "    ValueArray  length = " << _ax.dimension_0() << std::endl
//        << std::endl;

//     const ordinal_type w = 10;
//     if (_ap.size() && _aj.size() && _ax.size()) {
//       os << std::setw(w) <<  "Row" << "  "
//          << std::setw(w) <<  "Col" << "  "
//          << std::setw(w) <<  "Val" << std::endl;
//       for (ordinal_type i=0;i<_m;++i) {
//         size_type jbegin = _ap[i], jend = _ap[i+1];
//         for (size_type j=jbegin;j<jend;++j) {
//           SpT val = _ax[j];
//           os << std::setw(w) <<      i << "  "
//              << std::setw(w) << _aj[j] << "  "
//              << std::setw(w) <<    val << std::endl;
//         }
//       }
//     }
//     os.copyfmt(prec);

//     return os;
//   }

  // template<class Scalar>
  // size_t
  // OrientationTools<Scalar>::DenseMatrix::countNumNonZeros(const Scalar epsilon) const {
  //   size_t nnz = 0;
  //   for (ordinal_type j=0;j<NumCols();++j) {
  //     for (ordinal_type i=0;i<NumRows();++i) {
  //       const Scalar val = Value(i,j);
  //       nnz += ((val*val) > epsilon);
  //     }
  //   }
  //   return nnz;
  // }

  // template<class Scalar>
  // std::ostream&
  // OrientationTools<Scalar>::DenseMatrix::showMe(std::ostream &os) const {
  //   std::ofstream prec;
  //   prec.copyfmt(os);

  //   os.precision(3);

  //   os << " -- OrientationTools::DenseMatrix -- " << std::endl
  //      << "    # of Rows              = " << _m << std::endl
  //      << "    # of Cols              = " << _n << std::endl
  //      << "    Col Stride             = " << _cs << std::endl
  //      << "    Row Stride             = " << _rs << std::endl
  //      << std::endl
  //      << "    ValueArray dimensions  = " << _a.dimension_0() << std::endl
  //      << std::endl;

  //   const ordinal_type w = 10;
  //   if (_a.size()) {
  //     for (ordinal_type i=0;i<_m;++i) {
  //       for (ordinal_type j=0;j<_n;++j) {
  //         const Scalar val = this->Value(i,j);
  //         os << std::setw(w) << val << "  ";
  //       }
  //       os << std::endl;
  //     }
  //   }
  //   os.copyfmt(prec);

  //   return os;
  // }

