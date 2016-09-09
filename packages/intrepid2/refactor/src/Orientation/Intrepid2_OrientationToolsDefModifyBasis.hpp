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
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {
  
  template<typename SpT>
  template<typename ptViewType>
  KOKKOS_INLINE_FUNCTION
  bool OrientationTools<SpT>::
  isLeftHandedCell(const ptViewType pts) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( pts.rank() != 2,  // npts x ndim
                              ">>> ERROR (Intrepid::OrientationTools::isLeftHandedCell): " \
                              "Point array is supposed to have rank 2.");
#endif
    typedef typename ptViewType::value_type value_type;
    
    const auto dim = pts.dimension(1);
    value_type det = 0.0;
    switch (dim) {
    case 2: {
      // need 3 points (origin, x end point, y end point)
      const value_type v[2][2] = { { pts(1,0) - pts(0,0), pts(1,1) - pts(0,1) },
                                   { pts(2,0) - pts(0,0), pts(2,1) - pts(0,1) } };
      
      det = (v[0][0]*v[1][1] - v[1][0]*v[0][1]);
      break;
    }
    case 3: {
      // need 4 points (origin, x end point, y end point, z end point)
      const value_type v[3][3] = { { pts(1,0) - pts(0,0), pts(1,1) - pts(0,1), pts(1,2) - pts(0,2) },
                                   { pts(2,0) - pts(0,0), pts(2,1) - pts(0,1), pts(2,2) - pts(0,2) },
                                   { pts(3,0) - pts(0,0), pts(3,1) - pts(0,1), pts(3,2) - pts(0,2) } };
      
      det = (v[0][0] * v[1][1] * v[2][2] +
             v[0][1] * v[1][2] * v[2][0] +
             v[0][2] * v[1][0] * v[2][1] -
             v[0][2] * v[1][1] * v[2][0] -
             v[0][0] * v[1][2] * v[2][1] -
             v[0][1] * v[1][0] * v[2][2]);
      break;
    }
    default:{
      INTREPID2_TEST_FOR_ABORT( true, 
                                ">>> ERROR (Intrepid::Orientation::isLeftHandedCell): " \
                                "Dimension of points must be 2 or 3");
    }
    }
    return (det < 0.0);
  }
  
  template<typename SpT>
  template<typename elemOrtValueType, class ...elemOrtProperties,
           typename elemNodeValueType, class ...elemNodeProperties>
  void
  OrientationTools<SpT>::
  getOrientation(/**/  Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                 const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                 const shards::CellTopology cellTopo) {
    // small meta data modification and it uses shards; let's do this on host
    typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
    auto elemOrtsHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemOrts);
    auto elemNodesHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemNodes);
    
    const ordinal_type numCells = elemNodes.dimension(0);
    for (auto cell=0;cell<numCells;++cell) {
      const auto nodes = Kokkos::subview(elemNodesHost, cell, Kokkos::ALL());
      elemOrtsHost(cell) = Orientation::getOrientation(cellTopo, nodes);
    }
    
    Kokkos::deep_copy(elemOrts, elemOrtsHost);
  }
}

#endif


//   template<typename SpT>
//   template<class ArrayType>
//   void OrientationTools<Scalar>::applyCoeffMatrix(ArrayType &                                   outValues,
//                                                   const ArrayType &                             refValues,
//                                                   const OrientationTools<Scalar>::CoeffMatrix & C,
//                                                   const unsigned int                            offset,
//                                                   const unsigned int                            numDofs) {
//     const unsigned int numPts = refValues.dimension(1);
//     if (C.NumRows() && C.NumCols()) {
//       const int rank = refValues.rank();
//       switch (rank) {
//       case 2: {
//         for (auto i=0;i<numDofs;++i) {
//           const auto nnz      = C.NumNonZerosInRow(i);
//           const auto colsPtr  = C.ColsInRow(i);
//           const auto valsPtr  = C.ValuesInRow(i);
//           const auto ii = i + offset;

//           // sparse mat-vec
//           for (auto j=0;j<numPts;++j) {
//             Scalar temp = 0.0;
//             for (int p=0;p<nnz;++p)
//               temp += valsPtr[p]*refValues(colsPtr[p] + offset, j);
//             outValues(ii, j) = temp;
//           }
//         }
//         break;
//       }
//       case 3: {
//         const auto dimVal = refValues.dimension(2);
//         for (auto i=0;i<numDofs;++i) {
//           const auto nnz     = C.NumNonZerosInRow(i);
//           const auto colsPtr = C.ColsInRow(i);
//           const auto valsPtr = C.ValuesInRow(i);
//           const auto ii = i + offset;

//           // sparse mat-vec
//           for (auto j=0;j<numPts;++j)
//             for (auto k=0;k<dimVal;++k) {
//               Scalar temp = 0.0;
//               for (int p=0;p<nnz;++p)
//                 temp += valsPtr[p]*refValues(colsPtr[p] + offset, j, k);
//               outValues(ii, j, k) = temp;
//             }
//         }
//         break;
//       }
//       default: {
//         INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
//                                     ">>> ERROR (Intrepid::OrientationTools::applyCoeffMatrix): " \
//                                     "The rank of refValues is not 2 or 3.");
//       }
//       }
//     } else {
//       copyBasisValues(outValues,
//                       refValues,
//                       offset, numDofs,
//                       0,      numPts);
//     }
//   }

//   template<typename SpT>
//   template<class ArrayType>
//   void OrientationTools<Scalar>::copyBasisValues(ArrayType &        outValues,
//                                                  const ArrayType &  refValues,
//                                                  const unsigned int offRows, const unsigned int numRows,
//                                                  const unsigned int offCols, const unsigned int numCols) {
//     if (numRows && numCols) {
//       const int rank = refValues.rank();
//       switch (rank) {
//       case 2: {
//         for (auto i=0;i<numRows;++i)
//           for (auto j=0;j<numCols;++j)
//             outValues(i + offRows, j + offCols) = refValues(i + offRows, j + offCols);
//         break;
//       }
//       case 3: {
//         const int dimVal = refValues.dimension(2);
//         for (auto i=0;i<numRows;++i)
//           for (auto j=0;j<numCols;++j) {
//             const auto ii = i + offRows, jj = j + offCols;
//             for (auto k=0;k<dimVal;++k)
//               outValues(ii, jj, k) = refValues(ii, jj, k);
//           }
//         break;
//       }
//       default: {
//         INTREPID2_TEST_FOR_ABORT( true, std::invalid_argument,
//                                     ">>> ERROR (Intrepid::OrientationTools::copyBasis): " \
//                                     "The rank of refValues is not 2 or 3.");
//       }
//       }
//     }
//   }


//   // ------------------------------------------------------------------------------------
//   // Public interface
//   //
//   //
//   template<typename SpT>
//   template<class ArrayType>
//   void OrientationTools<Scalar>::getModifiedBasisFunctions(ArrayType &                        outValues,
//                                                            const ArrayType &                  refValues,
//                                                            const BasisSet<Scalar,ArrayType> & basisSet,
//                                                            const Orientation                  ort) {
//     const Basis<Scalar,ArrayType>& cellBasis = basisSet.getCellBasis();
//     const EFunctionSpace space = basisSet.getFunctionSpace();

//     const int numBasis = cellBasis.getCardinality();
//     const int numPts = refValues.dimension(1);


// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_ABORT( !( numBasis <= outValues.dimension(0) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
//                                 "Basis cardinality is bigger than outValues dimension(0).");

//     INTREPID2_TEST_FOR_ABORT( !( numBasis <= refValues.dimension(0) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions: " \
//                                 "Basis cardinality is bigger than refValues dimension(0).");

//     INTREPID2_TEST_FOR_ABORT( !( refValues.dimension(0) <= outValues.dimension(0) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions: " \
//                                 "Dimension(0) in outValues is less than the dimension(0) in refValues.");

//     INTREPID2_TEST_FOR_ABORT( !( numPts <= outValues.dimension(1) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
//                                 "Dimension(1) in refValues is greater than the number of points on outValues.");
// #endif
//     shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();

//     // early return after simply copied
//     const auto key = cellTopo.getBaseCellTopologyData()->key;
//     if (key   == shards::Line<>::key ||
//         space == FUNCTION_SPACE_HVOL) {
//       const unsigned int numRows = refValues.dimension(0);
//       const unsigned int numCols = refValues.dimension(1);
//       copyBasisValues(outValues,
//                       refValues,
//                       0, numRows,
//                       0, numCols);
//       return;
//     }

//     // topology structure,
//     // note that topology substructure dimension may be different from basis subdimension
//     const unsigned int cellDim = cellTopo.getDimension();
//     const unsigned int numVert = cellTopo.getVertexCount();
//     const unsigned int numEdge = cellTopo.getEdgeCount();
//     const unsigned int numFace = cellTopo.getFaceCount();
//     const unsigned int numIntr = 1;

//     // vertex copy
//     unsigned int offset = 0;
//     if (offset < numBasis) {
//       for (int i=0;i<numVert;++i) {
//         const int ord = cellBasis.getDofOrdinal(0, i, 0);
//         const unsigned int ndof= cellBasis.getDofTag(ord)[3];
//         if (ndof)
//           copyBasisValues(outValues,
//                           refValues,
//                           offset, ndof,
//                           0,      numPts);
//         offset += ndof;
//       }
//     }

//     // edge rotation
//     int ortEdge[12];
//     ort.getEdgeOrientation(ortEdge, numEdge);

//     if (offset < numBasis) {
//       for (int i=0;i<numEdge;++i) {
//         const int ord = cellBasis.getDofOrdinal(1, i, 0);
//         const unsigned int ndof = cellBasis.getDofTag(ord)[3];
//         const Basis<Scalar,ArrayType>& lineBasis = basisSet.getLineBasis();
//         if (ndof) {
//           typename OrientationTools<Scalar>::CoeffMatrix C;
//           switch (space) {
//           case FUNCTION_SPACE_HGRAD:
//             OrientationTools<Scalar>::getEdgeCoeffMatrix_HGRAD(C,
//                                                                lineBasis,
//                                                                cellBasis,
//                                                                i,
//                                                                ortEdge[i]);
//             break;
//           case FUNCTION_SPACE_HCURL:
//             break;
//           case FUNCTION_SPACE_HDIV:
//             break;
//           default:
//             INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
//                                         ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
//                                         "Functions space is invalid.");
//           }

//           if (verbose) {
//             *verboseStreamPtr << " -- Edge ID : " << i << ", ndof = " << ndof << ", " << EFunctionSpaceToString(space) << "\n";
//             C.showMe(*verboseStreamPtr);
//           }
//           applyCoeffMatrix(outValues,
//                            refValues,
//                            C, offset, ndof);

//         }
//         offset += ndof;
//       }
//     }

//     // face rotation
//     int ortFace[6];
//     ort.getFaceOrientation(ortFace, numFace);

//     if (offset < numBasis) {
//       for (int i=0;i<numFace;++i) {
//         const int ord = cellBasis.getDofOrdinal(2, i, 0);
//         const unsigned int ndof = cellBasis.getDofTag(ord)[3];
//         if (ndof) {
//           typename OrientationTools<Scalar>::CoeffMatrix C;
//           switch (space) {
//           case FUNCTION_SPACE_HGRAD: {
//             shards::CellTopology faceTopo(cellTopo.getCellTopologyData(2, i));
//             const unsigned int key = faceTopo.getBaseCellTopologyData()->key;
//             if        (key == shards::Triangle<>::key) {
//               const Basis<Scalar,ArrayType>& trigBasis = basisSet.getTriangleBasis();
//               OrientationTools<Scalar>::getTriangleCoeffMatrix_HGRAD(C,
//                                                                      trigBasis,
//                                                                      cellBasis,
//                                                                      i,
//                                                                      ortFace[i]);
//             } else if (key == shards::Quadrilateral<>::key) {
//               const Basis<Scalar,ArrayType>& quadBasis = basisSet.getQuadrilateralBasis();
//               OrientationTools<Scalar>::getQuadrilateralCoeffMatrix_HGRAD(C,
//                                                                           quadBasis,
//                                                                           cellBasis,
//                                                                           i,
//                                                                           ortFace[i]);
//             } else {
//               INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
//                                           ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
//                                           "Face topology is invalid.");
//             }
//             break;
//           }
//           case FUNCTION_SPACE_HCURL:
//             break;
//           case FUNCTION_SPACE_HDIV:
//             break;
//           default:
//             INTREPID2_TEST_FOR_ABORT( true, std::runtime_error,
//                                         ">>> ERROR (Intrepid::OrientationTools::getModifiedBasisFunctions): " \
//                                         "Functions space is invalid.");
//           }
//           if (verbose) {
//             *verboseStreamPtr << " -- Face ID : " << i << ", ndof = " << ndof << ", " << EFunctionSpaceToString(space) << "\n";
//             C.showMe(*verboseStreamPtr);
//           }
//           applyCoeffMatrix(outValues,
//                            refValues,
//                            C, offset, ndof);
//         }
//         offset += ndof;
//       }
//     }

//     // interior copy
//     if (offset < numBasis) {
//       for (int intrId=0;intrId<numIntr;++intrId) {
//         const int ord = cellBasis.getDofOrdinal(cellDim, intrId, 0);
//         const unsigned int ndof = cellBasis.getDofTag(ord)[3];
//         if (ndof)
//           copyBasisValues(outValues,
//                           refValues,
//                           offset, ndof,
//                           0,      numPts);
//         offset += ndof;
//       }
//     }
//   }

//   template<typename SpT>
//   bool OrientationTools<Scalar>::verbose = false;

//   // apply transposed map (reverse mapping for test only)
//   template<typename SpT>
//   bool OrientationTools<Scalar>::reverse = false;

//   template<typename SpT>
//   std::ostream* OrientationTools<Scalar>::verboseStreamPtr = &std::cout;





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
//     const int nrows = (transpose ? b.NumCols() : b.NumRows());
//     const int ncols = (transpose ? b.NumRows() : b.NumCols());
//     size_type nnz = b.countNumNonZeros(eps);
//     createInternalArrays(nrows, ncols, nnz);

//     // construct sparse array
//     nnz = 0;
//     for (int i=0;i<nrows;++i) {
//       _ap(i) = nnz;
//       for (int j=0;j<ncols;++j) {
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

//     const int w = 10;
//     if (_ap.size() && _aj.size() && _ax.size()) {
//       os << std::setw(w) <<  "Row" << "  "
//          << std::setw(w) <<  "Col" << "  "
//          << std::setw(w) <<  "Val" << std::endl;
//       for (int i=0;i<_m;++i) {
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
  //   for (int j=0;j<NumCols();++j) {
  //     for (int i=0;i<NumRows();++i) {
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

  //   const int w = 10;
  //   if (_a.size()) {
  //     for (int i=0;i<_m;++i) {
  //       for (int j=0;j<_n;++j) {
  //         const Scalar val = this->Value(i,j);
  //         os << std::setw(w) << val << "  ";
  //       }
  //       os << std::endl;
  //     }
  //   }
  //   os.copyfmt(prec);

  //   return os;
  // }

