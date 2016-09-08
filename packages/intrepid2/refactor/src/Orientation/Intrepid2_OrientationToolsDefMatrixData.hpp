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
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  template<typename SpT>
  void
  OrientationTools<SpT>::
  initQuadrilateral(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
                    const EFunctionSpace space,
                    const ordinal_type order) {

    switch (space) {
    case FUNCTION_SPACE_HGRAD: {
      Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
      Basis_HGRAD_QUAD_Cn_FEM<SpT> cellBasis(order);
      
      const ordinal_type numEdge = 4, numOrt = 2;
      for (auto edgeId=0;edgeId<numEdge;++edgeId)
        for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
          MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
          Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
        }
      break;
    }            
    case FUNCTION_SPACE_HCURL:
    case FUNCTION_SPACE_HDIV: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
                                    "Not yet implemented.");
      break;
    }
    case FUNCTION_SPACE_HVOL: {
      // do nothing
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
                                    "Invalid function space.");
      break;
    }
    }
  }

  template<typename SpT>
  void
  OrientationTools<SpT>::
  initTriangle(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
               const EFunctionSpace space,
               const ordinal_type order) {
    
    switch (space) {
    case FUNCTION_SPACE_HGRAD: {
      Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
      Basis_HGRAD_TRI_Cn_FEM <SpT> cellBasis(order);
      
      const ordinal_type numEdge = 3, numOrt = 2;
      for (auto edgeId=0;edgeId<numEdge;++edgeId)
        for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
          MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
          Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
        }
      break;
    }            
    case FUNCTION_SPACE_HCURL:
    case FUNCTION_SPACE_HDIV: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
                                    "Not yet implemented.");
      break;
    }
    case FUNCTION_SPACE_HVOL: {
      // do nothing
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
                                    "Invalid function space.");
      break;
    }
    }
  }

  template<typename SpT>
  void OrientationTools<SpT>::initialize(const shards::CellTopology cellTopo,
                                         const EFunctionSpace space,
                                         const ordinal_type order) {
    const auto key = cellTopo.getBaseCellTopologyData()->key;
    switch (key) {
    case shards::Triangle<>::key : {
      if (!trigMatrixData.span())
        trigMatrixData = Kokkos::View<CoeffMatrixType*****,SpT>("trigMatrixData", 
                                                                3,  // # of function space
                                                                Parameters::MaxOrder, // # of orders
                                                                1,  // # of subcell dims
                                                                3,  // # of edges
                                                                2); // # of orts
      auto MatrixData = Kokkos::subview(trigMatrixData, 
                                        space, order - 1, 
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()); 
      initTriangle(MatrixData, space, order);
      break;
    }
    case shards::Quadrilateral<>::key : {
      if (!quadMatrixData.span())
        quadMatrixData = Kokkos::View<CoeffMatrixType*****,SpT>("quadMatrixData", 
                                                                3,  // # of function space
                                                                Parameters::MaxOrder, // # of orders
                                                                1,  // # of subcell dims
                                                                4,  // # of edges
                                                                2); // # of orts
      auto MatrixData = Kokkos::subview(quadMatrixData, 
                                        space, order - 1, 
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()); 
      initQuadrilateral(MatrixData, space, order);
      break;
    }
    }
  }

  template<typename SpT>
  void OrientationTools<SpT>::finalize() {
    const Kokkos::View<CoeffMatrix***,ExecSpaceType> null;
    quadMatrixData = null;
    trigMatrixData = null;
  }

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





// template<typename SpT>
// void
// OrientationTools<SpT>::
// initHexahedron(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
//                const EFunctionSpace space,
//                const ordinal_type order) {
    
//   switch (space) {
//   case FUNCTION_SPACE_HGRAD: {
//     Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
//     Basis_HGRAD_QUAD_Cn_FEM<SpT> faceBasis(order);
//     Basis_HGRAD_HEXA_Cn_FEM<SpT> cellBasis(order);
      
//     {
//       const ordinal_type numEdge = 12, numOrt = 2;
//       for (auto edgeId=0;edgeId<numEdge;++edgeId)
//         for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
//           const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
//           MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
//         }
//     }
//     {
//       const ordinal_type numFace = 6, numOrt = 8;
//       for (auto faceId=0;faceId<numFace;++faceId)
//         for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
//           const auto C = Impl::OrientationTools::getQuadrilateralCoeffMatrix_HGRAD(faceBasis, cellBasis, faceId, faceOrt);
//           MatrixData(1, faceId, faceOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(1, faceId, faceOrt), C);
//         }
//     }
//     break;
//   }            
//   case FUNCTION_SPACE_HCURL:
//   case FUNCTION_SPACE_HDIV: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
//                                   "Not yet implemented.");
//     break;
//   }
//   case FUNCTION_SPACE_HVOL: {
//     // do nothing
//     break;
//   }
//   default: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
//                                   "Invalid function space.");
//     break;
//   }
//   }
// }

// template<typename SpT>
// void
// OrientationTools<SpT>::
// initTetrahedron(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
//                 const EFunctionSpace space,
//                 const ordinal_type order) {
    
//   switch (space) {
//   case FUNCTION_SPACE_HGRAD: {
//     Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
//     Basis_HGRAD_TRI_Cn_FEM <SpT> faceBasis(order);
//     Basis_HGRAD_TET_Cn_FEM <SpT> cellBasis(order);
      
//     {
//       const ordinal_type numEdge = 6, numOrt = 2;
//       for (auto edgeId=0;edgeId<numEdge;++edgeId)
//         for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
//           const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
//           MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
//         }
//     }
//     {
//       const ordinal_type numFace = 4, numOrt = 6;
//       for (auto faceId=0;faceId<numFace;++faceId)
//         for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
//           const auto C = Impl::OrientationTools::getTriangleCoeffMatrix_HGRAD(faceBasis, cellBasis, faceId, faceOrt);
//           MatrixData(1, faceId, faceOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(1, faceId, faceOrt), C);
//         }
//     }
//     break;
//   }            
//   case FUNCTION_SPACE_HCURL:
//   case FUNCTION_SPACE_HDIV: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
//                                   "Not yet implemented.");
//     break;
//   }
//   case FUNCTION_SPACE_HVOL: {
//     // do nothing
//     break;
//   }
//   default: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): " \
//                                   "Invalid function space.");
//     break;
//   }
//   }
// }
