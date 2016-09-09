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

/** \file   Intrepid_OrientationTools.hpp
    \brief  Header file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_ORIENTATIONTOOLS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_PointTools.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "Teuchos_LAPACK.hpp"

namespace Intrepid2 {

  namespace Impl {

    class OrientationTools {
    public:

      // -----------------------------------------------------------------------------
      // Point modification
      //
      //

      /** \brief  Computes modified point for line segment.
          
          \param  ot       [out] - modified point value
          \param  pt       [in]  - input point in [-1.0 , 1.0]
          \param  ort      [in]  - orientation number between 0 and 1
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedLinePoint(ValueType &ot,
                           const ValueType pt,
                           const ordinal_type ort);
      
      /** \brief  Computes modified point for triangle.
          
          \param  ot0      [out] - modified coordinate 0
          \param  ot1      [out] - modified coordinate 1
          \param  pt0      [out] - input coordinate 0
          \param  pt1      [out] - input coordinate 1
          \param  ort      [in]  - orientation number between 0 and 5
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedTrianglePoint(ValueType &ot0,
                               ValueType &ot1,
                               const ValueType pt0,
                               const ValueType pt1,
                               const ordinal_type ort);
      
      /** \brief  Computes modified point for quadrilateral.
          
          \param  ot0      [out] - modified coordinate 0
          \param  ot1      [out] - modified coordinate 1
          \param  pt0      [out] - input coordinate 0
          \param  pt1      [out] - input coordinate 1
          \param  ort      [in]  - orientation number between 0 and 7
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedQuadrilateralPoint(ValueType &ot0,
                                    ValueType &ot1,
                                    const ValueType pt0,
                                    const ValueType pt1,
                                    const ordinal_type ort);

      /** \brief  Computes modified parameterization maps of 1- and 2-subcells with orientation.
          
          \param  outPoints       [out] - rank-2 (P,D2) array with points in 1D or 2D modified domain with orientation
          \param  refPoints       [in]  - rank-2 (P,D2) array with points in 1D or 2D parameter domain
          \param  cellTopo        [in]  - cell topology of the parameterized domain (1- and 2-subcells)
          \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
      */
      template<typename outPointViewType,
               typename refPointViewType>
      inline
      static void 
      mapToModifiedReference(outPointViewType outPoints,
                             const refPointViewType refPoints,
                             const shards::CellTopology cellTopo,
                             const ordinal_type cellOrt = 0);

      // -----------------------------------------------------------------------------
      // Coefficient Matrix
      //
      //
      
      /** \brief  Compute coefficient matrix by collocating point values
          
          \param  lineBasis [in]  - line basis function 
          \param  cellBasis [in]  - cell basis function
          \param  edgeId    [in]  - edge Id in the cell topology
          \param  edgeOrt   [in]  - orientation number between 0 and 1

          \return rank 2 coefficient matrix
      */
      template<typename outputViewType,
               typename lineBasisType,
               typename cellBasisType>
      inline
      static void
      getEdgeCoeffMatrix_HGRAD(outputViewType &output,
                               const lineBasisType lineBasis,
                               const cellBasisType cellBasis,
                               const ordinal_type edgeId,
                               const ordinal_type edgeOrt);

      // /** \brief  Compute coefficient matrix by collocating point values
          
      //     \param  faceBasis [in]  - triangle face basis function 
      //     \param  cellBasis [in]  - cell basis function
      //     \param  faceId    [in]  - face Id in the cell topology
      //     \param  faceOrt   [in]  - orientation number between 0 and 5

      //     \return rank 2 coefficient matrix
      // */
      // template<typename ExecSpaceType,
      //          typename ValueType>
      // inline
      // static Kokkos::View<ValueType**,ExecSpaceType,Kokkos::LayoutStride>
      // getTriangleCoeffMatrix_HGRAD(const Basis<ExecSpaceType,ValueType,ValueType> faceBasis,
      //                              const Basis<ExecSpaceType,ValueType,ValueType> cellBasis,
      //                              const ordinal_type faceId,
      //                              const ordinal_type faceOrt);
      
      // /** \brief  Compute coefficient matrix by collocating point values
          
      //     \param  faceBasis [in]  - triangle face basis function 
      //     \param  cellBasis [in]  - cell basis function
      //     \param  faceId    [in]  - face Id in the cell topology
      //     \param  faceOrt   [in]  - orientation number between 0 and 7

      //     \return rank 2 coefficient matrix

      //     For simplicity, this one does not use tensor product space; 
      //     later It would be better if we use it. 
      // */
      // template<typename ExecSpaceType,
      //          typename ValueType>
      // inline
      // static Kokkos::View<ValueType**,ExecSpaceType,Kokkos::LayoutStride> 
      // getQuadrilateralCoeffMatrix_HGRAD(const Basis<ExecSpaceType,ValueType,ValueType> faceBasis,
      //                                   const Basis<ExecSpaceType,ValueType,ValueType> cellBasis,
      //                                   const ordinal_type faceId,
      //                                   const ordinal_type faceOrt);
    };
  }

  template<typename ExecSpaceType>
  class OrientationTools {
  public:
    // function space, order, subcell ordinal, orientation, matrix m x n
    typedef Kokkos::View<double******,ExecSpaceType> MatrixDataViewType;

    static MatrixDataViewType quadEdgeData;
    
  private:

    // space and order is fixed
    inline
    static void initQuadrilateral(Kokkos::View<double****,Kokkos::LayoutStride,ExecSpaceType> matData,
                                  const EFunctionSpace space,
                                  const ordinal_type order);

  public:
    inline 
    static void initialize(const shards::CellTopology cellTopo, 
                           const EFunctionSpace space,
                           const ordinal_type order);

    inline 
    static void finalize();

    // if an element is aligned left, it is an error.
    template<typename ptViewType>
    KOKKOS_INLINE_FUNCTION
    static bool 
    isLeftHandedCell(const ptViewType pts);

    // compute orientations of cells in a workset
    template<typename elemOrtValueType, class ...elemOrtProperties,
             typename elemNodeValueType, class ...elemNodeProperties>
    inline 
    static void
    getOrientation(/**/  Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                   const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                   const shards::CellTopology cellTopo);
        
    template<typename outValueValueType, class ...outValueProperties,
             typename refValueValueType, class ...refValueProperties,
             typename elemOrtValueType,  class ...elemOrtProperties,
             typename quadBasisType>
    inline
    static void 
    getModifiedHgradBasisQuadrilateral(/**/  Kokkos::DynRankView<outValueValueType,outValueProperties...> outValues,
                                       const Kokkos::DynRankView<refValueValueType,refValueProperties...> refValues,
                                       const Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                                       const quadBasisType quadBasis);
  };
  
  template<typename T> 
  typename OrientationTools<T>::MatrixDataViewType OrientationTools<T>::quadEdgeData;
}

// include templated function definitions
#include "Intrepid2_OrientationToolsDefModifyPoints.hpp"
#include "Intrepid2_OrientationToolsDefCoeffMatrix.hpp"
#include "Intrepid2_OrientationToolsDefMatrixData.hpp"
#include "Intrepid2_OrientationToolsDefModifyBasis.hpp"

#endif





  //   class CoeffMatrix {
  //   private:
  //     ordinal_type _m, _n;

  //     Kokkos::View<size_type*,ExecSpaceType> _ap;      //!< pointers to column index and values
  //     Kokkos::View<ordinal_type*,ExecSpaceType> _aj;   //!< column index compressed format
  //     Kokkos::View<double*,ExecSpaceType> _ax;         //!< values

  //     inline 
  //     void createInternalArrays(const ordinal_type m,
  //                               const ordinal_type n,
  //                               const size_type nnz) {
  //       _m = m;
  //       _n = n;
        
  //       _ap = Kokkos::View<size_type*,   SpT>("OrientationTools::CoeffMatrix::RowPtrArray", m+1);
  //       _aj = Kokkos::View<ordinal_type*,SpT>("OrientationTools::CoeffMatrix::ColsArray",   nnz);
  //       _ax = Kokkos::View<double*,      SpT>("OrientationTools::CoeffMAtrix::ValuesArray", nnz);
  //     }

  //   public:
  //     KOKKOS_INLINE_FUNCTION
  //     CoeffMatrix()
  //       : _m(0), _n(0), _ap(), _aj(), _ax() { }
      
  //     KOKKOS_INLINE_FUNCTION
  //     CoeffMatrix(const CoeffMatrix &b) = default;

  //     KOKKOS_INLINE_FUNCTION
  //     ordinal_type NumRows() const { 
  //       return _m; 
  //     }

  //     KOKKOS_INLINE_FUNCTION
  //     ordinal_type NumCols() const { 
  //       return _n; 
  //     }

  //     KOKKOS_INLINE_FUNCTION
  //     size_type RowPtr(const ordinal_type i) const { 
  //       return _ap(i); 
  //     }

  //     KOKKOS_INLINE_FUNCTION
  //     Kokkos::View<ordinal_type*,ExecSpaceType> ColsInRow(const ordinal_type i) const {
  //       return Kokkos::subview(_aj, Kokkos::pair<ordinal_type,ordinal_type>(_ap(i), _ap(i+1)));
  //     }

  //     KOKKOS_INLINE_FUNCTION
  //     Kokkos::View<double*,ExecSpaceType> ValuesInRow(const ordinal_type i) const {
  //       return Kokkos::subview(_ax, Kokkos::pair<ordinal_type,ordinal_type>(_ap(i), _ap(i+1)));
  //     }

  //     KOKKOS_INLINE_FUNCTION
  //     ordinal_type NumNonZerosInRow(const ordinal_type i) const {
  //       return (_ap(i+1) - _ap(i));
  //     }
  //   };

