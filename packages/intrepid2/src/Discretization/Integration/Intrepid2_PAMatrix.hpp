// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PAMatrix.hpp
    \brief  Header file for the Intrepid2::PAMatrix class; provides support for matrix partial assembly.
    \author Created by Nathan V. Roberts.
*/

#ifndef __INTREPID2_PAMATRIX_HPP__
#define __INTREPID2_PAMATRIX_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::PAMatrix
      \brief Provides support for structure-aware integration.
  */
  template<typename DeviceType = Kokkos::DefaultExecutionSpace::device_type,
           typename Scalar = double>
  class PAMatrix {
  public:
    using View1D = Kokkos::View<Scalar*,DeviceType>;
    struct OpSpec {
      View1D opView;
      int M;
      int N;
    };
    struct PointDataSpec {
      int C;
      int P;
      int a0;
      int b0;
      int aSpan;
      int bSpan;
      
      bool operator<(const PointDataSpec &otherSpec) const {
        if      (C < otherSpec.C) return true;
        else if (C > otherSpec.C) return false;
        
        if      (P < otherSpec.P) return true;
        else if (P > otherSpec.P) return false;
        
        if      (a0 < otherSpec.a0) return true;
        else if (a0 > otherSpec.a0) return false;
        
        if      (b0 < otherSpec.b0) return true;
        else if (b0 > otherSpec.b0) return false;
        
        if      (aSpan < otherSpec.aSpan) return true;
        else if (aSpan > otherSpec.aSpan) return false;
        
        if      (bSpan < otherSpec.bSpan) return true;
        else if (bSpan > otherSpec.bSpan) return false;
        
        // all equal
        return false;
      }
      
      bool operator==(const PointDataSpec &otherSpec) const {
        return         (C == otherSpec.C) && (P == otherSpec.P)
            &&       (a0 == otherSpec.a0) && (b0 == otherSpec.b0)
            && (aSpan == otherSpec.aSpan) && (bSpan == otherSpec.bSpan);
      }
    };
    Data<Scalar,DeviceType> _composedWeightedTransform; // (C,P[,D1[,D2]]), used for general case
    std::map<PointDataSpec,View1D> _pointDataCache; // copies of appropriate slices of _composedWeightedTransform; will be regenerated when recomputePointData() is called.
    TensorData<Scalar,DeviceType> _cellMeasures; // (C,P); used for separable case
    TransformedBasisValues<Scalar,DeviceType> _basisValuesLeft, _basisValuesRight;
    ScalarView<Orientation,DeviceType> _orientations;
    static constexpr bool layoutLeft_ = true; // BLAS expects this
    int maxIntermediateSize_ = 0;
    
    using ComponentSequence = std::tuple<std::vector<OpSpec>, PointDataSpec, std::vector<OpSpec>, int, int, int, int>; // left, pointData, right, left offset, left output span, right offset, right input span
    std::vector<ComponentSequence> componentIntegralsToSum_;
    
    bool _separable = false; // separable means that we can perform integrals in reference space, and separately in each tensorial component dimension. (currently unused)
    
    PAMatrix()
    {}

    void init(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
              const TensorData<Scalar,DeviceType> cellMeasures,
              const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
              const ScalarView<Orientation,DeviceType> orientations);
    /** \brief   Constructs a <b>PAMatrix</b>  representing the contraction of \a <b>basisValuesLeft</b> against \a <b>basisValuesRight</b> containers on
                 point and space dimensions, weighting each point according to <b>cellMeasures</b>.

        \param  basisValuesLeft      [in] - Left input container, with logical shape (C,F1,P,D)
        \param  cellMeasures             [in] - Point weight container, with logical shape (C,P)
        \param  basisValuesRight    [in] - Right input container with logical shape (C,F2,P,D)
        \param  orientations             [in] - orientations container, with shape (C) (optional)

        On construction, computes (if it will be needed) a <b>composedTransform</b> object with logical shape (C,P), (C,P,D), or (C,P,D,D) that stores (det J) M^T_L M_R, where det J represents <b>cellMeasures</b> and M_L and M_R represent the basis transformations for the left and right basis, respectively.
    */
    PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
             const TensorData<Scalar,DeviceType> cellMeasures,
             const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
             const ScalarView<Orientation,DeviceType> orientations = ScalarView<Orientation,DeviceType>());
    
    /** \brief   Constructs a <b>PAMatrix</b>  representing the contraction of \a <b>basisValues</b> against itself in
                 point and space dimensions, weighting each point according to <b>cellMeasures</b>.

        \param  basisValues               [in] - Transformed basis values input container, with logical shape (C,F,P,D)
        \param  cellMeasures             [in] - Point weight container, with logical shape (C,P)
        \param  orientations             [in] - orientations container, with shape (C)
        
        On construction, computes (if it will be needed) a <b>composedTransform</b> object with logical shape (C,P), (C,P,D), or (C,P,D,D) that stores (det J) M^T M, where det J represents <b>cellMeasures</b> and M represents the basis transformations for the reference-space basis.
    */
    PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValues,
             const TensorData<Scalar,DeviceType> cellMeasures,
             const ScalarView<Orientation,DeviceType> orientations);
    
    //! Recomputes point data according to updated basis transformations and cell measures.
    void recomputePointData(const TransformedBasisValues<Scalar,DeviceType> basisValuesWithUpdatedTransformations,
                            const TensorData<Scalar,DeviceType> updatedCellMeasures);
    
    //! Recomputes point data according to updated basis transformations and cell measures.
    void recomputePointData(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeftWithUpdatedTransformations,
                            const TensorData<Scalar,DeviceType> updatedCellMeasures,
                            const TransformedBasisValues<Scalar,DeviceType> basisValuesRightWithUpdatedTransformations);
    
    /** \brief   Allocates storage for a fully-assembled matrix.
        \return <b>integrals</b>, a container with logical shape (C,F1,F2), suitable for passing to assemble().
    */
    Data<Scalar,DeviceType> allocateMatrixStorage() const;
  
    /** \brief   Allocates and returns a container with shape (C,F1).
        \return  a container with logical shape (C,F1), suitable for passing to extractColumn().
    */
    Data<Scalar,DeviceType> allocateColumnStorage() const;
    
    /** \brief   Allocates and returns a container with shape (C,F), where F=min(F1,F2).
        \return  a container with logical shape (C,F), suitable for passing to extractDiagonal().
    */
    Data<Scalar,DeviceType> allocateDiagonalStorage() const;
    
    /** \brief   Allocates and returns a view with shape (C).
        \return  a container with logical shape (C), suitable for passing to extractEntry().
    */
    Data<Scalar,DeviceType> allocateEntryStorage() const;
    
    /** \brief   Allocates and returns a view with shape (C,F2).
        \return  a container with logical shape (C,F2), suitable for passing to extractRow().
    */
    Data<Scalar,DeviceType> allocateRowStorage() const;
    
    /** \brief   Allocates and returns a vector with shape (C,F2).
        \return  a container with logical shape (C,F2), suitable for passing to apply() as input.
    */
    ScalarView<Scalar,DeviceType> allocateInputVector() const;
    
    /** \brief   Allocates and returns a multi-vector with shape (C,F2,N), where N is the number of input vectors.
        \return  a container with logical shape (C,F2,N), suitable for passing to apply() as input.
    */
    ScalarView<Scalar,DeviceType> allocateInputMultiVector(const ordinal_type &n) const;
    
    /** \brief   Allocates and returns a vector with shape (C,F1).
        \return  a container with logical shape (C,F1), suitable for passing to apply() as output.
    */
    ScalarView<Scalar,DeviceType> allocateOutputVector() const;
    
    /** \brief   Allocates and returns a multi-vector with shape (C,F1,N), where N is the number of output vectors.
        \return  a container with logical shape (C,F1,N), suitable for passing to apply() as output.
    */
    ScalarView<Scalar,DeviceType> allocateOutputMultiVector(const ordinal_type &n) const;

    /** \brief   Allocates and returns workspace storage suitable for providing to apply() to an vector with a workset of size worksetSize.
        \return
    */
    Kokkos::View<Scalar*,DeviceType> allocateWorkspace(const ordinal_type &worksetSize) const;
    
    /** \brief   Allocates and returns workspace storage suitable for providing to apply() to an n-multivector with a workset of size worksetSize.
        \return
    */
    Kokkos::View<Scalar*,DeviceType> allocateWorkspace(const ordinal_type &worksetSize, const ordinal_type &n) const;
    
    /** \brief  Applies the matrix to <b>inputVector</b>, placing the result in <b>outputVector</b>, without explicit assembly and storage of the matrix itself.

        \param  outputVector             [out] - the result of applying the matrix to the input vector
        \param  inputVector               [in] - the vector to which the matrix will applied
        \param  workspace                   [in] - memory allocated by allocateWorkspace(), called with arguments (C,N) or, if worksetSize is nonzero, (worksetSize,N).
        \param  worksetSize               [in] - the number of cells to apply to at a time.  If unspecified, applies to all cells.
        
        <b>outputVector</b> and <b>inputVector</b> may have shapes (C,F1) and (C,F2), representing single vectors, or shapes (C,F1,N) and (C,F2,N), representing multi-vectors.
    */
    template <typename OutputViewType, typename InputViewType>
    void apply(const OutputViewType &outputVector,
               const  InputViewType & inputVector,
               const Kokkos::View<Scalar*,DeviceType> &workspace,
               const bool sumInto = false,
               const int worksetSizeIn = 0) const;
    
    /** \brief   Fully assembles the matrix.
        \param   integrals          [out] - Output matrix, with logical shape (C,F,F).  See allocateMatrixStorage().
    */
    void assemble(Data<Scalar,DeviceType> &integrals) const;
    
    /** \brief   Extracts the <b>j</b>th column of the matrix, placing it in <b>column</b>.
        \param   row          [out] - Output container with logical shape (C,F1).  See allocateColumnStorage().
    */
    void extractColumn(const Data<Scalar,DeviceType> &column, const ordinal_type &j) const;
    
    /** \brief   Extracts the diagonal of the matrix, placing it in <b>diagonal</b>.  If the matrix is not square, returns the portion of the matrix for which i==j.
        \param   diagonal          [out] - Output container with logical shape (C,F), F=min(F1,F2).  See allocateDiagonalStorage().
    */
    void extractDiagonal(const Data<Scalar,DeviceType> &diagonal) const;
    
    /** \brief   Extracts the <b>i</b>th row of the matrix, placing it in <b>row</b>.
        \param   row          [out] - Output view with logical shape (C,F2).  See allocateRowStorage().
    */
    void extractRow(const Data<Scalar,DeviceType> &row, const ordinal_type &i) const;
    
    /** \brief   Extracts matrix entries for each cell at (i,j).
        \param   entry [out] - Output container with logical shape (C).  See allocateEntryStorage().
        \param   i          [in] - row index.
        \param   j          [in] - column index.
     
     \note This method is asymptotically more expensive per entry than extracting diagonals, rows, and columns.  The cost of this evaluation scales with the number of quadrature points, generally O(p^d), with no possibility of reuse of intermediate sums from one row/column to another.  The diagonal, row, and column extraction methods, on the other hand, produce O(p^d) values in O(p^{d+1}) time.
    */
    void extractEntry(const Data<Scalar,DeviceType> &entry, const ordinal_type &i, const ordinal_type &j) const;
    
    
    
    //! Returns the cumulative estimated floating-point operation count for calls to gemm() from PAMatrix.
    static double gemmFlopCount();
    
    //! Accumulates into a static variable with cumulative flop count.  Returns its current value.
    static double recordGEMMFlops(const ordinal_type &M, const ordinal_type &N, const ordinal_type &K);
    
    //! Returns the cumulative time spent in calls to gemm() from PAMatrix.
    static double gemmTimeSeconds(); // returns cumulative time spent in gemm() calls
    
    //! if baseFlopCount and baseTimeSeconds are provided, these will be subtracted from the cumulative counts to provide a throughput since a previous measurement.
    static double gemmThroughputGFlops(const double baseFlopCount = 0, const double baseTimeSeconds = 0);

  }; // end PAMatrix class

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_PAMatrixDef.hpp>

#endif
