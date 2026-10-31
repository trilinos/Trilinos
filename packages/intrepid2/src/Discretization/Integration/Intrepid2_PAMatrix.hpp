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
    std::map<PointDataSpec,View1D> _pointDataCache; // copies of appropriate slices of _composedWeightedTransform
    TensorData<Scalar,DeviceType> _cellMeasures; // (C,P); used for separable case
    TransformedBasisValues<Scalar,DeviceType> _basisValuesLeft, _basisValuesRight;
    ScalarView<Orientation,DeviceType> _orientations;
    static constexpr bool layoutLeft_ = true; // BLAS expects this
    int maxIntermediateSize_ = 0;

    using ComponentSequence = std::tuple<std::vector<OpSpec>, PointDataSpec, std::vector<OpSpec>, int, int, int, int>; // left, pointData, right, left offset, left output span, right offset, right input span
    std::vector<ComponentSequence> componentIntegralsToSum_;
    std::vector<std::pair<int,int>> componentIntegralGroups_; // for each entry in componentIntegralsToSum_: (left group, right group), where a group is a (family, component) pair; used by assemble()
    // Ordinal filters (e.g., for Serendipity bases): matrix row/column i corresponds to full tensor-basis field
    // ordinal filter(i).  Empty when there is no filter.  componentIntegralsToSum_ refers to full tensor-basis ordinals.
    // Supported by assemble() (unit-column path); other operations reject filtered PAMatrix objects.
    Kokkos::View<ordinal_type*,DeviceType> leftOrdinalFilter_, rightOrdinalFilter_;
    int numUnfilteredFieldsLeft_ = 0, numUnfilteredFieldsRight_ = 0;
    bool hasOrdinalFilter() const { return (leftOrdinalFilter_.extent_int(0) > 0) || (rightOrdinalFilter_.extent_int(0) > 0); }

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

    /** \brief   Allocates storage for a fully-assembled matrix.
        \return <b>integrals</b>, a container with logical shape (C,F1,F2), suitable for passing to assemble().
    */
    Data<Scalar,DeviceType> allocateMatrixStorage() const;

    /** \brief   Allocates and returns a container with shape (C,F), where F=min(F1,F2).
        \return  a container with logical shape (C,F), suitable for passing to extractDiagonal().
    */
    Data<Scalar,DeviceType> allocateDiagonalStorage() const;

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

    /** \brief   Fully assembles the matrix, using sum factorization.
        \param   integrals          [out] - Output matrix, with logical shape (C,F,F).  See allocateMatrixStorage().

        For supported structures (scalar point weights; one-dimensional tensor components; no orientations), the matrix
        is assembled by direct sum factorization: the point weights are contracted against the per-component products of
        left and right 1D operators, one component at a time, by GEMMs over all cells of a workset.  Passes that share
        their last-component operators share the last contraction, and the output is written once.  With sparse ordinal
        filters (e.g., Serendipity bases), the last contraction computes only the kept entries.  Workset size is chosen
        from a memory budget (environment variable INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB, default 512).  Setting
        INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS selects the previous path, which computes blocks of columns by
        contraction of unit trial columns, for these structures too.  Other structures are assembled by apply() to
        blocked unit columns.  Ordinal filters are supported on the first two paths only.
    */
    void assemble(Data<Scalar,DeviceType> &integrals) const;

    /** \brief   Extracts the diagonal of the matrix, placing it in <b>diagonal</b>.  If the matrix is not square, returns the portion of the matrix for which i==j.
        \param   diagonal          [out] - Output container with logical shape (C,F), F=min(F1,F2).  See allocateDiagonalStorage().
    */
    void extractDiagonal(const Data<Scalar,DeviceType> &diagonal) const;

    //! The left basis values (with transform) from which the matrix was constructed.
    const TransformedBasisValues<Scalar,DeviceType> & leftBasisValues()  const { return _basisValuesLeft;  }

    //! The right basis values (with transform) from which the matrix was constructed.
    const TransformedBasisValues<Scalar,DeviceType> & rightBasisValues() const { return _basisValuesRight; }

    //! The cell measures (quadrature weights times Jacobian determinants) from which the matrix was constructed.
    const TensorData<Scalar,DeviceType> & cellMeasures() const { return _cellMeasures; }

    //! The cell orientations (empty if none were provided).
    const ScalarView<Orientation,DeviceType> & orientations() const { return _orientations; }

    /** \name GEMM statistics (opt-in)
        PAMatrix can record an estimated floating-point operation count and the time spent in the GEMM calls it issues
        (in apply() and related methods), for performance studies.  Recording is off by default; enable it with
        setGEMMStatisticsEnabled(true).

        Limitations:
        - The statistics are process-wide and cumulative: they are shared by all PAMatrix instances (for a given
          DeviceType and Scalar), and are never reset.  To measure an interval, record gemmFlopCount() and
          gemmTimeSeconds() at its start and pass them to gemmThroughputGFlops().
        - Recording is not thread-safe; do not issue PAMatrix operations from multiple host threads while it is enabled.
        - Only the GEMMs PAMatrix itself issues are counted.  Other work (e.g., point-data multiplication and the
          sum-factorized assembly kernels) is neither counted nor timed.
        - Flop counts are estimates (2K-1 operations per entry of an MxK by KxN product).
        - Times come from a Teuchos::TimeMonitor named "gemm", and include any time the GEMM spends waiting on the
          execution space; the timer itself adds per-call overhead.
    */
    //@{
    //! Enables or disables recording of GEMM statistics (off by default).
    static void setGEMMStatisticsEnabled(const bool enabled);

    //! Returns true if GEMM statistics are being recorded.
    static bool gemmStatisticsEnabled();

    //! Returns the cumulative estimated floating-point operation count for recorded GEMM calls.
    static double gemmFlopCount();

    //! Adds the estimated flop count for an MxK by KxN product to the cumulative count, if recording is enabled.  Returns the cumulative count.
    static double recordGEMMFlops(const ordinal_type &M, const ordinal_type &N, const ordinal_type &K);

    //! Returns the cumulative time spent in recorded GEMM calls.
    static double gemmTimeSeconds();

    //! Returns the GEMM throughput in GFlop/s; baseFlopCount and baseTimeSeconds, if provided, are subtracted from the cumulative counts, to give the throughput since a previous measurement.
    static double gemmThroughputGFlops(const double baseFlopCount = 0, const double baseTimeSeconds = 0);
    //@}

  }; // end PAMatrix class

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_PAMatrixDef.hpp>

#endif
