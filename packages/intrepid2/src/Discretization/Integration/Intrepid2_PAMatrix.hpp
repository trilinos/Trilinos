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

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"

#include "Kokkos_Core.hpp"

#include <map>
#include <tuple>
#include <vector>

namespace Intrepid2 {

  /** \brief Options for PAMatrix::assemble().  The defaults are appropriate in most cases; the other settings are
             primarily for testing and performance comparison.
  */
  struct PAMatrixAssemblyOptions
  {
    //! Algorithm used by assemble().
    enum class Algorithm
    {
      Automatic,     //!< direct sum factorization where the structure supports it, column-blocked unit-column contraction otherwise
      UnitColumns,   //!< column-blocked unit-column contraction where supported, apply() to unit columns otherwise
      ApplyToColumns //!< apply() to blocked unit columns (the most general path)
    };
    Algorithm algorithm = Algorithm::Automatic;

    //! Workspace budget, in MB, used to choose the workset size for assembly.
    double workspaceMB = 512.0;

    //! Use of the sparse last stage in the direct path when ordinal filters are present.
    enum class SparseLastStage
    {
      Automatic, //!< used when fewer than a quarter of the (row, column) pairs of the unfiltered matrix are kept
      Always,    //!< used whenever there is an ordinal filter
      Never
    };
    SparseLastStage sparseLastStage = SparseLastStage::Automatic;

    /** \brief Returns default options, overridden by any of the following environment variables (for testing and
               performance comparison without recompiling):
        - INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM: "automatic", "unit-columns", or "apply"
        - INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB: positive number
        - INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST: "automatic", "1" (always), or "0" (never)
        - INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS (deprecated; equivalent to INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM=unit-columns)
        The variables are read once, on the first call; invalid values throw std::invalid_argument.
    */
    static PAMatrixAssemblyOptions defaultsWithEnvironmentOverrides();
  };

  /** \class Intrepid2::PAMatrix
      \brief Partially assembled (matrix-free) representation of element matrices for integrals of tensor-product
             (hypercube) bases, supporting matrix-free application by sum factorization, sum-factorized full assembly,
             and diagonal extraction.

      The matrix is represented by the 1D basis operators (values or derivatives at quadrature points) of each tensor
      component, together with the pointwise data (transforms composed with cell measures) at quadrature points.
  */
  template<typename DeviceType = Kokkos::DefaultExecutionSpace::device_type,
           typename Scalar = double>
  class PAMatrix {
  public:
    using View1D = Kokkos::View<Scalar*,DeviceType>;
    using BasisPtrType = BasisPtr<DeviceType,Scalar,Scalar>;

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
        return std::tie(C, P, a0, b0, aSpan, bSpan) < std::tie(otherSpec.C, otherSpec.P, otherSpec.a0, otherSpec.b0, otherSpec.aSpan, otherSpec.bSpan);
      }

      bool operator==(const PointDataSpec &otherSpec) const {
        return std::tie(C, P, a0, b0, aSpan, bSpan) == std::tie(otherSpec.C, otherSpec.P, otherSpec.a0, otherSpec.b0, otherSpec.aSpan, otherSpec.bSpan);
      }
    };
    //! left, pointData, right, left offset, left output span, right offset, right input span
    using ComponentSequence = std::tuple<std::vector<OpSpec>, PointDataSpec, std::vector<OpSpec>, int, int, int, int>;

    /** \brief   Constructs a <b>PAMatrix</b>  representing the contraction of \a <b>basisValuesLeft</b> against \a <b>basisValuesRight</b> containers on
                 point and space dimensions, weighting each point according to <b>cellMeasures</b>.

        \param  basisValuesLeft      [in] - Left input container, with logical shape (C,F1,P,D)
        \param  cellMeasures         [in] - Point weight container, with logical shape (C,P)
        \param  basisValuesRight     [in] - Right input container with logical shape (C,F2,P,D)
        \param  orientations         [in] - orientations container, with shape (C) (optional)
        \param  basisLeft            [in] - the basis for the left values; required if any orientation is not aligned to the reference
        \param  basisRight           [in] - the basis for the right values; required if any orientation is not aligned to the reference

        With orientations, the matrix represented is the one for the oriented bases, O_L phi and O_R psi, where O_L and O_R
        are the operators applied by OrientationTools::modifyBasisByOrientation(); i.e., O_L M O_R^T, where M is the
        matrix for the unoriented bases.

        On construction, computes (if it will be needed) a <b>composedTransform</b> object with logical shape (C,P), (C,P,D), or (C,P,D,D) that stores (det J) M^T_L M_R, where det J represents <b>cellMeasures</b> and M_L and M_R represent the basis transformations for the left and right basis, respectively.
    */
    PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
             const TensorData<Scalar,DeviceType> cellMeasures,
             const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
             const ScalarView<Orientation,DeviceType> orientations = ScalarView<Orientation,DeviceType>(),
             const BasisPtrType basisLeft  = Teuchos::null,
             const BasisPtrType basisRight = Teuchos::null);

    /** \brief   Constructs a <b>PAMatrix</b>  representing the contraction of \a <b>basisValues</b> against itself in
                 point and space dimensions, weighting each point according to <b>cellMeasures</b>.

        \param  basisValues              [in] - Transformed basis values input container, with logical shape (C,F,P,D)
        \param  cellMeasures             [in] - Point weight container, with logical shape (C,P)
        \param  orientations             [in] - orientations container, with shape (C)
        \param  basis                    [in] - the basis; required if any orientation is not aligned to the reference

        On construction, computes (if it will be needed) a <b>composedTransform</b> object with logical shape (C,P), (C,P,D), or (C,P,D,D) that stores (det J) M^T M, where det J represents <b>cellMeasures</b> and M represents the basis transformations for the reference-space basis.
    */
    PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValues,
             const TensorData<Scalar,DeviceType> cellMeasures,
             const ScalarView<Orientation,DeviceType> orientations,
             const BasisPtrType basis = Teuchos::null);

    //! Constructs an empty PAMatrix, which may be assigned to from a constructed one.
    PAMatrix() = default;

    /** \brief   Allocates storage for a fully-assembled matrix.
        \return <b>integrals</b>, a container with logical shape (C,F1,F2), suitable for passing to assemble().
    */
    Data<Scalar,DeviceType> allocateMatrixStorage() const;

    /** \brief   Allocates and returns a container with shape (C,F), where F=min(F1,F2).
        \return  a container with logical shape (C,F), suitable for passing to extractDiagonal().
    */
    Data<Scalar,DeviceType> allocateDiagonalStorage() const;

    /** \brief   Allocates and returns workspace storage suitable for providing to apply() to a vector with a workset of size worksetSize.
    */
    Kokkos::View<Scalar*,DeviceType> allocateWorkspace(const ordinal_type &worksetSize) const;

    /** \brief   Allocates and returns workspace storage suitable for providing to apply() to an n-multivector with a workset of size worksetSize.
    */
    Kokkos::View<Scalar*,DeviceType> allocateWorkspace(const ordinal_type &worksetSize, const ordinal_type &n) const;

    /** \brief  Applies the matrix to <b>inputVector</b>, placing the result in <b>outputVector</b>, without explicit assembly and storage of the matrix itself.

        \param  outputVector             [out] - the result of applying the matrix to the input vector
        \param  inputVector               [in] - the vector to which the matrix will applied
        \param  workspace                 [in] - memory allocated by allocateWorkspace(), called with arguments (C,N) or, if worksetSize is nonzero, (worksetSize,N).
        \param  sumInto                   [in] - if true, the result is added to <b>outputVector</b>.
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
        \param   integrals          [out] - Output matrix, with logical shape (C,F1,F2).  See allocateMatrixStorage().
        \param   options             [in] - assembly options; see PAMatrixAssemblyOptions.

        For supported structures (scalar point weights and one-dimensional tensor components), the matrix is assembled by
        direct sum factorization: the point weights are contracted against the per-component products of left and right 1D
        operators, one component at a time, by GEMMs over all cells of a workset.  Passes that share their last-component
        operators share the last contraction, and the output is written once.  With sparse ordinal filters (e.g., Serendipity
        bases), the last contraction computes only the kept entries.  Other structures are assembled by contraction of
        blocked unit columns, or by apply() to them.  Workset sizes are chosen from the workspace budget in the options.
        Orientations, if any, are applied to the assembled matrix (as by OrientationTools::modifyMatrixByOrientation()).
    */
    void assemble(Data<Scalar,DeviceType> &integrals,
                  const PAMatrixAssemblyOptions &options = PAMatrixAssemblyOptions::defaultsWithEnvironmentOverrides()) const;

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

    //! The left basis (null if none was provided).
    const BasisPtrType & leftBasis()  const { return _basisLeft;  }

    //! The right basis (null if none was provided).
    const BasisPtrType & rightBasis() const { return _basisRight; }

    //! True if the left or right basis values have an ordinal filter (e.g., for Serendipity bases).
    bool hasOrdinalFilter() const { return (leftOrdinalFilter_.extent_int(0) > 0) || (rightOrdinalFilter_.extent_int(0) > 0); }

    //! True if orientations were provided and at least one cell is not aligned to the reference.
    bool hasNontrivialOrientations() const { return _hasNontrivialOrientations; }

    /** \name GEMM statistics (opt-in)
        PAMatrix can record an estimated floating-point operation count and the time spent in the GEMM calls it issues,
        for performance studies.  Recording is off by default; enable it with setGEMMStatisticsEnabled(true).

        Limitations:
        - The statistics are process-wide and cumulative: they are shared by all PAMatrix instances, and are never reset.
          To measure an interval, record gemmFlopCount() and gemmTimeSeconds() at its start and pass them to
          gemmThroughputGFlops().
        - Recording is not thread-safe; do not issue PAMatrix operations from multiple host threads while it is enabled.
        - Only the GEMMs PAMatrix itself issues are counted.  Other work (e.g., point-data multiplication and the
          sum-factorized assembly kernels) is neither counted nor timed.
        - Flop counts are estimates (2K-1 operations per entry of an MxK by KxN product).
        - Times come from a Teuchos::TimeMonitor named "Intrepid2::PAMatrix gemm"; each recorded GEMM is followed by a
          fence of its execution space, so that the times include the GEMM's execution.
    */
    //@{
    //! Enables or disables recording of GEMM statistics (off by default).
    static void setGEMMStatisticsEnabled(const bool enabled);

    //! Returns true if GEMM statistics are being recorded.
    static bool gemmStatisticsEnabled();

    //! Returns the cumulative estimated floating-point operation count for recorded GEMM calls.
    static double gemmFlopCount();

    //! Returns the cumulative time spent in recorded GEMM calls.
    static double gemmTimeSeconds();

    //! Returns the GEMM throughput in GFlop/s; baseFlopCount and baseTimeSeconds, if provided, are subtracted from the cumulative counts, to give the throughput since a previous measurement.
    static double gemmThroughputGFlops(const double baseFlopCount = 0, const double baseTimeSeconds = 0);
    //@}

    /** \brief True if assemble() can use direct sum factorization (or contraction of unit columns) for this matrix's
               structure: scalar point weights, one-dimensional tensor components, and equal left and right tensor component
               counts.  Otherwise, assemble() applies the matrix to blocked unit columns.
    */
    bool supportsDirectAssembly() const;

    /** \name Implementation details
        These are not part of the supported interface.  They are public only because they define KOKKOS_LAMBDAs, and
        CUDA does not allow extended __host__ __device__ lambdas in private or protected member functions.
    */
    //@{
    void init(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
              const TensorData<Scalar,DeviceType> cellMeasures,
              const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
              const ScalarView<Orientation,DeviceType> orientations);

    //! assembles the matrix for the unoriented bases
    void assembleUnoriented(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const;

    //! assembles by direct sum factorization (requires supportsDirectAssembly())
    void assembleDirect(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const;

    //! assembles by contraction of unit columns (requires supportsDirectAssembly())
    void assembleUnitColumns(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const;

    //! assembles by apply() to blocked unit columns
    void assembleByApply(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const;

    //! applies orientations in place to an assembled (C,F1,F2) matrix
    void orientAssembledMatrix(Data<Scalar,DeviceType> &integrals) const;

    //! the workspace size, in entries per cell and vector, required by apply()
    size_t applyWorkspacePerCell() const;

    /** \brief Applies the unoriented matrix to cells [startCell, startCell+Cw), accumulating into <b>output</b>.  Cell c of the
               workset is cell c+inputCellOffset of <b>input</b> (shape (.,F2[,N])) and cell c+outputCellOffset of
               <b>output</b> (shape (.,F1[,N])).  workspace1 and workspace2 must each hold maxIntermediateSize_ * Cw * N entries.
    */
    template <typename OutputViewType, typename InputViewType>
    void applyUnorientedWorkset(const OutputViewType &output, const ordinal_type outputCellOffset,
                                const InputViewType &input, const ordinal_type inputCellOffset,
                                const ordinal_type startCell, const ordinal_type Cw, const ordinal_type N,
                                Scalar *workspace1, Scalar *workspace2) const;
    //@}

  private:
    Data<Scalar,DeviceType> _composedWeightedTransform; // (C,P[,D1[,D2]]), used for general case
    std::map<PointDataSpec,View1D> _pointDataCache; // copies of appropriate slices of _composedWeightedTransform
    TensorData<Scalar,DeviceType> _cellMeasures; // (C,P)
    TransformedBasisValues<Scalar,DeviceType> _basisValuesLeft, _basisValuesRight;
    ScalarView<Orientation,DeviceType> _orientations;
    BasisPtrType _basisLeft, _basisRight;
    bool _hasNontrivialOrientations = false;
    size_t maxIntermediateSize_ = 0;

    std::vector<ComponentSequence> componentIntegralsToSum_;
    std::vector<std::pair<int,int>> componentIntegralGroups_; // for each entry in componentIntegralsToSum_: (left group, right group), where a group is a (family, component) pair; used by assemble()
    // Ordinal filters (e.g., for Serendipity bases): matrix row/column i corresponds to full tensor-basis field
    // ordinal filter(i).  Empty when there is no filter.  componentIntegralsToSum_ refers to full tensor-basis ordinals.
    // Supported by assemble(); other operations reject filtered PAMatrix objects.
    Kokkos::View<ordinal_type*,DeviceType> leftOrdinalFilter_, rightOrdinalFilter_;
    int numUnfilteredFieldsLeft_ = 0, numUnfilteredFieldsRight_ = 0;
  }; // end PAMatrix class

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_PAMatrixDef.hpp>

#endif
