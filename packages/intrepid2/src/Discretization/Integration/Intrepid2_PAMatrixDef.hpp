// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PAMatrixDef.hpp
    \brief  Header file for the Intrepid2::PAMatrix implementations; provides support for matrix partial assembly.
    \author Created by Nathan V. Roberts.
*/

#ifndef __INTREPID2_PAMATRIX_DEF_HPP__
#define __INTREPID2_PAMATRIX_DEF_HPP__

#include "Intrepid2_PAMatrix.hpp"

#include "Intrepid2_DataDimensionInfo.hpp"
#include "Intrepid2_DataTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include <Teuchos_BLAS.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_INTREPID2_KOKKOSKERNELS
#include <KokkosBlas3_gemm.hpp>
#endif

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>

namespace Intrepid2 {

namespace Impl
{
  //! process-wide GEMM statistics state for PAMatrix; see PAMatrix::setGEMMStatisticsEnabled()
  struct PAMatrixGEMMStatistics
  {
    static bool &enabled()           { static bool value = false; return value; }
    static double &cumulativeFlops() { static double value = 0;   return value; }
    static Teuchos::RCP<Teuchos::Time> timer()
    {
      static const auto value = Teuchos::TimeMonitor::getNewTimer("Intrepid2::PAMatrix gemm");
      return value;
    }
  };

  //! True if a DeviceType's execution space runs on the host and its memory is host-accessible, so that GEMMs can be
  //! issued directly to the host BLAS (Teuchos::BLAS).
  template<typename DeviceType>
  constexpr bool isHostDevice()
  {
    using ExecutionSpace = typename DeviceType::execution_space;
    using MemorySpace    = typename DeviceType::memory_space;
    return Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemorySpace>::accessible
        && Kokkos::SpaceAccessibility<ExecutionSpace, Kokkos::HostSpace>::accessible;
  }

  //! True if PAMatrix<DeviceType,...> issues its GEMMs through Teuchos::BLAS: host devices without KokkosKernels, and
  //! serial host devices with it.  (With KokkosKernels, threaded host execution spaces such as OpenMP use
  //! KokkosBlas::gemm, which calls the TPL BLAS where one is enabled.)
  template<typename DeviceType>
  constexpr bool gemmUsesTeuchosBLAS()
  {
#ifdef HAVE_INTREPID2_KOKKOSKERNELS
  #ifdef KOKKOS_ENABLE_SERIAL
    return isHostDevice<DeviceType>() && std::is_same<typename DeviceType::execution_space, Kokkos::Serial>::value;
  #else
    return false;
  #endif
#else
    return isHostDevice<DeviceType>();
#endif
  }

  //! True if PAMatrix<DeviceType,...> supports GEMMs, i.e., its operations can be used in this build.  Non-host devices
  //! (e.g., CUDA, HIP, SYCL) require Intrepid2 to be built with KokkosKernels.
  template<typename DeviceType>
  constexpr bool gemmSupported()
  {
#ifdef HAVE_INTREPID2_KOKKOSKERNELS
    return true;
#else
    return isHostDevice<DeviceType>();
#endif
  }

  /** \brief Column-major GEMM, C = alpha * op(A) * op(B) + beta * C, with op(A) of size M x K and op(B) of size K x N, on
             data accessible from DeviceType.  LDA is the leading dimension of A; B and C are packed (LDB = K or N, LDC = M).
             The backend is chosen from DeviceType: Teuchos::BLAS for host devices (without KokkosKernels, or serial with
             it), and KokkosBlas::gemm on DeviceType's execution space otherwise.  Returns once the result is available to
             subsequent kernels on DeviceType's execution space, and to the host if DeviceType is a host device.
  */
  template<typename DeviceType,typename Scalar>
  void gemm(const char transA, const char transB,
            const ordinal_type M, const ordinal_type N, const ordinal_type K,
            const Scalar alpha, const Scalar* A, const ordinal_type LDA,
            const Scalar* B, const Scalar beta, Scalar* C)
  {
    static_assert(gemmSupported<DeviceType>(), "PAMatrix on this device type requires Intrepid2 to be built with KokkosKernels");
    INTREPID2_TEST_FOR_EXCEPTION(LDA <= 0, std::invalid_argument, "LDA must be positive");
    if ((M == 0) || (N == 0)) return;

    using ExecutionSpace = typename DeviceType::execution_space;
    const bool recordStatistics = PAMatrixGEMMStatistics::enabled();
    Teuchos::RCP<Teuchos::TimeMonitor> gemmTimer;
    if (recordStatistics)
    {
      ExecutionSpace().fence(); // so that the timer measures only this GEMM
      gemmTimer = Teuchos::rcp(new Teuchos::TimeMonitor(*PAMatrixGEMMStatistics::timer()));
    }

    if constexpr (gemmUsesTeuchosBLAS<DeviceType>())
    {
      // any outstanding kernels that produce A or B must complete before the host BLAS reads them
      ExecutionSpace().fence();
      const Teuchos::ETransp trA = (transA == 'T') ? Teuchos::TRANS : (transA == 'C') ? Teuchos::CONJ_TRANS : Teuchos::NO_TRANS;
      const Teuchos::ETransp trB = (transB == 'T') ? Teuchos::TRANS : (transB == 'C') ? Teuchos::CONJ_TRANS : Teuchos::NO_TRANS;
      const ordinal_type LDB = std::max<ordinal_type>(1, (trB == Teuchos::NO_TRANS) ? K : N);
      const ordinal_type LDC = std::max<ordinal_type>(1, M);
      Teuchos::BLAS<int,Scalar> blas;
      blas.GEMM(trA, trB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
    else
    {
#ifdef HAVE_INTREPID2_KOKKOSKERNELS
      using ConstView2D = Kokkos::View<const Scalar**, Kokkos::LayoutLeft, DeviceType, Kokkos::MemoryUnmanaged>;
      using      View2D = Kokkos::View<      Scalar**, Kokkos::LayoutLeft, DeviceType, Kokkos::MemoryUnmanaged>;
      using Stride2D    = Kokkos::View<const Scalar**, Kokkos::LayoutStride, DeviceType, Kokkos::MemoryUnmanaged>;
      const ordinal_type rowsA = (transA == 'N') ? M : K;
      const ordinal_type colsA = (transA == 'N') ? K : M;
      Kokkos::LayoutStride layoutA(rowsA, 1, colsA, LDA);
      Stride2D AView(A, layoutA);
      ConstView2D BView = (transB != 'N') ? ConstView2D(B,N,K) : ConstView2D(B,K,N);
      View2D CView(C,M,N);
      const char tA[] = {transA, '\0'};
      const char tB[] = {transB, '\0'};
      KokkosBlas::gemm(ExecutionSpace(), tA, tB, alpha, AView, BView, beta, CView);
#endif
    }

    if (recordStatistics)
    {
      ExecutionSpace().fence(); // so that the timer measures the GEMM's execution, not just its launch
      gemmTimer = Teuchos::null;
      PAMatrixGEMMStatistics::cumulativeFlops() += double(M) * double(N) * double(2 * K - 1);
    }
  }

  //! Parses a positive number from an environment variable; throws if the value is not a positive number.
  inline double parsePositiveEnvironmentNumber(const char *name, const char *value)
  {
    char *end = nullptr;
    const double result = std::strtod(value, &end);
    INTREPID2_TEST_FOR_EXCEPTION((end == value) || (*end != '\0') || !(result > 0) || !(result < std::numeric_limits<double>::infinity()),
                                 std::invalid_argument, std::string(name) + " must be a positive number; got '" + value + "'");
    return result;
  }

  /** \brief For point data A(C,P,Da,Db) and vector-valued B(P[,Db],N,Cw), computes the pointwise matrix-vector product
             C(N,Cw,P[,Da]) for cells startCell..startCell+Cw-1 of A.  A's cell extent is numDataCells, which is either
             the total number of cells or 1 (for data that is constant in the cell).
  */
  template<typename DeviceType, class Scalar>
  void pointDataMultiply(const ordinal_type numDataCells, const ordinal_type startCell, const ordinal_type Cw,
                         const ordinal_type numPoints, const ordinal_type aSpan, const ordinal_type bSpan,
                         const Scalar* A, const Scalar *B, Scalar *C, const ordinal_type N)
  {
    // Limitation: only aSpan == 1 is supported.  Non-trivial Da would require packing the vector index into the point component
    // over which the left integral performs a dot product (and additional arguments describing the point components).
    INTREPID2_TEST_FOR_EXCEPTION(aSpan != 1, std::invalid_argument, "aSpan != 1 is not yet supported");

    using ExecutionSpace = typename DeviceType::execution_space;
    // Iterate in the output's memory order (n fastest, then cell, then point, then a), so that writes to C are contiguous;
    // MDRangePolicy's default iteration makes the last index fastest on host.
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{aSpan,numPoints,Cw,N});

    Kokkos::parallel_for("pointDataMultiply", policy,
                         KOKKOS_LAMBDA(const ordinal_type &a, const ordinal_type &point, const ordinal_type &cell, const ordinal_type &n)
    {
      // layout left for A, B, C
      const size_t dataCell = (numDataCells == 1) ? 0 : size_t(cell + startCell);
      const size_t strideA  = size_t(aSpan) * numPoints * numDataCells; // stride in b
      const Scalar* entryA = A + dataCell + (point + size_t(a) * numPoints) * numDataCells; // (C,P,Da,Db)
      const Scalar* entryB = B + point + (size_t(n) + size_t(cell) * N) * bSpan * numPoints; // (P[,Db],N,Cw)
      Scalar value = 0;
      for (int b=0; b<bSpan; b++)
      {
        value += *entryA * *entryB;
        entryA += strideA;
        entryB += numPoints;
      }
      C[n + (size_t(cell) + (point + size_t(a) * numPoints) * Cw) * N] = value; // (N,Cw,P[,Da])
    });
  }

  //! Given (C,P[,D,D]) transform and (C,P) pointwise weights, construct a suitable container for storing the pointwise weighted transform.
  template<typename DeviceType,class Scalar>
  Data<Scalar,DeviceType> allocateComposedWeightedTransform(const Data<Scalar,DeviceType> &composedTransform,
                                                            const TensorData<Scalar,DeviceType> &pointWeights)
  {
    auto cellDimInfo = composedTransform.getDimensionInfo(0); // cell dimension
    int numTensorComponents = pointWeights.numTensorComponents();
    if (pointWeights.separateFirstComponent())
    {
      auto cellDimInfo_points = pointWeights.getTensorComponent(0).getDimensionInfo(0);
      cellDimInfo = combinedDimensionInfo(cellDimInfo, cellDimInfo_points);
    }
    else
    {
      for (int r=0; r<numTensorComponents; r++)
      {
        auto cellDimInfo_r = pointWeights.getTensorComponent(r).getDimensionInfo(0);
        cellDimInfo = combinedDimensionInfo(cellDimInfo, cellDimInfo_r);
      }
    }

    int numPoints = composedTransform.extent_int(1);
    DimensionInfo pointDimInfo {numPoints,GENERAL,numPoints,numPoints,-1};

    if (composedTransform.rank() == 2)
    {
      return Data<Scalar,DeviceType>({cellDimInfo,pointDimInfo});
    }
    else if (composedTransform.rank() == 3)
    {
      auto D1DimInfo = composedTransform.getDimensionInfo(2);
      return Data<Scalar,DeviceType>({cellDimInfo,pointDimInfo,D1DimInfo});
    }
    else if (composedTransform.rank() == 4)
    {
      auto D1DimInfo = composedTransform.getDimensionInfo(2);
      auto D2DimInfo = composedTransform.getDimensionInfo(3);
      return Data<Scalar,DeviceType>({cellDimInfo,pointDimInfo,D1DimInfo,D2DimInfo});
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported rank for composedTransform");
    }
  }

  //! Packs a (F,P) or (F,P,D) 1D operator into a 1D view; if transposeToPoints, the layout is (P[,D],F) (layout left), else (F,P[,D]).
  template<typename DeviceType, class Scalar, class OpDataType>
  Kokkos::View<Scalar*,DeviceType> pack1DOperator(const OpDataType &opData, const bool transposeToPoints)
  {
    using ExecutionSpace = typename DeviceType::execution_space;
    const int opFields = opData.extent_int(0);
    const int opPoints = opData.extent_int(1);
    Kokkos::View<Scalar*,DeviceType> opView(transposeToPoints ? "rightOp 1D view" : "leftOp 1D view", opData.size());
    if (opData.rank() == 2) // (F,P)
    {
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{opFields,opPoints});
      Kokkos::parallel_for("pack 1D opView", policy,
      KOKKOS_LAMBDA(const int &field, const int &pt)
      {
        const int idx = transposeToPoints ? pt + field * opPoints : field + pt * opFields;
        opView(idx) = opData(field,pt);
      });
    }
    else if (opData.rank() == 3) // (F,P,D)
    {
      const int opDim = opData.extent_int(2);
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{opFields,opPoints,opDim});
      Kokkos::parallel_for("pack 1D opView", policy,
      KOKKOS_LAMBDA(const int &field, const int &pt, const int &d)
      {
        const int idx = transposeToPoints ? d + (pt + field * opPoints) * opDim : field + (pt + d * opPoints) * opFields;
        opView(idx) = opData(field,pt,d);
      });
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "PAMatrix: Unsupported component operator rank");
    }
    return opView;
  }
} // namespace Impl

inline PAMatrixAssemblyOptions PAMatrixAssemblyOptions::defaultsWithEnvironmentOverrides()
{
  static const PAMatrixAssemblyOptions options = []()
  {
    PAMatrixAssemblyOptions result;
    if (const char *value = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS"))
    {
      if (std::string(value) != "0") result.algorithm = Algorithm::UnitColumns;
    }
    if (const char *value = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM"))
    {
      const std::string algorithm(value);
      if      (algorithm == "automatic")    result.algorithm = Algorithm::Automatic;
      else if (algorithm == "unit-columns") result.algorithm = Algorithm::UnitColumns;
      else if (algorithm == "apply")        result.algorithm = Algorithm::ApplyToColumns;
      else INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                        "INTREPID2_PAMATRIX_ASSEMBLE_ALGORITHM must be automatic, unit-columns, or apply; got '" + algorithm + "'");
    }
    if (const char *value = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB"))
    {
      result.workspaceMB = Impl::parsePositiveEnvironmentNumber("INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB", value);
    }
    if (const char *value = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST"))
    {
      const std::string setting(value);
      if      (setting == "automatic") result.sparseLastStage = SparseLastStage::Automatic;
      else if (setting == "1")         result.sparseLastStage = SparseLastStage::Always;
      else if (setting == "0")         result.sparseLastStage = SparseLastStage::Never;
      else INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                        "INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST must be automatic, 1, or 0; got '" + setting + "'");
    }
    return result;
  }();
  return options;
}

//! the number of cells per workset that fit perCellEntries entries of Scalar per cell within budgetMB, in [1, numCells]
template<class Scalar>
ordinal_type pamatrixWorksetSize(const double budgetMB, const size_t perCellEntries, const ordinal_type numCells)
{
  INTREPID2_TEST_FOR_EXCEPTION(!(budgetMB > 0), std::invalid_argument, "PAMatrix workspace budget must be positive");
  const double budgetBytes = budgetMB * 1024.0 * 1024.0;
  const double cellsInBudget = budgetBytes / (double(sizeof(Scalar)) * double(std::max<size_t>(perCellEntries, 1)));
  const double cells = std::min(double(numCells), std::floor(cellsInBudget));
  return std::max<ordinal_type>(1, ordinal_type(cells));
}

template<typename DeviceType,class Scalar>
PAMatrix<DeviceType,Scalar>::PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
                                      const TensorData<Scalar,DeviceType> cellMeasures,
                                      const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
                                      const ScalarView<Orientation,DeviceType> orientations,
                                      const BasisPtrType basisLeft,
                                      const BasisPtrType basisRight)
:
_cellMeasures(cellMeasures),
_basisValuesLeft(basisValuesLeft),
_basisValuesRight(basisValuesRight),
_orientations(orientations),
_basisLeft(basisLeft),
_basisRight(basisRight)
{
  init(basisValuesLeft, cellMeasures, basisValuesRight, orientations);
} // PAMatrix()

template<typename DeviceType,class Scalar>
PAMatrix<DeviceType,Scalar>::PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValues,
                                      const TensorData<Scalar,DeviceType> cellMeasures,
                                      const ScalarView<Orientation,DeviceType> orientations,
                                      const BasisPtrType basis)
:
PAMatrix<DeviceType,Scalar>(basisValues,cellMeasures,basisValues,orientations,basis,basis)
{}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::init(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeftInput,
                                       const TensorData<Scalar,DeviceType> cellMeasures,
                                       const TransformedBasisValues<Scalar,DeviceType> basisValuesRightInput,
                                       const ScalarView<Orientation,DeviceType> orientations) {
  using ExecutionSpace = typename DeviceType::execution_space;

  // orientations matter only if some cell is not aligned to the reference
  _hasNontrivialOrientations = false;
  if (orientations.size() > 0)
  {
    INTREPID2_TEST_FOR_EXCEPTION(orientations.extent_int(0) != basisValuesLeftInput.numCells(), std::invalid_argument,
                                 "PAMatrix: orientations must have one entry per cell");
    int numNotAligned = 0;
    Kokkos::parallel_reduce("PAMatrix: count non-aligned orientations", Kokkos::RangePolicy<ExecutionSpace>(0, orientations.extent_int(0)),
    KOKKOS_LAMBDA(const int &cell, int &count) {
      if (!orientations(cell).isAlignedToReference()) count++;
    }, numNotAligned);
    _hasNontrivialOrientations = (numNotAligned > 0);
  }
  INTREPID2_TEST_FOR_EXCEPTION(_hasNontrivialOrientations && ((_basisLeft == Teuchos::null) || (_basisRight == Teuchos::null)),
                               std::invalid_argument, "PAMatrix: the left and right bases must be provided when orientations are not aligned to the reference");

  // Ordinal filters (used for Serendipity bases) select a subset of the fields of an underlying tensor basis.  The
  // component integrals below are set up for the unfiltered (full tensor) bases; the filters are applied when results are
  // written (see assemble()).
  leftOrdinalFilter_  = basisValuesLeftInput .basisValues().ordinalFilter();
  rightOrdinalFilter_ = basisValuesRightInput.basisValues().ordinalFilter();
  auto unfiltered = [](const TransformedBasisValues<Scalar,DeviceType> &values) {
    if (values.basisValues().ordinalFilter().extent_int(0) == 0) return values;
    BasisValues<Scalar,DeviceType> unfilteredValues = values.basisValues();
    unfilteredValues.setOrdinalFilter(Kokkos::View<ordinal_type*,DeviceType>());
    return values.transform().isValid() ? TransformedBasisValues<Scalar,DeviceType>(values.transform(), unfilteredValues)
                                        : TransformedBasisValues<Scalar,DeviceType>(values.numCells(), unfilteredValues);
  };
  const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft  = unfiltered(basisValuesLeftInput);
  const TransformedBasisValues<Scalar,DeviceType> basisValuesRight = unfiltered(basisValuesRightInput);
  numUnfilteredFieldsLeft_  = basisValuesLeft .numFields();
  numUnfilteredFieldsRight_ = basisValuesRight.numFields();
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter() && _hasNontrivialOrientations, std::invalid_argument,
                               "PAMatrix does not yet support ordinal filters together with orientations that are not aligned to the reference");

  const int spaceDim = basisValuesLeft.spaceDim();

  // MARK: checks for supported construction
  INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.spaceDim() != basisValuesRight.spaceDim(), std::invalid_argument, "basisValuesLeft and basisValuesRight must agree on the space dimension");
  INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.numCells() != basisValuesRight.numCells(), std::invalid_argument, "basisValuesLeft and basisValuesRight must agree on the number of cells");
  INTREPID2_TEST_FOR_EXCEPTION((cellMeasures.rank() != 2) || (cellMeasures.extent_int(0) != basisValuesLeft.numCells()) || (cellMeasures.extent_int(1) != basisValuesLeft.numPoints()),
                               std::invalid_argument, "cellMeasures must have shape (C,P), matching the basis values");
  if (_hasNontrivialOrientations)
  {
    INTREPID2_TEST_FOR_EXCEPTION(_basisLeft ->getCardinality() != basisValuesLeft .numFields(), std::invalid_argument, "PAMatrix: basisLeft cardinality does not match basisValuesLeft");
    INTREPID2_TEST_FOR_EXCEPTION(_basisRight->getCardinality() != basisValuesRight.numFields(), std::invalid_argument, "PAMatrix: basisRight cardinality does not match basisValuesRight");
  }

  const int leftFamilyCount  =  basisValuesLeft.basisValues().numFamilies();
  const int rightFamilyCount = basisValuesRight.basisValues().numFamilies();

  // we require that the number of tensor components in the vectors is the same for each vector entry
  // this is not strictly necessary, but it makes implementation easier, and we don't at present anticipate other use cases
  int numTensorComponentsLeft = -1;
  const bool leftIsVectorValued = basisValuesLeft.vectorData().isValid();

  if (leftIsVectorValued)
  {
    const auto &refVectorLeft   = basisValuesLeft.vectorData();
    int numFamiliesLeft         = refVectorLeft.numFamilies();
    int numVectorComponentsLeft = refVectorLeft.numComponents();
    for (int familyOrdinal=0; familyOrdinal<numFamiliesLeft; familyOrdinal++)
    {
      for (int vectorComponent=0; vectorComponent<numVectorComponentsLeft; vectorComponent++)
      {
        const TensorData<Scalar,DeviceType> &tensorData = refVectorLeft.getComponent(familyOrdinal,vectorComponent);
        if (tensorData.numTensorComponents() > 0)
        {
          if (numTensorComponentsLeft == -1)
          {
            numTensorComponentsLeft = tensorData.numTensorComponents();
          }
          INTREPID2_TEST_FOR_EXCEPTION(numTensorComponentsLeft != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in basisValuesLeft must have the same number of tensor components as every other");
        }
      }
    }
  }
  else
  {
    numTensorComponentsLeft = basisValuesLeft.basisValues().tensorData(0).numTensorComponents(); // family ordinal 0
    for (int familyOrdinal = 0; familyOrdinal < leftFamilyCount; familyOrdinal++)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.basisValues().tensorData(familyOrdinal).numTensorComponents() != numTensorComponentsLeft, std::invalid_argument, "All families must match in the number of tensor components");
    }
  }
  int numTensorComponentsRight = -1;
  const bool rightIsVectorValued = basisValuesRight.vectorData().isValid();

  if (rightIsVectorValued)
  {
    const auto &refVectorRight   = basisValuesRight.vectorData();
    int numFamiliesRight         = refVectorRight.numFamilies();
    int numVectorComponentsRight = refVectorRight.numComponents();
    for (int familyOrdinal=0; familyOrdinal<numFamiliesRight; familyOrdinal++)
    {
      for (int vectorComponent=0; vectorComponent<numVectorComponentsRight; vectorComponent++)
      {
        const auto &tensorData = refVectorRight.getComponent(familyOrdinal,vectorComponent);
        if (tensorData.numTensorComponents() > 0)
        {
          if (numTensorComponentsRight == -1)
          {
            numTensorComponentsRight = tensorData.numTensorComponents();
          }
          INTREPID2_TEST_FOR_EXCEPTION(numTensorComponentsRight != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in basisValuesRight must have the same number of tensor components as every other");
        }
      }
    }
    INTREPID2_TEST_FOR_EXCEPTION(numTensorComponentsRight != numTensorComponentsLeft, std::invalid_argument, "Right families must match left in the number of tensor components");
  }
  else
  {
    // check that right tensor component count agrees with left
    for (int familyOrdinal=0; familyOrdinal< rightFamilyCount; familyOrdinal++)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basisValuesRight.basisValues().tensorData(familyOrdinal).numTensorComponents() != numTensorComponentsLeft, std::invalid_argument, "Right families must match left in the number of tensor components");
    }
  }

  {
    // MARK: prepare composed transformation matrices
    const Data<Scalar,DeviceType> & leftTransform  = basisValuesLeft.transform();
    const Data<Scalar,DeviceType> & rightTransform = basisValuesRight.transform();
    const bool transposeLeft  = true;
    const bool transposeRight = false;
    // transforms can be matrices -- (C,P,D,D): rank 4 -- or scalar weights -- (C,P): rank 2 -- or vector weights -- (C,P,D): rank 3
    Data<Scalar,DeviceType> composedTransform;
    // invalid/empty transforms are used when the identity is intended.
    const int leftRank  = leftTransform.rank();
    const int rightRank = rightTransform.rank();

    if (leftTransform.isValid() && rightTransform.isValid())
    {
      const bool bothRank4 = (leftRank == 4) && (rightRank == 4);
      const bool bothRank3 = (leftRank == 3) && (rightRank == 3);
      const bool bothRank2 = (leftRank == 2) && (rightRank == 2);
      const bool ranks32   = ((leftRank == 3) && (rightRank == 2)) || ((leftRank == 2) && (rightRank == 3));
      const bool ranks42   = ((leftRank == 4) && (rightRank == 2)) || ((leftRank == 2) && (rightRank == 4));

      if (bothRank4) // (C,P,D,D)
      {
        composedTransform = Data<Scalar,DeviceType>::allocateMatMatResult(transposeLeft, leftTransform, transposeRight, rightTransform);
        composedTransform.storeMatMat(transposeLeft, leftTransform, transposeRight, rightTransform);
      }
      else if (bothRank3) // (C,P,D)
      {
        // re-cast leftTransform as a rank 4 (C,P,1,D) object -- a 1 x D matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = leftTransform.getExtents();
        auto variationTypes = leftTransform.getVariationTypes();
        extents[3]               = extents[2];
        extents[2]               = 1;
        variationTypes[3]        = variationTypes[2];
        variationTypes[2]        = CONSTANT;
        auto leftTransformMatrix = leftTransform.shallowCopy(newRank, extents, variationTypes);

        // re-cast rightTransform as a rank 4 (C,P,1,D) object -- a 1 x D matrix at each (C,P)
        extents                  = rightTransform.getExtents();
        variationTypes           = rightTransform.getVariationTypes();
        extents[3]               = extents[2];
        extents[2]               = 1;
        variationTypes[3]        = variationTypes[2];
        variationTypes[2]        = CONSTANT;
        auto rightTransformMatrix = rightTransform.shallowCopy(newRank, extents, variationTypes);

        composedTransform = Data<Scalar,DeviceType>::allocateMatMatResult(transposeLeft, leftTransformMatrix, transposeRight, rightTransformMatrix);
        composedTransform.storeMatMat(transposeLeft, leftTransformMatrix, transposeRight, rightTransformMatrix);
      }
      else if (bothRank2)
      {
        composedTransform = leftTransform.allocateInPlaceCombinationResult(leftTransform, rightTransform);
        composedTransform.storeInPlaceProduct(leftTransform, rightTransform);

        // re-cast composedTranform as a rank 4 (C,P,1,1) object -- a 1 x 1 matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = composedTransform.getExtents();
        auto variationTypes = composedTransform.getVariationTypes();
        composedTransform = composedTransform.shallowCopy(newRank, extents, variationTypes);
      }
      else if (ranks32) // rank 2 / rank 3 combination.
      {
        const auto & rank3Transform = (leftRank == 3) ? leftTransform : rightTransform;
        const auto & rank2Transform = (leftRank == 2) ? leftTransform : rightTransform;

        composedTransform = DataTools::multiplyByCPWeights(rank3Transform, rank2Transform);

        // re-cast composedTransform as a rank 4 object:
        // logically, the original rank-3 transform can be understood as a 1xD matrix.  The composed transform is leftTransform^T * rightTransform, so:
        // - if left  has the rank-3 transform, composedTransform should be a (C,P,D,1) object -- a D x 1 matrix at each (C,P).
        // - if right has the rank-3 transform, composedTransform should be a (C,P,1,D) object -- a 1 x D matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = composedTransform.getExtents();
        auto variationTypes = composedTransform.getVariationTypes();
        if (leftRank != 3)
        {
          // (if left has the rank-3 transform, extents[3] and variationTypes[3] will already be 1 and CONSTANT)
          extents[3]               = extents[2];
          extents[2]               = 1;
          variationTypes[3]        = variationTypes[2];
          variationTypes[2]        = CONSTANT;
        }
        composedTransform = composedTransform.shallowCopy(newRank, extents, variationTypes);
      }
      else if (ranks42) // rank 4 / rank 2 combination.
      {
        if (leftRank == 4)
        {
          // want to transpose left matrix, and multiply by the values from rightTransform
          // start with the multiplication:
          auto composedTransformTransposed = DataTools::multiplyByCPWeights(leftTransform, rightTransform);
          composedTransform = DataTools::transposeMatrix(composedTransformTransposed);
        }
        else // (leftRank == 2)
        {
          composedTransform = DataTools::multiplyByCPWeights(rightTransform, leftTransform);
        }
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else if (leftTransform.isValid())
    {
      // rightTransform is the identity
      switch (leftRank)
      {
        case 4: composedTransform = DataTools::transposeMatrix(leftTransform); break;
        case 3:
        {
          // - if left  has the rank-3 transform, composedTransform should be a (C,P,D,1) object -- a D x 1 matrix at each (C,P).
          const int newRank   = 4;
          auto extents        = leftTransform.getExtents();
          auto variationTypes = leftTransform.getVariationTypes();

          composedTransform = leftTransform.shallowCopy(newRank, extents, variationTypes);
        }
          break;
        case 2: composedTransform = leftTransform; break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else if (rightTransform.isValid())
    {
      // leftTransform is the identity
      switch (rightRank)
      {
        case 4: composedTransform = rightTransform; break;
        case 3:
        {
          // - if right has the rank-3 transform, composedTransform should be a (C,P,1,D) object -- a 1 x D matrix at each (C,P).
          const int newRank   = 4;
          auto extents        = rightTransform.getExtents();
          auto variationTypes = rightTransform.getVariationTypes();
          extents[3]          = extents[2];
          variationTypes[3]   = variationTypes[2];
          extents[2]          = 1;
          variationTypes[2]   = CONSTANT;

          composedTransform = rightTransform.shallowCopy(newRank, extents, variationTypes);
        }
          break;
        case 2: composedTransform = rightTransform; break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else
    {
      // both left and right transforms are identity
      Kokkos::Array<ordinal_type,4> extents {basisValuesLeft.numCells(),basisValuesLeft.numPoints(),spaceDim,spaceDim};
      Kokkos::Array<DataVariationType,4> variationTypes {CONSTANT,CONSTANT,BLOCK_PLUS_DIAGONAL,BLOCK_PLUS_DIAGONAL};

      Kokkos::View<Scalar*,DeviceType> identityUnderlyingView("Intrepid2::PAMatrix - identity view",spaceDim);
      Kokkos::deep_copy(identityUnderlyingView, 1.0);
      composedTransform = Data<Scalar,DeviceType>(identityUnderlyingView,extents,variationTypes);
    }
    // allocate weighted transform
    _composedWeightedTransform = Impl::allocateComposedWeightedTransform<DeviceType,Scalar>(composedTransform,cellMeasures);
    auto composedWeightedTransform = _composedWeightedTransform; // avoid implicit reference to this
    // MARK: fill weighted transform container
    int rank = composedWeightedTransform.rank();
    int cellDataExtent    = composedWeightedTransform.getDataExtent(0);
    int numPoints         = composedWeightedTransform.getDataExtent(1);
    int d1_dim            = composedWeightedTransform.getDataExtent(2);
    int d2_dim            = composedWeightedTransform.getDataExtent(3);
    auto d1_variationType = composedWeightedTransform.getVariationTypes()[2];

    if (rank == 2)
    {
      Kokkos::Array<int,2> lowerBounds {0,0};
      Kokkos::Array<int,2> upperBounds {cellDataExtent,numPoints};
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>(lowerBounds, upperBounds);

      Kokkos::parallel_for("compute weighted transform", policy,
                           KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &pointOrdinal) {
        const Scalar & w = cellMeasures(cellDataOrdinal, pointOrdinal);
        Scalar & result  = composedWeightedTransform.getWritableEntry(cellDataOrdinal,pointOrdinal);
        result = w * composedTransform(cellDataOrdinal,pointOrdinal);
      });
    }
    else if ((rank == 3) || ((rank == 4) && (d1_variationType == BLOCK_PLUS_DIAGONAL)))
    {
      Kokkos::Array<int,3> lowerBounds {0,0,0};
      Kokkos::Array<int,3> upperBounds {cellDataExtent,numPoints,d1_dim};
      bool passThroughMatrixDims = (d1_variationType == BLOCK_PLUS_DIAGONAL); // if BLOCK_PLUS_DIAGONAL, it's a matrix, but everything is packed into the D1 dimension, and we want to sidestep the logic that tries to compute the matrix entry index based on (d1,d2) arguments.
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>(lowerBounds, upperBounds);

      Kokkos::parallel_for("compute weighted transform", policy,
                           KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &pointOrdinal, const int &d1) {
        const Scalar & w = cellMeasures(cellDataOrdinal, pointOrdinal);
        Scalar & result  = composedWeightedTransform.getWritableEntryWithPassThroughOption(passThroughMatrixDims,cellDataOrdinal,pointOrdinal,d1);
        result = w * composedTransform.getEntryWithPassThroughOption(passThroughMatrixDims,cellDataOrdinal,pointOrdinal,d1);
      });
    }
    else if (rank == 4)
    {
      Kokkos::Array<int,4> lowerBounds {0,0,0,0};
      Kokkos::Array<int,4> upperBounds {cellDataExtent,numPoints,d1_dim,d2_dim};
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>(lowerBounds, upperBounds);

      Kokkos::parallel_for("compute weighted transform", policy,
                           KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &pointOrdinal, const int &d1, const int &d2) {
        const Scalar & w = cellMeasures(cellDataOrdinal, pointOrdinal);
        Scalar & result  = composedWeightedTransform.getWritableEntry(cellDataOrdinal,pointOrdinal,d1,d2);
        const Scalar & originalTransform = composedTransform(cellDataOrdinal,pointOrdinal,d1,d2);
        result = w * originalTransform;
      });
    }
  }

  // MARK: Set up component integrations
  const int leftComponentCount  =  leftIsVectorValued ? basisValuesLeft. vectorData().numComponents() : 1;
  const int rightComponentCount = rightIsVectorValued ? basisValuesRight.vectorData().numComponents() : 1;

  // The 1D operators are packed once per distinct underlying view, so that passes sharing a component operator share
  // the packed view (assemble() groups passes by operator identity).  Left operators are packed (F,P[,D]); right ones (P[,D],F).
  std::map<std::tuple<const void*, size_t, bool>, View1D> packedOperators;
  auto packedOperator = [&](const auto &opData, const bool isRight) -> View1D
  {
    const auto key = std::make_tuple(static_cast<const void*>(opData.data()), size_t(opData.size()), isRight);
    auto entry = packedOperators.find(key);
    if (entry != packedOperators.end()) return entry->second;
    View1D opView = Impl::pack1DOperator<DeviceType,Scalar>(opData, isRight);
    packedOperators[key] = opView;
    return opView;
  };

  const int numCells  = _composedWeightedTransform.extent_int(0);
  const int numPoints = _composedWeightedTransform.extent_int(1);

  int leftFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
  for (int leftFamilyOrdinal=0; leftFamilyOrdinal<leftFamilyCount; leftFamilyOrdinal++)
  {
    const int leftFieldSpan = leftIsVectorValued ? basisValuesLeft.vectorData() .numFieldsInFamily(leftFamilyOrdinal)
                                                 : basisValuesLeft.basisValues().numFieldsInFamily(leftFamilyOrdinal);

    // "a" keeps track of the spatial dimension over which we are integrating in the left vector.
    // Components are allowed to span several dimensions; we keep track of the offset for the component in a_offset
    int a_offset = 0;
    for (int leftComponentOrdinal=0; leftComponentOrdinal<leftComponentCount; leftComponentOrdinal++)
    {
      TensorData<Scalar,DeviceType> leftComponent = leftIsVectorValued ? basisValuesLeft.vectorData().getComponent(leftFamilyOrdinal, leftComponentOrdinal)
                                                                       : basisValuesLeft.basisValues().tensorData(leftFamilyOrdinal);
      if (!leftComponent.isValid())
      {
         // represents zero
        a_offset += basisValuesLeft.vectorData().numDimsForComponent(leftComponentOrdinal);
        continue;
      }
      // left operators contract in the point (and space) dimensions
      std::vector<OpSpec> leftOperators;
      for (int r=0; r<leftComponent.numTensorComponents(); r++)
      {
        const auto opData = leftComponent.getTensorComponent(r).getUnderlyingView();
        leftOperators.push_back({packedOperator(opData, false), opData.extent_int(0), opData.extent_int(1)});
      }

      int rightFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
      for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
      {
        const int rightFieldSpan = rightIsVectorValued ? basisValuesRight.vectorData() .numFieldsInFamily(rightFamilyOrdinal)
                                                       : basisValuesRight.basisValues().numFieldsInFamily(rightFamilyOrdinal);
        // "b" keeps track of the spatial dimension over which we are integrating in the right vector
        // components are allowed to span several dimensions; we keep track of the offset for the component in b_offset
        int b_offset = 0;
        for (int rightComponentOrdinal=0; rightComponentOrdinal<rightComponentCount; rightComponentOrdinal++)
        {
          TensorData<Scalar,DeviceType> rightComponent =
             rightIsVectorValued ? basisValuesRight.vectorData().getComponent(rightFamilyOrdinal, rightComponentOrdinal)
                                 : basisValuesRight.basisValues().tensorData(rightFamilyOrdinal);
          if (!rightComponent.isValid())
          {
             // represents zero
            b_offset += basisValuesRight.vectorData().numDimsForComponent(rightComponentOrdinal);
            continue;
          }

          INTREPID2_TEST_FOR_EXCEPTION(leftComponent.numTensorComponents() != rightComponent.numTensorComponents(), std::invalid_argument, "left TensorData and right TensorData have different number of tensor components.  This is not supported.");

          // right operators contract in the field dimension
          std::vector<OpSpec> rightOperators;
          for (int r=0; r<rightComponent.numTensorComponents(); r++)
          {
            const auto opData = rightComponent.getTensorComponent(r).getUnderlyingView();
            rightOperators.push_back({packedOperator(opData, true), opData.extent_int(1), opData.extent_int(0)});
          }

          const int aSpan =  leftComponent.extent_int(2);
          const int bSpan = rightComponent.extent_int(2);

          PointDataSpec pointDataSpec{numCells, numPoints, a_offset, b_offset, aSpan, bSpan};

          if (_pointDataCache.find(pointDataSpec) == _pointDataCache.end())
          {
            const size_t pointDataSize = size_t(numCells) * numPoints * aSpan * bSpan;
            View1D pointDataView("pointDataView", pointDataSize);

            auto composedWeightedTransform = _composedWeightedTransform;

            // layout left: (C,P[,Da[,Db]])
            if (_composedWeightedTransform.rank() == 2) // (C,P): pointwise weight
            {
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numPoints});
              Kokkos::parallel_for("pack 1D pointData", policy,
                                   KOKKOS_LAMBDA(const int &cell, const int &pt)
                                   {
                pointDataView(cell + size_t(pt) * numCells) = composedWeightedTransform(cell,pt);
              });
            }
            else if (_composedWeightedTransform.rank() == 3) // (C,P,D): contract in b or expand in a
            {
              const bool contraction = (bSpan > 1);
              const int dOffset = contraction ? b_offset : a_offset;
              const int dSpan   = contraction ?    bSpan : aSpan;
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numPoints,dSpan});
              Kokkos::parallel_for("pack 1D pointData", policy,
                                   KOKKOS_LAMBDA(const int &cell, const int &pt, const int &d)
                                   {
                pointDataView(cell + (pt + size_t(d) * numPoints) * numCells) = composedWeightedTransform(cell,pt,dOffset + d);
              });
            }
            else if (_composedWeightedTransform.rank() == 4) // (C,P,D,D): contract in b and expand in a at each point
            {
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{numCells,numPoints,aSpan,bSpan});
              Kokkos::parallel_for("pack 1D pointData", policy,
                                   KOKKOS_LAMBDA(const int &cell, const int &pt, const int &da, const int &db)
                                   {
                pointDataView(cell + (pt + (da + size_t(db) * aSpan) * numPoints) * numCells) = composedWeightedTransform(cell,pt,a_offset + da,b_offset + db);
              });
            }
            _pointDataCache[pointDataSpec] = pointDataView;
          }
          {
            const int leftGroup  = leftFamilyOrdinal  * leftComponentCount  + leftComponentOrdinal;
            const int rightGroup = rightFamilyOrdinal * rightComponentCount + rightComponentOrdinal;
            componentIntegralGroups_.push_back({leftGroup, rightGroup});
          }
          componentIntegralsToSum_.push_back({leftOperators,pointDataSpec,rightOperators,
                                              leftFieldOrdinalOffset,leftFieldSpan,rightFieldOrdinalOffset,rightFieldSpan});

          b_offset += rightIsVectorValued ? basisValuesRight.vectorData().numDimsForComponent(rightComponentOrdinal) : 1;
        }
        rightFieldOrdinalOffset += rightFieldSpan;
      }
      a_offset += leftIsVectorValued ? basisValuesLeft.vectorData().numDimsForComponent(leftComponentOrdinal) : 1;
    }
    leftFieldOrdinalOffset += leftFieldSpan;
  }
  ExecutionSpace().fence();

  // set maxIntermediateSize_: the per-cell size required for intermediate computations, which is used to size the workspaces
  const size_t F2 = basisValuesRight.numFields();
  const size_t F1 = basisValuesLeft.numFields();
  maxIntermediateSize_ = std::max(F1,F2);

  for (const auto &entry : componentIntegralsToSum_)
  {
    const auto & leftOps       = std::get<0>(entry);
    const auto & pointDataSpec = std::get<1>(entry);
    const auto & rightOps      = std::get<2>(entry);

    // per-cell sizes of the intermediates in apply(): the right input span, then each right contraction (from the last
    // component), the point weighting, and each left contraction
    size_t perCellSize = std::get<6>(entry);
    for (int r=int(rightOps.size())-1; r>=0; r--)
    {
      perCellSize = perCellSize / rightOps[r].N * rightOps[r].M;
      maxIntermediateSize_ = std::max(perCellSize,maxIntermediateSize_);
    }

    perCellSize = perCellSize / pointDataSpec.bSpan * pointDataSpec.aSpan;
    maxIntermediateSize_ = std::max(perCellSize,maxIntermediateSize_);

    for (int r=int(leftOps.size())-1; r>=0; r--)
    {
      perCellSize = perCellSize / leftOps[r].N * leftOps[r].M;
      maxIntermediateSize_ = std::max(perCellSize,maxIntermediateSize_);
    }
    INTREPID2_TEST_FOR_EXCEPTION(perCellSize != size_t(std::get<4>(entry)), std::logic_error,
                                 "PAMatrix: inconsistent component operator extents");
  }
}

template<typename DeviceType,class Scalar>
Data<Scalar,DeviceType> PAMatrix<DeviceType,Scalar>::allocateMatrixStorage() const
{
  // Allocates a (C,F,F) container for storing integral data

  // With ordinal filters (e.g., Serendipity), the stored basis values report the filtered field counts, which is the
  // size of the matrix that assemble() produces.

  // determine cellDataExtent and variation type.  We currently support CONSTANT, MODULAR, and GENERAL as possible output variation types, depending on the inputs.
  // If cellMeasures has non-trivial tensor structure, the rank-1 cell Data object is the first component.
  // If cellMeasures has trivial tensor structure, then the first and only component has the cell index in its first dimension.
  // I.e., either way the relevant Data object is cellMeasures.getTensorComponent(0)
  const int CELL_DIM = 0;
  const auto cellMeasureData = _cellMeasures.getTensorComponent(0);

  DimensionInfo combinedCellDimInfo = cellMeasureData.getDimensionInfo(CELL_DIM);
  // transforms may be invalid, indicating an identity transform.  If so, it will not constrain the output at all.
  if (_basisValuesLeft.transform().isValid())
  {
    combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo, _basisValuesLeft.transform().getDimensionInfo(CELL_DIM));
  }
  if (_basisValuesRight.transform().isValid())
  {
    combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo, _basisValuesRight.transform().getDimensionInfo(CELL_DIM));
  }

  DataVariationType cellVariationType = combinedCellDimInfo.variationType;
  int cellDataExtent                  = combinedCellDimInfo.dataExtent;

  const int numCells       = _basisValuesLeft.numCells();
  const int numFieldsLeft  = _basisValuesLeft.numFields();
  const int numFieldsRight = _basisValuesRight.numFields();

  if (_hasNontrivialOrientations)
  {
    // orientations make the matrix cell-dependent even when geometry and weights are otherwise constant across cells
    cellVariationType = GENERAL;
    cellDataExtent = numCells;
  }

  Kokkos::Array<int,3> extents {numCells, numFieldsLeft, numFieldsRight};
  Kokkos::Array<DataVariationType,3> variationTypes {cellVariationType,GENERAL,GENERAL};

  if (cellVariationType != CONSTANT)
  {
    Kokkos::View<Scalar***,DeviceType> data("Intrepid2::PAMatrix matrix storage",cellDataExtent,numFieldsLeft,numFieldsRight);
    return Data<Scalar,DeviceType>(data, extents, variationTypes);
  }
  else
  {
    Kokkos::View<Scalar**,DeviceType> data("Intrepid2::PAMatrix matrix storage",numFieldsLeft,numFieldsRight);
    return Data<Scalar,DeviceType>(data, extents, variationTypes);
  }
} // allocateMatrixStorage()

template<typename DeviceType,class Scalar>
Data<Scalar,DeviceType> PAMatrix<DeviceType,Scalar>::allocateDiagonalStorage() const
{
  const int CELL_DIM = 0;
  const auto cellMeasureData = _cellMeasures.getTensorComponent(0);

  DimensionInfo combinedCellDimInfo = cellMeasureData.getDimensionInfo(CELL_DIM);
  if (_basisValuesLeft.transform().isValid())
  {
    combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo,
                                                _basisValuesLeft.transform().getDimensionInfo(CELL_DIM));
  }
  if (_basisValuesRight.transform().isValid())
  {
    combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo,
                                                _basisValuesRight.transform().getDimensionInfo(CELL_DIM));
  }

  DataVariationType cellVariationType = combinedCellDimInfo.variationType;
  int cellDataExtent                  = combinedCellDimInfo.dataExtent;

  const int numCells       = _basisValuesLeft.numCells();
  const int numFieldsLeft  = _basisValuesLeft.numFields();
  const int numFieldsRight = _basisValuesRight.numFields();
  const int numDiagonalEntries = std::min(numFieldsLeft, numFieldsRight);

  if (_hasNontrivialOrientations)
  {
    cellVariationType = GENERAL;
    cellDataExtent = numCells;
  }

  Kokkos::Array<int,2> extents {numCells, numDiagonalEntries};
  Kokkos::Array<DataVariationType,2> variationTypes {cellVariationType, GENERAL};

  if (cellVariationType != CONSTANT)
  {
    Kokkos::View<Scalar**,DeviceType> data("Intrepid2::PAMatrix diagonal storage", cellDataExtent, numDiagonalEntries);
    return Data<Scalar,DeviceType>(data, extents, variationTypes);
  }
  else
  {
    Kokkos::View<Scalar*,DeviceType> data("Intrepid2::PAMatrix diagonal storage", numDiagonalEntries);
    return Data<Scalar,DeviceType>(data, extents, variationTypes);
  }
}

template<typename DeviceType,class Scalar>
bool PAMatrix<DeviceType,Scalar>::supportsDirectAssembly() const
{
  if (componentIntegralsToSum_.size() == 0) return false;
  // the grouped paths share point buffers across passes, so all passes must have the same point count, number of
  // tensor components, and cell count
  const int numPoints           = std::get<1>(componentIntegralsToSum_[0]).P;
  const int numTensorComponents = int(std::get<2>(componentIntegralsToSum_[0]).size());
  const int numCells            = _basisValuesLeft.numCells();
  for (const auto &integralTuple : componentIntegralsToSum_)
  {
    const auto &leftOps       = std::get<0>(integralTuple);
    const PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
    const auto &rightOps      = std::get<2>(integralTuple);
    const int leftSpan        = std::get<4>(integralTuple);
    const int rightSpan       = std::get<6>(integralTuple);

    if ((pointDataSpec.P != numPoints) || (int(rightOps.size()) != numTensorComponents) || (pointDataSpec.C != numCells)) return false;
    if ((leftOps.size() != rightOps.size()) || (rightOps.size() == 0) ||
        (rightOps.size() > 8) || (pointDataSpec.aSpan != 1) || (pointDataSpec.bSpan != 1))
    {
      return false;
    }

    size_t rightFieldProduct = 1, pointProduct = 1, leftFieldProduct = 1, leftPointProduct = 1;
    for (size_t r=0; r<rightOps.size(); r++)
    {
      rightFieldProduct *= rightOps[r].N;
      pointProduct      *= rightOps[r].M;
      leftFieldProduct  *= leftOps[r].M;
      leftPointProduct  *= leftOps[r].N;
    }
    if ((rightFieldProduct != size_t(rightSpan)) || (leftFieldProduct != size_t(leftSpan)) ||
        (pointProduct != size_t(pointDataSpec.P)) || (leftPointProduct != size_t(pointDataSpec.P)))
    {
      return false;
    }
  }
  return true;
}

//! Deprecated: use PAMatrix::supportsDirectAssembly().
template<typename DeviceType,class Scalar>
bool pamatrixAssembleUsesUnitColumnPath(const PAMatrix<DeviceType,Scalar> &matrix)
{
  return matrix.supportsDirectAssembly();
}

/** \brief Direct sum-factorized assembly (of the unoriented matrix).

    For one integration pass (scalar point weights W, left operators L_r of size F_r x Q_r, right operators R_r of size
    Q_r x G_r, r = 0..R-1, component 0 fastest in both field and point ordinals),
      M(c, I, J) = sum_p W(c,p) prod_r L_r(i_r, p_r) R_r(p_r, j_r).
    With the pair operators P_r(q, (i,j)) = L_r(i,q) R_r(q,j), stage r maps a tensor stored as (p_r, X) to (X, (i_r j_r))
    by the GEMM T' = T^T P_r.  Starting from W stored as (p_0, ..., p_{R-1}, c), after R-1 stages the tensor is
    (p_{R-1}, c, (i_0 j_0), ..., (i_{R-2} j_{R-2})).

    Passes are grouped by their last-component operator pair (L_{R-1}, R_{R-1}); within a group, the stage-(R-1) inputs
    are summed, so that each group needs one last-stage contraction.  The last stage either
    - (unfiltered, or dense filters) accumulates all groups' (c, ..., (i_{R-1} j_{R-1})) results by GEMM into one
      buffer, which a single permutation writes to the output; or
    - (sparse ordinal filters, e.g. Serendipity) computes each kept output entry directly as a length-Q_{R-1} dot product
      against the group's pair operator, skipping the full tensor-product intermediate.
    Cells are processed in worksets sized to the memory budget in the options. */
template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assembleDirect(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using OrdView1D = Kokkos::View<ordinal_type*,DeviceType>;

  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();
  const bool hasLeftFilter  = leftOrdinalFilter_ .extent_int(0) > 0;
  const bool hasRightFilter = rightOrdinalFilter_.extent_int(0) > 0;

  const int outputRank = integrals.getUnderlyingViewRank();
  INTREPID2_TEST_FOR_EXCEPTION((outputRank != 2) && (outputRank != 3), std::invalid_argument,
      "PAMatrix::assemble() supports only rank-2 or rank-3 underlying matrix storage");
  // only the cells for which the output stores data are computed (all, for rank-3 GENERAL storage; the first
  // getDataExtent(0) for MODULAR; one, for rank-2 storage, which is constant in the cell)
  const ordinal_type numCellsComputed = (outputRank == 3) ? integrals.getDataExtent(0) : 1;

  const auto &passes = componentIntegralsToSum_;
  const int numPasses = int(passes.size());
  const int R = int(std::get<2>(passes[0]).size());

  // full-ordinal row/column -> output row/column (or -1 if filtered out); with family offsets, a pass's full ordinal is
  // offset + (tensor ordinal within the family)
  auto hostLeftFilter  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), leftOrdinalFilter_);
  auto hostRightFilter = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rightOrdinalFilter_);
  const ordinal_type numFullLeft  = hasLeftFilter  ? numUnfilteredFieldsLeft_  : numFieldsLeft;
  const ordinal_type numFullRight = hasRightFilter ? numUnfilteredFieldsRight_ : numFieldsRight;
  std::vector<ordinal_type> rowForFull(numFullLeft, -1), colForFull(numFullRight, -1);
  for (ordinal_type row=0; row<numFieldsLeft;  row++) rowForFull[hasLeftFilter  ? hostLeftFilter(row)  : row] = row;
  for (ordinal_type col=0; col<numFieldsRight; col++) colForFull[hasRightFilter ? hostRightFilter(col) : col] = col;

  // group passes by (left family span/offset, right family span/offset, last-component operator pair): within a group,
  // the last-stage contraction is shared.  (init() packs each distinct 1D operator once, so operator identity can be
  // compared by pointer.)
  struct Group { std::vector<int> passes; };
  std::vector<Group> groups;
  {
    auto sameOp = [](const OpSpec &a, const OpSpec &b) { return (a.opView.data() == b.opView.data()) && (a.M == b.M) && (a.N == b.N); };
    for (int pass=0; pass<numPasses; pass++)
    {
      const auto &t = passes[pass];
      int found = -1;
      for (int g=0; g<int(groups.size()); g++)
      {
        const auto &u = passes[groups[g].passes[0]];
        if ((std::get<3>(t) == std::get<3>(u)) && (std::get<5>(t) == std::get<5>(u)) &&
            (std::get<4>(t) == std::get<4>(u)) && (std::get<6>(t) == std::get<6>(u)) &&
            sameOp(std::get<0>(t)[R-1], std::get<0>(u)[R-1]) && sameOp(std::get<2>(t)[R-1], std::get<2>(u)[R-1]))
        {
          bool sameShape = true;
          for (int r=0; r<R-1; r++)
          {
            sameShape = sameShape && (std::get<0>(t)[r].M == std::get<0>(u)[r].M) && (std::get<2>(t)[r].N == std::get<2>(u)[r].N);
          }
          if (sameShape) { found = g; break; }
        }
      }
      if (found < 0) { groups.push_back(Group()); found = int(groups.size()) - 1; }
      groups[found].passes.push_back(pass);
    }
  }
  bool anyGroupHasMultiplePasses = false;
  for (const auto &group : groups) anyGroupHasMultiplePasses = anyGroupHasMultiplePasses || (group.passes.size() > 1);

  // per-group geometry: extents, and row/column offset tables over the pair (i,j) index space of the first R-1
  // components and the last component; output row = rowForFull[leftOffset + I], I = I_prefix + i_{R-1} * strideLast
  struct GroupLayout
  {
    int leftOffset, rightOffset;
    std::vector<int> F, G, Q;      // per component
    size_t prefixPairs = 1;        // prod_{r<R-1} F_r G_r
    int prefixLeft = 1, prefixRight = 1;
    View1D pairOpLast;             // (Q_{R-1}, F_{R-1} G_{R-1})
    OrdView1D prefixRowOffset, prefixColOffset; // per prefix pair: I_prefix, J_prefix
    OrdView1D rowMap, colMap;      // full ordinal within the family -> output row/column (or -1)
  };
  auto pairOperator = [](const OpSpec &leftOp, const OpSpec &rightOp) -> View1D
  {
    const int F = leftOp.M, Q = leftOp.N, Gn = rightOp.N;
    View1D P(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble pair operator"), size_t(Q) * F * Gn);
    const auto Lop = leftOp.opView;
    const auto Rop = rightOp.opView;
    Kokkos::parallel_for("PAMatrix::assemble(): pair operator", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Q,F,Gn}),
    KOKKOS_LAMBDA(const ordinal_type &q, const ordinal_type &i, const ordinal_type &j) {
      P(q + (i + size_t(j) * F) * Q) = Lop(i + size_t(q) * F) * Rop(q + size_t(j) * Q);
    });
    return P;
  };

  std::vector<GroupLayout> layouts(groups.size());
  // groups that share family offsets and component shapes have identical (c, prefix pairs, (i_{R-1} j_{R-1})) layouts,
  // so their dense last-stage results can be accumulated (beta = 1) before a single write
  std::vector<int> writeFamily(groups.size(), -1);
  int numWriteFamilies = 0;
  for (int g=0; g<int(groups.size()); g++)
  {
    const auto &t = passes[groups[g].passes[0]];
    auto &L = layouts[g];
    L.leftOffset  = std::get<3>(t);
    L.rightOffset = std::get<5>(t);
    for (int r=0; r<R; r++)
    {
      L.F.push_back(std::get<0>(t)[r].M);
      L.Q.push_back(std::get<0>(t)[r].N);
      L.G.push_back(std::get<2>(t)[r].N);
    }
    for (int r=0; r<R-1; r++) { L.prefixPairs *= size_t(L.F[r]) * L.G[r]; L.prefixLeft *= L.F[r]; L.prefixRight *= L.G[r]; }
    OrdView1D rowOffsets("PAMatrix assemble prefix row offsets", L.prefixPairs), colOffsets("PAMatrix assemble prefix col offsets", L.prefixPairs);
    auto rowHost = Kokkos::create_mirror_view(rowOffsets);
    auto colHost = Kokkos::create_mirror_view(colOffsets);
    for (size_t pair=0; pair<L.prefixPairs; pair++)
    {
      size_t remainder = pair;
      int I = 0, J = 0, leftStride = 1, rightStride = 1;
      for (int r=0; r<R-1; r++)
      {
        const int ij = int(remainder % (size_t(L.F[r]) * L.G[r]));
        remainder /= (size_t(L.F[r]) * L.G[r]);
        I += (ij % L.F[r]) * leftStride;
        J += (ij / L.F[r]) * rightStride;
        leftStride *= L.F[r];
        rightStride *= L.G[r];
      }
      rowHost(pair) = I;
      colHost(pair) = J;
    }
    Kokkos::deep_copy(rowOffsets, rowHost);
    Kokkos::deep_copy(colOffsets, colHost);
    L.prefixRowOffset = rowOffsets;
    L.prefixColOffset = colOffsets;
    L.pairOpLast = pairOperator(std::get<0>(t)[R-1], std::get<2>(t)[R-1]);

    // device copies of the full->output maps, offset to this group's family
    L.rowMap = OrdView1D("PAMatrix assemble row map", size_t(L.prefixLeft) * L.F[R-1]);
    L.colMap = OrdView1D("PAMatrix assemble col map", size_t(L.prefixRight) * L.G[R-1]);
    auto rowMapHost = Kokkos::create_mirror_view(L.rowMap);
    auto colMapHost = Kokkos::create_mirror_view(L.colMap);
    for (int I=0; I<L.prefixLeft*L.F[R-1];  I++) rowMapHost(I) = rowForFull[L.leftOffset  + I];
    for (int J=0; J<L.prefixRight*L.G[R-1]; J++) colMapHost(J) = colForFull[L.rightOffset + J];
    Kokkos::deep_copy(L.rowMap, rowMapHost);
    Kokkos::deep_copy(L.colMap, colMapHost);
  }
  for (int g=0; g<int(groups.size()); g++)
  {
    for (int h=0; h<g; h++)
    {
      const auto &A = layouts[g], &B = layouts[h];
      if ((A.leftOffset == B.leftOffset) && (A.rightOffset == B.rightOffset) && (A.F == B.F) && (A.G == B.G))
      {
        writeFamily[g] = writeFamily[h];
        break;
      }
    }
    if (writeFamily[g] < 0) writeFamily[g] = numWriteFamilies++;
  }
  {
    // order the groups so that each write family's groups are contiguous: their dense last-stage GEMMs accumulate into a
    // single buffer, which another family's GEMM (beta = 0) must not overwrite before the family's result is written
    std::vector<int> order(groups.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&writeFamily](const int a, const int b) { return writeFamily[a] < writeFamily[b]; });
    std::vector<Group> sortedGroups;
    std::vector<GroupLayout> sortedLayouts;
    std::vector<int> sortedWriteFamily;
    for (const int g : order)
    {
      sortedGroups.push_back(groups[g]);
      sortedLayouts.push_back(layouts[g]);
      sortedWriteFamily.push_back(writeFamily[g]);
    }
    groups      = sortedGroups;
    layouts     = sortedLayouts;
    writeFamily = sortedWriteFamily;
  }

  // Sparse-filter path: when few (row, column) pairs survive, compute each kept entry directly at the last stage.
  bool useSparseLast = false;
  {
    const bool hasFilter = hasLeftFilter || hasRightFilter;
    const double keptFraction = double(numFieldsLeft) * numFieldsRight / (double(numFullLeft) * numFullRight);
    switch (options.sparseLastStage)
    {
      case PAMatrixAssemblyOptions::SparseLastStage::Automatic: useSparseLast = hasFilter && (keptFraction < 0.25); break;
      case PAMatrixAssemblyOptions::SparseLastStage::Always:    useSparseLast = hasFilter; break;
      case PAMatrixAssemblyOptions::SparseLastStage::Never:     useSparseLast = false; break;
    }
  }

  // kept (prefix pair, last-component (i,j)) -> output (row, col), per group, for the sparse last stage
  std::vector<OrdView1D> keptPairV(groups.size()), keptIJV(groups.size()), keptRowV(groups.size()), keptColV(groups.size());
  if (useSparseLast)
  {
    for (int g=0; g<int(groups.size()); g++)
    {
      const auto &L = layouts[g];
      const int FLast = L.F[R-1], GLast = L.G[R-1];
      std::vector<ordinal_type> keptPair, keptIJ, keptRow, keptCol;
      auto rowHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), L.rowMap);
      auto colHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), L.colMap);
      auto prefixRowHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), L.prefixRowOffset);
      auto prefixColHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), L.prefixColOffset);
      for (size_t pair=0; pair<L.prefixPairs; pair++)
      {
        for (int jLast=0; jLast<GLast; jLast++)
        {
          const ordinal_type col = colHost(prefixColHost(pair) + jLast * L.prefixRight);
          if (col < 0) continue;
          for (int iLast=0; iLast<FLast; iLast++)
          {
            const ordinal_type row = rowHost(prefixRowHost(pair) + iLast * L.prefixLeft);
            if (row < 0) continue;
            keptPair.push_back(ordinal_type(pair)); keptIJ.push_back(iLast + jLast * FLast); keptRow.push_back(row); keptCol.push_back(col);
          }
        }
      }
      const int numKept = int(keptPair.size());
      keptPairV[g] = OrdView1D("kept pair", numKept); keptIJV[g] = OrdView1D("kept ij", numKept);
      keptRowV[g]  = OrdView1D("kept row",  numKept); keptColV[g] = OrdView1D("kept col", numKept);
      auto a = Kokkos::create_mirror_view(keptPairV[g]); auto b = Kokkos::create_mirror_view(keptIJV[g]);
      auto c = Kokkos::create_mirror_view(keptRowV[g]);  auto d = Kokkos::create_mirror_view(keptColV[g]);
      for (int k=0; k<numKept; k++) { a(k) = keptPair[k]; b(k) = keptIJ[k]; c(k) = keptRow[k]; d(k) = keptCol[k]; }
      Kokkos::deep_copy(keptPairV[g], a); Kokkos::deep_copy(keptIJV[g], b); Kokkos::deep_copy(keptRowV[g], c); Kokkos::deep_copy(keptColV[g], d);
    }
  }

  // workset size from the memory budget: the largest per-cell intermediate is max over stages of (current tensor size)
  size_t perCellEntries = 0;
  for (const auto &L : layouts)
  {
    size_t size = 1;
    for (int r=0; r<R; r++) size *= L.Q[r];
    size_t maxSize = size;
    for (int r=0; r<R-1; r++)
    {
      size = size / L.Q[r] * L.F[r] * L.G[r];
      maxSize = std::max(maxSize, size);
    }
    if (!useSparseLast) maxSize = std::max(maxSize, size / L.Q[R-1] * L.F[R-1] * L.G[R-1]);
    perCellEntries = std::max(perCellEntries, maxSize);
  }
  // ping and pong buffers, plus the group accumulator (if any group has more than one pass) and the dense last-stage accumulator
  const int numBuffers = 2 + (anyGroupHasMultiplePasses ? 1 : 0) + (useSparseLast ? 0 : 1);
  const ordinal_type worksetSize = pamatrixWorksetSize<Scalar>(options.workspaceMB, numBuffers * perCellEntries, numCellsComputed);

  const size_t bufferSize = perCellEntries * size_t(worksetSize);
  View1D bufferA(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct A"), bufferSize);
  View1D bufferB(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct B"), bufferSize);
  View1D groupSumBuffer;
  if (anyGroupHasMultiplePasses) groupSumBuffer = View1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct group sum"), bufferSize);
  View1D lastAccumulator;
  if (!useSparseLast) lastAccumulator = View1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct last"), bufferSize);

  if (outputRank == 3) Kokkos::deep_copy(integrals.template getUnderlyingView<3>(), Scalar(0));
  else                 Kokkos::deep_copy(integrals.template getUnderlyingView<2>(), Scalar(0));

  // pair operators for the first R-1 components, per pass
  std::vector<std::vector<View1D>> pairOps(numPasses);
  for (int pass=0; pass<numPasses; pass++)
  {
    const auto &t = passes[pass];
    for (int r=0; r<R-1; r++) pairOps[pass].push_back(pairOperator(std::get<0>(t)[r], std::get<2>(t)[r]));
  }

  for (ordinal_type cellStart=0; cellStart<numCellsComputed; cellStart+=worksetSize)
  {
    const ordinal_type Cw = std::min(worksetSize, numCellsComputed - cellStart);

    for (int g=0; g<int(groups.size()); g++)
    {
      const auto &L = layouts[g];
      const int QLast = L.Q[R-1];
      const size_t lastInputSize = size_t(QLast) * Cw * L.prefixPairs; // (p_{R-1}, c, prefix pairs)
      // with a single pass, the stage-(R-1) input stays in the ping/pong buffers; otherwise, the passes' inputs are summed
      const bool singlePass = (groups[g].passes.size() == 1);
      View1D groupSum = singlePass ? View1D() : groupSumBuffer;

      bool firstPass = true;
      for (int pass : groups[g].passes)
      {
        const auto &t = passes[pass];
        const auto &pointDataSpec = std::get<1>(t);
        const int numPointsTotal = pointDataSpec.P;
        const int C = pointDataSpec.C;
        const auto pointData = _pointDataCache.find(pointDataSpec)->second;
        // stage-0 input (p_0, ..., p_{R-1}, c) into bufferA
        {
          auto dest = bufferA;
          Kokkos::parallel_for("PAMatrix::assemble(): transpose point weights",
                               Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{Cw,numPointsTotal}),
          KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &point) {
            const size_t dataCell = (C == 1) ? 0 : size_t(cell + cellStart);
            dest(point + size_t(cell) * numPointsTotal) = pointData(dataCell + size_t(point) * C);
          });
        }
        View1D current = bufferA, other = bufferB;
        size_t currentSize = size_t(numPointsTotal) * Cw;
        for (int r=0; r<R-1; r++)
        {
          const int Q = L.Q[r], FG = L.F[r] * L.G[r];
          const size_t X = currentSize / Q;
          INTREPID2_TEST_FOR_EXCEPTION(X > size_t(std::numeric_limits<ordinal_type>::max()), std::invalid_argument,
                                       "PAMatrix::assemble(): GEMM dimension overflow; reduce the workspace budget");
          const bool sumIntoGroup = (r == R-2) && !singlePass;
          View1D target = sumIntoGroup ? groupSum : other;
          const Scalar beta = (sumIntoGroup && !firstPass) ? Scalar(1) : Scalar(0);
          Impl::gemm<DeviceType,Scalar>('T', 'N', ordinal_type(X), FG, Q, Scalar(1), current.data(), Q,
                                        pairOps[pass][r].data(), beta, target.data());
          currentSize = X * FG;
          if (!sumIntoGroup) std::swap(current, other);
        }
        if (singlePass)
        {
          groupSum = current; // the single pass's stage-(R-1) input
        }
        else if (R == 1)
        {
          // no prefix stages: the group sum is the weights themselves
          auto src = bufferA; auto dst = groupSum;
          const bool accumulate = !firstPass;
          Kokkos::parallel_for("PAMatrix::assemble(): sum weights", Kokkos::RangePolicy<ExecutionSpace>(0, currentSize),
          KOKKOS_LAMBDA(const size_t &k) { dst(k) = accumulate ? dst(k) + src(k) : src(k); });
        }
        firstPass = false;
      }
      // groupSum: (p_{R-1}, c, prefix pairs), size lastInputSize

      const int FLast = L.F[R-1], GLast = L.G[R-1];
      const int leftStrideLast = L.prefixLeft, rightStrideLast = L.prefixRight;
      const size_t prefixPairs = L.prefixPairs;
      auto prefixRow = L.prefixRowOffset;
      auto prefixCol = L.prefixColOffset;
      auto rowMap = L.rowMap;
      auto colMap = L.colMap;

      if (!useSparseLast)
      {
        // last stage: (p_{R-1}, [c, prefix pairs]) -> ([c, prefix pairs], (i_{R-1} j_{R-1}))
        const size_t X = lastInputSize / QLast;
        INTREPID2_TEST_FOR_EXCEPTION(X > size_t(std::numeric_limits<ordinal_type>::max()), std::invalid_argument,
                                     "PAMatrix::assemble(): GEMM dimension overflow; reduce the workspace budget");
        bool firstInFamily = true, lastInFamily = true;
        for (int h=0; h<g; h++) if (writeFamily[h] == writeFamily[g]) firstInFamily = false;
        for (int h=g+1; h<int(groups.size()); h++) if (writeFamily[h] == writeFamily[g]) lastInFamily = false;
        Impl::gemm<DeviceType,Scalar>('T', 'N', ordinal_type(X), FLast * GLast, QLast, Scalar(1), groupSum.data(), QLast,
                                      L.pairOpLast.data(), firstInFamily ? Scalar(0) : Scalar(1), lastAccumulator.data());
        if (!lastInFamily) continue;
        auto source = lastAccumulator;
        // cell outermost, so that writes stay within one cell's matrix.  (Separate lambdas for the two output ranks:
        // CUDA extended lambdas cannot be defined in generic lambdas or capture in constexpr-if contexts.)
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw, FLast*GLast, ordinal_type(prefixPairs)});
        if (outputRank == 3)
        {
          auto matrixView = integrals.template getUnderlyingView<3>();
          Kokkos::parallel_for("PAMatrix::assemble(): write direct contraction", policy,
          KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &ijLast, const ordinal_type &pair) {
            const int iLast = ijLast % FLast, jLast = ijLast / FLast;
            const ordinal_type row = rowMap(prefixRow(pair) + iLast * leftStrideLast);
            const ordinal_type col = colMap(prefixCol(pair) + jLast * rightStrideLast);
            if ((row < 0) || (col < 0)) return;
            matrixView(cell + cellStart, row, col) += source(cell + (pair + size_t(ijLast) * prefixPairs) * Cw);
          });
        }
        else
        {
          auto matrixView = integrals.template getUnderlyingView<2>();
          Kokkos::parallel_for("PAMatrix::assemble(): write direct contraction", policy,
          KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &ijLast, const ordinal_type &pair) {
            const int iLast = ijLast % FLast, jLast = ijLast / FLast;
            const ordinal_type row = rowMap(prefixRow(pair) + iLast * leftStrideLast);
            const ordinal_type col = colMap(prefixCol(pair) + jLast * rightStrideLast);
            if ((row < 0) || (col < 0)) return;
            matrixView(row, col) += source(cell + (pair + size_t(ijLast) * prefixPairs) * Cw);
          });
        }
      }
      else
      {
        // sparse last stage: for each kept (row, col) of this group, value(c) = sum_q groupSum(q, c, prefixPair) P(q, ij)
        const int numKept = keptPairV[g].extent_int(0);
        if (numKept == 0) continue;
        auto keptPair = keptPairV[g]; auto keptIJ = keptIJV[g]; auto keptRow = keptRowV[g]; auto keptCol = keptColV[g];
        auto source = groupSum;
        auto P = L.pairOpLast;
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numKept, Cw});
        if (outputRank == 3)
        {
          auto matrixView = integrals.template getUnderlyingView<3>();
          Kokkos::parallel_for("PAMatrix::assemble(): sparse last stage", policy,
          KOKKOS_LAMBDA(const ordinal_type &k, const ordinal_type &cell) {
            const size_t sOffset = size_t(QLast) * (cell + size_t(keptPair(k)) * Cw);
            const size_t pOffset = size_t(QLast) * keptIJ(k);
            Scalar value = 0;
            for (int q=0; q<QLast; q++) value += source(sOffset + q) * P(pOffset + q);
            matrixView(cell + cellStart, keptRow(k), keptCol(k)) += value;
          });
        }
        else
        {
          auto matrixView = integrals.template getUnderlyingView<2>();
          Kokkos::parallel_for("PAMatrix::assemble(): sparse last stage", policy,
          KOKKOS_LAMBDA(const ordinal_type &k, const ordinal_type &cell) {
            const size_t sOffset = size_t(QLast) * (cell + size_t(keptPair(k)) * Cw);
            const size_t pOffset = size_t(QLast) * keptIJ(k);
            Scalar value = 0;
            for (int q=0; q<QLast; q++) value += source(sOffset + q) * P(pOffset + q);
            matrixView(keptRow(k), keptCol(k)) += value;
          });
        }
      }
    }
  }
  ExecutionSpace().fence();
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assembleUnitColumns(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using ScratchView = Kokkos::View<Scalar*,DeviceType,Kokkos::MemoryUnmanaged>;

  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();
  const int outputRank = integrals.getUnderlyingViewRank();
  const ordinal_type numCellsComputed = (outputRank == 3) ? integrals.getDataExtent(0) : 1;

  // with ordinal filters, block column j is full right ordinal rightOrdinalFilter_(columnStart+j), and the left
  // contraction produces all numUnfilteredFieldsLeft_ rows, of which rows leftOrdinalFilter_(i) are kept
  const bool hasLeftFilter  = leftOrdinalFilter_ .extent_int(0) > 0;
  const bool hasRightFilter = rightOrdinalFilter_.extent_int(0) > 0;
  const auto leftOrdinalFilter  = leftOrdinalFilter_;
  const auto rightOrdinalFilter = rightOrdinalFilter_;
  const ordinal_type numLeftRowsComputed = hasLeftFilter ? numUnfilteredFieldsLeft_ : numFieldsLeft;

  // Passes are grouped by their left and right (family, component).  Passes with the same right group share their
  // right operators, so the unit columns' point values are computed once per right group; passes with the same left
  // group share their left operators, so their weighted point values are summed before a single left contraction.
  const int numIntegrationPasses = int(componentIntegralsToSum_.size());
  const int numTensorComponents  = int(std::get<2>(componentIntegralsToSum_[0]).size());
  const int numPoints            = std::get<1>(componentIntegralsToSum_[0]).P;
  int numLeftGroups = 0, numRightGroups = 0;
  for (const auto &groups : componentIntegralGroups_)
  {
    numLeftGroups  = std::max(numLeftGroups,  groups.first  + 1);
    numRightGroups = std::max(numRightGroups, groups.second + 1);
  }

  const ordinal_type maxColumnBlockSize = 64;
  const ordinal_type columnBlockSize = std::max<ordinal_type>(1, std::min(numFieldsRight, maxColumnBlockSize));

  // per cell: unit-column point values for each right group, two contraction workspaces, and the output columns
  const size_t perCellEntries = size_t(columnBlockSize) * (size_t(numPoints) * numRightGroups + 2 * maxIntermediateSize_ + numLeftRowsComputed);
  const ordinal_type worksetSize = pamatrixWorksetSize<Scalar>(options.workspaceMB, perCellEntries, numCellsComputed);

  for (ordinal_type columnStart=0; columnStart<numFieldsRight; columnStart+=columnBlockSize)
  {
    const ordinal_type blockSize = std::min(columnBlockSize, numFieldsRight - columnStart);
    for (ordinal_type cellStart=0; cellStart<numCellsComputed; cellStart+=worksetSize)
    {
      const ordinal_type Cw = std::min(worksetSize, numCellsComputed - cellStart);
      Kokkos::View<Scalar***,DeviceType> outputColumns("PAMatrix assemble output columns", Cw, numLeftRowsComputed, blockSize);

      // point values of this block's unit columns, one buffer per right group, each with layout (P,N,Cw) as consumed by
      // the weighting below
      const size_t pointBufferSize = size_t(numPoints) * blockSize * Cw;
      View1D rightFills(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble unit-column point values"), pointBufferSize * numRightGroups);
      std::vector<bool> rightGroupFilled(numRightGroups, false);
      for (int integrationPass=0; integrationPass<numIntegrationPasses; integrationPass++)
      {
        const int rightGroup = componentIntegralGroups_[integrationPass].second;
        if (rightGroupFilled[rightGroup]) continue;
        rightGroupFilled[rightGroup] = true;

        const auto &integralTuple = componentIntegralsToSum_[integrationPass];
        const auto &rightOps      = std::get<2>(integralTuple);
        const int rightOrdinalOffset = std::get<5>(integralTuple);
        const int rightInputSpan     = std::get<6>(integralTuple);

        Kokkos::Array<OpSpec,8> rightOpsArray;
        Kokkos::Array<int,8> rightFieldExtents;
        Kokkos::Array<int,8> pointExtents;
        for (int r=0; r<numTensorComponents; r++)
        {
          rightOpsArray[r] = rightOps[r];
          rightFieldExtents[r] = rightOps[r].N;
          pointExtents[r]      = rightOps[r].M;
        }
        auto fill = Kokkos::subview(rightFills, Kokkos::make_pair(rightGroup * pointBufferSize, (rightGroup + 1) * pointBufferSize));

        // Point values of the unit trial columns: for column (f_0,...,f_{R-1}) and point (p_0,...,p_{R-1}), the value is
        // prod_r op_r(p_r, f_r).  Component 0 varies fastest in both column and point ordinals, so each thread handles one
        // (cell, column, p_1..p_{R-1}) and writes the contiguous run over p_0, after computing the product of the other
        // components once.  (The 1D operators do not depend on the cell.)
        const int numPoints0     = pointExtents[0];
        const int numPointsOuter = numPoints / numPoints0;
        auto unitColumnPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,blockSize,numPointsOuter});
        Kokkos::parallel_for("PAMatrix::assemble(): unit trial columns to points", unitColumnPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &blockColumn, const ordinal_type &pointOuter)
        {
          const int globalColumn = hasRightFilter ? rightOrdinalFilter(columnStart + blockColumn) : columnStart + blockColumn;
          const size_t dest = size_t(numPoints0) * pointOuter + (blockColumn + size_t(cell) * blockSize) * numPoints;
          if ((globalColumn < rightOrdinalOffset) || (globalColumn >= rightOrdinalOffset + rightInputSpan))
          {
            for (int p0=0; p0<numPoints0; p0++) fill(dest + p0) = Scalar(0);
            return;
          }
          int columnRemainder = globalColumn - rightOrdinalOffset;
          const int field0 = columnRemainder % rightFieldExtents[0];
          columnRemainder /= rightFieldExtents[0];
          int pointRemainder = pointOuter;
          Scalar outerProduct = Scalar(1);
          for (int r=1; r<numTensorComponents; r++)
          {
            const int field = columnRemainder % rightFieldExtents[r];
            const int point = pointRemainder % pointExtents[r];
            columnRemainder /= rightFieldExtents[r];
            pointRemainder  /= pointExtents[r];
            const auto op = rightOpsArray[r];
            outerProduct *= op.opView(point + size_t(field) * op.M);
          }
          const auto op0 = rightOpsArray[0];
          for (int p0=0; p0<numPoints0; p0++) fill(dest + p0) = outerProduct * op0.opView(p0 + size_t(field0) * op0.M);
        });
      }

      const size_t workspaceSize = maxIntermediateSize_ * size_t(Cw) * blockSize;
      View1D workspace(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble unit-column workspace"), 2 * workspaceSize);
      ScratchView workspace1(workspace.data(), workspaceSize);
      ScratchView workspace2(workspace.data() + workspaceSize, workspaceSize);
      // the left contraction starts from pointDataOut; keep the parity used by apply() so the GEMM chain alternates the same way
      ScratchView pointDataOut = (numTensorComponents % 2 == 0) ? workspace2 : workspace1;

      for (int leftGroup=0; leftGroup<numLeftGroups; leftGroup++)
      {
        // sum, over this left group's passes, of (point weights) * (unit-column point values), in layout (N,Cw,P)
        int representativePass = -1;
        Kokkos::deep_copy(Kokkos::subview(pointDataOut, Kokkos::make_pair(size_t(0), pointBufferSize)), Scalar(0));
        for (int integrationPass=0; integrationPass<numIntegrationPasses; integrationPass++)
        {
          if (componentIntegralGroups_[integrationPass].first != leftGroup) continue;
          representativePass = integrationPass;
          const PointDataSpec &pointDataSpec = std::get<1>(componentIntegralsToSum_[integrationPass]);
          const int rightGroup = componentIntegralGroups_[integrationPass].second;
          auto fill = Kokkos::subview(rightFills, Kokkos::make_pair(rightGroup * pointBufferSize, (rightGroup + 1) * pointBufferSize));
          const auto weights = _pointDataCache.find(pointDataSpec)->second;
          auto accumulated = pointDataOut;
          const int C = pointDataSpec.C;
          const int P = pointDataSpec.P;
          const int N = blockSize;
          // accumulated(n,cell,point) += weights(cell,point) * fill(point,n,cell); n fastest to match the output layout
          auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{P,Cw,N});
          Kokkos::parallel_for("PAMatrix::assemble(): accumulate weighted unit-column point values", policy,
          KOKKOS_LAMBDA(const ordinal_type &point, const ordinal_type &cell, const ordinal_type &n)
          {
            const size_t dataCell = (C == 1) ? 0 : size_t(cell + cellStart);
            accumulated(n + (cell + size_t(point) * Cw) * N) += weights(dataCell + size_t(point) * C) * fill(point + (n + size_t(cell) * N) * P);
          });
        }
        if (representativePass < 0) continue;

        const auto &integralTuple = componentIntegralsToSum_[representativePass];
        const auto &leftOps       = std::get<0>(integralTuple);
        const PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
        const int leftOrdinalOffset  = std::get<3>(integralTuple);
        const int leftOutputSpan     = std::get<4>(integralTuple);

        size_t Nr = size_t(pointDataSpec.P) * Cw * blockSize * pointDataSpec.aSpan;
        const int numLeftIntegrals = int(leftOps.size());
        for (int i=0; i<numLeftIntegrals; i++)
        {
          auto in  = ((i+numTensorComponents+1)%2 == 0) ? workspace1 : workspace2;
          auto out = ((i+numTensorComponents+1)%2 == 0) ? workspace2 : workspace1;

          const int r = numLeftIntegrals-i-1;
          const auto &op = leftOps[r];
          const ordinal_type Fr = op.M;
          const ordinal_type Pr = op.N;
          INTREPID2_TEST_FOR_EXCEPTION(Nr % Pr != 0, std::invalid_argument, "Error: Nr must be a multiple of Pr");
          INTREPID2_TEST_FOR_EXCEPTION(Nr / Pr > size_t(std::numeric_limits<ordinal_type>::max()), std::invalid_argument,
                                       "PAMatrix::assemble(): GEMM dimension overflow; reduce the workspace budget");
          Impl::gemm<DeviceType,Scalar>('N', 'T', Fr, ordinal_type(Nr/Pr), Pr, Scalar(1), op.opView.data(), Fr, in.data(), Scalar(0), out.data());
          Nr = (Nr / Pr) * Fr;
        }
        auto finalOut = ((numLeftIntegrals+numTensorComponents+1)%2 == 0) ? workspace1 : workspace2;
        auto sumPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,leftOutputSpan,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): sum unit-column final output", sumPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &field, const ordinal_type &blockColumn)
        {
          const size_t idx = field + (blockColumn + size_t(cell) * blockSize) * leftOutputSpan;
          outputColumns(cell,field+leftOrdinalOffset,blockColumn) += finalOut(idx);
        });
      }

      if (outputRank == 3)
      {
        auto matrix = integrals.template getUnderlyingView<3>();
        auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,numFieldsLeft,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): copy unit-column results", copyPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &row, const ordinal_type &blockColumn)
        {
          const ordinal_type computedRow = hasLeftFilter ? leftOrdinalFilter(row) : row;
          matrix(cell + cellStart,row,columnStart + blockColumn) = outputColumns(cell,computedRow,blockColumn);
        });
      }
      else
      {
        auto matrix = integrals.template getUnderlyingView<2>();
        auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numFieldsLeft,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): copy constant-cell unit-column results", copyPolicy,
        KOKKOS_LAMBDA(const ordinal_type &row, const ordinal_type &blockColumn)
        {
          const ordinal_type computedRow = hasLeftFilter ? leftOrdinalFilter(row) : row;
          matrix(row,columnStart + blockColumn) = outputColumns(0,computedRow,blockColumn);
        });
      }
      ExecutionSpace().fence();
    }
  }
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assembleByApply(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter(), std::invalid_argument,
      "PAMatrix::assemble() supports ordinal filters only for structures supported by direct assembly (see supportsDirectAssembly())");

  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();
  const int outputRank = integrals.getUnderlyingViewRank();
  const ordinal_type numCellsComputed = (outputRank == 3) ? integrals.getDataExtent(0) : 1;

  const ordinal_type maxColumnBlockSize = 64;
  const ordinal_type columnBlockSize = std::max<ordinal_type>(1, std::min(numFieldsRight, maxColumnBlockSize));

  // per cell: input and output columns, and two contraction workspaces
  const size_t perCellEntries = size_t(columnBlockSize) * (size_t(numFieldsRight) + numFieldsLeft + 2 * maxIntermediateSize_);
  const ordinal_type worksetSize = pamatrixWorksetSize<Scalar>(options.workspaceMB, perCellEntries, numCellsComputed);

  for (ordinal_type columnStart=0; columnStart<numFieldsRight; columnStart+=columnBlockSize)
  {
    const ordinal_type blockSize = std::min(columnBlockSize, numFieldsRight - columnStart);
    for (ordinal_type cellStart=0; cellStart<numCellsComputed; cellStart+=worksetSize)
    {
      const ordinal_type Cw = std::min(worksetSize, numCellsComputed - cellStart);
      Kokkos::View<Scalar***,DeviceType> inputColumns ("PAMatrix assemble input unit columns",  Cw, numFieldsRight, blockSize);
      Kokkos::View<Scalar***,DeviceType> outputColumns("PAMatrix assemble output columns",      Cw, numFieldsLeft,  blockSize);
      auto identityPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{Cw,blockSize});
      Kokkos::parallel_for("PAMatrix::assemble(): set unit multivector", identityPolicy,
      KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &blockColumn)
      {
        inputColumns(cell,columnStart + blockColumn,blockColumn) = Scalar(1);
      });
      const size_t workspaceSize = maxIntermediateSize_ * size_t(Cw) * blockSize;
      View1D workspace(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble apply workspace"), 2 * workspaceSize);
      applyUnorientedWorkset(outputColumns, 0, inputColumns, 0, cellStart, Cw, blockSize,
                             workspace.data(), workspace.data() + workspaceSize);

      if (outputRank == 3)
      {
        auto matrix = integrals.template getUnderlyingView<3>();
        auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,numFieldsLeft,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): copy general columns", copyPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &row, const ordinal_type &blockColumn)
        {
          matrix(cell + cellStart,row,columnStart + blockColumn) = outputColumns(cell,row,blockColumn);
        });
      }
      else
      {
        auto matrix = integrals.template getUnderlyingView<2>();
        auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numFieldsLeft,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): copy constant-cell columns", copyPolicy,
        KOKKOS_LAMBDA(const ordinal_type &row, const ordinal_type &blockColumn)
        {
          matrix(row,columnStart + blockColumn) = outputColumns(0,row,blockColumn);
        });
      }
      ExecutionSpace().fence();
    }
  }
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assembleUnoriented(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const
{
  const bool direct = supportsDirectAssembly();
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter() && !direct, std::invalid_argument,
      "PAMatrix::assemble() supports ordinal filters only for structures supported by direct assembly (see supportsDirectAssembly())");
  switch (options.algorithm)
  {
    case PAMatrixAssemblyOptions::Algorithm::Automatic:
      if (direct) assembleDirect(integrals, options);
      else        assembleByApply(integrals, options);
      break;
    case PAMatrixAssemblyOptions::Algorithm::UnitColumns:
      if (direct) assembleUnitColumns(integrals, options);
      else        assembleByApply(integrals, options);
      break;
    case PAMatrixAssemblyOptions::Algorithm::ApplyToColumns:
      assembleByApply(integrals, options);
      break;
  }
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::orientAssembledMatrix(Data<Scalar,DeviceType> &integrals) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  // modifyMatrixByOrientation() computes O_L M O_R^T
  auto matrixView = integrals.getUnderlyingView();
  Kokkos::DynRankView<Scalar,DeviceType> unoriented("PAMatrix unoriented matrix", matrixView.extent(0), matrixView.extent(1), matrixView.extent(2));
  Kokkos::deep_copy(unoriented, matrixView);
  OrientationTools<DeviceType>::modifyMatrixByOrientation(matrixView, unoriented, _orientations, _basisLeft.get(), _basisRight.get());
  ExecutionSpace().fence();
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assemble(Data<Scalar,DeviceType> &integrals, const PAMatrixAssemblyOptions &options) const
{
  const ordinal_type numCells       = _basisValuesLeft.numCells();
  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();

  INTREPID2_TEST_FOR_EXCEPTION(integrals.rank() != 3, std::invalid_argument,
      "PAMatrix::assemble() requires a rank-3 logical output container");
  INTREPID2_TEST_FOR_EXCEPTION(integrals.extent_int(0) != numCells ||
                               integrals.extent_int(1) != numFieldsLeft ||
                               integrals.extent_int(2) != numFieldsRight, std::invalid_argument,
      "PAMatrix::assemble() output container has incompatible extents");
  const int outputRank = integrals.getUnderlyingViewRank();
  INTREPID2_TEST_FOR_EXCEPTION((outputRank != 2) && (outputRank != 3), std::invalid_argument,
      "PAMatrix::assemble() supports only rank-2 or rank-3 underlying matrix storage");
  INTREPID2_TEST_FOR_EXCEPTION(_hasNontrivialOrientations && ((outputRank != 3) || (integrals.getDataExtent(0) != numCells)), std::invalid_argument,
      "PAMatrix::assemble() with orientations requires storage for every cell; see allocateMatrixStorage()");

  assembleUnoriented(integrals, options);
  if (_hasNontrivialOrientations) orientAssembledMatrix(integrals);
}
template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::extractDiagonal(const Data<Scalar,DeviceType> &diagonal) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter(), std::invalid_argument, "PAMatrix::extractDiagonal() does not yet support ordinal filters");

  const ordinal_type numCells       = _basisValuesLeft.numCells();
  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();
  const ordinal_type numDiagonalEntries = std::min(numFieldsLeft, numFieldsRight);

  INTREPID2_TEST_FOR_EXCEPTION(diagonal.rank() != 2, std::invalid_argument,
      "PAMatrix::extractDiagonal() requires a rank-2 logical output container");
  INTREPID2_TEST_FOR_EXCEPTION(diagonal.extent_int(0) != numCells ||
                               diagonal.extent_int(1) != numDiagonalEntries, std::invalid_argument,
      "PAMatrix::extractDiagonal() output container has incompatible extents");

  // fallback: assemble the full matrix (with orientations), and copy its diagonal.  (Assembly uses a workspace budget, but
  // the assembled temporary holds the matrices of all cells for which the diagonal stores data.)
  const auto assembleAndCopyDiagonal = [this, numCells, numFieldsLeft, numFieldsRight, numDiagonalEntries, diagonal]()
  {
    const ordinal_type numCellsStored = (diagonal.getUnderlyingViewRank() == 2) ? diagonal.getDataExtent(0) : 1;
    // the matrix needs storage only for the cells for which the diagonal stores data (all cells, with orientations)
    Kokkos::View<Scalar***,DeviceType> matrixView("PAMatrix extractDiagonal assembled temporary", numCellsStored, numFieldsLeft, numFieldsRight);
    Kokkos::Array<int,3> matrixExtents {numCells, numFieldsLeft, numFieldsRight};
    const DataVariationType cellVariation = (numCellsStored == numCells) ? GENERAL : MODULAR;
    Kokkos::Array<DataVariationType,3> matrixVariationTypes {cellVariation, GENERAL, GENERAL};
    Data<Scalar,DeviceType> matrixData(matrixView, matrixExtents, matrixVariationTypes);
    this->assemble(matrixData);

    if (diagonal.getUnderlyingViewRank() == 2)
    {
      auto diagonalView = diagonal.template getUnderlyingView<2>();
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCellsStored,numDiagonalEntries});
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): copy assembled diagonal", policy,
      KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &field)
      {
        diagonalView(cell,field) = matrixView(cell,field,field);
      });
    }
    else if (diagonal.getUnderlyingViewRank() == 1)
    {
      auto diagonalView = diagonal.template getUnderlyingView<1>();
      auto policy = Kokkos::RangePolicy<ExecutionSpace>(0,numDiagonalEntries);
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): copy assembled constant-cell diagonal", policy,
      KOKKOS_LAMBDA(const ordinal_type &field)
      {
        diagonalView(field) = matrixView(0,field,field);
      });
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          "PAMatrix::extractDiagonal() supports only rank-1 or rank-2 underlying diagonal storage");
    }
    ExecutionSpace().fence();
  };

  bool useSumFactorizedPath = !_hasNontrivialOrientations;
  for (const auto &integralTuple : componentIntegralsToSum_)
  {
    const auto &leftOps       = std::get<0>(integralTuple);
    const PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
    const auto &rightOps      = std::get<2>(integralTuple);
    const int leftOffset      = std::get<3>(integralTuple);
    const int leftSpan        = std::get<4>(integralTuple);
    const int rightOffset     = std::get<5>(integralTuple);
    const int rightSpan       = std::get<6>(integralTuple);

    if ((leftOps.size() != rightOps.size()) || (leftOps.size() == 0) ||
        (pointDataSpec.aSpan != 1) || (pointDataSpec.bSpan != 1) ||
        (leftOffset != rightOffset) || (leftSpan != rightSpan) ||
        (leftOffset + leftSpan > numDiagonalEntries))
    {
      useSumFactorizedPath = false;
      break;
    }

    int fieldProduct = 1;
    int pointProduct = 1;
    for (size_t r=0; r<leftOps.size(); r++)
    {
      if ((leftOps[r].M != rightOps[r].N) || (leftOps[r].N != rightOps[r].M))
      {
        useSumFactorizedPath = false;
        break;
      }
      fieldProduct *= leftOps[r].M;
      pointProduct *= leftOps[r].N;
    }
    if (!useSumFactorizedPath) break;
    if ((fieldProduct != leftSpan) || (pointProduct != pointDataSpec.P))
    {
      useSumFactorizedPath = false;
      break;
    }
  }

  if (!useSumFactorizedPath)
  {
    assembleAndCopyDiagonal();
    return;
  }

  diagonal.clear();
  Data<Scalar,DeviceType> diagonalData = diagonal;
  const int diagonalCellDataExtent = diagonal.getDataExtent(0);

  for (const auto &integralTuple : componentIntegralsToSum_)
  {
    const auto &leftOps       = std::get<0>(integralTuple);
    const PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
    const auto &rightOps      = std::get<2>(integralTuple);
    const int diagonalOffset  = std::get<3>(integralTuple);
    const int numTensorComponents = static_cast<int>(leftOps.size());

    std::vector<int> fieldExtents(numTensorComponents);
    std::vector<int> pointExtents(numTensorComponents);
    int fieldProduct = 1;
    int pointProduct = 1;
    int maxScratchEntriesPerCell = 1;
    for (int r=0; r<numTensorComponents; r++)
    {
      fieldExtents[r] = leftOps[r].M;
      pointExtents[r] = leftOps[r].N;
      fieldProduct *= fieldExtents[r];
      pointProduct *= pointExtents[r];
    }
    for (int r=1; r<numTensorComponents; r++)
    {
      int prefixPointProduct = 1;
      for (int s=0; s<r; s++) prefixPointProduct *= pointExtents[s];
      int suffixFieldProduct = 1;
      for (int s=r; s<numTensorComponents; s++) suffixFieldProduct *= fieldExtents[s];
      maxScratchEntriesPerCell = std::max(maxScratchEntriesPerCell, prefixPointProduct * suffixFieldProduct);
    }

    auto pointData = _pointDataCache.find(pointDataSpec)->second;

    if (numTensorComponents == 1)
    {
      auto leftOp  = leftOps[0];
      auto rightOp = rightOps[0];
      const int numFields = fieldExtents[0];
      const int numPoints = pointExtents[0];
      const int pointDataCellExtent = pointDataSpec.C;
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{diagonalCellDataExtent,numFields});
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): 1D sum-factorized diagonal", policy,
      KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &field)
      {
        const int pointDataCell = (pointDataCellExtent == 1) ? 0 : cell;
        Scalar value = 0;
        for (int point=0; point<numPoints; point++)
        {
          const Scalar L = leftOp.opView(field + point * leftOp.M);
          const Scalar R = rightOp.opView(point + field * rightOp.M);
          const Scalar W = pointData(pointDataCell + size_t(point) * pointDataCellExtent);
          value += L * R * W;
        }
        diagonalData.getWritableEntry(cell, diagonalOffset + field) += value;
      });
      ExecutionSpace().fence();
      continue;
    }

    const size_t scratchEntries = size_t(diagonalCellDataExtent) * maxScratchEntriesPerCell;
    Kokkos::View<Scalar*,DeviceType> scratchA("PAMatrix diagonal scratch A", scratchEntries);
    Kokkos::View<Scalar*,DeviceType> scratchB("PAMatrix diagonal scratch B", scratchEntries);

    const int lastDim = numTensorComponents - 1;
    {
      auto leftOp  = leftOps[lastDim];
      auto rightOp = rightOps[lastDim];
      int prefixPointProduct = 1;
      for (int s=0; s<lastDim; s++) prefixPointProduct *= pointExtents[s];
      const int numFieldsLast = fieldExtents[lastDim];
      const int numPointsLast = pointExtents[lastDim];
      const int pointDataCellExtent = pointDataSpec.C;
      const size_t totalEntries = size_t(diagonalCellDataExtent) * prefixPointProduct * numFieldsLast;
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): initialize last diagonal contraction", Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<size_t>>(0,totalEntries),
      KOKKOS_LAMBDA(const size_t &entryOrdinal)
      {
        const int cell = int(entryOrdinal % diagonalCellDataExtent);
        size_t remainder = entryOrdinal / diagonalCellDataExtent;
        const int prefixPointOrdinal = int(remainder % prefixPointProduct);
        const int field = int(remainder / prefixPointProduct);
        const int pointDataCell = (pointDataCellExtent == 1) ? 0 : cell;
        Scalar value = 0;
        for (int point=0; point<numPointsLast; point++)
        {
          const int pointOrdinal = prefixPointOrdinal + point * prefixPointProduct;
          const Scalar L = leftOp.opView(field + point * leftOp.M);
          const Scalar R = rightOp.opView(point + field * rightOp.M);
          const Scalar W = pointData(pointDataCell + size_t(pointOrdinal) * pointDataCellExtent);
          value += L * R * W;
        }
        scratchA(entryOrdinal) = value;
      });
      ExecutionSpace().fence();
    }

    bool scratchAIsCurrent = true;
    for (int r=numTensorComponents-2; r>=1; r--)
    {
      auto inScratch  = scratchAIsCurrent ? scratchA : scratchB;
      auto outScratch = scratchAIsCurrent ? scratchB : scratchA;
      auto leftOp  = leftOps[r];
      auto rightOp = rightOps[r];

      int prefixPointProduct = 1;
      for (int s=0; s<r; s++) prefixPointProduct *= pointExtents[s];
      const int pointExtent = pointExtents[r];
      const int previousPrefixPointProduct = prefixPointProduct * pointExtent;
      int suffixFieldProduct = 1;
      for (int s=r; s<numTensorComponents; s++) suffixFieldProduct *= fieldExtents[s];
      int suffixTailProduct = 1;
      for (int s=r+1; s<numTensorComponents; s++) suffixTailProduct *= fieldExtents[s];
      const size_t totalEntries = size_t(diagonalCellDataExtent) * prefixPointProduct * suffixFieldProduct;

      Kokkos::parallel_for("PAMatrix::extractDiagonal(): diagonal tensor contraction stage", Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<size_t>>(0,totalEntries),
      KOKKOS_LAMBDA(const size_t &entryOrdinal)
      {
        const int cell = int(entryOrdinal % diagonalCellDataExtent);
        size_t remainder = entryOrdinal / diagonalCellDataExtent;
        const int prefixPointOrdinal = int(remainder % prefixPointProduct);
        const int suffixFieldOrdinal = int(remainder / prefixPointProduct);
        const int field = suffixFieldOrdinal % leftOp.M;
        const int tailFieldOrdinal = suffixFieldOrdinal / leftOp.M;
        Scalar value = 0;
        for (int point=0; point<pointExtent; point++)
        {
          const int previousPrefixPointOrdinal = prefixPointOrdinal + point * prefixPointProduct;
          const size_t inputIndex = cell + (previousPrefixPointOrdinal + size_t(tailFieldOrdinal) * previousPrefixPointProduct) * diagonalCellDataExtent;
          const Scalar L = leftOp.opView(field + point * leftOp.M);
          const Scalar R = rightOp.opView(point + field * rightOp.M);
          value += L * R * inScratch(inputIndex);
        }
        outScratch(entryOrdinal) = value;
      });
      ExecutionSpace().fence();
      scratchAIsCurrent = !scratchAIsCurrent;
    }

    auto finalScratch = scratchAIsCurrent ? scratchA : scratchB;
    {
      auto leftOp  = leftOps[0];
      auto rightOp = rightOps[0];
      const int pointExtent = pointExtents[0];
      int suffixTailProduct = 1;
      for (int s=1; s<numTensorComponents; s++) suffixTailProduct *= fieldExtents[s];
      const size_t totalEntries = size_t(diagonalCellDataExtent) * fieldProduct;
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): final diagonal accumulation", Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<size_t>>(0,totalEntries),
      KOKKOS_LAMBDA(const size_t &entryOrdinal)
      {
        const int cell = int(entryOrdinal % diagonalCellDataExtent);
        const int fieldOrdinal = int(entryOrdinal / diagonalCellDataExtent);
        const int field0 = fieldOrdinal % leftOp.M;
        const int tailFieldOrdinal = fieldOrdinal / leftOp.M;
        Scalar value = 0;
        for (int point=0; point<pointExtent; point++)
        {
          const size_t inputIndex = cell + (point + size_t(tailFieldOrdinal) * pointExtent) * diagonalCellDataExtent;
          const Scalar L = leftOp.opView(field0 + point * leftOp.M);
          const Scalar R = rightOp.opView(point + field0 * rightOp.M);
          value += L * R * finalScratch(inputIndex);
        }
        diagonalData.getWritableEntry(cell, diagonalOffset + fieldOrdinal) += value;
      });
      ExecutionSpace().fence();
    }
  }
}


template<typename DeviceType,class Scalar>
size_t PAMatrix<DeviceType,Scalar>::applyWorkspacePerCell() const
{
  size_t size = maxIntermediateSize_ * 2;
  if (_hasNontrivialOrientations)
  {
    // unoriented and oriented input values, and unoriented and oriented output values
    size += 2 * size_t(_basisValuesRight.numFields()) + 2 * size_t(_basisValuesLeft.numFields());
  }
  return size;
}

template<typename DeviceType,class Scalar>
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize) const
{
  return allocateWorkspace(worksetSize, 1);
}

template<typename DeviceType,class Scalar>
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize,
                                                                                const ordinal_type &n) const
{
  return View1D("PAMatrix workspace", applyWorkspacePerCell() * size_t(worksetSize) * size_t(n));
}

template<typename DeviceType,class Scalar>
template<typename OutputViewType, typename InputViewType>
void PAMatrix<DeviceType,Scalar>::applyUnorientedWorkset(const OutputViewType &outputVector, const ordinal_type outputCellOffset,
                                                         const InputViewType &inputVector, const ordinal_type inputCellOffset,
                                                         const ordinal_type startCell, const ordinal_type Cw, const ordinal_type N,
                                                         Scalar *workspace1Data, Scalar *workspace2Data) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using ScratchView    = Kokkos::View<Scalar*, DeviceType, Kokkos::MemoryUnmanaged>;

  const size_t workspaceSize = maxIntermediateSize_ * size_t(Cw) * N;
  ScratchView workspace1(workspace1Data, workspaceSize);
  ScratchView workspace2(workspace2Data, workspaceSize);

  const Scalar alpha = 1.0;
  const Scalar beta  = 0.0;

  // For vector-dot-vector integrals (e.g.), we need to integrate left x components against right x components, etc., and sum.
  // Each of these is one pass, and we accumulate in outputVector.
  const int numIntegrationPasses = int(componentIntegralsToSum_.size());
  for (int integrationPass=0; integrationPass<numIntegrationPasses; integrationPass++)
  {
    const auto &integral_tuple = componentIntegralsToSum_[integrationPass];

    const auto & leftIntegrals = std::get<0>(integral_tuple);
    const PointDataSpec & pointDataSpec = std::get<1>(integral_tuple);
    // right integrals: replace F2j basis coefficients with evaluations at Pj
    const auto & rightIntegrals = std::get<2>(integral_tuple);
    const int numRightIntegrals = int(rightIntegrals.size());

    const auto  leftOrdinalOffset = std::get<3>(integral_tuple);
    const auto     leftOutputSpan = std::get<4>(integral_tuple);
    const auto rightOrdinalOffset = std::get<5>(integral_tuple);
    const auto     rightInputSpan = std::get<6>(integral_tuple);

    const int numLeftIntegrals = int(leftIntegrals.size());

    // set workspace1 to the input data in an appropriate order
    // (C,F,N) arguments, where F=F_0…F_d, and the tensor ordering of these has F_0 as the fastest-moving index.
    // We want to start with input basis coefficients that are in a tensor product ordering with shape
    // (N,C,F), where the fastest-moving indices are on the left (LayoutLeft ordering).
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,rightInputSpan,N});
    Kokkos::parallel_for("PAMatrix::apply(): copy inputVector into workspace", policy,
                         KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                         {
      workspace1(n + (c + size_t(f) * Cw) * N) = inputVector.access(c + inputCellOffset,f+rightOrdinalOffset,n);
    });

    // the right integrals are ordered in the natural dimension ordering: x integrals come first.
    // At the start, input has shape (N,C,F_0,...,F_d); we regard this as a matrix (NCF_0…F_{d-1} x F_d)
    // We left-multiply the transpose of this by a right-integral matrix A_d with shape (P_d,F_d), with result
    // of shape (P_d,N,C,F_0,...,F_{d-1}).  In the next iteration, we multiply by A_{d-1} to get
    // (P_{d-1},P_d,N,C,F_0,...,F_{d-2}), continuing until we have (P_0,…,P_d,N,C) = (P,N,C).
    // We then weight with point data; as we do, we reshape to a form (N,C,P).
    // We then apply the transpose of left-integral matrix A_d to produce (F_d,N,C,P_0,…,P_{d-1}),
    // ending with (F_0,...,F_d,N,C) = (F,N,C), where now the F's belong to the left basis.
    // Finally, we accumulate the result into outputVector as (C,F,N), in Kokkos's layout for outputVector.

    // define Nr to be the size of the most recent workspace data
    size_t Nr = size_t(Cw) * N * rightInputSpan;
    for (int j=0; j<numRightIntegrals; j++)
    {
      // we alternate whether we are placing intermediate results in workspace1 or workspace2
      auto  in = (j%2 == 0) ? workspace1 : workspace2;
      auto out = (j%2 == 0) ? workspace2 : workspace1;

      const int r = numRightIntegrals-j-1; // start with final component, work back to 0th component.
      const auto &op = rightIntegrals[r];

      // multiply P_r x F_r matrix with transpose of (… x F_r)
      const ordinal_type Pr = op.M;
      const ordinal_type Fr = op.N;
      const ordinal_type LDA = Pr; // will need to revise if we ever pad our operators (for byte alignment)
      INTREPID2_TEST_FOR_EXCEPTION(Nr % Fr != 0, std::invalid_argument, "Error: Nr must be a multiple of Fr");
      INTREPID2_TEST_FOR_EXCEPTION(Nr / Fr > size_t(std::numeric_limits<ordinal_type>::max()), std::invalid_argument,
                                   "PAMatrix::apply(): GEMM dimension overflow; use a smaller workset size");
      Impl::gemm<DeviceType,Scalar>('N', 'T', Pr, ordinal_type(Nr/Fr), Fr, alpha, op.opView.data(), LDA, in.data(), beta, out.data());
      Nr = (Nr / Fr) * Pr;
    }
    auto  pointDataIn = (numRightIntegrals%2 == 0) ? workspace1 : workspace2; // pointwise result from contractions so far
    auto pointDataOut = (numRightIntegrals%2 == 0) ? workspace2 : workspace1; // pointwise output from weighting with pointData

    const auto &pointData = _pointDataCache.find(pointDataSpec)->second; // pointwise weights with shape (C,P[,Da[,Db]])
    // pointDataIn has shape (P,N,Cw); pointDataOut will have shape (N,Cw,P) (layout left).
    Impl::pointDataMultiply<DeviceType,Scalar>(pointDataSpec.C, startCell, Cw, pointDataSpec.P, pointDataSpec.aSpan, pointDataSpec.bSpan,
                                               pointData.data(), pointDataIn.data(), pointDataOut.data(), N);

    Nr = Nr / pointDataSpec.bSpan * pointDataSpec.aSpan; // contracted in b, expanded in a

    for (int i=0; i<numLeftIntegrals; i++)
    {
      // we alternate whether we are placing intermediate results in workspace1 or workspace2
      auto  in = ((i+numRightIntegrals+1)%2 == 0) ? workspace1 : workspace2;
      auto out = ((i+numRightIntegrals+1)%2 == 0) ? workspace2 : workspace1;

      const int r = numLeftIntegrals-i-1; // start with final component, work back to 0th component.
      const auto &op = leftIntegrals[r];

      // multiply Fr x Pr matrix with transpose of (… x Pr): result is (Fr x …)
      const ordinal_type Fr = op.M;
      const ordinal_type Pr = op.N;
      const ordinal_type LDA = Fr; // will need to revise if we ever pad our operators (for byte alignment)
      INTREPID2_TEST_FOR_EXCEPTION(Nr % Pr != 0, std::invalid_argument, "Error: Nr must be a multiple of Pr");
      INTREPID2_TEST_FOR_EXCEPTION(Nr / Pr > size_t(std::numeric_limits<ordinal_type>::max()), std::invalid_argument,
                                   "PAMatrix::apply(): GEMM dimension overflow; use a smaller workset size");
      Impl::gemm<DeviceType,Scalar>('N', 'T', Fr, ordinal_type(Nr/Pr), Pr, alpha, op.opView.data(), LDA, in.data(), beta, out.data());
      Nr = (Nr / Pr) * Fr;
    }
    auto finalOut = ((numLeftIntegrals+numRightIntegrals+1)%2 == 0) ? workspace1 : workspace2;
    // Sum finalOut into outputVector
    policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,leftOutputSpan,N});
    Kokkos::parallel_for("PAMatrix::apply(): sum finalOut into outputVector", policy,
                         KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                         {
      const size_t idx = f + (n + size_t(c) * N) * leftOutputSpan; // (F,N,C) layout-left
      outputVector.access(c+outputCellOffset,f+leftOrdinalOffset,n) += finalOut(idx);
    });
  } // integrationPass for loop
}

template<typename DeviceType,class Scalar>
template<typename OutputViewType, typename InputViewType>
void PAMatrix<DeviceType,Scalar>::apply(const OutputViewType &outputVector,
                                        const  InputViewType & inputVector,
                                        const Kokkos::View<Scalar*,DeviceType> &workspace,
                                        const bool sumInto,
                                        const int worksetSizeIn) const
{
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter(), std::invalid_argument, "PAMatrix::apply() does not yet support ordinal filters");

  using ExecutionSpace = typename DeviceType::execution_space;

  INTREPID2_TEST_FOR_EXCEPTION((inputVector.rank() < 2) || (inputVector.rank() > 3) || (outputVector.rank() != inputVector.rank()),
                               std::invalid_argument, "PAMatrix::apply(): input and output must both have shape (C,F) or (C,F,N)");
  const ordinal_type C  = inputVector.extent_int(0); // C, F2, N
  const ordinal_type F2 = inputVector.extent_int(1);
  const ordinal_type N  = (inputVector.rank() > 2) ? inputVector.extent_int(2) : 1;
  const ordinal_type F1 = outputVector.extent_int(1); // C, F1, N
  const ordinal_type outputN = (outputVector.rank() > 2) ? outputVector.extent_int(2) : 1;

  INTREPID2_TEST_FOR_EXCEPTION(C != _basisValuesLeft.numCells(), std::invalid_argument, "PAMatrix::apply(): input cell count does not match the matrix");
  INTREPID2_TEST_FOR_EXCEPTION(outputVector.extent_int(0) != C, std::invalid_argument, "PAMatrix::apply(): output cell count does not match the input");
  INTREPID2_TEST_FOR_EXCEPTION(F2 != _basisValuesRight.numFields(), std::invalid_argument, "PAMatrix::apply(): input field count does not match the right basis");
  INTREPID2_TEST_FOR_EXCEPTION(F1 != _basisValuesLeft.numFields(),  std::invalid_argument, "PAMatrix::apply(): output field count does not match the left basis");
  INTREPID2_TEST_FOR_EXCEPTION(outputN != N, std::invalid_argument, "PAMatrix::apply(): input and output vector counts differ");

  const ordinal_type worksetSize = ((worksetSizeIn > 0) && (worksetSizeIn < C)) ? worksetSizeIn : C;
  INTREPID2_TEST_FOR_EXCEPTION(applyWorkspacePerCell() * size_t(worksetSize) * N > workspace.extent(0), std::invalid_argument,
                               "PAMatrix::apply(): allocated workspace size is not sufficient; see allocateWorkspace()");

  if (!sumInto)
    Kokkos::deep_copy(outputVector, Scalar(0));

  const size_t workspaceSize = maxIntermediateSize_ * size_t(worksetSize) * N;
  Scalar * const workspace1 = workspace.data();
  Scalar * const workspace2 = workspace.data() + workspaceSize;

  if (!_hasNontrivialOrientations)
  {
    for (ordinal_type startCell=0; startCell<C; startCell+=worksetSize)
    {
      const ordinal_type Cw = std::min(worksetSize, C - startCell);
      applyUnorientedWorkset(outputVector, startCell, inputVector, startCell, startCell, Cw, N, workspace1, workspace2);
    }
    ExecutionSpace().fence();
    return;
  }

  // With orientations, the matrix is O_L M O_R^T: we apply O_R^T to the input, apply M, and apply O_L to the result.
  // OrientationTools treats the trailing dimension of (Cw,F,N) views as a point dimension, so that all N vectors of a
  // workset are oriented by one call.
  using DynScratchView = Kokkos::DynRankView<Scalar, DeviceType, Kokkos::MemoryUnmanaged>;
  for (ordinal_type startCell=0; startCell<C; startCell+=worksetSize)
  {
    const ordinal_type Cw = std::min(worksetSize, C - startCell);
    const std::pair<int,int> cellRange = {startCell, startCell + Cw};
    size_t offset = 2 * workspaceSize;
    DynScratchView unorientedInput (workspace.data() + offset, Cw, F2, N); offset += size_t(Cw) * F2 * N;
    DynScratchView orientedInput   (workspace.data() + offset, Cw, F2, N); offset += size_t(Cw) * F2 * N;
    DynScratchView unorientedOutput(workspace.data() + offset, Cw, F1, N); offset += size_t(Cw) * F1 * N;
    DynScratchView orientedOutput  (workspace.data() + offset, Cw, F1, N);

    auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,F2,N});
    Kokkos::parallel_for("PAMatrix::apply(): copy input workset", copyPolicy,
    KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
    {
      unorientedInput(c,f,n) = inputVector.access(c + startCell, f, n);
    });

    auto orientationsWorkset = Kokkos::subview(_orientations, cellRange);
    OrientationTools<DeviceType>::modifyBasisByOrientationTranspose(orientedInput, unorientedInput, orientationsWorkset, _basisRight.get());

    Kokkos::deep_copy(unorientedOutput, Scalar(0));
    applyUnorientedWorkset(unorientedOutput, 0, orientedInput, 0, startCell, Cw, N, workspace1, workspace2);

    OrientationTools<DeviceType>::modifyBasisByOrientation(orientedOutput, unorientedOutput, orientationsWorkset, _basisLeft.get());

    auto sumPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,F1,N});
    Kokkos::parallel_for("PAMatrix::apply(): sum oriented output into outputVector", sumPolicy,
    KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
    {
      outputVector.access(c+startCell,f,n) += orientedOutput(c,f,n);
    });
  }
  ExecutionSpace().fence();
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::setGEMMStatisticsEnabled(const bool enabled)
{
  Impl::PAMatrixGEMMStatistics::enabled() = enabled;
}

template<typename DeviceType,class Scalar>
bool PAMatrix<DeviceType,Scalar>::gemmStatisticsEnabled()
{
  return Impl::PAMatrixGEMMStatistics::enabled();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmFlopCount()
{
  return Impl::PAMatrixGEMMStatistics::cumulativeFlops();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmTimeSeconds()
{
  // rank-local (no communication)
  return Impl::PAMatrixGEMMStatistics::timer()->totalElapsedTime();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmThroughputGFlops(const double baseFlopCount, const double baseTimeSeconds)
{
  const double timeInSeconds = gemmTimeSeconds() - baseTimeSeconds;
  if (timeInSeconds > 0)
  {
    const double approximateFlopCountTotal = gemmFlopCount() - baseFlopCount;
    return approximateFlopCountTotal / timeInSeconds / 1.0e9;
  }
  else
  {
    return 0;
  }
}

} // end namespace Intrepid2
#endif
