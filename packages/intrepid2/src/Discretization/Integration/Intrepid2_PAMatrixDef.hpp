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
#include "Intrepid2_IntegrationTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include <Teuchos_BLAS.hpp>

#ifdef HAVE_INTREPID2_KOKKOSKERNELS
#include <KokkosBlas.hpp>
#endif

#include <algorithm>
#include <cstdlib>
#include <type_traits>

namespace Intrepid2 {

namespace Impl
{
//! For matrix-valued A(C,P,Da,Db) and vector-valued (P[,Db],N,C), output C(N,C,P,Da) representing the pointwise matrix-vector product.
template<typename DeviceType, class Scalar>
void pointDataMultiply(const ordinal_type numCells, const ordinal_type numPoints, const ordinal_type aSpan, const ordinal_type bSpan,
                       const Scalar* A, const Scalar *B, Scalar *C, const ordinal_type N)
{
  // with layout left convention:
  // A has shape (C,P,Da,Db).
  // B has shape (P[,Db],N,C).
  // C has shape (N,C,P[,Da]).

  // Limitation: only aSpan == 1 is supported.  Non-trivial Da would require packing the vector index into the point component
  // over which the left integral performs a dot product (and additional arguments describing the point components).
  INTREPID2_TEST_FOR_EXCEPTION(aSpan != 1, std::invalid_argument, "aSpan != 1 is not yet supported");

  using ExecutionSpace = typename DeviceType::execution_space;
  // Iterate in the output's memory order (n fastest, then cell, then point, then a), so that writes to C are contiguous;
  // MDRangePolicy's default iteration makes the last index fastest on host.
  auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{aSpan,numPoints,numCells,N});

  Kokkos::parallel_for("pointDataMultiply", policy,
                       KOKKOS_LAMBDA(const ordinal_type &a, const ordinal_type &point, const ordinal_type &cell, const ordinal_type &n)
  {
    Scalar value = 0;
    // assume layout left for A,B,C
    const int b0 = 0; // placeholder to make formulas clear
    const Scalar* entryA = A + cell  + (point + (a +   b0 * aSpan ) * numPoints) * numCells; // shape (C,P,Da,Db).
    const Scalar* entryB = B + point + (   b0 + (n + cell *     N ) * bSpan    ) * numPoints; // (P[,Db],N,C)
    const ordinal_type strideA = aSpan * numPoints * numCells;
    const ordinal_type strideB = numPoints;
    for (int b=0; b<bSpan; b++)
    {
      value += *entryA * *entryB;
      entryA += strideA;
      entryB += strideB;
    }
    Scalar* entryC = C + n + (cell + (point + a * numPoints) * numCells) * N; // (N,C,P[,Da])
    *entryC = value;
  });
  ExecutionSpace().fence();
}

  template<typename DeviceType,typename Scalar>
  std::enable_if_t<std::is_same<typename DeviceType::execution_space, typename Kokkos::Serial::execution_space>::value>
  gemm(const char transA, const char transB,
       const ordinal_type &M, const ordinal_type &N, const ordinal_type &K,
       const Scalar &alpha, const Scalar* A, const ordinal_type &LDA,
       const Scalar *B, const Scalar &beta, Scalar *C)
  {
    Teuchos::RCP<Teuchos::TimeMonitor> gemmTimer;
    if (PAMatrix<DeviceType,Scalar>::gemmStatisticsEnabled())
    {
      static const auto timer = Teuchos::TimeMonitor::getNewTimer("gemm");
      gemmTimer = Teuchos::rcp(new Teuchos::TimeMonitor(*timer));
    }
    Teuchos::ETransp trA = (transA == 'T') ? Teuchos::TRANS : (transA == 'C') ? Teuchos::CONJ_TRANS : Teuchos::NO_TRANS;
    Teuchos::ETransp trB = (transB == 'T') ? Teuchos::TRANS : (transB == 'C') ? Teuchos::CONJ_TRANS : Teuchos::NO_TRANS;
    Teuchos::BLAS<int,Scalar> blas;
    const ordinal_type LDB = (trB == Teuchos::TRANS) ? N : K;
    const ordinal_type LDC = M;
    INTREPID2_TEST_FOR_EXCEPTION(LDA==0, std::invalid_argument, "LDA cannot be 0");
    INTREPID2_TEST_FOR_EXCEPTION(LDB==0, std::invalid_argument, "LDB cannot be 0");
    INTREPID2_TEST_FOR_EXCEPTION(LDC==0, std::invalid_argument, "LDC cannot be 0");
    blas.GEMM(trA, trB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    PAMatrix<DeviceType,Scalar>::recordGEMMFlops(M,N,K);
  }

#ifdef HAVE_INTREPID2_KOKKOSKERNELS
  template<typename DeviceType,typename Scalar>
  std::enable_if_t<!std::is_same<typename DeviceType::execution_space, typename Kokkos::Serial::execution_space>::value>
  gemm(const char transA, const char transB,
       const ordinal_type &M, const ordinal_type &N, const ordinal_type &K,
       const Scalar &alpha, const Scalar* A, const ordinal_type &LDA,
       const Scalar *B, const Scalar &beta, Scalar *C)
  {
    Teuchos::RCP<Teuchos::TimeMonitor> gemmTimer;
    if (PAMatrix<DeviceType,Scalar>::gemmStatisticsEnabled())
    {
      static const auto timer = Teuchos::TimeMonitor::getNewTimer("gemm");
      gemmTimer = Teuchos::rcp(new Teuchos::TimeMonitor(*timer));
    }
    using ConstView2D = Kokkos::View<const Scalar**, Kokkos::LayoutLeft, DeviceType, Kokkos::MemoryUnmanaged>;
    using      View2D = Kokkos::View<      Scalar**, Kokkos::LayoutLeft, DeviceType, Kokkos::MemoryUnmanaged>;
    ConstView2D AView = (transA != 'N') ? ConstView2D(A,K,M) : ConstView2D(A,M,K);
    ConstView2D BView = (transB != 'N') ? ConstView2D(B,N,K) : ConstView2D(B,K,N);
    View2D CView(C,M,N);

    typename DeviceType::execution_space exec_space;
    KokkosBlas::gemm(exec_space, &transA, &transB, alpha, AView, BView, beta, CView);
    PAMatrix<DeviceType,Scalar>::recordGEMMFlops(M,N,K);
  }
#else
  template<typename DeviceType,typename Scalar>
  std::enable_if_t<!std::is_same<typename DeviceType::execution_space, typename Kokkos::Serial::execution_space>::value>
  gemm(const char transA, const char transB,
       const ordinal_type &M, const ordinal_type &N, const ordinal_type &K,
       const Scalar &alpha, const Scalar* A, const ordinal_type &LDA,
       const Scalar *B, const Scalar &beta, Scalar *C)
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "To support gemm on this ExecutionSpace, please build Intrepid2 with KokkosKernels");
  }
#endif

// Define GemmDeviceType: use Kokkos-supported GPUs if enabled; otherwise use serial.  Note that on macOS if you use serial and are using Apple's BLAS on an M-series Mac, it will run on the M-series GPU (and will be very fast).
#if defined(KOKKOS_ENABLE_CUDA)
using GemmDeviceType = Kokkos::Cuda;
#elif defined(KOKKOS_ENABLE_HIP)
using GemmDeviceType = Kokkos::HIP;
#else
using GemmDeviceType = Kokkos::Serial;
#endif

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
} // namespace Impl

template<typename DeviceType,class Scalar>
PAMatrix<DeviceType,Scalar>::PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
                                      const TensorData<Scalar,DeviceType> cellMeasures,
                                      const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
                                      const ScalarView<Orientation,DeviceType> orientations)
:
_cellMeasures(cellMeasures),
_basisValuesLeft(basisValuesLeft),
_basisValuesRight(basisValuesRight),
_orientations(orientations)
{
  init(basisValuesLeft, cellMeasures, basisValuesRight, orientations);
} // PAMatrix()

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::init(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeftInput,
                                  const TensorData<Scalar,DeviceType> cellMeasures,
                                  const TransformedBasisValues<Scalar,DeviceType> basisValuesRightInput,
                                  const ScalarView<Orientation,DeviceType> orientations) {
  using ExecutionSpace = typename DeviceType::execution_space;

  const bool layoutLeft = layoutLeft_;

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
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter() && (orientations.size() > 0), std::invalid_argument,
                               "PAMatrix does not support ordinal filters together with orientations");

  const int spaceDim = basisValuesLeft.spaceDim();

  // MARK: checks for supported construction
  INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.spaceDim() != basisValuesRight.spaceDim(), std::invalid_argument, "basisValuesLeft and basisValuesRight must agree on the space dimension");

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
    Kokkos::Array<int,7> maxFieldsForComponentLeft  {0,0,0,0,0,0,0};
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
          INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsLeft != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in basisValuesLeft must have the same number of tensor components as every other");
          for (int r=0; r<numTensorComponentsLeft; r++)
          {
            maxFieldsForComponentLeft[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentLeft[r]);
          }
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
    Kokkos::Array<int,7> maxFieldsForComponentRight {0,0,0,0,0,0,0};
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
          INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsRight != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in basisValuesRight must have the same number of tensor components as every other");
          for (int r=0; r<numTensorComponentsRight; r++)
          {
            maxFieldsForComponentRight[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentRight[r]);
          }
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

  // for now, we don't check for separability: we always treat as non-separable.
  // we could gain some performance in the not-too-common case that the integrals are separable; I haven't yet figured out how much.
  {
    _separable = false;
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

        composedTransform = Data<Scalar,DeviceType>::allocateMatMatResult(transposeLeft, leftTransformMatrix, transposeRight, rightTransformMatrix); // false: don't transpose
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
        if (leftRank == 3)
        {
          // extents[3] and variationTypes[3] will already be 1 and CONSTANT, respectively
        }
        else
        {
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
      composedTransform = rightTransform;
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

      Kokkos::View<Scalar*,DeviceType> identityUnderlyingView("Intrepid2::FST::integrate() - identity view",spaceDim);
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
      // set up the individual operators as 1D views
      // left operators contract in the point (and space) dimensions
      std::vector<OpSpec> leftOperators;
      for (int r=0; r<leftComponent.numTensorComponents(); r++)
      {
        const auto opData  = leftComponent.getTensorComponent(r).getUnderlyingView();
        const int opFields = opData.extent_int(0);
        const int opPoints = opData.extent_int(1);
        View1D opView("leftOp 1D view", opData.size());
        if (opData.rank() == 2) // (F,P)
        {
          auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{opFields,opPoints});
          Kokkos::parallel_for("pack 1D opView", policy,
          KOKKOS_LAMBDA(const int &field, const int &pt)
          {
            const int idx = layoutLeft ? field + pt * opFields : pt + field * opPoints;
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
            const int idx = layoutLeft ? field + (pt + d * opPoints) * opFields : d + (pt + field * opPoints) * opDim;
            opView(idx) = opData(field,pt,d);
          });
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "PAMatrix: Unsupported component operator rank");
        }
        leftOperators.push_back({opView,opFields,opPoints});
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

          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(leftComponent.numTensorComponents() != rightComponent.numTensorComponents(), std::invalid_argument, "left TensorData and right TensorData have different number of tensor components.  This is not supported.");

          // right operators contract in the field dimension
          std::vector<OpSpec> rightOperators;
          for (int r=0; r<rightComponent.numTensorComponents(); r++)
          {
            const auto  opData = rightComponent.getTensorComponent(r).getUnderlyingView();
            const int opFields = opData.extent_int(0);
            const int opPoints = opData.extent_int(1);

            View1D opView("rightOp 1D view", opData.size());
            if (opData.rank() == 2) // (F,P), but will pack as (P,F)
            {
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{opFields,opPoints});
              Kokkos::parallel_for("pack 1D opView", policy,
              KOKKOS_LAMBDA(const int &field, const int &pt)
              {
                const int idx = layoutLeft ? pt + field * opPoints : field + pt * opFields;
                opView(idx) = opData(field,pt);
              });
            }
            else if (opData.rank() == 3) // (F,P,D), but will pack as (P,D,F) for contraction in F
            {
              const int opDim = opData.extent_int(2);
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{opFields,opPoints,opDim});
              Kokkos::parallel_for("pack 1D opView", policy,
              KOKKOS_LAMBDA(const int &field, const int &pt, const int &d)
              {
                const int idx = layoutLeft ? d + (pt + field * opPoints) * opDim : field + (pt + d * opPoints) * opFields;
                opView(idx) = opData(field,pt,d);
              });
            }
            else
            {
              INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "PAMatrix: Unsupported component operator rank");
            }
            rightOperators.push_back({opView,opPoints,opFields});
          }

          const int aSpan =  leftComponent.extent_int(2);
          const int bSpan = rightComponent.extent_int(2);

          const int numCells      = _composedWeightedTransform.extent_int(0);
          const int numPoints     = _composedWeightedTransform.extent_int(1);

          PointDataSpec pointDataSpec{numCells, numPoints, a_offset, b_offset, aSpan, bSpan};

          if (_pointDataCache.find(pointDataSpec) == _pointDataCache.end())
          {
            const int pointDataSize = numCells * numPoints * aSpan * bSpan;
            View1D pointDataView("pointDataView", pointDataSize);

            auto composedWeightedTransform = _composedWeightedTransform;

            if (_composedWeightedTransform.rank() == 2) // (C,P): pointwise weight
            {
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numPoints});
              Kokkos::parallel_for("pack 1D pointData", policy,
                                   KOKKOS_LAMBDA(const int &cell, const int &pt)
                                   {
                const int idx = layoutLeft ? cell + pt * numCells : pt + cell * numPoints ;
                pointDataView(idx) = composedWeightedTransform(cell,pt);
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
                const int idx = layoutLeft ? cell + (pt + d * numPoints) * numCells : d + (pt + cell * numPoints) * dSpan;
                pointDataView(idx) = composedWeightedTransform(cell,pt,dOffset + d);
              });
            }
            else if (_composedWeightedTransform.rank() == 4) // (C,P,D,D): contract in b and expand in a at each point
            {
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{numCells,numPoints,aSpan,bSpan});
              Kokkos::parallel_for("pack 1D pointData", policy,
                                   KOKKOS_LAMBDA(const int &cell, const int &pt, const int &da, const int &db)
                                   {
                const int idx = layoutLeft ? cell + (pt + (da + db * aSpan) * numPoints) * numCells
                                           : db + (da + (pt + cell * numPoints) * aSpan) * bSpan ;
                pointDataView(idx) = composedWeightedTransform(cell,pt,a_offset + da,b_offset + db);
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

  // set maxIntermediateSize_: the per-cell size required for intermediate computations, which is used to size the workspaces
  const int F2 = basisValuesRight.extent_int(1); // C,F,P,…
  const int F1 = basisValuesLeft.extent_int(1);
  maxIntermediateSize_ = std::max(F1,F2);

  for (const auto &entry : componentIntegralsToSum_)
  {
    const auto & leftOps       = std::get<0>(entry);
    const auto & pointDataSpec = std::get<1>(entry);
    const auto & rightOps      = std::get<2>(entry);

    int perCellSize = F2; // num basis coefficients in the vector we multiply
    // we start the contraction on the right
    for (const auto &rightOp : rightOps)
    {
      perCellSize /= rightOp.N;
      perCellSize *= rightOp.M;
      maxIntermediateSize_ = max(perCellSize,maxIntermediateSize_);
    }

    perCellSize /= pointDataSpec.bSpan;
    perCellSize *= pointDataSpec.aSpan;
    maxIntermediateSize_ = max(perCellSize,maxIntermediateSize_);

    for (const auto &leftOp : leftOps)
    {
      perCellSize /= leftOp.N;
      perCellSize *= leftOp.M;
      maxIntermediateSize_ = max(perCellSize,maxIntermediateSize_);
    }
  }

}


template<typename DeviceType,class Scalar>
PAMatrix<DeviceType,Scalar>::PAMatrix(const TransformedBasisValues<Scalar,DeviceType> basisValues,
                                      const TensorData<Scalar,DeviceType> cellMeasures,
                                      const ScalarView<Orientation,DeviceType> orientations)
:
PAMatrix<DeviceType,Scalar>(basisValues,cellMeasures,basisValues,orientations)
{}

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
  const auto leftTransform = _basisValuesLeft.transform();

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

  if (_orientations.size() > 0)
  {
    // Orientation application is cell-dependent even when geometry and weights
    // are otherwise constant across cells.  Keep oriented assembled matrices in
    // explicit (C,F1,F2) storage so OrientationTools never receives a rank-2
    // constant-cell view for a logically rank-3 matrix.
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

  if (_orientations.size() > 0)
  {
    // Orientations can make otherwise constant cell data cell-dependent.
    // Preserve per-cell diagonal storage so fallback assembly/copy and future
    // oriented diagonal paths do not collapse distinct oriented cells.
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
bool pamatrixAssembleUsesUnitColumnPath(const PAMatrix<DeviceType,Scalar> &matrix)
{
  bool useUnitColumnPath = (matrix._orientations.size() == 0) && (matrix.componentIntegralsToSum_.size() > 0)
                        && (matrix.componentIntegralGroups_.size() == matrix.componentIntegralsToSum_.size());
  if (!useUnitColumnPath) return false;
  // the grouped unit-column path shares point buffers across passes, so all passes must have the same point count,
  // number of tensor components, and cell count
  const int numPoints           = std::get<1>(matrix.componentIntegralsToSum_[0]).P;
  const int numTensorComponents = int(std::get<2>(matrix.componentIntegralsToSum_[0]).size());
  const int numCells            = matrix._basisValuesLeft.numCells();
  for (const auto &integralTuple : matrix.componentIntegralsToSum_)
  {
    if ((std::get<1>(integralTuple).P != numPoints) || (int(std::get<2>(integralTuple).size()) != numTensorComponents) ||
        (std::get<1>(integralTuple).C != numCells))
    {
      return false;
    }
  }
  for (const auto &integralTuple : matrix.componentIntegralsToSum_)
  {
    const auto &leftOps       = std::get<0>(integralTuple);
    const typename PAMatrix<DeviceType,Scalar>::PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
    const auto &rightOps      = std::get<2>(integralTuple);
    const int leftSpan        = std::get<4>(integralTuple);
    const int rightSpan       = std::get<6>(integralTuple);

    if ((leftOps.size() != rightOps.size()) || (rightOps.size() == 0) ||
        (rightOps.size() > 8) || (pointDataSpec.aSpan != 1) || (pointDataSpec.bSpan != 1))
    {
      useUnitColumnPath = false;
      break;
    }

    int rightFieldProduct = 1;
    int pointProduct = 1;
    for (size_t r=0; r<rightOps.size(); r++)
    {
      rightFieldProduct *= rightOps[r].N;
      pointProduct      *= rightOps[r].M;
    }
    int leftFieldProduct = 1;
    int leftPointProduct = 1;
    for (size_t r=0; r<leftOps.size(); r++)
    {
      leftFieldProduct *= leftOps[r].M;
      leftPointProduct *= leftOps[r].N;
    }
    if ((rightFieldProduct != rightSpan) || (leftFieldProduct != leftSpan) ||
        (pointProduct != pointDataSpec.P) || (leftPointProduct != pointDataSpec.P))
    {
      useUnitColumnPath = false;
      break;
    }
  }
  return useUnitColumnPath;
}

//! true unless INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS is set, in which case assemble() uses the unit-column path even
//! where the direct-contraction path applies (for verification and performance comparison)
inline bool pamatrixAssembleDirectContractionEnabled()
{
  static const bool enabled = (std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_UNIT_COLUMNS") == nullptr);
  return enabled;
}

/** \brief Direct sum-factorized assembly.

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
    Cells are processed in worksets sized to a memory budget (INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB, default 512).

    Requires the structure accepted by pamatrixAssembleUsesUnitColumnPath() (scalar point weights, one-dimensional
    tensor components, equal left and right component counts, no orientations). */
template<typename DeviceType,class Scalar>
void pamatrixAssembleDirect(const PAMatrix<DeviceType,Scalar> &matrix, Data<Scalar,DeviceType> &integrals)
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using View1D    = Kokkos::View<Scalar*,DeviceType>;
  using OrdView1D = Kokkos::View<ordinal_type*,DeviceType>;
  using OpSpec    = typename PAMatrix<DeviceType,Scalar>::OpSpec;

  const ordinal_type numCells       = matrix._basisValuesLeft.numCells();
  const ordinal_type numFieldsLeft  = matrix._basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = matrix._basisValuesRight.numFields();
  const bool hasLeftFilter  = matrix.leftOrdinalFilter_ .extent_int(0) > 0;
  const bool hasRightFilter = matrix.rightOrdinalFilter_.extent_int(0) > 0;

  const int outputRank = integrals.getUnderlyingViewRank();
  INTREPID2_TEST_FOR_EXCEPTION((outputRank != 2) && (outputRank != 3), std::invalid_argument,
      "PAMatrix::assemble() supports only rank-2 or rank-3 underlying matrix storage");
  // with rank-2 storage (constant in the cell), assemble cell 0 only
  const ordinal_type numCellsComputed = (outputRank == 3) ? numCells : 1;

  const auto &passes = matrix.componentIntegralsToSum_;
  const int numPasses = int(passes.size());
  const int R = int(std::get<2>(passes[0]).size());

  // full-ordinal row/column -> output row/column (or -1 if filtered out); with family offsets, a pass's full ordinal is
  // offset + (tensor ordinal within the family)
  auto hostLeftFilter  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), matrix.leftOrdinalFilter_);
  auto hostRightFilter = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), matrix.rightOrdinalFilter_);
  const ordinal_type numFullLeft  = hasLeftFilter  ? matrix.numUnfilteredFieldsLeft_  : numFieldsLeft;
  const ordinal_type numFullRight = hasRightFilter ? matrix.numUnfilteredFieldsRight_ : numFieldsRight;
  std::vector<ordinal_type> rowForFull(numFullLeft, -1), colForFull(numFullRight, -1);
  for (ordinal_type row=0; row<numFieldsLeft;  row++) rowForFull[hasLeftFilter  ? hostLeftFilter(row)  : row] = row;
  for (ordinal_type col=0; col<numFieldsRight; col++) colForFull[hasRightFilter ? hostRightFilter(col) : col] = col;

  // group passes by (left family span/offset, right family span/offset, last-component operator pair): within a group,
  // the last-stage contraction is shared
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

  // per-group geometry: extents, and host row/column offset tables over the pair (i,j) index space of the first R-1
  // components and the last component; output row = rowForFull[leftOffset + I], I = I_prefix + i_{R-1} * strideLast
  struct GroupLayout
  {
    int leftOffset, rightOffset;
    std::vector<int> F, G, Q;      // per component
    int prefixPairs = 1;           // prod_{r<R-1} F_r G_r
    int prefixLeft = 1, prefixRight = 1;
    View1D pairOpLast;             // (Q_{R-1}, F_{R-1} G_{R-1})
    OrdView1D prefixRowOffset, prefixColOffset; // per prefix pair: I_prefix, J_prefix
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
    for (int r=0; r<R-1; r++) { L.prefixPairs *= L.F[r] * L.G[r]; L.prefixLeft *= L.F[r]; L.prefixRight *= L.G[r]; }
    OrdView1D rowOffsets("PAMatrix assemble prefix row offsets", L.prefixPairs), colOffsets("PAMatrix assemble prefix col offsets", L.prefixPairs);
    auto rowHost = Kokkos::create_mirror_view(rowOffsets);
    auto colHost = Kokkos::create_mirror_view(colOffsets);
    for (int pair=0; pair<L.prefixPairs; pair++)
    {
      int remainder = pair, I = 0, J = 0, leftStride = 1, rightStride = 1;
      for (int r=0; r<R-1; r++)
      {
        const int ij = remainder % (L.F[r] * L.G[r]);
        remainder /= (L.F[r] * L.G[r]);
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
    // last-component pair operator
    {
      const int F = L.F[R-1], Q = L.Q[R-1], Gn = L.G[R-1];
      View1D P(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble pair operator"), size_t(Q) * F * Gn);
      const Scalar * const Lop = std::get<0>(t)[R-1].opView.data();
      const Scalar * const Rop = std::get<2>(t)[R-1].opView.data();
      Kokkos::parallel_for("PAMatrix::assemble(): last pair operator", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Q,F,Gn}),
      KOKKOS_LAMBDA(const ordinal_type &q, const ordinal_type &i, const ordinal_type &j) {
        P(q + (i + j * F) * size_t(Q)) = Lop[i + q * F] * Rop[q + j * Q];
      });
      L.pairOpLast = P;
    }
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

  // Sparse-filter path: when few (row, column) pairs survive, compute each kept entry directly at the last stage.
  // Kept pairs, per group: (prefix pair, last-component (i,j)) -> output (row, col)
  // (INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST=1 forces it whenever there is a filter, =0 disables it; for testing)
  bool useSparseLast = false;
  {
    const double keptFraction = double(numFieldsLeft) * numFieldsRight / (double(numFullLeft) * numFullRight);
    useSparseLast = (hasLeftFilter || hasRightFilter) && (keptFraction < 0.25);
    const char *sparseSetting = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_SPARSE_LAST");
    if (sparseSetting != nullptr) useSparseLast = (hasLeftFilter || hasRightFilter) && (std::atoi(sparseSetting) != 0);
  }

  // workset size from the memory budget: the largest per-cell intermediate is max over stages of (current tensor size)
  size_t perCellDoubles = 0;
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
    perCellDoubles = std::max(perCellDoubles, maxSize);
  }
  const char *budgetString = std::getenv("INTREPID2_PAMATRIX_ASSEMBLE_WORKSPACE_MB");
  const double budgetBytes = ((budgetString != nullptr) ? std::atof(budgetString) : 512.0) * 1024.0 * 1024.0;
  // three buffers of perCellDoubles (ping, pong, group accumulator), plus the last-stage accumulator
  const ordinal_type worksetSize = std::max<ordinal_type>(1, std::min<ordinal_type>(numCellsComputed,
                                   ordinal_type(budgetBytes / (4.0 * sizeof(Scalar) * double(perCellDoubles)))));

  View1D bufferA(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct A"), perCellDoubles * worksetSize);
  View1D bufferB(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct B"), perCellDoubles * worksetSize);
  View1D groupSum(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct group sum"), perCellDoubles * worksetSize);
  View1D lastAccumulator;
  if (!useSparseLast) lastAccumulator = View1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble direct last"), perCellDoubles * worksetSize);

  if (outputRank == 3) Kokkos::deep_copy(integrals.template getUnderlyingView<3>(), Scalar(0));
  else                 Kokkos::deep_copy(integrals.template getUnderlyingView<2>(), Scalar(0));

  // pair operators for the first R-1 components, per pass
  std::vector<std::vector<View1D>> pairOps(numPasses);
  for (int pass=0; pass<numPasses; pass++)
  {
    const auto &t = passes[pass];
    for (int r=0; r<R-1; r++)
    {
      const int F = std::get<0>(t)[r].M, Q = std::get<0>(t)[r].N, Gn = std::get<2>(t)[r].N;
      View1D P(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PAMatrix assemble pair operator"), size_t(Q) * F * Gn);
      const Scalar * const Lop = std::get<0>(t)[r].opView.data();
      const Scalar * const Rop = std::get<2>(t)[r].opView.data();
      Kokkos::parallel_for("PAMatrix::assemble(): pair operator", Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Q,F,Gn}),
      KOKKOS_LAMBDA(const ordinal_type &q, const ordinal_type &i, const ordinal_type &j) {
        P(q + (i + j * F) * size_t(Q)) = Lop[i + q * F] * Rop[q + j * Q];
      });
      pairOps[pass].push_back(P);
    }
  }

  for (ordinal_type cellStart=0; cellStart<numCellsComputed; cellStart+=worksetSize)
  {
    const ordinal_type Cw = std::min(worksetSize, numCellsComputed - cellStart);

    // for the dense last stage, the accumulator holds (Q-contracted) results for all groups in turn; we write each group's
    // result to the output after its last GEMM, so the accumulator is per group
    for (int g=0; g<int(groups.size()); g++)
    {
      const auto &L = layouts[g];
      const int QLast = L.Q[R-1];
      const size_t lastInputSize = size_t(QLast) * Cw * L.prefixPairs; // (p_{R-1}, c, prefix pairs)

      // sum the group's passes' stage-(R-1) inputs into groupSum
      bool firstPass = true;
      for (int pass : groups[g].passes)
      {
        const auto &t = passes[pass];
        const auto &pointDataSpec = std::get<1>(t);
        const int numPointsTotal = pointDataSpec.P;
        const int C = pointDataSpec.C;
        const Scalar * const pointData = matrix._pointDataCache.find(pointDataSpec)->second.data();
        // stage-0 input (p_0, ..., p_{R-1}, c) into bufferA
        {
          auto dest = bufferA;
          Kokkos::parallel_for("PAMatrix::assemble(): transpose point weights",
                               Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{Cw,numPointsTotal}),
          KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &point) {
            const ordinal_type dataCell = (C == 1) ? 0 : cell + cellStart;
            dest(point + size_t(cell) * numPointsTotal) = pointData[dataCell + point * C];
          });
        }
        View1D current = bufferA, other = bufferB;
        size_t currentSize = size_t(numPointsTotal) * Cw;
        for (int r=0; r<R-1; r++)
        {
          const int Q = L.Q[r], FG = L.F[r] * L.G[r];
          const size_t X = currentSize / Q;
          const bool lastPrefixStage = (r == R-2);
          // the group's passes share everything after this stage, so sum at the final prefix stage (beta = 1)
          View1D target = lastPrefixStage ? groupSum : other;
          const Scalar beta = (lastPrefixStage && !firstPass) ? Scalar(1) : Scalar(0);
          Impl::gemm<Impl::GemmDeviceType>('T', 'N', ordinal_type(X), FG, Q, Scalar(1.0), current.data(), Q,
                                           pairOps[pass][r].data(), beta, target.data());
          currentSize = X * FG;
          if (!lastPrefixStage) std::swap(current, other);
        }
        if (R == 1)
        {
          // no prefix stages: the group sum is the weights themselves
          auto src = bufferA; auto dst = groupSum;
          const bool accumulate = !firstPass;
          Kokkos::parallel_for("PAMatrix::assemble(): sum weights", Kokkos::RangePolicy<ExecutionSpace>(0, ordinal_type(currentSize)),
          KOKKOS_LAMBDA(const ordinal_type &k) { dst(k) = accumulate ? dst(k) + src(k) : src(k); });
        }
        firstPass = false;
      }
      // groupSum: (p_{R-1}, c, prefix pairs), size lastInputSize

      const int FLast = L.F[R-1], GLast = L.G[R-1];
      const int leftStrideLast = L.prefixLeft, rightStrideLast = L.prefixRight;
      const int leftOffset = L.leftOffset, rightOffset = L.rightOffset;
      const int prefixPairs = L.prefixPairs;
      auto prefixRow = L.prefixRowOffset;
      auto prefixCol = L.prefixColOffset;

      // device copies of the full->output maps, offset to this group's family
      OrdView1D rowMap("PAMatrix assemble row map", size_t(L.prefixLeft) * FLast);
      OrdView1D colMap("PAMatrix assemble col map", size_t(L.prefixRight) * GLast);
      {
        auto rowHost = Kokkos::create_mirror_view(rowMap);
        auto colHost = Kokkos::create_mirror_view(colMap);
        for (int I=0; I<L.prefixLeft*FLast;  I++) rowHost(I) = rowForFull[leftOffset  + I];
        for (int J=0; J<L.prefixRight*GLast; J++) colHost(J) = colForFull[rightOffset + J];
        Kokkos::deep_copy(rowMap, rowHost);
        Kokkos::deep_copy(colMap, colHost);
      }

      if (!useSparseLast)
      {
        // last stage: (p_{R-1}, [c, prefix pairs]) -> ([c, prefix pairs], (i_{R-1} j_{R-1}))
        const size_t X = lastInputSize / QLast;
        bool firstInFamily = true, lastInFamily = true;
        for (int h=0; h<g; h++) if (writeFamily[h] == writeFamily[g]) firstInFamily = false;
        for (int h=g+1; h<int(groups.size()); h++) if (writeFamily[h] == writeFamily[g]) lastInFamily = false;
        Impl::gemm<Impl::GemmDeviceType>('T', 'N', ordinal_type(X), FLast * GLast, QLast, Scalar(1.0), groupSum.data(), QLast,
                                         L.pairOpLast.data(), firstInFamily ? Scalar(0.0) : Scalar(1.0), lastAccumulator.data());
        if (!lastInFamily) continue;
        auto source = lastAccumulator;
        auto write = [&](auto matrixView, auto rankTag) {
          constexpr int rank = decltype(rankTag)::value;
          // cell outermost, so that writes stay within one cell's matrix; for each (cell, last-component pair), the prefix
          // pairs are contiguous in neither source nor destination, but the destination row is fixed by (iLast, pair) and the
          // writes for consecutive prefix pairs fall in a few rows (Kokkos MDRange iterates the last index fastest)
          Kokkos::parallel_for("PAMatrix::assemble(): write direct contraction",
                               Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw, FLast*GLast, prefixPairs}),
          KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &ijLast, const ordinal_type &pair) {
            const int iLast = ijLast % FLast, jLast = ijLast / FLast;
            const ordinal_type row = rowMap(prefixRow(pair) + iLast * leftStrideLast);
            const ordinal_type col = colMap(prefixCol(pair) + jLast * rightStrideLast);
            if ((row < 0) || (col < 0)) return;
            const Scalar value = source(cell + (pair + size_t(ijLast) * prefixPairs) * Cw);
            if constexpr (rank == 3) matrixView(cell + cellStart, row, col) += value;
            else                     matrixView(row, col) += value;
          });
        };
        if (outputRank == 3) write(integrals.template getUnderlyingView<3>(), std::integral_constant<int,3>());
        else                 write(integrals.template getUnderlyingView<2>(), std::integral_constant<int,2>());
      }
      else
      {
        // sparse last stage: for each kept (row, col) of this group, value(c) = sum_q groupSum(q, c, prefixPair) P(q, ij)
        std::vector<int> keptPair, keptIJ, keptRow, keptCol;
        auto rowHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowMap);
        auto colHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colMap);
        auto prefixRowHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), prefixRow);
        auto prefixColHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), prefixCol);
        for (int pair=0; pair<prefixPairs; pair++)
        {
          for (int jLast=0; jLast<GLast; jLast++)
          {
            const ordinal_type col = colHost(prefixColHost(pair) + jLast * rightStrideLast);
            if (col < 0) continue;
            for (int iLast=0; iLast<FLast; iLast++)
            {
              const ordinal_type row = rowHost(prefixRowHost(pair) + iLast * leftStrideLast);
              if (row < 0) continue;
              keptPair.push_back(pair); keptIJ.push_back(iLast + jLast * FLast); keptRow.push_back(row); keptCol.push_back(col);
            }
          }
        }
        const int numKept = int(keptPair.size());
        if (numKept == 0) continue;
        OrdView1D keptPairV("kept pair", numKept), keptIJV("kept ij", numKept), keptRowV("kept row", numKept), keptColV("kept col", numKept);
        {
          auto a = Kokkos::create_mirror_view(keptPairV); auto b = Kokkos::create_mirror_view(keptIJV);
          auto c = Kokkos::create_mirror_view(keptRowV);  auto d = Kokkos::create_mirror_view(keptColV);
          for (int k=0; k<numKept; k++) { a(k) = keptPair[k]; b(k) = keptIJ[k]; c(k) = keptRow[k]; d(k) = keptCol[k]; }
          Kokkos::deep_copy(keptPairV, a); Kokkos::deep_copy(keptIJV, b); Kokkos::deep_copy(keptRowV, c); Kokkos::deep_copy(keptColV, d);
        }
        auto source = groupSum;
        auto P = L.pairOpLast;
        auto write = [&](auto matrixView, auto rankTag) {
          constexpr int rank = decltype(rankTag)::value;
          Kokkos::parallel_for("PAMatrix::assemble(): sparse last stage",
                               Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numKept, Cw}),
          KOKKOS_LAMBDA(const ordinal_type &k, const ordinal_type &cell) {
            const Scalar * const s = &source(size_t(QLast) * (cell + size_t(keptPairV(k)) * Cw));
            const Scalar * const p = &P(size_t(QLast) * keptIJV(k));
            Scalar value = 0;
            for (int q=0; q<QLast; q++) value += s[q] * p[q];
            if constexpr (rank == 3) matrixView(cell + cellStart, keptRowV(k), keptColV(k)) += value;
            else                     matrixView(keptRowV(k), keptColV(k)) += value;
          });
        };
        if (outputRank == 3) write(integrals.template getUnderlyingView<3>(), std::integral_constant<int,3>());
        else                 write(integrals.template getUnderlyingView<2>(), std::integral_constant<int,2>());
      }
    }
  }
  ExecutionSpace().fence();
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assemble(Data<Scalar,DeviceType> &integrals) const
{
  using ExecutionSpace = typename DeviceType::execution_space;
  using MultiVector = ScalarView<Scalar,DeviceType>;
  using ScratchView = Kokkos::View<Scalar*,DeviceType,Kokkos::MemoryUnmanaged>;

  const ordinal_type numCells       = _basisValuesLeft.numCells();
  const ordinal_type numFieldsLeft  = _basisValuesLeft.numFields();
  const ordinal_type numFieldsRight = _basisValuesRight.numFields();

  INTREPID2_TEST_FOR_EXCEPTION(integrals.rank() != 3, std::invalid_argument,
      "PAMatrix::assemble() requires a rank-3 logical output container");
  INTREPID2_TEST_FOR_EXCEPTION(integrals.extent_int(0) != numCells ||
                               integrals.extent_int(1) != numFieldsLeft ||
                               integrals.extent_int(2) != numFieldsRight, std::invalid_argument,
      "PAMatrix::assemble() output container has incompatible extents");

  const bool useUnitColumnPath = pamatrixAssembleUsesUnitColumnPath(*this);
  if (useUnitColumnPath && pamatrixAssembleDirectContractionEnabled())
  {
    pamatrixAssembleDirect(*this, integrals);
    return;
  }
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter() && !useUnitColumnPath, std::invalid_argument,
      "PAMatrix::assemble() supports ordinal filters only for structures handled by its unit-column path");

  // with ordinal filters, block column j is full right ordinal rightOrdinalFilter_(columnStart+j), and the left
  // contraction produces all numUnfilteredFieldsLeft_ rows, of which rows leftOrdinalFilter_(i) are kept
  const bool hasLeftFilter  = leftOrdinalFilter_ .extent_int(0) > 0;
  const bool hasRightFilter = rightOrdinalFilter_.extent_int(0) > 0;
  const auto leftOrdinalFilter  = leftOrdinalFilter_;
  const auto rightOrdinalFilter = rightOrdinalFilter_;
  const ordinal_type numLeftRowsComputed = hasLeftFilter ? numUnfilteredFieldsLeft_ : numFieldsLeft;

  const ordinal_type maxColumnBlockSize = 64;
  const ordinal_type columnBlockSize = std::min(numFieldsRight, maxColumnBlockSize);

  for (ordinal_type columnStart=0; columnStart<numFieldsRight; columnStart+=columnBlockSize)
  {
    const ordinal_type blockSize = std::min(columnBlockSize, numFieldsRight - columnStart);
    MultiVector outputColumns("PAMatrix assemble output columns", numCells, numLeftRowsComputed, blockSize);

    auto workspace = this->allocateWorkspace(numCells, blockSize);

    if (useUnitColumnPath)
    {
      Kokkos::deep_copy(outputColumns, Scalar(0));

      // Passes are grouped by their left and right (family, component).  Passes with the same right group share their
      // right operators, so the unit columns' point values are computed once per right group; passes with the same left
      // group share their left operators, so their weighted point values are summed before a single left contraction.
      // This replaces (#left groups) x (#right groups) fill/contraction sequences with #right fills and #left
      // contractions (e.g., 4 + 4 instead of 16 for a full D x D gradient weight in 4D).
      const int numIntegrationPasses = int(componentIntegralsToSum_.size());
      const int numTensorComponents  = int(std::get<2>(componentIntegralsToSum_[0]).size());
      const int numPoints            = std::get<1>(componentIntegralsToSum_[0]).P;
      int numLeftGroups = 0, numRightGroups = 0;
      for (const auto &groups : componentIntegralGroups_)
      {
        numLeftGroups  = std::max(numLeftGroups,  groups.first  + 1);
        numRightGroups = std::max(numRightGroups, groups.second + 1);
      }

      // point values of this block's unit columns, one buffer per right group, each with layout (P,N,C) as consumed by
      // pointDataMultiply()
      const int pointBufferSize = numPoints * blockSize * numCells;
      Kokkos::View<Scalar*,DeviceType> rightFills("PAMatrix assemble unit-column point values", pointBufferSize * numRightGroups);
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

        Kokkos::Array<typename PAMatrix<DeviceType,Scalar>::OpSpec,8> rightOpsArray;
        Kokkos::Array<int,8> rightFieldExtents;
        Kokkos::Array<int,8> pointExtents;
        for (int r=0; r<numTensorComponents; r++)
        {
          rightOpsArray[r] = rightOps[r];
          rightFieldExtents[r] = rightOps[r].N;
          pointExtents[r]      = rightOps[r].M;
        }
        Scalar * const fill = rightFills.data() + rightGroup * pointBufferSize;

        // Point values of the unit trial columns: for column (f_0,...,f_{R-1}) and point (p_0,...,p_{R-1}), the value is
        // prod_r op_r(p_r, f_r).  Component 0 varies fastest in both column and point ordinals, so each thread handles one
        // (cell, column, p_1..p_{R-1}) and writes the contiguous run over p_0, after computing the product of the other
        // components once.  (The 1D operators do not depend on the cell.)
        const int numPoints0     = pointExtents[0];
        const int numPointsOuter = numPoints / numPoints0;
        auto unitColumnPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,blockSize,numPointsOuter});
        Kokkos::parallel_for("PAMatrix::assemble(): unit trial columns to points", unitColumnPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &blockColumn, const ordinal_type &pointOuter)
        {
          const int globalColumn = hasRightFilter ? rightOrdinalFilter(columnStart + blockColumn) : columnStart + blockColumn;
          Scalar * const dest = fill + numPoints0 * pointOuter + (blockColumn + cell * blockSize) * numPoints;
          if ((globalColumn < rightOrdinalOffset) || (globalColumn >= rightOrdinalOffset + rightInputSpan))
          {
            for (int p0=0; p0<numPoints0; p0++) dest[p0] = Scalar(0);
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
            outerProduct *= op.opView(point + field * op.M);
          }
          const auto op0 = rightOpsArray[0];
          const Scalar * const op0Column = op0.opView.data() + field0 * op0.M;
          for (int p0=0; p0<numPoints0; p0++) dest[p0] = outerProduct * op0Column[p0];
        });
      }
      ExecutionSpace().fence();

      const int workspace1Size = maxIntermediateSize_ * numCells * blockSize;
      const int workspace2Size = maxIntermediateSize_ * numCells * blockSize;
      ScratchView workspace1(workspace.data(), workspace1Size);
      ScratchView workspace2(workspace.data() + workspace1Size, workspace2Size);
      // the left contraction starts from pointDataOut; keep the parity used by apply() so the GEMM chain alternates the same way
      ScratchView pointDataOut = (numTensorComponents % 2 == 0) ? workspace2 : workspace1;

      for (int leftGroup=0; leftGroup<numLeftGroups; leftGroup++)
      {
        // sum, over this left group's passes, of (point weights) * (unit-column point values), in layout (N,C,P)
        int representativePass = -1;
        Kokkos::deep_copy(Kokkos::subview(pointDataOut, Kokkos::make_pair(0, pointBufferSize)), Scalar(0));
        for (int integrationPass=0; integrationPass<numIntegrationPasses; integrationPass++)
        {
          if (componentIntegralGroups_[integrationPass].first != leftGroup) continue;
          representativePass = integrationPass;
          const PointDataSpec &pointDataSpec = std::get<1>(componentIntegralsToSum_[integrationPass]);
          const int rightGroup = componentIntegralGroups_[integrationPass].second;
          const Scalar * const fill = rightFills.data() + rightGroup * pointBufferSize;
          const Scalar * const weights = _pointDataCache.find(pointDataSpec)->second.data();
          Scalar * const accumulated = pointDataOut.data();
          const int C = pointDataSpec.C;
          const int P = pointDataSpec.P;
          const int N = blockSize;
          // accumulated(n,cell,point) += weights(cell,point) * fill(point,n,cell); n fastest to match the output layout
          auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{P,C,N});
          Kokkos::parallel_for("PAMatrix::assemble(): accumulate weighted unit-column point values", policy,
          KOKKOS_LAMBDA(const ordinal_type &point, const ordinal_type &cell, const ordinal_type &n)
          {
            accumulated[n + (cell + point * C) * N] += weights[cell + point * C] * fill[point + (n + cell * N) * P];
          });
          ExecutionSpace().fence();
        }
        if (representativePass < 0) continue;

        const auto &integralTuple = componentIntegralsToSum_[representativePass];
        const auto &leftOps       = std::get<0>(integralTuple);
        const PointDataSpec &pointDataSpec = std::get<1>(integralTuple);
        const int leftOrdinalOffset  = std::get<3>(integralTuple);
        const int leftOutputSpan     = std::get<4>(integralTuple);

        Kokkos::Array<typename PAMatrix<DeviceType,Scalar>::OpSpec,8> leftOpsArray;
        for (int r=0; r<int(leftOps.size()); r++) leftOpsArray[r] = leftOps[r];

        int Nr = pointDataSpec.P * numCells * blockSize * pointDataSpec.aSpan;
        const int numLeftIntegrals = int(leftOps.size());
        for (int i=0; i<numLeftIntegrals; i++)
        {
          auto in  = ((i+numTensorComponents+1)%2 == 0) ? workspace1 : workspace2;
          auto out = ((i+numTensorComponents+1)%2 == 0) ? workspace2 : workspace1;

          const int r = numLeftIntegrals-i-1;
          auto op = leftOpsArray[r];
          const ordinal_type & Fr = op.M;
          const ordinal_type & Pr = op.N;

          const auto A = op.opView.data();
          const ordinal_type LDA = Fr;
          const auto B = in.data();
          auto C = out.data();
          INTREPID2_TEST_FOR_EXCEPTION(Nr % Pr != 0, std::invalid_argument, "Error: Nr must be a multiple of Pr");
          Impl::gemm<Impl::GemmDeviceType>('N', 'T', Fr, Nr/Pr, Pr, Scalar(1.0), A, LDA, B, Scalar(0.0), C);
          Nr = (Nr * Fr) / Pr;
        }
        auto finalOut = ((numLeftIntegrals+numTensorComponents+1)%2 == 0) ? workspace1 : workspace2;
        auto sumPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,leftOutputSpan,blockSize});
        Kokkos::parallel_for("PAMatrix::assemble(): sum unit-column final output", sumPolicy,
        KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &field, const ordinal_type &blockColumn)
        {
          const ordinal_type idx = field + (blockColumn + cell * blockSize) * leftOutputSpan;
          outputColumns(cell,field+leftOrdinalOffset,blockColumn) += finalOut(idx);
        });
        ExecutionSpace().fence();
      }
    }
    else
    {
      MultiVector inputColumns ("PAMatrix assemble input unit columns", numCells, numFieldsRight, blockSize);
      auto identityPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numFieldsRight,blockSize});
      Kokkos::parallel_for("PAMatrix::assemble(): set unit multivector", identityPolicy,
      KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &field, const ordinal_type &blockColumn)
      {
        inputColumns(cell,field,blockColumn) = (field == columnStart + blockColumn) ? Scalar(1) : Scalar(0);
      });
      this->apply(outputColumns, inputColumns, workspace, false, numCells);
    }

    if (integrals.getUnderlyingViewRank() == 3)
    {
      auto matrix = integrals.template getUnderlyingView<3>();
      auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numCells,numFieldsLeft,blockSize});
      Kokkos::parallel_for("PAMatrix::assemble(): copy general columns", copyPolicy,
      KOKKOS_LAMBDA(const ordinal_type &cell, const ordinal_type &row, const ordinal_type &blockColumn)
      {
        const ordinal_type computedRow = hasLeftFilter ? leftOrdinalFilter(row) : row;
        matrix(cell,row,columnStart + blockColumn) = outputColumns(cell,computedRow,blockColumn);
      });
    }
    else if (integrals.getUnderlyingViewRank() == 2)
    {
      auto matrix = integrals.template getUnderlyingView<2>();
      auto copyPolicy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numFieldsLeft,blockSize});
      Kokkos::parallel_for("PAMatrix::assemble(): copy constant-cell columns", copyPolicy,
      KOKKOS_LAMBDA(const ordinal_type &row, const ordinal_type &blockColumn)
      {
        const ordinal_type computedRow = hasLeftFilter ? leftOrdinalFilter(row) : row;
        matrix(row,columnStart + blockColumn) = outputColumns(0,computedRow,blockColumn);
      });
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          "PAMatrix::assemble() supports only rank-2 or rank-3 underlying matrix storage");
    }
  }
  ExecutionSpace().fence();
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

  const auto assembleAndCopyDiagonal = [this, numCells, numFieldsLeft, numFieldsRight, numDiagonalEntries, diagonal]()
  {
    Kokkos::View<Scalar***,DeviceType> matrixView("PAMatrix extractDiagonal assembled temporary",
                                                  numCells, numFieldsLeft, numFieldsRight);
    Kokkos::Array<int,3> matrixExtents {numCells, numFieldsLeft, numFieldsRight};
    Kokkos::Array<DataVariationType,3> matrixVariationTypes {GENERAL, GENERAL, GENERAL};
    Data<Scalar,DeviceType> matrixData(matrixView, matrixExtents, matrixVariationTypes);
    this->assemble(matrixData);

    if (diagonal.getUnderlyingViewRank() == 2)
    {
      auto diagonalView = diagonal.template getUnderlyingView<2>();
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numDiagonalEntries});
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

  bool useSumFactorizedPath = (_orientations.size() == 0);
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
          const Scalar W = pointData(pointDataCell + point * pointDataCellExtent);
          value += L * R * W;
        }
        diagonalData.getWritableEntry(cell, diagonalOffset + field) += value;
      });
      ExecutionSpace().fence();
      continue;
    }

    const int scratchEntries = diagonalCellDataExtent * maxScratchEntriesPerCell;
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
      const int totalEntries = diagonalCellDataExtent * prefixPointProduct * numFieldsLast;
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): initialize last diagonal contraction", Kokkos::RangePolicy<ExecutionSpace>(0,totalEntries),
      KOKKOS_LAMBDA(const ordinal_type &entryOrdinal)
      {
        const int cell = entryOrdinal % diagonalCellDataExtent;
        int remainder = entryOrdinal / diagonalCellDataExtent;
        const int prefixPointOrdinal = remainder % prefixPointProduct;
        const int field = remainder / prefixPointProduct;
        const int pointDataCell = (pointDataCellExtent == 1) ? 0 : cell;
        Scalar value = 0;
        for (int point=0; point<numPointsLast; point++)
        {
          const int pointOrdinal = prefixPointOrdinal + point * prefixPointProduct;
          const Scalar L = leftOp.opView(field + point * leftOp.M);
          const Scalar R = rightOp.opView(point + field * rightOp.M);
          const Scalar W = pointData(pointDataCell + pointOrdinal * pointDataCellExtent);
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
      const int totalEntries = diagonalCellDataExtent * prefixPointProduct * suffixFieldProduct;

      Kokkos::parallel_for("PAMatrix::extractDiagonal(): diagonal tensor contraction stage", Kokkos::RangePolicy<ExecutionSpace>(0,totalEntries),
      KOKKOS_LAMBDA(const ordinal_type &entryOrdinal)
      {
        const int cell = entryOrdinal % diagonalCellDataExtent;
        int remainder = entryOrdinal / diagonalCellDataExtent;
        const int prefixPointOrdinal = remainder % prefixPointProduct;
        const int suffixFieldOrdinal = remainder / prefixPointProduct;
        const int field = suffixFieldOrdinal % leftOp.M;
        const int tailFieldOrdinal = suffixFieldOrdinal / leftOp.M;
        Scalar value = 0;
        for (int point=0; point<pointExtent; point++)
        {
          const int previousPrefixPointOrdinal = prefixPointOrdinal + point * prefixPointProduct;
          const int inputIndex = cell + (previousPrefixPointOrdinal + tailFieldOrdinal * previousPrefixPointProduct) * diagonalCellDataExtent;
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
      const int totalEntries = diagonalCellDataExtent * fieldProduct;
      Kokkos::parallel_for("PAMatrix::extractDiagonal(): final diagonal accumulation", Kokkos::RangePolicy<ExecutionSpace>(0,totalEntries),
      KOKKOS_LAMBDA(const ordinal_type &entryOrdinal)
      {
        const int cell = entryOrdinal % diagonalCellDataExtent;
        const int fieldOrdinal = entryOrdinal / diagonalCellDataExtent;
        const int field0 = fieldOrdinal % leftOp.M;
        const int tailFieldOrdinal = fieldOrdinal / leftOp.M;
        Scalar value = 0;
        for (int point=0; point<pointExtent; point++)
        {
          const int inputIndex = cell + (point + tailFieldOrdinal * pointExtent) * diagonalCellDataExtent;
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
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize) const
{
  int size1D = maxIntermediateSize_ * worksetSize * 2;

  if (_orientations.size() > 0)
  {
    // then we need to apply orientations on the way in and on the way out, and we need workspace to do that
    size1D +=  _basisValuesLeft.numFields() * worksetSize * 2; // to support sumInto = true, we need intermediate storage for oriented left values (hence the factor of 2).
    size1D += _basisValuesRight.numFields() * worksetSize;
  }

  return View1D("PAMatrix workspace", size1D);
}

template<typename DeviceType,class Scalar>
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize,
                                                                                const ordinal_type &n) const
{
  int size1D = maxIntermediateSize_ * worksetSize * n * 2;

  if (_orientations.size() > 0)
  {
    // then we need to apply orientations on the way in and on the way out, and we need workspace to do that
    size1D +=  _basisValuesLeft.numFields() * worksetSize * n * 2; // to support sumInto = true, we need intermediate storage for oriented left values (hence the factor of 2).
    size1D += _basisValuesRight.numFields() * worksetSize * n;
  }


  return View1D("PAMatrix workspace", size1D);
}

template<typename DeviceType,class Scalar>
template<typename OutputViewType, typename InputViewType>
void PAMatrix<DeviceType,Scalar>::apply(const OutputViewType &outputVector,
                                        const  InputViewType & inputVector,
                                        const Kokkos::View<Scalar*,DeviceType> &workspace,
                                        const bool sumInto,
                                        const int worksetSizeIn) const
{
  // (Note: callers allocate the workspace for a given workset size; see allocateWorkspace().)
  INTREPID2_TEST_FOR_EXCEPTION(hasOrdinalFilter(), std::invalid_argument, "PAMatrix::apply() does not yet support ordinal filters");

  using ExecutionSpace = typename DeviceType::execution_space;

  const ordinal_type C  = inputVector.extent_int(0); // C, F2, N
  const ordinal_type F2 = inputVector.extent_int(1);
  const ordinal_type N  = inputVector.extent_int(2);
  const ordinal_type F1 = outputVector.extent_int(1); // C, F1, N

  const int worksetSize = (worksetSizeIn > 0) ? worksetSizeIn : C;

  using    ScratchView = Kokkos::View       <Scalar*, DeviceType, Kokkos::MemoryUnmanaged>;

  const double alpha = 1.0;
  const double beta  = 0.0;

  if (!sumInto)
    Kokkos::deep_copy(outputVector, 0.0);

  int startCell = 0;
  while (startCell < C)
  {
    const int Cw = (worksetSize + startCell <= C) ? worksetSize : C - startCell;

    int workspace1_size = maxIntermediateSize_ * Cw * N;
    int workspace2_size = maxIntermediateSize_ * Cw * N;
    ScratchView    workspace1   (workspace.data(),                   workspace1_size);
    ScratchView    workspace2   (workspace.data() + workspace1_size, workspace2_size);
    INTREPID2_TEST_FOR_EXCEPTION(workspace1_size + workspace2_size > workspace.extent_int(0), std::invalid_argument, "Allocated workspace size is not sufficient");

    std::pair<int,int> cellRange = {startCell, startCell + Cw};

    // For vector-dot-vector integrals (e.g.), we need to integrate left x components against right x components, etc., and sum.
    // Each of these is one pass, and we accumulate in outputVector.
    const int numIntegrationPasses = int(componentIntegralsToSum_.size());

    using DynScratchView = Kokkos::DynRankView<Scalar,  DeviceType, Kokkos::MemoryUnmanaged>;
    int workspace3_size = Cw * F2 * N;
    int workspace4_size = Cw * F1 * N;
    int workspace5_size = Cw * F1 * N;
    DynScratchView workspace3, workspace4, workspace5;
    if (_orientations.size() > 0)
    {
      Teuchos::TimeMonitor orientationTimer = *Teuchos::TimeMonitor::getNewTimer("PAMatrix apply orientations right");

      int work_offset = workspace1_size + workspace2_size;
      workspace3 = DynScratchView(workspace.data() + work_offset, Cw, F2, N);
      work_offset += workspace3_size;
      workspace4 = DynScratchView(workspace.data() + work_offset, Cw, F1, N);
      work_offset += workspace4_size;
      workspace5 = DynScratchView(workspace.data() + work_offset, Cw, F1, N);
      work_offset += workspace5_size;
      INTREPID2_TEST_FOR_EXCEPTION(work_offset > workspace.extent_int(0), std::invalid_argument, "Allocated workspace size is not sufficient");
      // we accumulate in workspace4, so clear first:
      Kokkos::deep_copy(workspace4, 0.0);
      auto orientationsWorkset = Kokkos::subview(_orientations, cellRange);
      for (int n=0; n<N; n++)
      {
        auto inputVectorWorkset_n = Kokkos::subview( inputVector,   cellRange, Kokkos::ALL, n);
        auto workspace3_n         = Kokkos::subview( workspace3,  Kokkos::ALL, Kokkos::ALL, n);
        OrientationTools<DeviceType>::modifyBasisByOrientation(workspace3_n, inputVectorWorkset_n, orientationsWorkset,
                                                               _basisValuesRight.basisValues().getBasis().get());
      }
    }

    for (int integrationPass=0; integrationPass<numIntegrationPasses; integrationPass++)
    {
      const auto &integral_tuple = componentIntegralsToSum_[integrationPass];

      const auto & leftIntegrals = std::get<0>(integral_tuple);
      const PointDataSpec & pointDataSpec = std::get<1>(integral_tuple);
      // right integrals: replace F2j basis coefficients with evaluations at Pj
      const auto & rightIntegrals = std::get<2>(integral_tuple);
      int numRightIntegrals = int(rightIntegrals.size());

      const auto  leftOrdinalOffset = std::get<3>(integral_tuple);
      const auto     leftOutputSpan = std::get<4>(integral_tuple);
      const auto rightOrdinalOffset = std::get<5>(integral_tuple);
      const auto     rightInputSpan = std::get<6>(integral_tuple);

      int numLeftIntegrals = int(leftIntegrals.size());

      // set workspace1 to the (oriented) input data in an appropriate order
      // (C,F,N) arguments, where F=F_0…F_d, and the tensor ordering of these has F_0 as the fastest-moving index.
      // We want to start with input basis coefficients that are in a tensor product ordering with shape
      // (N,C,F), where the fastest-moving indices are on the left (LayoutLeft ordering).
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,rightInputSpan,N});
      if (_orientations.size() > 0)
      {
        Kokkos::parallel_for("PAMatrix::apply(): copy inputVector into workspace", policy,
                             KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                             {
          const ordinal_type idx = n + (c + f * Cw) * N;
          workspace1(idx) = workspace3(c,f+rightOrdinalOffset,n);
        });
      }
      else
      {
        Kokkos::parallel_for("PAMatrix::apply(): copy inputVector into workspace", policy,
                             KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                             {
          const ordinal_type idx = n + (c + f * Cw) * N;
          workspace1(idx) = inputVector(c + startCell,f+rightOrdinalOffset,n);
        });
      }
      ExecutionSpace().fence();

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
      int Nr = Cw * N * rightInputSpan;
      for (int j=0; j<numRightIntegrals; j++)
      {
        // we alternate whether we are placing intermediate results in workspace1 or workspace2
        auto  in = (j%2 == 0) ? workspace1 : workspace2;
        auto out = (j%2 == 0) ? workspace2 : workspace1;

        const int r = numRightIntegrals-j-1; // start with final component, work back to 0th component.

        auto op = rightIntegrals[r];

        // multiply P_r x F_r matrix with transpose of (… x F_r)
        const ordinal_type & Pr = op.M;
        const ordinal_type & Fr = op.N;

        const auto A = op.opView.data();
        const ordinal_type LDA = Pr; // will need to revise if we ever pad our operators (for byte alignment)
        const auto B = in.data();
        auto CData = out.data();
        INTREPID2_TEST_FOR_EXCEPTION(Nr % Fr != 0, std::invalid_argument, "Error: Nr must be a multiple of Fr");
        Impl::gemm<Impl::GemmDeviceType>('N', 'T', Pr, Nr/Fr, Fr, alpha, A, LDA, B, beta, CData);
        Nr = (Nr * Pr) / Fr;

      }
      auto  pointDataIn = (numRightIntegrals%2 == 0) ? workspace1 : workspace2; // pointwise result from contractions so far
      auto pointDataOut = (numRightIntegrals%2 == 0) ? workspace2 : workspace1; // pointwise output from weighting with pointData

      auto pointData = _pointDataCache.find(pointDataSpec)->second; // pointwise weights with shape (P[,Da[,Db]])
      // pointDataIn has shape (P,N,C); pointDataOut will have shape (N,C,P) (layout left).
      Impl::pointDataMultiply<DeviceType,Scalar>(pointDataSpec.C, pointDataSpec.P, pointDataSpec.aSpan, pointDataSpec.bSpan,
                                                 pointData.data(), pointDataIn.data(), pointDataOut.data(), N);

      Nr = Nr * pointDataSpec.aSpan / pointDataSpec.bSpan; // contracted in b, expanded in a


      for (int i=0; i<numLeftIntegrals; i++)
      {
        // we alternate whether we are placing intermediate results in workspace1 or workspace2
        auto  in = ((i+numRightIntegrals+1)%2 == 0) ? workspace1 : workspace2;
        auto out = ((i+numRightIntegrals+1)%2 == 0) ? workspace2 : workspace1;

        const int r = numRightIntegrals-i-1; // start with final component, work back to 0th component.
        auto op = leftIntegrals[r];

        // multiply Fr x Pr matrix with transpose of (… x Pr): result is (Fr x …)
        const ordinal_type & Fr = op.M;
        const ordinal_type & Pr = op.N;

        const auto A = op.opView.data();
        const ordinal_type LDA = Fr; // will need to revise if we ever pad our operators (for byte alignment)
        const auto B = in.data();
        auto CData = out.data();

        INTREPID2_TEST_FOR_EXCEPTION(Nr % Pr != 0, std::invalid_argument, "Error: Nr must be a multiple of Pr");
        Impl::gemm<Impl::GemmDeviceType>('N', 'T', Fr, Nr/Pr, Pr, alpha, A, LDA, B, beta, CData);
        Nr = (Nr * Fr) / Pr;

      }
      auto finalOut = ((numLeftIntegrals+numRightIntegrals+1)%2 == 0) ? workspace1 : workspace2;
      // Sum finalOut into outputVector
      policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,leftOutputSpan,N});
      if (_orientations.size() > 0)
      {
        // orientations will get applied in workspace4
        Kokkos::parallel_for("PAMatrix::apply(): sum finalOut into outputVector", policy,
                             KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                             {
          const ordinal_type idx = f + (n + c * N) * leftOutputSpan; // (F,N,C) layout-left
          workspace4(c,f+leftOrdinalOffset,n) += finalOut(idx);
        });
      }
      else
      {
        // if no orientations to apply, then we can accumulate directly int outputVector
        Kokkos::parallel_for("PAMatrix::apply(): sum finalOut into outputVector", policy,
                             KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
                             {
          const ordinal_type idx = f + (n + c * N) * leftOutputSpan; // (F,N,C) layout-left
          outputVector(c+startCell,f+leftOrdinalOffset,n) += finalOut(idx);
        });
      }
      ExecutionSpace().fence();
    } // integrationPass for loop
    if (_orientations.size() > 0)
    {
      auto orientationsWorkset = Kokkos::subview(_orientations, cellRange);
      for (int n=0; n<N; n++)
      {
        Teuchos::TimeMonitor orientationTimer = *Teuchos::TimeMonitor::getNewTimer("PAMatrix apply orientations left");
        auto outputVectorWorkset_n = Kokkos::subview( outputVector,   cellRange, Kokkos::ALL, n);
        auto workspace4_n          = Kokkos::subview( workspace4,   Kokkos::ALL, Kokkos::ALL, n);
        auto workspace5_n          = Kokkos::subview( workspace5,   Kokkos::ALL, Kokkos::ALL, n);
        OrientationTools<DeviceType>::modifyBasisByOrientation(workspace5_n, workspace4_n, orientationsWorkset,
                                                               _basisValuesLeft.basisValues().getBasis().get());
      }

      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{Cw,F1,N});
      Kokkos::parallel_for("PAMatrix::apply(): sum oriented output into outputVector", policy,
      KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &f, const ordinal_type &n)
      {
        outputVector(c+startCell,f,n) += workspace5(c,f,n);
      });
    }
    startCell += Cw;
  } // while (startCell < C)
}

namespace Impl
{
  //! process-wide GEMM statistics state for PAMatrix; see PAMatrix::setGEMMStatisticsEnabled()
  template<typename DeviceType,class Scalar>
  struct PAMatrixGEMMStatistics
  {
    static bool &enabled()         { static bool value = false; return value; }
    static double &cumulativeFlops() { static double value = 0;   return value; }
  };
}

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::setGEMMStatisticsEnabled(const bool enabled)
{
  Impl::PAMatrixGEMMStatistics<DeviceType,Scalar>::enabled() = enabled;
}

template<typename DeviceType,class Scalar>
bool PAMatrix<DeviceType,Scalar>::gemmStatisticsEnabled()
{
  return Impl::PAMatrixGEMMStatistics<DeviceType,Scalar>::enabled();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::recordGEMMFlops(const ordinal_type &M, const ordinal_type &N, const ordinal_type &K)
{
  double &cumulativeCount = Impl::PAMatrixGEMMStatistics<DeviceType,Scalar>::cumulativeFlops();
  if (gemmStatisticsEnabled())
  {
    // the product of M x K with K x N: each entry in the final matrix costs K multiplies and K - 1 adds
    cumulativeCount += double(M) * double(N) * double(2 * K - 1);
  }
  return cumulativeCount;
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmFlopCount()
{
  return Impl::PAMatrixGEMMStatistics<DeviceType,Scalar>::cumulativeFlops();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmTimeSeconds()
{
  // rank-local (no communication)
  auto timer = Teuchos::TimeMonitor::lookupCounter("gemm");
  return (timer == Teuchos::null) ? 0.0 : timer->totalElapsedTime();
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmThroughputGFlops(const double baseFlopCount, const double baseTimeSeconds)
{
  auto timeInSeconds = gemmTimeSeconds() - baseTimeSeconds;

  if (timeInSeconds > 0)
  {
    auto approximateFlopCountTotal = gemmFlopCount() - baseFlopCount;
    return approximateFlopCountTotal / timeInSeconds / 1.0e9;
  }
  else
  {
    return 0;
  }
}


} // end namespace Intrepid2
#endif
