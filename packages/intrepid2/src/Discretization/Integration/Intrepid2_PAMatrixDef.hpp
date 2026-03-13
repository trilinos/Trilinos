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

//#ifdef __APPLE__
//#include <Accelerate/Accelerate.h>
//#else
//#include <cblas.h>
//#endif

#include <Teuchos_BLAS.hpp>

#ifdef HAVE_INTREPID2_KOKKOSKERNELS
#include <KokkosBlas.hpp>
#endif

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
  
  // TODO: For non-trivial Da (≠1), we need to pack in the vector info into the point component for which the left integral will perform a dot product.  This will require some additional indexing logic, and we'll need additional arguments telling us about the point components: the number on the "left" of the vector point component, and either the number of points in the vector point component, or the number of points to the "right" of the vector point component.  For now, we throw an exception if aSpan != 1
  INTREPID2_TEST_FOR_EXCEPTION(aSpan != 1, std::invalid_argument, "aSpan != 1 is not yet supported");
  
  using ExecutionSpace = typename DeviceType::execution_space;
  auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<4>>({0,0,0,0},{N,numCells,numPoints,aSpan});
  
  Kokkos::parallel_for("pointDataMultiply", policy,
                       KOKKOS_LAMBDA(const ordinal_type &n, const ordinal_type &cell, const ordinal_type &point, const ordinal_type &a)
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
    Teuchos::TimeMonitor gemmTimer = *Teuchos::TimeMonitor::getNewTimer("gemm");
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
    Teuchos::TimeMonitor gemmTimer = *Teuchos::TimeMonitor::getNewTimer("gemm");
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

/*!
 Given tensor data with shape (D1, D2, …, Dn), prepare for a contraction in the Dk dimension by reordering as
  (Di, D1, D2, …, D{k-1}, D{k+1}, …, Dn)
 */
  template<typename DeviceType,class Scalar>
  class TensorReorderForGemmFunctor
  {
  public:
    using ExecutionSpace = typename DeviceType::execution_space;
    using View1D = Kokkos::View<Scalar*,DeviceType>;
    
    View1D outputView_;
    View1D  inputView_;
    
    int  leftDims_ = 1; // product D1 * D2 … * D{k-1}
    int      kDim_ = 1; // Dk
    int rightDims_ = 1; // product D{k+1} * … * Dn
    
    static constexpr bool layoutLeft_ = true; // aka column-major (Fortran-style): columns are together
    
    //! outputView and inputView must each be large enough to accommodate leftDims * iDim * rightDims, but they may be oversized.
    TensorReorderForGemmFunctor(View1D outputView, View1D inputView,
                                int leftDims, int kDim, int rightDims)
    :
    outputView_(outputView),
    inputView_(inputView),
    leftDims_(leftDims),
    kDim_(kDim),
    rightDims_(rightDims)
    {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const int &i, const int &j, const int &k) const
    {
      // source has (i,k,j)
      // dest   has (k,i,j)
      
      // i should iterate over leftDims (flattened)
      // j should iterate over rightDims (flattened)
      // k should iterate over Dk
      if (layoutLeft_) // column-major
      {
        const int dest_idx = k + (i + j * leftDims_) * kDim_;
        const int  src_idx = i + (k + j * kDim_    ) * leftDims_;
        outputView_(dest_idx) = inputView_(src_idx);
      }
      else
      {
        const int dest_idx = j + (i + k * leftDims_) * rightDims_;
        const int  src_idx = j + (k + i * kDim_    ) * rightDims_;
        outputView_(dest_idx) = inputView_(src_idx);
      }
    }
    
    void run()
    {
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{leftDims_,rightDims_,kDim_});
      Kokkos::parallel_for("PAMatrix: tensor reorder for gemm", policy, *this);
    }
  };

  template<typename DeviceType,class Scalar>
  class GemmSequenceFunctor
  {
  public:
    using ExecutionSpace = typename DeviceType::execution_space;
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    int numCells_;
    int numPoints_;
    int D1_ = -1;    // first dimension of pointwise data, if it is vector- or matrix-valued
    int D2_ = -1;    // second dimension of pointwise data, if it is matrix-valued
    int numVectors_; // N, the number of vectors in the multi-vector input/output
    int inputSize_;  // total number of  input data entries per cell
    int outputSize_; // total number of output data entries per cell
    int maxIntermediateSize_; // the highest entry count we'll need per cell as we apply operators
    int fad_size_output_;
    
    using View2D = Kokkos::View<  Scalar**,DeviceType>;
    using View3D = Kokkos::View< Scalar***,DeviceType>;
    using View4D = Kokkos::View<Scalar****,DeviceType>;
    
    Kokkos::Array<int,8> rightOpRowDims_; // point dimension
    Kokkos::Array<int,8> rightOpColDims_; // field dimension
    
    Kokkos::Array<int,8>  leftOpRowDims_; // field dimension
    Kokkos::Array<int,8>  leftOpColDims_; // point dimension
    
    Kokkos::Array< View2D, 8> refSpaceOpsRight_; // 2D views, shape (Pj,F2j) (ith component); if operator is vector-valued or tensor-valued, the Pi dimension will be packed (so that it will have more entries than the nominal point count)
    Kokkos::Array< View2D, 8> refSpaceOpsLeft_;  // 2D views, shape (F1i,Pi)
    
    int numOpsRight_;
    int numOpsLeft_;
    
    int pointDataRank_          = -1; // 0 for scalar data, 1 for vector, 2 for matrix
    int pointExpansionFactor_   =  1; // > 1 for vector/matrix-valued point data when right vector evaluation produces a scalar
    int pointContractionFactor_ =  1; // > 1 for vector/matrix-valued point data when right vector evaluation produces a vector
    
    View3D  inputView_; // shape (C,F2,N), where F2 is the row dimension of the full operator, and N is the number of vectors
    View3D outputView_; // shape (C,F1,N), where F1 is the column dimension of the full operator
    
    //! bottleneck constructor.  Only one of scalarPointData, vectorPointData, matrixPointData may be non-empty.
    GemmSequenceFunctor(View3D outputView, View3D inputView,
                        std::vector<View2D> refSpaceOpsRight, std::vector<View2D> refSpaceOpsLeft,
                        View2D scalarPointData, View3D vectorPointData, View4D matrixPointData)
    :
    inputView_(inputView),
    outputView_(outputView)
    {
      const bool allocateFadStorage = !(std::is_standard_layout<Scalar>::value && std::is_trivial<Scalar>::value);
      if (allocateFadStorage)
      {
        fad_size_output_ = dimension_scalar(inputView_);
      }
      numOpsRight_ = int(refSpaceOpsRight.size());
      INTREPID2_TEST_FOR_EXCEPTION(numOpsRight_ > refSpaceOpsRight_.size(), std::invalid_argument, "Too many right ops");
      
      numVectors_ = inputView_.extent_int(2);
      INTREPID2_TEST_FOR_EXCEPTION(numVectors_ != outputView_.extent_int(2), std::invalid_argument, "inputView and outputView must agree on the number of vectors");
      
      inputSize_  =  inputView_.extent_int(1) * numVectors_;   // F2 * N
      outputSize_ = outputView_.extent_int(1) * numVectors_;   // F1 * N
      int maxSize = max(inputSize_,outputSize_);
      int currentSize = inputSize_;
      for (int rj=0; rj<numOpsRight_; rj++)
      {
        const auto & rightOp = refSpaceOpsRight[rj];
        refSpaceOpsRight_[rj] = rightOp;
        // right ops convert from F2j dims to Pj dims
        rightOpRowDims_[rj] = rightOp.extent_int(0); // Pj
        rightOpColDims_[rj] = rightOp.extent_int(1); // F2j
        
        currentSize = currentSize / rightOpColDims_[rj] * rightOpRowDims_[rj];
        maxSize = max(currentSize, maxSize);
      }
      if (scalarPointData.size() > 0)
      {
        pointDataRank_ = 0;
        numPoints_ = scalarPointData.extent_int(1); // (C,P)
      }
      else if (vectorPointData.size() > 0)
      {
        pointDataRank_ = 1;
        numPoints_ = vectorPointData.extent_int(1); // (C,P,D)
        D1_        = vectorPointData.extent_int(2); // (C,P,D)
        
        if (currentSize == numVectors_ * numPoints_)
        {
          pointExpansionFactor_ = D1_;
        }
        else if (currentSize == numVectors_ * numPoints_ * D1_)
        {
          pointContractionFactor_ = D1_;
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "incompatible size sequence");
        }
      }
      else if (matrixPointData.size() > 0)
      {
        pointDataRank_ = 2;
        numPoints_ = matrixPointData.extent_int(1); // (C,P,D,D)
        D1_        = matrixPointData.extent_int(2); // (C,P,D,D)
        D2_        = matrixPointData.extent_int(3); // (C,P,D,D)
        
        if (currentSize == numVectors_ * numPoints_)
        {
          pointExpansionFactor_ = D1_ * D2_;
        }
        else if (currentSize == numVectors_ * numPoints_ * D2_)
        {
          // pointwise mat-vec: contract by D2, expand by D1
          pointContractionFactor_ = D2_;
          pointExpansionFactor_   = D1_;
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "incompatible size sequence");
        }
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "must specify scalar, vector, or matrix-valued point data");
      }
      currentSize *= pointExpansionFactor_;
      currentSize /= pointContractionFactor_;
      maxSize = max(currentSize, maxSize);
      
      numOpsLeft_  = int(refSpaceOpsLeft.size());
      INTREPID2_TEST_FOR_EXCEPTION(numOpsLeft_ > refSpaceOpsLeft_.size(), std::invalid_argument, "Too many left ops");
      for (int ri=0; ri<numOpsLeft_; ri++)
      {
        const auto & leftOp = refSpaceOpsLeft[ri];
        refSpaceOpsLeft_[ri] = leftOp;
        // left ops convert from Pi dims to F2i dims
        leftOpRowDims_[ri] = leftOp.extent_int(0); // F1i
        leftOpColDims_[ri] = leftOp.extent_int(1); // Pi
        
        currentSize = currentSize / leftOpColDims_[ri] * leftOpRowDims_[ri];
        maxSize = max(currentSize, maxSize);
      }
      INTREPID2_TEST_FOR_EXCEPTION(currentSize != outputSize_, std::invalid_argument, "Incompatible dimensions");
    }
    
    GemmSequenceFunctor(std::vector<View2D> refSpaceOpsRight, View2D scalarPointData, std::vector<View2D> refSpaceOpsLeft)
    :
    GemmSequenceFunctor(refSpaceOpsRight, refSpaceOpsLeft, scalarPointData, View3D(), View4D())
    {}
    
    GemmSequenceFunctor(std::vector<View2D> refSpaceOpsRight, View3D vectorPointData, std::vector<View2D> refSpaceOpsLeft)
    :
    GemmSequenceFunctor(refSpaceOpsRight, refSpaceOpsLeft, View2D(), vectorPointData, View4D())
    {}
    
    GemmSequenceFunctor(std::vector<View2D> refSpaceOpsRight, View4D matrixPointData, std::vector<View2D> refSpaceOpsLeft)
    :
    GemmSequenceFunctor(refSpaceOpsRight, refSpaceOpsLeft, View2D(), View3D(), matrixPointData)
    {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      const int cellOrdinal  = teamMember.league_rank();
      const int threadNumber = teamMember.team_rank();
      const int numThreads   = teamMember.team_size(); // num threads
      
      using ScratchView = Kokkos::View<Scalar*, DeviceType, Kokkos::MemoryUnmanaged>;
      
      // we alternate between using these two workspaces as the destination for the gemms
      ScratchView workspace1;
      ScratchView workspace2;
      
      if (fad_size_output_ > 0) 
      {
        workspace1 = ScratchView(teamMember.team_shmem(), maxIntermediateSize_, fad_size_output_);
        workspace2 = ScratchView(teamMember.team_shmem(), maxIntermediateSize_, fad_size_output_);
      }
      else 
      {
        workspace1 = ScratchView(teamMember.team_shmem(), maxIntermediateSize_);
        workspace2 = ScratchView(teamMember.team_shmem(), maxIntermediateSize_);
      }
      
      // TODO: sequence of gemms.  We will need to either move memory to allow the gemms to take place on the whole structure,
      //       or to specify gemms on contiguous chunks; due to the slicing of the tensor, in general the input data does not
      //       have structure of a monolithic matrix, but it does have matrix blocks whose products can be placed into a result
      //       in a memory-contiguous fashion.  Splitting into smaller gemms *might* allow expression of more parallelism, but
      //       I doubt we can beat vendor-provided implementations.
      
//      Scalar result;
//      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember,0,maxFields), [&] (const int& fieldOrdinal, Scalar &contractionThusFar)
//      {
//        
//      }, result);
      
      // synchronize threads
      teamMember.team_barrier();
    }
    
  };

// blas.GEMM(trA, trB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

  //! take an M x K matrix A and contract with an N1 x K x N2 tensor B to produce a M x N1 x N2 output C.
  //! This is done in terms of a series of constituent gemms, iterating over the n2 dimension.  We launch these in a Kokkos::parallel_for on the
  //! *host* execution space.  This allows us to invoke a synchronous gemm call in an asynchronous way.  In particular, on macOS, Apple's Accelerate
  //! framework provides a gemm implementation that invokes the GPU on M-series processors, but this waits for completion before it returns, and
  //! typically does not saturate the GPU.  If Kokkos is built with OpenMP support, we can thus increase parallelism by however many OpenMP threads
  //! are available.  Similar considerations apply to KokkosKernels's gemm implementation on CUDA or HIP DeviceType.  We do need to be careful not to
  //! launch a KokkosKernels gemm under OpenMP with an OpenMP DispatchExecutionSpace: the basic rule here is that GemmDeviceType must
//! be different from DispatchExecutionSpace unless they are both Serial.
  template<typename GemmDeviceType, class Scalar, typename DispatchExecutionSpace=Kokkos::DefaultHostExecutionSpace>
  std::enable_if_t<
    !std::is_same<DispatchExecutionSpace, GemmDeviceType>::value ||
    (std::is_same<DispatchExecutionSpace, Kokkos::Serial>::value && std::is_same<GemmDeviceType, Kokkos::Serial>::value)
  >
  matrixTensorContractionLayoutLeft(const ordinal_type &M, const ordinal_type &N1, const ordinal_type &N2, const ordinal_type &K,
                                    const Scalar &alpha, const Scalar* A, const ordinal_type &LDA,
                                    const Scalar *B,
                                    const Scalar &beta, Scalar *C)
  {
    // TODO: once everything is working, assuming we have not found a use for this code, delete it and the tests against it.
    // We assume layout left, so that the B tensor (i,k,j) index flattens to i + (k + j * K) * N1.
    // This means that the slice B(:,:,j) is a N1 x K matrix, contiguous in memory, at offset j * K * N1.
    // C(m,i,j) -> m + (i + j * N1) * M
    // Similarly, the slice C(:,:,j) is a M x N1 matrix, contiguous in memory, at offset j * N1 * M.
    auto policy = Kokkos::RangePolicy<DispatchExecutionSpace>(0,N2);
    
    const ordinal_type KN1 = K * N1;
    const ordinal_type N1M = N1 * M;
    Kokkos::parallel_for("matrixTensorContractionLayoutLeft: GEMM dispatch", policy,
                         KOKKOS_LAMBDA(const ordinal_type &j)
    {
      const auto B_j = B + j * KN1;
      const auto C_j = C + j * N1M;
      gemm<GemmDeviceType>('N', 'T', M, N1, K, alpha, A, LDA, B_j, beta, C_j);
      
      using namespace std;
//      cout << "j = " << j << std::endl;
//      cout << "computed A * B^T = C:\n";
//      cout << "A:\n";
//      for (int m=0; m<M; m++)
//      {
//        cout << "[ ";
//        for (int k=0; k<K; k++)
//        {
//          cout << *(A + m + k * M) << " ";
//        }
//        cout << "]\n";
//      }
//      cout << "B:\n";
//      for (int n1=0; n1<N1; n1++)
//      {
//        cout << "[ ";
//        for (int k=0; k<K; k++)
//        {
//          cout << *(B_j + n1 + k * N1) << " ";
//        }
//        cout << "]\n";
//      }
//      
//      cout << "C:\n";
//      for (int n1=0; n1<N1; n1++)
//      {
//        cout << "[ ";
//        for (int m=0; m<M; m++)
//        {
//          cout << *(C_j + m + n1 * M) << " ";
//        }
//        cout << "]\n";
//      }
    });
    
    DispatchExecutionSpace().fence();
  }

  //! Given (C,P[,D,D]) transform and (C,P) pointwise weights, construct a suitable container for storing the pointwise weighted transform.
  template<typename DeviceType,class Scalar>
  Data<Scalar,DeviceType> allocateComposedWeightedTransform(const Data<Scalar,DeviceType> &composedTransform,
                                                            const TensorData<Scalar,DeviceType> &pointWeights)
  {
    auto cellDimInfo = composedTransform.getDimensionInfo(0); // cell dimension
    int numTensorComponents = pointWeights.numTensorComponents();
    const int & numLogicalCells = cellDimInfo.logicalExtent;
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
void PAMatrix<DeviceType,Scalar>::init(const TransformedBasisValues<Scalar,DeviceType> basisValuesLeft,
                                  const TensorData<Scalar,DeviceType> cellMeasures,
                                  const TransformedBasisValues<Scalar,DeviceType> basisValuesRight,
                                  const ScalarView<Orientation,DeviceType> orientations) {
  using ExecutionSpace = typename DeviceType::execution_space;

  const bool layoutLeft = layoutLeft_;
  
  const bool  leftHasOrdinalFilter =  basisValuesLeft.basisValues().ordinalFilter().extent_int(0) > 0;
  const bool rightHasOrdinalFilter = basisValuesRight.basisValues().ordinalFilter().extent_int(0) > 0;
  TEUCHOS_TEST_FOR_EXCEPTION(leftHasOrdinalFilter || rightHasOrdinalFilter, std::invalid_argument, "Ordinal filters for BasisValues are not yet supported by PAMatrix");
  
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
  const int numPointTensorComponents = cellMeasures.numTensorComponents() - 1;
    
  // for now, we don't check for separability: we always treat as non-separable.
  // we could gain some performance in the not-too-common case that the integrals are separable; I haven't yet figured out how much.
//  if ((numPointTensorComponents == numTensorComponentsLeft) && basisValuesLeft.axisAligned() && basisValuesRight.axisAligned())
//  {
//    _separable = true;
//  }
//  else // general case (not axis-aligned + affine tensor-product structure)
  {
    _separable = false;
    // MARK: prepare composed transformation matrices
    const Data<Scalar,DeviceType> & leftTransform  = basisValuesLeft.transform();
    const Data<Scalar,DeviceType> & rightTransform = basisValuesRight.transform();
    const bool transposeLeft  = true;
    const bool transposeRight = false;
    //    auto timer = Teuchos::TimeMonitor::getNewTimer("mat-mat");
    //    timer->start();
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
          // extents[3]               = 1;
          // variationTypes[3]        = CONSTANT;
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
    bool haveLaunchedContributionToCurrentFamilyLeft = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
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
      
      int rightFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families // TODO: figure out what a nonzero value means for matrix-free apply() implementation
      for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
      {
        const int rightFieldSpan = rightIsVectorValued ? basisValuesRight.vectorData() .numFieldsInFamily(rightFamilyOrdinal)
                                                       : basisValuesRight.basisValues().numFieldsInFamily(rightFamilyOrdinal);
        // "b" keeps track of the spatial dimension over which we are integrating in the right vector
        // components are allowed to span several dimensions; we keep track of the offset for the component in b_offset
        bool haveLaunchedContributionToCurrentFamilyRight = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
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
  
  // Ordinal filter is used for Serendipity basis; we don't yet support Serendipity for PAMatrix.
  // (When we do, the strategy will likely be to apply the right filter at the "middle" of the operator sequence, and the left filter at the end.  This does mean that the intermediate containers for the right operators will be sized for the unfiltered basis; the intermediate containers for the left operators will be sized like unfiltered left x filtered right.)
  const bool  leftHasOrdinalFilter =  _basisValuesLeft.basisValues().ordinalFilter().extent_int(0) > 0;
  const bool rightHasOrdinalFilter = _basisValuesRight.basisValues().ordinalFilter().extent_int(0) > 0;
  TEUCHOS_TEST_FOR_EXCEPTION(leftHasOrdinalFilter || rightHasOrdinalFilter, std::invalid_argument, "Ordinal filters for BasisValues are not yet supported by PAMatrix");
  
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
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize) const
{
  using View1D = Kokkos::View<Scalar*,DeviceType>;
  int size1D = maxIntermediateSize_ * worksetSize * 2;
  
  if (_orientations.size() > 0)
  {
    // then we need to apply orientations on the way in and on the way out, and we need workspace to do that
    size1D +=  _basisValuesLeft.numFields() * worksetSize * 2; // to support sumInto = true, we need intermediate storage for oriented left values (hence the factor of 2).
    size1D += _basisValuesRight.numFields() * worksetSize;
  }
//  std::cout  << "allocated " << size1D << " entries of workspace.\n";
    
  return View1D("PAMatrix workspace", size1D);
}

template<typename DeviceType,class Scalar>
Kokkos::View<Scalar*,DeviceType> PAMatrix<DeviceType,Scalar>::allocateWorkspace(const ordinal_type &worksetSize,
                                                                                const ordinal_type &n) const
{
  using View1D = Kokkos::View<Scalar*,DeviceType>;
  int size1D = maxIntermediateSize_ * worksetSize * n * 2;
  
  if (_orientations.size() > 0)
  {
    // then we need to apply orientations on the way in and on the way out, and we need workspace to do that
    size1D +=  _basisValuesLeft.numFields() * worksetSize * n * 2; // to support sumInto = true, we need intermediate storage for oriented left values (hence the factor of 2).
    size1D += _basisValuesRight.numFields() * worksetSize * n;
  }
  
//  std::cout  << "allocated " << size1D << " entries of workspace.\n";
  
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
  // TODO: revise to take a single workspace argument.  We should manage subdivision internally.
  
  using ExecutionSpace = typename DeviceType::execution_space;
  using View1D = Kokkos::View<Scalar*,DeviceType>;
  
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
      int work_offset = workspace1_size + workspace2_size;
      workspace3 = DynScratchView(workspace.data() + work_offset, Cw, F2, N);
      work_offset += workspace3_size;
      workspace4 = DynScratchView(workspace.data() + work_offset, Cw, F1, N);
      work_offset += workspace4_size;
      workspace5 = DynScratchView(workspace.data() + work_offset, Cw, F1, N);
      work_offset += workspace5_size;
//      std::cout  << "using " << work_offset << " entries of workspace.\n";
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
        auto C = out.data();
        INTREPID2_TEST_FOR_EXCEPTION(Nr % Fr != 0, std::invalid_argument, "Error: Nr must be a multiple of Fr");
        Impl::gemm<Impl::GemmDeviceType>('N', 'T', Pr, Nr/Fr, Fr, alpha, A, LDA, B, beta, C);
        //      Impl::matrixTensorContractionLayoutLeft<Impl::GemmDeviceType>(Fr, N1, N2, Pr, alpha, A, LDA, B, beta, C);
        Nr = (Nr * Pr) / Fr;
        
        //      {
        //        // DEBUGGING
        //        using namespace std;
        //        cout << "j=" << j << ", result: [";
        //        for (int i=0; i<Nr; i++)
        //        {
        //          cout << out(i) << " ";
        //        }
        //        cout << "]\n";
        //      }
      }
      auto  pointDataIn = (numRightIntegrals%2 == 0) ? workspace1 : workspace2; // pointwise result from contractions so far
      auto pointDataOut = (numRightIntegrals%2 == 0) ? workspace2 : workspace1; // pointwise output from weighting with pointData
      
      auto pointData = _pointDataCache.find(pointDataSpec)->second; // pointwise weights with shape (P[,Da[,Db]])
      // pointDataIn has shape (P,N,C); pointDataOut will have shape (N,C,P) (layout left).
      Impl::pointDataMultiply<DeviceType,Scalar>(pointDataSpec.C, pointDataSpec.P, pointDataSpec.aSpan, pointDataSpec.bSpan,
                                                 pointData.data(), pointDataIn.data(), pointDataOut.data(), N);
      
      Nr = Nr * pointDataSpec.aSpan / pointDataSpec.bSpan; // contracted in b, expanded in a
      
      //    {
      //      // DEBUGGING
      //      using namespace std;
      //      cout << "After pointData contraction, result: [";
      //      for (int j=0; j<Nr; j++)
      //      {
      //        cout << pointDataOut(j) << " ";
      //      }
      //      cout << "]\n";
      //    }
      
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
        auto C = out.data();
        //      {
        //        // DEBUGGING
        //        using namespace std;
        //        cout << "i=" << i << ", A: [";
        //        for (int j=0; j<Fr*Pr; j++)
        //        {
        //          cout << op.opView(j) << " ";
        //        }
        //        cout << "]\n";
        //        cout << "B: [";
        //        for (int j=0; j<Nr; j++)
        //        {
        //          cout << in(j) << " ";
        //        }
        //        cout << "]\n";
        //        cout << "gemm args: ('N', 'T', " << Fr << ", " << Nr/Pr << ", " << Pr << ", " << alpha << ", A, " << LDA << ", B, " << beta << ", C)\n";
        //      }
        
        INTREPID2_TEST_FOR_EXCEPTION(Nr % Pr != 0, std::invalid_argument, "Error: Nr must be a multiple of Pr");
        Impl::gemm<Impl::GemmDeviceType>('N', 'T', Fr, Nr/Pr, Pr, alpha, A, LDA, B, beta, C);
        //      Impl::matrixTensorContractionLayoutLeft<Impl::GemmDeviceType>(Fr, N1, N2, Pr, alpha, A, LDA, B, beta, C);
        Nr = (Nr * Fr) / Pr;
        
        //      {
        //        // DEBUGGING
        //        using namespace std;
        //        cout << "i=" << i << ", result: [";
        //        for (int j=0; j<Nr; j++)
        //        {
        //          cout << out(j) << " ";
        //        }
        //        cout << "]\n";
        //      }
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

template<typename DeviceType,class Scalar>
void PAMatrix<DeviceType,Scalar>::assemble(Data<Scalar,DeviceType> &integrals) const
{
  //placeholder implementation: just invoke IntegrationTools
  using ExecutionSpace = typename DeviceType::execution_space;
  using MemorySpace    = typename DeviceType::memory_space;
  
  bool sumInto = false;
  double approximateFlopCountIntegrate = 0;
  IntegrationTools<DeviceType>::integrate(integrals, _basisValuesLeft, _cellMeasures, _basisValuesRight, sumInto, &approximateFlopCountIntegrate);
  ExecutionSpace().fence();
  
  auto leftBasis  =  _basisValuesLeft.basisValues().getBasis();
  auto rightBasis = _basisValuesRight.basisValues().getBasis();
  
  if (_orientations.size() > 0)
  {
    // modify integrals by orientations -- we are NOT allowed to use the same view as source and result, so let's create a mirror view for source.
    auto unorientatedValues = Kokkos::create_mirror_view_and_copy(MemorySpace(), integrals.getUnderlyingView());
    OrientationTools<DeviceType>::modifyMatrixByOrientation(integrals.getUnderlyingView(), unorientatedValues,
                                                            _orientations, leftBasis.get(), rightBasis.get());
    ExecutionSpace().fence();
  }
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::recordGEMMFlops(const ordinal_type &M, const ordinal_type &N, const ordinal_type &K)
{
  static double cumulativeCount = 0;
  
  // compute the product of  M x K with K x N: each entry in the final matrix costs K multiplies and K - 1 adds
  const double approximateFlops = M * N * (2 * K - 1);
  cumulativeCount += approximateFlops;
  return cumulativeCount;
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmFlopCount()
{
  auto approximateFlopCountTotal = recordGEMMFlops(0,0,0);
  return approximateFlopCountTotal;
}

template<typename DeviceType,class Scalar>
double PAMatrix<DeviceType,Scalar>::gemmTimeSeconds() // returns cumulative time spent in gemm() calls
{
  Teuchos::stat_map_type statData;
  std::vector<std::string> statNames;
  Teuchos::TimeMonitor::computeGlobalTimerStatistics(statData, statNames, Teuchos::Intersection, "gemm");
  
  if (statData["gemm"].size() > 0)
  {
    const double timeInSeconds = statData["gemm"][0].first;
    return timeInSeconds;
  }
  else
  {
    return 0;
  }
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
