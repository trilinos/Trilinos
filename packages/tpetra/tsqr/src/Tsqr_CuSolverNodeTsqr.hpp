// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_CuSolverNodeTsqr.hpp
/// \brief Declaration and definition of CuSolverNodeTsqr.

#ifndef TSQR_CUSOLVERNODETSQR_HPP
#define TSQR_CUSOLVERNODETSQR_HPP

#include "TpetraTSQR_config.h"

#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
#include "Tsqr_NodeTsqr.hpp"
#include "Tsqr_Impl_CuBlas.hpp"
#include "Tsqr_Impl_CuSolver.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <memory>
#include <type_traits>

#define TSQR_IMPL_CATCH( message ) \
  catch (std::exception& e) { \
    threw = true; \
    err = std::unique_ptr<std::ostringstream> (new std::ostringstream); \
    *err << prefix << message << std::endl << e.what (); \
  } \
  TEUCHOS_TEST_FOR_EXCEPTION \
    (threw, std::runtime_error, \
     (err.get () == nullptr ? "Unknown error" : err->str ())); \
  do {} while (false)

#define TSQR_IMPL_CHECK_LAST_CUDA_ERROR( location ) \
  do { \
    cudaError_t errCode = cudaGetLastError (); \
    if (errCode != cudaSuccess ) { \
      const char* errorString = cudaGetErrorString (errCode); \
      TEUCHOS_TEST_FOR_EXCEPTION \
        (true, std::runtime_error, "At \"" << (location) << "\", " \
         "CUDA is in the following error state: " << errorString); \
    } \
  } while (false)

namespace TSQR {
  namespace Impl {

    using cusolver_memory_space = Kokkos::CudaSpace;
    using cusolver_execution_space = Kokkos::Cuda;
    using host_device_type = Kokkos::Device<
      Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;

    // Mapping from Scalar to Kokkos value type.
    // e.g., Scalar=std::complex<double> -> Kokkos::complex<double>.

    template<class Scalar>
    using non_const_kokkos_value_type = typename Kokkos::ArithTraits<
        typename std::remove_const<Scalar>::type
      >::val_type;

    template<class Scalar>
    using kokkos_view_value_type = typename std::conditional<
        std::is_const<Scalar>::value,
        const non_const_kokkos_value_type<Scalar>,
        non_const_kokkos_value_type<Scalar>
      >::type;

    // vector_type, device_vector_type, and host_vector_type

    template<class T, class MemorySpace>
    using vector_type = Kokkos::View<T*, MemorySpace>;

    template<class T>
    using device_vector_type = vector_type<T, cusolver_memory_space>;

    template<class T>
    using host_vector_type = vector_type<T, host_device_type>;

    template<class T>
    void
    reallocDeviceVectorIfNeeded (device_vector_type<T>& vec,
                                 const char label[],
                                 const size_t minSize)
    {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;

      if (size_t (vec.size ()) < minSize) {
        vec = device_vector_type<T> ();
        auto alloc = view_alloc (std::string (label), WithoutInitializing);
        vec = device_vector_type<T> (alloc, minSize);
      }
    }

    // vec_view_type & device_vec_view_type

    template<class T, class MemorySpace>
    using vec_view_type =
      Kokkos::View<T*, MemorySpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    template<class T>
    using device_vec_view_type = vec_view_type<T, cusolver_memory_space>;

    // matrix_type & device_matrix_type

    template<class T, class MemorySpace>
    using matrix_type = Kokkos::View<T**, Kokkos::LayoutLeft, MemorySpace>;

    template<class T>
    using device_matrix_type = matrix_type<T, cusolver_memory_space>;

    // mat_view_type, device_mat_view_type, & host_mat_view_type

    template<class T, class MemorySpace>
    using mat_view_type =
      Kokkos::View<T**, Kokkos::LayoutLeft, MemorySpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    template<class T>
    using device_mat_view_type =
      mat_view_type<T, cusolver_memory_space>;

    template<class T>
    using host_mat_view_type = mat_view_type<T, host_device_type>;

    // get_mat_view, get_host_mat_view, & get_device_mat_view

    template<class Scalar, class MemorySpace>
    static mat_view_type<kokkos_view_value_type<Scalar>, MemorySpace>
    get_mat_view (const size_t nrows,
                  const size_t ncols,
                  Scalar A[],
                  const size_t lda)
    {
      static_assert
        (! std::is_const<non_const_kokkos_value_type<Scalar> >::value,
         "non_const_kokkos_value_type is const.");
      using KVVT = kokkos_view_value_type<Scalar>; // preserves const
      static_assert
        ((std::is_const<Scalar>::value && std::is_const<KVVT>::value) ||
         (! std::is_const<Scalar>::value && ! std::is_const<KVVT>::value),
         "kokkos_view_value_type failed to preserve const-ness.");
      KVVT* A_raw = reinterpret_cast<KVVT*> (A);

      mat_view_type<KVVT, MemorySpace> A_full (A_raw, lda, ncols);
      const std::pair<size_t, size_t> rowRange (0, nrows);
      return Kokkos::subview (A_full, rowRange, Kokkos::ALL ());
    }

    template<class Scalar>
    static host_mat_view_type<kokkos_view_value_type<Scalar>>
    get_host_mat_view (const size_t nrows,
                       const size_t ncols,
                       Scalar A[],
                       const size_t lda)
    {
      return get_mat_view<Scalar, host_device_type>
        (nrows, ncols, A, lda);
    }

    template<class Scalar, class Ordinal>
    static host_mat_view_type<kokkos_view_value_type<Scalar>>
    get_host_mat_view (const MatView<Ordinal, Scalar>& A_host)
    {
      const size_t nrows (A_host.extent (0));
      const size_t ncols (A_host.extent (1));
      const size_t lda (A_host.stride (1));
      return get_mat_view<Scalar, host_device_type>
        (nrows, ncols, A_host.data (), lda);
    }

    template<class Scalar>
    static device_mat_view_type<kokkos_view_value_type<Scalar>>
    get_device_mat_view (const size_t nrows,
                         const size_t ncols,
                         Scalar A[],
                         const size_t lda)
    {
      return get_mat_view<Scalar, cusolver_memory_space> (nrows, ncols, A, lda);
    }

    /// \brief Given rank-1 backing storage, return a device matrix
    ///   view with the given dimensions (numRows by numCols), that
    ///   has contiguous storage.  Reallocate storage if needed.
    ///
    /// "Contiguous storage" means that if A is the matrix view
    /// result, then A.stride(1) == A.extent(0).
    template<class T>
    device_mat_view_type<T>
    get_contiguous_device_mat_view (device_vector_type<T>& storage,
                                    const size_t numRows,
                                    const size_t numCols)
    {
      const char prefix[] = "TSQR::Impl::get_contiguous_device_mat_view: ";

      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( prefix );

      const size_t currentStorageSize (storage.extent (0));
      const size_t requiredStorageSize = numRows * numCols;
      if (currentStorageSize < requiredStorageSize) {
        // It costs about as much to allocate 8B on device as 800B.
        constexpr size_t minStorageSize = 100;
        const size_t newStorageSize =
          std::max (minStorageSize, requiredStorageSize);

        // Free it first, so that two allocations won't coexist.
        storage = device_vector_type<T> ();
        using Kokkos::view_alloc;
        using Kokkos::WithoutInitializing;
        const char label[] = "matrixStorage";

        TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::Impl::get_contiguous_device_mat_view: Right before allocating" );

        try {
          storage = device_vector_type<T>
            (view_alloc (std::string (label), WithoutInitializing),
             newStorageSize);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, prefix << "Allocating rank-1 "
             "View of size " << newStorageSize << " to represent a "
             << numRows << " x " << numCols << " matrix threw: "
             << std::endl << e.what ());
        }
      }
      return device_mat_view_type<T> (storage.data (),
                                      numRows, numCols);
    }

    template<class T>
    host_mat_view_type<T>
    get_contiguous_host_mat_view (host_vector_type<T>& storage,
                                  const size_t numRows,
                                  const size_t numCols)
    {
      const char prefix[] = "TSQR::Impl::get_contiguous_host_mat_view: ";

      const size_t currentStorageSize (storage.extent (0));
      const size_t requiredStorageSize = numRows * numCols;
      if (currentStorageSize < requiredStorageSize) {
        // It costs about as much to allocate 8B on host as 800B.
        constexpr size_t minStorageSize = 100;
        const size_t newStorageSize =
          std::max (minStorageSize, requiredStorageSize);

        // Free it first, so that two allocations won't coexist.
        storage = host_vector_type<T> ();
        using Kokkos::view_alloc;
        using Kokkos::WithoutInitializing;
        const char label[] = "hostMatrixStorage";

        try {
          storage = host_vector_type<T>
            (view_alloc (std::string (label), WithoutInitializing),
             newStorageSize);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, prefix << "Allocating rank-1 "
             "host View of size " << newStorageSize << " to store a "
             << numRows << " x " << numCols << " matrix threw: "
             << std::endl << e.what ());
        }
      }
      return host_mat_view_type<T> (storage.data (),
                                    numRows, numCols);
    }

    // info_type & const_info_type

    using info_type = Kokkos::View<int, cusolver_memory_space>;
    using const_info_type = Kokkos::View<const int, cusolver_memory_space>;

    template<class LocalOrdinal, class Scalar>
    class CuSolverNodeFactorOutput :
      public NodeFactorOutput<LocalOrdinal, Scalar>
    {
    public:
      //using cuda_value_type = typename Impl::CudaValue<Scalar>::type;
      using kokkos_value_type = non_const_kokkos_value_type<Scalar>;
      using const_tau_type = device_vector_type<const kokkos_value_type>;
      using const_unmanaged_tau_type =
        device_vec_view_type<const kokkos_value_type>;

      CuSolverNodeFactorOutput (const const_tau_type& tau,
                                const const_info_type& info) :
        tau_ (tau), info_ (info)
      {}

      const_unmanaged_tau_type tau () const { return tau_; }

      int info () const {
        int info_h = 0;
        Kokkos::deep_copy (info_h, info_);
        return info_h;
      }

    private:
      const_tau_type tau_;
      const_info_type info_;
    };

    template<class ScalarType, class IndexType>
    class SetDiagonalEntriesToOne {
      static_assert (! std::is_const<ScalarType>::value,
        "SetDiagonalEntriesToOne requires a View of nonconst.");
    public:
      SetDiagonalEntriesToOne
        (const device_mat_view_type<ScalarType>& A) : A_ (A) {}
      KOKKOS_INLINE_FUNCTION void
      operator() (const IndexType j) const {
        A_(j,j) = ScalarType (1.0);
      }
    private:
      device_mat_view_type<ScalarType> A_;
    };

    template<class ScalarType>
    void
    set_diagonal_entries_to_one
      (const device_mat_view_type<ScalarType>& A)
    {
      static_assert (! std::is_const<ScalarType>::value,
        "set_diagonal_entries_to_one requires a View of nonconst.");
      using LO =
        typename std::make_signed<decltype (A.extent (1)) >::type;
      const LO ncols = std::min (A.extent (0), A.extent (1));
      using Kokkos::RangePolicy;
      RangePolicy<cusolver_execution_space, LO> range (0, ncols);
      Kokkos::parallel_for
        ("set_diagonal_entries_to_one", range,
         SetDiagonalEntriesToOne<ScalarType, LO> (A));
    }

  } // namespace Impl

  /// \class CuSolverNodeTsqr
  /// \brief NodeTsqr implementation based on cuSOLVER.
  /// \author Mark Hoemmen
  template<class LocalOrdinal, class Scalar>
  class CuSolverNodeTsqr : public NodeTsqr<LocalOrdinal, Scalar>
  {
  private:
    using base_type = NodeTsqr<LocalOrdinal, Scalar>;
    using my_factor_output_type =
      Impl::CuSolverNodeFactorOutput<LocalOrdinal, Scalar>;
    using kokkos_value_type =
      Impl::non_const_kokkos_value_type<Scalar>;

  public:
    using ordinal_type = typename base_type::ordinal_type;
    using scalar_type = typename base_type::scalar_type;
    using factor_output_type = typename base_type::factor_output_type;

    CuSolverNodeTsqr () = default;

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const override
    {
      return Teuchos::parameterList ("NodeTsqr");
    }

    void
    setParameterList
      (const Teuchos::RCP<Teuchos::ParameterList>&) override
    {}

    std::string description () const override {
      return "CuSolverNodeTsqr";
    }

    bool wants_device_memory () const override { return true; }

    bool ready () const override {
      return true;
    }

    bool
    QR_produces_R_factor_with_nonnegative_diagonal () const override
    {
      return false;
    }

    size_t cache_size_hint () const override {
      return 0;
    }

  private:
    using tau_type = Impl::device_vector_type<kokkos_value_type>;

    // must return owning, since we'll pass off to factor output
    tau_type
    get_tau (const LocalOrdinal numCols) const
    {
      using Impl::reallocDeviceVectorIfNeeded;
      Impl::reallocDeviceVectorIfNeeded (tau_, "tau", size_t (numCols));
      return tau_;
    }

    using work_type = Impl::device_vector_type<kokkos_value_type>;
    using nonowning_work_type =
      Impl::device_vec_view_type<kokkos_value_type>;

    nonowning_work_type
    get_work_for_geqrf (const LocalOrdinal numRows,
                        const LocalOrdinal numCols,
                        Scalar A[],
                        const LocalOrdinal lda) const
    {
      auto info = get_info ();
      TSQR::Impl::CuSolver<Scalar> solver (info.data ());
      const int lwork =
        solver.compute_QR_lwork (numRows, numCols, A, lda);
      // Avoid constant reallocation by setting a minimum lwork.
      constexpr int min_lwork = 128;
      const int new_lwork = lwork < min_lwork ? min_lwork : lwork;
      using Impl::reallocDeviceVectorIfNeeded;
      reallocDeviceVectorIfNeeded (work_, "work", new_lwork);
      return nonowning_work_type (work_);
    }

    nonowning_work_type
    get_work_for_apply_Q_factor (const ApplyType& apply_type,
                                 const LocalOrdinal nrows,
                                 const LocalOrdinal ncols_C,
                                 const LocalOrdinal ncols_Q,
                                 const Scalar A[],
                                 const LocalOrdinal lda,
                                 const Scalar tau[],
                                 Scalar C[],
                                 const LocalOrdinal ldc) const
    {
      auto info = get_info ();
      TSQR::Impl::CuSolver<Scalar> solver (info.data ());
      const char side = 'L';
      const char trans = apply_type.toString ()[0];
      const int lwork =
        solver.apply_Q_factor_lwork (side, trans,
                                     nrows, ncols_C, ncols_Q,
                                     A, lda, tau, C, ldc);
      // Avoid constant reallocation by setting a minimum lwork.
      constexpr int min_lwork = 128;
      const int new_lwork = lwork < min_lwork ? min_lwork : lwork;
      using Impl::reallocDeviceVectorIfNeeded;
      reallocDeviceVectorIfNeeded (work_, "work", new_lwork);
      return nonowning_work_type (work_);
    }

    // must return owning, since we'll pass off to factor output
    Impl::info_type
    get_info () const
    {
      if (info_.data () == nullptr) {
        info_ = Impl::info_type ("info");
      }
      // "get last error" model will avoid doing multiple info allocations.
      return info_;
    }

    Impl::device_mat_view_type<kokkos_value_type>
    get_Q_copy (const LocalOrdinal nrows,
                const LocalOrdinal ncols,
                const Scalar Q[], // DEVICE MEMORY
                const LocalOrdinal ldq) const
    {
      using Impl::get_contiguous_device_mat_view;
      auto Q_copy =
        get_contiguous_device_mat_view (matrixStorage_, nrows, ncols);
      auto Q_view = Impl::get_device_mat_view (nrows, ncols, Q, ldq);
      // NOTE (mfh 17 Dec 2019) We're copying device to device, so the
      // Kokkos::deep_copy noncontiguity problem does not apply.
      Kokkos::deep_copy (Q_copy, Q_view);
      return Q_copy;
    }

    Impl::device_mat_view_type<kokkos_value_type>
    get_B_copy (const LocalOrdinal nrows_and_ncols,
                const Scalar B[], // HOST MEMORY
                const LocalOrdinal ldb) const
    {
      auto B_copy =
        Impl::get_contiguous_device_mat_view (matrixStorageB_,
                                              nrows_and_ncols,
                                              nrows_and_ncols);
      // Use copy_from_host, which knows how to avoid the
      // Kokkos::deep_copy noncontiguity problem.
      Scalar* B_copy_raw = reinterpret_cast<Scalar*> (B_copy.data ());
      const LocalOrdinal B_copy_stride (B_copy.extent (1));
      MatView<LocalOrdinal, Scalar> B_copy_matview
        (nrows_and_ncols, nrows_and_ncols, B_copy_raw, B_copy_stride);
      MatView<LocalOrdinal, const Scalar> B_matview
        (nrows_and_ncols, nrows_and_ncols, B, ldb);
      this->copy_from_host (B_copy_matview, B_matview);
      return B_copy;
    }

    void
    extract_R (const LocalOrdinal nrows,
               const LocalOrdinal ncols,
               const Scalar A[], // DEVICE POINTER
               const LocalOrdinal lda,
               Scalar R[], // HOST POINTER
               const LocalOrdinal ldr,
               const bool /* contiguous_cache_blocks */) const
    {
      using std::endl;
      const char prefix[] = "TSQR::CuSolverNodeTsqr::extract_R: ";

      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "Top of TSQR::CuSolverNodeTsqr::extract_R" );

      std::unique_ptr<std::ostringstream> err;
      bool threw = false;

      using Impl::get_device_mat_view;
      using a_view_type = decltype (get_device_mat_view<const Scalar>
                                    (nrows, ncols, A, lda));
      a_view_type A_view;
      try {
        A_view = get_device_mat_view<const Scalar>
          (nrows, ncols, A, lda);
      }
      TSQR_IMPL_CATCH( "get_device_mat_view of A threw: " );

      auto R_view =
        Impl::get_host_mat_view<Scalar> (ncols, ncols, R, ldr);

      try {
        // Fill R (including lower triangle) with zeros.
        //Kokkos::deep_copy (R_view, kokkos_value_type {});

        // The above code throws the following exception, even though
        // R_view is most definitely a host View:
        //
        // TSQR::CuSolverNodeTsqr::extract_R:
        // Kokkos::deep_copy(R_view, 0) threw an exception:
        // cudaDeviceSynchronize() error( cudaErrorIllegalAddress): an
        // illegal memory access was encountered
        // .../kokkos/core/src/Cuda/Kokkos_Cuda_Instance.cpp:120

        MatView<LocalOrdinal, Scalar> R_mv (ncols, ncols, R, ldr);
        deep_copy (R_mv, Scalar {});
      }
      TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_view, 0.0) threw: " );

      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::extract_R, "
                                       "after deep_copy(R_mv, 0.0)" );

      // Copy out the upper triangle of the R factor from A into R.
      //
      // The following (pseudo)code often does not work:
      //
      // auto A_view_top = subview(A_view, {0, ncols}, ALL());
      // Kokkos::deep_copy(R_view, A_view_top);
      //
      // Kokkos throws an exception, claiming "no available copy
      // mechanism."  This is probably because A_view is not packed.
      // This means that cudaMemcpy won't work, so Kokkos must execute
      // a kernel to copy the data.  However, that kernel must be able
      // to access both Views.  In this case, it (thinks it) can't,
      // because R_view is a HostSpace View and A_view_top is a device
      // View (even though it may be a CudaUVMSpace View).

      using Kokkos::ALL;
      using Kokkos::subview;
      using LO = LocalOrdinal;
      const std::pair<LO, LO> rowRange (0, ncols);
      auto A_view_top = subview (A_view, rowRange, ALL ());

      if (size_t (A_view_top.stride (1)) == size_t (A_view_top.extent (0))) {
        try {
          Kokkos::deep_copy (R_view, A_view_top);
        }
        TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_view, A_view_top) "
                         "for contiguous A_view_top threw: ");
        TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::extract_R, "
                                         "after attempting "
                                         "Kokkos::deep_copy(R_view, A_view_top) "
                                         "with contiguous A_view_top" );
      }
      else { // A_view_top is NOT contiguous
        // Packed device version of R.
        Impl::device_mat_view_type<kokkos_value_type> R_contig_d;
        try {
          using Impl::get_contiguous_device_mat_view;
          R_contig_d = get_contiguous_device_mat_view (matrixStorage_,
                                                       ncols, ncols);
        }
        TSQR_IMPL_CATCH( "R_contig_d = get_contiguous_device_mat_view threw: " );

        TEUCHOS_ASSERT( size_t (R_contig_d.extent (0)) == size_t (ncols) );
        TEUCHOS_ASSERT( size_t (R_contig_d.extent (1)) == size_t (ncols) );
        TEUCHOS_ASSERT( size_t (R_contig_d.stride (1)) == size_t (ncols) );

        try {
          Kokkos::deep_copy (R_contig_d, A_view_top);
        }
        TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_contig_d, A_view_top) threw: ");

        if (R_view.extent (0) < R_view.stride (1)) {
          // R_view is not contiguous, so we can't deep_copy directly
          // from R_contig_d (device View) to R_view (host View).  We
          // need an intermediate contiguous host View, R_contig_h.
          auto R_contig_h =
            Impl::get_contiguous_host_mat_view (hostMatrixStorage_,
                                                ncols, ncols);
          TEUCHOS_ASSERT( size_t (R_contig_h.extent (0)) == size_t (ncols) );
          TEUCHOS_ASSERT( size_t (R_contig_h.extent (1)) == size_t (ncols) );
          TEUCHOS_ASSERT( size_t (R_contig_h.stride (1)) == size_t (ncols) );
          try {
            Kokkos::deep_copy (R_contig_h, R_contig_d);
          }
          TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_contig_h, R_contig_d) threw: ");
          try {
            Kokkos::deep_copy (R_view, R_contig_h);
          }
          TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_view, R_contig_h) threw: ");
        }
        else { // R_view is contiguous, so we can deep_copy directly
          try {
            Kokkos::deep_copy (R_view, R_contig_d);
          }
          TSQR_IMPL_CATCH( "Kokkos::deep_copy(R_view, R_contig_d) threw: ");
        }
      }

      try {
        for (LO j = 0; j < ncols; ++j) {
          auto R_j = subview (R_view, ALL (), j);
          for (LO i = j + LO(1); i < LO (R_j.extent(0)); ++i) {
            R_j(i) = kokkos_value_type {};
          }
        }
      }
      TSQR_IMPL_CATCH( "Filling lower triangle of R_view with zeros threw: ");
    }

  public:
    Teuchos::RCP<factor_output_type>
    factor (const LocalOrdinal nrows,
            const LocalOrdinal ncols,
            Scalar A[],
            const LocalOrdinal lda,
            Scalar R[],
            const LocalOrdinal ldr,
            const bool contigCacheBlocks) const override
    {
      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::factor (top)" );

      // It's a common case to call factor() again and again with the
      // same pointers.  In that case, it's wasteful for us to
      // allocate a new tau array each time, especially since most
      // users want explicit Q anyway (and thus will never see tau).
      auto tau = get_tau (ncols);
      // FIXME (mfh 11 Dec 2019) TSQR::Impl::CuBlas takes
      // std::complex, but Kokkos::View stores Kokkos::complex.  We're
      // assuming they have the same alignment here, but all of Tpetra
      // assumes that.
      Scalar* tau_raw = reinterpret_cast<Scalar*> (tau.data ());
      auto work = get_work_for_geqrf (nrows, ncols, A, lda);
      Scalar* work_raw = reinterpret_cast<Scalar*> (work.data ());
      const int lwork (work.extent (0));
      auto info = get_info ();

      TSQR::Impl::CuSolver<Scalar> solver (info.data ());
      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::factor, "
                                       "before solver.compute_QR" );
      try {
        solver.compute_QR (nrows, ncols, A, lda, tau_raw,
                           work_raw, lwork);
      }
      catch (std::exception& e) {
        std::ostringstream err;
        err << "TSQR::CuSolverNodeTsqr::factor: CuSolver::compute_QR "
          "threw an exception: " << std::endl << e.what ();
        throw std::runtime_error (err.str ());
      }
      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::factor, "
                                       "after solver.compute_QR, "
                                       "before extract_R" );
      try {
        this->extract_R (nrows, ncols, A, lda, R, ldr,
                         contigCacheBlocks);
      }
      catch (std::exception& e) {
        std::ostringstream err;
        err << "TSQR::CuSolverNodeTsqr::factor: extract_R "
          "threw an exception: " << std::endl << e.what ();
        throw std::runtime_error (err.str ());
      }

      TSQR_IMPL_CHECK_LAST_CUDA_ERROR( "TSQR::CuSolverNodeTsqr::factor, "
                                       "after extract_R" );
      return Teuchos::rcp (new my_factor_output_type (tau, info));
    }

  private:
    const my_factor_output_type&
    get_my_factor_output (const factor_output_type& factor_output) const
    {
      const char prefix[] = "TSQR::CuSolverNodeTsqr: ";

      const my_factor_output_type* output_ptr =
        dynamic_cast<const my_factor_output_type*> (&factor_output);
      if (output_ptr == nullptr) {
        const std::string this_name = Teuchos::typeName (*this);
        const std::string factor_output_type_name =
          Teuchos::TypeNameTraits<my_factor_output_type>::name ();
        const std::string dynamic_type_name =
          Teuchos::demangleName (typeid (factor_output).name ());
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, prefix << "Input "
           "factor_output_type object was not created by the same "
           "type of CuSolverNodeTsqr object as this one.  This "
           "object has type " << this_name << " and its subclass of "
           "factor_output_type has type " << factor_output_type_name
           << ", but the input factor_output_type object has dynamic "
           "type " << dynamic_type_name << ".");
      }
      return *output_ptr;
    }

  public:
    void
    apply (const ApplyType& apply_type,
           const LocalOrdinal nrows,
           const LocalOrdinal ncols_Q,
           const Scalar Q[],
           const LocalOrdinal ldq,
           const factor_output_type& factor_output,
           const LocalOrdinal ncols_C,
           Scalar C[],
           const LocalOrdinal ldc,
           const bool contigCacheBlocks) const override
    {
      const char prefix[] = "TSQR::CuSolverNodeTsqr::apply: ";

      // Quick exit and error tests
      if (ncols_Q == 0 || ncols_C == 0 || nrows == 0) {
        return;
      }
      else if (ldc < nrows) {
        std::ostringstream os;
        os << prefix << "ldc (= " << ldc << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }
      else if (ldq < nrows) {
        std::ostringstream os;
        os << prefix << "ldq (= " << ldq << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }

      const char side = 'L';
      const char trans = apply_type.toString ()[0];
      auto tau = get_my_factor_output (factor_output).tau ();
      // FIXME (mfh 11 Dec 2019) TSQR::Impl::CuBlas takes
      // std::complex, but Kokkos::View stores Kokkos::complex.  We're
      // assuming they have the same alignment here, but all of Tpetra
      // assumes that.
      const Scalar* tau_raw =
        reinterpret_cast<const Scalar*> (tau.data ());
      auto work =
        get_work_for_apply_Q_factor (apply_type,
                                     nrows, ncols_C, ncols_Q,
                                     Q, ldq, tau_raw, C, ldc);
      Scalar* work_raw = reinterpret_cast<Scalar*> (work.data ());
      const int lwork (work.extent (0));
      auto info = get_info ();

      TSQR::Impl::CuSolver<Scalar> solver (info.data ());
      solver.apply_Q_factor (side, trans,
                             nrows, ncols_C, ncols_Q,
                             Q, ldq, tau_raw,
                             C, ldc,
                             work_raw, lwork);
    }

    /// \brief Copy from a host matrix, to "native" NodeTsqr device
    ///   storage.
    virtual void
    copy_from_host (const MatView<LocalOrdinal, Scalar>& C_dev,
                    const MatView<LocalOrdinal, const Scalar>& C_host) const
    {
      const char prefix[] =
        "TSQR::CuSolverNodeTsqr::copy_from_host: ";

      const size_t nrows (C_dev.extent (0));
      const size_t ncols (C_dev.extent (1));
      TEUCHOS_ASSERT( nrows == size_t (C_host.extent (0)) );
      TEUCHOS_ASSERT( ncols == size_t (C_host.extent (1)) );

      auto C_dev_view = Impl::get_device_mat_view<Scalar>
        (nrows, ncols, C_dev.data (), C_dev.stride (1));
      auto C_host_view = Impl::get_host_mat_view<const Scalar>
        (nrows, ncols, C_host.data (), C_host.stride (1));

      // NOTE (mfh 17 Dec 2019) If C_host is contiguous, that is, if
      // C_host.stride(1) == C_host.extent(0), then we can
      // Kokkos::deep_copy directly.  Otherwise, Kokkos::deep_copy
      // will throw an exception, claiming "no available copy
      // mechanism."  This is because cudaMemcpy won't work, so Kokkos
      // must execute a kernel to copy the data.  (Kokkos doesn't seem
      // to exploit any of the various 2-D or 3-D array copying
      // functions that CUDA provides.)  That kernel must be able to
      // access both Views.  We deal with this with a fall-back path
      // that uses temporary contiguous storage.

      if (C_dev_view.stride (1) == C_dev_view.extent (0) &&
          C_host_view.stride (1) == C_host_view.extent (0)) {
        // Both Views are contiguous.
        try {
          Kokkos::deep_copy (C_dev_view, C_host_view);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, prefix <<
             "Kokkos::deep_copy(C_dev_view, C_host_view) (both "
             "contiguous) threw: " << e.what ());
        }
      }
      else {
        // We need to make a contiguous copy of host storage.
        auto C_host_copy = Impl::get_contiguous_host_mat_view
          (hostMatrixStorage_, nrows, ncols);
        TEUCHOS_ASSERT( C_host_copy.stride (1) ==
                        C_host_copy.extent (0) );
        try {
          Kokkos::deep_copy (C_host_copy, C_host_view);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, prefix <<
             "Kokkos::deep_copy(C_host_copy, C_host_view) threw: "
             << e.what ());
        }

        if (C_dev_view.stride (1) == C_dev_view.extent (0)) {
          try {
            Kokkos::deep_copy (C_dev_view, C_host_copy);
          }
          catch (std::exception& e) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (true, std::runtime_error, prefix <<
               "Kokkos::deep_copy(C_dev_view, C_host_copy) threw: "
               << e.what ());
          }
        }
        else {
          auto C_dev_copy = Impl::get_contiguous_device_mat_view
            (matrixStorage_, nrows, ncols);
          try {
            Kokkos::deep_copy (C_dev_copy, C_host_copy);
          }
          catch (std::exception& e) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (true, std::runtime_error, prefix <<
               "Kokkos::deep_copy(C_dev_copy, C_host_copy) threw: "
               << e.what ());
          }
          try {
            Kokkos::deep_copy (C_dev_view, C_dev_copy);
          }
          catch (std::exception& e) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (true, std::runtime_error, prefix <<
               "Kokkos::deep_copy(C_dev_view, C_dev_copy) threw: "
               << e.what ());
          }
        }
      }
    }

    /// \brief Copy from "native" NodeTsqr device storage, to a packed
    ///   host matrix.
    Matrix<LocalOrdinal, Scalar>
    copy_to_host
      (const MatView<LocalOrdinal, Scalar>& C) const override
    {
      using LO = LocalOrdinal;
      const LO nrows (C.extent (0));
      const LO ncols (C.extent (1));
      const LO ldc (C.stride (1));
      auto C_dev =
        Impl::get_device_mat_view<const Scalar> (nrows, ncols,
                                                 C.data (), ldc);
      Matrix<LO, Scalar> C_copy (nrows, ncols);
      auto C_host = Impl::get_host_mat_view (C_copy.view ());

      // NOTE (mfh 17 Dec 2019) Directly calling
      // Kokkos::deep_copy(C_host, C_dev) may not necessarily work,
      // since C_dev need not be contiguous.  In that case, Kokkos
      // would throw an exception, claiming "no available copy
      // mechanism."  The work-around is to create a packed device
      // View, copy C_dev into it, then copy the packed View to
      // C_host.
      try {
        Kokkos::deep_copy (C_host, C_dev);
      }
      catch (std::exception& /* e */) {
        auto C_dev_copy =
          Impl::get_contiguous_device_mat_view (matrixStorage_,
                                                nrows, ncols);
        Kokkos::deep_copy (C_dev_copy, C_dev);
        Kokkos::deep_copy (C_host, C_dev_copy);
      }
      return C_copy;
    }

    /// \brief Fill C (DEVICE MEMORY) with the first C.extent(1)
    ///   columns of the identity matrix.  Assume that C has already
    ///   been pre-filled with zeros.
    void
    set_diagonal_entries_to_one
      (const MatView<LocalOrdinal, Scalar>& C) const override
    {
      auto C_view =
        Impl::get_device_mat_view (C.extent (0), C.extent (1),
                                   C.data (), C.stride (1));
      Impl::set_diagonal_entries_to_one (C_view);
    }

    void
    explicit_Q (const LocalOrdinal nrows,
                const LocalOrdinal ncols_Q,
                const Scalar Q[], // DEVICE MEMORY
                const LocalOrdinal ldq,
                const factor_output_type& factor_output,
                const LocalOrdinal ncols_C,
                Scalar C[], // DEVICE MEMORY
                const LocalOrdinal ldc,
                const bool contigCacheBlocks) const override
    {
      using Impl::get_device_mat_view;
      auto C_view = get_device_mat_view (nrows, ncols_C, C, ldc);
      using IST = Impl::non_const_kokkos_value_type<Scalar>;
      deep_copy (C_view, IST {});
      Impl::set_diagonal_entries_to_one (C_view);
      apply (ApplyType::NoTranspose,
             nrows, ncols_Q, Q, ldq, factor_output,
             ncols_C, C, ldc, contigCacheBlocks);
    }

    void
    Q_times_B (const LocalOrdinal nrows,
               const LocalOrdinal ncols,
               Scalar Q[], // DEVICE MEMORY
               const LocalOrdinal ldq,
               const Scalar B[], // HOST MEMORY
               const LocalOrdinal ldb,
               const bool /* contigCacheBlocks */) const override
    {
      // Take the easy exit if available.
      if (ncols == 0 || nrows == 0) {
        return;
      }

      // _GEMM doesn't permit the in/out matrix to alias either of the
      // two input matrices, so we must make a copy.
      auto Q_copy = get_Q_copy (nrows, ncols, Q, ldq);

      // We assume that B is in host memory, so we need to copy it to
      // device before we can use cuBLAS.
      auto B_copy = get_B_copy (ncols, B, ldb);

      constexpr Scalar ZERO {};
      constexpr Scalar ONE (1.0);

      using Impl::CuBlas;
      CuBlas<Scalar> blas;

      const char transa = 'N';
      const char transb = 'N';
      // FIXME (mfh 11 Dec 2019) TSQR::Impl::CuBlas takes
      // std::complex, but Kokkos::View stores Kokkos::complex.  We're
      // assuming they have the same alignment here, but all of Tpetra
      // assumes that.
      const Scalar* Q_copy_raw =
        reinterpret_cast<const Scalar*> (Q_copy.data ());
      const int Q_copy_stride (Q_copy.stride (1));
      blas.gemm (transa, transb, nrows, ncols, ncols,
                 ONE, Q_copy_raw, Q_copy_stride,
                 B, ldb, ZERO, Q, ldq);
    }

    void
    cache_block (const LocalOrdinal /* nrows */,
                 const LocalOrdinal /* ncols */,
                 Scalar /* A_out */ [],
                 const Scalar /*A_in */ [],
                 const LocalOrdinal /* lda_in */) const override
    {}

    void
    un_cache_block (const LocalOrdinal /* nrows */,
                    const LocalOrdinal /* ncols */,
                    Scalar /* A_out */[],
                    const LocalOrdinal /* lda_out */,
                    const Scalar /* A_in */ []) const override
    {}

    void
    fill_with_zeros (const LocalOrdinal nrows,
                     const LocalOrdinal ncols,
                     Scalar A[],
                     const LocalOrdinal lda,
                     const bool /* contigCacheBlocks */) const override
    {
      auto A_view = Impl::get_device_mat_view (nrows, ncols, A, lda);
      Kokkos::deep_copy (A_view, kokkos_value_type {});
    }

  private:
    mutable tau_type tau_;
    mutable work_type work_;
    mutable Impl::info_type info_;
    mutable Impl::device_vector_type<kokkos_value_type> matrixStorage_;
    mutable Impl::device_vector_type<kokkos_value_type> matrixStorageB_;
    mutable Impl::host_vector_type<kokkos_value_type> hostMatrixStorage_;
  };

} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
#endif // TSQR_CUSOLVERNODETSQR_HPP
