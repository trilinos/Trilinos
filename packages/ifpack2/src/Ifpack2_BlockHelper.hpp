// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKHELPER_IMPL_HPP
#define IFPACK2_BLOCKHELPER_IMPL_HPP

#include "Ifpack2_BlockHelper_Timers.hpp"

namespace Ifpack2 {

  namespace BlockHelperDetails {

    namespace KB = KokkosBatched;

    ///
    /// view decorators for unmanaged and const memory
    ///
    using do_not_initialize_tag = Kokkos::ViewAllocateWithoutInitializing;

    template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
    using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::is_unmanaged |
                                              MemoryTraitsType::is_random_access |
                                              flag>;

    template <typename ViewType>
    using Unmanaged = Kokkos::View<typename ViewType::data_type,
                                   typename ViewType::array_layout,
                                   typename ViewType::device_type,
                                  MemoryTraits<typename ViewType::memory_traits,Kokkos::Unmanaged> >;
    template <typename ViewType>
    using Atomic = Kokkos::View<typename ViewType::data_type,
                                typename ViewType::array_layout,
                                typename ViewType::device_type,
                                MemoryTraits<typename ViewType::memory_traits,Kokkos::Atomic> >;
    template <typename ViewType>
    using Const = Kokkos::View<typename ViewType::const_data_type,
                               typename ViewType::array_layout,
                               typename ViewType::device_type,
                               typename ViewType::memory_traits>;
    template <typename ViewType>
    using ConstUnmanaged = Const<Unmanaged<ViewType> >;

    template <typename ViewType>
    using AtomicUnmanaged = Atomic<Unmanaged<ViewType> >;

    template <typename ViewType>
    using Unmanaged = Kokkos::View<typename ViewType::data_type,
                                   typename ViewType::array_layout,
                                   typename ViewType::device_type,
                                   MemoryTraits<typename ViewType::memory_traits,Kokkos::Unmanaged> >;


    template <typename ViewType>
    using Scratch = Kokkos::View<typename ViewType::data_type,
                                 typename ViewType::array_layout,
                                 typename ViewType::execution_space::scratch_memory_space,
                                 MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;

    /// 
    /// tpetra little block index
    ///
    template<typename LayoutType> struct TpetraLittleBlock;
    template<> struct TpetraLittleBlock<Kokkos::LayoutLeft> {
      template<typename T> KOKKOS_INLINE_FUNCTION
      static T getFlatIndex(const T i, const T j, const T blksize) { return i+j*blksize; }
    };
    template<> struct TpetraLittleBlock<Kokkos::LayoutRight> {
      template<typename T> KOKKOS_INLINE_FUNCTION
      static T getFlatIndex(const T i, const T j, const T blksize) { return i*blksize+j; }
    };

    ///
    /// block tridiag scalar type
    ///
    template<typename T> struct BlockTridiagScalarType { typedef T type; };
#if defined(IFPACK2_BLOCKHELPER_USE_SMALL_SCALAR_FOR_BLOCKTRIDIAG)
    template<> struct BlockTridiagScalarType<double> { typedef float type; };
    //template<> struct SmallScalarType<Kokkos::complex<double> > { typedef Kokkos::complex<float> type; };
#endif

    ///
    /// cuda specialization
    ///
    template<typename T> struct is_cuda                 { enum : bool { value = false }; };
#if defined(KOKKOS_ENABLE_CUDA)
    template<> struct is_cuda<Kokkos::Cuda>             { enum : bool { value = true  }; };
#endif

    ///
    /// hip specialization
    ///
    template<typename T> struct is_hip                  { enum : bool { value = false }; };
#if defined(KOKKOS_ENABLE_HIP)
    template<> struct is_hip<Kokkos::HIP> { enum : bool { value = true  }; };
#endif

    ///
    /// sycl specialization
    ///
    template<typename T> struct is_sycl                  { enum : bool { value = false }; };
#if defined(KOKKOS_ENABLE_SYCL)
    template<> struct is_sycl<Kokkos::Experimental::SYCL> { enum : bool { value = true  }; };
#endif

    template<typename T> struct is_device                  { enum : bool { value = is_cuda<T>::value || is_hip<T>::value || is_sycl<T>::value }; };

    
    ///
    /// execution space instance
    ///
    template<typename T>
    struct ExecutionSpaceFactory {
      static void createInstance(T &exec_instance) {
        exec_instance = T();
      }
#if defined(KOKKOS_ENABLE_CUDA)
      static void createInstance(const cudaStream_t &s, T &exec_instance) {
        exec_instance = T();
      }
#endif
    };

#if defined(KOKKOS_ENABLE_CUDA)
    template<>
    struct ExecutionSpaceFactory<Kokkos::Cuda> {
      static void createInstance(Kokkos::Cuda &exec_instance) {
        exec_instance = Kokkos::Cuda();
      }
      static void createInstance(const cudaStream_t &s, Kokkos::Cuda &exec_instance) {
        exec_instance = Kokkos::Cuda(s);
      }      
    };
#endif

#if defined(KOKKOS_ENABLE_HIP)
    template<>
    struct ExecutionSpaceFactory<Kokkos::HIP> {
      static void createInstance(Kokkos::HIP &exec_instance) {
	exec_instance = Kokkos::HIP();
      }
    };
#endif

#if defined(KOKKOS_ENABLE_SYCL)
    template<>
    struct ExecutionSpaceFactory<Kokkos::Experimental::SYCL> {
      static void createInstance(Kokkos::Experimental::SYCL &exec_instance) {
	exec_instance = Kokkos::Experimental::SYCL();
      }
    };
#endif

#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKHELPER_ENABLE_PROFILE)
#define IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN \
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaProfilerStart());

#define IFPACK2_BLOCKHELPER_PROFILER_REGION_END \
    { KOKKOS_IMPL_CUDA_SAFE_CALL( cudaProfilerStop() ); }
#else
    /// later put vtune profiler region
#define IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN
#define IFPACK2_BLOCKHELPER_PROFILER_REGION_END
#endif

    
    ///
    /// utility functions
    ///
    template<typename CommPtrType>
    std::string get_msg_prefix (const CommPtrType &comm) {
      const auto rank = comm->getRank();
      const auto nranks = comm->getSize();
      std::stringstream ss;
      ss << "Rank " << rank << " of " << nranks << ": ";
      return ss.str();
    }

    ///
    /// custom multiple varilable reduce and scan
    ///
    template<typename T, int N>
    struct ArrayValueType {
      T v[N];
      KOKKOS_INLINE_FUNCTION
      ArrayValueType() {
        for (int i=0;i<N;++i)
          this->v[i] = 0;
      }
      KOKKOS_INLINE_FUNCTION
      ArrayValueType(const ArrayValueType &b) {
        for (int i=0;i<N;++i)
          this->v[i] = b.v[i];
      }
    };
    template<typename T, int N>
    static
    KOKKOS_INLINE_FUNCTION
    void
    operator+=(ArrayValueType<T,N> &a,
               const ArrayValueType<T,N> &b) {
      for (int i=0;i<N;++i)
        a.v[i] += b.v[i];
    }

    ///
    /// custom reducer functor for compile time array variable
    ///
    template<typename T, int N, typename ExecSpace>
    struct SumReducer {
      typedef SumReducer reducer;
      typedef ArrayValueType<T,N> value_type;
      typedef Kokkos::View<value_type,ExecSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
      value_type *value;

      KOKKOS_INLINE_FUNCTION
      SumReducer(value_type &val) : value(&val) {}

      KOKKOS_INLINE_FUNCTION
      void join(value_type &dst, value_type const &src) const {
        for (int i=0;i<N;++i)
          dst.v[i] += src.v[i];
      }
      KOKKOS_INLINE_FUNCTION
      void init(value_type &val) const {
        for (int i=0;i<N;++i)
          val.v[i] = Kokkos::reduction_identity<T>::sum();
      }
      KOKKOS_INLINE_FUNCTION
      value_type& reference() {
        return *value;
      }
      KOKKOS_INLINE_FUNCTION
      result_view_type view() const {
        return result_view_type(value);
      }
    };


    ///
    /// implementation typedefs
    ///
    template <typename MatrixType>
    struct ImplType {
      ///
      /// matrix type derived types
      ///
      typedef size_t size_type;
      typedef typename MatrixType::scalar_type scalar_type;
      typedef typename MatrixType::local_ordinal_type local_ordinal_type;
      typedef typename MatrixType::global_ordinal_type global_ordinal_type;
      typedef typename MatrixType::node_type node_type;

      ///
      /// kokkos arithmetic traits of scalar_type
      ///
      typedef typename Kokkos::Details::ArithTraits<scalar_type>::val_type impl_scalar_type;
      typedef typename Kokkos::ArithTraits<impl_scalar_type>::mag_type magnitude_type;

      typedef typename BlockTridiagScalarType<impl_scalar_type>::type btdm_scalar_type;
      typedef typename Kokkos::ArithTraits<btdm_scalar_type>::mag_type btdm_magnitude_type;

      ///
      /// default host execution space
      ///
      typedef Kokkos::DefaultHostExecutionSpace host_execution_space;

      ///
      /// tpetra types
      ///
      typedef typename node_type::device_type node_device_type;
      typedef typename node_device_type::execution_space node_execution_space;
      typedef typename node_device_type::memory_space node_memory_space;

#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKHELPER_USE_CUDA_SPACE)
      /// force to use cuda space instead uvm space
      typedef node_execution_space execution_space;
      typedef typename std::conditional<std::is_same<node_memory_space,Kokkos::CudaUVMSpace>::value,
                                        Kokkos::CudaSpace,
                                        node_memory_space>::type memory_space;
      typedef Kokkos::Device<execution_space,memory_space> device_type;
#else
      typedef node_device_type device_type;
      typedef node_execution_space execution_space;
      typedef node_memory_space memory_space;
#endif

      typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_multivector_type;
      typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> tpetra_map_type;
      typedef Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> tpetra_import_type;
      typedef Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_row_matrix_type;
      typedef Tpetra::CrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_crs_matrix_type;
      typedef Tpetra::CrsGraph<local_ordinal_type,global_ordinal_type,node_type> tpetra_crs_graph_type;
      typedef Tpetra::BlockCrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_crs_matrix_type;
      typedef typename tpetra_block_crs_matrix_type::little_block_type tpetra_block_access_view_type;
      typedef Tpetra::BlockMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_multivector_type;
      typedef typename tpetra_block_crs_matrix_type::crs_graph_type::local_graph_device_type local_crs_graph_type;

      ///
      /// simd vectorization
      ///
      template<typename T, int l> using Vector = KB::Vector<T,l>;
      template<typename T> using SIMD = KB::SIMD<T>;
      template<typename T, typename M> using DefaultVectorLength = KB::DefaultVectorLength<T,M>;
      template<typename T, typename M> using DefaultInternalVectorLength = KB::DefaultInternalVectorLength<T,M>;

      static constexpr int vector_length = DefaultVectorLength<btdm_scalar_type,memory_space>::value;
      static constexpr int internal_vector_length = DefaultInternalVectorLength<btdm_scalar_type,memory_space>::value;
      typedef Vector<SIMD<btdm_scalar_type>,vector_length> vector_type;
      typedef Vector<SIMD<btdm_scalar_type>,internal_vector_length> internal_vector_type;

      ///
      /// commonly used view types
      ///
      typedef Kokkos::View<size_type*,device_type> size_type_1d_view;
      typedef Kokkos::View<size_type**,device_type> size_type_2d_view;
      typedef Kokkos::View<local_ordinal_type*,device_type> local_ordinal_type_1d_view;
      typedef Kokkos::View<local_ordinal_type**,device_type> local_ordinal_type_2d_view;
      // tpetra block crs values
      typedef Kokkos::View<impl_scalar_type*,device_type> impl_scalar_type_1d_view;
      typedef Kokkos::View<impl_scalar_type*,node_device_type> impl_scalar_type_1d_view_tpetra;

      // tpetra multivector values (layout left): may need to change the typename more explicitly
      typedef Kokkos::View<impl_scalar_type**,Kokkos::LayoutLeft,device_type> impl_scalar_type_2d_view;
      typedef Kokkos::View<impl_scalar_type**,Kokkos::LayoutLeft,node_device_type> impl_scalar_type_2d_view_tpetra;

      // packed data always use layout right
      typedef Kokkos::View<vector_type*,device_type> vector_type_1d_view;
      typedef Kokkos::View<vector_type***,Kokkos::LayoutRight,device_type> vector_type_3d_view;
      typedef Kokkos::View<vector_type****,Kokkos::LayoutRight,device_type> vector_type_4d_view;
      typedef Kokkos::View<internal_vector_type***,Kokkos::LayoutRight,device_type> internal_vector_type_3d_view;
      typedef Kokkos::View<internal_vector_type****,Kokkos::LayoutRight,device_type> internal_vector_type_4d_view;
      typedef Kokkos::View<internal_vector_type*****,Kokkos::LayoutRight,device_type> internal_vector_type_5d_view;
      typedef Kokkos::View<btdm_scalar_type**,Kokkos::LayoutRight,device_type> btdm_scalar_type_2d_view;
      typedef Kokkos::View<btdm_scalar_type***,Kokkos::LayoutRight,device_type> btdm_scalar_type_3d_view;
      typedef Kokkos::View<btdm_scalar_type****,Kokkos::LayoutRight,device_type> btdm_scalar_type_4d_view;
      typedef Kokkos::View<btdm_scalar_type*****,Kokkos::LayoutRight,device_type> btdm_scalar_type_5d_view;
    };


    ///
    /// Manage the distributed part of the computation of residual norms.
    ///
    template<typename MatrixType>
    struct NormManager {
    public:
      using impl_type = ImplType<MatrixType>;
      using host_execution_space = typename impl_type::host_execution_space;
      using magnitude_type = typename impl_type::magnitude_type;

    private:
      bool collective_;
      int sweep_step_, sweep_step_upper_bound_;
#ifdef HAVE_IFPACK2_MPI
      MPI_Request mpi_request_;
      MPI_Comm comm_;
#endif
      magnitude_type work_[3];

    public:
      NormManager() = default;
      NormManager(const NormManager &b) = default;
      NormManager(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
        sweep_step_ = 1;
        sweep_step_upper_bound_ = 1;
        collective_ = comm->getSize() > 1;
        if (collective_) {
#ifdef HAVE_IFPACK2_MPI
          const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
          TEUCHOS_ASSERT( ! mpi_comm.is_null());
          comm_ = *mpi_comm->getRawMpiComm();
#endif
        }
        const magnitude_type zero(0), minus_one(-1);
        work_[0] = zero;
        work_[1] = zero;
        work_[2] = minus_one;
      }

      // Check the norm every sweep_step sweeps.
      void setCheckFrequency(const int sweep_step) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(sweep_step < 1, "sweep step must be >= 1");
        sweep_step_upper_bound_ = sweep_step;
        sweep_step_ = 1;
      }

      // Get the buffer into which to store rank-local squared norms.
      magnitude_type* getBuffer() { return &work_[0]; }

      // Call MPI_Iallreduce to find the global squared norms.
      void ireduce(const int sweep, const bool force = false) {
        if ( ! force && sweep % sweep_step_) return;

        IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::NormManager::Ireduce", Ireduce);

        work_[1] = work_[0];
#ifdef HAVE_IFPACK2_MPI
        auto send_data = &work_[1];
        auto recv_data = &work_[0];
        if (collective_) {
# if defined(IFPACK2_BLOCKTRIDICONTAINER_USE_MPI_3)
          MPI_Iallreduce(send_data, recv_data, 1,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_, &mpi_request_);
# else
          MPI_Allreduce (send_data, recv_data, 1,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_);
# endif
        }
#endif
      }

      // Check if the norm-based termination criterion is met. tol2 is the
      // tolerance squared. Sweep is the sweep index. If not every iteration is
      // being checked, this function immediately returns false. If a check must
      // be done at this iteration, it waits for the reduction triggered by
      // ireduce to complete, then checks the global norm against the tolerance.
      bool checkDone (const int sweep, const magnitude_type tol2, const bool force = false) {
        // early return
        if (sweep <= 0) return false;

        IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::NormManager::CheckDone", CheckDone);

        TEUCHOS_ASSERT(sweep >= 1);
        if ( ! force && (sweep - 1) % sweep_step_) return false;
        if (collective_) {
#ifdef HAVE_IFPACK2_MPI
# if defined(IFPACK2_BLOCKTRIDICONTAINER_USE_MPI_3)
          MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
# else
          // Do nothing.
# endif
#endif
        }
        bool r_val = false;
        if (sweep == 1) {
          work_[2] = work_[0];
        } else {
          r_val = (work_[0] < tol2*work_[2]);
        }

        // adjust sweep step
        const auto adjusted_sweep_step = 2*sweep_step_;
        if (adjusted_sweep_step < sweep_step_upper_bound_) {
          sweep_step_ = adjusted_sweep_step;
        } else {
          sweep_step_ = sweep_step_upper_bound_;
        }
        return r_val;
      }

      // After termination has occurred, finalize the norms for use in
      // get_norms{0,final}.
      void finalize () {
        work_[0] = std::sqrt(work_[0]); // after converged
        if (work_[2] >= 0)
          work_[2] = std::sqrt(work_[2]); // first norm
        // if work_[2] is minus one, then norm is not requested.
      }

      // Report norms to the caller.
      const magnitude_type getNorms0 () const { return work_[2]; }
      const magnitude_type getNormsFinal () const { return work_[0]; }
    };

    template<typename MatrixType>
    void reduceVector(const ConstUnmanaged<typename BlockHelperDetails::ImplType<MatrixType>::impl_scalar_type_1d_view> zz,
                      /* */ typename BlockHelperDetails::ImplType<MatrixType>::magnitude_type *vals) {
      IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::ReduceVector", ReduceVector);

      using impl_type = BlockHelperDetails::ImplType<MatrixType>;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
#if 0
      const auto norm2 = KokkosBlas::nrm1(zz);
#else
      impl_scalar_type norm2(0);
      Kokkos::parallel_reduce
        ("ReduceMultiVector::Device",
         Kokkos::RangePolicy<typename impl_type::execution_space>(0,zz.extent(0)),
         KOKKOS_LAMBDA(const local_ordinal_type &i, impl_scalar_type &update) {
          update += zz(i);
        }, norm2);
#endif
      vals[0] = Kokkos::ArithTraits<impl_scalar_type>::abs(norm2);

      IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename ImplType<MatrixType>::execution_space)
    }

  } // namespace BlockHelperDetails

} // namespace Ifpack2

#endif
