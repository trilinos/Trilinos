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
                               MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;
template <typename ViewType>
using Atomic = Kokkos::View<typename ViewType::data_type,
                            typename ViewType::array_layout,
                            typename ViewType::device_type,
                            MemoryTraits<typename ViewType::memory_traits, Kokkos::Atomic> >;
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
                               MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;

template <typename ViewType>
using Scratch = Kokkos::View<typename ViewType::data_type,
                             typename ViewType::array_layout,
                             typename ViewType::execution_space::scratch_memory_space,
                             MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;

///
/// tpetra little block index
///
template <typename LayoutType>
struct TpetraLittleBlock;
template <>
struct TpetraLittleBlock<Kokkos::LayoutLeft> {
  template <typename T>
  KOKKOS_INLINE_FUNCTION static T getFlatIndex(const T i, const T j, const T blksize) { return i + j * blksize; }
};
template <>
struct TpetraLittleBlock<Kokkos::LayoutRight> {
  template <typename T>
  KOKKOS_INLINE_FUNCTION static T getFlatIndex(const T i, const T j, const T blksize) { return i * blksize + j; }
};

///
/// block tridiag scalar type
///
template <typename T>
struct BlockTridiagScalarType { typedef T type; };
#if defined(IFPACK2_BLOCKHELPER_USE_SMALL_SCALAR_FOR_BLOCKTRIDIAG)
template <>
struct BlockTridiagScalarType<double> { typedef float type; };
// template<> struct SmallScalarType<Kokkos::complex<double> > { typedef Kokkos::complex<float> type; };
#endif

///
/// cuda specialization
///
template <typename T>
struct is_cuda {
  enum : bool { value = false };
};
#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct is_cuda<Kokkos::Cuda> {
  enum : bool { value = true };
};
#endif

///
/// hip specialization
///
template <typename T>
struct is_hip {
  enum : bool { value = false };
};
#if defined(KOKKOS_ENABLE_HIP)
template <>
struct is_hip<Kokkos::HIP> {
  enum : bool { value = true };
};
#endif

///
/// sycl specialization
///
template <typename T>
struct is_sycl {
  enum : bool { value = false };
};
#if defined(KOKKOS_ENABLE_SYCL)
template <>
struct is_sycl<Kokkos::Experimental::SYCL> {
  enum : bool { value = true };
};
#endif

template <typename T>
struct is_device {
  enum : bool { value = is_cuda<T>::value || is_hip<T>::value || is_sycl<T>::value };
};

///
/// execution space instance
///
template <typename T>
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
template <>
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
template <>
struct ExecutionSpaceFactory<Kokkos::HIP> {
  static void createInstance(Kokkos::HIP &exec_instance) {
    exec_instance = Kokkos::HIP();
  }
};
#endif

#if defined(KOKKOS_ENABLE_SYCL)
template <>
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
  { KOKKOS_IMPL_CUDA_SAFE_CALL(cudaProfilerStop()); }
#else
/// later put vtune profiler region
#define IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN
#define IFPACK2_BLOCKHELPER_PROFILER_REGION_END
#endif

///
/// utility functions
///
template <typename CommPtrType>
std::string get_msg_prefix(const CommPtrType &comm) {
  const auto rank   = comm->getRank();
  const auto nranks = comm->getSize();
  std::stringstream ss;
  ss << "Rank " << rank << " of " << nranks << ": ";
  return ss.str();
}

///
/// custom multiple varilable reduce and scan
///
template <typename T, int N>
struct ArrayValueType {
  T v[N];
  KOKKOS_INLINE_FUNCTION
  ArrayValueType() {
    for (int i = 0; i < N; ++i)
      this->v[i] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  ArrayValueType(const ArrayValueType &b) {
    for (int i = 0; i < N; ++i)
      this->v[i] = b.v[i];
  }
};
template <typename T, int N>
static KOKKOS_INLINE_FUNCTION void
operator+=(ArrayValueType<T, N> &a,
           const ArrayValueType<T, N> &b) {
  for (int i = 0; i < N; ++i)
    a.v[i] += b.v[i];
}

///
/// custom reducer functor for compile time array variable
///
template <typename T, int N, typename ExecSpace>
struct SumReducer {
  typedef SumReducer reducer;
  typedef ArrayValueType<T, N> value_type;
  typedef Kokkos::View<value_type, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
  value_type *value;

  KOKKOS_INLINE_FUNCTION
  SumReducer(value_type &val)
    : value(&val) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type &dst, value_type const &src) const {
    for (int i = 0; i < N; ++i)
      dst.v[i] += src.v[i];
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &val) const {
    for (int i = 0; i < N; ++i)
      val.v[i] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() {
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
  typedef MatrixType matrix_type;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  ///
  /// kokkos arithmetic traits of scalar_type
  ///
  typedef typename KokkosKernels::ArithTraits<scalar_type>::val_type impl_scalar_type;
  typedef typename KokkosKernels::ArithTraits<impl_scalar_type>::mag_type magnitude_type;

  typedef typename BlockTridiagScalarType<impl_scalar_type>::type btdm_scalar_type;
  typedef typename KokkosKernels::ArithTraits<btdm_scalar_type>::mag_type btdm_magnitude_type;

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
  typedef typename std::conditional<std::is_same<node_memory_space, Kokkos::CudaUVMSpace>::value,
                                    Kokkos::CudaSpace,
                                    node_memory_space>::type memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;
#else
  typedef node_device_type device_type;
  typedef node_execution_space execution_space;
  typedef node_memory_space memory_space;
#endif

  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> tpetra_multivector_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> tpetra_map_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> tpetra_import_type;
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> tpetra_row_matrix_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> tpetra_crs_matrix_type;
  typedef Tpetra::CrsGraph<local_ordinal_type, global_ordinal_type, node_type> tpetra_crs_graph_type;
  typedef Tpetra::BlockCrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> tpetra_block_crs_matrix_type;
  typedef typename tpetra_block_crs_matrix_type::little_block_type tpetra_block_access_view_type;
  typedef Tpetra::BlockMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> tpetra_block_multivector_type;
  typedef typename tpetra_block_crs_matrix_type::crs_graph_type::local_graph_device_type local_crs_graph_type;

  ///
  /// simd vectorization
  ///
  template <typename T, int l>
  using Vector = KB::Vector<T, l>;
  template <typename T>
  using SIMD = KB::SIMD<T>;
  template <typename T, typename M>
  using DefaultVectorLength = KB::DefaultVectorLength<T, M>;
  template <typename T, typename M>
  using DefaultInternalVectorLength = KB::DefaultInternalVectorLength<T, M>;

  static constexpr int vector_length          = DefaultVectorLength<btdm_scalar_type, memory_space>::value;
  static constexpr int internal_vector_length = DefaultInternalVectorLength<btdm_scalar_type, memory_space>::value;
  static constexpr int half_vector_length     = (vector_length > 1) ? (vector_length / 2) : 1;
  typedef Vector<SIMD<btdm_scalar_type>, vector_length> vector_type;
  typedef Vector<SIMD<btdm_scalar_type>, internal_vector_length> internal_vector_type;

  ///
  /// commonly used view types
  ///
  typedef Kokkos::View<size_type *, device_type> size_type_1d_view;
  typedef Kokkos::View<size_type **, device_type> size_type_2d_view;
  typedef Kokkos::View<int64_t ***, Kokkos::LayoutRight, device_type> i64_3d_view;
  typedef Kokkos::View<local_ordinal_type *, device_type> local_ordinal_type_1d_view;
  typedef Kokkos::View<local_ordinal_type **, device_type> local_ordinal_type_2d_view;
  // tpetra block crs values
  typedef Kokkos::View<impl_scalar_type *, device_type> impl_scalar_type_1d_view;
  typedef Kokkos::View<impl_scalar_type *, node_device_type> impl_scalar_type_1d_view_tpetra;

  // tpetra multivector values (layout left): may need to change the typename more explicitly
  typedef Kokkos::View<impl_scalar_type **, Kokkos::LayoutLeft, device_type> impl_scalar_type_2d_view;
  typedef Kokkos::View<impl_scalar_type **, Kokkos::LayoutLeft, node_device_type> impl_scalar_type_2d_view_tpetra;
  typedef Kokkos::View<const impl_scalar_type **, Kokkos::LayoutLeft, node_device_type> const_impl_scalar_type_2d_view_tpetra;

  // packed data always use layout right
  typedef Kokkos::View<vector_type *, device_type> vector_type_1d_view;
  typedef Kokkos::View<vector_type ***, Kokkos::LayoutRight, device_type> vector_type_3d_view;
  typedef Kokkos::View<vector_type ****, Kokkos::LayoutRight, device_type> vector_type_4d_view;
  typedef Kokkos::View<internal_vector_type ***, Kokkos::LayoutRight, device_type> internal_vector_type_3d_view;
  typedef Kokkos::View<internal_vector_type ****, Kokkos::LayoutRight, device_type> internal_vector_type_4d_view;
  typedef Kokkos::View<internal_vector_type *****, Kokkos::LayoutRight, device_type> internal_vector_type_5d_view;
  typedef Kokkos::View<btdm_scalar_type **, Kokkos::LayoutRight, device_type> btdm_scalar_type_2d_view;
  typedef Kokkos::View<btdm_scalar_type ***, Kokkos::LayoutRight, device_type> btdm_scalar_type_3d_view;
  typedef Kokkos::View<btdm_scalar_type ****, Kokkos::LayoutRight, device_type> btdm_scalar_type_4d_view;
  typedef Kokkos::View<btdm_scalar_type *****, Kokkos::LayoutRight, device_type> btdm_scalar_type_5d_view;

  // maximum supported block size for BTD/BJacobi solvers (some codepaths may allow larger)
  enum : int { max_blocksize = 32 };
};

///
/// A - Tridiags(A), i.e., R in the splitting A = D + R.
///
template <typename MatrixType>
struct AmD {
  using impl_type                       = BlockHelperDetails::ImplType<MatrixType>;
  using local_ordinal_type_1d_view      = typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view               = typename impl_type::size_type_1d_view;
  using i64_3d_view                     = typename impl_type::i64_3d_view;
  using impl_scalar_type_1d_view_tpetra = Unmanaged<typename impl_type::impl_scalar_type_1d_view_tpetra>;
  // rowptr points to the start of each row of A_colindsub.
  size_type_1d_view rowptr, rowptr_remote;
  // Indices into A's rows giving the blocks to extract. rowptr(i) points to
  // the i'th row. Thus, g.entries(A_colindsub(rowptr(row) : rowptr(row+1))),
  // where g is A's graph, are the columns AmD uses. If seq_method_, then
  // A_colindsub contains all the LIDs and A_colindsub_remote is empty. If !
  // seq_method_, then A_colindsub contains owned LIDs and A_colindsub_remote
  // contains the remote ones.
  local_ordinal_type_1d_view A_colindsub, A_colindsub_remote;
  // Precomputed direct offsets to A,x blocks, for owned entries (OverlapTag case) or all entries (AsyncTag case)
  i64_3d_view A_x_offsets;
  // Precomputed direct offsets to A,x blocks, for non-owned entries (OverlapTag case). For AsyncTag case this is left empty.
  i64_3d_view A_x_offsets_remote;

  // Currently always true.
  bool is_tpetra_block_crs;

  // If is_tpetra_block_crs, then this is a pointer to A_'s value data.
  impl_scalar_type_1d_view_tpetra tpetra_values;

  AmD()             = default;
  AmD(const AmD &b) = default;
};

template <typename MatrixType>
struct PartInterface {
  using local_ordinal_type         = typename BlockHelperDetails::ImplType<MatrixType>::local_ordinal_type;
  using local_ordinal_type_1d_view = typename BlockHelperDetails::ImplType<MatrixType>::local_ordinal_type_1d_view;
  using local_ordinal_type_2d_view = typename BlockHelperDetails::ImplType<MatrixType>::local_ordinal_type_2d_view;

  PartInterface()                       = default;
  PartInterface(const PartInterface &b) = default;

  // Some terms:
  //   The matrix A is split as A = D + R, where D is the matrix of tridiag
  // blocks and R is the remainder.
  //   A part is roughly a synonym for a tridiag. The distinction is that a part
  // is the set of rows belonging to one tridiag and, equivalently, the off-diag
  // rows in R associated with that tridiag. In contrast, the term tridiag is
  // used to refer specifically to tridiag data, such as the pointer into the
  // tridiag data array.
  //   Local (lcl) row are the LIDs. lclrow lists the LIDs belonging to each
  // tridiag, and partptr points to the beginning of each tridiag. This is the
  // LID space.
  //   Row index (idx) is the ordinal in the tridiag ordering. lclrow is indexed
  // by this ordinal. This is the 'index' space.
  //   A flat index is the mathematical index into an array. A pack index
  // accounts for SIMD packing.

  // Local row LIDs. Permutation from caller's index space to tridiag index
  // space.
  local_ordinal_type_1d_view lclrow;
  // partptr_ is the pointer array into lclrow_.
  local_ordinal_type_1d_view partptr;  // np+1
  local_ordinal_type_2d_view partptr_sub;
  local_ordinal_type_1d_view partptr_schur;
  // packptr_(i), for i the pack index, indexes partptr_. partptr_(packptr_(i))
  // is the start of the i'th pack.
  local_ordinal_type_1d_view packptr;  // npack+1
  local_ordinal_type_1d_view packptr_sub;
  local_ordinal_type_1d_view packindices_sub;
  local_ordinal_type_2d_view packindices_schur;
  // part2rowidx0_(i) is the flat row index of the start of the i'th part. It's
  // an alias of partptr_ in the case of no overlap.
  local_ordinal_type_1d_view part2rowidx0;  // np+1
  local_ordinal_type_1d_view part2rowidx0_sub;
  // part2packrowidx0_(i) is the packed row index. If vector_length is 1, then
  // it's the same as part2rowidx0_; if it's > 1, then the value is combined
  // with i % vector_length to get the location in the packed data.
  local_ordinal_type_1d_view part2packrowidx0;  // np+1
  local_ordinal_type_2d_view part2packrowidx0_sub;
  local_ordinal_type part2packrowidx0_back;  // So we don't need to grab the array from the GPU.
  // rowidx2part_ maps the row index to the part index.
  local_ordinal_type_1d_view rowidx2part;  // nr
  local_ordinal_type_1d_view rowidx2part_sub;
  // True if lcl{row|col} is at most a constant away from row{idx|col}. In
  // practice, this knowledge is not particularly useful, as packing for batched
  // processing is done at the same time as the permutation from LID to index
  // space. But it's easy to detect, so it's recorded in case an optimization
  // can be made based on it.
  bool row_contiguous;

  local_ordinal_type max_partsz;
  local_ordinal_type max_subpartsz;
  local_ordinal_type n_subparts_per_part;
  local_ordinal_type nparts;
};

///
/// Manage the distributed part of the computation of residual norms.
///
template <typename MatrixType>
struct NormManager {
 public:
  using impl_type            = ImplType<MatrixType>;
  using host_execution_space = typename impl_type::host_execution_space;
  using magnitude_type       = typename impl_type::magnitude_type;

 private:
  bool collective_;
  int sweep_step_, sweep_step_upper_bound_;
#ifdef HAVE_IFPACK2_MPI
  MPI_Request mpi_request_;
  MPI_Comm comm_;
#endif
  magnitude_type work_[3];

 public:
  NormManager()                     = default;
  NormManager(const NormManager &b) = default;
  NormManager(const Teuchos::RCP<const Teuchos::Comm<int> > &comm) {
    sweep_step_             = 1;
    sweep_step_upper_bound_ = 1;
    collective_             = comm->getSize() > 1;
    if (collective_) {
#ifdef HAVE_IFPACK2_MPI
      const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
      TEUCHOS_ASSERT(!mpi_comm.is_null());
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
    sweep_step_             = 1;
  }

  // Get the buffer into which to store rank-local squared norms.
  magnitude_type *getBuffer() { return &work_[0]; }

  // Call MPI_Iallreduce to find the global squared norms.
  void ireduce(const int sweep, const bool force = false) {
    if (!force && sweep % sweep_step_) return;

    IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::NormManager::Ireduce", Ireduce);

    work_[1] = work_[0];
#ifdef HAVE_IFPACK2_MPI
    auto send_data = &work_[1];
    auto recv_data = &work_[0];
    if (collective_) {
#if defined(IFPACK2_BLOCKTRIDICONTAINER_USE_MPI_3)
      MPI_Iallreduce(send_data, recv_data, 1,
                     Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                     MPI_SUM, comm_, &mpi_request_);
#else
      MPI_Allreduce(send_data, recv_data, 1,
                    Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                    MPI_SUM, comm_);
#endif
    }
#endif
  }

  // Check if the norm-based termination criterion is met. tol2 is the
  // tolerance squared. Sweep is the sweep index. If not every iteration is
  // being checked, this function immediately returns false. If a check must
  // be done at this iteration, it waits for the reduction triggered by
  // ireduce to complete, then checks the global norm against the tolerance.
  bool checkDone(const int sweep, const magnitude_type tol2, const bool force = false) {
    // early return
    if (sweep <= 0) return false;

    IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::NormManager::CheckDone", CheckDone);

    TEUCHOS_ASSERT(sweep >= 1);
    if (!force && (sweep - 1) % sweep_step_) return false;
    if (collective_) {
#ifdef HAVE_IFPACK2_MPI
#if defined(IFPACK2_BLOCKTRIDICONTAINER_USE_MPI_3)
      MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
#else
      // Do nothing.
#endif
#endif
    }
    bool r_val = false;
    if (sweep == 1) {
      work_[2] = work_[0];
    } else {
      r_val = (work_[0] < tol2 * work_[2]);
    }

    // adjust sweep step
    const auto adjusted_sweep_step = 2 * sweep_step_;
    if (adjusted_sweep_step < sweep_step_upper_bound_) {
      sweep_step_ = adjusted_sweep_step;
    } else {
      sweep_step_ = sweep_step_upper_bound_;
    }
    return r_val;
  }

  // After termination has occurred, finalize the norms for use in
  // get_norms{0,final}.
  void finalize() {
    work_[0] = std::sqrt(work_[0]);  // after converged
    if (work_[2] >= 0)
      work_[2] = std::sqrt(work_[2]);  // first norm
    // if work_[2] is minus one, then norm is not requested.
  }

  // Report norms to the caller.
  const magnitude_type getNorms0() const { return work_[2]; }
  const magnitude_type getNormsFinal() const { return work_[0]; }
};

template <typename MatrixType>
void reduceVector(const ConstUnmanaged<typename BlockHelperDetails::ImplType<MatrixType>::impl_scalar_type_1d_view> zz,
                  /* */ typename BlockHelperDetails::ImplType<MatrixType>::magnitude_type *vals) {
  IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
  IFPACK2_BLOCKHELPER_TIMER("BlockTriDi::ReduceVector", ReduceVector);

  using impl_type          = BlockHelperDetails::ImplType<MatrixType>;
  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
#if 0
      const auto norm2 = KokkosBlas::nrm1(zz);
#else
  impl_scalar_type norm2(0);
  Kokkos::parallel_reduce(
      "ReduceMultiVector::Device",
      Kokkos::RangePolicy<typename impl_type::execution_space>(0, zz.extent(0)),
      KOKKOS_LAMBDA(const local_ordinal_type &i, impl_scalar_type &update) {
        update += zz(i);
      },
      norm2);
#endif
  vals[0] = KokkosKernels::ArithTraits<impl_scalar_type>::abs(norm2);

  IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
  IFPACK2_BLOCKHELPER_TIMER_FENCE(typename ImplType<MatrixType>::execution_space)
}

}  // namespace BlockHelperDetails

}  // namespace Ifpack2

#endif
