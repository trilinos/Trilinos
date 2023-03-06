#ifndef TPETRA_DETAILS_SPMV_HPP
#define TPETRA_DETAILS_SPMV_HPP

#include "KokkosBlas1_scal.hpp"

#include "Tpetra_Details_debug_cwp.hpp"

namespace Tpetra {
namespace Details {

/*! \brief Construct a view of on-rank entries of a row of a matrix
    \param A The matrix
    \param aOff The off-rank offsets of A
    \param i The row of a to produce a view of
    \tparam CrsMatrix the type of A
    \tparam RowOffsetView the type of aOff
*/
template<typename CrsMatrix,typename RowOffsetView>
struct OnRankRowViewer {
  KOKKOS_INLINE_FUNCTION static KokkosSparse::SparseRowViewConst<CrsMatrix> view(
      const CrsMatrix &A,
      const RowOffsetView &aOff,
      const typename CrsMatrix::ordinal_type i
  ) {
    if (size_t(i) >= A.graph.row_map.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.graph.row_map.extent(0) %d i\n",
      __FILE__, __LINE__, int(A.graph.row_map.extent(0)), int(i));
    }

    const typename CrsMatrix::size_type start = A.graph.row_map(i);
    const typename CrsMatrix::ordinal_type stride = 1;

    if (size_t(i) >= aOff.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d aOff.extent(0) %d i\n",
      __FILE__, __LINE__, int(aOff.extent(0)), int(i));
    }
    if (start > aOff(i)) {
      CWP_PRINTF("%s:%d BAD SUB %d aOff(i) - %d start\n",
      __FILE__, __LINE__, int(aOff(i)), int(start));
    }

    const typename CrsMatrix::ordinal_type length = aOff(i) - start;

    if (size_t(start) >= A.values.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.values.extent(0) %d start\n",
      __FILE__, __LINE__, int(A.values.extent(0)), int(start));
    }
    if (size_t(start) >= A.graph.entries.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.graph.entries.extent(0) %d start\n",
      __FILE__, __LINE__, int(A.graph.entries.extent(0)), int(start));
    }

    if (size_t(start+length) > A.values.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.values.extent(0) %d start+length\n",
      __FILE__, __LINE__, int(A.values.extent(0)), int(start+length));
    }
    if (size_t(start+length) > A.graph.entries.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.graph.entries.extent(0) %d start+length\n",
      __FILE__, __LINE__, int(A.graph.entries.extent(0)), int(start+length));
    }

    return KokkosSparse::SparseRowViewConst<CrsMatrix> (
      &A.values(start),
      &A.graph.entries(start),
      stride,
      length
    );
  }
};

/*! \brief Construct a view of on-rank entries of a row of a matrix
    \param A The matrix
    \param aOff The off-rank offsets of A
    \param i The row of a to produce a view of
    \tparam CrsMatrix the type of A
    \tparam RowOffsetView the type of aOff
*/
template<typename CrsMatrix, typename RowOffsetView>
struct OffRankRowViewer {
  KOKKOS_INLINE_FUNCTION static KokkosSparse::SparseRowViewConst<CrsMatrix> view(
      const CrsMatrix &A,
      const RowOffsetView &aOff,
      const typename CrsMatrix::ordinal_type i
  ) {

    if (size_t(i) >= aOff.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d aOff.extent(0) %d i\n",
      __FILE__, __LINE__, int(aOff.extent(0)), int(i));
    }
    if (size_t(i+1) >= A.graph.row_map.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.graph.row_map.extent(0) %d i+1\n",
      __FILE__, __LINE__, int(A.graph.row_map.extent(0)), int(i+1));
    }

    const typename CrsMatrix::size_type start = aOff(i);

    if (start > A.graph.row_map(i+1)) {
      CWP_PRINTF("%s:%d BAD SUB %d A.graph.row_map(i+1) - %d start\n",
      __FILE__, __LINE__, int(A.graph.row_map(i+1)), int(start));
    }


    const typename CrsMatrix::ordinal_type stride = 1;
    const typename CrsMatrix::ordinal_type length = A.graph.row_map(i+1) - start;

    if (size_t(start) >= A.values.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.values.extent(0) %d start\n",
      __FILE__, __LINE__, int(A.values.extent(0)), int(start));
    }
    if (size_t(start) >= A.graph.entries.extent(0)) {
      CWP_PRINTF("%s:%d BAD ACCESS %d A.graph.entries.extent(0) %d start\n",
      __FILE__, __LINE__, int(A.graph.entries.extent(0)), int(start));
    }

    return KokkosSparse::SparseRowViewConst<CrsMatrix> (
      &A.values(start),
      &A.graph.entries(start),
      stride,
      length
    );
  }
};

// KokkosKernels SPMV_Functor
template <
typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector,
typename OffsetDeviceViewType, 
typename RowViewer, // OnRankRowView or OffRankRowView
bool CONJ
>
struct SpmvFunctor {

  using y_value_type = typename YVector::non_const_value_type;
  using x_value_type = typename XVector::non_const_value_type;

  using value_type = typename AMatrix::non_const_value_type;
  using ordinal_type = typename AMatrix::non_const_ordinal_type; 
  using size_type = typename AMatrix::non_const_size_type; 

  using execution_space = typename AMatrix::execution_space;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::Details::ArithTraits<value_type>;

  const Alpha alpha_;
  AMatrix A_;
  XVector x_;
  const Beta beta_;
  YVector y_;
  OffsetDeviceViewType offRankOffsets_;
  const ordinal_type rowsPerTeam_;

  SpmvFunctor(const Alpha alpha, const AMatrix A, const XVector x,
              const Beta beta, const YVector y, const OffsetDeviceViewType &offRankOffsets,
              const int rowsPerTeam)
      : alpha_(alpha),
        A_(A),
        x_(x),
        beta_(beta),
        y_(y),
        offRankOffsets_(offRankOffsets),
        rowsPerTeam_(rowsPerTeam) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iRow) const {
    if (iRow >= A_.numRows()) {
      return;
    }

    const auto row = RowViewer::view(A_, offRankOffsets_, iRow);

    const ordinal_type row_length = static_cast<ordinal_type>(row.length);

    y_value_type sum = 0;
    for (ordinal_type iEntry = 0; iEntry < row_length; iEntry++) {
      const value_type val =
          CONJ ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
      sum += val * x_(row.colidx(iEntry));
    }

    sum *= alpha_;

    if (0 == beta_) {
      y_(iRow) = sum;
    } else {
      y_(iRow) = beta_ * y_(iRow) + sum;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& dev) const {
    if (dev.team_rank() == 0 && dev.league_rank() == 0) {
      CWP_PRINTF("%s:%d inside SpmvFunctor\n",
      __FILE__, __LINE__);
    }
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, 0, rowsPerTeam_),
        [&](const ordinal_type& loop) {
          const ordinal_type iRow =
              static_cast<ordinal_type>(dev.league_rank()) * rowsPerTeam_ +
              loop;
          if (iRow >= A_.numRows()) {
            return;
          }
          const auto row = RowViewer::view(A_, offRankOffsets_, iRow);
          CWP_PRINTF("%s:%d %d,%d got view\n", __FILE__, __LINE__, dev.league_rank(), dev.team_rank());
          const ordinal_type row_length = static_cast<ordinal_type>(row.length);

          y_value_type sum = 0;

          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(dev, row_length),
              [&](const ordinal_type& iEntry, y_value_type& lsum) {
                if (size_t(iEntry) >= row.length) {
                  CWP_PRINTF("%s:%d BAD ACCESS %d row.length %d iEntry\n",
                  __FILE__, __LINE__, int(row.length), int(iEntry));
                }
                const value_type val = CONJ ? ATV::conj(row.value(iEntry))
                                                : row.value(iEntry);

                const ordinal_type xi = row.colidx(iEntry);
                if (size_t(xi) >= x_.extent(0)) {
                  CWP_PRINTF("%s:%d BAD ACCESS %d x_.extent(0) %d xi\n",
                  __FILE__, __LINE__, int(x_.extent(0)), int(xi));
                }
                if (xi < 0) {
                  CWP_PRINTF("%s:%d BAD ACCESS %d x_.extent(0) %d xi\n",
                  __FILE__, __LINE__, int(x_.extent(0)), int(xi));
                }
                lsum += val * x_(xi);
              },
              sum);

          CWP_PRINTF("%s:%d %d,%d did pps\n", __FILE__, __LINE__, dev.league_rank(), dev.team_rank());

          Kokkos::single(Kokkos::PerThread(dev), [&]() {
            sum *= alpha_;

            CWP_PRINTF("%s:%d %d,%d did *=\n", __FILE__, __LINE__, dev.league_rank(), dev.team_rank());

            if (size_t(iRow) >= y_.extent(0)) {
              CWP_PRINTF("%s:%d BAD ACCESS %d y_.extent(0_) %d iRow\n", 
              __FILE__, __LINE__, int(y_.extent(0)), int(iRow));
              return;
            }
            if (iRow < 0) {
              CWP_PRINTF("%s:%d BAD ACCESS %d y_.extent(0_) %d iRow\n", 
              __FILE__, __LINE__, int(y_.extent(0)), int(iRow));
              return;
            }
            CWP_PRINTF("%s:%d %d,%d iRow=%d y_.extent(0)=%d\n", 
            __FILE__, __LINE__, int(dev.league_rank()), int(dev.team_rank()),
            int(iRow), int(y_.extent(0)));

            if (0 == beta_) {
              y_(iRow) = sum;
            } else {
              y_(iRow) = beta_ * y_(iRow) + sum;
            }
          });

          CWP_PRINTF("%s:%d %d,%d did loop\n", __FILE__, __LINE__, dev.league_rank(), dev.team_rank());
        });
  }


  struct LaunchParams {
    ordinal_type rowsPerThread;
    ordinal_type rowsPerTeam;
    int teamSize;
    int vectorLength;
  };

  template <class execution_space>
  static LaunchParams launch_parameters(ordinal_type numRows, int64_t nnz) {

    TEUCHOS_TEST_FOR_EXCEPTION(nnz <= int64_t(0), std::out_of_range, 
    "can't determine launch parameters for nnz=" << nnz);

    TEUCHOS_TEST_FOR_EXCEPTION(numRows <= ordinal_type(0), std::out_of_range, 
    "can't determine launch parameters for numRows=" << nnz);

    LaunchParams ret;

    int64_t nnzPerRow = nnz / numRows;
    if (nnzPerRow < 1) nnzPerRow = 1;

    int maxVectorLength = 1;
  #ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<execution_space, Kokkos::Cuda>::value)
      maxVectorLength = 32;
  #endif
  #ifdef KOKKOS_ENABLE_HIP
    if (std::is_same<execution_space, Kokkos::Experimental::HIP>::value)
      maxVectorLength = 64;
  #endif
    ret.vectorLength = 1;
    while (ret.vectorLength < maxVectorLength && ret.vectorLength * 6 < nnzPerRow) {
      ret.vectorLength *= 2;
    }

    // Determine rows per thread
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>())
      ret.rowsPerThread = 1;
    else {
      if (nnzPerRow < 20 && nnz > 5000000) {
        ret.rowsPerThread = 256;
      } else
        ret.rowsPerThread = 64;
    }

    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
      ret.teamSize = 256 / ret.vectorLength;
    } else {
      ret.teamSize = 1;
    }

    ret.rowsPerTeam = ret.rowsPerThread * ret.teamSize;

    // never happens?
    if (ret.rowsPerTeam < 0) {
      int64_t nnzPerTeam = 4096;
      int64_t conc         = execution_space::concurrency();
      while ((conc * nnzPerTeam * 4 > nnz) && (nnzPerTeam > 256))
        nnzPerTeam /= 2;
      ret.rowsPerTeam = (nnzPerTeam + nnzPerRow - 1) / nnzPerRow;
    }
    
    return ret;
  }
  
  // a kind of crummy scal implementation that we can put in an execution space instance
  struct ScalTag{};
  KOKKOS_INLINE_FUNCTION void operator()(const ScalTag&, const ordinal_type i) const {
    if (i >= ordinal_type(y_.extent(0))) {return;}
    if (0 == beta_) {
      y_(i) = 0;
    } else {
      y_(i) *= beta_;
    }
  }

  static void launch(const Alpha &alpha, 
    const AMatrix &A, 
    XVector &X, 
    const Beta &beta, 
    YVector &Y,
    const OffsetDeviceViewType &offRankOffsets,
    const execution_space &space) {

    CWP_CERR(__FILE__ << ":" << __LINE__
        << " A.numRows()=" << A.numRows()
        << " A.nnz()=" << A.nnz()
        << "\n");

    // if there's no matrix, still need to scale beta
    if (A.numRows() <= 0 || A.nnz() <= 0) {
      if (1 != beta) {
        SpmvFunctor op(alpha, A, X, beta, Y, offRankOffsets, /*unused*/ 0);
        Kokkos::parallel_for("SpmvFunctor scal",
                            Kokkos::RangePolicy<ScalTag, execution_space>(
                                space, 0, Y.extent(0)),
                            op);
      }
      return;
    }



    if constexpr (Details::Spaces::is_gpu_exec_space<execution_space>()) {

      LaunchParams params = launch_parameters<execution_space>(A.numRows(), A.nnz());

      CWP_CERR(__FILE__ << ":" << __LINE__ << ": rowsPerTeam=" << params.rowsPerTeam
                << " vectorLength=" << params.vectorLength
                << " teamSize=" << params.teamSize
                << "\n");

      const bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
      const bool use_static_schedule  = false;  // Forces the use of a static schedule
      if (params.rowsPerTeam < 1) {
        throw std::runtime_error("rowsPerTeam was < 1");
      }
      int64_t worksets = (Y.extent(0) + params.rowsPerTeam - 1) / params.rowsPerTeam;

      SpmvFunctor func(alpha, A, X, beta, Y, offRankOffsets, params.rowsPerTeam);

      if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
        Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
            policy(space, 1, 1);
        if (params.teamSize < 0)
          policy = Kokkos::TeamPolicy<execution_space,
                                      Kokkos::Schedule<Kokkos::Dynamic>>(
              space, worksets, Kokkos::AUTO, params.vectorLength);
        else
          policy = Kokkos::TeamPolicy<execution_space,
                                      Kokkos::Schedule<Kokkos::Dynamic>>(
              space, worksets, params.teamSize, params.vectorLength);
        Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Dynamic>", policy,
                            func);
      } else {
        Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>
            policy(space, 1, 1);
        if (params.teamSize < 0)
          policy =
              Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
                  space, worksets, Kokkos::AUTO, params.vectorLength);
        else
          policy =
              Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
                  space, worksets, params.teamSize, params.vectorLength);
        Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Static>", policy,
                            func);
      }
    } else {
      SpmvFunctor func(alpha, A, X, beta, Y, offRankOffsets, 1);

      const bool useDynamicSchedule = false;
      const bool useStaticSchedule = false;
      if (((A.nnz() > 10000000) || useDynamicSchedule) && !useStaticSchedule) {

        Kokkos::parallel_for(
            "KokkosSparse::spmv<NoTranspose,Dynamic>",
            Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
                space, 0, A.numRows()),
            func);
      } else {
        Kokkos::parallel_for(
            "KokkosSparse::spmv<NoTranspose,Static>",
            Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
                space, 0, A.numRows()),
            func);
      }
    }
  }
}; // SpmvFunctor

// cwp 21 Sep 2022
// A functor that does on-rank or off-rank part of a local Sparse-matrix multivector product
// KokkosKernels does not currently have a 4-array CSR
// SPMV_MV_LayoutLeft_Functor
template<
typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector,
typename OffsetDeviceViewType, 
typename RowViewer, // OnRankRowView or OffRankRowView
bool CONJ
> class SpmvMvFunctor {
public:
  using y_value_type = typename YVector::non_const_value_type;
  using x_value_type = typename XVector::non_const_value_type;
  using value_type = typename AMatrix::non_const_value_type;
  using ordinal_type = typename AMatrix::non_const_ordinal_type; 
  using size_type = typename AMatrix::non_const_size_type; 

  using execution_space = typename AMatrix::execution_space;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::Details::ArithTraits<value_type>;


  Alpha alpha_;
  AMatrix A_;
  XVector X_;
  Beta beta_;
  YVector Y_;
  OffsetDeviceViewType offRankOffsets_;
  ordinal_type rowsPerThread_;
  int vectorLength_;

public:
  SpmvMvFunctor(const Alpha &alpha, 
  const AMatrix &A, 
  XVector &X, 
  const Beta &beta, 
  YVector &Y,
  const OffsetDeviceViewType &offRankOffsets,
  ordinal_type rowsPerThread,
  int vectorLength) 
    : alpha_(alpha), A_(A), X_(X), beta_(beta), Y_(Y),
      offRankOffsets_(offRankOffsets),
      rowsPerThread_(rowsPerThread), vectorLength_(vectorLength) {}

  
  // a kind of crummy scal implementation that we can put in an execution space instance
  struct ScalTag{};
  KOKKOS_INLINE_FUNCTION void operator()(const ScalTag&, const ordinal_type i) const {
    if (i >= ordinal_type(Y_.extent(0))) {return;}
    for (ordinal_type k = 0; k < ordinal_type(Y_.extent(1)); ++k) {
      if (0 == beta_) {
        Y_(i, k) = 0;
      } else {
        Y_(i, k) *= beta_;
      }
    }
  }

  struct SimpleTag{};

  /* simplified version brought from Kokkos Kernels
  */
  template <int UNROLL>
  KOKKOS_INLINE_FUNCTION void strip_mine(const team_member& dev,
                                          const ordinal_type& iRow,
                                          const ordinal_type& kk) const {
    y_value_type sum[UNROLL];

    for (int k = 0; k < UNROLL; ++k) {
      sum[k] = Kokkos::Details::ArithTraits<y_value_type>::zero();
    }

    const auto row = RowViewer::view(A_, offRankOffsets_, iRow);

    // split loop indices 0..row.length over the vector lanes of the calling thread
    // each lane will have a subset of the partial products
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(dev, row.length), [&](ordinal_type iEntry) {
          const value_type val =
              CONJ ? Kokkos::Details::ArithTraits<value_type>::conj(
                              row.value(iEntry))
                        : row.value(iEntry);
          const ordinal_type ind = row.colidx(iEntry);
          for (int k = 0; k < UNROLL; ++k) {
            sum[k] += val * X_(ind, kk + k);
          }
        }
    );

#if 0
    for (int ii = 0; ii < UNROLL; ++ii) {
      int tr = dev.team_rank();
      int lr = dev.league_rank();
      printf("pre:         <%d,%d> sum(%d, %d) = %d\n",
          lr, tr, int(iRow), int(ii), int(sum[ii]));
    }
#endif

    for (int ii = 0; ii < UNROLL; ++ii) {

        // in CUDA, each lane is actually a sub-warp (thread) so each thread's sum is not correct
        // in serial, each lane is a true vector lane so sum[k] is already correct
        if (Details::Spaces::is_gpu_exec_space<execution_space>()) {
          // try to sum up the sum[ii] in each vector lane
          y_value_type sumt;
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(dev, vectorLength_),
              [&](ordinal_type, y_value_type& lsum) {
                lsum += sum[ii];
              },
              sumt);
          // sumt is broadcast to every lane of the thread

          // now every lane should have the complete product
          // of the row this thread was working on
          sum[ii] = sumt * alpha_;
        } else {
          sum[ii] *= alpha_;
        }
    }

#if 0
    for (int ii = 0; ii < UNROLL; ++ii) {
      int tr = dev.team_rank();
      int lr = dev.league_rank();
      printf("all-reduced: <%d,%d> sum(%d, %d) = %d\n",
          lr, tr, int(iRow), int(ii), int(sum[ii]));
    }
#endif

    // split 0..UNROLL over the vector lanes of the calling thread
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(dev, UNROLL), [&](ordinal_type k) {
          if (0 == beta_) {
            Y_(iRow, kk + k) = sum[k];
          } else {
            Y_(iRow, kk + k) = beta_ * Y_(iRow, kk + k) + sum[k];
          }
          
        }
    );
    
  }

  /* simplified version brought from Kokkos Kernels
  */
  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    for (ordinal_type loop = 0; loop < rowsPerThread_; ++loop) {
      const ordinal_type iRow =
          (dev.league_rank() * dev.team_size() + dev.team_rank()) *
              rowsPerThread_ + loop;
      if (iRow >= A_.numRows()) {
        return;
      }

      // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
      // needs to have the same type as n.
      ordinal_type kk = 0;

      const ordinal_type n = X_.extent(1);

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for (; kk + 8 <= n; kk += 8) strip_mine<8>(dev, iRow, kk);
      if (kk < n) {
        switch (n - kk) {
#else   // NOT a GPU
      if ((n > 16) && (n % 16 == 1)) {
        strip_mine<17>(dev, iRow, kk);
        kk += 17;
      }

      for (; kk + 16 <= n; kk += 16) {
        strip_mine<16>(dev, iRow, kk);
      }
      if (kk < n) {
        switch (n - kk) {
          case 15: strip_mine<15>(dev, iRow, kk); break;

          case 14: strip_mine<14>(dev, iRow, kk); break;

          case 13: strip_mine<13>(dev, iRow, kk); break;

          case 12: strip_mine<12>(dev, iRow, kk); break;

          case 11: strip_mine<11>(dev, iRow, kk); break;

          case 10: strip_mine<10>(dev, iRow, kk); break;

          case 9: strip_mine<9>(dev, iRow, kk); break;

          case 8: strip_mine<8>(dev, iRow, kk); break;
#endif  // if/else: __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__
          case 7: strip_mine<7>(dev, iRow, kk); break;

          case 6: strip_mine<6>(dev, iRow, kk); break;

          case 5: strip_mine<5>(dev, iRow, kk); break;

          case 4: strip_mine<4>(dev, iRow, kk); break;

          case 3: strip_mine<3>(dev, iRow, kk); break;

          case 2: strip_mine<2>(dev, iRow, kk); break;

          case 1: strip_mine<1>(dev, iRow, kk); break; // was strip_mine_1
        }
      }
    }
  }


  KOKKOS_INLINE_FUNCTION void operator()(SimpleTag, const ordinal_type i) const {

    if (i >= A_.numRows()) {
      return;
    }

    const KokkosSparse::SparseRowViewConst<AMatrix> row = 
        RowViewer::view(A_, offRankOffsets_, i);

    // this implementation may be best for a single vector (not multivector)
    for (ordinal_type k = 0; k < Y_.extent(1); ++k) {
      y_value_type sum = 0;

      for (ordinal_type ri = 0; ri < row.length; ++ri) {
        value_type A_ij = row.value(ri);
        ordinal_type j = row.colidx(ri);
        value_type X_jk = X_(j, k);
        // if (j >= X_.extent(0)) {
        //   printf("X row violation %d >= %d\n", int(j), int(X_.extent(0)));
        // }
        // if (X_jk > 30 || X_jk != X_jk) {
        //   printf("XXXXXXX %d %d %d %e\n", int(i), int(j), int(k), double(X_jk));
        // }
        // if (A_ij > 30 || A_ij != A_ij) {
        //   printf("AAAAAAA %d %d %d %e\n", int(i), int(j), int(k), double(A_ij));
        // }
        sum += A_ij * X_jk;
      }
      sum *= alpha_;

      if (0 == beta_) {
        Y_(i,k) = sum;
      } else {
        Y_(i,k) = beta_ * Y_(i,k) + sum;
      } 
    }
  }

  static void launch(const Alpha &alpha, 
    const AMatrix &A, 
    XVector &X, 
    const Beta &beta, 
    YVector &Y,
    const OffsetDeviceViewType &offRankOffsets,
    const execution_space &space) {

#if 0
    std::cerr << __FILE__<<":"<<__LINE__<<": " << Tpetra::getDefaultComm()->getRank() << ","
              << " y[" << Y.extent(0) << "] = A[" << A.numRows() <<"," <<A.numCols() << "] * x[" << X.extent(0) << "]\n";
#endif
#if 0
    int xNans;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<execution_space>(space, 0, X.extent(0)),
        [&](size_t i, int &lnans) {
          for (size_t k = 0; k < X.extent(1); ++k) {
            lnans += std::isnan(X(i,k)) || std::isinf(X(i,k));
          }
        }
    ,xNans);

  int yNans;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<execution_space>(space, 0, Y.extent(0)),
        [&](size_t i, int &lnans) {
          for (size_t k = 0; k < Y.extent(1); ++k) {
            lnans += std::isnan(Y(i,k)) || std::isinf(Y(i,k));
          }
        }
    ,yNans);

    if (std::is_same<RowViewer, OnRankRowViewer<OffsetDeviceViewType>>::value) {
      std::cerr << __FILE__<<":"<<__LINE__<<": " << Tpetra::getDefaultComm()->getRank() << ","
                << " SpmvMvFunctor::on-rank  alpha=" << alpha << " beta=" << beta << " xNans=" << xNans << " yNans=" << yNans << "\n";
    } else if (std::is_same<RowViewer, OffRankRowViewer<OffsetDeviceViewType>>::value) {
      std::cerr << __FILE__<<":"<<__LINE__<<": " << Tpetra::getDefaultComm()->getRank() << ","
                << " SpmvMvFunctor::off-rank alpha=" << alpha << " beta=" << beta << " xNans=" << xNans << " yNans=" << yNans << "\n";
    }
#endif


    // still need to scal A.numRows() == 0
    if (0 == A.numRows() || 0 == alpha) {
      // TODO: define a common scal operator
      SpmvMvFunctor op(alpha, A, X, beta, Y, offRankOffsets, 0/*unused*/, 0/*unused*/);
      Kokkos::parallel_for("SpmvMvFunctor scal",
                          Kokkos::RangePolicy<ScalTag, execution_space>(
                              space, 0, Y.extent(0)),
                          op);
      return;
    }


#if 0
#if defined(KOKKOS_ENABLE_SERIAL)
    if (std::is_same<execution_space, Kokkos::Serial>::value) {
      std::cerr << __FILE__<<":"<<__LINE__<<": SpmvMvFunctor::launch simple only\n";
      SpmvMvFunctor op(alpha, A, X, beta, Y, offRankOffsets, 0/*unused*/, 0/*unused*/);
      Kokkos::parallel_for("SpmvMvFunctor simple",
        Kokkos::RangePolicy<SimpleTag, execution_space>(
            space, 0, A.numRows()),
        op);
      return;
    }
#endif
#endif

    if (A.numRows() < 1) {
      throw std::runtime_error("A.numRows() < 1");
    }
    const ordinal_type NNZPerRow = A.nnz() / A.numRows();

    ordinal_type vectorLength = 1;
    while ((vectorLength * 2 * 3 <= NNZPerRow) && (vectorLength < 8)) {
      vectorLength *= 2;
    }
    

    const ordinal_type rowsPerThread =
        KokkosSparse::RowsPerThread<execution_space>(NNZPerRow);

    SpmvMvFunctor op(alpha, A, X, beta, Y, offRankOffsets,
              rowsPerThread, vectorLength);

    const ordinal_type teamSize =
        Kokkos::TeamPolicy<execution_space>(space, 
            rowsPerThread, Kokkos::AUTO, vectorLength)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    CWP_CERR("teamSize=" << teamSize << " vectorLength=" << vectorLength << "\n");
    const ordinal_type rowsPerTeam = rowsPerThread * teamSize;
    if (rowsPerTeam < 1) {
      throw std::runtime_error("rowsPerTeam < 1");
    }
    const size_type nteams = (A.numRows() + rowsPerTeam - 1) / rowsPerTeam;
    Kokkos::parallel_for("SpmvMvFunctor",
                        Kokkos::TeamPolicy<execution_space>(
                            space, nteams, teamSize, vectorLength),
                        op);

#if 0
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<execution_space>(space, 0, X.extent(0)),
        [&](size_t i, int &lnans) {
          for (size_t k = 0; k < X.extent(1); ++k) {
            lnans += std::isnan(X(i,k)) || std::isinf(X(i,k));
          }
        }
    ,xNans);

    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<execution_space>(space, 0, Y.extent(0)),
        [&](size_t i, int &lnans) {
          for (size_t k = 0; k < Y.extent(1); ++k) {
            lnans += std::isnan(Y(i,k)) || std::isinf(Y(i,k));
          }
        }
    ,yNans);

    if (std::is_same<RowViewer, OnRankRowViewer<OffsetDeviceViewType>>::value) {
      std::cerr << __FILE__<<":"<<__LINE__<<": " << Tpetra::getDefaultComm()->getRank() << ",";
      std::cerr << " SpmvMvFunctor::on-rank  xNans=" << xNans << " yNans=" << yNans << "\n";
    } else if (std::is_same<RowViewer, OffRankRowViewer<OffsetDeviceViewType>>::value) {
      std::cerr << __FILE__<<":"<<__LINE__<<": " << Tpetra::getDefaultComm()->getRank() << ",";
      std::cerr << " SpmvMvFunctor::off-rank xNans=" << xNans << " yNans=" << yNans << "\n";
    }
#endif
  }
}; // SpmvMvFunctor

/* Simplified SpmvMvTransFunctor from Kokkos Kernels
*/
template<
typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector,
typename OffsetDeviceViewType, 
typename RowViewer, // OnRankRowView or OffRankRowView
bool CONJ
> struct SpmvMvTransFunctor {
public:
  using y_value_type = typename YVector::non_const_value_type;
  using x_value_type = typename XVector::non_const_value_type;
  using value_type = typename AMatrix::non_const_value_type;
  using ordinal_type = typename AMatrix::non_const_ordinal_type; 
  using size_type = typename AMatrix::non_const_size_type; 

  using execution_space = typename AMatrix::execution_space;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::Details::ArithTraits<value_type>;

private:

  Alpha alpha_;
  AMatrix A_;
  XVector X_;
  Beta beta_;
  YVector Y_;
  OffsetDeviceViewType offRankOffsets_;
public:
  ordinal_type rowsPerTeam_ = 0;
private:

  SpmvMvTransFunctor(const Alpha& alpha, const AMatrix& A,
                     const XVector& X, const Beta& beta,
                     const YVector& Y, const OffsetDeviceViewType &offRankOffsets)
      : alpha_(alpha),
        A_(A),
        X_(X),
        beta_(beta),
        Y_(Y),
        offRankOffsets_(offRankOffsets) {}

public:
  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {

    const ordinal_type n = X_.extent(1);

    const ordinal_type teamWork = dev.league_rank() * rowsPerTeam_;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, rowsPerTeam_), [&](ordinal_type loop) {
          // iRow represents a row of the matrix, so its correct type is
          // ordinal_type.
          const ordinal_type iRow = teamWork + loop;
          if (iRow >= A_.numRows()) {
            return;
          }

          const auto row = RowViewer::view(A_, offRankOffsets_, iRow);

          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(dev, row.length),
              [&](ordinal_type iEntry) {
                const value_type val =
                    CONJ
                        ? Kokkos::Details::ArithTraits<value_type>::conj(
                              row.value(iEntry))
                        : row.value(iEntry);
                const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
                for (ordinal_type k = 0; k < n; ++k) {
                  Kokkos::atomic_add(
                      &Y_(ind, k),
                      static_cast<y_value_type>(alpha_ * val * X_(iRow, k)));
                }
              });
        });
  }

  // a kind of crummy scal implementation that we can put in an execution space instance
  struct ScalTag{};
  KOKKOS_INLINE_FUNCTION void operator()(const ScalTag&, const ordinal_type i) const {
    if (i >= ordinal_type(Y_.extent(0))) {return;}
    for (ordinal_type k = 0; k < ordinal_type(Y_.extent(1)); ++k) {
      if (0 == beta_) {
        Y_(i, k) = 0;
      } else {
        Y_(i, k) *= beta_;
      }
    }
  }

  static void launch(const Alpha &alpha, 
    const AMatrix &A, 
    XVector &X, 
    const Beta &beta, 
    YVector &Y,
    const OffsetDeviceViewType &offRankOffsets,
    const execution_space &space) {

    SpmvMvTransFunctor op(alpha, A, X, beta, Y, offRankOffsets);
    Kokkos::parallel_for("SpmvMvTransFunctor scal",
                        Kokkos::RangePolicy<ScalTag, execution_space>(
                            space, 0, Y.extent(0)),
                        op);

    // nothing to do
    if (0 == A.numRows()) {
      return;
    }

    const ordinal_type NNZPerRow = A.nnz() / A.numRows();

    ordinal_type vectorLength = 1;
    while ((vectorLength * 2 * 3 <= NNZPerRow) && (vectorLength < 8)) {
      vectorLength *= 2;
    }

    ordinal_type nrow = A.numRows();

    const ordinal_type rowsPerThread =
        KokkosSparse::RowsPerThread<execution_space>(NNZPerRow);
    const ordinal_type teamSize =
        Kokkos::TeamPolicy<execution_space>(space, 
            rowsPerThread, Kokkos::AUTO, vectorLength)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    const ordinal_type rowsPerTeam = rowsPerThread * teamSize;
    if (rowsPerTeam < 1) {
      throw std::runtime_error("rowsPerTeam < 1");
    }
    op.rowsPerTeam_ = rowsPerTeam;
    const size_type nteams = (nrow + rowsPerTeam - 1) / rowsPerTeam;
    Kokkos::parallel_for("SpmvMvTransFunctor",
                        Kokkos::TeamPolicy<execution_space>(
                            space, nteams, teamSize, vectorLength),
                        op);
  }

}; // SpmvMvTransFunctor

// SPMV_Transpose_Functor
// This TransposeFunctor is functional, but not necessarily performant.
template <typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector,
typename OffsetDeviceViewType, 
typename RowViewer, // OnRankRowView or OffRankRowView
bool CONJ>
struct SpmvTransFunctor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::Details::ArithTraits<value_type> ATV;
  typedef typename YVector::non_const_value_type y_value_type;

  const Alpha alpha;
  AMatrix A_;
  XVector x_;
  YVector y_;
  OffsetDeviceViewType offRankOffsets_;
  ordinal_type rows_per_team;

  SpmvTransFunctor(const Alpha& alpha_, const AMatrix&A,
                         const XVector& x, const YVector& y, const OffsetDeviceViewType &offRankOffsets)
      : alpha(alpha_), A_(A), x_(x), y_(y), offRankOffsets_(offRankOffsets) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type iRow) const {
    const auto row = RowViewer::view(A_, offRankOffsets_, iRow);
    const ordinal_type row_length = row.length;
    for (ordinal_type iEntry = 0; iEntry < row_length; iEntry++) {
      const value_type val =
          CONJ ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);
      Kokkos::atomic_add(&y_(ind),
                         static_cast<y_value_type>(alpha * val * x_(iRow)));
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    const ordinal_type teamWork = dev.league_rank() * rows_per_team;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, rows_per_team), [&](ordinal_type loop) {
          // iRow represents a row of the matrix, so its correct type is
          // ordinal_type.
          const ordinal_type iRow = teamWork + loop;
          if (iRow >= A_.numRows()) {
            return;
          }

          const auto row = RowViewer::view(A_, offRankOffsets_, iRow);
          const ordinal_type row_length = row.length;
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(dev, row_length),
              [&](ordinal_type iEntry) {
                const value_type val = CONJ ? ATV::conj(row.value(iEntry))
                                                 : row.value(iEntry);
                const ordinal_type ind = row.colidx(iEntry);
                Kokkos::atomic_add(&y_(ind), static_cast<y_value_type>(
                                                  alpha * val * x_(iRow)));
              });
        });
  }

static void launch(const Alpha &alpha, 
    const AMatrix &A, 
    XVector &X, 
    const Beta &beta, 
    YVector &Y,
    const OffsetDeviceViewType &offRankOffsets,
    const execution_space &space) {

    if (1 != beta) {
      KokkosBlas::scal(Y, beta, Y);
    }

    if (A.numRows() <= static_cast<ordinal_type>(0)) {
      return;
    }

    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {

      // Assuming that no row contains duplicate entries, NNZPerRow
      // cannot be more than the number of columns of the matrix.  Thus,
      // the appropriate type is ordinal_type.
      const ordinal_type NNZPerRow = A.nnz() / A.numRows();

      int vector_length     = 1;
      int max_vector_length = 1;
  #ifdef KOKKOS_ENABLE_CUDA
      if (std::is_same<execution_space, Kokkos::Cuda>::value)
        max_vector_length = 32;
  #endif
  #ifdef KOKKOS_ENABLE_HIP
      if (std::is_same<execution_space, Kokkos::Experimental::HIP>::value)
        max_vector_length = 64;
  #endif
      while ((vector_length * 2 * 3 <= NNZPerRow) &&
            (vector_length < max_vector_length))
        vector_length *= 2;

      typename AMatrix::const_ordinal_type nrow = A.numRows();

      SpmvTransFunctor op(alpha, A, X, Y, offRankOffsets);

      const ordinal_type rows_per_thread =
          KokkosSparse::RowsPerThread<execution_space>(NNZPerRow);
      const ordinal_type team_size =
          Kokkos::TeamPolicy<execution_space>(rows_per_thread, Kokkos::AUTO,
                                              vector_length)
              .team_size_recommended(op, Kokkos::ParallelForTag());
      const ordinal_type rows_per_team = rows_per_thread * team_size;
      if (rows_per_team <= 0) {
        throw std::runtime_error("rows_per_team <= 0");
      }
      op.rows_per_team                 = rows_per_team;
      const ordinal_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
      Kokkos::parallel_for(
          "KokkosSparse::spmv<Transpose>",
          Kokkos::TeamPolicy<execution_space>(space, nteams, team_size, vector_length),
          op);
    } else {
        Kokkos::parallel_for("KokkosSparse::spmv<Transpose>",
                            Kokkos::RangePolicy<execution_space>(space, 0, A.numRows()),
                            SpmvTransFunctor(alpha, A, X, Y, offRankOffsets));
    }
  }

};


/*! \brief
*/
template <typename RowViewer, typename Alpha, typename AMatrix, typename XVector, 
typename Beta, typename YVector, typename RowOffsetView>
void spmv(const typename AMatrix::execution_space &execSpace, const Alpha &alpha, const AMatrix &A, const XVector &X, const Beta &beta, const YVector &Y,
const Teuchos::ETransp &mode, const RowOffsetView &offRankOffsets) {

  using execution_space = typename AMatrix::execution_space;

  static_assert(Kokkos::SpaceAccessibility<execution_space, typename AMatrix::memory_space>::accessible,
  "execution space must be able to access A");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
  "execution space must be able to access X");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YVector::memory_space>::accessible,
  "execution space must be able to access Y");

  // unmanaged versions
  using UX = Kokkos::View<typename XVector::data_type, typename XVector::array_layout,
        typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UY = Kokkos::View<typename YVector::data_type, typename YVector::array_layout,
        typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UA = KokkosSparse::CrsMatrix<typename AMatrix::value_type,
    typename AMatrix::ordinal_type,
    typename AMatrix::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename AMatrix::size_type>;
  using UR = Kokkos::View<typename RowOffsetView::data_type, typename RowOffsetView::array_layout,
        typename RowOffsetView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  UA uA(A);
  UX uX(X);
  UY uY(Y);
  UR uOffRankOffsets(offRankOffsets);

  CWP_CERR(__FILE__ << ":" << __LINE__
           << " "   << Y.extent(0) << "," << Y.extent(1)
           << " = " << A.numRows() << "," << A.numCols()
           << " x " << X.extent(0) << "," << X.extent(1)
           << "\n");

  // TEUCHOS_TEST_FOR_EXCEPTION(Y.extent(1) != X.extent(1), std::logic_error, 
  // "Y cols " << Y.extent(1) << " != " << " X cols " << X.extent(1));

  // launch conjugate / transpose / on-rank / off-rank / vector / multivector SpMV
  switch(mode) {
    case Teuchos::ETransp::TRANS: {
      if (1 == X.extent(1) && 1 == Y.extent(1)) {
        auto X0 = Kokkos::subview(uX, Kokkos::ALL, 0);
        auto Y0 = Kokkos::subview(uY, Kokkos::ALL, 0);
        using Op = SpmvTransFunctor<Alpha, AMatrix, decltype(X0), Beta, decltype(Y0), RowOffsetView, RowViewer, false>;
        Op::launch(alpha, uA, X0, beta, Y0, uOffRankOffsets, execSpace);
      } else {
        using Op = SpmvMvTransFunctor<Alpha, AMatrix, XVector, Beta, YVector, RowOffsetView, RowViewer, false>;
        Op::launch(alpha, uA, uX, beta, uY, uOffRankOffsets, execSpace);
        return;
      }
    }
    case Teuchos::ETransp::NO_TRANS: {

      // TEUCHOS_TEST_FOR_EXCEPTION(int64_t(Y.extent(0)) != int64_t(A.numRows()), std::logic_error, 
      // "y rows " << Y.extent(0) << " != " << " A rows " << A.numRows());

      // TEUCHOS_TEST_FOR_EXCEPTION(int64_t(A.numCols()) != int64_t(X.extent(0)), std::logic_error, 
      // "A cols " << A.numCols() << " != " << " X rows " << X.extent(0));

       if(int64_t(Y.extent(0)) != int64_t(A.numRows())) {
        CWP_CERR(__FILE__ <<":" << __LINE__ << " y rows " << Y.extent(0) << " != " << " A rows " << A.numRows() << "\n");
       }

       if(int64_t(A.numCols()) != int64_t(X.extent(0))) {
        CWP_CERR(__FILE__ <<":" << __LINE__ << " A cols " << A.numCols() << " != " << " X rows " << X.extent(0) << "\n");
       }

      if (1 == X.extent(1) && 1 == Y.extent(1)) {
        CWP_CERR(__FILE__ << ":" << __LINE__ << "\n");
        auto X0 = Kokkos::subview(uX, Kokkos::ALL, 0);
        auto Y0 = Kokkos::subview(uY, Kokkos::ALL, 0);
        using Op = SpmvFunctor<Alpha, AMatrix, decltype(X0), Beta, decltype(Y0), RowOffsetView, RowViewer, false>;
        CWP_CERR(__FILE__ << ":" << __LINE__ << "\n");
        Op::launch(alpha, uA, X0, beta, Y0, uOffRankOffsets, execSpace);
        CWP_CERR(__FILE__ << ":" << __LINE__ << "\n");
      } else {
        using Op = SpmvMvFunctor<Alpha, AMatrix, XVector, Beta, YVector, RowOffsetView, RowViewer, false>;
        Op::launch(alpha, uA, uX, beta, uY, uOffRankOffsets, execSpace);
      }
      return;
    }
    case Teuchos::ETransp::CONJ_TRANS: {
      if (1 == X.extent(1) && 1 == Y.extent(1)) {
        auto X0 = Kokkos::subview(uX, Kokkos::ALL, 0);
        auto Y0 = Kokkos::subview(uY, Kokkos::ALL, 0);
        using Op = SpmvTransFunctor<Alpha, AMatrix, decltype(X0), Beta, decltype(Y0), RowOffsetView, RowViewer, true>;
        Op::launch(alpha, uA, X0, beta, Y0, uOffRankOffsets, execSpace);
      } else {
        using Op = SpmvMvTransFunctor<Alpha, AMatrix, XVector, Beta, YVector, RowOffsetView, RowViewer, true>;
        Op::launch(alpha, uA, uX, beta, uY, uOffRankOffsets, execSpace);
        return;
      }
    }
    default:
      throw std::runtime_error("unexpected Teuchos::ETransp mode in on-rank SpMV");
  }

} // spmv



} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_SPMV_HPP
