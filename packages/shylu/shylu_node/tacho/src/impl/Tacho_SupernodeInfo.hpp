// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_SUPERNODE_INFO_HPP__
#define __TACHO_SUPERNODE_INFO_HPP__

#include "Tacho_Util.hpp"
#if defined(KOKKOS_ENABLE_CUDA)
 #include <cusparse_v2.h>
#elif defined(KOKKOS_ENABLE_HIP)
 #if __has_include(<rocm-core/rocm_version.h>)
  #include <rocm-core/rocm_version.h>
 #else
  #include <rocm_version.h>
 #endif
 #include <rocsparse/rocsparse.h>
 #define ROCM_VERSION ROCM_VERSION_MAJOR * 10000 + ROCM_VERSION_MINOR * 100 + ROCM_VERSION_PATCH
#endif

/// \file Tacho_SupernodeInfo.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

struct SuperNodeInfoInitReducer {
  using reducer = SuperNodeInfoInitReducer;
  struct ValueType {
    ordinal_type max_nchildren;
    ordinal_type max_supernode_size;
    ordinal_type max_num_cols;
    ordinal_type max_schur_size;
    size_type nnz;
  };
  using value_type = struct ValueType;

  using result_view_type = Kokkos::View<value_type, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  value_type *value;

  KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer() : value() {}                                         //= default;
  KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer(const SuperNodeInfoInitReducer &b) : value(b.value) {} //= default;
  KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer(value_type &val) : value(&val) {}

  KOKKOS_INLINE_FUNCTION void join(value_type &dst, value_type &src) const {
    dst.max_nchildren = (src.max_nchildren > dst.max_nchildren ? src.max_nchildren : dst.max_nchildren);
    dst.max_supernode_size =
        (src.max_supernode_size > dst.max_supernode_size ? src.max_supernode_size : dst.max_supernode_size);
    dst.max_num_cols =
        (src.max_num_cols > dst.max_num_cols ? src.max_num_cols : dst.max_num_cols);
    dst.max_schur_size = (src.max_schur_size > dst.max_schur_size ? src.max_schur_size : dst.max_schur_size);
    dst.nnz += src.nnz;
  }

  KOKKOS_INLINE_FUNCTION void join(value_type &dst, const value_type &src) const {
    dst.max_nchildren = (src.max_nchildren > dst.max_nchildren ? src.max_nchildren : dst.max_nchildren);
    dst.max_supernode_size =
        (src.max_supernode_size > dst.max_supernode_size ? src.max_supernode_size : dst.max_supernode_size);
    dst.max_num_cols =
        (src.max_num_cols > dst.max_num_cols ? src.max_num_cols : dst.max_num_cols);
    dst.max_schur_size = (src.max_schur_size > dst.max_schur_size ? src.max_schur_size : dst.max_schur_size);
    dst.nnz += src.nnz;
  }

  KOKKOS_INLINE_FUNCTION void init(value_type &val) const {
    val.max_nchildren = Kokkos::reduction_identity<ordinal_type>::max();
    val.max_supernode_size = Kokkos::reduction_identity<ordinal_type>::max();
    val.max_num_cols = Kokkos::reduction_identity<ordinal_type>::max();
    val.max_schur_size = Kokkos::reduction_identity<ordinal_type>::max();
    val.nnz = Kokkos::reduction_identity<size_type>::sum();
  }

  KOKKOS_INLINE_FUNCTION
  value_type &reference() { return *value; }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return result_view_type(value); }
};

template <typename ValueType, typename DeviceType> struct SupernodeInfo {
  using value_type = ValueType;

  using device_type = DeviceType;
  using exec_space = typename device_type::execution_space;

  using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
  using host_memory_space = typename host_device_type::memory_space;

  using crs_matrix_type = CrsMatrixBase<value_type, device_type>;

  using ordinal_type_array = Kokkos::View<ordinal_type *, device_type>;
  using size_type_array = Kokkos::View<size_type *, device_type>;
  using value_type_array = Kokkos::View<value_type *, device_type>;
  using int_type_array = Kokkos::View<int*, Kokkos::LayoutLeft, device_type>;

  using ordinal_pair_type = Kokkos::pair<ordinal_type, ordinal_type>;
  using ordinal_pair_type_array = Kokkos::View<ordinal_pair_type *, device_type>;
  using value_type_matrix = Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type>;
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  struct Supernode {
    mutable int32_t lock;

    ordinal_type row_begin; // beginning row
    ordinal_type m, n;      // panel dimension

    // column connectivity (gid - dof, sid - supernode)
    ordinal_type gid_col_begin, gid_col_end, sid_col_begin, sid_col_end;
    ordinal_type nchildren, *children; // children[MaxDependenceSize]; // hierarchy

    ordinal_type max_decendant_schur_size;     // workspace
    ordinal_type max_decendant_supernode_size; // workspace

    value_type *l_buf, *u_buf;

    bool do_not_apply_pivots;

    // for using SpMV
    size_t nnzU;
    int* rowptrU;
    int* colindU;
    value_type* nzvalsU;

    size_t nnzL;
    int* rowptrL;
    int* colindL;
    value_type* nzvalsL;

    bool spmv_explicit_transpose;
#if defined(KOKKOS_ENABLE_CUDA)
    cusparseSpMatDescr_t U_cusparse;
    cusparseSpMatDescr_t L_cusparse;
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_spmat_descr descrU;
    rocsparse_spmat_descr descrL;
#endif

    KOKKOS_INLINE_FUNCTION
    Supernode()
        : lock(0), row_begin(0), m(0), n(0), gid_col_begin(0), gid_col_end(0), sid_col_begin(0), sid_col_end(0),
          nchildren(0), children(NULL), max_decendant_schur_size(0), max_decendant_supernode_size(0), l_buf(NULL),
          u_buf(NULL), do_not_apply_pivots(false) {
      // for (ordinal_type i=0;i<MaxDependenceSize;++i) children[i] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    Supernode(const Supernode &b)
        : lock(0), row_begin(b.row_begin), m(b.m), n(b.n), gid_col_begin(b.gid_col_begin), gid_col_end(b.gid_col_end),
          sid_col_begin(b.sid_col_begin), sid_col_end(b.sid_col_end), nchildren(b.nchildren), children(b.children),
          max_decendant_schur_size(b.max_decendant_schur_size),
          max_decendant_supernode_size(b.max_decendant_supernode_size), l_buf(b.l_buf), u_buf(b.u_buf),
          do_not_apply_pivots(b.do_not_apply_pivots) {
      // for (ordinal_type i=0;i<b.nchildren;++i) children[i] = b.children[i];
    }

    KOKKOS_INLINE_FUNCTION
    ~Supernode() {}
  };
  using supernode_type = struct Supernode;
  using supernode_type_array = Kokkos::View<supernode_type *, device_type>;

  ///
  /// info for symbolic
  ///
  // ConstUnmanagedViewType<supernode_type_array> supernodes;
  UnmanagedViewType<supernode_type_array> supernodes;

  /// dof mapping to sparse matrix
  UnmanagedViewType<ordinal_type_array> gid_colidx;

  /// supernode map and panel size configuration
  /// first - sid, second - blk , blk_superpanel_colidx;
  /// the last sid is dummy but last blk is ending point of the block
  UnmanagedViewType<ordinal_pair_type_array> sid_block_colidx;

  ///
  /// max parameter
  ///
  ordinal_type max_nchildren, max_supernode_size, max_num_cols, max_schur_size;

  ///
  /// frontal matrix subassembly mode and serialization parameter
  ///
  short front_update_mode, serial_thres_size; // 0 - lock, 1 - atomic

  ///
  /// info for solve (rhs multivector)
  UnmanagedViewType<value_type_matrix> x;

  KOKKOS_INLINE_FUNCTION
  SupernodeInfo()
      : supernodes(), gid_colidx(), sid_block_colidx(), max_nchildren(), max_supernode_size(), max_num_cols(), max_schur_size(),
        front_update_mode(), serial_thres_size(), x() {}
  //= default;

  KOKKOS_INLINE_FUNCTION
  SupernodeInfo(const SupernodeInfo &b)
      : supernodes(b.supernodes), gid_colidx(b.gid_colidx), sid_block_colidx(b.sid_block_colidx),
        max_nchildren(b.max_nchildren), max_supernode_size(b.max_supernode_size), max_num_cols(b.max_num_cols), max_schur_size(b.max_schur_size),
        front_update_mode(b.front_update_mode), serial_thres_size(b.serial_thres_size), x(b.x) {}
  //= default;

  static inline void initialize(/* */ SupernodeInfo &self,
                                /* */ supernode_type_array &supernodes_,
                                /* */ ordinal_pair_type_array &sid_block_colidx_,
                                /* */ value_type_array &superpanel_buf_,
                                /// constrol
                                bool allocate_l_buf_,
                                // symbolic input
                                const ordinal_type_array &snodes_, const size_type_array &gid_ptr_,
                                const ordinal_type_array &gid_colidx_, const size_type_array &sid_ptr_,
                                const ordinal_type_array &sid_colidx_, const ordinal_type_array &blk_colidx_,
                                // tree hierarchy
                                const ordinal_type_array &stree_parent_, const size_type_array &stree_ptr_,
                                const ordinal_type_array &stree_children_) {
    const ordinal_type nsupernodes = snodes_.extent(0) - 1;

    /// allocate and assign supernodes
    supernodes_ = supernode_type_array("supernodes", nsupernodes); // managed view

    sid_block_colidx_ = ordinal_pair_type_array("sid_block_colidx", sid_colidx_.span());

    // by default, update mode is atomic: 0 - mutex lock, 1 - atomic
    self.front_update_mode = 1;
    self.serial_thres_size = 0;

    /// workspace parameter initialization
    self.max_nchildren = 0;
    self.max_supernode_size = 0;
    self.max_num_cols = 0;
    self.max_schur_size = 0;

    Kokkos::RangePolicy<exec_space> supernodes_range_policy(0, nsupernodes);
    SuperNodeInfoInitReducer::value_type init_reduce_val;
    Kokkos::parallel_reduce(
        supernodes_range_policy,
        KOKKOS_LAMBDA(const ordinal_type &sid, SuperNodeInfoInitReducer::value_type &update) {
          auto &s = supernodes_(sid);

          s.row_begin = snodes_(sid);
          s.m = snodes_(sid + 1) - snodes_(sid);
          s.n = blk_colidx_(sid_ptr_(sid + 1) - 1);

          s.gid_col_begin = gid_ptr_(sid);
          s.gid_col_end = gid_ptr_(sid + 1);
          s.sid_col_begin = sid_ptr_(sid);
          s.sid_col_end = sid_ptr_(sid + 1);

          for (ordinal_type i = s.sid_col_begin; i < s.sid_col_end; ++i) {
            sid_block_colidx_(i).first = sid_colidx_(i);
            sid_block_colidx_(i).second = blk_colidx_(i);
          }

          s.nchildren = stree_ptr_(sid + 1) - stree_ptr_(sid);
          s.children = &stree_children_(stree_ptr_(sid));
          // const ordinal_type offset = stree_ptr_(sid);
          // for (ordinal_type i=0;i<s.nchildren;++i)
          //   s.children[i] = stree_children_(offset + i);

          update.max_nchildren = max(update.max_nchildren, s.nchildren);
          update.max_supernode_size = max(update.max_supernode_size, s.m);
          update.max_num_cols = max(update.max_num_cols, s.n);
          update.max_schur_size = max(update.max_schur_size, s.n - s.m);

          update.nnz += (s.m * s.n);                         /// upper
          update.nnz += (allocate_l_buf_ ? (s.m * s.n) : 0); /// lower

          s.max_decendant_supernode_size = s.m;
          s.max_decendant_schur_size = s.n - s.m;
        },
        SuperNodeInfoInitReducer(init_reduce_val));

    {
      // finding max decendant supernode information is sequential
      auto h_supernodes = Kokkos::create_mirror_view(host_memory_space(), supernodes_);
      auto h_stree_parent = Kokkos::create_mirror_view(host_memory_space(), stree_parent_);
      Kokkos::deep_copy(h_supernodes, supernodes_);
      Kokkos::deep_copy(h_stree_parent, stree_parent_);
      for (ordinal_type sid = 0; sid < nsupernodes; ++sid) {
        auto &s = h_supernodes(sid);
        const ordinal_type sidpar = h_stree_parent(sid);
        if (sidpar != -1) {
          auto &spar = h_supernodes(sidpar);
          spar.max_decendant_supernode_size = max(s.max_decendant_supernode_size, spar.max_decendant_supernode_size);
          spar.max_decendant_schur_size = max(s.max_decendant_schur_size, spar.max_decendant_schur_size);
        }
      }
      Kokkos::deep_copy(supernodes_, h_supernodes);
    }

    self.max_nchildren = init_reduce_val.max_nchildren;
    self.max_supernode_size = init_reduce_val.max_supernode_size;
    self.max_num_cols = init_reduce_val.max_num_cols;
    self.max_schur_size = init_reduce_val.max_schur_size;

    // supernodal factor array; data is held outside with a managed view
    // supernode does not include this view
    // for the case that the same data structure is reused, the buffer will
    // be zero'ed for each numeric factorization.
    superpanel_buf_ = value_type_array(do_not_initialize_tag("superpanel_buf"), init_reduce_val.nnz);
    Kokkos::parallel_scan(
        supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid, size_type &update, const bool &final) {
          auto &s = supernodes_(sid);
          const ordinal_type u_buf_size = s.m * s.n;
          const ordinal_type l_buf_size = allocate_l_buf_ ? (s.m * s.n) : 0;
          if (final) {
            s.u_buf = &superpanel_buf_(update);
            s.l_buf = (allocate_l_buf_ ? s.u_buf + u_buf_size : NULL);
          }
          update += (u_buf_size + l_buf_size);
        });

    self.supernodes = supernodes_; // unmanaged view, data is held outside
    self.gid_colidx = gid_colidx_;
    self.sid_block_colidx = sid_block_colidx_;
  }

  inline void initialize(/* */ supernode_type_array &supernodes_,
                         /* */ ordinal_pair_type_array &sid_block_colidx_,
                         /* */ value_type_array &superpanel_buf_,
                         /// control
                         const bool allocate_l_buf_,
                         // symbolic input
                         const ordinal_type_array &snodes_, const size_type_array &gid_ptr_,
                         const ordinal_type_array &gid_colidx_, const size_type_array &sid_ptr_,
                         const ordinal_type_array &sid_colidx_, const ordinal_type_array &blk_colidx_,
                         // tree hierarchy
                         const ordinal_type_array &stree_parent_, const size_type_array &stree_ptr_,
                         const ordinal_type_array &stree_children_) {
    initialize(*this,
               /// output
               supernodes_, sid_block_colidx_, superpanel_buf_,
               /// control
               allocate_l_buf_,
               /// super node input
               snodes_, gid_ptr_, gid_colidx_, sid_ptr_, sid_colidx_, blk_colidx_, stree_parent_, stree_ptr_,
               stree_children_);
  }

  static inline void copySparseToSuperpanels(SupernodeInfo &self,
                                             /// control
                                             const bool copy_to_l_buf,
                                             /// input from sparse matrix
                                             const size_type_array &ap, const ordinal_type_array &aj,
                                             const value_type_array &ax, const ordinal_type_array &perm,
                                             const ordinal_type_array &peri) {
    const ordinal_type nsupernodes = self.supernodes.extent(0);
    using policy_type = Kokkos::TeamPolicy<exec_space, Kokkos::Schedule<Kokkos::Static>>;

    value_type_array axt;
    if (copy_to_l_buf) {
      axt = value_type_array("axt", ax.extent(0));
      policy_type policy(ap.extent(0) - 1, Kokkos::AUTO());
      Kokkos::parallel_for(
          policy, KOKKOS_LAMBDA(const typename policy_type::member_type &member) {
            const ordinal_type row = member.league_rank();
            const ordinal_type kbeg = ap(row), kend = ap(row + 1);

            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, kbeg, kend), [&](const ordinal_type &k) {
              const ordinal_type i = aj(k), j = row;
              {
                const ordinal_type lbeg = ap(i), lend = ap(i + 1);
                ordinal_type *first = aj.data() + lbeg;
                ordinal_type *last  = aj.data() + lend;
                ordinal_type *loc =
                    lower_bound(first, last, j, [](ordinal_type left, ordinal_type right) { return left < right; });
                TACHO_TEST_FOR_ABORT(*loc != j, "transpose fail");
                axt(lbeg + loc - first) = ax(k);
              }
            });
          });
    }

    {
      policy_type policy(nsupernodes, Kokkos::AUTO()); // team and vector sizes are AUTO selected.
      Kokkos::parallel_for(
          policy, KOKKOS_LAMBDA(const typename policy_type::member_type &member) {
            const ordinal_type sid = member.league_rank();
            const auto s = self.supernodes(sid);
            /// copy to upper triangular
            {
              UnmanagedViewType<value_type_matrix> tgt_u(s.u_buf, s.m, s.n);
              UnmanagedViewType<value_type_matrix> tgt_lp(s.l_buf, s.n, s.m);
              const auto tgt_l = Kokkos::subview(tgt_lp, range_type(s.m, s.n), Kokkos::ALL());

              // row major access to sparse src
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, s.m), [&](const ordinal_type &i) {
                const ordinal_type ii = i + s.row_begin,                // row in U
                    row = perm(ii), kbeg = ap(row), kend = ap(row + 1); // row in A

                const ordinal_type jjbeg = (copy_to_l_buf ? s.row_begin : ii);
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, kbeg, kend),
                                     [&, jjbeg](const ordinal_type &k) { // Value capture is a workaround for cuda +
                                                                         // gcc-7.2 compiler bug w/c++14
                                       const ordinal_type jj = peri(aj(k) /* col in A */); // col in U
                                       if (jjbeg <= jj) {
                                         ordinal_type *first = self.gid_colidx.data() + s.gid_col_begin;
                                         ordinal_type *last = self.gid_colidx.data() + s.gid_col_end;
                                         ordinal_type *loc =
                                             lower_bound(first, last, jj, [](ordinal_type left, ordinal_type right) {
                                               return left < right;
                                             });
                                         TACHO_TEST_FOR_ABORT(*loc != jj, "copy is wrong");
                                         const ordinal_type j = loc - first;
                                         tgt_u(i, j) = ax(k);
                                         if (j >= s.m && copy_to_l_buf) {
                                           tgt_l(j - s.m, i) = axt(k);
                                         }
                                       }
                                     });
              });
            }
          });
    }
  }

  inline void copySparseToSuperpanels( // input from sparse matrix
      const bool copy_to_l_buf, const size_type_array &ap, const ordinal_type_array &aj, const value_type_array &ax,
      const ordinal_type_array &perm, const ordinal_type_array &peri) {
    copySparseToSuperpanels(*this, copy_to_l_buf, ap, aj, ax, perm, peri);
  }

  static inline void createCrsMatrix(SupernodeInfo &self, crs_matrix_type &A,
                                     const bool replace_value_with_one = false) {
    // // count m, n, nnz
    // const ordinal_type nsupernodes = self.supernodes.extent(0);

    // auto d_last = Kokkos::subview(self.supernodes, nsupernodes - 1);
    // auto h_last = Kokkos::create_mirror_view(host_memory_space(), d_last);
    // Kokkos::deep_copy(h_last, d_last);
    // auto last = h_last();
    // //auto &last = supernodes(nsupernodes - 1);

    // const ordinal_type mm = last.row_begin + last.m, nn = mm;

    // Kokkos::RangePolicy<exec_space> supernodes_range_policy(0,nsupernodes);

    // // parallel for/scan version
    // size_type_array ap_tmp("ap", mm+1);
    // Kokkos::parallel_for
    //   (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid) {
    //     // row major access to sparse src
    //     const auto &s = self.supernodes(sid);
    //     const ordinal_type soffset = s.row_begin;
    //     for (ordinal_type i=0;i<s.m;++i) {
    //       const ordinal_type row = i+soffset; // row in sparse matrix
    //       ap_tmp(row) = (s.n - i); // upper triangular only
    //     }
    //   });

    // size_type_array ap("ap", mm+1);
    // Kokkos::parallel_scan
    //   (Kokkos::RangePolicy<exec_space>(0,mm+1),
    //    KOKKOS_LAMBDA(const ordinal_type &i, size_type &update, const bool &final) {
    //     if (final)
    //       ap(i) = update;
    //     update += ap_tmp(i);
    //   });

    // // fill the matrix
    // auto d_nnz = Kokkos::subview(ap, mm);
    // auto h_nnz = Kokkos::create_mirror_view(host_memory_space(), d_nnz);
    // Kokkos::deep_copy(h_nnz, d_nnz);

    // const auto nnz = h_nnz();
    // ordinal_type_array aj("aj", nnz);
    // value_type_array ax("ax", nnz);

    // Kokkos::parallel_for
    //   (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid) {
    //     const auto &s = self.supernodes(sid);

    //     UnmanagedViewType<value_type_matrix> src(s.u_buf, s.m, s.n);

    //     // row major access to sparse src
    //     const ordinal_type
    //       soffset = s.row_begin,
    //       goffset = s.gid_col_begin;

    //     for (ordinal_type i=0;i<s.m;++i) {
    //       const size_type beg = ap(i+soffset);
    //       for (ordinal_type j=i,k=0;j<s.n;++j,++k) {
    //         const ordinal_type col = self.gid_colidx(j+goffset);
    //         aj(beg+k) = col;
    //         ax(beg+k) = replace_value_with_one ? 1.0 : src(i, j);
    //       }
    //     }
    //   });

    // // set triple to crs matrix
    // A.clear();
    // A.setExternalMatrix(mm, nn, nnz, ap, aj, ax);
  }

  inline void createCrsMatrix(crs_matrix_type &A, const bool replace_value_with_one = false) {
    createCrsMatrix(*this, A, replace_value_with_one);
  }
};

} // namespace Tacho

#endif
