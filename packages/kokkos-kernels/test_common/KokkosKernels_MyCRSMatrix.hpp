//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include "KokkosKernels_Utils.hpp"
namespace MyKokkosSparse {

template <typename OrdinalType, typename Device, typename MemoryTraits, typename SizeType>
class StaticCrsGraph {
 public:
  typedef OrdinalType data_type;
  typedef typename Device::execution_space execution_space;
  typedef Device device_type;
  typedef SizeType size_type;

  typedef Kokkos::View<const size_type*, device_type> row_map_type;
  typedef Kokkos::View<data_type*, device_type> entries_type;

  entries_type entries;
  row_map_type row_map;
  OrdinalType num_cols;

  //! Construct an empty view.
  StaticCrsGraph() : entries(), row_map(), num_cols() {}

  //! Copy constructor (shallow copy).
  StaticCrsGraph(const StaticCrsGraph& rhs) : entries(rhs.entries), row_map(rhs.row_map), num_cols(rhs.num_cols) {}

  template <class EntriesType, class RowMapType>
  StaticCrsGraph(const EntriesType& entries_, const RowMapType& row_map_) : entries(entries_), row_map(row_map_) {}
  template <class EntriesType, class RowMapType>
  StaticCrsGraph(const EntriesType& entries_, const RowMapType& row_map_, OrdinalType numCols_)
      : entries(entries_), row_map(row_map_), num_cols(numCols_) {}
  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  StaticCrsGraph& operator=(const StaticCrsGraph& rhs) {
    entries = rhs.entries;
    row_map = rhs.row_map;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  data_type numCols() const { return num_cols; }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~StaticCrsGraph() {}
  KOKKOS_INLINE_FUNCTION
  data_type numRows() const {
    return (row_map.extent(0) != 0) ? row_map.extent(0) - static_cast<size_type>(1) : static_cast<size_type>(0);
  }
};

template <typename ScalarType, typename OrdinalType, typename Device, typename MemoryTraits, typename SizeType>
class CrsMatrix {
 public:
  typedef typename Kokkos::ViewTraits<ScalarType*, Device, void, void>::host_mirror_space host_mirror_space;

  typedef typename Device::execution_space execution_space;
  typedef typename Device::memory_space memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;
  typedef ScalarType value_type;
  typedef OrdinalType ordinal_type;
  typedef MemoryTraits memory_traits;
  typedef SizeType size_type;

  typedef StaticCrsGraph<OrdinalType, Device, MemoryTraits, SizeType> StaticCrsGraphType;
  typedef typename StaticCrsGraphType::entries_type index_type;
  typedef typename index_type::non_const_value_type const_ordinal_type;
  typedef typename index_type::non_const_value_type non_const_ordinal_type;
  typedef typename StaticCrsGraphType::row_map_type row_map_type;
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits> values_type;
  typedef CrsMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits, SizeType> HostMirror;

  StaticCrsGraphType graph;
  values_type values;
  CrsMatrix() : numCols_(0) {}
  CrsMatrix(const std::string& label, const OrdinalType& ncols, const values_type& vals,
            const StaticCrsGraphType& graph_)
      : graph(graph_), values(vals), numCols_(ncols) {}

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const { return graph.numRows(); }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const { return numCols_; }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz() const { return graph.entries.extent(0); }
  ordinal_type numCols_;
};

template <typename myExecSpace, typename crsMat_t>
crsMat_t get_crsmat(typename crsMat_t::row_map_type::non_const_type::value_type* xadj,
                    typename crsMat_t::index_type::non_const_type::value_type* adj,
                    typename crsMat_t::values_type::non_const_type::value_type* ew,
                    typename crsMat_t::row_map_type::non_const_type::value_type ne,
                    typename crsMat_t::index_type::non_const_type::value_type nv, int is_one_based) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
  typedef typename crsMat_t::index_type::non_const_type cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;

  typedef typename row_map_view_t::value_type size_type;
  typedef typename cols_view_t::value_type lno_t;
  typedef typename values_view_t::value_type scalar_t;

  row_map_view_t rowmap_view("rowmap_view", nv + 1);
  cols_view_t columns_view("colsmap_view", ne);
  values_view_t values_view("values_view", ne);

  KokkosKernels::Impl::copy_vector<scalar_t*, values_view_t, myExecSpace>(ne, ew, values_view);
  KokkosKernels::Impl::copy_vector<lno_t*, cols_view_t, myExecSpace>(ne, adj, columns_view);
  KokkosKernels::Impl::copy_vector<size_type*, row_map_view_t, myExecSpace>(nv + 1, xadj, rowmap_view);

  size_type ncols = 0;
  KokkosKernels::Impl::view_reduce_max<cols_view_t, myExecSpace>(ne, columns_view, ncols);
  ncols += 1;

  if (is_one_based) {
    // if algorithm is mkl_csrmultcsr convert to 1 base so that we dont
    // dublicate the memory at the experiments/
    KokkosKernels::Impl::kk_a_times_x_plus_b<row_map_view_t, row_map_view_t, int, int, myExecSpace>(nv + 1, rowmap_view,
                                                                                                    rowmap_view, 1, 1);
    KokkosKernels::Impl::kk_a_times_x_plus_b<cols_view_t, cols_view_t, int, int, myExecSpace>(ne, columns_view,
                                                                                              columns_view, 1, 1);
  }

  graph_t static_graph(columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
  return crsmat;
}

template <typename myExecSpace, typename in_crsMat_t, typename out_crsMat_t>
out_crsMat_t copy_crsmat(in_crsMat_t inputMat) {
  /*
      typedef typename out_crsMat_t::StaticCrsGraphType graph_t;
      typedef typename out_crsMat_t::row_map_type::non_const_type
     row_map_view_t; typedef typename out_crsMat_t::index_type::non_const_type
     cols_view_t; typedef typename out_crsMat_t::values_type::non_const_type
     values_view_t;


      typedef typename in_crsMat_t::StaticCrsGraphType in_graph_t;
      typedef typename in_crsMat_t::row_map_type::const_type in_row_map_view_t;
      typedef typename in_crsMat_t::index_type::const_type   in_cols_view_t;
      typedef typename in_crsMat_t::values_type::const_type in_values_view_t;

      typedef typename row_map_view_t::value_type size_type;
      typedef typename cols_view_t::value_type   lno_t;
      typedef typename values_view_t::value_type scalar_t;



      const size_type nv = inputMat.numRows();
      const size_type ne = inputMat.graph.entries.extent(0);

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);

      KokkosKernels::Impl::copy_vector<in_values_view_t , values_view_t,
     myExecSpace>(ne, inputMat.values, values_view);
      KokkosKernels::Impl::copy_vector<in_cols_view_t , cols_view_t,
     myExecSpace>(ne, inputMat.graph.entries, columns_view);
      KokkosKernels::Impl::copy_vector<in_row_map_view_t , row_map_view_t,
     myExecSpace>(nv+1, inputMat.graph.row_map, rowmap_view);

      size_type ncols = 0;
      KokkosKernels::Impl::view_reduce_max<cols_view_t, myExecSpace>(ne,
     columns_view, ncols); ncols += 1;

      graph_t static_graph (columns_view, rowmap_view);
      out_crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
      return crsmat;
      */
}
}  // namespace MyKokkosSparse
