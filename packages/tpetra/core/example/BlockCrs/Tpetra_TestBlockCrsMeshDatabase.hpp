// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TPETRA_TEST_BLOCKCRS_MESHDATABASE_HPP__
#define __TPETRA_TEST_BLOCKCRS_MESHDATABASE_HPP__

#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <type_traits>

#include <Teuchos_Comm.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockMultiVector.hpp>

namespace BlockCrsTest {
  // typedefs
  typedef double value_type;
  typedef double magnitude_type;

  typedef Tpetra::Map<> map_type;
  typedef typename map_type::local_ordinal_type local_ordinal_type;
  typedef typename map_type::global_ordinal_type global_ordinal_type;
  typedef typename map_type::node_type node_type;

  typedef Tpetra::Import<> tpetra_import_type;
  typedef Tpetra::Export<> tpetra_export_type;
  typedef Tpetra::MultiVector<value_type> tpetra_multivector_type;
  typedef Tpetra::MultiVector<magnitude_type> tpetra_multivector_magnitude_type;

  typedef Tpetra::CrsGraph<> tpetra_crs_graph_type;
  typedef Tpetra::CrsMatrix<value_type> tpetra_crs_matrix_type;
  typedef Tpetra::RowMatrix<value_type> tpetra_rowmatrix_type;
  typedef Tpetra::BlockCrsMatrix<value_type> tpetra_blockcrs_matrix_type;

  typedef typename node_type::device_type device_type;
  typedef typename node_type::execution_space exec_space;
  typedef typename node_type::memory_space mem_space;
  typedef Kokkos::DefaultHostExecutionSpace host_space;

  typedef Kokkos::pair<local_ordinal_type,local_ordinal_type> local_ordinal_range_type;
  typedef Kokkos::pair<global_ordinal_type,global_ordinal_type> global_ordinal_range_type;

  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;

  // utils
  template<typename T1, typename T2, typename CompareType>
  KOKKOS_INLINE_FUNCTION
  static T1* lower_bound(T1* first, T1* last, const T2& val,
                         CompareType compare) {
    T1 *it;
    local_ordinal_type step = 0, count = last - first;
    while (count > 0) {
      it = first; step = count/2; it += step;
      if (compare(*it,val)) {
        first = ++it;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    return first;
  }

  template<typename T1, typename T2>
  KOKKOS_FORCEINLINE_FUNCTION
  static T1* lower_bound(T1* first, T1* last, const T2& val) {
    return lower_bound(first, last, val, [](T1 left, T2 right) { return left < right; });
  }

  template<typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  static void heapify(T1 *v, T2 n, T2 i) {
    while (true) {
      T2 largest = i;
      T2 l = 2*i + 1;
      T2 r = 2*i + 2;

      if (l < n && v[l] > v[largest]) largest = l;
      if (r < n && v[r] > v[largest]) largest = r;
      if (largest == i)
        break;
      // swap
      T1 tmp = v[i]; v[i] = v[largest]; v[largest] = tmp;
      i = largest;
    }
  }

  template<typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  static void heap_sort(T1 *v, T2 n) {
    for (T2 i=n/2-1;i>=0;--i) heapify(v, n, i);
    for (T2 i=n-1;i>=0;--i) {
      T1 tmp = v[0]; v[0] = v[i]; v[i] = tmp;
      heapify(v, i, 0);
    }
  }

  // mesh database
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  static T
  get_block_crs_entry(const LO i, const LO j, const LO k,
                      const LO diff_i, const LO diff_j, const LO diff_k,
                      const LO l0, const LO l1) {
    // any reasonable deterministic number to fill diagonal blocks
    // inside of the block, put some weight on diagonals
    // make off diagonal blocks a bit smaller than the diagonals
    // add epsilon to remove clean zeros
    return ( ( (i%13)/14.0 + (j%17)/12.0 + (k%19)/10.0 ) +  // any number to fill diagonal blocks
             ( l0 == l1 ? 10.0*(l0%10)/10.0 : ((l0+l1)%10)/10.0) +
             ( diff_i*((i-j)%4)/18.0 + diff_j*((i-j)%7)/24.0 + diff_k*((i-j)%5)/16.0) +
             ( 1.0e-7 ) );
  }

  template<typename T>
  KOKKOS_INLINE_FUNCTION
  static T
  get_multi_vector_entry(const GO gid, const LO j) {
    return ((gid + j)%100)/100.0 + 1.0e-7;
  }

  struct MeshDatabase {
  public:
    // global perspective to the mesh structure (finite volume interior node only)
    struct StructuredProcGrid {
    public:
      typedef int process_rank_type;

      process_rank_type _num_procs_i;
      process_rank_type _num_procs_j;
      process_rank_type _num_procs_k;
      process_rank_type _num_procs_jk;
      process_rank_type _max_num_procs;
      process_rank_type _rank;
      process_rank_type _proc_i;
      process_rank_type _proc_j;
      process_rank_type _proc_k;

      void
      init (const process_rank_type num_procs_i,
            const process_rank_type num_procs_j,
            const process_rank_type num_procs_k)
      {
        _num_procs_i = num_procs_i;
        _num_procs_j = num_procs_j;
        _num_procs_k = num_procs_k;
        _num_procs_jk = num_procs_j*num_procs_k;

        const process_rank_type bigger_ij = num_procs_i > num_procs_j ? num_procs_i : num_procs_j;
        const process_rank_type bigger_jk = num_procs_j > num_procs_k ? num_procs_j : num_procs_k;
        _max_num_procs = bigger_ij > bigger_jk ? bigger_ij : bigger_jk;
      }

      void setRank (const process_rank_type rank) {
        _rank = rank;
        _proc_i = _rank / _num_procs_jk;
        _proc_k = _rank % _num_procs_jk;
        _proc_j = _proc_k / _num_procs_k;
        _proc_k = _proc_k % _num_procs_k;
      }

      StructuredProcGrid() = default;
      StructuredProcGrid(const StructuredProcGrid &b) = default;
      StructuredProcGrid (const process_rank_type num_procs_i,
                          const process_rank_type num_procs_j,
                          const process_rank_type num_procs_k)
      {
        init (num_procs_i,
              num_procs_j,
              num_procs_k);
      }

      StructuredProcGrid (const process_rank_type rank,
                          const process_rank_type num_procs_i,
                          const process_rank_type num_procs_j,
                          const process_rank_type num_procs_k)
      {
        init (num_procs_i,
              num_procs_j,
              num_procs_k);
        setRank (rank);
      }
    };

    struct StructuredBlock {
    public:
      GO _num_global_elements_i;
      GO _num_global_elements_j;
      GO _num_global_elements_k;
      GO _num_global_elements_jk;
      StructuredBlock() = default;
      StructuredBlock(const StructuredBlock &b) = default;
      StructuredBlock(const GO num_global_elements_i,
                      const GO num_global_elements_j,
                      const GO num_global_elements_k)
        : _num_global_elements_i (num_global_elements_i),
          _num_global_elements_j (num_global_elements_j),
          _num_global_elements_k (num_global_elements_k),
          _num_global_elements_jk (num_global_elements_j * num_global_elements_k)
      {}

      template<class IndexType>
      KOKKOS_INLINE_FUNCTION GO
      ijk_to_idx (const IndexType i, const IndexType j, const IndexType k) const {
        static_assert (std::is_integral<IndexType>::value,
                       "IndexType must be a built-in integer type.");
        return (i*_num_global_elements_j + j)*_num_global_elements_k + k;
      }

      template<class IndexType>
      KOKKOS_INLINE_FUNCTION void
      idx_to_ijk (const GO idx, IndexType& i, IndexType& j, IndexType& k) const {
        static_assert (std::is_integral<IndexType>::value,
                       "IndexType must be a built-in integer type.");
        i = idx / _num_global_elements_jk;
        k = idx % _num_global_elements_jk;
        j = k   / _num_global_elements_k;
        k = k   % _num_global_elements_k;
      }

      KOKKOS_INLINE_FUNCTION GO getNumElements() const {
        return _num_global_elements_i*_num_global_elements_jk;
      }
    };

    struct StructuredBlockPart {
      local_ordinal_range_type _range_i, _range_j, _range_k;
      StructuredBlockPart() = default;
      StructuredBlockPart(const StructuredBlockPart &b) = default;
      StructuredBlockPart(const LO range_i_beg, const LO range_i_end,
                          const LO range_j_beg, const LO range_j_end,
                          const LO range_k_beg, const LO range_k_end)
        : _range_i(range_i_beg, range_i_end),
          _range_j(range_j_beg, range_j_end),
          _range_k(range_k_beg, range_k_end) {}

      inline void getOwnedRange(local_ordinal_range_type &range_i,
                                local_ordinal_range_type &range_j,
                                local_ordinal_range_type &range_k) const {
        range_i = _range_i;
        range_j = _range_j;
        range_k = _range_k;
      }

      inline void getRemoteRange(const StructuredBlock &sb,
                                 local_ordinal_range_type &range_i,
                                 local_ordinal_range_type &range_j,
                                 local_ordinal_range_type &range_k) const {
        range_i.first  = _range_i.first  - (_range_i.first > 0);
        range_i.second = _range_i.second + (_range_i.second < sb._num_global_elements_i);

        range_j.first  = _range_j.first  - (_range_j.first > 0);
        range_j.second = _range_j.second + (_range_j.second < sb._num_global_elements_j);

        range_k.first  = _range_k.first  - (_range_k.first > 0);
        range_k.second = _range_k.second + (_range_k.second < sb._num_global_elements_k);
      }

      inline GO getNumElements() const {
        return ( (_range_i.second - _range_i.first)*
                 (_range_j.second - _range_j.first)*
                 (_range_k.second - _range_k.first) );
      }
    };

  public:
    typedef StructuredProcGrid::process_rank_type process_rank_type;

    Teuchos::RCP<const Teuchos::Comm<int> > _comm;
    StructuredBlock _sb;
    StructuredProcGrid _grid;
    StructuredBlockPart _owned;

    typedef Kokkos::View<GO*,device_type> global_ordinal_view_type;
    typedef Kokkos::View<GO*,host_space> global_ordinal_view_host_type;

    global_ordinal_view_host_type _element_gids;
    global_ordinal_view_host_type _owned_element_gids;
    global_ordinal_view_host_type _remote_element_gids;

    MeshDatabase() = default;
    MeshDatabase(const MeshDatabase &b) = default;

    MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm,
                 LO num_global_elements_i,
                 LO num_global_elements_j,
                 LO num_global_elements_k,
                 const process_rank_type num_procs_i,
                 const process_rank_type num_procs_j,
                 const process_rank_type num_procs_k)
      : _comm(comm),
        _sb(num_global_elements_i, num_global_elements_j, num_global_elements_k),
        _grid(_comm->getRank(), num_procs_i, num_procs_j, num_procs_k) {

      // uniform partitions on the structured block
      const LO iparts = _sb._num_global_elements_i / _grid._num_procs_i;
      const LO jparts = _sb._num_global_elements_j / _grid._num_procs_j;
      const LO kparts = _sb._num_global_elements_k / _grid._num_procs_k;

      const LO ibeg = _grid._proc_i * iparts;
      const LO jbeg = _grid._proc_j * jparts;
      const LO kbeg = _grid._proc_k * kparts;

      const LO itmp = ibeg + iparts;
      const LO jtmp = jbeg + jparts;
      const LO ktmp = kbeg + kparts;

      const LO iend = itmp < _sb._num_global_elements_i ? itmp : _sb._num_global_elements_i;
      const LO jend = jtmp < _sb._num_global_elements_j ? jtmp : _sb._num_global_elements_j;
      const LO kend = ktmp < _sb._num_global_elements_k ? ktmp : _sb._num_global_elements_k;

#if 0
      if (_grid._rank == 0) {
        printf(" gridsi= %d %d %d\n", _grid._num_procs_i, _grid._num_procs_j, _grid._num_procs_k);
        printf(" gridid= %d %d %d\n", _grid._proc_i, _grid._proc_j, _grid._proc_k);
        printf(" parts = %d %d %d\n", iparts, jparts, kparts);
        printf(" beg   = %d %d %d\n", ibeg, jbeg, kbeg);
        printf(" tmp   = %d %d %d\n", itmp, jtmp, ktmp);
        printf(" end   = %d %d %d\n", iend, jend, kend);
      }
#endif
      // local elements owned by this proc
      _owned = StructuredBlockPart(ibeg, iend,
                                   jbeg, jend,
                                   kbeg, kend);

      // count the number of owned elements
      const GO num_owned_elements = _owned.getNumElements();

      // remote elements ids
      local_ordinal_range_type remote_range_i, remote_range_j, remote_range_k;
      _owned.getRemoteRange(_sb, remote_range_i, remote_range_j, remote_range_k);

      // remote elements possibly exist for all six faces
      const local_ordinal_range_type face[6][3]
        = { { local_ordinal_range_type( remote_range_i.first, _owned._range_i.first), _owned._range_j, _owned._range_k },
            { local_ordinal_range_type(_owned._range_i.second, remote_range_i.second), _owned._range_j, _owned._range_k },
            { _owned._range_i, local_ordinal_range_type( remote_range_j.first, _owned._range_j.first), _owned._range_k },
            { _owned._range_i, local_ordinal_range_type(_owned._range_j.second, remote_range_j.second), _owned._range_k },
            { _owned._range_i, _owned._range_j, local_ordinal_range_type( remote_range_k.first, _owned._range_k.first) },
            { _owned._range_i, _owned._range_j, local_ordinal_range_type(_owned._range_k.second, remote_range_k.second) } };

      // count the number of remote elements
      GO num_remote_elements = 0;
      for (LO f=0;f<6;++f)
        num_remote_elements += ( (face[f][0].second - face[f][0].first) *
                                 (face[f][1].second - face[f][1].first) *
                                 (face[f][2].second - face[f][2].first) );

      const GO num_elements = num_owned_elements + num_remote_elements;
      _element_gids = global_ordinal_view_host_type("element gids", num_elements);
      _owned_element_gids  = Kokkos::subview(_element_gids, global_ordinal_range_type(0,num_owned_elements));
      _remote_element_gids = Kokkos::subview(_element_gids, global_ordinal_range_type(num_owned_elements,num_elements));

      {
        LO l = 0;
        for (LO i=_owned._range_i.first;i<_owned._range_i.second;++i)
          for (LO j=_owned._range_j.first;j<_owned._range_j.second;++j)
            for (LO k=_owned._range_k.first;k<_owned._range_k.second;++k)
              _owned_element_gids(l++) = _sb.ijk_to_idx(i,j,k);

        GO *beg = _owned_element_gids.data(), *end = beg + l;
        std::sort(beg, end);
      }
      {
        LO l = 0;
        for (LO f=0;f<6;++f)
          for (LO i=face[f][0].first;i<face[f][0].second;++i)
            for (LO j=face[f][1].first;j<face[f][1].second;++j)
              for (LO k=face[f][2].first;k<face[f][2].second;++k)
                _remote_element_gids(l++) = _sb.ijk_to_idx(i,j,k);
        GO *beg = _remote_element_gids.data(), *end = beg + l;
        std::sort(beg, end);
      }
    }

    ~MeshDatabase() = default;

    StructuredBlock getStructuredBlock() const { return _sb; }
    StructuredBlockPart getStructuredBlockPart() const { return _owned; }

    size_t getNumOwnedElements() const  { return _owned_element_gids.extent(0); }
    size_t getNumRemoteElements() const { return _remote_element_gids.extent(0); }
    size_t getNumElements() const {return _element_gids.extent(0);}

    global_ordinal_view_host_type getOwnedElementGlobalIDs() {return _owned_element_gids;}
    global_ordinal_view_host_type getRemoteElementGlobalIDs() {return _remote_element_gids;}

    global_ordinal_view_host_type getElementGlobalIDs() {return _element_gids;}

    // Debugging output
    void print(std::ostream & oss) {
      std::ostringstream ss;
      ss << "[" << _grid._rank << "::"
         << _grid._proc_i << "," << _grid._proc_j << "," << _grid._proc_k << "]";

      oss << ss.str() << " Global Elements = ["
          << _sb._num_global_elements_i << "x"
          << _sb._num_global_elements_j << "x"
          << _sb._num_global_elements_k <<"]\n";

      oss << ss.str() <<" Stop/Start Elements   = "
          << "[" << _owned._range_i.first << "," << _owned._range_i.second << ")x"
          << "[" << _owned._range_j.first << "," << _owned._range_j.second << ")x"
          << "[" << _owned._range_k.first << "," << _owned._range_k.second << ")\n";

      oss << ss.str()<<" Owned Global Elements = ";
      for(size_t i=0;i<_owned_element_gids.extent(0);++i) {
        oss << _owned_element_gids(i) << " ";
      }

      oss<<"\n"<<ss.str()<<" Remote Global Elements = ";
      for(size_t i=0;i<_remote_element_gids.extent(0);++i) {
        oss << _remote_element_gids[i] << " ";
      }

      oss << std::endl;
    }
  };

  struct LocalGraphConstruction {
  private:
    MeshDatabase::StructuredBlock _sb;
    MeshDatabase::global_ordinal_view_type _owned_gids;

    local_ordinal_range_type _remote_range_i;
    local_ordinal_range_type _remote_range_j;
    local_ordinal_range_type _remote_range_k;

    typedef typename tpetra_crs_graph_type::local_graph_device_type::row_map_type::non_const_type rowptr_view_type;
    rowptr_view_type _rowptr;

    typedef typename rowptr_view_type::non_const_value_type scan_value_type;

  public:
    LocalGraphConstruction(const MeshDatabase::StructuredBlock &sb,
                           const MeshDatabase::global_ordinal_view_type &owned_gids,
                           const local_ordinal_range_type &remote_range_i,
                           const local_ordinal_range_type &remote_range_j,
                           const local_ordinal_range_type &remote_range_k,
                           const rowptr_view_type &rowptr)
      : _sb(sb),
        _owned_gids(owned_gids),
        _remote_range_i(remote_range_i),
        _remote_range_j(remote_range_j),
        _remote_range_k(remote_range_k),
        _rowptr(rowptr) {};

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO &local_idx,
                    scan_value_type &update,
                    const bool final) const {
      LO cnt = 0;
      if (local_idx < static_cast<LO> (_owned_gids.extent (0))) {
        LO i, j, k;
        const GO global_idx = _owned_gids(local_idx);
        _sb.idx_to_ijk(global_idx, i, j, k);

        cnt = 1; // self

        cnt += ((i-1) >= _remote_range_i.first );
        cnt += ((i+1) <  _remote_range_i.second);

        cnt += ((j-1) >= _remote_range_j.first );
        cnt += ((j+1) <  _remote_range_j.second);

        cnt += ((k-1) >= _remote_range_k.first );
        cnt += ((k+1) <  _remote_range_k.second);
      }

      _rowptr(local_idx) = cnt;
      if (final)
        _rowptr(local_idx) = update;
      update += cnt;
    }

    KOKKOS_INLINE_FUNCTION
    void init(scan_value_type &update) const {
      update = 0;
    }

    KOKKOS_INLINE_FUNCTION
    void join(scan_value_type &update,
              const scan_value_type &input) const {
      update += input;
    }

    inline
    void run() {
      Kokkos::parallel_scan(_owned_gids.extent(0)+1, *this);
    }
  };

  struct LocalGraphFill {
  private:
    typedef typename tpetra_crs_graph_type::local_graph_device_type::row_map_type::non_const_type rowptr_view_type;
    typedef typename tpetra_crs_graph_type::local_graph_device_type::entries_type colidx_view_type;

    MeshDatabase::StructuredBlock _sb;
    MeshDatabase::global_ordinal_view_type _owned_gids;
    MeshDatabase::global_ordinal_view_type _remote_gids;

    rowptr_view_type _rowptr;
    colidx_view_type _colidx;

    local_ordinal_range_type _owned_range_i;
    local_ordinal_range_type _owned_range_j;
    local_ordinal_range_type _owned_range_k;
    local_ordinal_range_type _remote_range_i;
    local_ordinal_range_type _remote_range_j;
    local_ordinal_range_type _remote_range_k;

  public:
    LocalGraphFill(const MeshDatabase::StructuredBlock& sb,
                   const MeshDatabase::global_ordinal_view_type& owned_gids,
                   const MeshDatabase::global_ordinal_view_type& remote_gids,
                   const local_ordinal_range_type& owned_range_i,
                   const local_ordinal_range_type& owned_range_j,
                   const local_ordinal_range_type& owned_range_k,
                   const local_ordinal_range_type& remote_range_i,
                   const local_ordinal_range_type& remote_range_j,
                   const local_ordinal_range_type& remote_range_k,
                   const rowptr_view_type& rowptr,
                   const colidx_view_type& colidx)
      : _sb(sb),
        _owned_gids(owned_gids),
        _remote_gids(remote_gids),
        _rowptr(rowptr),
        _colidx(colidx),
        _owned_range_i(owned_range_i),
        _owned_range_j(owned_range_j),
        _owned_range_k(owned_range_k),
        _remote_range_i(remote_range_i),
        _remote_range_j(remote_range_j),
        _remote_range_k(remote_range_k)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const LO& local_idx) const {
      LO i, j, k;
      const GO global_idx = _owned_gids(local_idx);
      _sb.idx_to_ijk(global_idx, i, j, k);

      const LO lbeg = _rowptr(local_idx);
      LO lcnt = lbeg;

      // self
      _colidx(lcnt++) = local_idx;

      // owned and remote gids are separately sorted
      const auto num_owned_elements = _owned_gids.extent(0);
      const auto num_remote_elements = _remote_gids.extent(0);

      const auto owned_first= _owned_gids.data();
      const auto owned_last = owned_first + (num_owned_elements - 1);
      const auto remote_first = _remote_gids.data();
      const auto remote_last = remote_first + (num_remote_elements-1);

      // sides on i
      {
        const auto i_minus_one = i-1;
        if (i_minus_one >= _owned_range_i.first) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i_minus_one,j,k)) -
                              owned_first );
        } else if (i_minus_one >= _remote_range_i.first) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i_minus_one,j,k)) -
                              remote_first +
                              num_owned_elements );
        }
        const auto i_plus_one = i+1;
        if (i_plus_one < _owned_range_i.second) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i_plus_one,j,k)) -
                              owned_first );
        } else if (i_plus_one < _remote_range_i.second) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i_plus_one,j,k)) -
                              remote_first +
                              num_owned_elements );
        }
      }

      // sides on j
      {
        const auto j_minus_one = j-1;
        if (j_minus_one >= _owned_range_j.first) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i,j_minus_one,k)) -
                              owned_first );
        } else if (j_minus_one >= _remote_range_j.first) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i,j_minus_one,k)) -
                              remote_first +
                              num_owned_elements );
        }
        const auto j_plus_one = j+1;
        if (j_plus_one < _owned_range_j.second) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i,j_plus_one,k)) -
                              owned_first );
        } else if (j_plus_one < _remote_range_j.second) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i,j_plus_one,k)) -
                              remote_first +
                              num_owned_elements );
        }
      }

      // sides on k
      {
        const auto k_minus_one = k-1;
        if (k_minus_one >= _owned_range_k.first) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i,j,k_minus_one)) -
                              owned_first );
        } else if (k_minus_one >= _remote_range_k.first) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i,j,k_minus_one)) -
                              remote_first +
                              num_owned_elements );
        }
        const auto k_plus_one = k+1;
        if (k_plus_one < _owned_range_k.second) {
          _colidx(lcnt++) = ( lower_bound(owned_first, owned_last, _sb.ijk_to_idx(i,j,k_plus_one)) -
                              owned_first );
        } else if (k_plus_one < _remote_range_k.second) {
          _colidx(lcnt++) = ( lower_bound(remote_first, remote_last, _sb.ijk_to_idx(i,j,k_plus_one)) -
                              remote_first +
                              num_owned_elements );
        }
      }

      // sort
      heap_sort(&_colidx(lbeg), lcnt-lbeg);
    }

    inline
    void run() {
      Kokkos::parallel_for(_owned_gids.extent(0), *this);
    }

  };

}

#endif


