// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef __TPETRA_TEST_BLOCKCRS_MESHDATABASE_HPP__
#define __TPETRA_TEST_BLOCKCRS_MESHDATABASE_HPP__

#include <iostream>
#include <sstream>

#include <Teuchos_Comm.hpp>
#include <Tpetra_DefaultPlatform.hpp>
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

  typedef Kokkos::DefaultExecutionSpace exec_space;
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
    T2 largest = i;
    T2 l = 2*i + 1;
    T2 r = 2*i + 2;

    if (l < n && v[l] > v[largest]) largest = l;
    if (r < n && v[r] > v[largest]) largest = r;
    if (largest != i) {
      // swap
      T1 tmp = v[i]; v[i] = v[largest]; v[largest] = tmp;
      heapify(v, n, largest);
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
      LO _num_procs_i, _num_procs_j, _num_procs_k, _num_procs_jk, _max_num_procs;
      LO _rank, _proc_i, _proc_j, _proc_k;
      StructuredProcGrid() = default;
      StructuredProcGrid(const StructuredProcGrid &b) = default;
      StructuredProcGrid(const LO num_procs_i, 
                         const LO num_procs_j, 
                         const LO num_procs_k) 
        : _num_procs_i(num_procs_i), 
          _num_procs_j(num_procs_j), 
          _num_procs_k(num_procs_k), 
          _num_procs_jk(num_procs_j*num_procs_k) {
        const LO bigger_ij = num_procs_i > num_procs_j ? num_procs_i : num_procs_j;
        const LO bigger_jk = num_procs_j > num_procs_k ? num_procs_j : num_procs_k;
        _max_num_procs = bigger_ij > bigger_jk ? bigger_ij : bigger_jk;
      }

      void setRank(const LO rank) {
        _rank = rank;
        _proc_i = _rank / _num_procs_jk;
        _proc_k = _rank % _num_procs_jk;
        _proc_j = _proc_k / _num_procs_k;
        _proc_k = _proc_k % _num_procs_k;
      }
    };
    
    struct StructuredBlock {
    public:
      LO _num_global_elements_i, 
        _num_global_elements_j, 
        _num_global_elements_k, 
        _num_global_elements_jk;
      StructuredBlock() = default;
      StructuredBlock(const StructuredBlock &b) = default;
      StructuredBlock(const LO num_global_elements_i,
                      const LO num_global_elements_j,
                      const LO num_global_elements_k) 
        : _num_global_elements_i(num_global_elements_i),
          _num_global_elements_j(num_global_elements_j),
          _num_global_elements_k(num_global_elements_k),
          _num_global_elements_jk(num_global_elements_j*num_global_elements_k){}
      
      inline GO ijk_to_idx(const LO i, const LO j, const LO k) const {
        return (i*_num_global_elements_j + j)*_num_global_elements_k + k;
      }
    
      inline void idx_to_ijk(const GO idx, LO &i, LO &j, LO &k) const {
        i = idx / _num_global_elements_jk;
        k = idx % _num_global_elements_jk;
        j = k   / _num_global_elements_k;
        k = k   % _num_global_elements_k;
      }
    
      inline GO getNumElements() const { 
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
    Teuchos::RCP<const Teuchos::Comm<int> > _comm;
    StructuredBlock _sb;
    StructuredProcGrid _grid;
    StructuredBlockPart _owned;

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
                 LO num_procs_i, 
                 LO num_procs_j,
                 LO num_procs_k) 
      : _comm(comm), 
        _sb(num_global_elements_i, num_global_elements_j, num_global_elements_k),
        _grid(num_procs_i, num_procs_j, num_procs_k) {
      _grid.setRank(_comm->getRank());
    
      // uniform partitoins on the structured block
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

}

#endif


