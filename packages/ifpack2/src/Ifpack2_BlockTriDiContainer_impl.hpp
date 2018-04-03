/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_BLOCKTRIDICONTAINER_IMPL_HPP
#define IFPACK2_BLOCKTRIDICONTAINER_IMPL_HPP

#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_BlockMultiVector.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>

#include <memory>

namespace Ifpack2 {

  namespace BlockTriDiContainerDetails {
    ///
    /// view decorators for unmanaged and const memory
    ///
    template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
    using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
                                              MemoryTraitsType::RandomAccess |
                                              MemoryTraitsType::Atomic |
                                              flag>;
    
    template <typename ViewType>
    using Unmanaged = Kokkos::View<typename ViewType::data_type,
                                   typename ViewType::array_layout,
                                   typename ViewType::device_type,
                                   MemoryTraits<typename ViewType::memory_traits,
                                                Kokkos::Unmanaged> >;
    template <typename ViewType>
    using Const = Kokkos::View<typename ViewType::const_data_type, 
                               typename ViewType::array_layout,
                               typename ViewType::device_type, 
                               typename ViewType::memory_traits>;
    template <typename ViewType>
    using ConstUnmanaged = Const<Unmanaged<ViewType> >;    
    
    ///
    /// utility functions
    ///
    template <typename ViewType>
    typename ViewType::HostMirror 
    create_host_mirror_view_and_sync(const ViewType& v) {
      const auto hv = Kokkos::create_mirror_view(v);
      Kokkos::deep_copy(hv, v);
      return hv;
    }
    
    // Consolidate error output prefix.
    template<typename CommPtrType>
    std::string get_msg_prefix (const CommPtrType &comm) {
      const auto rank = comm->getRank();
      const auto nranks = comm->getSize();
      std::stringstream ss;
      ss << "Rank " << rank << " of " << nranks << ": ";
      return ss.str();
    }

    // this is used for host only
    template<typename T, int N>
    struct ArrayValueType {
      T v[N];
      inline ArrayValueType() = default;
      inline ArrayValueType(const ArrayValueType &b) {
        for (int i=0;i<N;++i) 
          this->v[i] = b.v[i];
      }      
    };
    template<typename T, int N>
    inline static volatile ArrayValueType<T,N>& 
    operator+=(volatile ArrayValueType<T,N> &a, 
               volatile const ArrayValueType<T,N> &b) {
      for (int i=0;i<N;++i) 
        a.v[i] += b.v[i];
      return a;
    }
    template<typename T, int N>
    inline static ArrayValueType<T,N>& 
    operator+=(ArrayValueType<T,N> &a, 
               const ArrayValueType<T,N> &b) {
      for (int i=0;i<N;++i) 
        a.v[i] += b.v[i];
      return a;
    }
    
    template<typename T, int N>
    struct SumReducer {
      typedef SumReducer reducer;
      typedef ArrayValueType<T,N> value_type;
      typedef Kokkos::View<value_type,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
      value_type *value;
      inline SumReducer() = default;
      inline SumReducer(const SumReducer &b) = default;
      inline SumReducer(value_type &val) : value(&val) {}
      inline void join(value_type &dst, value_type &src) const {
        for (int i=0;i<N;++i) 
          dst.v[i] += src.v[i];
      }
      inline void join(volatile value_type &dst, const volatile value_type &src) const {
        for (int i=0;i<N;++i) 
          dst.v[i] += src.v[i];
      }          
      inline void init(value_type &val) const {
        for (int i=0;i<N;++i)         
          val.v[i] = Kokkos::reduction_identity<T>::sum();
      }
      inline value_type& reference() {
        return *value;
      }
      inline result_view_type view() const {
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

      ///
      /// tpetra types
      ///
      typedef typename node_type::device_type device_type;      
      typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_multivector_type;
      typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> tpetra_map_type;
      typedef Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> tpetra_import_type;
      typedef Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_row_matrix_type;
      typedef Tpetra::BlockCrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_crs_matrix_type;
      typedef Tpetra::BlockMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_multivector_type;
      typedef typename tpetra_block_crs_matrix_type::crs_graph_type::local_graph_type local_crs_graph_type;

      ///
      /// simd vectorization
      ///
      template<typename T, int l> using Vector = KokkosBatched::Experimental::Vector<T,l>;
      template<typename T> using SIMD = KokkosBatched::Experimental::SIMD<T>;
      template<typename T, typename M> using DefaultVectorLength = KokkosBatched::Experimental::DefaultVectorLength<T,M>;

      enum : int { vector_length = DefaultVectorLength<impl_scalar_type,typename device_type::memory_space>::value };
      typedef Vector<SIMD<impl_scalar_type>,vector_length> vector_type;

      ///
      /// commonly used view types 
      ///
      typedef Kokkos::View<size_type*,device_type> size_type_1d_view;
      typedef Kokkos::View<local_ordinal_type*,device_type> local_ordinal_type_1d_view;
      // tpetra block crs values
      typedef Kokkos::View<impl_scalar_type*,device_type> impl_scalar_type_1d_view;
      // tpetra multivector values (layout left): may need to change the typename more explicitly 
      typedef Kokkos::View<impl_scalar_type**,Kokkos::LayoutLeft,device_type> impl_scalar_type_2d_view;
      // packed data always use layout right
      typedef Kokkos::View<vector_type*,device_type> vector_type_1d_view;
      typedef Kokkos::View<vector_type***,Kokkos::LayoutRight,device_type> vector_type_3d_view;
      typedef Kokkos::View<impl_scalar_type****,Kokkos::LayoutRight,device_type> impl_scalar_type_4d_view;


      ///
      /// cuda specialization
      ///
      template<typename ExecSpace>
      struct is_cuda {
        enum : bool { value = false }; 
      };
#if defined(HAVE_IFPACK2_CUDA) && defined (KOKKOS_ENABLE_CUDA)
      template<>
      struct is_cuda<Kokkos::Cuda> {
        enum : bool { value = true };         
      };
#endif
      
    };

    ///
    /// setup sequential importer
    ///
    template<typename MatrixType>
    typename Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_import_type> 
    createBlockCrsTpetraImporter
    (const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A) {
      using map_type = typename ImplType<MatrixType>::tpetra_map_type;
      using mv_type = typename ImplType<MatrixType>::tpetra_block_multivector_type;
      
      const auto g = A->getCrsGraph();  // tpetra crs graph object
      const auto blocksize = A->getBlockSize();
      const auto src = Teuchos::rcp(new map_type(mv_type::makePointMap(*g.getDomainMap(), blocksize)));
      const auto tgt = Teuchos::rcp(new map_type(mv_type::makePointMap(*g.getColMap()   , blocksize)));

      return Teuchos::rcp(new typename ImplType<MatrixType>::tpetra_import_type(src, tgt));
    }

    template<typename MatrixType>
    struct PartInterface {
      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;

      PartInterface() = default;
      PartInterface(const PartInterface &b) = default;

      // Some terms:
      //   The matrix A is split as A = D + R, where D is the matrix of tridiag
      // blocks and R is the remainder.
      //   A part is roughly a synonym for a tridiag. The distinction is that a part
      // is the set of rows belonging to one tridiag and, equivalently, the off-diag
      // rows in R associated with that tridiag. In contrast, the term tridiag is
      // used to refer specifically to tridiag data, such as the pointer into the
      // tridiag data array.
      //   Local (lcl) row arge the LIDs. lclrow lists the LIDs belonging to each
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
      local_ordinal_type_1d_view partptr; // np+1
      // packptr_(i), for i the pack index, indexes partptr_. partptr_(packptr_(i))
      // is the start of the i'th pack.
      local_ordinal_type_1d_view packptr; // npack+1
      // part2rowidx0_(i) is the flat row index of the start of the i'th part. It's
      // an alias of partptr_ in the case of no overlap.
      local_ordinal_type_1d_view part2rowidx0; // np+1
      // part2packrowidx0_(i) is the packed row index. If vector_length is 1, then
      // it's the same as part2rowidx0_; if it's > 1, then the value is combined
      // with i % vector_length to get the location in the packed data.
      local_ordinal_type_1d_view part2packrowidx0; // np+1
      local_ordinal_type part2packrowidx0_back; // So we don't need to grab the array from the GPU.
      // rowidx2part_ maps the row index to the part index.
      local_ordinal_type_1d_view rowidx2part; // nr
      // If needed, permute owned domain main LIDs to owned column map LIDs.
      local_ordinal_type_1d_view dm2cm;
      // True if lcl{row|col} is at most a constant away from row{idx|col}. In
      // practice, this knowledge is not particularly useful, as packing for batched
      // processing is done at the same time as the permutation from LID to index
      // space. But it's easy to detect, so it's recorded in case an optimization
      // can be made based on it.
      bool row_contiguous;
    };

    ///
    /// setup part interface using the container partitions array
    ///
    template<typename MatrixType>
    PartInterface<MatrixType> 
    createPartInterface(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
                        const Teuchos::Array<Teuchos::Array<typename ImplType<MatrixType>::local_ordinal_type> > &partitions) {
      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;

      enum : int { vector_length = ImplType<MatrixType>::vector_length };

      const auto comm = A->getRowMap()->getComm();
      
      PartInterface<MatrixType> interf;
      
      const bool jacobi = partitions.size() == 0;
      const local_ordinal_type A_n_lclrows = A->getNodeNumRows();
      const local_ordinal_type nparts = jacobi ? A_n_lclrows : partitions.size();

#if defined(BLOCKTRIDICONTAINER_DEBUG)              
      local_ordinal_type nrows = 0;
      if (jacobi)       nrows = nparts;
      else              for (local_ordinal_type i=0;i<nparts;++i) nrows += partitions[i].size();

      TEUCHOS_TEST_FOR_EXCEPT_MSG
        (nrows != A_n_lclrows, get_msg_prefix(comm) << "The #rows implied by the local partition is not "
         << "the same as getNodeNumRows: " << nrows << " vs " << A_n_lclrows);
#endif
      
      // permutation vector
      std::vector<local_ordinal_type> p;
      if (!jacobi) {
        // reorder parts to maximize simd packing efficiency
        const local_ordinal_type nparts = partitions.size();
        p.resize(nparts);
        
        typedef std::pair<local_ordinal_type,local_ordinal_type> size_idx_pair_type;
        std::vector<size_idx_pair_type> partsz(nparts);
        for (local_ordinal_type i=0;i<nparts;++i)
          partsz[i] = size_idx_pair_type(partitions[i].size(), i);
        std::sort(partsz.begin(), partsz.end(),
                  [] (const size_idx_pair_type& x, const size_idx_pair_type& y) { 
                    return x.first > y.first; 
                  });
        for (local_ordinal_type i=0;i<nparts;++i)
          p[i] = partsz[i].second;
      }

      // allocate parts
      interf.partptr = local_ordinal_type_1d_view("partptr", nparts + 1);
      interf.lclrow = local_ordinal_type_1d_view("lclrow", A_n_lclrows);
      interf.part2rowidx0 = local_ordinal_type_1d_view("part2rowidx0", nparts + 1);
      interf.part2packrowidx0 = local_ordinal_type_1d_view("part2packrowidx0", nparts + 1);
      interf.rowidx2part = local_ordinal_type_1d_view("rowidx2part", A_n_lclrows);

      // mirror to host and compute on host execution space
      const auto partptr = Kokkos::create_mirror_view(interf.partptr);
      const auto lclrow = Kokkos::create_mirror_view(interf.lclrow);
      const auto part2rowidx0 = Kokkos::create_mirror_view(interf.part2rowidx0);
      const auto part2packrowidx0 = Kokkos::create_mirror_view(interf.part2packrowidx0);
      const auto rowidx2part = Kokkos::create_mirror_view(interf.rowidx2part);
      
      // Determine parts.
      interf.row_contiguous = true;
      partptr(0) = 0;
      part2rowidx0(0) = 0;
      part2packrowidx0(0) = 0;
      local_ordinal_type pack_nrows = 0;
      for (local_ordinal_type ip=0;ip<nparts;++ip) {
        const auto* part = jacobi ? NULL : &partitions[p[ip]];
        const local_ordinal_type ipnrows = jacobi ? 1 : part->size();
        TEUCHOS_ASSERT(ip == 0 || (jacobi || ipnrows <= static_cast<local_ordinal_type>(partitions[p[ip-1]].size())));
        TEUCHOS_TEST_FOR_EXCEPT_MSG(ipnrows == 0, 
                                    get_msg_prefix(comm) 
                                    << "partition " << p[ip]
                                    << " is empty, which is not allowed.");
        //assume No overlap.
        part2rowidx0(ip+1) = part2rowidx0(ip) + ipnrows;
        // Since parts are ordered in nonincreasing size, the size of the first
        // part in a pack is the size for all parts in the pack.
        if (ip % vector_length == 0) pack_nrows = ipnrows;
        part2packrowidx0(ip+1) = part2packrowidx0(ip) + ((ip+1) % vector_length == 0 || ip+1 == nparts ? pack_nrows : 0);
        const local_ordinal_type os = partptr(ip);
        for (local_ordinal_type i=0;i<ipnrows;++i) {
          const auto lcl_row = jacobi ? ip : (*part)[i];
          TEUCHOS_TEST_FOR_EXCEPT_MSG(lcl_row < 0 || lcl_row >= A_n_lclrows, 
                                      get_msg_prefix(comm) 
                                      << "partitions[" << p[ip] << "]["
                                      << i << "] = " << lcl_row 
                                      << " but input matrix implies limits of [0, " << A_n_lclrows-1
                                      << "].");
          lclrow(os+i) = lcl_row;
          rowidx2part(os+i) = ip;
          if (interf.row_contiguous && os+i > 0 && lclrow((os+i)-1) + 1 != lcl_row)
            interf.row_contiguous = false;
        }
        partptr(ip+1) = os + ipnrows;
      }
#if defined(BLOCKTRIDICONTAINER_DEBUG)              
      TEUCHOS_ASSERT(partptr(nparts) == nrows);
#endif
      if (lclrow(0) != 0) interf.row_contiguous = false;

      Kokkos::deep_copy(interf.partptr, partptr);
      Kokkos::deep_copy(interf.lclrow, lclrow);

      //assume No overlap. Thus:
      interf.part2rowidx0 = interf.partptr;
      Kokkos::deep_copy(interf.part2packrowidx0, part2packrowidx0);

      interf.part2packrowidx0_back = part2packrowidx0(part2packrowidx0.extent(0) - 1);
      Kokkos::deep_copy(interf.rowidx2part, rowidx2part);

      { // Fill packptr.
        local_ordinal_type npacks = 0;
        for (local_ordinal_type ip=1;ip<=nparts;++ip)
          if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
            ++npacks;
        interf.packptr = local_ordinal_type_1d_view("packptr", npacks + 1);
        const auto packptr = Kokkos::create_mirror_view(interf.packptr);
        packptr(0) = 0;
        for (local_ordinal_type ip=1,k=1;ip<=nparts;++ip)
          if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
            packptr(k++) = ip;
        Kokkos::deep_copy(interf.packptr, packptr);
      }

      return interf;
    }

    ///
    /// block tridiagonals
    ///
    template <typename MatrixType>
    struct BlockTridiags {
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view;
      using vector_type_3d_view = typename ImplType<MatrixType>::vector_type_3d_view;

      // flat_td_ptr(i) is the index into flat-array values of the start of the
      // i'th tridiag. pack_td_ptr is the same, but for packs. If vector_length ==
      // 1, pack_td_ptr is the same as flat_td_ptr; if vector_length > 1, then i %
      // vector_length is the position in the pack.
      size_type_1d_view flat_td_ptr, pack_td_ptr;
      // List of local column indices into A from which to grab
      // data. flat_td_ptr(i) points to the start of the i'th tridiag's data.
      local_ordinal_type_1d_view A_colindsub;
      // Tridiag block values. pack_td_ptr(i) points to the start of the i'th
      // tridiag's pack, and i % vector_length gives the position in the pack.
      vector_type_3d_view values;

      // Index into row-major block of a tridiag.
      template <typename idx_type> 
      static KOKKOS_FORCEINLINE_FUNCTION
      idx_type IndexToRow (const idx_type& ind) { return (ind + 1) / 3; }
      // Given a row of a row-major tridiag, return the index of the first block
      // in that row.
      template <typename idx_type> 
      static KOKKOS_FORCEINLINE_FUNCTION
      idx_type RowToIndex (const idx_type& row) { return row > 0 ? 3*row - 1 : 0; }
      // Number of blocks in a tridiag having a given number of rows.
      template <typename idx_type> 
      static KOKKOS_FORCEINLINE_FUNCTION
      idx_type NumBlocks (const idx_type& nrows) { return nrows > 0 ? 3*nrows - 2 : 0; }
    };

    
    ///
    /// block tridiags initialization from part interface
    /// 
    template<typename MatrixType>
    BlockTridiags<MatrixType> 
    createBlockTridiags(const PartInterface<MatrixType> &interf) {
      using device_type = typename ImplType<MatrixType>::device_type;
      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using size_type = typename ImplType<MatrixType>::size_type;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view;

      enum : int { vector_length = ImplType<MatrixType>::vector_length };

      BlockTridiags<MatrixType> btdm;

      const auto partptr = create_host_mirror_view_and_sync(interf.partptr);
      const local_ordinal_type ntridiags = interf.partptr.extent(0) - 1;
      
      { // construct the flat index pointers into the tridiag values array.
        btdm.flat_td_ptr = size_type_1d_view("btdm.flat_td_ptr", ntridiags + 1);
        const auto flat_td_ptr = Kokkos::create_mirror_view(btdm.flat_td_ptr);
        flat_td_ptr(0) = 0;
        for (local_ordinal_type ti = 1; ti <= ntridiags; ++ti) {
          const local_ordinal_type nrows = partptr(ti) - partptr(ti-1);
          flat_td_ptr(ti) = flat_td_ptr(ti-1) + (3*nrows - 2);
        }
#if defined(BLOCKTRIDICONTAINER_DEBUG)        
        {
          const size_type nnz = 3*partptr(partptr.size() - 1) - 2*ntridiags;
          TEUCHOS_ASSERT(flat_td_ptr(ntridiags) == nnz);
        }
#endif
        Kokkos::deep_copy(btdm.flat_td_ptr, flat_td_ptr);
      }
      
      // And the packed index pointers.
      if (vector_length == 1) {
        btdm.pack_td_ptr = btdm.flat_td_ptr;
      } else {
        const auto packptr = create_host_mirror_view_and_sync(interf.packptr);
        const local_ordinal_type npacks = packptr.extent(0) - 1;
        btdm.pack_td_ptr = size_type_1d_view("btdm.pack_td_ptr", ntridiags + 1);
        const auto pack_td_ptr = Kokkos::create_mirror_view(btdm.pack_td_ptr);
        size_type nblks = 0;
        for (local_ordinal_type pki = 0; pki < npacks; ++pki) {
          const local_ordinal_type parti = packptr(pki);
          for (local_ordinal_type pti = parti; pti < packptr(pki+1); ++pti)
            pack_td_ptr(pti) = nblks;
          nblks += BlockTridiags<MatrixType>::NumBlocks(partptr(parti+1) - partptr(parti));
        }
        pack_td_ptr(ntridiags) = nblks;
        Kokkos::deep_copy(btdm.pack_td_ptr, pack_td_ptr);
      }
      return btdm;
    }

    // Set the tridiags to be I to the full pack block size. That way, if a
    // tridiag within a pack is shorter than the longest one, the extra blocks are
    // processed in a safe way. Similarly, in the solve phase, if the extra blocks
    // in the packed multvector are 0, and the tridiag LU reflects the extra I
    // blocks, then the solve proceeds as though the extra blocks aren't
    // present. Since this extra work is part of the SIMD calls, it's not actually
    // extra work. Instead, it means we don't have to put checks or masks in, or
    // quiet NaNs. This functor has to be called just once, in the symbolic phase,
    // since the numeric phase fills in only the used entries, leaving these I
    // blocks intact.
    template<typename MatrixType>
    void setTridiagsToIdentity(const BlockTridiags<MatrixType> &btdm, 
                               const typename ImplType<MatrixType>::local_ordinal_type_1d_view &packptr) {
      using device_type = typename ImplType<MatrixType>::device_type;
      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using size_type = typename ImplType<MatrixType>::size_type;
      
      const local_ordinal_type blocksize = btdm.values.extent(1);
      Kokkos::parallel_for
        (Kokkos::RangePolicy<typename device_type::execution_space>(0, packptr.extent(0) - 1),
         KOKKOS_LAMBDA(const local_ordinal_type k) {
          for (size_type i=btdm.pack_td_ptr(packptr(k)),iend=btdm.pack_td_ptr(packptr(k+1));i<iend;i += 3)
            for (local_ordinal_type j=0;j<blocksize;++j)
              btdm.values(i,j,j) = 1;
        });
    }
    
    ///
    /// A - Tridiags(A), i.e., R in the splitting A = D + R.
    ///
    template <typename MatrixType>
    struct AmD {
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view;
      using impl_scalar_type_1d_view = typename ImplType<MatrixType>::impl_scalar_type_1d_view;

      // rowptr points to the start of each row of A_colindsub.
      size_type_1d_view rowptr, rowptr_remote;
      // Indices into A's rows giving the blocks to extract. rowptr(i) points to
      // the i'th row. Thus, g.entries(A_colindsub(rowptr(row) : rowptr(row+1))),
      // where g is A's graph, are the columns AmD uses. If seq_method_, then
      // A_colindsub contains all the LIDs and A_colindsub_remote is empty. If !
      // seq_method_, then A_colindsub contains owned LIDs and A_colindsub_remote
      // contains the remote ones.
      local_ordinal_type_1d_view A_colindsub, A_colindsub_remote;
      
      // Currently always true.
      bool is_tpetra_block_crs;

      // If is_tpetra_block_crs, then this is a pointer to A_'s value data.
      impl_scalar_type_1d_view tpetra_values;
    };

    ///
    /// symbolic phase, on host : create R = A - D, pack D
    ///
    template<typename MatrixType>
    void
    performSymbolicPhase
    (const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
     const PartInterface<MatrixType> &interf,
     BlockTridiags<MatrixType> &btdm,
     AmD<MatrixType> &amd,
     const bool overlap_comm) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::symbolic");
#endif
      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using global_ordinal_type = typename ImplType<MatrixType>::global_ordinal_type;
      using size_type = typename ImplType<MatrixType>::size_type;
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view;
      using vector_type_3d_view = typename ImplType<MatrixType>::vector_type_3d_view;
      using block_crs_matrix_type = typename ImplType<MatrixType>::tpetra_block_crs_matrix_type;

      using device_type = typename ImplType<MatrixType>::device_type;
      using host_execution_space = Kokkos::DefaultHostExecutionSpace;

      enum : int { vector_length = ImplType<MatrixType>::vector_length };

      const auto comm = A->getRowMap()->getComm();
      const auto& g = A->getCrsGraph();
      const auto blocksize = A->getBlockSize();      

      // mirroring to host
      const auto partptr = create_host_mirror_view_and_sync(interf.partptr);
      const auto lclrow = create_host_mirror_view_and_sync(interf.lclrow);
      const auto rowidx2part = create_host_mirror_view_and_sync(interf.rowidx2part);
      const auto part2rowidx0 = create_host_mirror_view_and_sync(interf.part2rowidx0);
      const auto packptr = create_host_mirror_view_and_sync(interf.packptr);

      const local_ordinal_type nrows = partptr(partptr.extent(0) - 1);

      // find column to row map on host
      Kokkos::View<local_ordinal_type*,host_execution_space> col2row("col2row", A->getNodeNumCols());
      Kokkos::deep_copy(col2row, Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
      {
        const auto rowmap = g.getRowMap();
        const auto colmap = g.getColMap();
        const auto dommap = g.getDomainMap();
        TEUCHOS_ASSERT( !(rowmap.is_null() || colmap.is_null() || dommap.is_null()));

        Kokkos::parallel_for
          (Kokkos::RangePolicy<host_execution_space>(0,nrows),
           KOKKOS_LAMBDA(const local_ordinal_type &lr) {
            const global_ordinal_type gid = rowmap->getGlobalElement(lr);
            TEUCHOS_ASSERT(gid != Teuchos::OrdinalTraits<global_ordinal_type>::invalid());
            if (dommap->isNodeGlobalElement(gid)) {
              const local_ordinal_type lc = colmap->getLocalElement(gid);
#if defined(BLOCKTRIDICONTAINER_DEBUG)              
              TEUCHOS_TEST_FOR_EXCEPT_MSG(lc == Teuchos::OrdinalTraits<local_ordinal_type>::invalid(), 
                                          get_msg_prefix(comm) << "GID " << gid
                                          << " gives an invalid local column.");
#endif
              col2row(lc) = lr;
            }
          });
      }
      
      // construct the D and R graphs in A = D + R.
      { 
        const auto& local_graph = g.getLocalGraph();
        const auto& local_graph_rowptr = local_graph.row_map;
        TEUCHOS_ASSERT(local_graph_rowptr.size() == static_cast<size_t>(nrows + 1));
        const auto& local_graph_colidx = local_graph.entries;

        //assume no overlap.
        Kokkos::View<local_ordinal_type*,host_execution_space> lclrow2idx("lclrow2idx", nrows);
        Kokkos::parallel_for
          (Kokkos::RangePolicy<host_execution_space>(0,nrows),
           KOKKOS_LAMBDA(const local_ordinal_type &i) {
            lclrow2idx[lclrow(i)] = i;
          });

        // count (block) nnzs in D and R.
        typedef SumReducer<size_type,3> sum_reducer_type;
        typename sum_reducer_type::value_type sum_reducer_value;
        Kokkos::parallel_reduce
          (Kokkos::RangePolicy<host_execution_space>(0,nrows),
           KOKKOS_LAMBDA(const local_ordinal_type &lr, typename sum_reducer_type::value_type &update) {
            // LID -> index.
            const local_ordinal_type ri0 = lclrow2idx[lr];
            const local_ordinal_type pi0 = rowidx2part(ri0);
            for (size_type j=local_graph_rowptr(lr);j<local_graph_rowptr(lr+1);++j) {
              const local_ordinal_type lc = local_graph_colidx(j);
              const local_ordinal_type lc2r = col2row[lc];
              bool incr_R = false;
              do { // breakable
                if (lc2r == Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) {
                  incr_R = true;
                  break;
                }
                const local_ordinal_type ri = lclrow2idx[lc2r];
                const local_ordinal_type pi = rowidx2part(ri);
                if (pi != pi0) {
                  incr_R = true;
                  break;
                }
                // Test for being in the tridiag. This is done in index space. In
                // LID space, tridiag LIDs in a row are not necessarily related by
                // {-1, 0, 1}.
                if (ri0 + 1 >= ri && ri0 <= ri + 1)
                  ++update.v[0]; // D_nnz
                else
                  incr_R = true;
              } while (0);
              if (incr_R) {
                if (lc < nrows) ++update.v[1]; // R_nnz_owned
                else ++update.v[2]; // R_nnz_remote
              }
            }
          }, sum_reducer_type(sum_reducer_value));
        
        size_type D_nnz = sum_reducer_value.v[0];
        size_type R_nnz_owned = sum_reducer_value.v[1];
        size_type R_nnz_remote = sum_reducer_value.v[2];
        
        if (!overlap_comm) {
          R_nnz_owned += R_nnz_remote;
          R_nnz_remote = 0;
        }

        // construct the D graph.
        { 
          const auto flat_td_ptr = create_host_mirror_view_and_sync(btdm.flat_td_ptr);

          btdm.A_colindsub = local_ordinal_type_1d_view("btdm.A_colindsub", D_nnz);
          const auto D_A_colindsub = Kokkos::create_mirror_view(btdm.A_colindsub);

#if defined(BLOCKTRIDICONTAINER_DEBUG)              
          Kokkos::deep_copy(D_A_colindsub, Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
#endif

          const local_ordinal_type nparts = partptr.extent(0) - 1;
          Kokkos::parallel_for
            (Kokkos::RangePolicy<host_execution_space>(0, nparts),
             KOKKOS_LAMBDA(const local_ordinal_type &pi0) {
              const local_ordinal_type part_ri0 = part2rowidx0(pi0);
              local_ordinal_type offset = 0;
              for (local_ordinal_type ri0=partptr(pi0);ri0<partptr(pi0+1);++ri0) {
                const local_ordinal_type td_row_os 
                  = BlockTridiags<MatrixType>::RowToIndex(ri0 - part_ri0) + offset;
                offset = 1;
                const local_ordinal_type lr0 = lclrow(ri0);
                const size_type j0 = local_graph_rowptr(lr0);
                for (size_type j=j0;j<local_graph_rowptr(lr0+1);++j) {
                  const local_ordinal_type lc = local_graph_colidx(j);
                  const local_ordinal_type lc2r = col2row[lc];
                  if (lc2r == Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) continue;
                  const local_ordinal_type ri = lclrow2idx[lc2r];
                  const local_ordinal_type pi = rowidx2part(ri);
                  if (pi != pi0) continue;
                  if (ri + 1 < ri0 || ri > ri0 + 1) continue;
                  const local_ordinal_type row_entry = j - j0;
                  D_A_colindsub(flat_td_ptr(pi0) + ((td_row_os + ri) - ri0)) = row_entry;
                }
              }
            });
          
#if defined(BLOCKTRIDICONTAINER_DEBUG)              
          for (size_t i=0;i<D_A_colindsub.extent(0);++i)
            TEUCHOS_ASSERT(D_A_colindsub(i) != Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
#endif
          Kokkos::deep_copy(btdm.A_colindsub, D_A_colindsub);

          // Allocate values.
          { 
            const local_ordinal_type npacks = packptr.extent(0) - 1;
            local_ordinal_type nblks = 0; // Number of tridiag blocks, accounting for packing.
            Kokkos::parallel_reduce
              (Kokkos::RangePolicy<host_execution_space>(0, npacks),
               KOKKOS_LAMBDA(const local_ordinal_type pai, local_ordinal_type &update) {
                const local_ordinal_type pti = packptr(pai);
                const local_ordinal_type inrows = partptr(pti+1) - partptr(pti);
                update += BlockTridiags<MatrixType>::RowToIndex(inrows);
              }, nblks);
            
            btdm.values = vector_type_3d_view("btdm.values", nblks, blocksize, blocksize);
            if (vector_length > 1) setTridiagsToIdentity(btdm, interf.packptr);
          }
        }

        // Construct the R graph.        
        { 
          amd.rowptr = size_type_1d_view("amd.rowptr", nrows + 1);
          amd.A_colindsub = local_ordinal_type_1d_view("amd.A_colindsub", R_nnz_owned);
          const auto R_rowptr = Kokkos::create_mirror_view(amd.rowptr);
          const auto R_A_colindsub = Kokkos::create_mirror_view(amd.A_colindsub);
          R_rowptr(0) = 0;
          if (overlap_comm) {
            amd.rowptr_remote = size_type_1d_view("amd.rowptr_remote", nrows + 1);
            amd.A_colindsub_remote = local_ordinal_type_1d_view("amd.A_colindsub_remote", R_nnz_remote);
          }
          const auto R_rowptr_remote = Kokkos::create_mirror_view(amd.rowptr_remote);
          const auto R_A_colindsub_remote = Kokkos::create_mirror_view(amd.A_colindsub_remote);
          if (overlap_comm) R_rowptr_remote(0) = 0;

          Kokkos::parallel_for
            (Kokkos::RangePolicy<host_execution_space>(0,nrows),
             KOKKOS_LAMBDA(const local_ordinal_type &lr) {
              const local_ordinal_type ri0 = lclrow2idx[lr];
              const local_ordinal_type pi0 = rowidx2part(ri0);
              const size_type j0 = local_graph_rowptr(lr);
              for (size_type j=j0;j<local_graph_rowptr(lr+1);++j) {
                const local_ordinal_type lc = local_graph_colidx(j);
                const local_ordinal_type lc2r = col2row[lc];
                if (lc2r != Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) {
                  const local_ordinal_type ri = lclrow2idx[lc2r];
                  const local_ordinal_type pi = rowidx2part(ri);
                  if (pi == pi0 && ri + 1 >= ri0 && ri <= ri0 + 1)
                    continue;
                }
                const local_ordinal_type row_entry = j - j0;
                if (!overlap_comm || lc < nrows) 
                  ++R_rowptr(lr);
                else 
                  ++R_rowptr_remote(lr);
              }
            });
          
          // typedef ArrayValueType<size_type,2> update_type;
          // Kokkos::parallel_scan
          //   (Kokkos::RangePolicy<host_execution_space>(0,nrows+1),
          //    KOKKOS_LAMBDA(const local_ordinal_type &lr, 
          //                  update_type &update, 
          //                  const bool &final) {
          //     if (final) {
          //       R_rowptr(lr) += update.v[0];
          //       R_rowptr_remote(lr) += update.v[1];

          //       if (lr < nrows) {
          //         const local_ordinal_type ri0 = lclrow2idx[lr];
          //         const local_ordinal_type pi0 = rowidx2part(ri0);
                  
          //         size_type cnt_rowptr = R_rowptr(lr);
          //         size_type cnt_rowptr_remote = R_rowptr_remote(lr);

          //         const size_type j0 = local_graph_rowptr(lr);
          //         for (size_type j=j0;j<local_graph_rowptr(lr+1);++j) {
          //           const local_ordinal_type lc = local_graph_colidx(j);
          //           const local_ordinal_type lc2r = col2row[lc];
          //           if (lc2r != Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) {
          //             const local_ordinal_type ri = lclrow2idx[lc2r];
          //             const local_ordinal_type pi = rowidx2part(ri);
          //             if (pi == pi0 && ri + 1 >= ri0 && ri <= ri0 + 1)
          //               continue;
          //           }
          //           const local_ordinal_type row_entry = j - j0;
          //           if ( !overlap_comm || lc < nrows) 
          //             R_A_colindsub(cnt_rowptr++) = row_entry;
          //           else 
          //             R_A_colindsub_remote(cnt_rowptr_remote++) = row_entry;
          //         }
          //       }
          //     }
          //     update.v[0] += R_rowptr(lr);
          //     update.v[1] += R_rowptr_remote(lr);
          //   });

          TEUCHOS_ASSERT(R_rowptr(nrows) == R_nnz_owned);
          Kokkos::deep_copy(amd.rowptr, R_rowptr);
          Kokkos::deep_copy(amd.A_colindsub, R_A_colindsub);
          if (overlap_comm) {
            TEUCHOS_ASSERT(R_rowptr_remote(nrows) == R_nnz_remote);
            Kokkos::deep_copy(amd.rowptr_remote, R_rowptr_remote);
            Kokkos::deep_copy(amd.A_colindsub_remote, R_A_colindsub_remote);
          }

          // Allocate or view values.
          amd.tpetra_values = (const_cast<block_crs_matrix_type*>(A.get())->
                               template getValues<typename device_type::memory_space>());
        }
      }
    }

    ///
    /// numeric phase, initialize the preconditioner
    ///
    template<typename MatrixType>
    class ExtractAndFactorizeTridiags {
    public:
      using device_type = typename ImplType<MatrixType>::device_type;

      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using size_type = typename ImplType<MatrixType>::size_type;
      using impl_scalar_type = typename ImplType<MatrixType>::impl_scalar_type;
      using magnitude_type = typename ImplType<MatrixType>::magnitude_type;
      /// tpetra interface
      using block_crs_matrix_type = typename ImplType<MatrixType>::tpetra_block_crs_matrix_type;
      /// views
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view; 
      using impl_scalar_type_1d_view = typename ImplType<MatrixType>::impl_scalar_type_1d_view; 
      /// vectorization 
      using vector_type_3d_view = typename ImplType<MatrixType>::vector_type_3d_view;
      enum : int { vector_length = ImplType<MatrixType>::vector_length };
      /// flat view to vector
      //using impl_scalar_type_4d_view = typename ImplType<MatrixType>::impl_scalar_type_4d_view;
      /// team policy member type (used in cuda)
      using member_type = typename Kokkos::TeamPolicy<device_type>::member_type;      
 
    private:
      // part interface
      ConstUnmanaged<local_ordinal_type_1d_view> partptr, lclrow, packptr;
      // block crs matrix
      ConstUnmanaged<size_type_1d_view> A_rowptr;
      ConstUnmanaged<impl_scalar_type_1d_view> A_values;
      // block tridiags 
      ConstUnmanaged<size_type_1d_view> D_pack_td_ptr, D_flat_td_ptr;
      ConstUnmanaged<local_ordinal_type_1d_view> D_colindsub;
      // block tridiags values
      Unmanaged<vector_type_3d_view> D_vector_values;
      //Unmanaged<impl_scalar_type_4d_view> D_scalar_values;
      // diagonal safety
      const magnitude_type tiny;
      // shared information
      const local_ordinal_type blocksize, blocksize_square;

    public:
      ExtractAndFactorizeTridiags(const BlockTridiags<MatrixType> &btdm_, 
                                  const PartInterface<MatrixType> &interf_,
                                  const Teuchos::RCP<block_crs_matrix_type> &A_,
                                  const magnitude_type& tiny_) : 
        // interface
        partptr(interf_.partptr), 
        lclrow(interf_.lclrow), 
        packptr(interf_.packptr),
        // block crs matrix
        A_rowptr(A_->getLocalGraph().row_map), 
        A_values(A_->template getValues<typename device_type::memory_space>()),
        // block tridiags 
        D_flat_td_ptr(btdm_.flat_td_ptr), 
        D_pack_td_ptr(btdm_.pack_td_ptr), 
        D_vector_values(btdm_.values),
        // D_scalar_values(impl_scalar_type_4d_view(btdm_.values.data(), 
        //                                          btdm_.values.extent(0),
        //                                          btdm_.values.extent(1),
        //                                          btdm_.values.extent(2),
        //                                          vector_length)),
        blocksize(btdm_.values.extent(1)),
        blocksize_square(blocksize*blocksize),
        tiny(tiny_) {}
      
      KOKKOS_INLINE_FUNCTION 
      void 
      extract(const local_ordinal_type &packidx) const {
        local_ordinal_type partidx = packptr(packidx);
        local_ordinal_type npacks = packptr(packidx+1) - partidx;

        const size_type kps = pack_td_ptr(partidx);
        local_ordinal_type kfs[vector_length] = {};
        local_ordinal_type ri0[vector_length] = {};
        local_ordinal_type nrows[vector_length] = {};

        for (local_ordinal_type vi=0;vi<npacks;++vi,++partidx) {
          kfs[vi] = D_flat_td_ptr(partidx);
          ri0[vi] = partptr(partidx);
          nrows[vi] = partptr(partidx+1) - ri0[vi];
        }
        for (local_ordinal_type tr=0,j=0;tr<nrows[0];++tr) {
          for (local_ordinal_type e=0;e<3;++e) {
            const impl_scalar_type* block[vector_length] = {};
            for (local_ordinal_type vi=0;vi<npacks;++vi) {
              const size_type Aj = A_rowptr(lclrow(ri0[vi] + tr)) + D_A_colindsub(kfs[vi] + j);
              block[vi] = &A_values(Aj*blocksize_square);
            }
            const size_type pi = kps + j;
            ++j;
            for (local_ordinal_type ii=0;ii<blocksize;++ii) {
              for (local_ordinal_type jj=0;jj<blocksize;++jj) {
                const auto idx = ii*blocksize + jj;
                auto& v = D_vector_values(pi, ii, jj);
                for (local_ordinal_type vi=0;vi<npacks;++vi)
                  v[vi] = block[vi][idx];
              }
            }

            if (nrows[0] == 1) break;
            if (e == 1 && (tr == 0 || tr+1 == nrows[0])) break;
            for (local_ordinal_type vi=1;vi<npacks;++vi) {
              if ((e == 0 && nrows[vi] == 1) || (e == 1 && tr+1 == nrows[vi])) {
                npacks = vi;
                break;
              }
            }
          }
        }
      }
      
      KOKKOS_INLINE_FUNCTION 
      void 
      factorize (const local_ordinal_type& packidx) const {
        using namespace KokkosBatched::Experimental;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        
        // subview pattern
        auto A = Kokkos::subview(D_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        auto B = A;
        auto C = A;
        
        auto i0 = pack_td_ptr(packptr(packidx));
        const local_ordinal_type nrows 
          = BlockTridiags<MatrixType>::IndexToRow(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;
        
        A.assign_data( &D_vector_values(i0,0,0) );
        SerialLU<Algo::LU::Unblocked>::invoke(A, tiny);
        for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
          B.assign_data( &D_vector_values(i0+1,0,0) );
          SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsm::Blocked>
            ::invoke(one, A, B);
          C.assign_data( &D_vector_values(i0+2,0,0) );
          SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsm::Blocked>
            ::invoke(one, A, C);
          A.assign_data( &D_vector_values(i0+3,0,0) );
          SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>
            ::invoke(-one, C, B, one, A);
          SerialLU<Algo::LU::Unblocked>::invoke(A, tiny);
        }
      }
      
    public:
      ///
      /// host serial (vector intrinsic) vectorization
      ///
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const local_ordinal_type& packidx) const {
        extract(packidx);
        factorize(packidx);
      }
      
      ///
      /// cuda team parallel vectorization
      ///
      // KOKKOS_INLINE_FUNCTION 
      // void 
      // operator() (const member_type &member) const {
      //   const local_ordinal_type packidx = member.league_rank();
      //   Kokkos::parallel_For
      //     (Kokkos::ThreadVectorRange(member, vector_length),
      //      [&](const int idx) {
      //       extract(member, idx, packidx);
      //       factorize(member, idx, packidx);
      //     });
      // }
      
      void run() {
        // #if defined(HAVE_IFPACK2_CUDA) && defined (KOKKOS_ENABLE_CUDA)
        //         Kokkos::TeamPolicy<device_type> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vector_length);
        // #else
        Kokkos::RangePolicy<device_type> policy(0, packptr.extent(0) - 1);
        //#endif
        Kokkos::parallel_for(policy, *this);
      }
    }; 
    
    ///
    /// top level numeric interface
    ///
    template<typename MatrixType>
    void
    performNumericPhase
    (const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
     const PartInterface<MatrixType> &interf,
     BlockTridiags<MatrixType> &btdm,
     const typename ImplType<MatrixType>::magnitude_type tiny) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::numeric");
#endif
      ExtractAndFactorizeTridiags<MatrixType> function(btdm, interf, A, tiny);
      function.run();
    }

    ///
    /// pack multivector
    ///
    template<typename MatrixType>
    struct MultiVectorConverter {
    public:
      using impl_type = ImplType<MatrixType>;
      using device_type = typename impl_type::device_type;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using tpetra_multivector_type = typename impl_type::tpetra_multivector_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using impl_scalar_type_2d_view = typename impl_type::impl_scalar_type_2d_view;
      enum : int { vector_length = ImplType<MatrixType>::vector_length };

    private:
      // part interface
      ConstUnmanaged<local_ordinal_type_1d_view> partptr;
      ConstUnmanaged<local_ordinal_type_1d_view> packptr;
      ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
      ConstUnmanaged<local_ordinal_type_1d_view> part2rowidx0;
      ConstUnmanaged<local_ordinal_type_1d_view> lclrow;
      local_ordinal_type blocksize;
      local_ordinal_type num_vectors;

    public:
      // packed multivector output (or input)
      Unmanaged<vector_type_3d_view> packed_multivector;
      Unmanaged<impl_scalar_type_2d_view> scalar_multivector;
      
      MultiVectorConverter(const PartInterface<MatrixType> &interf,                    
                           const vector_type_3d_view &pmv) 
        : partptr(interf.partptr),
          packptr(interf.packptr),
          part2packrowidx0(interf.part2packrowidx0),
          part2rowidx0(interf.part2rowidx0),
          lclrow(interf.lclrow),
          packed_multivector(pmv),
          blocksize(pmv.extent(1)),
          num_vectors(pmv.extent(2)) {}
      
      KOKKOS_INLINE_FUNCTION
      void copy_to_packed_multivector(const local_ordinal_type &j, 
                                      const local_ordinal_type &vi, 
                                      const local_ordinal_type &pri, 
                                      const local_ordinal_type &nrow,  
                                      const local_ordinal_type &ri0) const {
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          for (local_ordinal_type i=0;i<blocksize;++i)
            packed_multivector(pri, i, col)[vi] = scalar_multivector(blocksize*lclrow(ri0+j)+i,col);
        }
      }
      
      KOKKOS_INLINE_FUNCTION
      void 
      operator() (const local_ordinal_type &packidx) const {
        local_ordinal_type partidx = packptr(packidx);
        local_ordinal_type npacks = packptr(packidx+1) - partidx;
        const local_ordinal_type pri0 = part2packrowidx0(partidx);

        local_ordinal_type ri0[vector_length] = {};
        local_ordinal_type nrows[vector_length] = {};
        for (local_ordinal_type vi=0;vi<npacks;++vi,++partidx) {
          ri0[vi] = part2rowidx0(partidx);
          nrows[vi] = part2rowidx0(partidx+1) - ri0[vi];
        }
        for (local_ordinal_type j=0;j<nrows[0];++j) {          
          local_ordinal_type cnt = 1;
          for (;cnt<npacks && j!= nrows[cnt];++cnt);
          npacks = cnt;
          const local_ordinal_type pri = pri0 + j;
          // team parallel later
          for (local_ordinal_type vi=0;vi<npacks;++vi) 
            copy_to_packed_multivector(j, vi, pri, nrows[vi], ri0[vi]);
        }
      }

      template<typename TpetraLocalViewType>
      void to_packed_multivector(const TpetraLocalViewType &scalar_multivector_) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::PermuteAndRepack");
#endif
        scalar_multivector = scalar_multivector_;
        const Kokkos::RangePolicy<typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
        Kokkos::parallel_for(policy, *this);
      }
    };


    ///
    /// solve tridiags
    ///
    template<typename MatrixType>
    struct SolveTridiags {
    public:
      using device_type = typename ImplType<MatrixType>::device_type;

      using local_ordinal_type = typename ImplType<MatrixType>::local_ordinal_type;
      using size_type = typename ImplType<MatrixType>::size_type;
      //using impl_scalar_type = typename ImplType<MatrixType>::impl_scalar_type;
      using magnitude_type = typename ImplType<MatrixType>::magnitude_type;
      /// views
      using local_ordinal_type_1d_view = typename ImplType<MatrixType>::local_ordinal_type_1d_view;
      using size_type_1d_view = typename ImplType<MatrixType>::size_type_1d_view; 
      /// vectorization 
      using vector_type_3d_view = typename ImplType<MatrixType>::vector_type_3d_view;
      //enum : int { vector_length = ImplType<MatrixType>::vector_length };
      /// flat view to vector
      //using impl_scalar_type_4d_view = typename ImplType<MatrixType>::impl_scalar_type_4d_view;
      /// team policy member type (used in cuda)
      //using member_type = typename Kokkos::TeamPolicy<device_type>::member_type;      
 
    private:
      // part interface
      ConstUnmanaged<local_ordinal_type_1d_view> partptr;
      ConstUnmanaged<local_ordinal_type_1d_view> packptr;
      ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
      // block tridiags 
      ConstUnmanaged<size_type_1d_view> pack_td_ptr;
      // block tridiags values
      ConstUnmanaged<vector_type_3d_view> D_vector_values;
      Unmanaged<vector_type_3d_view> X_vector_values;
      //Unmanaged<impl_scalar_type_4d_view> X_scalar_values;
      //ConstUnmanaged<impl_scalar_type_4d_view> D_scalar_values;

    public:
      SolveTridiags() = default;
      SolveTridiags(const PartInterface<MatrixType> &interf,                    
                    const BlockTridiags<MatrixType> &btdm, 
                    const vector_type_3d_view &pmv) 
        :
        // interface
        partptr(interf.partptr), 
        packptr(interf.packptr),
        part2packrowidx0(interf.part2packrowidx0),
        // block tridiags
        pack_td_ptr(btdm.pack_td_ptr), 
        D_vector_values(btdm.values),
        // block multivector
        X_vector_values(pmv)
      {}

    public:
      ///
      /// host serial (vector intrinsic) vectorization
      ///
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const local_ordinal_type& packidx) const {
        // using namespace KokkosBatched::Experimental;
        
        // // constant
        // const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        
        // // subview pattern
        // auto A = Kokkos::subview(D_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        // auto B = A;
        // auto X1 = Kokkos::subview(X_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        // auto X2 = X1;

        // // index counting
        // const local_ordinal_type partidx = packptr(packidx);
        // size_type i0 = pack_td_ptr(partidx);
        // local_ordinal_type r0 = part2packrowidx0(partidx);
        // const local_ordinal_type nrows = part2packrowidx0(packptr(packidx+1)) - r0;

        // // solve Lx = x
        // A.assign_data( &D_vector_values(i0,0,0) );
        // X1.assign_data( &X_vector_values(r0,0,0) );
        // SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsv::Blocked>
        //   ::invoke(one, A, X1);
        // for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
        //   B.assign_data( &D_vector_values(i0+2,0,0) );
        //   X2.assign_data( &X_vector_values(++r0,0,0) );
        //   SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>
        //     ::invoke(-one, B, X1, one, X2);
        //   A.assign_data( &D_vector_values(i0+3,0,0) );
        //   SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsv::Blocked>
        //     ::invoke(one, A, X1);
        //   X1.assign_data( X2.data() );
        // }

        // // solve Ux = x
        // SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsv::Blocked>
        //   ::invoke(one, A, X1);
        // for (local_ordinal_type i=nrows;i>1;--i) {
        //   i0 -= 3;
        //   B.assign_data( &D_vector_values(i0+1,0,0) );
        //   X2.assign_data( &X_vector_values(--r0,0,0) );          
        //   SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>
        //     ::invoke(-one, B, X1, one, X2); 
        //   A.assign_data( &D_vector_values(i0,0,0) );          
        //   SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsv::Blocked>
        //     ::invoke(one, A, X2);
        //   X1.assign_data( X2.data() );
        // }
      }
      
      // void run() {
      //   // #if defined(HAVE_IFPACK2_CUDA) && defined (KOKKOS_ENABLE_CUDA)
      //   //         Kokkos::TeamPolicy<device_type> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vector_length);
      //   // #else
      //   Kokkos::RangePolicy<device_type> policy(local_ordinal_type(0), 
      //                                           local_ordinal_type(packptr.extent(0) - 1));
      //   //#endif
      //   //Kokkos::parallel_for(policy, *this);
      //   Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const local_ordinal_type &i) { printf("0\n"); });
      // }
    }; 

    ///
    /// top level apply interface
    ///
    template<typename MatrixType>
    int 
    applyInverseJacobi 
    (// tpetra importer
     const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_import_type> &importer,
     // tpetra interface
     const typename ImplType<MatrixType>::tpetra_multivector_type &X,  // tpetra interface
     /* */ typename ImplType<MatrixType>::tpetra_multivector_type &Y,  // tpetra interface
     /* */ typename ImplType<MatrixType>::tpetra_multivector_type &Z,  // temporary tpetra interface
     // local object interface
     const PartInterface<MatrixType> &interf, // mesh interface
     const BlockTridiags<MatrixType> &btdm, // packed block tridiagonal matrices
     /* */ typename ImplType<MatrixType>::vector_type_1d_view &work, // workspace for packed multivector of right hand side 
     //const NormManager<MatrixType> &norm_manager,
     // preconditioner parameters
     const typename ImplType<MatrixType>::impl_scalar_type &damping_factor, 
     /* */ bool is_y_zero,
     const int max_num_sweeps, 
     const typename ImplType<MatrixType>::magnitude_type tol) { // if tol is zero, then we do not compute error norm

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::applyInverseJacobi");
#endif
      using impl_type = ImplType<MatrixType>;
      using device_type = typename impl_type::device_type;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using magnitude_type = typename impl_type::magnitude_type;
      using vector_type_1d_view = typename impl_type::vector_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using tpetra_multivector_type = typename impl_type::tpetra_multivector_type;

      // input check
      TEUCHOS_TEST_FOR_EXCEPT_MSG(max_num_sweeps <= 0, "Maximum number of sweeps must be >= 1.");

      // const parameters
      const magnitude_type tolerance = tol*tol;
      const local_ordinal_type blocksize = btdm.values.extent(1);
      const local_ordinal_type num_vectors = Y.getNumVectors();
      const local_ordinal_type num_blockrows = interf.part2packrowidx0_back;

      // if workspace is needed more, resize it
      const size_type work_span_required = num_blockrows*num_vectors*blocksize;
      if (work.span() < work_span_required) 
        work = vector_type_1d_view("vector workspace 1d view", work_span_required);

      // construct copy of Y again if num vectors are different
      if (Z.getNumVectors() != num_vectors) 
        Z = tpetra_multivector_type(importer->getTargetMap(), num_vectors, false);
    
      // wrap the workspace with 3d view
      vector_type_3d_view pmv(work.data(), num_blockrows, blocksize, num_vectors);
      const auto XX = X.template getLocalView<typename device_type::memory_space>();
      const auto YY = Y.template getLocalView<device_type>();
      //const auto ZZ = Z.getLocalView<device_type>();

      MultiVectorConverter<MatrixType> multivector_converter(interf, pmv);
      SolveTridiags<MatrixType> solve_tridiags(interf, btdm, pmv);
      
      // iterate
      int sweep;
      for (sweep=0;sweep<max_num_sweeps;++sweep) {
        if (is_y_zero) {
          // pmv := x(lclrow)
          multivector_converter.to_packed_multivector(XX);
        } else {
          // y := x - R y.
          Z.doImport(Y, *importer, Tpetra::REPLACE);
          //computeBMinusRX(YY, XX, Amd, ZZ);
          // pmv := y(lclrow).
          multivector_converter.to_packed_multivector(YY);
        }

        // pmv := inv(D) pmv.
        // Kokkos::parallel_for(10, KOKKOS_LAMBDA(const local_ordinal_type i) { 
        //     solve_tridiags(i);
        //   });
        {
          Kokkos::RangePolicy<typename device_type::execution_space> policy(0,10);
          Kokkos::parallel_for(policy, solve_tridiags);
        }
        //solve_tridiags.run();

        // y(lclrow) = (b - a) y(lclrow) + a pmv, with b = 1 always.
        //packMultiVectorAndComputeNorm();
        // PermuteAndRepack(pmv_, Y, damping_factor, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_,
        //                  y_is_zero, do_norm ? norm_mgr_->get_buffer() : nullptr).run();

        // if (norm_manager.is_active) {
        //   if (sweep + 1 == max_num_sweeps) {
        //     norm_manager.ireduce(sweep, true);
        //     norm_manager.checkDone(sweep + 1, tolerance, true);
        //   } else {
        //     norm_manager.ireduce(sweep);
        //   }
        // }

        is_y_zero = false;
      }

      // sqrt the norms for the caller's use.
      //if (norm_manager.is_active) norm_manager.finalize();

      return sweep;
    }




























    

    //     // Manage the distributed part of the computation of residual norms.
    //     template<typename MatrixType>
    //     struct NormManager {
    //     public: 
    //       using ImplType<MatrixType>::magnitude_type;
      
    //     private:
    //       bool collective_;
    //       int sweep_step_, bsz_, nvec_;
    //       std::vector<magnitude_type> wrk_;
    //       magnitude_type n0_;
    // #ifdef HAVE_MPI
    //       MPI_Request mpi_request_;
    //       MPI_Comm comm_;
    // #endif

    //     public:
    //       NormManager (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
    //         sweep_step_ = 1;
    //         n0_ = 0;
    //         collective_ = comm->getSize() > 1;
    //         if (collective_) {
    // #ifdef HAVE_MPI
    //           const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    //           TEUCHOS_ASSERT( ! mpi_comm.is_null());
    //           comm_ = *mpi_comm->getRawMpiComm();
    // #endif
    //         }
    //       }

    //       // Resize the buffer to accommodate nvec vectors in the multivector, for a
    //       // matrix having block size block_size.
    //       void resize (const int& block_size, const int& nvec) {
    //         bsz_ = block_size;
    //         nvec_ = nvec;
    //         wrk_.resize((2*block_size + 1)*nvec);
    //       }
      
    //       // Check the norm every sweep_step sweeps.
    //       void setCheckFrequency (const int& sweep_step) {
    //         TEUCHOS_TEST_FOR_EXCEPT_MSG(sweep_step < 1, "sweep step must be >= 1");
    //         sweep_step_ = sweep_step;
    //       }
      
    //       // Get the buffer into which to store rank-local squared norms.
    //       magnitude_type* get_buffer () { return wrk_.data(); }

    //       // Call MPI_Iallreduce to find the global squared norms.
    //       void ireduce (const int& sweep, const bool force = false) {
    //         if ( ! force && sweep % sweep_step_) return;
    //         const int n = bsz_*nvec_;
    //         if (collective_) {
    //           std::copy(wrk_.begin(), wrk_.begin() + n, wrk_.begin() + n);
    // #ifdef HAVE_MPI
    // #if MPI_VERSION >= 3
    //           MPI_Iallreduce(wrk_.data() + n, wrk_.data(), n,
    //                          Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
    //                          MPI_SUM, comm_, &mpi_request_);
    // #else
    //           MPI_Allreduce (wrk_.data() + n, wrk_.data(), n,
    //                          Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
    //                          MPI_SUM, comm_);
    // #endif
    // #endif
    //         }
    //       }
      
    //       // Check if the norm-based termination criterion is met. tol2 is the
    //       // tolerance squared. Sweep is the sweep index. If not every iteration is
    //       // being checked, this function immediately returns false. If a check must
    //       // be done at this iteration, it waits for the reduction triggered by
    //       // ireduce to complete, then checks the global norm against the tolerance.
    //       bool checkDone (const int& sweep, const magnitude_type& tol2, const bool force = false) {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //         TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NormManager::check_done");
    // #endif
    //         TEUCHOS_ASSERT(sweep >= 1);
    //         if ( ! force && (sweep - 1) % sweep_step_) return false;
    //         if (collective_) {
    // #ifdef HAVE_MPI
    // # if MPI_VERSION >= 3
    //           MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
    // # else
    //           // Do nothing.
    // # endif
    // #endif
    //         }
    //         const auto n = bsz_*nvec_;
    //         if (sweep == 1) {
    //           magnitude_type* const n0 = wrk_.data() + 2*n;
    //           for (int v = 0; v < nvec_; ++v) {
    //             const magnitude_type* const dn0 = wrk_.data() + v*bsz_;
    //             magnitude_type mdn0 = 0;
    //             for (int i = 0; i < bsz_; ++i)
    //               mdn0 = std::max(mdn0, dn0[i]);
    //             n0[v] = mdn0;
    //           }
    //           return false;
    //         } else {
    //           const auto n0 = wrk_.data() + 2*n;
    //           bool done = true;
    //           for (int v = 0; v < nvec_; ++v) {
    //             const magnitude_type* const dnf = wrk_.data() + v*bsz_;
    //             magnitude_type mdnf = 0;
    //             for (int i = 0; i < bsz_; ++i)
    //               mdnf = std::max(mdnf, dnf[i]);
    //             if (mdnf > tol2*n0[v]) {
    //               done = false;
    //               break;
    //             }
    //           }
    //           return done;
    //         }
    //       }
      
    //       // After termination has occurred, finalize the norms for use in
    //       // get_norms{0,final}.
    //       void finalize () {
    //         for (int v = 0; v < nvec_; ++v) {
    //           const magnitude_type* const dnf = wrk_.data() + v*bsz_;
    //           magnitude_type mdnf = 0;
    //           for (int i = 0; i < bsz_; ++i)
    //             mdnf = std::max(mdnf, dnf[i]);
    //           // This overwrites the receive buffer, but that's OK; at the time of
    //           // this write, we no longer need the data in this slot.
    //           wrk_[v] = mdnf;
    //         }
    //         for (int i = 0; i < nvec_; ++i)
    //           wrk_[i] = std::sqrt(wrk_[i]);
    //         magnitude_type* const nf = wrk_.data() + 2*bsz_*nvec_;
    //         for (int v = 0; v < nvec_; ++v)
    //           nf[v] = std::sqrt(nf[v]);
    //       }
      
    //       // Report norms to the caller.
    //       const magnitude_type* getNorms0 () const { return wrk_.data() + 2*bsz_*nvec_; }
    //       const magnitude_type* getNormsFinal () const { return wrk_.data(); }
    //     };

    
    //     // Repack from the BlockCrsMatrix to the tridiagonal storage, and factorize
    //     // the data in the latter.
    //     template<typename MatrixType> 

    //       // Solve D X = Y, where A + D + R. This general implementation uses subviews
    //       // and works on any platform.
    //       template <typename Layout, typename PartialSpecialization>
    //       class Solve {
    //       public: // for Cuda
    //         // Whether RHS is a vector or a multivector.
    //         struct VecTag {};
    //         struct MatTag {};

    //       private:
    //         ConstUnmanaged<SizeList> pack_td_ptr;
    //         ConstUnmanaged<typename BlockTridiags::Values> values;
    //         ConstUnmanaged<LOList> packptr, part2packrowidx0;
    //         Unmanaged<typename BlockTridiags::PackedMultiVector> X;

    //         template <typename Tag> void run () {
    //           Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
    //           Kokkos::parallel_for(rp, *this);
    //         }

    //         // X := a T \ X
    //         template <typename Tag, typename UploType, typename DiagType,
    //                   typename Scalar, typename MatA, typename MatX>
    //         KOKKOS_FORCEINLINE_FUNCTION static void
    //         trsmv (const Scalar& a, const MatA& A, MatX& X,
    //                typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsv<UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsv::Unblocked>
    //             ::invoke(a, A, X);
    //         }

    //         template <typename Tag, typename UploType, typename DiagType,
    //                   typename Scalar, typename MatA, typename MatX>
    //         KOKKOS_FORCEINLINE_FUNCTION static void
    //         trsmv (const Scalar& a, const MatA& A, MatX& X,
    //                typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsm<kbe::Side::Left, UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsm::Blocked>
    //             ::invoke(a, A, X);
    //         }

    //         // C := b C + a A B
    //         template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
    //         KOKKOS_FORCEINLINE_FUNCTION static void
    //         gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
    //                typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialGemv<kbe::Trans::NoTranspose, kbe::Algo::Gemv::Unblocked>::invoke(a, A, B, b, C);
    //         }

    //         template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
    //         KOKKOS_FORCEINLINE_FUNCTION static void
    //         gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
    //                typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialGemm<kbe::Trans::NoTranspose, kbe::Trans::NoTranspose, kbe::Algo::Gemm::Blocked>
    //             ::invoke(a, A, B, b, C);
    //         }

    //       public:
    //         Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
    //                const typename BlockTridiags::PackedMultiVector& X_)
    //           : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
    //             packptr(packptr_), part2packrowidx0(part2packrowidx0_),
    //             X(X_)
    //         {}

    //         template <typename Tag>
    //         KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
    //           using Kokkos::subview;
    //           using Kokkos::ALL;
    //           namespace kbe = KokkosBatched::Experimental;

    //           const auto one = Kokkos::ArithTraits<magnitude_type>::one();

    //           const auto partidx = packptr(packidx);
    //           Size i0 = pack_td_ptr(partidx);
    //           LO r0 = part2packrowidx0(partidx);
    //           const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

    //           // Solve L x = x.
    //           auto A = subview(values, i0, ALL(), ALL());
    //           auto X1 = subview(X, r0, ALL(), ALL());
    //           trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X1);
    //           for (LO i = 1; i < nrows; ++i) {
    //             const auto B = subview(values, i0+2, ALL(), ALL());
    //             r0++;
    //             const auto X2 = subview(X, r0, ALL(), ALL());
    //             gemmv<Tag>(-one, B, X1, one, X2);
    //             i0 += 3;
    //             A = subview(values, i0, ALL(), ALL());
    //             trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X2);
    //             X1 = X2;
    //           }

    //           // Solve U x = x.
    //           trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X1);
    //           for (LO i = nrows; i > 1; --i) {
    //             i0 -= 3;
    //             const auto B = subview(values, i0+1, ALL(), ALL());
    //             r0--;
    //             const auto X2 = subview(X, r0, ALL(), ALL());
    //             gemmv<Tag>(-one, B, X1, one, X2);
    //             A = subview(values, i0, ALL(), ALL());
    //             trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X2);
    //             X1 = X2;
    //           }
    //         }

    //         void run () {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //           TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
    // #endif
    //           if (X.extent(2) == 1)
    //             run<VecTag>();
    //           else
    //             run<MatTag>();
    //         }
    //       };

    //       // This specialization is for CPU and KNL. On KNL, it speeds up the solve by
    //       // ~20% for block size 5, more for smaller, less for bigger. The importance of
    //       // the apply phase of a preconditioner justifies the extra complexity of this
    //       // specialization.
    //       template <typename PartialSpecialization>
    //       class Solve<Kokkos::LayoutRight, PartialSpecialization> {
    //       public: // for Cuda
    //         // Whether RHS is a vector or a multivector.
    //         struct VecTag {};
    //         struct MatTag {};

    //       private:
    //         ConstUnmanaged<SizeList> pack_td_ptr;
    //         ConstUnmanaged<typename BlockTridiags::Values> values;
    //         ConstUnmanaged<LOList> packptr, part2packrowidx0;
    //         Unmanaged<typename BlockTridiags::PackedMultiVector> X;
    //         const int bs2, xsz;

    //         template <typename Tag> void run () {
    //           Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
    //           Kokkos::parallel_for(rp, *this);
    //         }

    //         // For speed in this important, O(nnz) operation, we don't use subviews, but
    //         // instead go directly to the raw-array interface.

    //         // X := a T \ X

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
    //                  typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsvInternalLower<kbe::Algo::Trsv::Unblocked>::invoke(
    //                                                                            true,
    //                                                                            values.dimension_1(),
    //                                                                            a,
    //                                                                            values.data() + i0*bs2, values.stride_1(), values.stride_2(),
    //                                                                            X.data() + r0*xsz, X.stride_1());
    //         }

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
    //                  typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsmInternalLeftLower<kbe::Algo::Trsm::Blocked>::invoke(
    //                                                                              true,
    //                                                                              values.dimension_1(), X.dimension_2(),
    //                                                                              a,
    //                                                                              values.data() + i0*bs2, values.stride_1(), values.stride_2(),
    //                                                                              X.data() + r0*xsz, X.stride_1(), X.stride_2());
    //         }

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
    //                  typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsvInternalUpper<kbe::Algo::Trsv::Unblocked>::invoke(
    //                                                                            false,
    //                                                                            values.dimension_1(),
    //                                                                            a,
    //                                                                            values.data() + i0*bs2, values.stride_1(), values.stride_2(),
    //                                                                            X.data() + r0*xsz, X.stride_1());
    //         }

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
    //                  typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialTrsmInternalLeftUpper<kbe::Algo::Trsm::Blocked>::invoke(
    //                                                                              false,
    //                                                                              values.dimension_1(), X.dimension_2(),
    //                                                                              a,
    //                                                                              values.data() + i0*bs2, values.stride_1(), values.stride_2(),
    //                                                                              X.data() + r0*xsz, X.stride_1(), X.stride_2());
    //         }

    //         // C := b C + a A B

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
    //                typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialGemvInternal<kbe::Algo::Gemv::Unblocked>::invoke(
    //                                                                       values.dimension_1(), values.dimension_2(),
    //                                                                       a,
    //                                                                       values.data() + a0*bs2, values.stride_1(), values.stride_2(),
    //                                                                       X.data() + b0*xsz, X.stride_1(),
    //                                                                       b,
    //                                                                       X.data() + c0*xsz, X.stride_1());
    //         }

    //         template <typename Tag>
    //         KOKKOS_FORCEINLINE_FUNCTION void
    //         gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
    //                typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
    //           namespace kbe = KokkosBatched::Experimental;
    //           kbe::SerialGemmInternal<kbe::Algo::Gemm::Blocked>::invoke(
    //                                                                     X.dimension_1(), X.dimension_2(), values.dimension_2(),
    //                                                                     a,
    //                                                                     values.data() + a0*bs2, values.stride_1(), values.stride_2(),
    //                                                                     X.data() + b0*xsz, X.stride_1(), X.stride_2(),
    //                                                                     b,
    //                                                                     X.data() + c0*xsz, X.stride_1(), X.stride_2());
    //         }

    //       public:
    //         Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
    //                const typename BlockTridiags::PackedMultiVector& X_)
    //           : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
    //             packptr(packptr_), part2packrowidx0(part2packrowidx0_),
    //             X(X_), bs2(values.dimension_1()*values.dimension_1()),
    //             xsz(values.dimension_1()*X.dimension_2())
    //         {}

    //         template <typename Tag>
    //         KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
    //           using Kokkos::subview;
    //           using Kokkos::ALL;
    //           namespace kbe = KokkosBatched::Experimental;

    //           const auto one = Kokkos::ArithTraits<magnitude_type>::one();

    //           const auto partidx = packptr(packidx);
    //           Size i0 = pack_td_ptr(partidx);
    //           LO r0 = part2packrowidx0(partidx);
    //           const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

    //           // Solve L x = x.
    //           trsmvlo<Tag>(one, i0, r0);
    //           for (LO i = 1; i < nrows; ++i) {
    //             r0++;
    //             gemmv<Tag>(-one, one, i0+2, r0-1, r0);
    //             i0 += 3;
    //             trsmvlo<Tag>(one, i0, r0);
    //           }

    //           // Solve U x = x.
    //           trsmvup<Tag>(one, i0, r0);
    //           for (LO i = nrows; i > 1; --i) {
    //             i0 -= 3;
    //             r0--;
    //             gemmv<Tag>(-one, one, i0+1, r0+1, r0);
    //             trsmvup<Tag>(one, i0, r0);
    //           }
    //         }

    //         void run () {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //           TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
    // #endif
    //           if (X.extent(2) == 1)
    //             run<VecTag>();
    //           else
    //             run<MatTag>();
    //         }
    //       };

    //     private:

    //       template <typename T> static KOKKOS_INLINE_FUNCTION constexpr T square (const T& x) { return x*x; }


    //       // Top-level Impl object initialization.
    //       void init (const Teuchos::RCP<const row_matrix_type>& matrix,
    //                  const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
    //                  const Teuchos::RCP<const import_type>& importer, int overlapLevel,
    //                  bool overlapCommAndComp, bool useSeqMethod) {
    //         seq_method_ = useSeqMethod;
    //         overlap_comm_ = overlapCommAndComp;
    //         validate_ = true;

    //         TEUCHOS_TEST_FOR_EXCEPT_MSG(
    //                                     overlapLevel != 0,
    //                                     "BlockTriDiContainer does not curently support OverlapLevel != 0; user provided "
    //                                     << overlapLevel);

    //         A_bcrs_ = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(matrix);
    //         TEUCHOS_TEST_FOR_EXCEPT_MSG(
    //                                     A_bcrs_.is_null(), "BlockTriDiContainer currently supports Tpetra::BlockCrsMatrix only.");

    //         init_importer(importer);
    //         init_parts(partitions);
    //         init_btdm();
    //       }

    //       }




    //       // Wrappers to the compute_b_minus_Rx implementations.
    //       // Y := B - R X, where A = D + R.
    //       void compute_b_minus_Rx (const mv_type& B, const mv_type& X, mv_type& Y) const {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //         TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
    // #endif
    //         const auto Bv = B.template getLocalView<device_type>();
    //         const auto Xv = X.template getLocalView<device_type>();
    //         const auto Yv = Y.template getLocalView<device_type>();
    //         if (amd_.is_tpetra_block_crs) {
    //           const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
    //           BlockTriDiContainerDetails::compute_b_minus_Rx(Bv, Xv, Yv, A_bcrs_->getBlockSize(),
    //                                                          amd_.rowptr, amd_.A_colindsub,
    //                                                          g.row_map, g.entries, amd_.tpetra_values);
    //         }
    //       }

    //       template <typename MvView>
    //       void compute_b_minus_Rx (const mv_type& B, const MvView& Xv,
    //                                typename BlockTridiags::PackedMultiVector& Y,
    //                                const bool first_apply) const {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //         TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
    // #endif
    //         const auto Bv = B.template getLocalView<device_type>();
    //         if (amd_.is_tpetra_block_crs) {
    //           const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
    //           BlockTriDiContainerDetails::compute_b_minus_Rx(
    //                                                          first_apply ? amd_.rowptr : amd_.rowptr_remote,
    //                                                          first_apply ? amd_.A_colindsub : amd_.A_colindsub_remote,
    //                                                          g.row_map, g.entries, amd_.tpetra_values,
    //                                                          Bv, Xv, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_, first_apply);
    //         }
    //       }

    //       template <typename MvView>
    //       void compute_b_minus_Rx (const mv_type& B, const MvView& X_owned, const MvView& X_remote,
    //                                typename BlockTridiags::PackedMultiVector& Y) const {
    // #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    //         TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
    // #endif
    //         const auto Bv = B.template getLocalView<device_type>();
    //         if (amd_.is_tpetra_block_crs) {
    //           const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
    //           BlockTriDiContainerDetails::compute_b_minus_Rx(
    //                                                          amd_.rowptr, amd_.A_colindsub,
    //                                                          g.row_map, g.entries, amd_.tpetra_values,
    //                                                          Bv, X_owned, X_remote, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_);
    //         }
    //       }

    //     public:
    //       Impl (BlockTriDiContainer<MatrixType, LocalScalarType>& container,
    //             const Teuchos::RCP<const row_matrix_type>& matrix,
    //             const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
    //             const Teuchos::RCP<const import_type>& importer, int overlapLevel,
    //             bool overlapCommAndComp = false, bool useSeqMethod = false)
    //         : container_(container)
    //       {
    //         init(matrix, partitions, importer, overlapLevel, overlapCommAndComp, useSeqMethod);
    //       }

    //       std::string describe () const {
    //         std::stringstream ss;
    //         ss << "seq_method " << seq_method_
    //            << " overlap_comm " << overlap_comm_
    //            << " dm2cm " << (dm2cm_.data() ? true : false);
    //         return ss.str();
    //       }


    //       static void debug_print (const Tridiags& t) {
    //         const auto& v = t.values;
    //         std::stringstream ss;
    //         ss << "v = [";
    //         for (size_t pi = 0; pi < t.pack_td_ptr.size() - 1; ++pi)
    //           for (LO i1 = 0; i1 < vectorization_traits::vector_length; ++i1)
    //             for (Size ind = t.pack_td_ptr(pi); ind < t.pack_td_ptr(pi+1); ++ind) {
    //               const auto i = Kokkos::make_pair(ind,i1);
    //               for (LO j = 0; j < v.extent_int(1); ++j)
    //                 for (LO k = 0; k < v.extent_int(2); ++k)
    //                   ss << " " << Details::Batched::idx(v,i,j,k);
    //             }
    //         ss << "]\n";
    //         std::cout << ss.str();
    //       }


    //     // Base class for any unimplemented typename combinations.
    //     template <typename MatrixType, typename LocalScalarType, typename ExeSpace>
    //     struct UnImpl {
    //       typedef typename MatrixType::scalar_type scalar_type;
    //       typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
    //       typedef typename MatrixType::local_ordinal_type local_ordinal_type;
    //       typedef typename MatrixType::global_ordinal_type global_ordinal_type;
    //       typedef typename MatrixType::node_type node_type;
    //       typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
    //       typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
    //       typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

    //       std::string describe () const { return ""; }
    //       int get_block_size () const { return 0; }
    //       int get_nvec () const { return 0; }
    //       const magnitude_type* get_norms0 () const { return nullptr; }
    //       const magnitude_type* get_norms_final () const { return nullptr; }
    //       void symbolic () {}
    //       void numeric (const magnitude_type add_to_diag = 0) {}
    //       int applyInverseJacobi (const mv_type& X, mv_type& Y, const scalar_type damping_factor,
    //                               bool y_is_zero, const int max_num_sweeps, magnitude_type tolerance,
    //                               const int check_tol_every) const
    //       { return 0; }
    //     };

    // #if defined HAVE_STOKHOS_IFPACK2
    //     // BlockTriDiContainer is built on KokkosBatch's SIMD type, which doesn't work
    //     // directly with Stokhos composite types. Hence, catch these with UnImpl.
    // #define IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Type)                \
    //     template <typename MatrixType, typename T, typename ExeSpace>       \
    //     struct ImplType<MatrixType, Type<T>, ExeSpace> : UnImplType<MatrixType, Type<T>, ExeSpace> { \
    //       typedef UnImplType<MatrixType, Type<T>, ExeSpace> ui;                 \
    //       Impl (BlockTriDiContainer<MatrixType, typename ui::scalar_type>& container, \
    //             const Teuchos::RCP<const typename ui::row_matrix_type>& matrix, \
    //             const Teuchos::Array<Teuchos::Array<typename ui::local_ordinal_type> >& partitions, \
    //             const Teuchos::RCP<const typename ui::import_type>& importer, int overlapLevel, \
    //             bool overlapCommAndComp = false, bool useSeqMethod = false) { \
    //         TEUCHOS_TEST_FOR_EXCEPT_MSG(                                    \
    //                                     true, "BlockTriDiContainer is not currently supported for Stokhos composite types."); \
    //       }                                                                 \
    //     };
    //     IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::MP::Vector)
    //     IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::UQ::PCE)
    // #undef IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS
    // #endif
  
} // namespace BlockTriDiContainerDetails

} // namespace Ifpack2

#endif
