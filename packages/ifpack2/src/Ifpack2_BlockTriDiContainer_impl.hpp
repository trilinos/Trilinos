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

#include "Ifpack2_BlockTriDiContainer_decl.hpp"
#include "Ifpack2_BlockTriDiContainer_impl.hpp"

#include <memory>

namespace Ifpack2 {

  namespace BlockTriDiContainerDetails {

    /// 
    /// implementation typedefs
    ///
    template <typename MatrixType>
    struct Impl {
      ///
      /// matrix type derived types
      ///
      typedef size_t size_type;
      typedef typename MatrixType::scalar_type scalar_type;
      typedef typename MatrixType::local_ordinal_type local_ordinal_type;
      typedef typename MatrixType::global_ordinal_type global_ordinal_type;
      typedef typename MatrixType::node_type node_type;

      ///
      /// kokkos arithmetic traits
      ///
      typedef typename Kokkos::Details::ArithTraits<scalar_type>::val_type impl_scalar_type;
      typedef typename Kokkos::ArithTraits<impl_scalar_type>::mag_type magnitude_type;

      ///
      /// tpetra types
      ///
      typedef typename node_type::device_type tpetra_device_type;      
      typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_multivector_type;
      typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> tpetra_map_type;
      typedef Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> tpetra_import_type;
      typedef Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_row_matrix_type;
      typedef Tpetra::BlockCrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_crs_matrix_type;
      typedef Tpetra::BlockMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_multivector_type;
      typedef typename tpetra_block_crs_matrix_type::crs_graph_type::local_graph_type local_crs_graph_type;
      //typedef Kokkos::View<impl_scalar_type*,tpetra_device_type> local_blockcrs_value_array_type;
      typedef typename Kokkos::View<double***,tpetra_device_type>::array_layout tpetra_layout_type;
      
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
      
      ///
      /// simd vectorization
      ///
      using KokkosBatched::Experimental::Vector;
      using KokkosBatched::Experimental::SIMD;
      using KokkosBatched::Experimental::DefaultVectorLength;

      enum : int { vector_length = DefaultVectorLength<impl_scalar_type,tpetra_device_type::memory_space>::value };
      typedef typename Vector<SIMD<impl_scalar_type>,vector_length> vector_type;

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
      using Const = Kokkos::View<typename ViewType::const_data_type, typename ViewType::array_layout,
                                 typename ViewType::device_type, typename ViewType::memory_traits>;
      template <typename ViewType>
      using ConstUnmanaged = Const<Unmanaged<ViewType> >;

      ///
      /// view types that are often used
      ///
      typedef Kokkos::View<size_type*,device_type> size_type_1d_view;
      typedef Kokkos::View<local_ordinal_type*,device_type> local_ordinal_type_1d_view;
      // tpetra block crs values
      typedef Kokkos::View<impl_scalar_type*,device_type> impl_scalar_type_1d_view;
      // packed data always use layout right
      typedef Kokkos::View<vector_type***,Kokkos::LayoutRight,device_type> vector_type_3d_view;
      typedef Kokkos::View<impl_scalar_type****,Kokkos::LayoutRight,Kokkos::device_type> impl_scalar_type_4d_view;

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
      std::string get_msg_prefix (const CommPtrType &comm) const {
        const auto rank = comm->getRank();
        const auto nranks = comm->getSize();
        std::stringstream ss;
        ss << "Rank " << rank << " of " << nranks << ": ";
        return ss.str();
      }

    };

    template<typename MatrixType>
    struct PartInterface {
      using Impl<MatrixType>::tpetra_device_type device_type;
      using Impl<MatrixType>::local_ordinal_type local_ordinal_type;

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
      Kokkos::View<local_ordinal_type*,device_type> lclrow;
      // partptr_ is the pointer array into lclrow_.
      Kokkos::View<local_ordinal_type*,device_type> partptr; // np+1
      // packptr_(i), for i the pack index, indexes partptr_. partptr_(packptr_(i))
      // is the start of the i'th pack.
      Kokkos::View<local_ordinal_type*,device_type> packptr; // npack+1
      // part2rowidx0_(i) is the flat row index of the start of the i'th part. It's
      // an alias of partptr_ in the case of no overlap.
      Kokkos::View<local_ordinal_type*,device_type> part2rowidx0; // np+1
      // part2packrowidx0_(i) is the packed row index. If vector_length is 1, then
      // it's the same as part2rowidx0_; if it's > 1, then the value is combined
      // with i % vector_length to get the location in the packed data.
      Kokkos::View<local_ordinal_type*,device_type> part2packrowidx0; // np+1
      local_ordinal_type part2packrowidx0_back; // So we don't need to grab the array from the GPU.
      // rowidx2part_ maps the row index to the part index.
      Kokkos::View<local_ordinal_type*,device_type> rowidx2part; // nr
      // If needed, permute owned domain main LIDs to owned column map LIDs.
      Kokkos::View<local_ordinal_type*,device_type> dm2cm;
      // True if lcl{row|col} is at most a constant away from row{idx|col}. In
      // practice, this knowledge is not particularly useful, as packing for batched
      // processing is done at the same time as the permutation from LID to index
      // space. But it's easy to detect, so it's recorded in case an optimization
      // can be made based on it.
      bool row_contiguous, col_contiguous;
    };

    ///
    /// block tridiagonals
    ///
    template <typename MatrixType>
    struct BlockTriDiags {
      using Impl<MatrixType>::tpetra_device_type device_type;
      using Impl<MatrixType>::local_ordinal_type local_ordinal_type;
      using Impl<MatrixType>::size_type size_type;
      using Impl<MatrixType>::vector_type vector_type;

      typedef Kokkos::View<size_type,device_type> size_type_view;
      typedef Kokkos::View<local_ordinal_type,device_type> local_ordinal_type_view;
      //typedef Kokkos::View<vector_type,device_type> values_type_view;
      
      // flat_td_ptr(i) is the index into flat-array values of the start of the
      // i'th tridiag. pack_td_ptr is the same, but for packs. If vector_length ==
      // 1, pack_td_ptr is the same as flat_td_ptr; if vector_length > 1, then i %
      // vector_length is the position in the pack.
      Kokkos::View<size_type,device_type> flat_td_ptr, pack_td_ptr;
      // List of local column indices into A from which to grab
      // data. flat_td_ptr(i) points to the start of the i'th tridiag's data.
      Kokkos::View<local_ordinal_type,device_type> A_colindsub;
      // Tridiag block values. pack_td_ptr(i) points to the start of the i'th
      // tridiag's pack, and i % vector_length gives the position in the pack.
      Kokkos::View<vector_type***,device_type> values;

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
    /// A - Tridiags(A), i.e., R in the splitting A = D + R.
    ///
    template <typename MatrixType>
    struct AmD {
      using Impl<MatrixType>::tpetra_device_type device_type;
      using Impl<MatrixType>::local_ordinal_type local_ordinal_type;
      using Impl<MatrixType>::size_type size_type;
      using Impl<MatrixType>::impl_scalar_type impl_scalar_type;

      // rowptr points to the start of each row of A_colindsub.
      Kokkos::View<size_type,device_type> rowptr, rowptr_remote;
      // Indices into A's rows giving the blocks to extract. rowptr(i) points to
      // the i'th row. Thus, g.entries(A_colindsub(rowptr(row) : rowptr(row+1))),
      // where g is A's graph, are the columns AmD uses. If seq_method_, then
      // A_colindsub contains all the LIDs and A_colindsub_remote is empty. If !
      // seq_method_, then A_colindsub contains owned LIDs and A_colindsub_remote
      // contains the remote ones.
      Kokkos::View<local_ordinal_type,device_type> A_colindsub, A_colindsub_remote;
      
      // Currently always true.
      bool is_tpetra_block_crs;

      // If is_tpetra_block_crs, then this is a pointer to A_'s value data.
      Kokkos::View<impl_scalar_type*,device_type> tpetra_values;
    };

    ///
    /// setup sequential importer
    ///
    template<typename MatrixType>
    static typename Impl<MatrixType>::tpetra_import_type
    createBlockCrsSequentialImporter(const typename Impl<MatrixType>::tpetra_block_crs_matrix_type &A) {
      using Impl<MatrixType>::tpetra_map_type;
      using Impl<MatrixType>::tpetra_block_multivector_type::makePointMap;
      
      const auto blocksize = A.getBlockSize(); // local ordinal
      const auto crs_graph = A.getCrsGraph(); // tpetra crs graph object

      const auto domain_map = crs_graph.getDomainMap(); //rcp
      const auto column_map = crs_graph.getColMap(); // rcp

      // create point map
      return Impl<MatrixType>::tpetra_import_type(Teuchos::rcp(new tpetra_map_type(makePointMap(domain_map))),
                                                  Teuchos::rcp(new tpetra_map_type(makePointMap(column_map))));
    }

//     ///
//     /// setup part interface using the container partitions array
//     ///
//     template<typename MatrixType>
//     static PartInterface<MatrixType> 
//     void createPartInterface(const typename Impl<MatrixType>::tpetra_block_crs_matrix_type &A,
//                              const Teuchos::Array<Teuchos::Array<typename Impl<MatrixType>::local_ordinal_type> > &partitions) {
//       using Impl<MatrixType>::tpetra_device_type device_type;
//       using Impl<MatrixType>::local_ordinal_type local_ordinal_type;

//       const auto comm = A.getRowMap()->getCopmm();
      
//       //PartInterface<MatrixType> interf;
      
//       const bool jacobi = partitions.extent(0) == 0;
//       const local_ordinal_type A_n_lclrows = A.getNodeNumRows();
//       const local_ordinal_type nparts = jacobi ? A_n_lclrows : partitions.extent(0);

// #if defined(BLOCKTRIDICONTAINER_DEBUG)              
//       local_ordinal_type nrows = 0;
//       if (jacobi)       nrows = nparts;
//       else              for (local_ordinal_type i = 0; i < nparts; ++i) nrows += partitions[i].size();

//       TEUCHOS_TEST_FOR_EXCEPT_MSG
//         (nrows != A_n_lclrows, get_msg_prefix(comm) << "The #rows implied by the local partition is not "
//          << "the same as getNodeNumRows: " << nrows << " vs " << A_n_lclrows);
// #endif

//       std::vector<local_ordinal_type> p;
//       if ( !jacobi) {
//         const LO nparts = partitions.size();
//         p.resize(nparts);
        
//         if (vectorization_traits::vector_length == 1) {
//           for (LO i = 0; i < nparts; ++i)
//             p[i] = i;
//         } else {
//           typedef std::pair<LO,LO> SzIdx;
//           std::vector<SzIdx> partsz(nparts);
//           for (LO i = 0; i < nparts; ++i)
//             partsz[i] = SzIdx(partitions[i].size(), i);
//           std::sort(partsz.begin(), partsz.end(),
//                     [=] (const SzIdx& x, const SzIdx& y) { return x.first > y.first; });
//           for (LO i = 0; i < nparts; ++i)
//             p[i] = partsz[i].second;
//         }
//       }
      
      
//     }
//     partptr_ = local_ordinal_typeList("partptr_", nparts + 1);
//         lclrow_ = local_ordinal_typeList("lclrow_", A_n_lclrows);
//         part2rowidx0_ = local_ordinal_typeList("part2rowidx0_", nparts + 1);
//         part2packrowidx0_ = local_ordinal_typeList("part2packrowidx0_", nparts + 1);
//         rowidx2part_ = local_ordinal_typeList("rowidx2part_", A_n_lclrows);
//         const auto partptr = Kokkos::create_mirror_view(partptr_);
//         const auto lclrow = Kokkos::create_mirror_view(lclrow_);
//         const auto part2rowidx0 = Kokkos::create_mirror_view(part2rowidx0_);
//         const auto part2packrowidx0 = Kokkos::create_mirror_view(part2packrowidx0_);
//         const auto rowidx2part = Kokkos::create_mirror_view(rowidx2part_);

//         // Determine parts.
//         row_contiguous_ = true;
//         partptr(0) = 0;
//         part2rowidx0(0) = 0;
//         part2packrowidx0(0) = 0;
//         LO pack_nrows = 0;
//         for (LO ip = 0; ip < nparts; ++ip) {
//           const auto* part = jacobi ? nullptr : &partitions[pp[ip]];
//           const LO ipnrows = jacobi ? 1 : part->size();
//           TEUCHOS_ASSERT(ip == 0 || vectorization_traits::vector_length == 1 ||
//                          (jacobi || ipnrows <= static_cast<LO>(partitions[pp[ip-1]].size())));
//           TEUCHOS_TEST_FOR_EXCEPT_MSG(
//                                       ipnrows == 0, get_msg_prefix() << "partition " << pp[ip]
//                                       << " is empty, which is not allowed.");
//           //assume No overlap.
//           part2rowidx0(ip+1) = part2rowidx0(ip) + ipnrows;
//           // Since parts are ordered in nonincreasing size, the size of the first
//           // part in a pack is the size for all parts in the pack.
//           if (ip % vectorization_traits::vector_length == 0)
//             pack_nrows = ipnrows;
//           part2packrowidx0(ip+1) = part2packrowidx0(ip) +
//             ((ip+1) % vectorization_traits::vector_length == 0 || ip+1 == nparts ? pack_nrows : 0);
//           const LO os = partptr(ip);
//           for (LO i = 0; i < ipnrows; ++i) {
//             const auto lcl_row = jacobi ? ip : (*part)[i];
//             TEUCHOS_TEST_FOR_EXCEPT_MSG(
//                                         lcl_row < 0 || lcl_row >= A_n_lclrows, get_msg_prefix() << "partitions[" << pp[ip] << "]["
//                                         << i << "] = " << lcl_row << " but input matrix implies limits of [0, " << A_n_lclrows-1
//                                         << "].");
//             lclrow(os+i) = lcl_row;
//             rowidx2part(os+i) = ip;
//             if (row_contiguous_ && os+i > 0 && lclrow((os+i)-1) + 1 != lcl_row)
//               row_contiguous_ = false;
//           }
//           partptr(ip+1) = os + ipnrows;
//         }
//         TEUCHOS_ASSERT(partptr(nparts) == nrows);
//         if (lclrow(0) != 0) row_contiguous_ = false;

//         Kokkos::deep_copy(partptr_, partptr);
//         Kokkos::deep_copy(lclrow_, lclrow);
//         //assume No overlap. Thus:
//         part2rowidx0_ = partptr_;
//         if (vectorization_traits::vector_length > 1)
//           Kokkos::deep_copy(part2packrowidx0_, part2packrowidx0);
//         else
//           part2packrowidx0_ = part2rowidx0_;
//         part2packrowidx0_back = part2packrowidx0(part2packrowidx0.size() - 1);
//         Kokkos::deep_copy(rowidx2part_, rowidx2part);

//         { // Fill packptr.
//           LO npacks = 0;
//           for (LO ip = 1; ip <= nparts; ++ip)
//             if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
//               ++npacks;
//           packptr_ = LOList("packptr_", npacks + 1);
//           const auto packptr = Kokkos::create_mirror_view(packptr_);
//           packptr(0) = 0;
//           for (LO ip = 1, k = 1; ip <= nparts; ++ip)
//             if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
//               packptr(k++) = ip;
//           Kokkos::deep_copy(packptr_, packptr);
//         }
//       }
    
    ///
    /// block tridiags initialization from part interface
    ///
    template<typename MatrixType>
    static 
    BlockTridiags<MatrixType> 
    createBlockTridiags(const PartInterface<MatrixType> &interf) {
      using Impl<MatrixType>::local_ordinal_type;
      using Impl<MatrixType>::size_type;
      using Impl<MatrixType>::vector_length;

      BlockTridiags<MatrixType> btdm;

      const auto partptr = create_host_mirror_view_and_sync(interf.partptr);
      const local_ordinal_type ntridiags = partptr.extent(0) - 1;
      
      { // Construct the flat index pointers into the tridiag values array.
        btdm.flat_td_ptr = Kokkos::View<size_type,device_type>("btdm.flat_td_ptr", ntridiags + 1);
        Kokkos::parallel_scan
          (Kokkos::RangePolicy<device_type>(0,ntridiags+1), 
           KOKKOS_LAMBDA(const local_ordinal_type i,
                         size_type &update,
                         const bool &final) {
            const local_ordinal_type cnt = i < ntridiags ? (3*(interf.partptr(i+1) - interf.partptr(i)) - 2) : 0;
            btdm.flat_td_ptr(i) = cnt;
            if (final) 
              btdm.flat_td_ptr(i) = update;
            update += cnt;
          });

#if defined(BLOCKTRIDICONTAINER_DEBUG)        
        const auto last_value_in_parptr = create_host_mirror_view_and_sync(Kokkos::subview(interf.partptr, ntridiags));        
        const size_type nnz = 3*last_value_in_parptr() - 2*ntridiags;
        const auto last_value_in_flat_td_ptr = create_host_mirror_view_and_sync(Kokkos::subview(btdm.flat_td_ptr, ntridiags));        
        TEUCHOS_ASSERT(last_value_in_flat_td_ptr() == nnz);
#endif
      }
      
      // And the packed index pointers.
      if (vector_length == 1) {
        btdm_.pack_td_ptr = btdm_.flat_td_ptr;
      } else {
        const local_ordinal_type npacks = interf.packptr.extent(0) - 1;
        Kokkos::View<size_type,device_type> work("work", npacks+1);        
        btdm.pack_td_ptr = Kokkos::View<size_type,device_type>("btdm.pack_td_ptr", ntridiags + 1);
        Kokkos::parallel_scan
          (Kokkos::RangePolicy<device_type>(0,npacks+1), 
           KOKKOS_LAMBDA(const local_ordinal_type i,
                         size_type &update,
                         const bool &final) {
            const auto ii = packptr(i);
            const local_ordinal_type cnt = i < npacks ? BlockTridiags::NumBlocks(interf.partptr(ii+1) - partptr(ii)) : 0;
            work(i) = cnt;
            if (final) 
              work(i) = update;
            update += cnt;
          });
        
        Kokkos::parallel_for
          (Kokkos::RangePolicy<device_type>(0,npacks), 
           KOKKOS_LAMBDA(const local_ordinal_type i) {
            const auto ii = packptr(i);
            for (local_ordinal_type ii = packptr(i); ii < packptr(i+1); ++ii)            
              pack_td_ptr(ii) = work(i);
            // last value
            if (i == 0) 
              pack_td_ptr(ntridiag) = work(npacks);
          });
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
                               const typename Impl<MatrixType>::local_ordinal_type_1d_view &packptr) {
      using Impl<MatrixType>::tpetra_device_type device_type;
      using Impl<MatrixType>::local_ordinal_type local_ordinal_type;
      using Impl<MatrixType>::size_type size_type;
      using Impl<MatrixType>::impl_scalar_type impl_scalar_type;

      const local_ordinal_type blocksize = btdm.values.extent(1);
      Kokkos::parallel_for
        (Kokkos::RangePolicy<device_type>(0, packptr.extent(0) - 1),
         KOKKOS_LAMBDA(const local_ordinal_type k) {
          for (size_type i=btdm.pack_td_ptr(packptr(k)),iend=btdm.pack_td_ptr(packptr(k+1));i<iend;i += 3)
            for (local_ordinal_type j=0;j<blocksize;++j)
              btdm.values(i,j,j) = 1;
        });
    }

    // Manage the distributed part of the computation of residual norms.
    template<typename MatrixType>
    struct NormManager {
    public: 
      using Impl<MatrixType>::magnitude_type;
      
    private:
      bool collective_;
      int sweep_step_, bsz_, nvec_;
      std::vector<magnitude_type> wrk_;
      magnitude_type n0_;
#ifdef HAVE_MPI
      MPI_Request mpi_request_;
      MPI_Comm comm_;
#endif

    public:
      NormManager (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
        sweep_step_ = 1;
        n0_ = 0;
        collective_ = comm->getSize() > 1;
        if (collective_) {
#ifdef HAVE_MPI
          const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
          TEUCHOS_ASSERT( ! mpi_comm.is_null());
          comm_ = *mpi_comm->getRawMpiComm();
#endif
        }
      }

      // Resize the buffer to accommodate nvec vectors in the multivector, for a
      // matrix having block size block_size.
      void resize (const int& block_size, const int& nvec) {
        bsz_ = block_size;
        nvec_ = nvec;
        wrk_.resize((2*block_size + 1)*nvec);
      }
      
      // Check the norm every sweep_step sweeps.
      void setCheckFrequency (const int& sweep_step) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(sweep_step < 1, "sweep step must be >= 1");
        sweep_step_ = sweep_step;
      }
      
      // Get the buffer into which to store rank-local squared norms.
      magnitude_type* get_buffer () { return wrk_.data(); }

      // Call MPI_Iallreduce to find the global squared norms.
      void ireduce (const int& sweep, const bool force = false) {
        if ( ! force && sweep % sweep_step_) return;
        const int n = bsz_*nvec_;
        if (collective_) {
          std::copy(wrk_.begin(), wrk_.begin() + n, wrk_.begin() + n);
#ifdef HAVE_MPI
#if MPI_VERSION >= 3
          MPI_Iallreduce(wrk_.data() + n, wrk_.data(), n,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_, &mpi_request_);
#else
          MPI_Allreduce (wrk_.data() + n, wrk_.data(), n,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_);
#endif
#endif
        }
      }
      
      // Check if the norm-based termination criterion is met. tol2 is the
      // tolerance squared. Sweep is the sweep index. If not every iteration is
      // being checked, this function immediately returns false. If a check must
      // be done at this iteration, it waits for the reduction triggered by
      // ireduce to complete, then checks the global norm against the tolerance.
      bool checkDone (const int& sweep, const magnitude_type& tol2, const bool force = false) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NormManager::check_done");
#endif
        TEUCHOS_ASSERT(sweep >= 1);
        if ( ! force && (sweep - 1) % sweep_step_) return false;
        if (collective_) {
#ifdef HAVE_MPI
# if MPI_VERSION >= 3
          MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
# else
          // Do nothing.
# endif
#endif
        }
        const auto n = bsz_*nvec_;
        if (sweep == 1) {
          magnitude_type* const n0 = wrk_.data() + 2*n;
          for (int v = 0; v < nvec_; ++v) {
            const magnitude_type* const dn0 = wrk_.data() + v*bsz_;
            magnitude_type mdn0 = 0;
            for (int i = 0; i < bsz_; ++i)
              mdn0 = std::max(mdn0, dn0[i]);
            n0[v] = mdn0;
          }
          return false;
        } else {
          const auto n0 = wrk_.data() + 2*n;
          bool done = true;
          for (int v = 0; v < nvec_; ++v) {
            const magnitude_type* const dnf = wrk_.data() + v*bsz_;
            magnitude_type mdnf = 0;
            for (int i = 0; i < bsz_; ++i)
              mdnf = std::max(mdnf, dnf[i]);
            if (mdnf > tol2*n0[v]) {
              done = false;
              break;
            }
          }
          return done;
        }
      }
      
      // After termination has occurred, finalize the norms for use in
      // get_norms{0,final}.
      void finalize () {
        for (int v = 0; v < nvec_; ++v) {
          const magnitude_type* const dnf = wrk_.data() + v*bsz_;
          magnitude_type mdnf = 0;
          for (int i = 0; i < bsz_; ++i)
            mdnf = std::max(mdnf, dnf[i]);
          // This overwrites the receive buffer, but that's OK; at the time of
          // this write, we no longer need the data in this slot.
          wrk_[v] = mdnf;
        }
        for (int i = 0; i < nvec_; ++i)
          wrk_[i] = std::sqrt(wrk_[i]);
        magnitude_type* const nf = wrk_.data() + 2*bsz_*nvec_;
        for (int v = 0; v < nvec_; ++v)
          nf[v] = std::sqrt(nf[v]);
      }
      
      // Report norms to the caller.
      const magnitude_type* getNorms0 () const { return wrk_.data() + 2*bsz_*nvec_; }
      const magnitude_type* getNormsFinal () const { return wrk_.data(); }
    };

//     // Map SIMD-packed multivectors to Tpetra::MultiVectors, accounting for
//     // preconditioner partitions. At the same time (for efficiency), compute the
//     // rank-local squared norms, if requested.
//     class PermuteAndRepack {
//     private:
//       enum Enum { left, right };
      
//       typename mv_type::dual_view_type::t_dev mv;
//       typename BlockTridiags::PackedMultiVector pmv;
//       Enum side;
//       impl_scalar_type damping_factor;
//       ConstUnmanaged<LOList> part2rowidx0, part2packrowidx0, lclrow, packptr;
//       bool y_is_zero;
//       magnitude_type* norm2sqr;
      
//       void init (const mv_type& mv_, const typename BlockTridiags::PackedMultiVector& pmv_,
//                  const Enum& side_, const impl_scalar_type& damping_factor_, const LOList& part2rowidx0_,
//                  const LOList& part2packrowidx0_, const LOList& lclrow_, const LOList& packptr_,
//                  const bool y_is_zero_ = false, magnitude_type* norm2sqr_ = nullptr) {
//         mv = mv_.template getLocalView<tpetra_device_type>();
//         pmv = pmv_;
//         value_count = pmv.extent(1)*pmv.extent(2);
//         side = side_;
//         damping_factor = damping_factor_;
//         part2rowidx0 = part2rowidx0_;
//         part2packrowidx0 = part2packrowidx0_;
//         lclrow = lclrow_;
//         packptr = packptr_;
//         y_is_zero = y_is_zero_;
//         norm2sqr = norm2sqr_;
//         TEUCHOS_ASSERT(lclrow.extent(0)*pmv.extent(1) == mv.extent(0));
//       }
      
//       template <typename Tag>
//       void run () {
//         const Kokkos::RangePolicy<Tag, typename device_type::execution_space>
//           rp(0, packptr.size() - 1);
//         if (norm2sqr) {
//           for (int i = 0; i < value_count; ++i)
//             norm2sqr[i] = Kokkos::ArithTraits<magnitude_type>::zero();
//           Kokkos::parallel_reduce(rp, *this, norm2sqr);
//         } else
//           Kokkos::parallel_for(rp, *this);
//       }
      
//     public:
//       struct LeftTag {};
//       struct LeftDFTag {};
//       struct RightTag {};
      
//       // For || reduce.
//       typedef magnitude_type value_type[];
//       int value_count;
      
//       // Packed multivector <- Tpetra::MultiVector.
//       PermuteAndRepack (const mv_type& mv_, typename BlockTridiags::PackedMultiVector& pmv_,
//                         const LOList& part2rowidx0_, const LOList& part2packrowidx0_,
//                         const LOList& lclrow_, const LOList& packptr_)
//       {
//         const auto one = Kokkos::ArithTraits<magnitude_type>::one();
//         init(mv_, pmv_, right, one, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_);
//       }
      
//       // Tpetra::MultiVector <- packed multivector.
//       PermuteAndRepack (const typename BlockTridiags::PackedMultiVector& pmv_, mv_type& mv_,
//                         const impl_scalar_type& damping_factor_, const LOList& part2rowidx0_,
//                         const LOList& part2packrowidx0_, const LOList& lclrow_, const LOList& packptr_,
//                         const bool y_is_zero_, magnitude_type* norm2sqr_ = nullptr)
//       {
//         init(mv_, pmv_, left, damping_factor_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_,
//              y_is_zero_, norm2sqr_);
//       }
      
//       KOKKOS_INLINE_FUNCTION void init (magnitude_type* dst) const {
//         for (int i = 0; i < value_count; ++i) dst[i] = 0;
//       }
      
//       KOKKOS_INLINE_FUNCTION
//       void join (volatile magnitude_type* dst, const volatile magnitude_type* src) const {
//         for (int i = 0; i < value_count; ++i) dst[i] += src[i];
//       }
      
//       template <typename Tag> KOKKOS_INLINE_FUNCTION
//       void operator() (const Tag&, const LO& packidx, magnitude_type* const nr = nullptr) const {
//         constexpr auto vl = vectorization_traits::vector_length;
//         const LO bs = pmv.extent(1), nrhs = pmv.extent(2);
//         LO partidx = packptr(packidx);
//         LO npacks = packptr(packidx+1) - partidx;
//         const LO pri0 = part2packrowidx0(partidx);
//         LO ri0[vl] = {0}, nrows[vl] = {0};
//         for (LO vi = 0; vi < npacks; ++vi, ++partidx) {
//           ri0[vi] = part2rowidx0(partidx);
//           nrows[vi] = part2rowidx0(partidx+1) - ri0[vi];
//         }
//         for (LO j = 0; j < nrows[0]; ++j) {
//           for (LO vi = 1; vi < npacks; ++vi)
//             if (j == nrows[vi]) {
//               npacks = vi;
//               break;
//             }
//           const LO pri = pri0 + j;
//           for (LO col = 0; col < nrhs; ++col) {
//             // In both cases, the written memory is traversed as contiguously as
//             // possible. The if statements should compile away since for a given
//             // Tag, the condition is always true or false at compile time.
//             if (std::is_same<Tag, RightTag>::value) {
//               for (LO i = 0; i < bs; ++i)
//                 for (LO vi = 0; vi < npacks; ++vi)
//                   Details::Batched::fastidx(pmv, Kokkos::make_pair(pri, vi), i, col) =
//                     mv(bs*lclrow(ri0[vi] + j) + i, col);
//             } else {
//               if (nr) {
//                 for (LO vi = 0; vi < npacks; ++vi) {
//                   const LO lr0 = bs*lclrow(ri0[vi] + j);
//                   const auto pair = Kokkos::make_pair(pri, vi);
//                   for (LO i = 0; i < bs; ++i) {
//                     auto& y = mv(lr0 + i, col);
//                     const auto& yc = Details::Batched::fastidx(pmv, pair, i, col);
//                     auto d = yc;
//                     if ( ! y_is_zero) d -= y;
//                     nr[bs*col + i] += BlockTriDiContainerDetails::abs2(d);
//                     if (std::is_same<Tag, LeftTag>::value)
//                       y = yc;
//                     else if (std::is_same<Tag, LeftDFTag>::value) {
//                       if (y_is_zero)
//                         y = damping_factor * yc;
//                       else
//                         y += damping_factor * d;
//                     }
//                   }
//                 }
//               } else {
//                 for (LO vi = 0; vi < npacks; ++vi) {
//                   const LO lr0 = bs*lclrow(ri0[vi] + j);
//                   const auto pair = Kokkos::make_pair(pri, vi);
//                   if (std::is_same<Tag, LeftTag>::value) {
//                     for (LO i = 0; i < bs; ++i)
//                       mv(lr0 + i, col) = Details::Batched::fastidx(pmv, pair, i, col);
//                   } else {
//                     for (LO i = 0; i < bs; ++i) {
//                       auto& y = mv(lr0 + i, col);
//                       const auto& yc = Details::Batched::fastidx(pmv, pair, i, col);
//                       if (y_is_zero)
//                         y = damping_factor * yc;
//                       else
//                         y += damping_factor * (yc - y);
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
      
//       void run () {
// #ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
//         TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::PermuteAndRepack");
// #endif
//         if (side == right)
//           run<RightTag>();
//         else if (damping_factor == Kokkos::ArithTraits<impl_scalar_type>::one())
//           run<LeftTag>();
//         else
//           run<LeftDFTag>();
//       }
//     };
    
    // Repack from the BlockCrsMatrix to the tridiagonal storage, and factorize
    // the data in the latter.
    template<typename MatrixType> 
    class ExtractAndFactorizeTridiags {
    public:
      using Impl<MatrixType>::local_ordinal_type;
      /// views
      using Impl<MatrixType>::local_ordinal_type_1d_view;
      using Impl<MatrixType>::size_type_1d_view; // tpetra row map view
      using Impl<MatrixType>::impl_scalar_type_1d_viewl; // tpetra value view
      /// vectorization 
      using Impl<MatrixType>::vector_type;
      using Impl<MatrixType>::vector_length;
      /// team policy member type (used in cuda)
      using typename Kokkos::TeamPolicy<device_type>::member_type;      
 
   private:
      // part interface
      ConstUnmanaged<local_ordinal_type_1d_view> partptr, lclrow, packptr;
      // block crs matrix
      ConstUnmanaged<size_type_1d_view> A_rowptr;
      ConstUnmanaged<impl_scalar_type_1d_viewl> A_values;
      // block tridiags 
      ConstUnmanaged<size_type_1d_view> D_pack_td_ptr, D_flat_td_ptr;
      ConstUnmanaged<local_ordinal_type_1d_view> D_colindsub;
      // block tridiags values
      Unmanaged<vector_type_3d_view> D_vector_values;
      Unmanaged<impl_scalar_type_4d_view> D_scalar_values;
      // diagonal safety
      const magnitude_type tiny;
      // shared information
      const local_ordinal_type blocksize, blocksize_square;

    public:
      ExtractAndFactorizeTridiags(const BlockTridiags<MatrixType> &btdm_, 
                                  const PartInterface<MatrixType> &interf_,
                                  const tpetra_block_crs_matrix_type &A_,
                                  const magnitude_type& tiny_) : 
        // interface
        partptr(interf_.partptr), 
        lclrow(interf_.lclrow), 
        packptr(interf_.packptr),
        // block crs matrix
        A_rowptr(A_.getLocalGraph().row_map), 
        A_values(A_.template getValues<typename device_type::memory_space>()),
        // block tridiags 
        D_flat_td_ptr(btdm_.flat_td_ptr), 
        D_pack_td_ptr(btdm_.pack_td_ptr), 
        D_vector_values(btdm_.values),
        D_scalar_values(impl_scalar_type_4d_view(btdm_.values.data(), 
                                                 btdm_.values.extent(0),
                                                 btdm_.values.extent(1),
                                                 btdm_.values.extent(2),
                                                 vector_length)),
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
          for (int e=0;e<3;++e) {
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
      extract(const member_type &member,  const vi, const local_ordinal_type &packidx) const {
        local_ordinal_type partidx = packptr(packidx);
        local_ordinal_type npacks = packptr(packidx+1) - partidx;

        const size_type kps = pack_td_ptr(partidx);
        local_ordinal_type kfs, ri0, nrows;

        if (vi < npacks) {
          kfs = D_flat_td_ptr(partidx+vi);
          ri0 = partptr(partidx+vi);
          nrows = partptr(partidx+1) - ri0;

          for (local_ordinal_type tr=0,j=0;tr<nrows[0];++tr) {
            for (int e=0;e<3;++e) {
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
      factorize (const LO& packidx) const {
        using namespace KokkosBatched::Experimental;

        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();

        // subview pattern
        auto A = Kokkos::subview(values, 0, Kokkkos::ALL(), Kokkos::ALL());
        auto B = A;
        auto C = A;

        auto i0 = pack_td_ptr(packptr(packidx));
        const local_ordinal_type nrows 
          = BlockTridiags::IndexToRow(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;

        A.assign_data( &values(i0,0,0) );
        SerialLU<Algo::LU::Unblocked>::invoke(A, diag_safety);
        for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
          B.assign_data( &values(i0+1,0,0) );
          SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsm::Blocked>
            ::invoke(one, A, B);
          C.assign_data( &values(i0+2,0,0) );
          SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsm::Blocked>
            ::invoke(one, A, C);
          A.assign_data( &values(i0+3,0,0) );
          SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>
            ::invoke(-one, C, B, one, A);
          SerialLU<Algo::LU::Unblocked>::invoke(A, diag_safety);
        }
      }


      KOKKOS_INLINE_FUNCTION void factorize (const typename Kokkos::TeamPolicy<device_type>::member_type &member,
                                             const local_ordinal_type &idx,
                                             const local_ordinal_type &packidx) const {
        using namespace KokkosBatched::Experimental;

        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();

        // subview pattern (layout right always; thus, right most index is the vector index)
        auto A = Kokkos::subview(flat_values, 0, Kokkkos::ALL(), Kokkos::ALL(), 0);
        auto B = A;
        auto C = A;
        
        auto i0 = pack_td_ptr(packptr(packidx));
        const local_ordinal_type nrows 
          = BlockTridiags::IndexToRow(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;

        A.assign_data( &flat_values(i0,0,0,idx) );
        
        TeamLU<Algo::LU::Unblocked>::invoke(member, A, diag_safety);
        for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
          B.assign_data( &flat_values(i0+1,0,0,idx) );
          TeamTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsm::Blocked>
            ::invoke(member, one, A, B);
          C.assign_data( &flat_values(i0+2,0,0,idx) );
          TeamTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsm::Blocked>
            ::invoke(member, one, A, C);
          A.assign_data( &flat_values(i0+3,0,0,idx) );
          TeamGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>
            ::invoke(member, -one, C, B, one, A);
          TeamLU<Algo::LU::Unblocked>::invoke(member, A, diag_safety);
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
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const typename Kokkos::TeamPolicy<device_type>::member_type &member) const {
        const local_ordinal_type packidx = member.league_rank();
        Kokkos::parallel_For
          (Kokkos::ThreadVectorRange(member, vector_length),
           [&](const int idx) {
            extract(member, idx, packidx);
            factorize(member, idx, packidx);
          });
      }
      
      void run() {
#if defined(HAVE_IFPACK2_CUDA) && defined (KOKKOS_ENABLE_CUDA)
        Kokkos::TeamPolicy<device_type> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vector_length);
#else
        Kokkos::RangePolicy<device_type> policy(0, packptr.extent(0) - 1);
#endif
        Kokkos::parallel_for(policy, *this);
      }
    };

      // Solve D X = Y, where A + D + R. This general implementation uses subviews
      // and works on any platform.
      template <typename Layout, typename PartialSpecialization>
      class Solve {
      public: // for Cuda
        // Whether RHS is a vector or a multivector.
        struct VecTag {};
        struct MatTag {};

      private:
        ConstUnmanaged<SizeList> pack_td_ptr;
        ConstUnmanaged<typename BlockTridiags::Values> values;
        ConstUnmanaged<LOList> packptr, part2packrowidx0;
        Unmanaged<typename BlockTridiags::PackedMultiVector> X;

        template <typename Tag> void run () {
          Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
          Kokkos::parallel_for(rp, *this);
        }

        // X := a T \ X
        template <typename Tag, typename UploType, typename DiagType,
                  typename Scalar, typename MatA, typename MatX>
        KOKKOS_FORCEINLINE_FUNCTION static void
        trsmv (const Scalar& a, const MatA& A, MatX& X,
               typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsv<UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsv::Unblocked>
            ::invoke(a, A, X);
        }

        template <typename Tag, typename UploType, typename DiagType,
                  typename Scalar, typename MatA, typename MatX>
        KOKKOS_FORCEINLINE_FUNCTION static void
        trsmv (const Scalar& a, const MatA& A, MatX& X,
               typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsm<kbe::Side::Left, UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsm::Blocked>
            ::invoke(a, A, X);
        }

        // C := b C + a A B
        template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
        KOKKOS_FORCEINLINE_FUNCTION static void
        gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
               typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialGemv<kbe::Trans::NoTranspose, kbe::Algo::Gemv::Unblocked>::invoke(a, A, B, b, C);
        }

        template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
        KOKKOS_FORCEINLINE_FUNCTION static void
        gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
               typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialGemm<kbe::Trans::NoTranspose, kbe::Trans::NoTranspose, kbe::Algo::Gemm::Blocked>
            ::invoke(a, A, B, b, C);
        }

      public:
        Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
               const typename BlockTridiags::PackedMultiVector& X_)
          : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
            packptr(packptr_), part2packrowidx0(part2packrowidx0_),
            X(X_)
        {}

        template <typename Tag>
        KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
          using Kokkos::subview;
          using Kokkos::ALL;
          namespace kbe = KokkosBatched::Experimental;

          const auto one = Kokkos::ArithTraits<magnitude_type>::one();

          const auto partidx = packptr(packidx);
          Size i0 = pack_td_ptr(partidx);
          LO r0 = part2packrowidx0(partidx);
          const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

          // Solve L x = x.
          auto A = subview(values, i0, ALL(), ALL());
          auto X1 = subview(X, r0, ALL(), ALL());
          trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X1);
          for (LO i = 1; i < nrows; ++i) {
            const auto B = subview(values, i0+2, ALL(), ALL());
            r0++;
            const auto X2 = subview(X, r0, ALL(), ALL());
            gemmv<Tag>(-one, B, X1, one, X2);
            i0 += 3;
            A = subview(values, i0, ALL(), ALL());
            trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X2);
            X1 = X2;
          }

          // Solve U x = x.
          trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X1);
          for (LO i = nrows; i > 1; --i) {
            i0 -= 3;
            const auto B = subview(values, i0+1, ALL(), ALL());
            r0--;
            const auto X2 = subview(X, r0, ALL(), ALL());
            gemmv<Tag>(-one, B, X1, one, X2);
            A = subview(values, i0, ALL(), ALL());
            trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X2);
            X1 = X2;
          }
        }

        void run () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
          TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
#endif
          if (X.extent(2) == 1)
            run<VecTag>();
          else
            run<MatTag>();
        }
      };

      // This specialization is for CPU and KNL. On KNL, it speeds up the solve by
      // ~20% for block size 5, more for smaller, less for bigger. The importance of
      // the apply phase of a preconditioner justifies the extra complexity of this
      // specialization.
      template <typename PartialSpecialization>
      class Solve<Kokkos::LayoutRight, PartialSpecialization> {
      public: // for Cuda
        // Whether RHS is a vector or a multivector.
        struct VecTag {};
        struct MatTag {};

      private:
        ConstUnmanaged<SizeList> pack_td_ptr;
        ConstUnmanaged<typename BlockTridiags::Values> values;
        ConstUnmanaged<LOList> packptr, part2packrowidx0;
        Unmanaged<typename BlockTridiags::PackedMultiVector> X;
        const int bs2, xsz;

        template <typename Tag> void run () {
          Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
          Kokkos::parallel_for(rp, *this);
        }

        // For speed in this important, O(nnz) operation, we don't use subviews, but
        // instead go directly to the raw-array interface.

        // X := a T \ X

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
                 typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsvInternalLower<kbe::Algo::Trsv::Unblocked>::invoke(
                                                                           true,
                                                                           values.dimension_1(),
                                                                           a,
                                                                           values.data() + i0*bs2, values.stride_1(), values.stride_2(),
                                                                           X.data() + r0*xsz, X.stride_1());
        }

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
                 typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsmInternalLeftLower<kbe::Algo::Trsm::Blocked>::invoke(
                                                                             true,
                                                                             values.dimension_1(), X.dimension_2(),
                                                                             a,
                                                                             values.data() + i0*bs2, values.stride_1(), values.stride_2(),
                                                                             X.data() + r0*xsz, X.stride_1(), X.stride_2());
        }

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
                 typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsvInternalUpper<kbe::Algo::Trsv::Unblocked>::invoke(
                                                                           false,
                                                                           values.dimension_1(),
                                                                           a,
                                                                           values.data() + i0*bs2, values.stride_1(), values.stride_2(),
                                                                           X.data() + r0*xsz, X.stride_1());
        }

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
                 typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialTrsmInternalLeftUpper<kbe::Algo::Trsm::Blocked>::invoke(
                                                                             false,
                                                                             values.dimension_1(), X.dimension_2(),
                                                                             a,
                                                                             values.data() + i0*bs2, values.stride_1(), values.stride_2(),
                                                                             X.data() + r0*xsz, X.stride_1(), X.stride_2());
        }

        // C := b C + a A B

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
               typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialGemvInternal<kbe::Algo::Gemv::Unblocked>::invoke(
                                                                      values.dimension_1(), values.dimension_2(),
                                                                      a,
                                                                      values.data() + a0*bs2, values.stride_1(), values.stride_2(),
                                                                      X.data() + b0*xsz, X.stride_1(),
                                                                      b,
                                                                      X.data() + c0*xsz, X.stride_1());
        }

        template <typename Tag>
        KOKKOS_FORCEINLINE_FUNCTION void
        gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
               typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
          namespace kbe = KokkosBatched::Experimental;
          kbe::SerialGemmInternal<kbe::Algo::Gemm::Blocked>::invoke(
                                                                    X.dimension_1(), X.dimension_2(), values.dimension_2(),
                                                                    a,
                                                                    values.data() + a0*bs2, values.stride_1(), values.stride_2(),
                                                                    X.data() + b0*xsz, X.stride_1(), X.stride_2(),
                                                                    b,
                                                                    X.data() + c0*xsz, X.stride_1(), X.stride_2());
        }

      public:
        Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
               const typename BlockTridiags::PackedMultiVector& X_)
          : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
            packptr(packptr_), part2packrowidx0(part2packrowidx0_),
            X(X_), bs2(values.dimension_1()*values.dimension_1()),
            xsz(values.dimension_1()*X.dimension_2())
        {}

        template <typename Tag>
        KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
          using Kokkos::subview;
          using Kokkos::ALL;
          namespace kbe = KokkosBatched::Experimental;

          const auto one = Kokkos::ArithTraits<magnitude_type>::one();

          const auto partidx = packptr(packidx);
          Size i0 = pack_td_ptr(partidx);
          LO r0 = part2packrowidx0(partidx);
          const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

          // Solve L x = x.
          trsmvlo<Tag>(one, i0, r0);
          for (LO i = 1; i < nrows; ++i) {
            r0++;
            gemmv<Tag>(-one, one, i0+2, r0-1, r0);
            i0 += 3;
            trsmvlo<Tag>(one, i0, r0);
          }

          // Solve U x = x.
          trsmvup<Tag>(one, i0, r0);
          for (LO i = nrows; i > 1; --i) {
            i0 -= 3;
            r0--;
            gemmv<Tag>(-one, one, i0+1, r0+1, r0);
            trsmvup<Tag>(one, i0, r0);
          }
        }

        void run () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
          TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
#endif
          if (X.extent(2) == 1)
            run<VecTag>();
          else
            run<MatTag>();
        }
      };

    private:

      template <typename T> static KOKKOS_INLINE_FUNCTION constexpr T square (const T& x) { return x*x; }


      // Top-level Impl object initialization.
      void init (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                 const Teuchos::RCP<const import_type>& importer, int overlapLevel,
                 bool overlapCommAndComp, bool useSeqMethod) {
        seq_method_ = useSeqMethod;
        overlap_comm_ = overlapCommAndComp;
        validate_ = true;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(
                                    overlapLevel != 0,
                                    "BlockTriDiContainer does not curently support OverlapLevel != 0; user provided "
                                    << overlapLevel);

        A_bcrs_ = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(matrix);
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
                                    A_bcrs_.is_null(), "BlockTriDiContainer currently supports Tpetra::BlockCrsMatrix only.");

        init_importer(importer);
        init_parts(partitions);
        init_btdm();
      }

      }



      void find_col2row (std::vector<LO>& col2row) {
        Teuchos::RCP<const map_type> rowmap, colmap, dommap; {
          const auto& g = A_bcrs_->getCrsGraph();
          rowmap = g.getRowMap();
          colmap = g.getColMap();
          dommap = g.getDomainMap();
          TEUCHOS_ASSERT( ! (rowmap.is_null() || colmap.is_null(), dommap.is_null()));
        }
        const LO nrows = partptr_(partptr_.size() - 1);
        col2row.resize(A_bcrs_->getNodeNumCols(), Teuchos::OrdinalTraits<LO>::invalid());
        for (LO lr = 0; lr < nrows; ++lr) {
          const GO gid = rowmap->getGlobalElement(lr);
          TEUCHOS_ASSERT(gid != Teuchos::OrdinalTraits<GO>::invalid());
          if ( ! dommap->isNodeGlobalElement(gid)) continue;
          const LO lc = colmap->getLocalElement(gid);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(
                                      lc == Teuchos::OrdinalTraits<LO>::invalid(), get_msg_prefix() << "GID " << gid
                                      << " gives an invalid local column.");
          col2row[lc] = lr;
        }
      }

      // Wrappers to the compute_b_minus_Rx implementations.
      // Y := B - R X, where A = D + R.
      void compute_b_minus_Rx (const mv_type& B, const mv_type& X, mv_type& Y) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
        const auto Bv = B.template getLocalView<tpetra_device_type>();
        const auto Xv = X.template getLocalView<tpetra_device_type>();
        const auto Yv = Y.template getLocalView<tpetra_device_type>();
        if (amd_.is_tpetra_block_crs) {
          const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
          BlockTriDiContainerDetails::compute_b_minus_Rx(Bv, Xv, Yv, A_bcrs_->getBlockSize(),
                                                         amd_.rowptr, amd_.A_colindsub,
                                                         g.row_map, g.entries, amd_.tpetra_values);
        }
      }

      template <typename MvView>
      void compute_b_minus_Rx (const mv_type& B, const MvView& Xv,
                               typename BlockTridiags::PackedMultiVector& Y,
                               const bool first_apply) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
        const auto Bv = B.template getLocalView<tpetra_device_type>();
        if (amd_.is_tpetra_block_crs) {
          const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
          BlockTriDiContainerDetails::compute_b_minus_Rx(
                                                         first_apply ? amd_.rowptr : amd_.rowptr_remote,
                                                         first_apply ? amd_.A_colindsub : amd_.A_colindsub_remote,
                                                         g.row_map, g.entries, amd_.tpetra_values,
                                                         Bv, Xv, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_, first_apply);
        }
      }

      template <typename MvView>
      void compute_b_minus_Rx (const mv_type& B, const MvView& X_owned, const MvView& X_remote,
                               typename BlockTridiags::PackedMultiVector& Y) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
        const auto Bv = B.template getLocalView<tpetra_device_type>();
        if (amd_.is_tpetra_block_crs) {
          const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
          BlockTriDiContainerDetails::compute_b_minus_Rx(
                                                         amd_.rowptr, amd_.A_colindsub,
                                                         g.row_map, g.entries, amd_.tpetra_values,
                                                         Bv, X_owned, X_remote, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_);
        }
      }

    public:
      Impl (BlockTriDiContainer<MatrixType, LocalScalarType>& container,
            const Teuchos::RCP<const row_matrix_type>& matrix,
            const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
            const Teuchos::RCP<const import_type>& importer, int overlapLevel,
            bool overlapCommAndComp = false, bool useSeqMethod = false)
        : container_(container)
      {
        init(matrix, partitions, importer, overlapLevel, overlapCommAndComp, useSeqMethod);
      }

      std::string describe () const {
        std::stringstream ss;
        ss << "seq_method " << seq_method_
           << " overlap_comm " << overlap_comm_
           << " dm2cm " << (dm2cm_.data() ? true : false);
        return ss.str();
      }

      // Top-level symbolic phase.
      void symbolic () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::symbolic");
#endif

        const LO nrows = partptr_(partptr_.size() - 1);
        const auto lclrow = create_host_mirror_view_and_sync(lclrow_);

        std::vector<LO> col2row;
        find_col2row(col2row);

        col_contiguous_ = col2row[0] == 0;
        if (col_contiguous_) {
          for (LO lr = 1; lr < nrows; ++lr)
            if (lclrow(col2row[lr-1]) + 1 != lclrow(col2row[lr])) {
              col_contiguous_ = false;
              break;
            }
        }

        { // Construct the D and R graphs in A = D + R.
          const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
          const auto& g_rowptr = g.row_map;
          TEUCHOS_ASSERT(g_rowptr.size() == static_cast<size_t>(nrows + 1));
          const auto& g_colind = g.entries;

          const auto rowidx2part = create_host_mirror_view_and_sync(rowidx2part_);
          const auto part2rowidx0 = create_host_mirror_view_and_sync(part2rowidx0_);

          //assume No overlap.
          std::vector<LO> lclrow2idx(nrows);
          for (LO i = 0; i < nrows; ++i)
            lclrow2idx[lclrow(i)] = i;

          // Count (block) nnzs in D and R.
          Size D_nnz = 0, R_nnz_owned = 0, R_nnz_remote = 0;
          for (LO lr = 0; lr < nrows; ++lr) {
            // LID -> index.
            const LO ri0 = lclrow2idx[lr];
            const LO pi0 = rowidx2part(ri0);
            for (Size j = g_rowptr(lr); j < g_rowptr(lr+1); ++j) {
              const LO lc = g_colind(j);
              const LO lc2r = col2row[lc];
              bool incr_R = false;
              do { // breakable
                if (lc2r == Teuchos::OrdinalTraits<LO>::invalid()) {
                  incr_R = true;
                  break;
                }
                const LO ri = lclrow2idx[lc2r];
                const LO pi = rowidx2part(ri);
                if (pi != pi0) {
                  incr_R = true;
                  break;
                }
                // Test for being in the tridiag. This is done in index space. In
                // LID space, tridiag LIDs in a row are not necessarily related by
                // {-1, 0, 1}.
                if (ri0 + 1 >= ri && ri0 <= ri + 1)
                  ++D_nnz;
                else
                  incr_R = true;
              } while (0);
              if (incr_R) {
                if (lc < nrows) ++R_nnz_owned;
                else ++R_nnz_remote;
              }
            }
          }
          if (! overlap_comm_) {
            R_nnz_owned += R_nnz_remote;
            R_nnz_remote = 0;
          }

          { // Construct the D graph.
            btdm_.A_colindsub = LOList("btdm.A_colindsub", D_nnz);
            const auto D_A_colindsub = Kokkos::create_mirror_view(btdm_.A_colindsub);
            if (validate_)
              Kokkos::deep_copy(D_A_colindsub, Teuchos::OrdinalTraits<LO>::invalid());
            const auto partptr = create_host_mirror_view_and_sync(partptr_);
            const auto flat_td_ptr = create_host_mirror_view_and_sync(btdm_.flat_td_ptr);
            const LO nparts = partptr.size() - 1;
            for (LO pi0 = 0; pi0 < nparts; ++pi0) {
              const LO part_ri0 = part2rowidx0_(pi0);
              LO offset = 0;
              for (LO ri0 = partptr(pi0); ri0 < partptr(pi0+1); ++ri0) {
                const LO td_row_os = BlockTridiags::RowToIndex(ri0 - part_ri0) + offset;
                offset = 1;
                const LO lr0 = lclrow(ri0);
                const Size j0 = g_rowptr(lr0);
                for (Size j = j0; j < g_rowptr(lr0+1); ++j) {
                  const LO lc = g_colind(j);
                  const LO lc2r = col2row[lc];
                  if (lc2r == Teuchos::OrdinalTraits<LO>::invalid()) continue;
                  const LO ri = lclrow2idx[lc2r];
                  const LO pi = rowidx2part(ri);
                  if (pi != pi0) continue;
                  if (ri + 1 < ri0 || ri > ri0 + 1) continue;
                  const LO row_entry = j - j0;
                  D_A_colindsub(flat_td_ptr(pi0) + ((td_row_os + ri) - ri0)) = row_entry;
                }
              }
            }
            if (validate_)
              for (size_t i = 0; i < D_A_colindsub.size(); ++i)
                TEUCHOS_ASSERT(D_A_colindsub(i) != Teuchos::OrdinalTraits<LO>::invalid());
            Kokkos::deep_copy(btdm_.A_colindsub, D_A_colindsub);
            { // Allocate values.
              const auto packptr = create_host_mirror_view_and_sync(packptr_);
              const LO npacks = packptr.size() - 1;
              LO nblks = 0; // Number of tridiag blocks, accounting for packing.
              for (LO pai = 0; pai < npacks; ++pai) {
                const LO pti = packptr(pai);
                const LO inrows = partptr(pti+1) - partptr(pti);
                nblks += BlockTridiags::RowToIndex(inrows);
              }
              const auto bs = A_bcrs_->getBlockSize();
              btdm_.values = typename BlockTridiags::Values("btdm.values", nblks, bs, bs);
              if (vectorization_traits::vector_length > 1)
                SetTridiagsToI(btdm_, packptr_).run();
            }
          }

          { // Construct the R graph.
            amd_.is_tpetra_block_crs = true;
            amd_.rowptr = SizeList("amd.rowptr", nrows + 1);
            amd_.A_colindsub = LOList("amd.A_colindsub", R_nnz_owned);
            const auto R_rowptr = Kokkos::create_mirror_view(amd_.rowptr);
            const auto R_A_colindsub = Kokkos::create_mirror_view(amd_.A_colindsub);
            R_rowptr(0) = 0;
            if (overlap_comm_) {
              amd_.rowptr_remote = SizeList("amd.rowptr_remote", nrows + 1);
              amd_.A_colindsub_remote = LOList("amd.A_colindsub_remote", R_nnz_remote);
            }
            const auto R_rowptr_remote = Kokkos::create_mirror_view(amd_.rowptr_remote);
            const auto R_A_colindsub_remote = Kokkos::create_mirror_view(amd_.A_colindsub_remote);
            if (overlap_comm_) R_rowptr_remote(0) = 0;
            for (LO lr = 0; lr < nrows; ++lr) {
              const LO ri0 = lclrow2idx[lr];
              const LO pi0 = rowidx2part(ri0);
              R_rowptr(lr+1) = R_rowptr(lr);
              if (overlap_comm_) R_rowptr_remote(lr+1) = R_rowptr_remote(lr);
              const Size j0 = g_rowptr(lr);
              for (Size j = j0; j < g_rowptr(lr+1); ++j) {
                const LO lc = g_colind(j);
                const LO lc2r = col2row[lc];
                if (lc2r != Teuchos::OrdinalTraits<LO>::invalid()) {
                  const LO ri = lclrow2idx[lc2r];
                  const LO pi = rowidx2part(ri);
                  if (pi == pi0 && ri + 1 >= ri0 && ri <= ri0 + 1)
                    continue;
                }
                const LO row_entry = j - j0;
                if ( ! overlap_comm_ || lc < nrows) {
                  R_A_colindsub(R_rowptr(lr+1)) = row_entry;
                  ++R_rowptr(lr+1);
                } else {
                  R_A_colindsub_remote(R_rowptr_remote(lr+1)) = row_entry;
                  ++R_rowptr_remote(lr+1);
                }
              }
            }
            TEUCHOS_ASSERT(R_rowptr(nrows) == R_nnz_owned);
            Kokkos::deep_copy(amd_.rowptr, R_rowptr);
            Kokkos::deep_copy(amd_.A_colindsub, R_A_colindsub);
            if (overlap_comm_) {
              TEUCHOS_ASSERT(R_rowptr_remote(nrows) == R_nnz_remote);
              Kokkos::deep_copy(amd_.rowptr_remote, R_rowptr_remote);
              Kokkos::deep_copy(amd_.A_colindsub_remote, R_A_colindsub_remote);
            }
            // Allocate or view values.
            if (amd_.is_tpetra_block_crs)
              amd_.tpetra_values = (const_cast<block_crs_matrix_type*>(A_bcrs_.get())->
                                    template getValues<typename tpetra_device_type::memory_space>());
          }
        }
      }

      static void debug_print (const Tridiags& t) {
        const auto& v = t.values;
        std::stringstream ss;
        ss << "v = [";
        for (size_t pi = 0; pi < t.pack_td_ptr.size() - 1; ++pi)
          for (LO i1 = 0; i1 < vectorization_traits::vector_length; ++i1)
            for (Size ind = t.pack_td_ptr(pi); ind < t.pack_td_ptr(pi+1); ++ind) {
              const auto i = Kokkos::make_pair(ind,i1);
              for (LO j = 0; j < v.extent_int(1); ++j)
                for (LO k = 0; k < v.extent_int(2); ++k)
                  ss << " " << Details::Batched::idx(v,i,j,k);
            }
        ss << "]\n";
        std::cout << ss.str();
      }

      // Top-level numeric phase.
      void numeric (const magnitude_type add_to_diag = 0) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::numeric");
#endif
        TEUCHOS_ASSERT( ! A_bcrs_.is_null());
        const auto A_rowptr = A_bcrs_->getCrsGraph().getLocalGraph().row_map;
        ExtractAndFactorizeTridiags(btdm_, partptr_, lclrow_, packptr_, A_rowptr,
                                    amd_.tpetra_values, add_to_diag).run();
        TEUCHOS_ASSERT(amd_.is_tpetra_block_crs);
      }

      // Top-level apply phase.
      int applyInverseJacobi (const mv_type& X, mv_type& Y, const impl_scalar_type& damping_factor,
                              bool y_is_zero, const int max_num_sweeps, magnitude_type tolerance,
                              const int check_tol_every) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::applyInverseJacobi");
#endif

        const bool do_norm = tolerance > Kokkos::ArithTraits<magnitude_type>::zero();
        if (do_norm) tolerance *= tolerance;
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
                                    do_norm && seq_method_,
                                    "The seq method for applyInverseJacobi, which in any case is for developer use only, " <<
                                    "does not support norm-based termination.");
        TEUCHOS_TEST_FOR_EXCEPT_MSG(max_num_sweeps <= 0, "Maximum number of sweeps must be >= 1.");

        // Set up work space if needed.
        if (pmv_.size() == 0 || Y.getNumVectors() != pmv_.extent(2)) {
          //todo Instead of reallocating, could get a subview if nvec is smaller
          // than the capacity.
          const auto nvec = Y.getNumVectors();
          const auto bs = A_bcrs_->getBlockSize();
          if (seq_method_)
            Y_remote_ = Teuchos::rcp(new mv_type(importer_->getTargetMap(), nvec));
          const auto nblks = part2packrowidx0_back;
          pmv_ = typename BlockTridiags::PackedMultiVector("pmv_" , nblks, bs, nvec);
        }

        have_norms_ = do_norm;
        if (do_norm) {
          if ( ! norm_mgr_) norm_mgr_ = std::make_shared<NormManager>(Y.getMap()->getComm());
          norm_mgr_->resize(A_bcrs_->getBlockSize(), Y.getNumVectors());
          norm_mgr_->set_check_frequency(check_tol_every);
        }

        int sweep;
        for (sweep = 0; sweep < max_num_sweeps; ++sweep) {
          if (y_is_zero) {
            // pmv := x(lclrow)
            PermuteAndRepack(X, pmv_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_).run();

          } else {
            if (seq_method_) {
              // y := x - R y.
              Y_remote_->doImport(Y, *importer_, Tpetra::REPLACE);
              compute_b_minus_Rx(X, *Y_remote_, Y);
              // pmv := y(lclrow).
              PermuteAndRepack(Y, pmv_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_).run();

            } else {

              // Fused y := x - R y and pmv := y(lclrow).
              if (overlap_comm_ || ! async_import_) {
                // Do R y_owned followed by R y_remote.
                if (async_import_)
                  async_import_->async_setup(Y);
                compute_b_minus_Rx(X, Y.template getLocalView<tpetra_device_type>(), pmv_, true);
                if (do_norm && sweep > 0 && norm_mgr_->check_done(sweep, tolerance)) {
                  if (async_import_) async_import_->cancel();
                  break;
                }
                if (async_import_) {
                  async_import_->sync_receive();
                  compute_b_minus_Rx(X, async_import_->get_mv_remote(), pmv_, false);
                }
              } else {
                // Use our custom import, but don't overlap comm with R y_owned.
                async_import_->sync_exchange(Y);
                if (do_norm && sweep > 0 && norm_mgr_->check_done(sweep, tolerance)) break;
                compute_b_minus_Rx(X, Y.template getLocalView<tpetra_device_type>(),
                                   async_import_->get_mv_remote(), pmv_);
              }

            }
          }

          // pmv := inv(D) pmv.
          Solve<layout_type, void>(btdm_, packptr_, part2packrowidx0_, pmv_).run();

          // y(lclrow) = (b - a) y(lclrow) + a pmv, with b = 1 always.
          PermuteAndRepack(pmv_, Y, damping_factor, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_,
                           y_is_zero, do_norm ? norm_mgr_->get_buffer() : nullptr).run();

          if (do_norm) {
            if (sweep + 1 == max_num_sweeps) {
              norm_mgr_->ireduce(sweep, true);
              norm_mgr_->check_done(sweep + 1, tolerance, true);
            } else {
              norm_mgr_->ireduce(sweep);
            }
          }

          y_is_zero = false;
        }

        // Sqrt the norms for the caller's use.
        if (do_norm) norm_mgr_->finalize();

        return sweep;
      }

      // Report norms to caller.
      const magnitude_type* get_norms0 () const {
        if ( ! have_norms_) return nullptr;
        return norm_mgr_->get_norms0();
      }

      const magnitude_type* get_norms_final () const {
        if ( ! have_norms_) return nullptr;
        return norm_mgr_->get_norms_final();
      }

      int get_block_size () const { return pmv_.extent(1); }
      int get_nvec () const { return pmv_.extent(2); }
    };

    // Base class for any unimplemented typename combinations.
    template <typename MatrixType, typename LocalScalarType, typename ExeSpace>
    struct UnImpl {
      typedef typename MatrixType::scalar_type scalar_type;
      typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
      typedef typename MatrixType::local_ordinal_type local_ordinal_type;
      typedef typename MatrixType::global_ordinal_type global_ordinal_type;
      typedef typename MatrixType::node_type node_type;
      typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
      typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
      typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

      std::string describe () const { return ""; }
      int get_block_size () const { return 0; }
      int get_nvec () const { return 0; }
      const magnitude_type* get_norms0 () const { return nullptr; }
      const magnitude_type* get_norms_final () const { return nullptr; }
      void symbolic () {}
      void numeric (const magnitude_type add_to_diag = 0) {}
      int applyInverseJacobi (const mv_type& X, mv_type& Y, const scalar_type damping_factor,
                              bool y_is_zero, const int max_num_sweeps, magnitude_type tolerance,
                              const int check_tol_every) const
      { return 0; }
    };

#if defined HAVE_STOKHOS_IFPACK2
    // BlockTriDiContainer is built on KokkosBatch's SIMD type, which doesn't work
    // directly with Stokhos composite types. Hence, catch these with UnImpl.
#define IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Type)                \
    template <typename MatrixType, typename T, typename ExeSpace>       \
    struct Impl<MatrixType, Type<T>, ExeSpace> : UnImpl<MatrixType, Type<T>, ExeSpace> { \
      typedef UnImpl<MatrixType, Type<T>, ExeSpace> ui;                 \
      Impl (BlockTriDiContainer<MatrixType, typename ui::scalar_type>& container, \
            const Teuchos::RCP<const typename ui::row_matrix_type>& matrix, \
            const Teuchos::Array<Teuchos::Array<typename ui::local_ordinal_type> >& partitions, \
            const Teuchos::RCP<const typename ui::import_type>& importer, int overlapLevel, \
            bool overlapCommAndComp = false, bool useSeqMethod = false) { \
        TEUCHOS_TEST_FOR_EXCEPT_MSG(                                    \
                                    true, "BlockTriDiContainer is not currently supported for Stokhos composite types."); \
      }                                                                 \
    };
    IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::MP::Vector)
    IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::UQ::PCE)
#undef IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS
#endif

  } // namespace BlockTriDiContainerDetails

} // namespace Ifpack2

#endif
