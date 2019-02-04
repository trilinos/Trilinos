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
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemm_Team_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Team_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Team_Impl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Team_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>

#include <memory>

// need to interface this into cmake variable (or only use this flag when it is necessary)
//#define IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE
#undef  IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
#include "cuda_profiler_api.h"
#endif

namespace Ifpack2 {

  namespace BlockTriDiContainerDetails {

    ///
    /// view decorators for unmanaged and const memory
    ///
    using do_not_initialize_tag = Kokkos::ViewAllocateWithoutInitializing;

    template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
    using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
                                              MemoryTraitsType::RandomAccess |
                                              flag>;
    
    template <typename ViewType>
    using Unmanaged = Kokkos::View<typename ViewType::data_type,
                                   typename ViewType::array_layout,
                                   typename ViewType::device_type,
                                   MemoryTraits<typename ViewType::memory_traits,Kokkos::Unmanaged> >;
    template <typename ViewType>
    using Const = Kokkos::View<typename ViewType::const_data_type, 
                               typename ViewType::array_layout,
                               typename ViewType::device_type, 
                               typename ViewType::memory_traits>;
    template <typename ViewType>
    using ConstUnmanaged = Const<Unmanaged<ViewType> >;    

    ///
    /// cuda specialization
    ///
    template<typename T> struct is_cuda        { enum : bool { value = false }; };
#if defined(KOKKOS_ENABLE_CUDA)
    template<> struct is_cuda<Kokkos::Cuda>    { enum : bool { value = true  }; };
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
    operator+=(volatile ArrayValueType<T,N> &a, 
               volatile const ArrayValueType<T,N> &b) {
      for (int i=0;i<N;++i) 
        a.v[i] += b.v[i];
    }
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
      void join(value_type &dst, value_type &src) const {
        for (int i=0;i<N;++i) 
          dst.v[i] += src.v[i];
      }
      KOKKOS_INLINE_FUNCTION 
      void join(volatile value_type &dst, const volatile value_type &src) const {
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

      ///
      /// default host execution space
      ///
      typedef Kokkos::DefaultHostExecutionSpace host_execution_space;
      
      ///
      /// tpetra types
      ///
      typedef typename node_type::device_type device_type;      
      typedef typename device_type::execution_space execution_space;
      typedef typename device_type::memory_space memory_space;
      typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_multivector_type;
      typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> tpetra_map_type;
      typedef Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> tpetra_import_type;
      typedef Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_row_matrix_type;
      typedef Tpetra::BlockCrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_crs_matrix_type;
      typedef typename tpetra_block_crs_matrix_type::little_block_type tpetra_block_access_view_type;
      typedef Tpetra::BlockMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> tpetra_block_multivector_type;
      typedef typename tpetra_block_crs_matrix_type::crs_graph_type::local_graph_type local_crs_graph_type;

      ///
      /// simd vectorization
      ///
      template<typename T, int l> using Vector = KokkosBatched::Experimental::Vector<T,l>;
      template<typename T> using SIMD = KokkosBatched::Experimental::SIMD<T>;
      template<typename T, typename M> using DefaultVectorLength = KokkosBatched::Experimental::DefaultVectorLength<T,M>;

      // Kyungjoo: hansen enum does not work (don't know why)
      // enum : int { vector_length = DefaultVectorLength<impl_scalar_type,memory_space>::value };
      static constexpr int vector_length = DefaultVectorLength<impl_scalar_type,memory_space>::value;
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
    };
    
    ///
    /// setup sequential importer
    ///
    template<typename MatrixType>
    typename Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_import_type> 
    createBlockCrsTpetraImporter(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A) {
      using impl_type = ImplType<MatrixType>;
      using tpetra_map_type = typename impl_type::tpetra_map_type;
      using tpetra_mv_type = typename impl_type::tpetra_block_multivector_type;
      using tpetra_import_type = typename impl_type::tpetra_import_type;

      const auto g = A->getCrsGraph();  // tpetra crs graph object
      const auto blocksize = A->getBlockSize();
      const auto src = Teuchos::rcp(new tpetra_map_type(tpetra_mv_type::makePointMap(*g.getDomainMap(), blocksize)));
      const auto tgt = Teuchos::rcp(new tpetra_map_type(tpetra_mv_type::makePointMap(*g.getColMap()   , blocksize)));

      return Teuchos::rcp(new tpetra_import_type(src, tgt));
    }

    // Partial replacement for forward-mode MultiVector::doImport. 
    // Permits overlapped communication and computation, but also supports sync'ed. 
    // I'm finding that overlapped comm/comp can give quite poor performance on some
    // platforms, so we can't just use it straightforwardly always.

    template<typename MatrixType>
    struct AsyncableImport {
    public:
      using impl_type = ImplType<MatrixType>;

    private:
      ///
      /// MPI wrapper
      ///
#if !defined(HAVE_IFPACK2_MPI)
      typedef int MPI_Request;
      typedef int MPI_Comm;
#endif
      /// teuchos mpi type traits does not recorgnize kokkos::complex (impl_scalar_type)
      /// use scalar_type for communication data type enum
      using scalar_type = typename impl_type::scalar_type;

      template <typename T>
      static int isend(const MPI_Comm comm, const T* buf, int count, int dest, int tag, MPI_Request* ireq) {
#ifdef HAVE_IFPACK2_MPI
        MPI_Request ureq;
        const auto dt = Teuchos::Details::MpiTypeTraits<scalar_type>::getType();
        int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, comm, ireq == NULL ? &ureq : ireq);
        if (ireq == NULL) MPI_Request_free(&ureq);
        return ret;
#else
        return 0;
#endif
      }

      template <typename T>
      static int irecv(const MPI_Comm comm, T* buf, int count, int src, int tag, MPI_Request* ireq) {
#ifdef HAVE_IFPACK2_MPI
        MPI_Request ureq;
        const auto dt = Teuchos::Details::MpiTypeTraits<scalar_type>::getType();
        int ret = MPI_Irecv(buf, count, dt, src, tag, comm, ireq == NULL ? &ureq : ireq);
        if (ireq == NULL) MPI_Request_free(&ureq);
        return ret;
#else
        return 0;
#endif
      }

      static int waitany(int count, MPI_Request* reqs, int* index) {
#ifdef HAVE_IFPACK2_MPI
        return MPI_Waitany(count, reqs, index, MPI_STATUS_IGNORE);
#else
        return 0;
#endif
      }

      static int waitall(int count, MPI_Request* reqs) {
#ifdef HAVE_IFPACK2_MPI
        return MPI_Waitall(count, reqs, MPI_STATUS_IGNORE);
#else
        return 0;
#endif
      }
      
    public:
      using tpetra_map_type = typename impl_type::tpetra_map_type;
      using tpetra_import_type = typename impl_type::tpetra_import_type;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using global_ordinal_type = typename impl_type::global_ordinal_type;
      using size_type = typename impl_type::size_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;

      using host_execution_space = typename impl_type::host_execution_space;
      using int_1d_view_host = Kokkos::View<int*,host_execution_space>;
      using local_ordinal_type_1d_view_host = typename impl_type::local_ordinal_type_1d_view::HostMirror;

      using execution_space = typename impl_type::execution_space;
      using memory_space = typename impl_type::memory_space;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view;
      using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
      using impl_scalar_type_2d_view = typename impl_type::impl_scalar_type_2d_view;

#ifdef HAVE_IFPACK2_MPI
      MPI_Comm comm;
#endif
      impl_scalar_type_2d_view remote_multivector;
      local_ordinal_type blocksize;
      
      template<typename T>
      struct SendRecvPair {
        T send, recv;
      };

      // (s)end and (r)eceive data:
      SendRecvPair<int_1d_view_host> pids;           // mpi ranks
      SendRecvPair<std::vector<MPI_Request> > reqs;  // MPI_Request is pointer, cannot use kokkos view
      SendRecvPair<size_type_1d_view> offset;        // offsets to local id list and data buffer
      SendRecvPair<local_ordinal_type_1d_view> lids; // local id list
      SendRecvPair<impl_scalar_type_1d_view> buffer; // data buffer

      local_ordinal_type_1d_view dm2cm; // permutation

      // for cuda
    public:
      void setOffsetValues(const Teuchos::ArrayView<const size_t> &lens, 
                           const size_type_1d_view &offs) { 
        // wrap lens to kokkos view and deep copy to device
        Kokkos::View<size_t*,host_execution_space> lens_host(const_cast<size_t*>(lens.getRawPtr()), lens.size());
        const auto lens_device = Kokkos::create_mirror_view(memory_space(), lens_host);
        Kokkos::deep_copy(lens_device, lens_host);
        
        // exclusive scan
        const Kokkos::RangePolicy<execution_space> policy(0,offs.extent(0));
        const local_ordinal_type lens_size = lens_device.extent(0);
        Kokkos::parallel_scan
          ("AsyncableImport::RangePolicy::setOffsetValues", 
           policy, KOKKOS_LAMBDA(const local_ordinal_type &i, size_type &update, const bool &final) {
            if (final) 
              offs(i) = update;            
            update += (i < lens_size ? lens_device[i] : 0);
          });
      }

    private:
      void createMpiRequests(const tpetra_import_type &import) {
        Tpetra::Distributor &distributor = import.getDistributor();

        // copy pids from distributor
        const auto pids_from = distributor.getProcsFrom();
        pids.recv = int_1d_view_host(do_not_initialize_tag("pids recv"), pids_from.size());
        memcpy(pids.recv.data(), pids_from.getRawPtr(), sizeof(int)*pids.recv.extent(0));

        const auto pids_to = distributor.getProcsTo();
        pids.send = int_1d_view_host(do_not_initialize_tag("pids send"), pids_to.size());
        memcpy(pids.send.data(), pids_to.getRawPtr(), sizeof(int)*pids.send.extent(0));

        // mpi requests 
        reqs.recv.resize(pids.recv.extent(0)); memset(reqs.recv.data(), 0, reqs.recv.size()*sizeof(MPI_Request));
        reqs.send.resize(pids.send.extent(0)); memset(reqs.send.data(), 0, reqs.send.size()*sizeof(MPI_Request));
        
        // construct offsets
        const auto lengths_to = distributor.getLengthsTo();
        offset.send = size_type_1d_view(do_not_initialize_tag("offset send"), lengths_to.size() + 1);
        
        const auto lengths_from = distributor.getLengthsFrom();
        offset.recv = size_type_1d_view(do_not_initialize_tag("offset recv"), lengths_from.size() + 1);

        setOffsetValues(lengths_to,   offset.send);
        setOffsetValues(lengths_from, offset.recv);
      }

      void createSendRecvIDs(const tpetra_import_type &import) {
        // For each remote PID, the list of LIDs to receive.
        const auto remote_lids = import.getRemoteLIDs();
        const local_ordinal_type_1d_view_host
          remote_lids_view_host(const_cast<local_ordinal_type*>(remote_lids.getRawPtr()), remote_lids.size());
        lids.recv = local_ordinal_type_1d_view(do_not_initialize_tag("lids recv"), remote_lids.size());
        Kokkos::deep_copy(lids.recv, remote_lids_view_host);

        // For each export PID, the list of LIDs to send.
        auto epids = import.getExportPIDs();
        auto elids = import.getExportLIDs();
        TEUCHOS_ASSERT(epids.size() == elids.size());
        lids.send = local_ordinal_type_1d_view(do_not_initialize_tag("lids send"), elids.size());
        auto lids_send_host = Kokkos::create_mirror_view(lids.send);

        // naive search (not sure if pids or epids are sorted)
        for (local_ordinal_type cnt=0,i=0,iend=pids.send.extent(0);i<iend;++i) {
          const auto pid_send_value = pids.send[i];
          for (local_ordinal_type j=0,jend=epids.size();j<jend;++j)
            if (epids[j] == pid_send_value) lids_send_host[cnt++] = elids[j];
#if !defined(__CUDA_ARCH__)
          TEUCHOS_ASSERT(static_cast<size_t>(cnt) == offset.send[i+1]);
#endif
        }
        Kokkos::deep_copy(lids.send, lids_send_host);
      }

    public:
      // for cuda, all tag types are public
      struct ToBuffer {};
      struct ToMultiVector {};

      // for cuda, kernel launch should be public too.
      template<typename PackTag>
      static 
      void copy(const local_ordinal_type_1d_view &lids_,
                const impl_scalar_type_1d_view &buffer_,
                const local_ordinal_type &ibeg_, 
                const local_ordinal_type &iend_,
                const impl_scalar_type_2d_view &multivector_,
                const local_ordinal_type blocksize_) {
        const local_ordinal_type num_vectors = multivector_.extent(1);
        const local_ordinal_type mv_blocksize = blocksize_*num_vectors;
        const local_ordinal_type idiff = iend_ - ibeg_;
        const auto abase = buffer_.data() + mv_blocksize*ibeg_;
        
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          using team_policy_type = Kokkos::TeamPolicy<execution_space>;
          const team_policy_type policy(idiff, num_vectors == 1 ? 1 : 2);
          Kokkos::parallel_for
            ("AsyncableImport::TeamPolicy::copy", 
             policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {          
              const local_ordinal_type i = member.league_rank();
              Kokkos::parallel_for
                (Kokkos::TeamThreadRange(member,num_vectors),[&](const local_ordinal_type &j) {
                  auto aptr = abase + blocksize_*(i + idiff*j);
                  auto bptr = &multivector_(blocksize_*lids_(i + ibeg_), j);
                  if (std::is_same<PackTag,ToBuffer>::value) 
                    Kokkos::parallel_for
                      (Kokkos::ThreadVectorRange(member,blocksize_),[&](const local_ordinal_type &k) {
                        aptr[k] = bptr[k];
                      });
                  else
                    Kokkos::parallel_for
                      (Kokkos::ThreadVectorRange(member,blocksize_),[&](const local_ordinal_type &k) {
                        bptr[k] = aptr[k];
                      });
                });
            });
#endif
        } else {
#if defined(__CUDA_ARCH__)
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else
          if (num_vectors == 1) {
            const Kokkos::RangePolicy<execution_space> policy(ibeg_, iend_); 
            Kokkos::parallel_for
              ("AsyncableImport::RangePolicy::copy",
               policy, KOKKOS_LAMBDA(const local_ordinal_type &i) {
                auto aptr = buffer_.data() + blocksize_*i;
                auto bptr = multivector_.data() + blocksize_*lids_(i);
                if (std::is_same<PackTag,ToBuffer>::value)
                  memcpy(aptr, bptr, sizeof(impl_scalar_type)*blocksize_);
                else
                  memcpy(bptr, aptr, sizeof(impl_scalar_type)*blocksize_);
              });

	  } else {
            Kokkos::MDRangePolicy
              <execution_space, Kokkos::Rank<2>, Kokkos::IndexType<local_ordinal_type> > 
              policy( { 0, 0 }, { idiff, num_vectors } );
            
            Kokkos::parallel_for
              ("AsyncableImport::MDRangePolicy::copy", 
               policy, KOKKOS_LAMBDA(const local_ordinal_type &i,
                                     const local_ordinal_type &j) { 
                auto aptr = abase + blocksize_*(i + idiff*j);
                auto bptr = &multivector_(blocksize_*lids_(i + ibeg_), j);
                if (std::is_same<PackTag,ToBuffer>::value) 
                  for (local_ordinal_type k=0;k<blocksize_;++k) aptr[k] = bptr[k];
                else
                  for (local_ordinal_type k=0;k<blocksize_;++k) bptr[k] = aptr[k];
              });
          }
#endif
        } 
        Kokkos::fence();
      }

      void createDataBuffer(const local_ordinal_type &num_vectors) {
        const size_type extent_0 = lids.recv.extent(0)*blocksize;
        const size_type extent_1 = num_vectors;
        if (remote_multivector.extent(0) == extent_0 &&
            remote_multivector.extent(1) == extent_1) {
          // skip
        } else {
          remote_multivector = 
            impl_scalar_type_2d_view(do_not_initialize_tag("remote multivector"), extent_0, extent_1);

          const auto send_buffer_size = offset.send[offset.send.extent(0)-1]*blocksize*num_vectors;
          buffer.send = impl_scalar_type_1d_view(do_not_initialize_tag("buffer send"), send_buffer_size);

          const auto recv_buffer_size = offset.recv[offset.recv.extent(0)-1]*blocksize*num_vectors;
          buffer.recv = impl_scalar_type_1d_view(do_not_initialize_tag("buffer recv"), recv_buffer_size);
        }
      }

      AsyncableImport (const Teuchos::RCP<const tpetra_map_type>& src_map,
                       const Teuchos::RCP<const tpetra_map_type>& tgt_map, 
                       const local_ordinal_type blocksize_,
                       const local_ordinal_type_1d_view dm2cm_) {
        blocksize = blocksize_;
        dm2cm = dm2cm_;

#ifdef HAVE_IFPACK2_MPI
        const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(tgt_map->getComm());
        TEUCHOS_ASSERT(!mpi_comm.is_null());
        comm = *mpi_comm->getRawMpiComm();
#endif
        const tpetra_import_type import(src_map, tgt_map);

        createMpiRequests(import);
        createSendRecvIDs(import);
      }

      void asyncSendRecv(const impl_scalar_type_2d_view &mv) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Setup::AsyncSendRecv");
#endif
#ifdef HAVE_IFPACK2_MPI
        // constants and reallocate data buffers if necessary
        const local_ordinal_type num_vectors = mv.extent(1);
        const local_ordinal_type mv_blocksize = blocksize*num_vectors;

        // send async
        for (local_ordinal_type i=0,iend=pids.send.extent(0);i<iend;++i) {
          copy<ToBuffer>(lids.send, buffer.send, offset.send(i), offset.send(i+1), 
                         mv, blocksize);
          isend(comm, 
                buffer.send.data() + offset.send[i]*mv_blocksize,
                (offset.send[i+1] - offset.send[i])*mv_blocksize,
                pids.send[i], 
                42,
                &reqs.send[i]);
        }

        // receive async
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i) {
          irecv(comm, 
                buffer.recv.data() + offset.recv[i]*mv_blocksize,
                (offset.recv[i+1] - offset.recv[i])*mv_blocksize,
                pids.recv[i],
                42,
                &reqs.recv[i]);
        }

        // I find that issuing an Iprobe seems to nudge some MPIs into action,
        // which helps with overlapped comm/comp performance.
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i) {
          int flag;
          MPI_Status stat;
          MPI_Iprobe(pids.recv[i], 42, comm, &flag, &stat);
        }
#endif
      }

      void cancel () {
#ifdef HAVE_IFPACK2_MPI
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i)
          MPI_Cancel(&reqs.recv[i]);
#endif
      }

      void syncRecv() {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Setup::SyncRecv");
#endif
#ifdef HAVE_IFPACK2_MPI
        // receive async.
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i) {
          local_ordinal_type idx = i;
          waitany(pids.recv.extent(0), reqs.recv.data(), &idx);
          copy<ToMultiVector>(lids.recv, buffer.recv, offset.recv(idx), offset.recv(idx+1),
                              remote_multivector, blocksize);
        }
        // wait on the sends to match all Isends with a cleanup operation.
        waitall(reqs.send.size(), reqs.send.data());
#endif
      }

      void syncExchange(const impl_scalar_type_2d_view &mv) {
        asyncSendRecv(mv);
        syncRecv();
      }

      impl_scalar_type_2d_view getRemoteMultiVectorLocalView() const { return remote_multivector; }
    };

    ///
    /// setup async importer
    ///
    template<typename MatrixType>
    Teuchos::RCP<AsyncableImport<MatrixType> > 
    createBlockCrsAsyncImporter(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A) {
      using impl_type = ImplType<MatrixType>;
      using tpetra_map_type = typename impl_type::tpetra_map_type;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using global_ordinal_type = typename impl_type::global_ordinal_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;

      const auto g = A->getCrsGraph();  // tpetra crs graph object
      const auto blocksize = A->getBlockSize();
      const auto domain_map = g.getDomainMap();
      const auto column_map = g.getColMap();
      
      std::vector<global_ordinal_type> gids;
      bool separate_remotes = true, found_first = false, need_owned_permutation = false;      
      for (size_t i=0;i<column_map->getNodeNumElements();++i) {
        const global_ordinal_type gid = column_map->getGlobalElement(i);
        if (!domain_map->isNodeGlobalElement(gid)) {
          found_first = true;
          gids.push_back(gid);
        } else if (found_first) {
          separate_remotes = false;
          break;
        }
        if (!need_owned_permutation && 
            domain_map->getLocalElement(gid) != static_cast<local_ordinal_type>(i)) {
          // The owned part of the domain and column maps are different
          // orderings. We *could* do a super efficient impl of this case in the
          // num_sweeps > 1 case by adding complexity to PermuteAndRepack. But,
          // really, if a caller cares about speed, they wouldn't make different
          // local permutations like this. So we punt on the best impl and go for
          // a pretty good one: the permutation is done in place in
          // compute_b_minus_Rx for the pure-owned part of the MVP. The only cost
          // is the presumably worse memory access pattern of the input vector.
          need_owned_permutation = true;
        }
      }
      
      if (separate_remotes) {
        const auto invalid = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
        const auto parsimonious_col_map 
          = Teuchos::rcp(new tpetra_map_type(invalid, gids.data(), gids.size(), 0, domain_map->getComm()));
        if (parsimonious_col_map->getGlobalNumElements() > 0) {
          // make the importer only if needed.
          local_ordinal_type_1d_view dm2cm;
          if (need_owned_permutation) {
            dm2cm = local_ordinal_type_1d_view("dm2cm", domain_map->getNodeNumElements());
            const auto dm2cm_host = Kokkos::create_mirror_view(dm2cm);
            for (size_t i=0;i<domain_map->getNodeNumElements();++i)
              dm2cm_host(i) = domain_map->getLocalElement(column_map->getGlobalElement(i));          
            Kokkos::deep_copy(dm2cm, dm2cm_host);
          }
          return Teuchos::rcp(new AsyncableImport<MatrixType>(domain_map, parsimonious_col_map, blocksize, dm2cm));
        }
      }
      return Teuchos::null;
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
      using impl_type = ImplType<MatrixType>;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;

      constexpr int vector_length = impl_type::vector_length;

      const auto comm = A->getRowMap()->getComm();
      
      PartInterface<MatrixType> interf;
      
      const bool jacobi = partitions.size() == 0;
      const local_ordinal_type A_n_lclrows = A->getNodeNumRows();
      const local_ordinal_type nparts = jacobi ? A_n_lclrows : partitions.size();

#if defined(BLOCKTRIDICONTAINER_DEBUG)              
      local_ordinal_type nrows = 0;
      if (jacobi)       
        nrows = nparts;
      else              
        for (local_ordinal_type i=0;i<nparts;++i) nrows += partitions[i].size();

      TEUCHOS_TEST_FOR_EXCEPT_MSG
        (nrows != A_n_lclrows, get_msg_prefix(comm) << "The #rows implied by the local partition is not "
         << "the same as getNodeNumRows: " << nrows << " vs " << A_n_lclrows);
#endif
      
      // permutation vector
      std::vector<local_ordinal_type> p;
      if (!jacobi) {
        // reorder parts to maximize simd packing efficiency
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
      interf.partptr = local_ordinal_type_1d_view(do_not_initialize_tag("partptr"), nparts + 1);
      interf.lclrow = local_ordinal_type_1d_view(do_not_initialize_tag("lclrow"), A_n_lclrows);
      interf.part2rowidx0 = local_ordinal_type_1d_view(do_not_initialize_tag("part2rowidx0"), nparts + 1);
      interf.part2packrowidx0 = local_ordinal_type_1d_view(do_not_initialize_tag("part2packrowidx0"), nparts + 1);
      interf.rowidx2part = local_ordinal_type_1d_view(do_not_initialize_tag("rowidx2part"), A_n_lclrows);

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
      using impl_type = ImplType<MatrixType>;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;

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

      bool is_diagonal_only;

      BlockTridiags() = default;
      BlockTridiags(const BlockTridiags &b) = default;

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
      using impl_type = ImplType<MatrixType>;
      using execution_space = typename impl_type::execution_space;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using size_type_1d_view = typename impl_type::size_type_1d_view;

      constexpr int vector_length = impl_type::vector_length;

      BlockTridiags<MatrixType> btdm;

      const local_ordinal_type ntridiags = interf.partptr.extent(0) - 1;
      
      { // construct the flat index pointers into the tridiag values array.
        btdm.flat_td_ptr = size_type_1d_view(do_not_initialize_tag("btdm.flat_td_ptr"), ntridiags + 1);
        const Kokkos::RangePolicy<execution_space> policy(0,ntridiags + 1);
        Kokkos::parallel_scan
          ("createBlockTridiags::RangePolicy::flat_td_ptr", 
           policy, KOKKOS_LAMBDA(const local_ordinal_type &i, size_type &update, const bool &final) {
            if (final) 
              btdm.flat_td_ptr(i) = update;       
            if (i < ntridiags) {
              const local_ordinal_type nrows = interf.partptr(i+1) - interf.partptr(i);
              update += btdm.NumBlocks(nrows);
            } 
          });
        
        const auto nblocks = Kokkos::create_mirror_view_and_copy
          (Kokkos::HostSpace(), Kokkos::subview(btdm.flat_td_ptr, ntridiags));
        btdm.is_diagonal_only = (static_cast<local_ordinal_type>(nblocks()) == ntridiags);
      }
      
      // And the packed index pointers.
      if (vector_length == 1) {
        btdm.pack_td_ptr = btdm.flat_td_ptr;
      } else {
        const local_ordinal_type npacks = interf.packptr.extent(0) - 1;
        btdm.pack_td_ptr = size_type_1d_view(do_not_initialize_tag("btdm.pack_td_ptr"), ntridiags + 1);
        const Kokkos::RangePolicy<execution_space> policy(0,npacks);
        Kokkos::parallel_scan
          ("createBlockTridiags::RangePolicy::pack_td_ptr", 
           policy, KOKKOS_LAMBDA(const local_ordinal_type &i, size_type &update, const bool &final) {
            const local_ordinal_type parti      = interf.packptr(i);
            const local_ordinal_type parti_next = interf.packptr(i+1);
            if (final) {
              const size_type nblks = update;
              for (local_ordinal_type pti=parti;pti<parti_next;++pti)
                btdm.pack_td_ptr(pti) = nblks;
              const local_ordinal_type nrows = interf.partptr(parti+1) - interf.partptr(parti);
              // last one
              if (i == npacks-1) 
                btdm.pack_td_ptr(ntridiags) = nblks + btdm.NumBlocks(nrows);
            }
            {
              const local_ordinal_type nrows = interf.partptr(parti+1) - interf.partptr(parti);
              update += btdm.NumBlocks(nrows);
            }
          });        
      }

      // values and A_colindsub are created in the symbolic phase

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
      using impl_type = ImplType<MatrixType>;
      using execution_space = typename impl_type::execution_space;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;

      using size_type_1d_view = typename impl_type::size_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;

      const ConstUnmanaged<size_type_1d_view> pack_td_ptr(btdm.pack_td_ptr);
      const local_ordinal_type blocksize = btdm.values.extent(1);
      
      if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
        constexpr int vector_length = impl_type::vector_length;
        using impl_scalar_type = typename impl_type::impl_scalar_type;
        using impl_scalar_type_4d_view = typename impl_type::impl_scalar_type_4d_view;
        using team_policy_type = Kokkos::TeamPolicy<execution_space>;
        const impl_scalar_type_4d_view values((impl_scalar_type*)btdm.values.data(), 
                                              btdm.values.extent(0), 
                                              btdm.values.extent(1),
                                              btdm.values.extent(2),
                                              vector_length);
	// vector_lengh (enum or constexpr) sometimes is not captured by device lambda
	// reference should not be used but some compilers interpret vector_length as 
	// a reference
	const local_ordinal_type vector_length_value = vector_length;
        const team_policy_type policy(packptr.extent(0)-1, Kokkos::AUTO(), vector_length); 
        Kokkos::parallel_for
          ("setTridiagsToIdentity::TeamPolicy", 
           policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {          
            const local_ordinal_type k = member.league_rank();
            const local_ordinal_type ibeg = pack_td_ptr(packptr(k));
            const local_ordinal_type iend = pack_td_ptr(packptr(k+1));
            const local_ordinal_type diff = iend - ibeg;
            const local_ordinal_type icount = diff/3 + (diff%3 > 0);
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,icount),[&](const local_ordinal_type &ii) {
		Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, vector_length_value),[&](const int &v) {
                    const local_ordinal_type i = ibeg + ii*3;
                    for (local_ordinal_type j=0;j<blocksize;++j) 
                      values(i,j,j,v) = 1;
                  });
              });
          });
#endif
      } else {
        // exclude from cuda to remove warning to compile host code on device
#if defined(__CUDA_ARCH__)
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else
        const Unmanaged<vector_type_3d_view> values = btdm.values;
        const Kokkos::RangePolicy<execution_space> policy(0, packptr.extent(0) - 1);
        Kokkos::parallel_for
          ("setTridiagsToIdentity::RangePolicy", 
           policy, KOKKOS_LAMBDA(const local_ordinal_type k) {
            for (size_type i=pack_td_ptr(packptr(k)),iend=pack_td_ptr(packptr(k+1));i<iend;i+=3)
              for (local_ordinal_type j=0;j<blocksize;++j)
                values(i,j,j) = 1;
          });
#endif
      }
    }
    
    ///
    /// A - Tridiags(A), i.e., R in the splitting A = D + R.
    ///
    template <typename MatrixType>
    struct AmD {
      using impl_type = ImplType<MatrixType>;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view;
      using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;

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

      AmD() = default;
      AmD(const AmD &b) = default;
    };

    ///
    /// symbolic phase, on host : create R = A - D, pack D
    ///
    template<typename MatrixType>
    void
    performSymbolicPhase(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
                         const PartInterface<MatrixType> &interf,
                         BlockTridiags<MatrixType> &btdm,
                         AmD<MatrixType> &amd,
                         const bool overlap_communication_and_computation) {

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::SymbolicPhase");
#endif
      using impl_type = ImplType<MatrixType>;
      using memory_space = typename impl_type::memory_space;
      using host_execution_space = typename impl_type::host_execution_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using global_ordinal_type = typename impl_type::global_ordinal_type;
      using size_type = typename impl_type::size_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using block_crs_matrix_type = typename impl_type::tpetra_block_crs_matrix_type;

      constexpr int vector_length = impl_type::vector_length;

      const auto comm = A->getRowMap()->getComm();
      const auto& g = A->getCrsGraph();
      const auto blocksize = A->getBlockSize();      

      // mirroring to host
      const auto partptr = Kokkos::create_mirror_view_and_copy     (Kokkos::HostSpace(), interf.partptr);
      const auto lclrow = Kokkos::create_mirror_view_and_copy      (Kokkos::HostSpace(), interf.lclrow);
      const auto rowidx2part = Kokkos::create_mirror_view_and_copy (Kokkos::HostSpace(), interf.rowidx2part);
      const auto part2rowidx0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interf.part2rowidx0);
      const auto packptr = Kokkos::create_mirror_view_and_copy     (Kokkos::HostSpace(), interf.packptr);

      const local_ordinal_type nrows = partptr(partptr.extent(0) - 1);

      // find column to row map on host
      Kokkos::View<local_ordinal_type*,host_execution_space> col2row("col2row", A->getNodeNumCols());
      Kokkos::deep_copy(col2row, Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
      {
        const auto rowmap = g.getRowMap();
        const auto colmap = g.getColMap();
        const auto dommap = g.getDomainMap();
        TEUCHOS_ASSERT( !(rowmap.is_null() || colmap.is_null() || dommap.is_null()));

#if !defined(__CUDA_ARCH__)        
        const Kokkos::RangePolicy<host_execution_space> policy(0,nrows);
        Kokkos::parallel_for
          ("performSymbolicPhase::RangePolicy::col2row",
           policy, KOKKOS_LAMBDA(const local_ordinal_type &lr) {
            const global_ordinal_type gid = rowmap->getGlobalElement(lr);
            TEUCHOS_ASSERT(gid != Teuchos::OrdinalTraits<global_ordinal_type>::invalid());
            if (dommap->isNodeGlobalElement(gid)) {
              const local_ordinal_type lc = colmap->getLocalElement(gid);
#  if defined(BLOCKTRIDICONTAINER_DEBUG)              
              TEUCHOS_TEST_FOR_EXCEPT_MSG(lc == Teuchos::OrdinalTraits<local_ordinal_type>::invalid(), 
                                          get_msg_prefix(comm) << "GID " << gid
                                          << " gives an invalid local column.");
#  endif
              col2row(lc) = lr;
            }
          });
#endif
      }

      // construct the D and R graphs in A = D + R.
      { 
        const auto& local_graph = g.getLocalGraph();
        const auto& local_graph_rowptr = local_graph.row_map;
        TEUCHOS_ASSERT(local_graph_rowptr.size() == static_cast<size_t>(nrows + 1));
        const auto& local_graph_colidx = local_graph.entries;

        //assume no overlap.

        Kokkos::View<local_ordinal_type*,host_execution_space> lclrow2idx("lclrow2idx", nrows);
        {
          const Kokkos::RangePolicy<host_execution_space> policy(0,nrows);
          Kokkos::parallel_for
            ("performSymbolicPhase::RangePolicy::lclrow2idx",
             policy, KOKKOS_LAMBDA(const local_ordinal_type &i) {
              lclrow2idx[lclrow(i)] = i;
            });
        }

        // count (block) nnzs in D and R.
        typedef SumReducer<size_type,3,host_execution_space> sum_reducer_type;
        typename sum_reducer_type::value_type sum_reducer_value;
        {
          const Kokkos::RangePolicy<host_execution_space> policy(0,nrows);
          Kokkos::parallel_reduce
            // profiling interface does not work
            (//"performSymbolicPhase::RangePolicy::count_nnz", 
             policy, KOKKOS_LAMBDA(const local_ordinal_type &lr, typename sum_reducer_type::value_type &update) {
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
        }
        size_type D_nnz = sum_reducer_value.v[0];
        size_type R_nnz_owned = sum_reducer_value.v[1];
        size_type R_nnz_remote = sum_reducer_value.v[2];

        if (!overlap_communication_and_computation) {
          R_nnz_owned += R_nnz_remote;
          R_nnz_remote = 0;
        }

        // construct the D graph.
        { 
          const auto flat_td_ptr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), btdm.flat_td_ptr);

          btdm.A_colindsub = local_ordinal_type_1d_view("btdm.A_colindsub", D_nnz);
          const auto D_A_colindsub = Kokkos::create_mirror_view(btdm.A_colindsub);

#if defined(BLOCKTRIDICONTAINER_DEBUG)              
          Kokkos::deep_copy(D_A_colindsub, Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
#endif

          const local_ordinal_type nparts = partptr.extent(0) - 1;
          {
            const Kokkos::RangePolicy<host_execution_space> policy(0, nparts);
            Kokkos::parallel_for
              ("performSymbolicPhase::RangePolicy<host_execution_space>::D_graph",
               policy, KOKKOS_LAMBDA(const local_ordinal_type &pi0) {
                const local_ordinal_type part_ri0 = part2rowidx0(pi0);
                local_ordinal_type offset = 0;
                for (local_ordinal_type ri0=partptr(pi0);ri0<partptr(pi0+1);++ri0) {
                  const local_ordinal_type td_row_os = btdm.RowToIndex(ri0 - part_ri0) + offset;
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
          }
#if defined(BLOCKTRIDICONTAINER_DEBUG)              
          for (size_t i=0;i<D_A_colindsub.extent(0);++i)
            TEUCHOS_ASSERT(D_A_colindsub(i) != Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
#endif
          Kokkos::deep_copy(btdm.A_colindsub, D_A_colindsub);
          
          // Allocate values.
          { 
            const local_ordinal_type npacks = packptr.extent(0) - 1;
            const auto pack_td_ptr_last = Kokkos::subview(btdm.pack_td_ptr, nparts);
            const auto num_packed_blocks = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pack_td_ptr_last);
            btdm.values = vector_type_3d_view("btdm.values", num_packed_blocks(), blocksize, blocksize);
            if (vector_length > 1) setTridiagsToIdentity(btdm, interf.packptr);
          }
        }
        
        // Construct the R graph.        
        { 
          amd.rowptr = size_type_1d_view("amd.rowptr", nrows + 1);
          amd.A_colindsub = local_ordinal_type_1d_view(Kokkos::ViewAllocateWithoutInitializing("amd.A_colindsub"), R_nnz_owned);

          const auto R_rowptr = Kokkos::create_mirror_view(amd.rowptr); 
          const auto R_A_colindsub = Kokkos::create_mirror_view(amd.A_colindsub);
          
          amd.rowptr_remote = size_type_1d_view("amd.rowptr_remote", overlap_communication_and_computation ? nrows + 1 : 0);
          amd.A_colindsub_remote = local_ordinal_type_1d_view(Kokkos::ViewAllocateWithoutInitializing("amd.A_colindsub_remote"), R_nnz_remote);
          
          const auto R_rowptr_remote = Kokkos::create_mirror_view(amd.rowptr_remote);
          const auto R_A_colindsub_remote = Kokkos::create_mirror_view(amd.A_colindsub_remote);

          {
            const Kokkos::RangePolicy<host_execution_space> policy(0,nrows);
            Kokkos::parallel_for
              ("performSymbolicPhase::RangePolicy<host_execution_space>::R_graph_count",
               policy, KOKKOS_LAMBDA(const local_ordinal_type &lr) {
                const local_ordinal_type ri0 = lclrow2idx[lr];
                const local_ordinal_type pi0 = rowidx2part(ri0);
                const size_type j0 = local_graph_rowptr(lr);
                for (size_type j=j0;j<local_graph_rowptr(lr+1);++j) {
                  const local_ordinal_type lc = local_graph_colidx(j);
                  const local_ordinal_type lc2r = col2row[lc];
                  if (lc2r != Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) {
                    const local_ordinal_type ri = lclrow2idx[lc2r];
                    const local_ordinal_type pi = rowidx2part(ri);
                    if (pi == pi0 && ri + 1 >= ri0 && ri <= ri0 + 1) {
                      continue;
                    }
                  }
                  // exclusive scan will be performed later
                  if (!overlap_communication_and_computation || lc < nrows) {
                    ++R_rowptr(lr);
                  } else {
                    ++R_rowptr_remote(lr);
                  }
                }
              });
          }

          // exclusive scan
          typedef ArrayValueType<size_type,2> update_type;
          {
            Kokkos::RangePolicy<host_execution_space> policy(0,nrows+1);
            Kokkos::parallel_scan
              ("performSymbolicPhase::RangePolicy<host_execution_space>::R_graph_fill",
               policy, KOKKOS_LAMBDA(const local_ordinal_type &lr, 
                                     update_type &update, 
                                     const bool &final) {
                update_type val;
                val.v[0] = R_rowptr(lr);
                if (overlap_communication_and_computation)
                  val.v[1] = R_rowptr_remote(lr);
                
                if (final) {
                  R_rowptr(lr) = update.v[0];
                  if (overlap_communication_and_computation)
                    R_rowptr_remote(lr) = update.v[1];
                  
                  if (lr < nrows) {
                    const local_ordinal_type ri0 = lclrow2idx[lr];
                    const local_ordinal_type pi0 = rowidx2part(ri0);
                    
                    size_type cnt_rowptr = R_rowptr(lr);
                    size_type cnt_rowptr_remote = overlap_communication_and_computation ? R_rowptr_remote(lr) : 0; // when not overlap_communication_and_computation, this value is garbage
                    
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
                      if (!overlap_communication_and_computation || lc < nrows) 
                        R_A_colindsub(cnt_rowptr++) = row_entry;
                      else 
                        R_A_colindsub_remote(cnt_rowptr_remote++) = row_entry;
                    }
                  }
                } 
                update += val;
              });
          }
          TEUCHOS_ASSERT(R_rowptr(nrows) == R_nnz_owned);
          Kokkos::deep_copy(amd.rowptr, R_rowptr);
          Kokkos::deep_copy(amd.A_colindsub, R_A_colindsub);
          if (overlap_communication_and_computation) {
            TEUCHOS_ASSERT(R_rowptr_remote(nrows) == R_nnz_remote);
            Kokkos::deep_copy(amd.rowptr_remote, R_rowptr_remote);
            Kokkos::deep_copy(amd.A_colindsub_remote, R_A_colindsub_remote);
          }
          
          // Allocate or view values.
          amd.tpetra_values = (const_cast<block_crs_matrix_type*>(A.get())->
                               template getValues<memory_space>());
        }
      }
    }
    
    ///
    /// numeric phase, initialize the preconditioner
    ///
    template<typename MatrixType>
    struct ExtractAndFactorizeTridiags {
    public:
      using impl_type = ImplType<MatrixType>;
      // a functor cannot have both device_type and execution_space; specialization error in kokkos
      using execution_space = typename impl_type::execution_space;
      using memory_space = typename impl_type::memory_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
      using magnitude_type = typename impl_type::magnitude_type;
      /// tpetra interface
      using block_crs_matrix_type = typename impl_type::tpetra_block_crs_matrix_type;
      /// views
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view; 
      using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view; 
      /// vectorization 
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using impl_scalar_type_4d_view = typename impl_type::impl_scalar_type_4d_view;
      static constexpr int vector_length = impl_type::vector_length;

      /// team policy member type (used in cuda)
      using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;      

    private:
      // part interface
      const ConstUnmanaged<local_ordinal_type_1d_view> partptr, lclrow, packptr;
      // block crs matrix (it could be Kokkos::UVMSpace::size_type, which is int)
      using a_rowptr_value_type = typename Kokkos::ViewTraits<local_ordinal_type*,typename impl_type::device_type>::size_type; 
      const ConstUnmanaged<Kokkos::View<a_rowptr_value_type*,typename impl_type::device_type> > A_rowptr;
      const ConstUnmanaged<impl_scalar_type_1d_view> A_values;
      // block tridiags 
      const ConstUnmanaged<size_type_1d_view> pack_td_ptr, flat_td_ptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> A_colindsub;
      const Unmanaged<vector_type_3d_view> vector_values;
      const Unmanaged<impl_scalar_type_4d_view> scalar_values;
      // shared information
      const local_ordinal_type blocksize, blocksize_square;
      // diagonal safety
      const magnitude_type tiny;
      const local_ordinal_type vector_length_value;

    public:
      ExtractAndFactorizeTridiags(const BlockTridiags<MatrixType> &btdm_, 
                                  const PartInterface<MatrixType> &interf_,
                                  const Teuchos::RCP<const block_crs_matrix_type> &A_,
                                  const magnitude_type& tiny_) : 
        // interface
        partptr(interf_.partptr), 
        lclrow(interf_.lclrow), 
        packptr(interf_.packptr),
        // block crs matrix
        A_rowptr(A_->getCrsGraph().getLocalGraph().row_map), 
        A_values(const_cast<block_crs_matrix_type*>(A_.get())->template getValues<memory_space>()),
        // block tridiags 
        pack_td_ptr(btdm_.pack_td_ptr), 
        flat_td_ptr(btdm_.flat_td_ptr), 
        A_colindsub(btdm_.A_colindsub),
        vector_values(btdm_.values),
        scalar_values((impl_scalar_type*)vector_values.data(),
                      vector_values.extent(0),
                      vector_values.extent(1),
                      vector_values.extent(2),
                      vector_length),
        blocksize(btdm_.values.extent(1)),
        blocksize_square(blocksize*blocksize),
        // diagonal weight to avoid zero pivots
        tiny(tiny_),
	vector_length_value(vector_length) {}

    private:

      inline
      void 
      extract(const local_ordinal_type &packidx) const {
        local_ordinal_type partidx = packptr(packidx);
        local_ordinal_type npacks = packptr(packidx+1) - partidx;

        const size_type kps = pack_td_ptr(partidx);
        local_ordinal_type kfs[vector_length] = {};
        local_ordinal_type ri0[vector_length] = {};
        local_ordinal_type nrows[vector_length] = {};

        for (local_ordinal_type vi=0;vi<npacks;++vi,++partidx) {
          kfs[vi] = flat_td_ptr(partidx);
          ri0[vi] = partptr(partidx);
          nrows[vi] = partptr(partidx+1) - ri0[vi];
        }
        for (local_ordinal_type tr=0,j=0;tr<nrows[0];++tr) {
          for (local_ordinal_type e=0;e<3;++e) {
            const impl_scalar_type* block[vector_length] = {};
            for (local_ordinal_type vi=0;vi<npacks;++vi) {
              const size_type Aj = A_rowptr(lclrow(ri0[vi] + tr)) + A_colindsub(kfs[vi] + j);
              block[vi] = &A_values(Aj*blocksize_square);
            }
            const size_type pi = kps + j;
            ++j;
            for (local_ordinal_type ii=0;ii<blocksize;++ii) {
              for (local_ordinal_type jj=0;jj<blocksize;++jj) {
                const auto idx = ii*blocksize + jj;
                auto& v = vector_values(pi, ii, jj);
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
      
      inline
      void 
      factorize (const local_ordinal_type& packidx) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level3::Blocked;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
                
        auto i0 = pack_td_ptr(packptr(packidx));
        const local_ordinal_type nrows 
          = BlockTridiags<MatrixType>::IndexToRow(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;

        // subview pattern
        auto A = Kokkos::subview(vector_values, i0, Kokkos::ALL(), Kokkos::ALL());
        KB::SerialLU<KB::Algo::LU::Unblocked>::invoke(A, tiny);

        if (nrows > 1) {
          auto B = A;
          auto C = A;
          for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
            B.assign_data( &vector_values(i0+1,0,0) );
            KB::SerialTrsm<KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>
              ::invoke(one, A, B);
            C.assign_data( &vector_values(i0+2,0,0) );
            KB::SerialTrsm<KB::Side::Right,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>
              ::invoke(one, A, C);
            A.assign_data( &vector_values(i0+3,0,0) );
            KB::SerialGemm<KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
              ::invoke(-one, C, B, one, A);
            KB::SerialLU<KB::Algo::LU::Unblocked>::invoke(A, tiny);
          }
        }
      }

      KOKKOS_INLINE_FUNCTION 
      void 
      extract(const member_type &member, 
              const local_ordinal_type &partidx,
	      const local_ordinal_type &v) const {
        const size_type kps = pack_td_ptr(partidx);
        const local_ordinal_type kfs = flat_td_ptr(partidx);
        const local_ordinal_type ri0 = partptr(partidx);
        const local_ordinal_type nrows = partptr(partidx+1) - ri0;
        for (local_ordinal_type tr=0,j=0;tr<nrows;++tr) {          
          local_ordinal_type lbeg = (tr == 0         ? 1 : 0);
          local_ordinal_type lend = (tr == nrows - 1 ? 2 : 3);
          for (local_ordinal_type l=lbeg;l<lend;++l,++j) { // l == 1, diagonal
            const size_type Aj = A_rowptr(lclrow(ri0 + tr)) + A_colindsub(kfs + j);
            const impl_scalar_type* block = &A_values(Aj*blocksize_square);
            const size_type pi = kps + j;
            Kokkos::parallel_for
              (Kokkos::TeamThreadRange(member,blocksize_square), [&](const local_ordinal_type &ij) {
                const local_ordinal_type ii = ij/blocksize;
                const local_ordinal_type jj = ij%blocksize;
                scalar_values(pi, ii, jj, v) = block[ii*blocksize + jj];
              });
          }
        }
      }

      KOKKOS_INLINE_FUNCTION 
      void 
      factorize(const member_type &member, 
		const local_ordinal_type &i0,
		const local_ordinal_type &nrows,
                const local_ordinal_type &v) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level3::Unblocked;

        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();

        // subview pattern
        auto A = Kokkos::subview(scalar_values, i0, Kokkos::ALL(), Kokkos::ALL(), 0);
        A.assign_data( &scalar_values(i0,0,0,v) );
        KB::TeamLU<member_type,KB::Algo::LU::Unblocked>::invoke(member, A, tiny);
        if (nrows > 1) {
          auto B = A;
          auto C = A;
          local_ordinal_type i = i0;
          for (local_ordinal_type tr=1;tr<nrows;++tr,i+=3) {
            B.assign_data( &scalar_values(i+1,0,0,v) );
            KB::TeamTrsm<member_type,KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>
              ::invoke(member, one, A, B);
            C.assign_data( &scalar_values(i+2,0,0,v) );
            KB::TeamTrsm<member_type,KB::Side::Right,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>
              ::invoke(member, one, A, C);
            A.assign_data( &scalar_values(i+3,0,0,v) );
            KB::TeamGemm<member_type,KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
              ::invoke(member, -one, C, B, one, A);
            KB::TeamLU<member_type,KB::Algo::LU::Unblocked>::invoke(member, A, tiny);
          }
        }
      }

    public:
      ///
      /// host serial (vector intrinsic) vectorization
      ///
      inline
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
      operator() (const member_type &member) const {
	// btdm is packed and sorted from largest one 
	const local_ordinal_type packidx = member.league_rank();
	const local_ordinal_type partidx = packptr(packidx);
	const local_ordinal_type npacks = packptr(packidx+1) - partidx;
	const local_ordinal_type i0 = pack_td_ptr(partidx);
	const local_ordinal_type nrows = partptr(partidx+1) - partptr(partidx);
	
	// every cuda threads should participate in computation (do not reduce by npacks)
	// Kyungjoo Why NVCC does not understand "constexpr"
	//   if I do set enum, some compilers take the enum value as ZERO and creates silent erros
	//   for other compilers, simply it generates constexpr is undefined on device. 
	//   seriously, better compiler support is necessary from NVIDIA
	Kokkos::parallel_for
	  (Kokkos::ThreadVectorRange(member, vector_length_value), [&](const local_ordinal_type &v) {
	    // extract team should not include any team barrier as it is masked out by npacks
	    if (v < npacks) extract(member, partidx + v, v);
	    member.team_barrier();
	    // factorize should be invoked by all threads as team kernels invokes barrier inside
	    factorize(member, i0, nrows, v);
	  });
      }

      void run() {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStart();
#endif

        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
	  const local_ordinal_type vl = vector_length;
          const Kokkos::TeamPolicy<execution_space> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl);
          Kokkos::parallel_for("ExtractAndFactorize::TeamPolicy::run", policy, *this);
#endif
        } else {
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else          
          const Kokkos::RangePolicy<execution_space> policy(0, packptr.extent(0) - 1);
          Kokkos::parallel_for("ExtractAndFactorize::RangePolicy::run", policy, *this);
#endif
        }
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStop();
#endif
      }
    }; 
    
    ///
    /// top level numeric interface
    ///
    template<typename MatrixType>
    void
    performNumericPhase(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
                        const PartInterface<MatrixType> &interf,
                        BlockTridiags<MatrixType> &btdm,
                        const typename ImplType<MatrixType>::magnitude_type tiny) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NumericPhase");
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
      using execution_space = typename impl_type::execution_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
      using magnitude_type = typename impl_type::magnitude_type;
      using tpetra_multivector_type = typename impl_type::tpetra_multivector_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using impl_scalar_type_2d_view = typename impl_type::impl_scalar_type_2d_view;
      static constexpr int vector_length = impl_type::vector_length;

      using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;
      
    private:
      // part interface
      const ConstUnmanaged<local_ordinal_type_1d_view> partptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> packptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
      const ConstUnmanaged<local_ordinal_type_1d_view> part2rowidx0;
      const ConstUnmanaged<local_ordinal_type_1d_view> lclrow;
      const local_ordinal_type blocksize;
      const local_ordinal_type num_vectors;

      // packed multivector output (or input)
      Unmanaged<vector_type_3d_view> packed_multivector;
      Unmanaged<impl_scalar_type_2d_view> scalar_multivector;
      impl_scalar_type damping_factor;

      template<typename TagType>
      KOKKOS_INLINE_FUNCTION
      void copy_multivectors(const local_ordinal_type &j, 
                             const local_ordinal_type &vi, 
                             const local_ordinal_type &pri, 
                             const local_ordinal_type &/* nrow */,  
                             const local_ordinal_type &ri0) const {
        if (TagType::id == 0) { // ToPackedMultiVectorTag
          for (local_ordinal_type col=0;col<num_vectors;++col) 
            for (local_ordinal_type i=0;i<blocksize;++i)
              packed_multivector(pri, i, col)[vi] = scalar_multivector(blocksize*lclrow(ri0+j)+i,col);
        } else if (TagType::id > 0) { //ToScalarMultiVector
          const impl_scalar_type df = damping_factor;
          for (local_ordinal_type col=0;col<num_vectors;++col) 
            for (local_ordinal_type i=0;i<blocksize;++i) {
              impl_scalar_type &y = scalar_multivector(blocksize*lclrow(ri0+j)+i,col);
              const impl_scalar_type yc = packed_multivector(pri, i, col)[vi];
              if (TagType::id == 1) y  = df*yc;
              else                  y += df*(yc - y); 
            }
        }
      }

      template<typename TagType>
      KOKKOS_INLINE_FUNCTION
      void copy_multivectors_with_norm(const local_ordinal_type &j, 
                                       const local_ordinal_type &vi, 
                                       const local_ordinal_type &pri, 
                                       const local_ordinal_type &/* nrow */,  
                                       const local_ordinal_type &ri0,
                                       /* */ magnitude_type *norm) const {
        if (TagType::id > 0) { //ToScalarMultiVector
          const impl_scalar_type df = damping_factor;
          for (local_ordinal_type col=0;col<num_vectors;++col) {
            const local_ordinal_type offset = col*blocksize;
            for (local_ordinal_type i=0;i<blocksize;++i) {
              impl_scalar_type &y = scalar_multivector(blocksize*lclrow(ri0+j)+i,col);
              const impl_scalar_type yc = packed_multivector(pri, i, col)[vi];
              const impl_scalar_type yd = TagType::id == 1 ? yc : yc - y;
              const magnitude_type abs_yd = Kokkos::ArithTraits<impl_scalar_type>::abs(yd);
              norm[offset+i] += abs_yd*abs_yd;
              if (TagType::id == 1) y  = df*yc;
              else                  y += df*yd; 
            }
          }
        }
      }

    public:

      // local reduction of norms with runtime array 
      // this value type and value_count is required for Kokkos
      typedef magnitude_type value_type[];
      int value_count;

      KOKKOS_INLINE_FUNCTION
      void 
      init(magnitude_type *dst) const {
        for (int i=0;i<value_count;++i)         
          dst[i] = Kokkos::reduction_identity<magnitude_type>::sum();
      }

      KOKKOS_INLINE_FUNCTION
      void 
      join(volatile magnitude_type *dst, const volatile magnitude_type *src) const {
        for (int i=0;i<value_count;++i) 
          dst[i] += src[i];
      }          

      MultiVectorConverter() = delete;
      MultiVectorConverter(const MultiVectorConverter &b) = default;
      MultiVectorConverter(const PartInterface<MatrixType> &interf,                    
                           const vector_type_3d_view &pmv) 
        : partptr(interf.partptr),
          packptr(interf.packptr),
          part2packrowidx0(interf.part2packrowidx0),
          part2rowidx0(interf.part2rowidx0),
          lclrow(interf.lclrow),
          blocksize(pmv.extent(1)),
          num_vectors(pmv.extent(2)),
          packed_multivector(pmv) {}

      // TODO:: modify this routine similar to the team level functions
      template<typename TagType>
      inline
      void 
      operator() (const TagType&, const local_ordinal_type &packidx, magnitude_type* const norm = NULL) const {
        local_ordinal_type partidx = packptr(packidx);
        local_ordinal_type npacks = packptr(packidx+1) - partidx;
        const local_ordinal_type pri0 = part2packrowidx0(partidx);

        local_ordinal_type ri0[vector_length] = {};
        local_ordinal_type nrows[vector_length] = {};
        for (local_ordinal_type v=0;v<npacks;++v,++partidx) {
          ri0[v] = part2rowidx0(partidx);
          nrows[v] = part2rowidx0(partidx+1) - ri0[v];
        }
        for (local_ordinal_type j=0;j<nrows[0];++j) {          
          local_ordinal_type cnt = 1;
          for (;cnt<npacks && j!= nrows[cnt];++cnt);
          npacks = cnt;
          const local_ordinal_type pri = pri0 + j;
          for (local_ordinal_type v=0;v<npacks;++v) {
            if (norm == NULL) copy_multivectors<TagType>(j, v, pri, nrows[v], ri0[v]);
            else              copy_multivectors_with_norm<TagType>(j, v, pri, nrows[v], ri0[v], norm);
          }
        }
      }

      template<typename TagType>
      KOKKOS_INLINE_FUNCTION
      void 
      operator() (const TagType&, const member_type &member, magnitude_type* const norm = NULL) const {
        const local_ordinal_type packidx = member.league_rank();
        const local_ordinal_type partidx_begin = packptr(packidx);
        const local_ordinal_type npacks = packptr(packidx+1) - partidx_begin;
        const local_ordinal_type pri0 = part2packrowidx0(partidx_begin);
	Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, npacks), [&](const local_ordinal_type &v) {
	    const local_ordinal_type partidx = partidx_begin + v;
	    const local_ordinal_type ri0 = part2rowidx0(partidx);
	    const local_ordinal_type nrows = part2rowidx0(partidx+1) - ri0;
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, nrows), [&](const local_ordinal_type &j) {
		const local_ordinal_type pri = pri0 + j;
                if (norm == NULL) copy_multivectors<TagType>(j, v, pri, nrows, ri0);
                else              copy_multivectors_with_norm<TagType>(j, v, pri, nrows, ri0, norm);
              });
          });
      }
      
      struct ToPackedMultiVectorTag       { enum : int { id = 0 }; };
      struct ToScalarMultiVectorFirstTag  { enum : int { id = 1 }; };
      struct ToScalarMultiVectorSecondTag { enum : int { id = 2 }; };
      
      template<typename TpetraLocalViewType>
      void to_packed_multivector(const TpetraLocalViewType &scalar_multivector_) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::MultiVectorConverter::ToPackedMultiVector");
#endif
        value_count = 0;
        scalar_multivector = scalar_multivector_;
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
	  const local_ordinal_type vl = vector_length;
          const Kokkos::TeamPolicy<execution_space,ToPackedMultiVectorTag> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl);
          Kokkos::parallel_for
            ("MultiVectorConverter::TeamPolicy::to_packed_multivector",
             policy, *this);
#endif
        } else {
#if defined(__CUDA_ARCH__)
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else
          const Kokkos::RangePolicy<execution_space,ToPackedMultiVectorTag> policy(0, packptr.extent(0) - 1);
          Kokkos::parallel_for
            ("MultiVectorConverter::RangePolicy::to_packed_multivector",
             policy, *this);
#endif
        }
      }
      
      template<typename TpetraLocalViewType>
      void to_scalar_multivector(const TpetraLocalViewType &scalar_multivector_, 
                                 const impl_scalar_type &damping_factor_,
                                 const bool &is_vectors_zero,
                                 /* */ magnitude_type *norm = NULL) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::MultiVectorConverter::ToScalarMultiVector");
#endif
        scalar_multivector = scalar_multivector_;
        damping_factor = damping_factor_;
        if (norm != NULL) {
          value_count = blocksize*num_vectors;
          for (int i=0;i<value_count;++i)
          norm[i] = Kokkos::ArithTraits<magnitude_type>::zero();
        } else {
          value_count = 0;
        }
        
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          // dynamic reduce does not support for vl > 1 
	  const local_ordinal_type vl = norm == NULL ? vector_length : 1;

          if (is_vectors_zero) {
            const Kokkos::TeamPolicy
              <execution_space,ToScalarMultiVectorFirstTag> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl);
            if (norm == NULL)  
              Kokkos::parallel_for
                ("MultiVectorConverter::TeamPolicy::to_scalar_multivector::fist", policy, *this);
            else               
              Kokkos::parallel_reduce
                ("MultiVectorConverter::TeamPolicy::to_scalar_multivector::fist_w_norm", policy, *this, norm);
          } else {
            const Kokkos::TeamPolicy
              <execution_space,ToScalarMultiVectorSecondTag> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl);
            if (norm == NULL)  
              Kokkos::parallel_for
                ("MultiVectorConverter::TeamPolicy::to_scalar_multivector::second", policy, *this);
            else               
              Kokkos::parallel_reduce
                ("MultiVectorConverter::TeamPolicy::to_scalar_multivector::second_w_norm", policy, *this, norm);              
          }
#endif
        } else {
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else
          if (is_vectors_zero) {
            const Kokkos::RangePolicy
              <execution_space,ToScalarMultiVectorFirstTag> policy(0, packptr.extent(0) - 1);
            if (norm == NULL)  
              Kokkos::parallel_for
                ("MultiVectorConverter::RangePolicy::to_scalar_multivector::fist", policy, *this);
            else               
              Kokkos::parallel_reduce
                ("MultiVectorConverter::RangePolicy::to_scalar_multivector::fist_w_norm", policy, *this, norm);              
          } else {
            const Kokkos::RangePolicy
              <execution_space,ToScalarMultiVectorSecondTag> policy(0, packptr.extent(0) - 1);
            if (norm == NULL)  
              Kokkos::parallel_for
                ("MultiVectorConverter::RangePolicy::to_scalar_multivector::second", policy, *this);
            else               
              Kokkos::parallel_reduce
                ("MultiVectorConverter::RangePolicy::to_scalar_multivector::second_w_norm", policy, *this, norm);              
          }
#endif
        }
      }
    };

    ///
    /// solve tridiags
    ///
    template<typename MatrixType>
    struct SolveTridiags {
    public:
      using impl_type = ImplType<MatrixType>;
      using execution_space = typename impl_type::execution_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
      using magnitude_type = typename impl_type::magnitude_type;
      /// views
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view; 
      /// vectorization 
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using impl_scalar_type_4d_view = typename impl_type::impl_scalar_type_4d_view;
      static constexpr int vector_length = impl_type::vector_length;

      /// team policy member type (used in cuda)
      using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;      
 
    private:
      // part interface
      const ConstUnmanaged<local_ordinal_type_1d_view> partptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> packptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
      // block tridiags 
      const ConstUnmanaged<size_type_1d_view> pack_td_ptr;
      // block tridiags values
      const ConstUnmanaged<vector_type_3d_view> D_vector_values;
      const Unmanaged<vector_type_3d_view> X_vector_values;
      // scalar tridiag values
      const ConstUnmanaged<impl_scalar_type_4d_view> D_scalar_values;
      const Unmanaged<impl_scalar_type_4d_view> X_scalar_values;
      
      const local_ordinal_type vector_length_value;

    public:
      SolveTridiags(const PartInterface<MatrixType> &interf,                    
                    const BlockTridiags<MatrixType> &btdm, 
                    const vector_type_3d_view &pmv) 
        :
        // interface
        partptr(interf.partptr), 
        packptr(interf.packptr),
        part2packrowidx0(interf.part2packrowidx0),
        // block tridiags and  multivector
        pack_td_ptr(btdm.pack_td_ptr), 
        D_vector_values(btdm.values),
        X_vector_values(pmv),
        // scalar tridiags and  multivector
        D_scalar_values((impl_scalar_type*)D_vector_values.data(),
                        D_vector_values.extent(0), 
                        D_vector_values.extent(1), 
                        D_vector_values.extent(2), 
                        vector_length),
        X_scalar_values((impl_scalar_type*)X_vector_values.data(),
                        X_vector_values.extent(0), 
                        X_vector_values.extent(1), 
                        X_vector_values.extent(2), 
                        vector_length),
	vector_length_value(vector_length)
      {}

    public:

      ///
      /// host serial (vector intrinsic) vectorization
      ///
      inline
      void 
      serialSolveSingleVector(const local_ordinal_type &blocksize, 
                              const local_ordinal_type &packidx) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level2::Unblocked;
        
        // base pointers
        auto A = D_vector_values.data();
        auto X = X_vector_values.data();

        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        // const local_ordinal_type blocksize = D_vector_values.extent(1);
        const local_ordinal_type astep = D_vector_values.stride_0();
        const local_ordinal_type as0 = blocksize; //D_vector_values.stride_1();
        const local_ordinal_type as1 = 1; //D_vector_values.stride_2();
        const local_ordinal_type xstep = X_vector_values.stride_0();
        const local_ordinal_type xs0 = 1; //X_vector_values.stride_1();

        // index counting
        const local_ordinal_type partidx = packptr(packidx);
        const size_type i0 = pack_td_ptr(partidx);
        const local_ordinal_type r0 = part2packrowidx0(partidx);
        const local_ordinal_type nrows = part2packrowidx0(packptr(packidx+1)) - r0;

        // move to starting point
        A += i0*astep; 
        X += r0*xstep;

        // solve Lx = x
        KOKKOSBATCHED_SERIAL_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE
          (AlgoType,KB::Diag::Unit,
           blocksize,blocksize,
           one, 
           A, as0, as1,
           X, xs0);

        for (local_ordinal_type tr=1;tr<nrows;++tr) {
          KOKKOSBATCHED_SERIAL_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,
             blocksize, blocksize,
             -one,
             A+2*astep, as0, as1,
             X, xs0,
             one,
             X+1*xstep, xs0);

          KOKKOSBATCHED_SERIAL_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,KB::Diag::Unit,
             blocksize,blocksize,
             one, 
             A+3*astep, as0, as1,
             X+1*xstep, xs0);

          A += 3*astep;
          X += 1*xstep;
        }
        
        // solve Ux = x
        KOKKOSBATCHED_SERIAL_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE
          (AlgoType,KB::Diag::NonUnit,
           blocksize, blocksize,
           one, 
           A, as0, as1,
           X, xs0);

        for (local_ordinal_type tr=nrows;tr>1;--tr) {
          A -= 3*astep;
          KOKKOSBATCHED_SERIAL_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,
             blocksize, blocksize,
             -one,
             A+1*astep, as0, as1,
             X, xs0,
             one,
             X-1*xstep, xs0);

          KOKKOSBATCHED_SERIAL_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,KB::Diag::NonUnit,
             blocksize, blocksize,
             one, 
             A, as0, as1,
             X-1*xstep,xs0);

          X -= 1*xstep;
        }
      }

      inline
      void 
      serialSolveMultiVector(const local_ordinal_type &/* blocksize */, 
                             const local_ordinal_type &packidx) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level3::Blocked;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        //const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);

        // subview pattern
        auto A = Kokkos::subview(D_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        auto X1 = Kokkos::subview(X_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        auto X2 = X1;

        // index counting
        const local_ordinal_type partidx = packptr(packidx);
        size_type i = pack_td_ptr(partidx);
        local_ordinal_type r = part2packrowidx0(partidx);
        const local_ordinal_type nrows = part2packrowidx0(packptr(packidx+1)) - r;

        // solve Lx = x
        A.assign_data( &D_vector_values(i,0,0) );
        X1.assign_data( &X_vector_values(r,0,0) );
        KB::SerialTrsm<KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>
          ::invoke(one, A, X1);
        for (local_ordinal_type tr=1;tr<nrows;++tr,i+=3) {
          A.assign_data( &D_vector_values(i+2,0,0) );
          X2.assign_data( &X_vector_values(++r,0,0) );
          KB::SerialGemm<KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
            ::invoke(-one, A, X1, one, X2);
          A.assign_data( &D_vector_values(i+3,0,0) );
          KB::SerialTrsm<KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>          
            ::invoke(one, A, X2);
          X1.assign_data( X2.data() );
        }
        
        // solve Ux = x
        KB::SerialTrsm<KB::Side::Left,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>          
          ::invoke(one, A, X1);
        for (local_ordinal_type tr=nrows;tr>1;--tr) {
          i -= 3;
          A.assign_data( &D_vector_values(i+1,0,0) );
          X2.assign_data( &X_vector_values(--r,0,0) );          
          KB::SerialGemm<KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
            ::invoke(-one, A, X1, one, X2); 
          A.assign_data( &D_vector_values(i,0,0) );          
          KB::SerialTrsm<KB::Side::Left,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>          
            ::invoke(one, A, X2);
          X1.assign_data( X2.data() );
        }
      }

      ///
      /// cuda team vectorization
      ///
      KOKKOS_INLINE_FUNCTION 
      void 
      teamSolveSingleVector(const member_type &member, 
                            const local_ordinal_type &blocksize,
                            const local_ordinal_type &i0,
                            const local_ordinal_type &r0,
                            const local_ordinal_type &nrows,
                            const local_ordinal_type &v) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level2::Unblocked;

        // base pointers
        auto A = D_scalar_values.data();
        auto X = X_scalar_values.data();

        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        //const local_ordinal_type num_vectors = X_scalar_values.extent(2);

        // const local_ordinal_type blocksize = D_scalar_values.extent(1);
        const local_ordinal_type astep = D_scalar_values.stride_0();
        const local_ordinal_type as0 = blocksize*vector_length; //D_scalar_values.stride_1();
        const local_ordinal_type as1 = vector_length; //D_scalar_values.stride_2();
        const local_ordinal_type xstep = X_scalar_values.stride_0();
        const local_ordinal_type xs0 = vector_length; //X_scalar_values.stride_1();

        // for multiple rhs
        //const local_ordinal_type xs0 = num_vectors*vector_length; //X_scalar_values.stride_1();
        //const local_ordinal_type xs1 = vector_length; //X_scalar_values.stride_2();

        // move to starting point
        A += i0*astep + v;
        X += r0*xstep + v;
        
        //for (local_ordinal_type col=0;col<num_vectors;++col) 
        {
          // solve Lx = x
          KOKKOSBATCHED_TEAM_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,member,KB::Diag::Unit,
             blocksize,blocksize,
             one, 
             A, as0, as1,
             X, xs0);

          for (local_ordinal_type tr=1;tr<nrows;++tr) {
            KOKKOSBATCHED_TEAM_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE
              (AlgoType,member,
               blocksize, blocksize,
               -one,
               A+2*astep, as0, as1,
               X, xs0,
               one,
               X+1*xstep, xs0);

            KOKKOSBATCHED_TEAM_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE
              (AlgoType,member,KB::Diag::Unit,
               blocksize,blocksize,
               one, 
               A+3*astep, as0, as1,
               X+1*xstep, xs0);
            
            A += 3*astep;
            X += 1*xstep;
          }
          
          // solve Ux = x
          KOKKOSBATCHED_TEAM_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE
            (AlgoType,member,KB::Diag::NonUnit,
             blocksize, blocksize,
             one, 
             A, as0, as1,
             X, xs0);

          for (local_ordinal_type tr=nrows;tr>1;--tr) {
            A -= 3*astep;
            KOKKOSBATCHED_TEAM_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE
              (AlgoType,member,
               blocksize, blocksize,
               -one,
               A+1*astep, as0, as1,
               X, xs0,
               one,
               X-1*xstep, xs0);
            
            KOKKOSBATCHED_TEAM_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE
              (AlgoType,member,KB::Diag::NonUnit,
               blocksize, blocksize,
               one, 
               A, as0, as1,
               X-1*xstep,xs0);
            
            X -= 1*xstep;
          }
          // for multiple rhs
          //X += xs1;
        }
      }

      KOKKOS_INLINE_FUNCTION 
      void 
      teamSolveMultiVector(const member_type &member, 
                           const local_ordinal_type &blocksize,
                           const local_ordinal_type &i0,
                           const local_ordinal_type &r0,
                           const local_ordinal_type &nrows,
                           const local_ordinal_type &v) const {
        namespace KB = KokkosBatched::Experimental;
        using AlgoType = KB::Algo::Level3::Blocked;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();

        // subview pattern
        auto A = Kokkos::subview(D_scalar_values, 0, Kokkos::ALL(), Kokkos::ALL(), 0);
        auto X1 = Kokkos::subview(X_scalar_values, 0, Kokkos::ALL(), Kokkos::ALL(), 0);
        auto X2 = X1;

        local_ordinal_type i = i0, r = r0;

        // solve Lx = x
        A.assign_data( &D_scalar_values(i,0,0,v) );
        X1.assign_data( &X_scalar_values(r,0,0,v) );
        KB::TeamTrsm<member_type,KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>
          ::invoke(member, one, A, X1);
        for (local_ordinal_type tr=1;tr<nrows;++tr,i+=3) {
          A.assign_data( &D_scalar_values(i+2,0,0,v) );
          X2.assign_data( &X_scalar_values(++r,0,0,v) );
          KB::TeamGemm<member_type,KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
            ::invoke(member, -one, A, X1, one, X2);
          A.assign_data( &D_scalar_values(i+3,0,0,v) );
          KB::TeamTrsm<member_type,KB::Side::Left,KB::Uplo::Lower,KB::Trans::NoTranspose,KB::Diag::Unit,AlgoType>          
            ::invoke(member, one, A, X2);
          X1.assign_data( X2.data() );
        }
        
        // solve Ux = x
        KB::TeamTrsm<member_type,KB::Side::Left,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>          
          ::invoke(member, one, A, X1);
        for (local_ordinal_type tr=nrows;tr>1;--tr) {
          i -= 3;
          A.assign_data( &D_scalar_values(i+1,0,0,v) );
          X2.assign_data( &X_scalar_values(--r,0,0,v) );          
          KB::TeamGemm<member_type,KB::Trans::NoTranspose,KB::Trans::NoTranspose,AlgoType>
            ::invoke(member, -one, A, X1, one, X2); 
          A.assign_data( &D_scalar_values(i,0,0,v) );          
          KB::TeamTrsm<member_type,KB::Side::Left,KB::Uplo::Upper,KB::Trans::NoTranspose,KB::Diag::NonUnit,AlgoType>          
            ::invoke(member, one, A, X2);
          X1.assign_data( X2.data() );
        }
      }


      template<int B> struct SingleVectorTag {};
      template<int B> struct MultiVectorTag {};
      
      template<int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const SingleVectorTag<B> &, const local_ordinal_type& packidx) const {
        const local_ordinal_type blocksize = (B == 0 ? D_vector_values.extent(1) : B);
        serialSolveSingleVector(blocksize, packidx);
      }      

      template<int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const MultiVectorTag<B> &, const local_ordinal_type& packidx) const {
        // this needs level 3 operations internal access (not yet)
        const local_ordinal_type blocksize = (B == 0 ? D_vector_values.extent(1) : B);
        serialSolveMultiVector(blocksize, packidx);
      }      

      template<int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const SingleVectorTag<B> &, const member_type &member) const {
	const local_ordinal_type packidx = member.league_rank();
	const local_ordinal_type partidx = packptr(packidx);	
	const local_ordinal_type i0 = pack_td_ptr(partidx);
	const local_ordinal_type r0 = part2packrowidx0(partidx);
        const local_ordinal_type nrows = partptr(partidx+1) - partptr(partidx);
        const local_ordinal_type blocksize = (B == 0 ? D_scalar_values.extent(1) : B);      
	Kokkos::parallel_for
	  (Kokkos::ThreadVectorRange(member, vector_length_value),[&](const int &v) {
	    teamSolveSingleVector(member, blocksize, i0, r0, nrows, v);
	  });
      }      

      template<int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const MultiVectorTag<B> &, const member_type &member) const {
	const local_ordinal_type packidx = member.league_rank();
	const local_ordinal_type partidx = packptr(packidx);	
	const local_ordinal_type i0 = pack_td_ptr(partidx);
	const local_ordinal_type r0 = part2packrowidx0(partidx);
        const local_ordinal_type nrows = partptr(partidx+1) - partptr(partidx);
        const local_ordinal_type blocksize = (B == 0 ? D_scalar_values.extent(1) : B);
	Kokkos::parallel_for
	  (Kokkos::ThreadVectorRange(member, vector_length_value),[&](const int &v) {
	    teamSolveMultiVector(member, blocksize, i0, r0, nrows, v);
	  });
      }      

      void run() {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStart();
#endif

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::SolveTridiags::Run");
#endif   
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          const local_ordinal_type num_vectors = X_scalar_values.extent(2);
	  const local_ordinal_type vl = vector_length;
#define BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(B)    \
          if (num_vectors == 1) {                                       \
            const Kokkos::TeamPolicy<execution_space,SingleVectorTag<B> > policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl); \
            Kokkos::parallel_for                                        \
              ("SolveTridiags::TeamPolicy::run<SingleVector>", policy, *this); \
          } else {                                                      \
            const Kokkos::TeamPolicy<execution_space,MultiVectorTag<B> > policy(packptr.extent(0) - 1, Kokkos::AUTO(), vl); \
            Kokkos::parallel_for                                        \
              ("SolveTridiags::TeamPolicy::run<MultiVector>", policy, *this); \
          } break

          const local_ordinal_type blocksize = D_scalar_values.extent(1);
          switch (blocksize) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(18);
          default : BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS

#endif
        } else {
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else
          const local_ordinal_type num_vectors = X_vector_values.extent(2);
#define BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(B)    \
          if (num_vectors == 1) {                                       \
            const Kokkos::RangePolicy<execution_space,SingleVectorTag<B> > policy(0, packptr.extent(0) - 1); \
            Kokkos::parallel_for                                        \
              ("SolveTridiags::RangePolicy::run<SingleVector>", policy, *this); \
          } else {                                                      \
            const Kokkos::RangePolicy<execution_space,MultiVectorTag<B> > policy(0, packptr.extent(0) - 1); \
            Kokkos::parallel_for                                        \
              ("SolveTridiags::RangePolicy::run<MultiVector>", policy, *this); \
          } break
          
          const local_ordinal_type blocksize = D_vector_values.extent(1);
          switch (blocksize) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS(18);
          default : BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_SOLVETRIDIAGS
#endif
        }

#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStop();
#endif
      }
    }; 
    
    ///
    /// compute local residula vector y = b - R x 
    ///
    template<typename MatrixType>
    struct ComputeResidualVector {
    public:
      using impl_type = ImplType<MatrixType>;
      using execution_space = typename impl_type::execution_space;
      using memory_space = typename impl_type::memory_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
      using magnitude_type = typename impl_type::magnitude_type;
      /// views
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view; 
      using tpetra_block_access_view_type = typename impl_type::tpetra_block_access_view_type; // block crs (layout right)
      using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view; 
      using impl_scalar_type_2d_view = typename impl_type::impl_scalar_type_2d_view; // block multivector (layout left)
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using impl_scalar_type_4d_view = typename impl_type::impl_scalar_type_4d_view;
      static constexpr int vector_length = impl_type::vector_length;

      /// team policy member type (used in cuda)
      using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;      

      // enum for max blocksize and vector length
      enum : int { max_blocksize = 32 };

    private:
      ConstUnmanaged<impl_scalar_type_2d_view> b;
      ConstUnmanaged<impl_scalar_type_2d_view> x; // x_owned
      ConstUnmanaged<impl_scalar_type_2d_view> x_remote;
      /* */Unmanaged<impl_scalar_type_2d_view> y;
      /* */Unmanaged<vector_type_3d_view> y_packed;
      /* */Unmanaged<impl_scalar_type_4d_view> y_packed_scalar;

      // AmD information
      const ConstUnmanaged<size_type_1d_view> rowptr, rowptr_remote;
      const ConstUnmanaged<local_ordinal_type_1d_view> colindsub, colindsub_remote;
      const ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

      // block crs graph information
      // for cuda (kokkos crs graph uses a different size_type from size_t)
      using a_rowptr_value_type = typename Kokkos::ViewTraits<local_ordinal_type*,typename impl_type::device_type>::size_type; 
      const ConstUnmanaged<Kokkos::View<a_rowptr_value_type*,typename impl_type::device_type> > A_rowptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> A_colind;

      // blocksize
      const local_ordinal_type blocksize_requested;

      // part interface
      const ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
      const ConstUnmanaged<local_ordinal_type_1d_view> part2rowidx0;
      const ConstUnmanaged<local_ordinal_type_1d_view> rowidx2part;
      const ConstUnmanaged<local_ordinal_type_1d_view> partptr;
      const ConstUnmanaged<local_ordinal_type_1d_view> lclrow;     
      const ConstUnmanaged<local_ordinal_type_1d_view> dm2cm;
      const bool is_dm2cm_active;

    public:
      template<typename LocalCrsGraphType>
      ComputeResidualVector(const AmD<MatrixType> &amd,
                            const LocalCrsGraphType &graph,
                            const local_ordinal_type &blocksize_requested_,
                            const PartInterface<MatrixType> &interf,
                            const local_ordinal_type_1d_view &dm2cm_) 
        : rowptr(amd.rowptr), rowptr_remote(amd.rowptr_remote),
          colindsub(amd.A_colindsub), colindsub_remote(amd.A_colindsub_remote),
          tpetra_values(amd.tpetra_values),
          A_rowptr(graph.row_map),
          A_colind(graph.entries),
          blocksize_requested(blocksize_requested_),
          part2packrowidx0(interf.part2packrowidx0),
          part2rowidx0(interf.part2rowidx0),
          rowidx2part(interf.rowidx2part),
          partptr(interf.partptr),
          lclrow(interf.lclrow),
          dm2cm(dm2cm_),
          is_dm2cm_active(dm2cm_.span() > 0)
      {}


      inline
      void 
      serialGemv(const local_ordinal_type &blocksize,
                 const impl_scalar_type * const __restrict__ AA, 
                 const impl_scalar_type * const __restrict__ xx,
                 /* */ impl_scalar_type * __restrict__ yy) const {
        for (local_ordinal_type k0=0;k0<blocksize;++k0) {
          impl_scalar_type val = 0;
          const local_ordinal_type offset = k0*blocksize;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#   pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#   pragma unroll
#endif
          for (local_ordinal_type k1=0;k1<blocksize;++k1) 
            val += AA[offset+k1]*xx[k1];
          yy[k0] -= val;
        }
      }

      template<typename bbViewType, typename yyViewType>
      KOKKOS_INLINE_FUNCTION       
      void 
      vectorCopy(const member_type &member,
                 const local_ordinal_type &blocksize,
                 const bbViewType &bb, 
                 const yyViewType &yy) const {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize), [&](const local_ordinal_type &k0)  {
            yy(k0) = bb(k0);
          });
      }

      template<typename AAViewType, typename xxViewType, typename yyViewType>
      KOKKOS_INLINE_FUNCTION       
      void 
      teamGemv(const member_type &member, 
               const local_ordinal_type &blocksize,
               const AAViewType &AA, 
               const xxViewType &xx, 
               const yyViewType &yy) const { 
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member, blocksize), 
           [&](const local_ordinal_type &k0) {
            impl_scalar_type val = 0;
            Kokkos::parallel_reduce
              (Kokkos::ThreadVectorRange(member, blocksize), 
               [&](const local_ordinal_type &k1, impl_scalar_type &update) {
                update += AA(k0,k1)*xx(k1);
              }, val);
            Kokkos::single(Kokkos::PerThread(member), [&]() {
                yy(k0) -= val;
              });
          });
      }

      struct SeqTag {};

      inline
      void 
      operator() (const SeqTag &, const local_ordinal_type& i) const {
        const local_ordinal_type blocksize = blocksize_requested;
        const local_ordinal_type blocksize_square = blocksize*blocksize;
        
        // constants
        const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);
        const local_ordinal_type num_vectors = y.extent(1);
        const local_ordinal_type row = i*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // y := b
          impl_scalar_type *yy = &y(row, col);
          const impl_scalar_type * const bb = &b(row, col);
          memcpy(yy, bb, sizeof(impl_scalar_type)*blocksize);         
        
          // y -= Rx
          const size_type A_k0 = A_rowptr[i];
          for (size_type k=rowptr[i];k<rowptr[i+1];++k) {
            const size_type j = A_k0 + colindsub[k];
            const impl_scalar_type * const AA = &tpetra_values(j*blocksize_square);
            const impl_scalar_type * const xx = &x(A_colind[j]*blocksize, col);
            serialGemv(blocksize,AA,xx,yy);
          }
        }
      }

      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const SeqTag &, const member_type &member) const {
        namespace KB = KokkosBatched::Experimental;
        
        // constants
        const local_ordinal_type blocksize = blocksize_requested;
        const local_ordinal_type blocksize_square = blocksize*blocksize;

        const local_ordinal_type i = member.league_rank();
        const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);
        const local_ordinal_type num_vectors = y.extent(1);

        // subview pattern
        auto bb = Kokkos::subview(b, block_range, 0);
        auto xx = bb;
        auto yy = Kokkos::subview(y, block_range, 0);
        auto A_block = ConstUnmanaged<tpetra_block_access_view_type>(NULL, blocksize, blocksize); 
        
        const local_ordinal_type row = i*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // y := b
          yy.assign_data(&y(row, col));
          bb.assign_data(&b(row, col));
	  if (member.team_rank() == 0) 
            vectorCopy(member, blocksize, bb, yy);
          member.team_barrier();

          // y -= Rx
          const size_type A_k0 = A_rowptr[i];              
          for (size_type k=rowptr[i];k<rowptr[i+1];++k) {
            const size_type j = A_k0 + colindsub[k];
            A_block.assign_data( &tpetra_values(j*blocksize_square) );
            xx.assign_data( &x(A_colind[j]*blocksize, col) );
            teamGemv(member, blocksize, A_block, xx, yy);
          }
        }
      }

      template<int B>
      struct AsyncTag {};

      template<int B>
      inline
      void 
      operator() (const AsyncTag<B> &, const local_ordinal_type &rowidx) const {
        const local_ordinal_type blocksize = (B == 0 ? blocksize_requested : B);
        const local_ordinal_type blocksize_square = blocksize*blocksize;

        // constants        
        const local_ordinal_type partidx = rowidx2part(rowidx);
        const local_ordinal_type pri = part2packrowidx0(partidx) + (rowidx - partptr(partidx));
        const local_ordinal_type v = partidx % vector_length;

        const local_ordinal_type num_vectors = y_packed.extent(2);
        const local_ordinal_type num_local_rows = lclrow.extent(0);

        // temporary buffer for y flat
        impl_scalar_type yy[B == 0 ? max_blocksize : B] = {};

        const local_ordinal_type lr = lclrow(rowidx);
        const local_ordinal_type row = lr*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // y := b
          memcpy(yy, &b(row, col), sizeof(impl_scalar_type)*blocksize);
          
          // y -= Rx
          const size_type A_k0 = A_rowptr[lr];
          for (size_type k=rowptr[lr];k<rowptr[lr+1];++k) {
            const size_type j = A_k0 + colindsub[k];
            const impl_scalar_type * const AA = &tpetra_values(j*blocksize_square);
            const local_ordinal_type A_colind_at_j = A_colind[j];
            if (A_colind_at_j < num_local_rows) {
              const auto loc = is_dm2cm_active ? dm2cm[A_colind_at_j] : A_colind_at_j;
              const impl_scalar_type * const xx = &x(loc*blocksize, col);
              serialGemv(blocksize, AA,xx,yy);
            } else {
              const auto loc = A_colind_at_j - num_local_rows;
              const impl_scalar_type * const xx_remote = &x_remote(loc*blocksize, col);
              serialGemv(blocksize, AA,xx_remote,yy);                
            }
          }
          // move yy to y_packed
          for (local_ordinal_type k=0;k<blocksize;++k) 
            y_packed(pri, k, col)[v] = yy[k];
        }
      }

      template<int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const AsyncTag<B> &, const member_type &member) const {
        const local_ordinal_type blocksize = (B == 0 ? blocksize_requested : B);
        const local_ordinal_type blocksize_square = blocksize*blocksize;

        // constants        
        const local_ordinal_type rowidx = member.league_rank();
        const local_ordinal_type partidx = rowidx2part(rowidx);
        const local_ordinal_type pri = part2packrowidx0(partidx) + (rowidx - partptr(partidx));
        const local_ordinal_type v = partidx % vector_length;

        const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);
        const local_ordinal_type num_vectors = y_packed_scalar.extent(2);
        const local_ordinal_type num_local_rows = lclrow.extent(0);

        // subview pattern
        auto bb = Kokkos::subview(b, block_range, 0);
        auto xx = Kokkos::subview(x, block_range, 0);
        auto xx_remote = Kokkos::subview(x_remote, block_range, 0);
        auto yy = Kokkos::subview(y_packed_scalar, 0, block_range, 0, 0);
        auto A_block = ConstUnmanaged<tpetra_block_access_view_type>(NULL, blocksize, blocksize); 

        const local_ordinal_type lr = lclrow(rowidx);
        const local_ordinal_type row = lr*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // y := b
          bb.assign_data(&b(row, col));
          yy.assign_data(&y_packed_scalar(pri, 0, col, v));
          if (member.team_rank() == 0) 
            vectorCopy(member, blocksize, bb, yy);
          member.team_barrier();

          // y -= Rx
          const size_type A_k0 = A_rowptr[lr];
          for (size_type k=rowptr[lr];k<rowptr[lr+1];++k) {
            const size_type j = A_k0 + colindsub[k];
            A_block.assign_data( &tpetra_values(j*blocksize_square) );
              
            const local_ordinal_type A_colind_at_j = A_colind[j];
            if (A_colind_at_j < num_local_rows) {
              const auto loc = is_dm2cm_active ? dm2cm[A_colind_at_j] : A_colind_at_j;
              xx.assign_data( &x(loc*blocksize, col) );
              teamGemv(member, blocksize, A_block, xx, yy);
            } else {
              const auto loc = A_colind_at_j - num_local_rows;
              xx_remote.assign_data( &x_remote(loc*blocksize, col) );
              teamGemv(member, blocksize, A_block, xx_remote, yy);
            }
          }
        }
      }
      
      template <int P, int B> struct OverlapTag {};

      template<int P, int B>
      inline
      void 
      operator() (const OverlapTag<P,B> &, const local_ordinal_type& rowidx) const {
        const local_ordinal_type blocksize = (B == 0 ? blocksize_requested : B);
        const local_ordinal_type blocksize_square = blocksize*blocksize;

        // constants        
        const local_ordinal_type partidx = rowidx2part(rowidx);
        const local_ordinal_type pri = part2packrowidx0(partidx) + (rowidx - partptr(partidx));
        const local_ordinal_type v = partidx % vector_length;

        const local_ordinal_type num_vectors = y_packed.extent(2);
        const local_ordinal_type num_local_rows = lclrow.extent(0);

        // temporary buffer for y flat
        impl_scalar_type yy[max_blocksize] = {};

        auto colindsub_used = (P == 0 ? colindsub : colindsub_remote); 
        auto rowptr_used = (P == 0 ? rowptr : rowptr_remote);

        const local_ordinal_type lr = lclrow(rowidx);
        const local_ordinal_type row = lr*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          if (P == 0) { 
            // y := b
            memcpy(yy, &b(row, col), sizeof(impl_scalar_type)*blocksize);
          } else {
            // y (temporary) := 0
            memset(yy, 0, sizeof(impl_scalar_type)*blocksize); 
          }
          
          // y -= Rx
          const size_type A_k0 = A_rowptr[lr];
          for (size_type k=rowptr_used[lr];k<rowptr_used[lr+1];++k) {
            const size_type j = A_k0 + colindsub_used[k];
            const impl_scalar_type * const AA = &tpetra_values(j*blocksize_square);
            const local_ordinal_type A_colind_at_j = A_colind[j];
            if (P == 0) {
              const auto loc = is_dm2cm_active ? dm2cm[A_colind_at_j] : A_colind_at_j;
              const impl_scalar_type * const xx = &x(loc*blocksize, col);
              serialGemv(blocksize,AA,xx,yy);
            } else {
              const auto loc = A_colind_at_j - num_local_rows;
              const impl_scalar_type * const xx_remote = &x_remote(loc*blocksize, col);
              serialGemv(blocksize,AA,xx_remote,yy);                
            }
          }
          // move yy to y_packed
          if (P == 0) {
            for (local_ordinal_type k=0;k<blocksize;++k) 
              y_packed(pri, k, col)[v] = yy[k];
          } else {
            for (local_ordinal_type k=0;k<blocksize;++k) 
              y_packed(pri, k, col)[v] += yy[k];
          }
        }
      }

      template<int P, int B>
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const OverlapTag<P,B> &, const member_type &member) const {
        const local_ordinal_type blocksize = (B == 0 ? blocksize_requested : B);
        const local_ordinal_type blocksize_square = blocksize*blocksize;

        // constants        
        const local_ordinal_type rowidx = member.league_rank();
        const local_ordinal_type partidx = rowidx2part(rowidx);
        const local_ordinal_type pri = part2packrowidx0(partidx) + (rowidx - partptr(partidx));
        const local_ordinal_type v = partidx % vector_length;

        const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);
        const local_ordinal_type num_vectors = y_packed_scalar.extent(2);
        const local_ordinal_type num_local_rows = lclrow.extent(0);

        // subview pattern
        auto bb = Kokkos::subview(b, block_range, 0);
        auto xx = Kokkos::subview(x, block_range, 0);
        auto xx_remote = Kokkos::subview(x_remote, block_range, 0);
        auto yy = Kokkos::subview(y_packed_scalar, 0, block_range, 0, 0);
        auto A_block = ConstUnmanaged<tpetra_block_access_view_type>(NULL, blocksize, blocksize);
        auto colindsub_used = (P == 0 ? colindsub : colindsub_remote); 
        auto rowptr_used = (P == 0 ? rowptr : rowptr_remote);

        const local_ordinal_type lr = lclrow(rowidx);
        const local_ordinal_type row = lr*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          yy.assign_data(&y_packed_scalar(pri, 0, col, v));
          if (P == 0) { 
            // y := b
            bb.assign_data(&b(row, col));
            if (member.team_rank() == 0) 
              vectorCopy(member, blocksize, bb, yy);
            member.team_barrier();            
          } 
          
          // y -= Rx
          const size_type A_k0 = A_rowptr[lr];
          for (size_type k=rowptr_used[lr];k<rowptr_used[lr+1];++k) {
            const size_type j = A_k0 + colindsub_used[k];
            A_block.assign_data( &tpetra_values(j*blocksize_square) );

            const local_ordinal_type A_colind_at_j = A_colind[j];
            if (P == 0) {
              const auto loc = is_dm2cm_active ? dm2cm[A_colind_at_j] : A_colind_at_j;
              xx.assign_data( &x(loc*blocksize, col) );
              teamGemv(member, blocksize, A_block, xx, yy);
            } else {
              const auto loc = A_colind_at_j - num_local_rows;
              xx_remote.assign_data( &x_remote(loc*blocksize, col) );
              teamGemv(member, blocksize, A_block, xx_remote, yy);
            }
          }
        }
      }

      // y = b - Rx; seq method
      template<typename MultiVectorLocalViewTypeY,
               typename MultiVectorLocalViewTypeB,
               typename MultiVectorLocalViewTypeX>
      void run(const MultiVectorLocalViewTypeY &y_, 
               const MultiVectorLocalViewTypeB &b_, 
               const MultiVectorLocalViewTypeX &x_) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStart();
#endif

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::ComputeResidual::Run<SeqTag>");
#endif
        y = y_; b = b_; x = x_; 
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          const Kokkos::TeamPolicy<execution_space,SeqTag> policy(rowptr.extent(0) - 1, Kokkos::AUTO());
          Kokkos::parallel_for
            ("ComputeResidual::TeamPolicy::run<SeqTag>", policy, *this);
#endif
        } else {
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else          
          const Kokkos::RangePolicy<execution_space,SeqTag> policy(0, rowptr.extent(0) - 1);
          Kokkos::parallel_for
            ("ComputeResidual::RangePolicy::run<SeqTag>", policy, *this);
#endif
        }

#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStop();
#endif
      }

      // y = b - R (x , x_remote)
      template<typename MultiVectorLocalViewTypeB,
               typename MultiVectorLocalViewTypeX,
               typename MultiVectorLocalViewTypeX_Remote>
      void run(const vector_type_3d_view &y_packed_, 
               const MultiVectorLocalViewTypeB &b_, 
               const MultiVectorLocalViewTypeX &x_,
               const MultiVectorLocalViewTypeX_Remote &x_remote_) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStart();
#endif

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::ComputeResidual::Run<AsyncTag>");
#endif
        b = b_; x = x_; x_remote = x_remote_;
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          y_packed_scalar = impl_scalar_type_4d_view((impl_scalar_type*)y_packed_.data(),
                                                     y_packed_.extent(0),
                                                     y_packed_.extent(1),
                                                     y_packed_.extent(2),
                                                     vector_length);
#endif
        } else {
          y_packed = y_packed_;
        }

        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          local_ordinal_type vl_power_of_two = 1;
          for (;vl_power_of_two<=blocksize_requested;vl_power_of_two*=2);
          vl_power_of_two *= (vl_power_of_two < blocksize_requested ? 2 : 1);
	  const local_ordinal_type vl = vl_power_of_two > vector_length ? vector_length : vl_power_of_two;
#define BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(B) {                \
            const Kokkos::TeamPolicy<execution_space,AsyncTag<B> > policy(rowidx2part.extent(0), Kokkos::AUTO(), vl); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::TeamPolicy::run<AsyncTag>",            \
               policy, *this); } break
          switch (blocksize_requested) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(18);
          default : BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL
#endif
        } else {          
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else          
#define BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(B) {                \
            const Kokkos::RangePolicy<execution_space,AsyncTag<B> > policy(0, rowidx2part.extent(0)); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::RangePolicy::run<AsyncTag>",           \
               policy, *this); } break
          switch (blocksize_requested) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(18);
          default : BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL
#endif
        }
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStop();
#endif
      }
      
      // y = b - R (y , y_remote)
      template<typename MultiVectorLocalViewTypeB,
               typename MultiVectorLocalViewTypeX,
               typename MultiVectorLocalViewTypeX_Remote>
      void run(const vector_type_3d_view &y_packed_, 
               const MultiVectorLocalViewTypeB &b_, 
               const MultiVectorLocalViewTypeX &x_,
               const MultiVectorLocalViewTypeX_Remote &x_remote_,
               const bool compute_owned) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStart();
#endif

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::ComputeResidual::Run<OverlapTag>");
#endif
        b = b_; x = x_; x_remote = x_remote_;
        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          y_packed_scalar = impl_scalar_type_4d_view((impl_scalar_type*)y_packed_.data(),
                                                     y_packed_.extent(0),
                                                     y_packed_.extent(1),
                                                     y_packed_.extent(2),
                                                     vector_length);
#endif
        } else {
          y_packed = y_packed_; 
        }

        if (is_cuda<execution_space>::value) {
#if defined(KOKKOS_ENABLE_CUDA)
          local_ordinal_type vl_power_of_two = 1;
          for (;vl_power_of_two<=blocksize_requested;vl_power_of_two*=2);
          vl_power_of_two *= (vl_power_of_two < blocksize_requested ? 2 : 1);
	  const local_ordinal_type vl = vl_power_of_two > vector_length ? vector_length : vl_power_of_two;
#define BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(B)  \
          if (compute_owned) {                                          \
            const Kokkos::TeamPolicy<execution_space,OverlapTag<0,B> > policy(rowidx2part.extent(0), Kokkos::AUTO(), vl); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::TeamPolicy::run<OverlapTag<0> >", policy, *this); \
          } else {                                                      \
            const Kokkos::TeamPolicy<execution_space,OverlapTag<1,B> > policy(rowidx2part.extent(0), Kokkos::AUTO(), vl); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::TeamPolicy::run<OverlapTag<1> >", policy, *this); \
          } break
          switch (blocksize_requested) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(18);
          default : BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL          
#endif
        } else {
#if defined(__CUDA_ARCH__)        
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: CUDA should not see this code"); 
#else          
#define BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(B)  \
          if (compute_owned) {                                          \
            const Kokkos::RangePolicy<execution_space,OverlapTag<0,B> > policy(0, rowidx2part.extent(0)); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::RangePolicy::run<OverlapTag<0> >", policy, *this); \
          } else {                                                      \
            const Kokkos::RangePolicy<execution_space,OverlapTag<1,B> > policy(0, rowidx2part.extent(0)); \
            Kokkos::parallel_for                                        \
              ("ComputeResidual::RangePolicy::run<OverlapTag<1> >", policy, *this); \
          } break

          switch (blocksize_requested) {
          case   3: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 3);
          case   5: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 5);
          case   7: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 7);
          case   9: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 9);
          case  10: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(10);
          case  11: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(11);
          case  16: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(16);
          case  17: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(17);
          case  18: BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL(18);
          default : BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL( 0);            
          }
#undef BLOCKTRIDICONTAINER_DETAILS_COMPUTERESIDUAL
#endif
        }
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
        cudaProfilerStop();
#endif
      }
    }; 

    ///
    /// Manage the distributed part of the computation of residual norms.
    ///
    template<typename MatrixType>
    struct NormManager {
    public: 
      using impl_type = ImplType<MatrixType>;
      using magnitude_type = typename impl_type::magnitude_type;

    private:
      bool collective_;
      int sweep_step_, blocksize_, num_vectors_;
#ifdef HAVE_IFPACK2_MPI
      MPI_Request mpi_request_;
      MPI_Comm comm_;
#endif
      std::vector<magnitude_type> work_;

    public:
      NormManager() = default;
      NormManager(const NormManager &b) = default;
      NormManager(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
        sweep_step_ = 1;
        blocksize_ = 0;
        num_vectors_ = 0;

        collective_ = comm->getSize() > 1;
        if (collective_) {
#ifdef HAVE_IFPACK2_MPI
          const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
          TEUCHOS_ASSERT( ! mpi_comm.is_null());
          comm_ = *mpi_comm->getRawMpiComm();
#endif
        }
      }

      int getBlocksize() const { return blocksize_; }
      int getNumVectors() const { return num_vectors_; }

      // Resize the buffer to accommodate nvec vectors in the multivector, for a
      // matrix having block size block_size.
      void resize(const int& blocksize, const int& num_vectors) {
        blocksize_ = blocksize;
        num_vectors_ = num_vectors;
        work_.resize((2*blocksize_ + 1)*num_vectors_);
      }
      
      // Check the norm every sweep_step sweeps.
      void setCheckFrequency(const int& sweep_step) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(sweep_step < 1, "sweep step must be >= 1");
        sweep_step_ = sweep_step;
      }
      
      // Get the buffer into which to store rank-local squared norms.
      magnitude_type* getBuffer() { return work_.data(); }

      // Call MPI_Iallreduce to find the global squared norms.
      void ireduce(const int& sweep, const bool force = false) {
        if ( ! force && sweep % sweep_step_) return;
        const int n = blocksize_*num_vectors_;
        if (collective_) {
          std::copy(work_.begin(), work_.begin() + n, work_.begin() + n);
#ifdef HAVE_IFPACK2_MPI
# if MPI_VERSION >= 3
          MPI_Iallreduce(work_.data() + n, work_.data(), n,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_, &mpi_request_);
# else
          MPI_Allreduce (work_.data() + n, work_.data(), n,
                         Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                         MPI_SUM, comm_);
# endif
#endif
        }
      }
      
      // Check if the norm-based termination criterion is met. tol2 is the
      // tolerance squared. Sweep is the sweep index. If not every iteration is
      // being checked, this function immediately returns false. If a check must
      // be done at this iteration, it waits for the reduction triggered by
      // ireduce to complete, then checks the global norm against the tolerance.
      bool checkDone (const int& sweep, const magnitude_type& tol2, const bool force = false) {
        // early return 
        if (sweep <= 0) return false;

#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NormManager::CheckDone");
#endif
        TEUCHOS_ASSERT(sweep >= 1);
        if ( ! force && (sweep - 1) % sweep_step_) return false;
        if (collective_) {
#ifdef HAVE_IFPACK2_MPI
# if MPI_VERSION >= 3
          MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
# else
          // Do nothing.
# endif
#endif
        }
        const auto n = blocksize_*num_vectors_;
        if (sweep == 1) {
          magnitude_type* const n0 = work_.data() + 2*n;
          for (int v = 0; v < num_vectors_; ++v) {
            const magnitude_type* const dn0 = work_.data() + v*blocksize_;
            magnitude_type mdn0 = 0;
            for (int i = 0; i < blocksize_; ++i)
              mdn0 = std::max(mdn0, dn0[i]);
            n0[v] = mdn0;
          }
          return false;
        } else {
          const auto n0 = work_.data() + 2*n;
          bool done = true;
          for (int v = 0; v < num_vectors_; ++v) {
            const magnitude_type* const dnf = work_.data() + v*blocksize_;
            magnitude_type mdnf = 0;
            for (int i = 0; i < blocksize_; ++i)
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
        for (int v = 0; v < num_vectors_; ++v) {
          const magnitude_type* const dnf = work_.data() + v*blocksize_;
          magnitude_type mdnf = 0;
          for (int i = 0; i < blocksize_; ++i)
            mdnf = std::max(mdnf, dnf[i]);
          // This overwrites the receive buffer, but that's OK; at the time of
          // this write, we no longer need the data in this slot.
          work_[v] = mdnf;
        }
        for (int i = 0; i < num_vectors_; ++i)
          work_[i] = std::sqrt(work_[i]);
        magnitude_type* const nf = work_.data() + 2*blocksize_*num_vectors_;
        for (int v = 0; v < num_vectors_; ++v)
          nf[v] = std::sqrt(nf[v]);
      }
      
      // Report norms to the caller.
      const magnitude_type* getNorms0 () const { return work_.data() + 2*blocksize_*num_vectors_; }
      const magnitude_type* getNormsFinal () const { return work_.data(); }
    };

    ///
    /// top level apply interface
    ///
    template<typename MatrixType>
    int 
    applyInverseJacobi(// importer
                       const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
                       const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_import_type> &tpetra_importer,
                       const Teuchos::RCP<AsyncableImport<MatrixType> > &async_importer,
                       const bool overlap_communication_and_computation,
                       // tpetra interface
                       const typename ImplType<MatrixType>::tpetra_multivector_type &X,  // tpetra interface
                       /* */ typename ImplType<MatrixType>::tpetra_multivector_type &Y,  // tpetra interface
                       /* */ typename ImplType<MatrixType>::tpetra_multivector_type &Z,  // temporary tpetra interface
                       // local object interface
                       const PartInterface<MatrixType> &interf, // mesh interface
                       const BlockTridiags<MatrixType> &btdm, // packed block tridiagonal matrices
                       const AmD<MatrixType> &amd, // R = A - D
                       /* */ typename ImplType<MatrixType>::vector_type_1d_view &work, // workspace for packed multivector of right hand side 
                       /* */ NormManager<MatrixType> &norm_manager,
                       // preconditioner parameters
                       const typename ImplType<MatrixType>::impl_scalar_type &damping_factor, 
                       /* */ bool is_y_zero,
                       const int max_num_sweeps, 
                       const typename ImplType<MatrixType>::magnitude_type tol,
                       const int check_tol_every) { 
      
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::ApplyInverseJacobi");
#endif
      using impl_type = ImplType<MatrixType>;
      using memory_space = typename impl_type::memory_space;

      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using magnitude_type = typename impl_type::magnitude_type;
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using vector_type_1d_view = typename impl_type::vector_type_1d_view;
      using vector_type_3d_view = typename impl_type::vector_type_3d_view;
      using tpetra_multivector_type = typename impl_type::tpetra_multivector_type;

      // either tpetra importer or async importer must be active
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!tpetra_importer.is_null() && !async_importer.is_null(), 
                                  "Neither Tpetra importer nor Async importer is null.");
      // max number of sweeps should be positive number
      TEUCHOS_TEST_FOR_EXCEPT_MSG(max_num_sweeps <= 0, 
                                  "Maximum number of sweeps must be >= 1.");

      // const parameters
      const bool is_norm_manager_active = tol > Kokkos::ArithTraits<magnitude_type>::zero();
      const bool is_seq_method_requested = !tpetra_importer.is_null();
      const bool is_async_importer_active = !async_importer.is_null();
      const magnitude_type tolerance = tol*tol;
      const local_ordinal_type blocksize = btdm.values.extent(1);
      const local_ordinal_type num_vectors = Y.getNumVectors();
      const local_ordinal_type num_blockrows = interf.part2packrowidx0_back;

      TEUCHOS_TEST_FOR_EXCEPT_MSG(is_norm_manager_active && is_seq_method_requested,
                                  "The seq method for applyInverseJacobi, " << 
                                  "which in any case is for developer use only, " <<
                                  "does not support norm-based termination.");

      // if workspace is needed more, resize it
      const size_type work_span_required = num_blockrows*num_vectors*blocksize;
      if (work.span() < work_span_required) 
        work = vector_type_1d_view("vector workspace 1d view", work_span_required);
     
      typename AsyncableImport<MatrixType>::impl_scalar_type_2d_view remote_multivector;
      if (is_seq_method_requested) {
        // construct copy of Y again if num vectors are different
        if (static_cast<local_ordinal_type>(Z.getNumVectors()) != num_vectors) 
          Z = tpetra_multivector_type(tpetra_importer->getTargetMap(), num_vectors, false);
      } else {
        if (is_async_importer_active) {
          // create comm data buffer and keep it here
          async_importer->createDataBuffer(num_vectors);
          remote_multivector = async_importer->getRemoteMultiVectorLocalView();
        }
      }

      // wrap the workspace with 3d view
      vector_type_3d_view pmv(work.data(), num_blockrows, blocksize, num_vectors);
      const auto XX = X.template getLocalView<memory_space>();
      const auto YY = Y.template getLocalView<memory_space>();
      const auto ZZ = Z.template getLocalView<memory_space>();

      MultiVectorConverter<MatrixType> multivector_converter(interf, pmv);
      SolveTridiags<MatrixType> solve_tridiags(interf, btdm, pmv);

      const local_ordinal_type_1d_view dummy_local_ordinal_type_1d_view;
      ComputeResidualVector<MatrixType> 
        compute_residual_vector(amd, A->getCrsGraph().getLocalGraph(), blocksize, interf, 
                                is_async_importer_active ? async_importer->dm2cm : dummy_local_ordinal_type_1d_view);
      
      // norm manager workspace resize
      if (is_norm_manager_active) {
        norm_manager.resize(blocksize, num_vectors);
        norm_manager.setCheckFrequency(check_tol_every);
      }

      // // iterate
      int sweep = 0;
      for (;sweep<max_num_sweeps;++sweep) {
        if (is_y_zero) {
          // pmv := x(lclrow)
          multivector_converter.to_packed_multivector(XX);
        } else {
          if (is_seq_method_requested) {
            // y := x - R y
            Z.doImport(Y, *tpetra_importer, Tpetra::REPLACE);
            compute_residual_vector.run(YY, XX, ZZ);

            // pmv := y(lclrow).
            multivector_converter.to_packed_multivector(YY);
          } else {
            // fused y := x - R y and pmv := y(lclrow); 
            if (overlap_communication_and_computation || !is_async_importer_active) {
              if (is_async_importer_active) async_importer->asyncSendRecv(YY);
              compute_residual_vector.run(pmv, XX, YY, remote_multivector, true);
              if (is_norm_manager_active && norm_manager.checkDone(sweep, tolerance)) {
                if (is_async_importer_active) async_importer->cancel();
                break;
              }
              if (is_async_importer_active) {
                async_importer->syncRecv();
                compute_residual_vector.run(pmv, XX, YY, remote_multivector, false);
              }
            } else {
              if (is_async_importer_active)
                async_importer->syncExchange(YY);
              if (is_norm_manager_active && norm_manager.checkDone(sweep, tolerance)) break;
              compute_residual_vector.run(pmv, XX, YY, remote_multivector);
            }
          }
        }
        
        // pmv := inv(D) pmv.
        solve_tridiags.run();
        
        // y(lclrow) = (b - a) y(lclrow) + a pmv, with b = 1 always.
        multivector_converter.to_scalar_multivector(YY, damping_factor, is_y_zero,
                                                    is_norm_manager_active ? norm_manager.getBuffer() : NULL);
        
        if (is_norm_manager_active) {
          if (sweep + 1 == max_num_sweeps) {
            norm_manager.ireduce(sweep, true);
            norm_manager.checkDone(sweep + 1, tolerance, true);
          } else {
            norm_manager.ireduce(sweep);
          }
        }

        is_y_zero = false;
      }

      //sqrt the norms for the caller's use.
      if (is_norm_manager_active) norm_manager.finalize();
      
      return sweep;
    }


    template<typename MatrixType>
    struct ImplObject { 
      using impl_type = ImplType<MatrixType>;
      using part_interface_type = PartInterface<MatrixType>;
      using block_tridiags_type = BlockTridiags<MatrixType>;
      using amd_type = AmD<MatrixType>;
      using norm_manager_type = NormManager<MatrixType>;
      using async_import_type = AsyncableImport<MatrixType>;
      
      // distructed objects
      Teuchos::RCP<const typename impl_type::tpetra_block_crs_matrix_type> A;
      Teuchos::RCP<const typename impl_type::tpetra_import_type> tpetra_importer;
      Teuchos::RCP<async_import_type> async_importer;
      bool overlap_communication_and_computation;
      
      // copy of Y (mutable to penentrate const)
      mutable typename impl_type::tpetra_multivector_type Z;
      
      // local objects
      part_interface_type part_interface;
      block_tridiags_type block_tridiags; // D
      amd_type a_minus_d; // R = A - D
      mutable typename impl_type::vector_type_1d_view work; // right hand side workspace
      mutable norm_manager_type norm_manager;      
    };
  
  } // namespace BlockTriDiContainerDetails

} // namespace Ifpack2

#endif
