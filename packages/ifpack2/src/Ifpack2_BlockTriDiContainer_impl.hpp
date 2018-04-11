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
      ArrayValueType() = default;
      ArrayValueType(const ArrayValueType &b) {
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
      /// default host execution space
      ///
      typedef Kokkos::DefaultHostExecutionSpace host_execution_space;
      
      ///
      /// tpetra types
      ///
      typedef typename node_type::device_type device_type;      
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
    createBlockCrsTpetraImporter(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A) {
      using impl_type = ImplType<MatrixType>;
      using map_type = typename impl_type::tpetra_map_type;
      using mv_type = typename impl_type::tpetra_block_multivector_type;
      
      const auto g = A->getCrsGraph();  // tpetra crs graph object
      const auto blocksize = A->getBlockSize();
      const auto src = Teuchos::rcp(new map_type(mv_type::makePointMap(*g.getDomainMap(), blocksize)));
      const auto tgt = Teuchos::rcp(new map_type(mv_type::makePointMap(*g.getColMap()   , blocksize)));

      return Teuchos::rcp(new typename ImplType<MatrixType>::tpetra_import_type(src, tgt));
    }

    // Partial replacement for forward-mode MultiVector::doImport. 
    // Permits overlapped communication and computation, but also supports sync'ed. 
    // I'm finding that overlapped comm/comp can give quite poor performance on some
    // platforms, so we can't just use it straightforwardly always.

    // we first try to communicate everything explicitly without relying on
    // gpu aware mpi, which means that data is moved to host explicitly and 
    // communicate among mpi ranks. I am not fully convinced that gpu aware mpi 
    // communication gains efficiency in sparse data
    template<typename MatrixType>
    struct AsyncableImport {
    public:
      using impl_type = ImplType<MatrixType>;

    private:
      ///
      /// MPI wrapper
      ///
#if !defined(HAVE_MPI)
      typedef int MPI_Request;
      typedef int MPI_Comm;
#endif
      /// teuchos mpi type traits does not recorgnize kokkos::complex (impl_scalar_type)
      /// use scalar_type for communication data type enum
      using scalar_type = typename impl_type::scalar_type;

      template <typename T>
      static int isend(const MPI_Comm comm, const T* buf, int count, int dest, int tag, MPI_Request* ireq) {
#ifdef HAVE_MPI
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
#ifdef HAVE_MPI
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
#ifdef HAVE_MPI
        return MPI_Waitany(count, reqs, index, MPI_STATUS_IGNORE);
#else
        return 0;
#endif
      }

      static int waitall(int count, MPI_Request* reqs) {
#ifdef HAVE_MPI
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

      // here we use everything on host space (device memory should be moved into host for communication)
      using host_execution_space = typename impl_type::host_execution_space;
      using int_1d_view_host = Kokkos::View<int*,host_execution_space>;
      using local_ordinal_type_1d_view_host = typename impl_type::local_ordinal_type_1d_view::HostMirror;
      using size_type_1d_view_host = typename impl_type::size_type_1d_view::HostMirror;
      using impl_scalar_type_1d_view_host = typename impl_type::impl_scalar_type_1d_view::HostMirror;
      using impl_scalar_type_2d_view_host = typename impl_type::impl_scalar_type_2d_view::HostMirror;

      // dm2cm permutation should happen on device
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using do_not_initialize_tag = Kokkos::ViewAllocateWithoutInitializing;

      MPI_Comm comm;
      impl_scalar_type_2d_view_host remote_multivector;
      local_ordinal_type blocksize;
      
      template<typename T>
      struct SendRecvPair {
        T send, recv;
      };

      // (s)end and (r)eceive data:
      SendRecvPair<int_1d_view_host> pids;                // mpi ranks
      SendRecvPair<std::vector<MPI_Request> > reqs;       // MPI_Request is pointer, cannot use kokkos view
      SendRecvPair<size_type_1d_view_host> offset;        // offsets to local id list and data buffer
      SendRecvPair<local_ordinal_type_1d_view_host> lids; // local id list
      SendRecvPair<impl_scalar_type_1d_view_host> buffer; // data buffer

      local_ordinal_type_1d_view dm2cm; // permutation

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
        reqs.recv.resize(pids.recv.extent(0));
        reqs.send.resize(pids.send.extent(0));
        
        // construct offsets
        const auto lengths_to = distributor.getLengthsTo();
        offset.send = size_type_1d_view_host(do_not_initialize_tag("offset send"), lengths_to.size() + 1);
        
        const auto lengths_from = distributor.getLengthsFrom();
        offset.recv = size_type_1d_view_host(do_not_initialize_tag("offset recv"), lengths_from.size() + 1);
        
        auto set_offset_values = [](const Teuchos::ArrayView<const size_t> &lens, 
                                    const size_type_1d_view_host &offs) { 
          const Kokkos::RangePolicy<host_execution_space> policy(0,lens.size());
          Kokkos::parallel_scan
          (policy,
           KOKKOS_LAMBDA(const local_ordinal_type &i, size_type &update, const bool &final) {
            if (final) 
              offs(i) = update;
            update += lens[i];
          });
        };
        set_offset_values(lengths_to,   offset.send);
        set_offset_values(lengths_from, offset.recv);
      }

      void createSendRecvIDs(const tpetra_import_type &import) {
        // For each remote PID, the list of LIDs to receive.
        const auto remote_lids = import.getRemoteLIDs();
        lids.recv = local_ordinal_type_1d_view_host(do_not_initialize_tag("lids recv"), remote_lids.size());
        memcpy(lids.recv.data(), remote_lids.getRawPtr(), sizeof(local_ordinal_type)*lids.recv.extent(0));

        // For each export PID, the list of LIDs to send.
        auto epids = import.getExportPIDs();
        auto elids = import.getExportLIDs();
        TEUCHOS_ASSERT(epids.size() == elids.size());
        lids.send = local_ordinal_type_1d_view_host(do_not_initialize_tag("lids send"), elids.size());

        // naive search (not sure if pids or epids are sorted)
        for (local_ordinal_type cnt=0,i=0,iend=pids.send.extent(0);i<iend;++i) {
          const auto pid_send_value = pids.send[i];
          for (local_ordinal_type j=0,jend=epids.size();j<jend;++j)
            if (epids[j] == pid_send_value) lids.send[cnt++] = elids[j];
        }
      }

      void createDataBuffer(const local_ordinal_type &num_vectors) {
        const size_type extent_0 = lids.recv.extent(0)*blocksize;
        const size_type extent_1 = num_vectors;
        if (remote_multivector.extent(0) == extent_0 &&
            remote_multivector.extent(1) == extent_1) {
          // skip
        } else {
          remote_multivector = 
            impl_scalar_type_2d_view_host(do_not_initialize_tag("remote multivector"), extent_0, extent_1);

          const auto send_buffer_size = offset.send[offset.send.extent(0)-1]*blocksize*num_vectors;
          buffer.send = impl_scalar_type_1d_view_host(do_not_initialize_tag("buffer send"), send_buffer_size);

          const auto recv_buffer_size = offset.recv[offset.recv.extent(0)-1]*blocksize*num_vectors;
          buffer.recv = impl_scalar_type_1d_view_host(do_not_initialize_tag("buffer recv"), recv_buffer_size);
        }
      }

      struct ToBuffer {};
      struct ToMultiVector {};

      template<typename PackTag>
      static 
      void copy(const local_ordinal_type_1d_view_host &lids_,
                const impl_scalar_type_1d_view_host &buffer_,
                const local_ordinal_type &ibeg_, 
                const local_ordinal_type &iend_,
                const impl_scalar_type_2d_view_host &multivector_,
                const local_ordinal_type blocksize_) {
        const local_ordinal_type num_vectors = multivector_.extent(1);
        const size_type datasize = blocksize_*sizeof(impl_scalar_type);
        if (num_vectors == 1) {
          const Kokkos::RangePolicy<host_execution_space> policy(ibeg_, iend_);
          Kokkos::parallel_for
            (policy,
             KOKKOS_LAMBDA(const local_ordinal_type &i) {
              auto aptr = buffer_.data() + blocksize_*i;
              auto bptr = multivector_.data() + blocksize_*lids_(i);
              if (std::is_same<PackTag,ToBuffer>::value) memcpy(aptr, bptr, datasize);
              else                                       memcpy(bptr, aptr, datasize);
            });
        } else { 
          const local_ordinal_type diff = iend_ - ibeg_;
          Kokkos::MDRangePolicy
            <host_execution_space, Kokkos::Rank<2>, Kokkos::IndexType<local_ordinal_type> > 
            policy( { ibeg_, 0 }, { iend_, num_vectors } );
          Kokkos::parallel_for
            (policy,
             KOKKOS_LAMBDA(const local_ordinal_type &i,
                           const local_ordinal_type &j) { 
              auto aptr = buffer_.data() + blocksize_*(i + diff*j);
              auto bptr = &multivector_(blocksize_*lids_(i), j);
              if (std::is_same<PackTag,ToBuffer>::value) memcpy(aptr, bptr, datasize);
              else                                       memcpy(bptr, aptr, datasize);
            });
        }
      }

    public:
      AsyncableImport (const Teuchos::RCP<const tpetra_map_type>& src_map,
                       const Teuchos::RCP<const tpetra_map_type>& tgt_map, 
                       const local_ordinal_type blocksize_,
                       const local_ordinal_type_1d_view dm2cm_) {
        blocksize = blocksize_;
        dm2cm = dm2cm_;

        const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(tgt_map->getComm());
        TEUCHOS_ASSERT(!mpi_comm.is_null());
        comm = *mpi_comm->getRawMpiComm();

        const tpetra_import_type import(src_map, tgt_map);

        createMpiRequests(import);
        createSendRecvIDs(import);
        createDataBuffer(1); // Optimistically set up for 1-vector case.
      }

      void asyncSendRecv(const impl_scalar_type_2d_view_host &mv) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Async_Setup::async_setup");
#endif
#ifdef HAVE_MPI
        // constants and reallocate data buffers if necessary
        const local_ordinal_type num_vectors = mv.extent(1);
        const local_ordinal_type mv_blocksize = blocksize*num_vectors;
        createDataBuffer(num_vectors);

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
#ifdef HAVE_MPI
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i)
          MPI_Cancel(&reqs.recv[i]);
#endif
      }

      void syncRecv() {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Async_Setup::sync_receive");
#endif
#ifdef HAVE_MPI
        // receive async.
        for (local_ordinal_type i=0,iend=pids.recv.extent(0);i<iend;++i) {
          local_ordinal_type idx;
          waitany(pids.recv.extent(0), reqs.recv.data(), &idx);
          copy<ToMultiVector>(lids.recv, buffer.recv, offset.recv(idx), offset.recv(idx+1),
                              remote_multivector, blocksize);
        }
        // wait on the sends to match all Isends with a cleanup operation.
        waitall(reqs.send.size(), reqs.send.data());
#endif
      }

      void syncExchange(const impl_scalar_type_2d_view_host &mv) {
        asyncSendRecv(mv);
        syncRecv();
      }
    };

    ///
    /// setup async importer
    ///
    template<typename MatrixType>
    Teuchos::RCP<const AsyncableImport<MatrixType> > 
    createBlockCrsAsyncImporter(const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A) {
      using impl_type = ImplType<MatrixType>;
      using map_type = typename impl_type::tpetra_map_type;
      using mv_type = typename impl_type::tpetra_block_multivector_type;
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
          = Teuchos::rcp(new map_type(invalid, gids.data(), gids.size(), 0, domain_map->getComm()));
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
                         const bool overlap_comm) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::SymbolicPhase");
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
          amd.A_colindsub = local_ordinal_type_1d_view(Kokkos::ViewAllocateWithoutInitializing("amd.A_colindsub"), R_nnz_owned);

          const auto R_rowptr = Kokkos::create_mirror_view(amd.rowptr); 
          const auto R_A_colindsub = Kokkos::create_mirror_view(amd.A_colindsub);
          
          amd.rowptr_remote = size_type_1d_view("amd.rowptr_remote", overlap_comm ? nrows + 1 : 0);
          amd.A_colindsub_remote = local_ordinal_type_1d_view(Kokkos::ViewAllocateWithoutInitializing("amd.A_colindsub_remote"), R_nnz_remote);
          
          const auto R_rowptr_remote = Kokkos::create_mirror_view(amd.rowptr_remote);
          const auto R_A_colindsub_remote = Kokkos::create_mirror_view(amd.A_colindsub_remote);
          
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
                // exclusive scan will be performed later
                if (!overlap_comm || lc < nrows) 
                  ++R_rowptr(lr);
                else 
                  ++R_rowptr_remote(lr);
              }
            });
          
          // exclusive scan
          typedef ArrayValueType<size_type,2> update_type;
          Kokkos::parallel_scan
            (Kokkos::RangePolicy<host_execution_space>(0,nrows+1),
             KOKKOS_LAMBDA(const local_ordinal_type &lr, 
                           update_type &update, 
                           const bool &final) {
              update_type val;
              val.v[0] = R_rowptr(lr);
              if (overlap_comm)
                val.v[1] = R_rowptr_remote(lr);
              
              if (final) {
                R_rowptr(lr) = update.v[0];
                if (overlap_comm)
                  R_rowptr_remote(lr) = update.v[1];
                
                if (lr < nrows) {
                  const local_ordinal_type ri0 = lclrow2idx[lr];
                  const local_ordinal_type pi0 = rowidx2part(ri0);
                  
                  size_type cnt_rowptr = R_rowptr(lr);
                  size_type cnt_rowptr_remote = overlap_comm ? R_rowptr_remote(lr) : 0; // when not overlap_comm, this value is garbage

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
                      R_A_colindsub(cnt_rowptr++) = row_entry;
                    else 
                      R_A_colindsub_remote(cnt_rowptr_remote++) = row_entry;
                  }
                }
              } 
              update += val;
            });

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
      ConstUnmanaged<size_type_1d_view> pack_td_ptr, flat_td_ptr;
      ConstUnmanaged<local_ordinal_type_1d_view> A_colindsub;
      Unmanaged<vector_type_3d_view> vector_values;
      // diagonal safety
      const magnitude_type tiny;
      // shared information
      const local_ordinal_type blocksize, blocksize_square;

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
        A_values(const_cast<block_crs_matrix_type*>(A_.get())->
                 template getValues<typename device_type::memory_space>()),
        // block tridiags 
        flat_td_ptr(btdm_.flat_td_ptr), 
        pack_td_ptr(btdm_.pack_td_ptr), 
        A_colindsub(btdm_.A_colindsub),
        vector_values(btdm_.values),
        blocksize(btdm_.values.extent(1)),
        blocksize_square(blocksize*blocksize),
        // diagonal weight to avoid zero pivots
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
      
      KOKKOS_INLINE_FUNCTION 
      void 
      factorize (const local_ordinal_type& packidx) const {
        using namespace KokkosBatched::Experimental;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        
        // subview pattern
        auto A = Kokkos::subview(vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        auto B = A;
        auto C = A;
        
        auto i0 = pack_td_ptr(packptr(packidx));
        const local_ordinal_type nrows 
          = BlockTridiags<MatrixType>::IndexToRow(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;
        
        A.assign_data( &vector_values(i0,0,0) );
        SerialLU<Algo::LU::Unblocked>::invoke(A, tiny);
        for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
          B.assign_data( &vector_values(i0+1,0,0) );
          SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsm::Blocked>
            ::invoke(one, A, B);
          C.assign_data( &vector_values(i0+2,0,0) );
          SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsm::Blocked>
            ::invoke(one, A, C);
          A.assign_data( &vector_values(i0+3,0,0) );
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
        Kokkos::RangePolicy<typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
        //#endif
        Kokkos::parallel_for(policy, *this);
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
      using device_type = typename impl_type::device_type;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using impl_scalar_type = typename impl_type::impl_scalar_type;
      using magnitude_type = typename impl_type::magnitude_type;
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
      impl_scalar_type damping_factor;

      // packed multivector output (or input)
      Unmanaged<vector_type_3d_view> packed_multivector;
      Unmanaged<impl_scalar_type_2d_view> scalar_multivector;

      struct ToPackedMultiVectorTag       { enum : int { id = 0 }; };
      struct ToScalarMultiVectorFirstTag  { enum : int { id = 1 }; };
      struct ToScalarMultiVectorSecondTag { enum : int { id = 2 }; };

      template<typename TagType>
      KOKKOS_INLINE_FUNCTION
      void copy_multivectors(const local_ordinal_type &j, 
                             const local_ordinal_type &vi, 
                             const local_ordinal_type &pri, 
                             const local_ordinal_type &nrow,  
                             const local_ordinal_type &ri0) const {
        if (TagType::id == 0) { // ToPackedMultiVectorTag
          for (local_ordinal_type col=0;col<num_vectors;++col) 
            for (local_ordinal_type i=0;i<blocksize;++i)
              packed_multivector(pri, i, col)[vi] = scalar_multivector(blocksize*lclrow(ri0+j)+i,col);
        } else if (TagType::id > 0) { //ToScalarMultiVector
          const impl_scalar_type zero(0), df = damping_factor;
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
                                       const local_ordinal_type &nrow,  
                                       const local_ordinal_type &ri0,
                                       /* */ magnitude_type *norm) const {
        if (TagType::id > 0) { //ToScalarMultiVector
          const impl_scalar_type zero(0), df = damping_factor;
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

      MultiVectorConverter() = default;
      MultiVectorConverter(const MultiVectorConverter &b) = default;
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

      template<typename TagType>
      KOKKOS_INLINE_FUNCTION
      void 
      operator() (const TagType&, const local_ordinal_type &packidx, magnitude_type* const norm = NULL) const {
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
          for (local_ordinal_type vi=0;vi<npacks;++vi) {
            if (norm == NULL) copy_multivectors<TagType>(j, vi, pri, nrows[vi], ri0[vi]);
            else              copy_multivectors_with_norm<TagType>(j, vi, pri, nrows[vi], ri0[vi], norm);
          }
        }
      }
      
      template<typename TpetraLocalViewType>
      void to_packed_multivector(const TpetraLocalViewType &scalar_multivector_) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::MultiVectorConverter::ToPackedMultiVector");
#endif
        value_count = 0;
        scalar_multivector = scalar_multivector_;
        const Kokkos::RangePolicy
          <ToPackedMultiVectorTag,typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
        Kokkos::parallel_for(policy, *this);
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

        if (is_vectors_zero) {
          const Kokkos::RangePolicy
            <ToScalarMultiVectorFirstTag,typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
          if (norm == NULL)  Kokkos::parallel_for(policy, *this);
          else               Kokkos::parallel_reduce(policy, *this, norm);              
        } else {
          const Kokkos::RangePolicy
            <ToScalarMultiVectorSecondTag,typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
          if (norm == NULL)  Kokkos::parallel_for(policy, *this);
          else               Kokkos::parallel_reduce(policy, *this, norm);              
        }
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
        using namespace KokkosBatched::Experimental;
        
        // constant
        const auto one = Kokkos::ArithTraits<magnitude_type>::one();
        const local_ordinal_type num_vectors = X_vector_values.extent(2);

        // subview pattern
        auto A = Kokkos::subview(D_vector_values, 0, Kokkos::ALL(), Kokkos::ALL());
        auto B = A;
        auto X1 = Kokkos::subview(X_vector_values, 0, Kokkos::ALL(), 0);
        auto X2 = X1;

        // index counting
        const local_ordinal_type partidx = packptr(packidx);
        size_type i0 = pack_td_ptr(partidx);
        local_ordinal_type r0 = part2packrowidx0(partidx);
        const local_ordinal_type nrows = part2packrowidx0(packptr(packidx+1)) - r0;

        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // solve Lx = x
          A.assign_data( &D_vector_values(i0,0,0) );
          X1.assign_data( &X_vector_values(r0,0,col) );
          SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsv::Blocked>
            ::invoke(one, A, X1);
          for (local_ordinal_type i=1;i<nrows;++i,i0+=3) {
            B.assign_data( &D_vector_values(i0+2,0,0) );
            X2.assign_data( &X_vector_values(++r0,0,col) );
            SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>
              ::invoke(-one, B, X1, one, X2);
            A.assign_data( &D_vector_values(i0+3,0,0) );
            SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Algo::Trsv::Blocked>
              ::invoke(one, A, X2);
            X1.assign_data( X2.data() );
          }
          
          // solve Ux = x
          SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsv::Blocked>
            ::invoke(one, A, X1);
          for (local_ordinal_type i=nrows;i>1;--i) {
            i0 -= 3;
            B.assign_data( &D_vector_values(i0+1,0,0) );
            X2.assign_data( &X_vector_values(--r0,0,col) );          
            SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>
              ::invoke(-one, B, X1, one, X2); 
            A.assign_data( &D_vector_values(i0,0,0) );          
            SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::Trsv::Blocked>
              ::invoke(one, A, X2);
            X1.assign_data( X2.data() );
          }
        }
      }
      
      void run() {
        // #if defined(HAVE_IFPACK2_CUDA) && defined (KOKKOS_ENABLE_CUDA)
        //         Kokkos::TeamPolicy<device_type> policy(packptr.extent(0) - 1, Kokkos::AUTO(), vector_length);
        // #else
        Kokkos::RangePolicy<typename device_type::execution_space> policy(0, packptr.extent(0) - 1);
        //#endif
        Kokkos::parallel_for(policy, *this);
      }
    }; 

    ///
    /// compute local residula vector y = b - R x 
    ///
    template<typename MatrixType>
    struct ComputeResidualVector {
    public:
      using impl_type = ImplType<MatrixType>;
      using device_type = typename impl_type::device_type;
      using local_ordinal_type = typename impl_type::local_ordinal_type;
      using size_type = typename impl_type::size_type;
      using magnitude_type = typename impl_type::magnitude_type;
      /// views
      using local_ordinal_type_1d_view = typename impl_type::local_ordinal_type_1d_view;
      using size_type_1d_view = typename impl_type::size_type_1d_view; 
      using tpetra_block_access_view_type = typename impl_type::tpetra_block_access_view_type; // block crs (layout right)
      using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view; 
      using impl_scalar_type_2d_view = typename impl_type::impl_scalar_type_2d_view; // block multivector (layout left)

      /// team policy member type (used in cuda)
      //using member_type = typename Kokkos::TeamPolicy<device_type>::member_type;      
 
    private:
      ConstUnmanaged<impl_scalar_type_2d_view> b;
      ConstUnmanaged<impl_scalar_type_2d_view> x;
      /* */Unmanaged<impl_scalar_type_2d_view> y;

      // AmD information
      ConstUnmanaged<size_type_1d_view> rowptr;
      ConstUnmanaged<local_ordinal_type_1d_view> colindsub;
      ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

      // block crs graph information
      ConstUnmanaged<size_type_1d_view> A_rowptr;
      ConstUnmanaged<local_ordinal_type_1d_view> A_colind;
      const local_ordinal_type blocksize, blocksize_square;

      // block access
      ConstUnmanaged<tpetra_block_access_view_type> A_null_block;
      
    public:
      template<typename LocalCrsGraphType>
      ComputeResidualVector(const AmD<MatrixType> &amd,
                            const LocalCrsGraphType &graph,
                            const local_ordinal_type blocksize_) 
        : rowptr(amd.rowptr),
          colindsub(amd.A_colindsub),
          tpetra_values(amd.tpetra_values),
          A_rowptr(graph.row_map),
          A_colind(graph.entries),
          blocksize(blocksize_), 
          blocksize_square(blocksize_*blocksize_),
          A_null_block(NULL, blocksize_, blocksize_) {}
      
      ///
      /// host serial 
      ///
      KOKKOS_INLINE_FUNCTION 
      void 
      operator() (const local_ordinal_type& i) const {
        using namespace KokkosBatched::Experimental;

        // constants
        const magnitude_type one(1);
        const Kokkos::pair<local_ordinal_type,local_ordinal_type> block_range(0, blocksize);
        const local_ordinal_type num_vectors = y.extent(1);

        // subview pattern
        auto bb = Kokkos::subview(b, block_range, 0);
        auto xx = bb;
        auto yy = Kokkos::subview(y, block_range, 0);
        auto A_block = A_null_block; // A_block is assumed to be compact (a dimension is equal to stride)
        
        const local_ordinal_type row = i*blocksize;
        for (local_ordinal_type col=0;col<num_vectors;++col) {
          // y := b
          yy.assign_data(&y(row, col));
          bb.assign_data(&b(row, col));
          SerialCopy<Trans::NoTranspose>::invoke(bb, yy);
        
          // y -= Rx
          const size_type A_k0 = A_rowptr[i];
          for (size_type k=rowptr[i];k<rowptr[i+1];++k) {
            const size_type j = A_k0 + colindsub[k];
            A_block.assign_data( &tpetra_values(j*blocksize_square) );
            xx.assign_data( &x(A_colind[j]*blocksize, col) );
            SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>::invoke(-one, A_block, xx, one, yy);  
          }
        }
      }
      
      // y = b - Rx;
      template<typename MultiVectorLocalViewTypeY,
               typename MultiVectorLocalViewTypeB,
               typename MultiVectorLocalViewTypeX>
      void run(const MultiVectorLocalViewTypeY &y_, 
               const MultiVectorLocalViewTypeB &b_, 
               const MultiVectorLocalViewTypeX &x_) {
        y = y_; b = b_; x = x_; 
        Kokkos::RangePolicy<typename device_type::execution_space> policy(0, rowptr.extent(0) - 1);
        Kokkos::parallel_for(policy, *this);
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
#ifdef HAVE_MPI
      MPI_Request mpi_request_;
      MPI_Comm comm_;
#endif
      std::vector<magnitude_type> work_;

    public:
      NormManager() = default;
      NormManager(const NormManager &b) = default;
      NormManager(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
        sweep_step_ = 1;
        collective_ = comm->getSize() > 1;
        if (collective_) {
#ifdef HAVE_MPI
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
#ifdef HAVE_MPI
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
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NormManager::checkDone");
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
    applyInverseJacobi(// tpetra importer
                       const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_block_crs_matrix_type> &A,
                       const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_import_type> &importer,
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
      const bool is_norm_manager_active = tol > Kokkos::ArithTraits<magnitude_type>::zero();
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
      const auto YY = Y.template getLocalView<typename device_type::memory_space>();
      const auto ZZ = Z.template getLocalView<typename device_type::memory_space>();

      MultiVectorConverter<MatrixType> multivector_converter(interf, pmv);
      SolveTridiags<MatrixType> solve_tridiags(interf, btdm, pmv);

      ComputeResidualVector<MatrixType> 
        compute_residual_vector(amd, A->getCrsGraph().getLocalGraph(), blocksize);
      
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
          // y := x - R y; in this notation, zz = yy, yy = xx - amd zz
          Z.doImport(Y, *importer, Tpetra::REPLACE);
          compute_residual_vector.run(YY, XX, ZZ);
          // pmv := y(lclrow).
          multivector_converter.to_packed_multivector(YY);
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
      
      printf("sweep = %d, X extents %d, %d, Y extents %d %d\n", 
             sweep, XX.extent(0), XX.extent(1), YY.extent(0), YY.extent(1));
      
      for (int j=0;j<XX.extent(1);++j) 
        for (int i=0;i<XX.extent(0);++i) 
          std::cout << " (i,j) = " << i << " " << j << " X = " << XX(i,j) << " Y = " << YY(i,j) << "\n";
      
      return sweep;
    }
  
  } // namespace BlockTriDiContainerDetails

} // namespace Ifpack2

#endif
