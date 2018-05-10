#include <cassert>
#include <vector>
#include <algorithm>

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosBatched_Util.hpp"

#define TEST_ASSERT(m, success)                                 \
  if ( !(m)) {                                                  \
    success = false;                                            \
    printf("FAILED: %s, at %d, %s\n", #m, __LINE__, __FILE__);  \
  }

namespace KokkosBatched {
  using namespace Experimental;

  namespace Test {

    typedef int ordinal_type;
    typedef int size_type;
    typedef double scalar_type;
#define BLOCKCRS_MAX_BLOCKSIZE 32
#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
    
    double LU_FlopCount(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      if (m > n)
        return (FLOP_MUL*(0.5*m*n*n-(1.0/6.0)*n*n*n+0.5*m*n-0.5*n*n+(2.0/3.0)*n) +
                FLOP_ADD*(0.5*m*n*n-(1.0/6.0)*n*n*n-0.5*m*n+        (1.0/6.0)*n));
      else
        return (FLOP_MUL*(0.5*n*m*m-(1.0/6.0)*m*m*m+0.5*n*m-0.5*m*m+(2.0/3.0)*m) +
                FLOP_ADD*(0.5*n*m*m-(1.0/6.0)*m*m*m-0.5*n*m+        (1.0/6.0)*m));
    }

    double Trsm_Lower_FlopCountLower(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    double Trsm_Upper_FlopCountUpper(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    double Gemm_FlopCount(int mm, int nn, int kk) {
      double m = (double)mm;    double n = (double)nn;    double k = (double)kk;
      return (FLOP_MUL*(m*n*k) +
              FLOP_ADD*(m*n*k));
    }

    template <typename aViewType, typename bViewType>
    double compute_relative_diff(const aViewType a, 
                                 const bViewType b) {
      // Bring the vectors to the host. This is just a correctness checker.
      auto aa = Kokkos::create_mirror_view(a); Kokkos::deep_copy(aa, a);
      auto bb = Kokkos::create_mirror_view(b); Kokkos::deep_copy(bb, b);

      double diff2 = 0, norm2 = 0;
      for (ordinal_type i=0,iend=aa.extent(0);i<iend;++i)
        for (ordinal_type j=0,jend=aa.extent(1);j<jend;++j)
          for (ordinal_type k=0,kend=aa.extent(2);k<kend;++k)
            for (ordinal_type l=0,lend=aa.extent(3);l<lend;++l) {
              const double 
                val  = aa(i,j,k,l),
                diff = aa(i,j,k,l) - bb(i,j,k,l);
              diff2 += diff*diff;
              norm2 += val*val;
            }

      return std::sqrt(diff2/norm2);
    }

    // Representation of a structured block mesh. The fastest index is k.
    struct StencilShape { 
      enum Enum { cross }; 
    };

    struct StructuredBlock {
      const ordinal_type ni, nj, nk;

      StructuredBlock (const ordinal_type ni_, 
                       const ordinal_type nj_, 
                       const ordinal_type nk_)
        : ni(ni_), nj(nj_), nk(nk_), _njnk(nj_*nk_) {}

      KOKKOS_INLINE_FUNCTION 
      size_type size () const { return ni*nj*nk; }

      KOKKOS_INLINE_FUNCTION 
      size_type ij2id (const ordinal_type i, 
                       const ordinal_type j) const { 
        return i*nj + j;
      }

      KOKKOS_INLINE_FUNCTION 
      void id2ij (const size_type id, 
                  ordinal_type& i, 
                  ordinal_type& j) const {
        i = id / nj;
        j = id % nj;
      }

      KOKKOS_INLINE_FUNCTION 
      size_type ijk2id (const ordinal_type i, 
                        const ordinal_type j, 
                        const ordinal_type k) const { 
        return (i*nj + j)*nk + k; 
      }
      
      KOKKOS_INLINE_FUNCTION 
      void id2ijk (const size_type id, 
                   ordinal_type& i, 
                   ordinal_type& j, 
                   ordinal_type& k) const {
        i = id / _njnk;
        k = id % _njnk;
        j = k / nk;
        k = k % nk;
      }
    private:
      const ordinal_type _njnk;
    };
    
    template <typename ExecSpace, typename ArrayLayout>
    struct CrsGraph {
      typedef ExecSpace exec_space;
      typedef ArrayLayout array_layout;      

      typedef Kokkos::View<size_type*,   array_layout,exec_space> row_ptr_type;
      typedef Kokkos::View<ordinal_type*,array_layout,exec_space> row_idx_type;
      typedef Kokkos::View<ordinal_type*,array_layout,exec_space> col_idx_type;

      row_ptr_type rowptr;
      row_idx_type rowidx;
      col_idx_type colidx;
      
      CrsGraph () 
        : rowptr("rowptr", 1), 
          rowidx("rowidx", 0),
          colidx("colidx", 0) {}

      KOKKOS_INLINE_FUNCTION
      bool isEmpty() const { 
        return (rowptr.extent(0) <= 1 || colidx.extent(0) == 0 || rowidx.extent(0) == 0); 
      }
      
      KOKKOS_INLINE_FUNCTION
      ordinal_type NumRows() const { 
        return (isEmpty() ? 0 : static_cast<ordinal_type>(rowptr.extent(0)) - 1); 
      }
      
      KOKKOS_INLINE_FUNCTION
      size_type NumNonZeros() const { 
        return (isEmpty() ? 0 : static_cast<size_type>(colidx.extent(0))); 
      }
    };

    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    CrsGraph<DstSpace,ArrayLayout>
    create_mirror(const CrsGraph<SrcSpace,ArrayLayout> src) {
      CrsGraph<DstSpace,ArrayLayout> dst;

      dst.rowptr = Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.rowptr);
      dst.rowidx = Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.rowidx);
      dst.colidx = Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.colidx);

      return dst;
    }

    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    void
    deep_copy(const CrsGraph<DstSpace,ArrayLayout> dst, const CrsGraph<SrcSpace,ArrayLayout> src) {
      Kokkos::deep_copy(dst.rowptr, src.rowptr);
      Kokkos::deep_copy(dst.rowidx, src.rowidx);
      Kokkos::deep_copy(dst.colidx, src.colidx);
    }
    
    // Given a structured block and a stencil (at present, just a 3D 1-hop cross),
    // construct a corresponding CRS graph.
    template<typename ArrayLayout>
    CrsGraph<Kokkos::DefaultHostExecutionSpace,ArrayLayout>
    create_graph_host_for_structured_block(const StructuredBlock mesh, 
                                           const StencilShape::Enum shape) {
      CrsGraph<Kokkos::DefaultHostExecutionSpace,ArrayLayout> graph;

      Kokkos::resize(graph.rowptr, mesh.size()+1);
      graph.rowptr[0] = 0;
      
      std::vector<ordinal_type> colidx, rowidx;
      switch (shape) {
      case StencilShape::cross:
        for (ordinal_type c=0;c<mesh.size();++c) {
          ordinal_type i, j, k, n = 0;

          mesh.id2ijk(c, i, j, k);
          
          rowidx.push_back(c); colidx.push_back(c); ++n;
          if (i > 0)         { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i-1, j, k  )); ++n; }
          if (i+1 < mesh.ni) { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i+1, j, k  )); ++n; }
          if (j > 0)         { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i, j-1, k  )); ++n; }
          if (j+1 < mesh.nj) { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i, j+1, k  )); ++n; }
          if (k > 0)         { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i, j,   k-1)); ++n; }
          if (k+1 < mesh.nk) { rowidx.push_back(c); colidx.push_back(mesh.ijk2id(i, j,   k+1)); ++n; }
          graph.rowptr[c+1] = graph.rowptr[c] + n;
        }
        break;
      }
      assert(graph.rowptr[mesh.size()] == static_cast<size_type>(colidx.size()));
      assert(graph.rowptr[mesh.size()] == static_cast<size_type>(rowidx.size()));

      for (ordinal_type c=0;c<mesh.size();++c)
        std::sort(colidx.begin() + graph.rowptr[c], colidx.begin() + graph.rowptr[c+1]);

      const ordinal_type nnz = graph.rowptr[mesh.size()];
      Kokkos::resize(graph.colidx, nnz);
      Kokkos::resize(graph.rowidx, nnz);
      for (ordinal_type c=0;c<nnz;++c) {
        graph.colidx[c] = colidx[c];
        graph.rowidx[c] = rowidx[c];
      }
      return graph;
    }

    template <typename ExeSpace, typename ArrayLayout>
    class BlockCrsMatrix {
    public:
      typedef ExeSpace exec_space;
      typedef ArrayLayout array_layout;
      typedef Test::CrsGraph<exec_space,array_layout> crs_graph_type;

      typedef scalar_type value_type;
      typedef Kokkos::View<scalar_type***,array_layout,exec_space> value_array_type;

    private:
      crs_graph_type _graph;
      ordinal_type _blocksize;
      value_array_type _values;

    public:
      BlockCrsMatrix() 
        : _graph(), 
          _blocksize(),
          _values() {} 

      BlockCrsMatrix(const BlockCrsMatrix &b) 
        : _graph(b._graph), 
          _blocksize(b._blocksize),
          _values(b._values) {} 

      BlockCrsMatrix (const crs_graph_type graph, 
                      const ordinal_type blocksize )
        : _graph(graph), 
          _blocksize(blocksize),
          _values("BlockCrsMatrix::_values", 
                  _graph.NumNonZeros(), _blocksize, _blocksize) {}

      BlockCrsMatrix (const crs_graph_type graph, 
                      const ordinal_type blocksize,
                      const value_array_type values)
        : _graph(graph), 
          _blocksize(blocksize),
          _values(values) {}

      ordinal_type BlockSize() const { return _blocksize; }
      crs_graph_type CrsGraph() const { return _graph; }
      value_array_type Values() const { return _values; }
    };

    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    BlockCrsMatrix<DstSpace,ArrayLayout>
    create_mirror(const BlockCrsMatrix<SrcSpace,ArrayLayout> src) {
      const auto graph = create_mirror<DstSpace>(src.CrsGraph());
      const auto blocksize = src.BlockSize();
      const auto values = Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.Values());
      return BlockCrsMatrix<DstSpace,ArrayLayout>(graph, blocksize, values);
    }

    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    void
    deep_copy(const BlockCrsMatrix<DstSpace,ArrayLayout> dst, 
              const BlockCrsMatrix<SrcSpace,ArrayLayout> src) {
      deep_copy(dst.CrsGraph(), src.CrsGraph());
      Kokkos::deep_copy(dst.Values(), src.Values());
    }
    
    template<typename ArrayLayout>
    void fill_block_crs_matrix_host(BlockCrsMatrix<Kokkos::DefaultHostExecutionSpace,ArrayLayout> A) {
      // extract graph and blocksizes
      const auto graph = A.CrsGraph();
      const auto values = A.Values();
      const ordinal_type blocksize = A.BlockSize();

      scalar_type 
        tmp[BLOCKCRS_MAX_BLOCKSIZE*BLOCKCRS_MAX_BLOCKSIZE], //[blocksize*blocksize], 
        diag_block[BLOCKCRS_MAX_BLOCKSIZE][BLOCKCRS_MAX_BLOCKSIZE], //[blocksize][blocksize], 
        offdiag_block[BLOCKCRS_MAX_BLOCKSIZE][BLOCKCRS_MAX_BLOCKSIZE]; //[blocksize][blocksize];
      
      Random<scalar_type> random;

      // for diagonal block, make spd
      {
        const ordinal_type iend = blocksize*blocksize;
        for (ordinal_type i=0;i<iend;++i) 
          tmp[i] = 2*(random.value() - 0.5);
        
        for (ordinal_type i=0;i<blocksize;++i) 
          for (ordinal_type j=i;j<blocksize;++j) {
            diag_block[i][j] = 0;
            for (ordinal_type k=0;k<blocksize;++k) 
              diag_block[i][j] += tmp[i*blocksize+k]*tmp[j*blocksize+k];
            if (i != j) diag_block[j][i]  = diag_block[i][j];    // symmetrize
            else        diag_block[i][j] *= 0.5*blocksize; // improve condition
          }
      } 
      
      {
        // for off diagonal; down-weight off-diag blocks to improve conditioning.
        for (ordinal_type i=0;i<blocksize;++i)
          for (ordinal_type j=0;j<blocksize;++j) 
            offdiag_block[i][j] = 0.1 * 2*(random.value() - 0.5);
      }
      
      for (ordinal_type r=0;r<graph.NumRows();++r) {
        // random number generator (-1, 1)
        const ordinal_type cbegin = graph.rowptr(r), cend = graph.rowptr(r+1);
        for (ordinal_type c=cbegin;c<cend;++c) {
          auto block = Kokkos::subview(values, c, Kokkos::ALL(), Kokkos::ALL());
          
          if (graph.colidx(c) == r) {
            for (ordinal_type i=0;i<blocksize;++i) 
              for (ordinal_type j=i;j<blocksize;++j) 
                block(i,j) = diag_block[i][j];
          } else {
            // for off diagonal; down-weight off-diag blocks to improve conditioning.
            for (ordinal_type i=0;i<blocksize;++i)
              for (ordinal_type j=0;j<blocksize;++j) 
                block(i,j) = offdiag_block[i][j];
          }
          
        }
      }
    }
    
    // nrhs should go after blocksize to match matrix dimensions consistently
    template <typename ExeSpace, typename ArrayLayout>
    class BlockMultiVector {
    public:
      typedef ExeSpace exec_space;
      typedef ArrayLayout array_layout;

      typedef scalar_type value_type;
      typedef Kokkos::View<scalar_type***,array_layout,exec_space> value_array_type;

    private:
      value_array_type _values;

    public:
      BlockMultiVector(const ordinal_type nvecs,
                       const ordinal_type nrows,                       
                       const ordinal_type blocksize )
        : _values("BlockMultiVector::_values", nvecs, nrows, blocksize) {}

      BlockMultiVector(const value_array_type values) 
        : _values(values) {}

      ordinal_type NumVectors() const { return _values.extent(0); }
      ordinal_type NumRows() const { return _values.extent(1); }
      ordinal_type BlockSize() const { return _values.extent(2); }

      value_array_type Values() const { return _values; }
    };
    
    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    BlockMultiVector<DstSpace,ArrayLayout>
    create_mirror(const BlockMultiVector<SrcSpace,ArrayLayout> src) {
      return BlockMultiVector<DstSpace,ArrayLayout>
        (Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.Values()));
    }
    
    template<typename DstSpace, typename SrcSpace, typename ArrayLayout>
    inline 
    void
    deep_copy(const BlockMultiVector<DstSpace,ArrayLayout> dst, 
              const BlockMultiVector<SrcSpace,ArrayLayout> src) {
      Kokkos::deep_copy(dst.Values(), src.Values());
    }
    
    template<typename ArrayLayout>
    void fill_block_multi_vector_host(BlockMultiVector<Kokkos::DefaultHostExecutionSpace,ArrayLayout> B) {
      const ordinal_type 
        jend = B.NumVectors(), 
        iend = B.NumRows(), 
        kend = B.BlockSize();
      
      auto B_val = B.Values();
      
      for (ordinal_type j=0;j<jend;++j) 
        for (ordinal_type i=0;i<iend;++i) 
          for (ordinal_type k=0;k<kend;++k) 
            B_val(j, i, k) = static_cast<double>((i+j+k)%7) - 3;
    }
    
    template <typename ExecSpace, typename ValueType, typename ArrayLayout>
    class BlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef Kokkos::View<value_type****,array_layout,exec_space> value_array_type;
      
    private:
      const ordinal_type _ntridiags, _nrows, _blocksize;
      // A B
      // C
      value_array_type _A, _B, _C;
      
    public:

      BlockTridiagMatrices (const ordinal_type ntridiags,
                            const ordinal_type nrows,
                            const ordinal_type blocksize)
        : _ntridiags(ntridiags), 
          _nrows(nrows),
          _blocksize(blocksize), 
          _A("BlockTridiagMatrix::_A", _ntridiags, _nrows,   _blocksize, _blocksize),
          _B("BlockTridiagMatrix::_B", _ntridiags, _nrows-1, _blocksize, _blocksize),
          _C("BlockTridiagMatrix::_C", _ntridiags, _nrows-1, _blocksize, _blocksize) {}

      BlockTridiagMatrices (const ordinal_type ntridiags,
                            const ordinal_type nrows,
                            const ordinal_type blocksize,
                            const value_array_type A,
                            const value_array_type B,
                            const value_array_type C)
        : _ntridiags(ntridiags), 
          _nrows(nrows),
          _blocksize(blocksize), 
          _A(A),
          _B(B),
          _C(C) {}

      value_array_type A() const { return _A; }
      value_array_type B() const { return _B; }
      value_array_type C() const { return _C; }
      
      ordinal_type BlockSize() const { return _blocksize; }
      ordinal_type NumRows() const { return _nrows; }
      ordinal_type NumTridiagMatrices() const { return _ntridiags; }
    };

    template<typename ExecSpace, typename ValueType, typename ArrayLayout>
    BlockTridiagMatrices<ExecSpace,ValueType,ArrayLayout> 
    create_block_tridiag_matrices(const ordinal_type ntridiags, 
                                  const ordinal_type nrows,
                                  const ordinal_type blocksize) {
      return BlockTridiagMatrices<ExecSpace,ValueType,ArrayLayout>
        (adjustDimension<ValueType>(ntridiags), nrows, blocksize);
    }

    template<typename DstSpace, typename SrcSpace, typename ValueType, typename ArrayLayout>
    inline 
    BlockTridiagMatrices<DstSpace,ValueType,ArrayLayout>
    create_mirror(const BlockTridiagMatrices<SrcSpace,ValueType,ArrayLayout> src) {
      return BlockTridiagMatrices<DstSpace,ValueType,ArrayLayout>
        (src.NumTridiagMatrices(),
         src.NumRows(),
         src.BlockSize(),
         Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.A()),
         Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.B()),
         Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.C()));
    }
    
    template<typename DstSpace, typename SrcSpace, typename ValueType, typename ArrayLayout>
    inline 
    void
    deep_copy(const BlockTridiagMatrices<DstSpace,ValueType,ArrayLayout> dst, 
              const BlockTridiagMatrices<SrcSpace,ValueType,ArrayLayout> src) {
      Kokkos::deep_copy(dst.A(), src.A());
      Kokkos::deep_copy(dst.B(), src.B());
      Kokkos::deep_copy(dst.C(), src.C());
    }

    template<typename ViewType>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if< std::is_same<typename ViewType::value_type,scalar_type>::value, scalar_type&>::type
    tdiag_val(const ViewType &A, 
              const ordinal_type &t, 
              const ordinal_type &i, 
              const ordinal_type &ii, 
              const ordinal_type &jj) { 
      return A(t, i, ii, jj);
    }

    template<typename ViewType>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if< !std::is_same<typename ViewType::value_type,scalar_type>::value, scalar_type&>::type
    tdiag_val(const ViewType &A,
              const ordinal_type &t, 
              const ordinal_type &i, 
              const ordinal_type &ii, 
              const ordinal_type &jj) { 
      typedef typename ViewType::value_type value_type;
      return A(t/value_type::vector_length, i, ii, jj)[t%value_type::vector_length];
    }

    template <typename ExeSpace, typename ValueType, typename ArrayLayout>
    class PartitionedBlockMultiVector {
    public:
      typedef ExeSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef Kokkos::View<value_type****,array_layout,exec_space> value_array_type;

    private:
      value_array_type _values;

    public:
      PartitionedBlockMultiVector(const ordinal_type nparts,
                                  const ordinal_type nvectors,
                                  const ordinal_type nrows,
                                  const ordinal_type blocksize)
        : _values("BlockMultiVector::_values", nparts, nvectors, nrows, blocksize) {}

      PartitionedBlockMultiVector(const value_array_type values) 
        : _values(values) {}
      
      ordinal_type NumPartitions() const { return _values.extent(0); }
      ordinal_type NumVectors() const { return _values.extent(1); }      
      ordinal_type NumRows() const { return _values.extent(2); }
      ordinal_type BlockSize() const { return _values.extent(3); }

      value_array_type Values() const { return _values; }
    };

    template<typename ExecSpace, typename ValueType, typename ArrayLayout>
    PartitionedBlockMultiVector<ExecSpace,ValueType,ArrayLayout> 
    create_partitioned_block_multi_vector(const ordinal_type nparts,
                                          const ordinal_type nvectors,
                                          const ordinal_type nrows,
                                          const ordinal_type blocksize) {
      return PartitionedBlockMultiVector
        <ExecSpace,ValueType,ArrayLayout>(adjustDimension<ValueType>(nparts), 
                                          nvectors,
                                          nrows, 
                                          blocksize);
    }
    
    template<typename DstSpace, typename SrcSpace, typename ValueType, typename ArrayLayout>
    inline 
    PartitionedBlockMultiVector<DstSpace,ValueType,ArrayLayout>
    create_mirror(const PartitionedBlockMultiVector<SrcSpace,ValueType,ArrayLayout> src) {
      return PartitionedBlockMultiVector<DstSpace,ValueType,ArrayLayout>
        (Kokkos::create_mirror_view(typename DstSpace::memory_space(), src.Values()));
    }
    
    template<typename DstSpace, typename SrcSpace, typename ValueType, typename ArrayLayout>
    inline 
    void
    deep_copy(const PartitionedBlockMultiVector<DstSpace,ValueType,ArrayLayout> dst, 
              const PartitionedBlockMultiVector<SrcSpace,ValueType,ArrayLayout> src) {
      Kokkos::deep_copy(dst.Values(), src.Values());
    }
    
    template<typename ValueType, typename ArrayLayout>
    void fill_partitioned_block_multi_vector_host
    (PartitionedBlockMultiVector<Kokkos::DefaultHostExecutionSpace,ValueType,ArrayLayout> B,
     const ordinal_type ninj) {
      const ordinal_type 
        iend = ninj, // B.NumPartitions(),
        jend = B.NumVectors(),
        kend = B.NumRows(), 
        lend = B.BlockSize();
      
      auto B_val = B.Values();
      for (ordinal_type i=0;i<iend;++i) 
        for (ordinal_type j=0;j<jend;++j) 
          for (ordinal_type k=0;k<kend;++k) 
            for (ordinal_type l=0;l<lend;++l)
              tdiag_val(B_val, i,j,k,l) = static_cast<double>((i+j+k+l)%7) - 3;
    }
    

    template <typename ExecSpace, typename ArrayLayout>
    class BlockCrsMatrixVectorProductByRow {
    public:
      typedef BlockCrsMatrix<ExecSpace,ArrayLayout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockMultiVector<ExecSpace,ArrayLayout> block_multi_vector_type;
      
    private:
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;

      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      ConstUnmanagedViewType<typename block_multi_vector_type::value_array_type> _x;
      /**/ UnmanagedViewType<typename block_multi_vector_type::value_array_type> _y;

      ordinal_type _blocksize;

    public:
      // A thread maps to a point row of the matrix.
      // loop = blksize*m
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type idx) const {
        // index of blockrow and row in a block
        const ordinal_type i  = idx/_blocksize;
        const ordinal_type ii = idx%_blocksize;
        
        // loop over multivectors
        const ordinal_type jend = _y.extent(0);
        for (ordinal_type j=0;j<jend;++j) {
          scalar_type tmp = 0;
          
          // block row 
          const ordinal_type 
            cbegin = _rowptr(i), cend = _rowptr(i+1);

          for (ordinal_type c=cbegin;c<cend;++c) {
            const ordinal_type col = _colidx(c);
            for (ordinal_type jj=0;jj<_blocksize;++jj) 
              tmp += _A(col,ii,jj)*_x(j, col, jj);
          }
          _y(j, i, ii) = tmp;
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_multi_vector_type x, 
               const block_multi_vector_type y) {
        _rowptr = A.CrsGraph().rowptr;
        _colidx = A.CrsGraph().colidx;

        _blocksize = A.BlockSize();

        _A = A.Values();
        _x = x.Values();
        _y = y.Values();

        Kokkos::RangePolicy<ExecSpace> policy(0, _x.extent(1)*_blocksize);
        Kokkos::parallel_for(policy, *this);
      }
    };

    template <typename ExecSpace, typename ArrayLayout>
    class BlockCrsMatrixVectorProductByBlockRow {
    public:
      typedef BlockCrsMatrix<ExecSpace,ArrayLayout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockMultiVector<ExecSpace,ArrayLayout> block_multi_vector_type;

    private:
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;

      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      ConstUnmanagedViewType<typename block_multi_vector_type::value_array_type> _x;
      /**/ UnmanagedViewType<typename block_multi_vector_type::value_array_type> _y;

      ordinal_type _blocksize;

    public:
      
      // A thread maps to a row block of the matrix.
      // loop = m
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type i) const {
        // loop over multivector colums
        const ordinal_type jend = _y.extent(0);
        for (ordinal_type j=0;j<jend;++j) {
          // set zero
          for (ordinal_type ii=0;ii<_blocksize;++ii) 
            _y(j, i, ii) = 0;
          
          // block row 
          const ordinal_type 
            cbegin = _rowptr(i), cend = _rowptr(i+1);
          
          for (ordinal_type c=cbegin;c<cend;++c) {
            const ordinal_type col = _colidx(c);
            for (ordinal_type ii=0;ii<_blocksize;++ii) {
              scalar_type tmp = 0;
              for (ordinal_type jj=0;jj<_blocksize;++jj) 
                tmp += _A(col,ii,jj)*_x(j, col, jj);
              _y(j, i, ii) += tmp;
            }
          }
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_multi_vector_type x, 
               const block_multi_vector_type y) {
        _rowptr = A.CrsGraph().rowptr;
        _colidx = A.CrsGraph().colidx;

        _blocksize = A.BlockSize();

        _A = A.Values();
        _x = x.Values();
        _y = y.Values();
        
        Kokkos::RangePolicy<ExecSpace> policy(0, _x.extent(1));
        Kokkos::parallel_for(policy, *this);
      }
    };

    template <typename ExecSpace, typename ValueType, typename ArrayLayout>
    class ExtractBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef StructuredBlock structured_block_mesh_type;
      typedef BlockCrsMatrix<exec_space,array_layout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;

    private:
      structured_block_mesh_type _mesh;
      ordinal_type _blocksize;
      
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::row_idx_type> _rowidx;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;
      
      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      /**/ UnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;
      
    public:
      ExtractBlockTridiagMatrices(const structured_block_mesh_type mesh) 
        : _mesh(mesh) {}

      template<typename TViewType,
               typename AViewType>
      KOKKOS_INLINE_FUNCTION
      void 
      elementwise_copy(const TViewType &T,
                       const AViewType &A,
                       const ordinal_type ij, 
                       const ordinal_type k,
                       const ordinal_type c,
                       const ordinal_type blocksize) const {
        for (ordinal_type ii=0;ii<blocksize;++ii)
          for (ordinal_type jj=0;jj<blocksize;++jj) 
            tdiag_val(T, ij, k, ii, jj) = A(c, ii, jj);
      }
      
      // A thread maps nonzero blocks
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type c) const {
        const ordinal_type row = _rowidx[c], col = _colidx[c];

        ordinal_type ri, rj, rk, ci, cj, ck;
        _mesh.id2ijk(row, ri, rj, rk);
        _mesh.id2ijk(col, ci, cj, ck);
  
        if (ri == ci && rj == cj) {
          const ordinal_type ij = _mesh.ij2id(ri, rj);
          // consider connectivity to k-direction
          switch (rk - ck) {
          case  1: elementwise_copy(_TC, _A, ij, ck, c, _blocksize); break;
          case  0: elementwise_copy(_TA, _A, ij, rk, c, _blocksize); break;
          case -1: elementwise_copy(_TB, _A, ij, rk, c, _blocksize); break;
          }
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_tridiag_matrices_type T) { 
        _rowptr = A.CrsGraph().rowptr;
        _rowidx = A.CrsGraph().rowidx;
        _colidx = A.CrsGraph().colidx;

        _A = A.Values();

        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();

        _blocksize = A.BlockSize();
        Kokkos::RangePolicy<ExecSpace> policy(0, _A.extent(0));
        Kokkos::parallel_for(policy, *this);
      }

      template<typename TViewType,
               typename AViewType>
      bool elementwise_check(const TViewType &T,
                             const AViewType &A,
                             const ordinal_type ij, 
                             const ordinal_type k,
                             const ordinal_type c,
                             const ordinal_type blocksize) const {
        const auto eps = 1e2*std::numeric_limits<scalar_type>::epsilon();
        for (ordinal_type ii=0;ii<blocksize;++ii)
          for (ordinal_type jj=0;jj<blocksize;++jj) 
            if ( std::abs(tdiag_val(T, ij, k, ii, jj) - A(c, ii, jj)) >= eps ) return false; 
        return true;
      }
      
      bool check() const {        
        auto rowptr = Kokkos::create_mirror_view(_rowptr); Kokkos::deep_copy(rowptr, _rowptr);
        auto colidx = Kokkos::create_mirror_view(_colidx); Kokkos::deep_copy(colidx, _colidx);
        auto TA     = Kokkos::create_mirror_view(_TA);     Kokkos::deep_copy(TA, _TA);
        auto TB     = Kokkos::create_mirror_view(_TB);     Kokkos::deep_copy(TB, _TB);
        auto TC     = Kokkos::create_mirror_view(_TC);     Kokkos::deep_copy(TC, _TC);
        auto A      = Kokkos::create_mirror_view(_A);      Kokkos::deep_copy(A, _A);

        const ordinal_type 
          ijend = adjustDimension<value_type>(_mesh.ni*_mesh.nj),
          kend = _mesh.nk;

        assert(ijend == ordinal_type(TA.extent(0))); assert((kend - 0) == ordinal_type(TA.extent(1))); 
        assert(ijend == ordinal_type(TB.extent(0))); assert((kend - 1) == ordinal_type(TB.extent(1)));
        assert(ijend == ordinal_type(TC.extent(0))); assert((kend - 1) == ordinal_type(TC.extent(1)));

        for (ordinal_type ij=0;ij<ijend;++ij) {
          ordinal_type i, j;
          _mesh.id2ij(ij, i, j);

          for (ordinal_type k=0;k<kend;++k) {
            const ordinal_type row = _mesh.ijk2id(i, j, k),
              idx_begin = rowptr[row], 
              idx_end = rowptr[row+1];

            // check
            bool found[3] = {}, same[3] = {};
            for (ordinal_type idx=idx_begin;idx<idx_end;++idx) {
              switch (row - colidx[idx]) {
              case  1: same[2] = elementwise_check(TC, A, ij, k-1, idx, _blocksize); found[2] = true; break;
              case  0: same[0] = elementwise_check(TA, A, ij, k,   idx, _blocksize); found[0] = true; break;
              case -1: same[1] = elementwise_check(TB, A, ij, k,   idx, _blocksize); found[1] = true; break;
              }
            }
            if      (k == 0)         assert(found[0] & same[0] && found[1] & same[1]); 
            else if (k == (kend-1))  assert(found[0] & same[0] && found[2] & same[2]); 
            else                     assert(found[0] & same[0] && found[1] & same[1] && found[2] & same[2]); 
          }
        }            
        return true;
      }
    };

    inline bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
      return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
              a == std::string("-") + std::string(b1));
    }

    // Command-line argument parser and holder.
    template<typename ExecSpace>
    struct Input {
      bool quiet, check;
      ordinal_type ni, nj, nk;
      ordinal_type bs; // block size
      ordinal_type nrhs; // #vectors in multivector
      ordinal_type opf, ops;
      StencilShape::Enum stencil_shape;
      
      Input (int argc, char** argv) {
        quiet = false;
        check = false;
        ni = nj = nk = 10;
        bs = 5;
        nrhs = 1;
        if (std::is_same<typename ExecSpace::memory_space,Kokkos::HostSpace>::value) {
          opf = 0; ops = 0; // range policy default
        } else {
          opf = 1; ops = 1; // team is default
        }
        stencil_shape = StencilShape::cross;

        for (ordinal_type i=1;i<argc;++i) {
          const std::string& token = argv[i];
          if (eq(token, "-nijk")) ni = nj = nk = std::atoi(argv[++i]);
          else if (eq(token, "-ni")) ni = std::atoi(argv[++i]);
          else if (eq(token, "-nj")) nj = std::atoi(argv[++i]);
          else if (eq(token, "-nk")) nk = std::atoi(argv[++i]);
          else if (eq(token, "-bs")) bs = std::atoi(argv[++i]);
          else if (eq(token, "-nrhs")) nrhs = std::atoi(argv[++i]);
          else if (eq(token, "-opf")) opf = std::atoi(argv[++i]);
          else if (eq(token, "-ops")) ops = std::atoi(argv[++i]);
          else if (eq(token, "-c", "-check")) check = true;
        }
        if (nk <= 1)
          throw std::runtime_error("k dimension is <= 1; must be >= 2.");
        if ( ! quiet) print(std::cout);
      }
      
      void print (std::ostream& os) const {
        os << "<I> ni " << ni << " nj " << nj << " nk " << nk
           << " bs " << bs
           << " nrhs " << nrhs
           << " opf " << opf << " ops " << ops
           << " sc " << stencil_shape << "\n";
      }
    };








  }
}
