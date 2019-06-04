#ifndef __TACHO_SOLVER_HPP__
#define __TACHO_SOLVER_HPP__

/// \file Tacho_Solver.hpp
/// \brief solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho.hpp"

namespace Tacho {

  template<typename ValueType,
           typename ExecSpace>
  struct Solver {
  public:
    typedef ValueType value_type;
    typedef ExecSpace exec_space;
    typedef typename ExecSpace::memory_space exec_memory_space;
    typedef Kokkos::DefaultHostExecutionSpace host_space;

    typedef Kokkos::View<size_type*,exec_space> size_type_array;
    typedef Kokkos::View<ordinal_type*,exec_space> ordinal_type_array;
    typedef Kokkos::View<value_type*,exec_space> value_type_array;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,exec_space> value_type_matrix; 

    typedef Kokkos::View<size_type*,host_space> size_type_array_host;
    typedef Kokkos::View<ordinal_type*,host_space> ordinal_type_array_host;
    typedef Kokkos::View<value_type*,host_space> value_type_array_host;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space> value_type_matrix_host; 

#if defined(TACHO_HAVE_METIS)
    typedef GraphTools_Metis graph_tools_type;
#else
    typedef GraphTools_CAMD graph_tools_type;
#endif

    typedef SymbolicTools symbolic_tools_type;
    typedef NumericTools<value_type,exec_space> numeric_tools_type;

  private:
    enum : int { Cholesky = 1,
                 UtDU = 2,
                 SymLU = 3,
                 LU = 4 };

    // ** solver mode
    bool _transpose;
    ordinal_type _mode;

    // ** problem
    ordinal_type _m;
    size_type _nnz;

    size_type_array    _ap; size_type_array_host    _h_ap;
    ordinal_type_array _aj; ordinal_type_array_host _h_aj;

    ordinal_type_array _perm; ordinal_type_array_host _h_perm;
    ordinal_type_array _peri; ordinal_type_array_host _h_peri;

    // ** symbolic factorization output
    // supernodes output
    ordinal_type _nsupernodes;
    ordinal_type_array _supernodes;

    // dof mapping to sparse matrix
    size_type_array _gid_super_panel_ptr;
    ordinal_type_array _gid_super_panel_colidx;

    // supernode map and panel size configuration
    size_type_array _sid_super_panel_ptr;
    ordinal_type_array _sid_super_panel_colidx, _blk_super_panel_colidx;

    // supernode elimination tree (parent - children)
    size_type_array _stree_ptr;
    ordinal_type_array _stree_children;

    // supernode elimination tree (child - parent)
    ordinal_type_array _stree_parent;

    // roots of supernodes
    ordinal_type_array_host _stree_roots;

    // ** numeric factorization output
    numeric_tools_type _N;

    // small dense matrix
    value_type_matrix_host _U;

    // ** options
    ordinal_type _verbose;              // print
    ordinal_type _small_problem_thres;  // smaller than this, use lapack
    ordinal_type _serial_thres_size;    // serialization threshold size
    ordinal_type _mb;                   // block size for byblocks algorithms
    ordinal_type _nb;                   // panel size for panel algorithms
    ordinal_type _front_update_mode;    // front update mode 0 - lock, 1 - atomic

    // parallelism and memory constraint is made via this parameter
    ordinal_type _max_num_superblocks;  // # of superblocks in the memoyrpool
    
  public:
    Solver()
      : _m(0), _nnz(0),
        _ap(), _h_ap(),_aj(), _h_aj(),
        _perm(), _h_perm(), _peri(), _h_peri(),
        _verbose(0),
        _small_problem_thres(1024),
        _serial_thres_size(-1),
        _mb(-1),
        _nb(-1),
        _front_update_mode(-1) {}

    Solver(const Solver &b) = default;

    void setVerbose(const ordinal_type verbose = 1) {
      _verbose = verbose;
    }
    void setSmallProblemThresholdsize(const ordinal_type small_problem_thres = 1024) {
      _small_problem_thres = small_problem_thres;
    }
    void setSerialThresholdsize(const ordinal_type serial_thres_size = -1) {
      _serial_thres_size = serial_thres_size;
    }
    void setBlocksize(const ordinal_type mb = -1) {
      _mb = mb;
    }
    void setPanelsize(const ordinal_type nb = -1) {
      _nb = nb;
    }
    void setFrontUpdateMode(const ordinal_type front_update_mode = 1) {
      _front_update_mode = front_update_mode;
    }
    void setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks = -1) {
      _max_num_superblocks = max_num_superblocks;
    }
    void setTransposeSolve(const bool transpose) {
      _transpose = transpose; // this option is not used yet
    }
    void setMatrixType(const int symmetric, // 0 - unsymmetric, 1 - structure sym, 2 - symmetric, 3 - hermitian
                       const bool is_positive_definite) { 
      switch (symmetric) {
      case 0: { _mode = LU; break; }
      case 1: { _mode = SymLU; break; }
      case 2: {
        if (is_positive_definite) {
          if (std::is_same<value_type,double>::value ||
              std::is_same<value_type,float>::value) { // real symmetric posdef
            _mode = Cholesky;
          } else { // complex symmetric posdef - does not exist; should be hermitian if the matrix is posdef
            TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                     "symmetric positive definite with complex values are not right combination: try hermitian positive definite");
          }
        } else { // real or complex symmetric indef
          _mode = UtDU;
        }
        break;
      }
      case 3: {
        if (is_positive_definite) { // real or complex hermitian posdef
          _mode = Cholesky;
        } else {
          if (std::is_same<value_type,double>::value ||
              std::is_same<value_type,float>::value) { // real symmetric indef 
            _mode = UtDU;
          } else { // complex hermitian indef
            _mode = SymLU; 
          }
        }
        break;
      }
      default: {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "symmetric argument is wrong");
      }
      }
      TACHO_TEST_FOR_EXCEPTION(_mode != Cholesky, std::logic_error, "Cholesky is only supported now");
    }

    // internal only
    int analyze(const ordinal_type m,
                const size_type_array_host &ap,
                const ordinal_type_array_host &aj,
                const ordinal_type_array_host &perm,
                const ordinal_type_array_host &peri) {
      if (_verbose) {
        printf("TachoSolver: Analyze\n");
        printf("====================\n");
      }

      if (_m == 0) {
        _m = m;
        
        _ap   = Kokkos::create_mirror_view(exec_memory_space(), ap); Kokkos::deep_copy(_ap, ap);
        _aj   = Kokkos::create_mirror_view(exec_memory_space(), aj); Kokkos::deep_copy(_aj, aj);

        _h_ap = ap;
        _h_aj = aj;

        _nnz = _ap(m);

        _h_perm = perm;
        _h_peri = peri;
      }

      if (_m < _small_problem_thres) {
        // for small problems, use lapack and no analysis
        if (_verbose) {
          printf("  Linear system A\n");
          printf("             number of equations:                             %10d\n", _m);
          printf("\n");
          printf("  A is a small problem ( < %d ) and LAPACK is used\n", _small_problem_thres);
          printf("\n");
        }
        return 0;
      } else {
        SymbolicTools S(_m, _h_ap, _h_aj, _h_perm, _h_peri);
        S.symbolicFactorize(_verbose);

        _nsupernodes = S.NumSupernodes();
        _stree_roots = S.SupernodesTreeRoots();

        _supernodes             = Kokkos::create_mirror_view(exec_memory_space(), S.Supernodes());            
        _gid_super_panel_ptr    = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelPtr());      
        _gid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.gidSuperPanelColIdx());   
        _sid_super_panel_ptr    = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelPtr());      
        _sid_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.sidSuperPanelColIdx());   
        _blk_super_panel_colidx = Kokkos::create_mirror_view(exec_memory_space(), S.blkSuperPanelColIdx());   
        _stree_parent           = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeParent());  
        _stree_ptr              = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreePtr());     
        _stree_children         = Kokkos::create_mirror_view(exec_memory_space(), S.SupernodesTreeChildren());

        Kokkos::deep_copy(_supernodes             , S.Supernodes());
        Kokkos::deep_copy(_gid_super_panel_ptr    , S.gidSuperPanelPtr());
        Kokkos::deep_copy(_gid_super_panel_colidx , S.gidSuperPanelColIdx());
        Kokkos::deep_copy(_sid_super_panel_ptr    , S.sidSuperPanelPtr());
        Kokkos::deep_copy(_sid_super_panel_colidx , S.sidSuperPanelColIdx());
        Kokkos::deep_copy(_blk_super_panel_colidx , S.blkSuperPanelColIdx());
        Kokkos::deep_copy(_stree_parent           , S.SupernodesTreeParent());
        Kokkos::deep_copy(_stree_ptr              , S.SupernodesTreePtr());
        Kokkos::deep_copy(_stree_children         , S.SupernodesTreeChildren());

        // perm and peri is updated during symbolic factorization
        _perm = Kokkos::create_mirror_view(exec_memory_space(), _h_perm); 
        _peri = Kokkos::create_mirror_view(exec_memory_space(), _h_peri); 

        Kokkos::deep_copy(_perm, _h_perm);
        Kokkos::deep_copy(_peri, _h_peri);
      }
      return 0;
    }

    template<typename arg_size_type_array,
             typename arg_ordinal_type_array>
    int analyze(const ordinal_type m,
                const arg_size_type_array &ap,
                const arg_ordinal_type_array &aj) {
      _m = m;

      _ap   = Kokkos::create_mirror_view(exec_memory_space(), ap); Kokkos::deep_copy(_ap, ap);
      _aj   = Kokkos::create_mirror_view(exec_memory_space(), aj); Kokkos::deep_copy(_aj, aj);

      _h_ap = Kokkos::create_mirror_view(ap); Kokkos::deep_copy(_h_ap, ap);
      _h_aj = Kokkos::create_mirror_view(aj); Kokkos::deep_copy(_h_aj, aj);

      _nnz = _h_ap(m);

      Graph graph(_m, _nnz, _h_ap, _h_aj);
      graph_tools_type G(graph);
      G.reorder(_verbose);

      _h_perm = G.PermVector(); 
      _h_peri = G.InvPermVector();

      // invoke a private function
      return analyze(_m, _h_ap, _h_aj, _h_perm, _h_peri);
    }

    int factorize(const value_type_array &ax) {
      if (_verbose) {
        printf("TachoSolver: Factorize\n");
        printf("======================\n");
      }

      if (_m < _small_problem_thres) {
        Kokkos::Impl::Timer timer;

        timer.reset();
        _U = value_type_matrix_host("U", _m, _m);
        auto h_ax = Kokkos::create_mirror_view(ax); Kokkos::deep_copy(h_ax, ax);
        for (ordinal_type i=0;i<_m;++i) {
          const size_type jbeg = _h_ap(i), jend = _h_ap(i+1);
          for (size_type j=jbeg;j<jend;++j) {
            const ordinal_type col = _h_aj(j);
            if (i <= col) 
              _U(i, col) = h_ax(j);
          }
        }
        const double t_copy = timer.seconds();
        
        timer.reset();
        Chol<Uplo::Upper,Algo::External>::invoke(_U);
        const double t_factor = timer.seconds();

        if (_verbose) {
          printf("Summary: NumericTools (SmallDenseFactorization)\n");
          printf("===============================================\n");
          printf("  Time\n");
          printf("             time for copying A into U:                       %10.6f s\n", t_copy);
          printf("             time for numeric factorization:                  %10.6f s\n", t_factor);
          printf("             total time spent:                                %10.6f s\n", (t_copy+t_factor));
          printf("\n");
        }
      } else {
        _N = numeric_tools_type(_m, _ap, _aj,
                                _perm, _peri,
                                _nsupernodes, _supernodes,
                                _gid_super_panel_ptr, _gid_super_panel_colidx,
                                _sid_super_panel_ptr, _sid_super_panel_colidx, _blk_super_panel_colidx,
                                _stree_parent, _stree_ptr, _stree_children,
                                _stree_roots);
        
        if (_serial_thres_size < 0) { // set default values
          _serial_thres_size = 64;
        }
        _N.setSerialThresholdSize(_serial_thres_size);

        if (_max_num_superblocks < 0) { // set default values
          _max_num_superblocks = 16;
        }
        _N.setMaxNumberOfSuperblocks(_max_num_superblocks);

        if (_front_update_mode < 0) { // set default values
          _front_update_mode = 1; // atomic is default
        }
        _N.setFrontUpdateMode(_front_update_mode);

#if !defined (KOKKOS_ENABLE_CUDA)
  #ifdef KOKKOS_ENABLE_DEPRECATED_CODE
        const ordinal_type nthreads = host_space::thread_pool_size(0);
  #else
        const ordinal_type nthreads = host_space::impl_thread_pool_size(0);
  #endif
#endif
        
        if (false) {
          // do nothing
        } 
#if !defined (KOKKOS_ENABLE_CUDA)
        else if (nthreads == 1) {
          if (_nb < 0) 
            _N.factorizeCholesky_Serial(ax, _verbose);
          else 
            _N.factorizeCholesky_SerialPanel(ax, _nb, _verbose);
        }
#endif
        else {
          const ordinal_type max_dense_size = max(_N.getMaxSupernodeSize(),_N.getMaxSchurSize());
          if (std::is_same<exec_memory_space,Kokkos::HostSpace>::value) {
            if (_nb < 0) { 
              _nb = 64;
              // if      (max_dense_size < 256)  _nb =  -1;
              // else if (max_dense_size < 512)  _nb =  64;
              // else if (max_dense_size < 1024) _nb = 128;
              // else if (max_dense_size < 8192) _nb = 256;
              // else                            _nb = 256;
            }
            if (_mb < 0) {
              if      (max_dense_size < 256)  _mb =  -1;
              else if (max_dense_size < 512)  _mb =  64;
              else if (max_dense_size < 1024) _mb = 128;
              else if (max_dense_size < 8192) _mb = 256;
              else                            _mb = 256;
            }
          } else {
            if (_nb < 0) { 
              _nb = 40;
              // if      (max_dense_size < 256)  _nb =  -1;
              // else if (max_dense_size < 512)  _nb =  64;
              // else if (max_dense_size < 1024) _nb = 128;
              // else if (max_dense_size < 8192) _nb = 256;
              // else                            _nb = 256;
            }
            if (_mb < 0) {
              if      (max_dense_size < 256)  _mb =  -1;
              else if (max_dense_size < 512)  _mb =  80;
              else if (max_dense_size < 1024) _mb = 120;
              else if (max_dense_size < 8192) _mb = 160;
              else                            _mb = 160;
            }            
          }
          
          if (_nb <= 0) 
            if (_mb > 0)
              _N.factorizeCholesky_ParallelByBlocks(ax, _mb, _verbose);
            else
              _N.factorizeCholesky_Parallel(ax, _verbose);
          else 
            if (_mb > 0) 
              _N.factorizeCholesky_ParallelByBlocksPanel(ax, _mb, _nb, _verbose);
            else
              _N.factorizeCholesky_ParallelPanel(ax, _nb, _verbose);
        }
      }
      return 0;
    }

    int solve(const value_type_matrix &x,
              const value_type_matrix &b,
              const value_type_matrix &t) {
      if (_verbose) {
        printf("TachoSolver: Solve\n");
        printf("==================\n");
      }

      if (_m < _small_problem_thres) {
        Kokkos::Impl::Timer timer;

        timer.reset();
        Kokkos::deep_copy(x, b);
        const double t_copy = timer.seconds();

        timer.reset();
        auto h_x = Kokkos::create_mirror_view(x); 
        Kokkos::deep_copy(h_x, x);
        Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
          ::invoke(Diag::NonUnit(), 1.0, _U, h_x);
        Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
          ::invoke(Diag::NonUnit(), 1.0, _U, h_x);
        Kokkos::deep_copy(x, h_x);         
        const double t_solve = timer.seconds();

        if (_verbose) {
          printf("Summary: NumericTools (SmallDenseSolve)\n");
          printf("=======================================\n");
          printf("  Time\n");
          printf("             time for extra work e.g.,copy rhs:               %10.6f s\n", t_copy);
          printf("             time for numeric solve:                          %10.6f s\n", t_solve);
          printf("             total time spent:                                %10.6f s\n", (t_solve+t_copy));
          printf("\n");
        }
      } else {
#if !defined (KOKKOS_ENABLE_CUDA)
  #ifdef KOKKOS_ENABLE_DEPRECATED_CODE
        const ordinal_type nthreads = host_space::thread_pool_size(0);
  #else
        const ordinal_type nthreads = host_space::impl_thread_pool_size(0);
  #endif
#endif
        TACHO_TEST_FOR_EXCEPTION(t.extent(0) < x.extent(0) ||
                                 t.extent(1) < x.extent(1), std::logic_error, "Temporary rhs vector t is smaller than x");
        auto tt = Kokkos::subview(t, 
                                  Kokkos::pair<ordinal_type,ordinal_type>(0, x.extent(0)),
                                  Kokkos::pair<ordinal_type,ordinal_type>(0, x.extent(1)));
        if (false) {
        } 
#if !defined (KOKKOS_ENABLE_CUDA)
        else if (nthreads == 1) {
          _N.solveCholesky_Serial(x, b, tt, _verbose);
        }
#endif
        else {
          _N.solveCholesky_Parallel(x, b, tt, _verbose);
        }
      }
      return 0;
    }

    double computeRelativeResidual(const value_type_array &ax,
                                   const value_type_matrix &x,
                                   const value_type_matrix &b) {
      CrsMatrixBase<value_type,exec_space> A; 
      A.setExternalMatrix(_m, _m, _nnz, _ap, _aj, ax);

      return _N.computeRelativeResidual(A, x, b);
    }

    int release() {
      if (_verbose) {
        printf("TachoSolver: Release\n");
        printf("====================\n");
      }
      _N.release(_verbose);

      {
        _transpose = false;
        _mode = 0;
        
        _m = 0;
        _nnz = 0;
        
        _ap = size_type_array();      _h_ap = size_type_array_host();
        _aj = ordinal_type_array();   _h_aj = ordinal_type_array_host();
        
        _perm = ordinal_type_array(); _h_perm = ordinal_type_array_host();
        _peri = ordinal_type_array(); _h_peri = ordinal_type_array_host();
        
        _nsupernodes = 0;
        _supernodes = ordinal_type_array();
        
        _gid_super_panel_ptr = size_type_array();
        _gid_super_panel_colidx = ordinal_type_array();
        
        _sid_super_panel_ptr = size_type_array();
        
        _sid_super_panel_colidx = ordinal_type_array();  
        _blk_super_panel_colidx = ordinal_type_array();
        
        _stree_ptr = size_type_array();
        _stree_children = ordinal_type_array();
        
        _stree_parent = ordinal_type_array();
        _stree_roots = ordinal_type_array_host();
        
        _N = numeric_tools_type();
        _U =value_type_matrix_host();
        
        _verbose = 0;
        _small_problem_thres = 1024;
        _serial_thres_size = -1; 
        _mb = -1;                
        _nb = -1;                
        _front_update_mode = -1; 
        
        _max_num_superblocks = -1;
      }
      return 0;
    }

  };

}

#endif
