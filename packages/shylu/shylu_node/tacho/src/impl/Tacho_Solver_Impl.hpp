#ifndef __TACHO_SOLVER_IMPL_HPP__
#define __TACHO_SOLVER_IMPL_HPP__

/// \file Tacho_Solver_Impl.hpp
/// \brief solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Internal.hpp"
#include "Tacho_Solver.hpp"

namespace Tacho {
  
  template<typename VT, typename ST>
  Solver<VT,ST>
  ::Solver()
    : _m(0), _nnz(0),
      _ap(), _h_ap(),_aj(), _h_aj(),
      _perm(), _h_perm(), _peri(), _h_peri(),
      _N(nullptr),
      _L0(nullptr),
      _L1(nullptr),
      _L2(nullptr),
      _verbose(0),
      _small_problem_thres(1024),
      _serial_thres_size(-1),
      _mb(-1),
      _nb(-1),
      _front_update_mode(-1),
      _levelset(0),
      _device_level_cut(0),
      _device_factor_thres(64),
      _device_solve_thres(128),
      _variant(2),
      _nstreams(16),
      _max_num_superblocks(-1) {}

  /// deleted
  // template<typename VT, typename ST>
  // Solver<VT,ST>
  // ::Solver(const Solver &b) = default;
  
  ///
  /// common options
  ///
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setVerbose(const ordinal_type verbose) {
    _verbose = verbose;
  }
  
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>  
  ::setSmallProblemThresholdsize(const ordinal_type small_problem_thres) {
    _small_problem_thres = small_problem_thres;
  }
  
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setTransposeSolve(const bool transpose) {
    _transpose = transpose; // this option is not used yet
  }
  
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setMatrixType(const int symmetric, // 0 - unsymmetric, 1 - structure sym, 2 - symmetric, 3 - hermitian
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

  ///
  /// tasking options
  ///
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setSerialThresholdsize(const ordinal_type serial_thres_size) {
    _serial_thres_size = serial_thres_size;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setBlocksize(const ordinal_type mb) {
    _mb = mb;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setPanelsize(const ordinal_type nb) {
    _nb = nb;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setFrontUpdateMode(const ordinal_type front_update_mode) {
    _front_update_mode = front_update_mode;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks) {
    _max_num_superblocks = max_num_superblocks;
  }

  ///
  /// Level set tools options
  ///
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setLevelSetScheduling(const bool levelset) {
    _levelset = levelset; 
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setLevelSetOptionDeviceLevelCut(const ordinal_type device_level_cut) {
    _device_level_cut = device_level_cut;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setLevelSetOptionDeviceFunctionThreshold(const ordinal_type device_factor_thres,
                                             const ordinal_type device_solve_thres) {
    _device_factor_thres = device_factor_thres;
    _device_solve_thres = device_solve_thres;
  }

  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setLevelSetOptionAlgorithmVariant(const ordinal_type variant) {
    if (variant > 2 || variant < 0) {
      std::logic_error("levelset algorithm variants range from 0 to 2");
    }
    _variant = variant;
  }
  
  template<typename VT, typename ST>  
  void 
  Solver<VT,ST>
  ::setLevelSetOptionNumStreams(const ordinal_type nstreams) {
    _nstreams = nstreams;
  }

  ///
  /// get interface
  ///
  template<typename VT, typename ST>  
  ordinal_type
  Solver<VT,ST>
  ::getNumSupernodes() const { 
    return _nsupernodes;
  } 

  template<typename VT, typename ST>  
  typename Solver<VT,ST>::ordinal_type_array
  Solver<VT,ST>
  ::getSupernodes() const { 
    return _supernodes;
  } 

  template<typename VT, typename ST>  
  typename Solver<VT,ST>::ordinal_type_array
  Solver<VT,ST>
  ::getPermutationVector() const { 
    return _perm;
  } 

  template<typename VT, typename ST>  
  typename Solver<VT,ST>::ordinal_type_array
  Solver<VT,ST>
  ::getInversePermutationVector() const { 
    return _peri; 
  }
  

  // internal only
  template<typename VT, typename ST>  
  int
  Solver<VT,ST>  
  ::analyze() {
    Graph graph(_m, _nnz, _h_ap, _h_aj);
    graph_tools_type G(graph);
    G.reorder(_verbose);

    _h_perm = G.PermVector(); 
    _h_peri = G.InvPermVector();

    // invoke a private function
    return analyze(_m, _h_ap, _h_aj, _h_perm, _h_peri);
  }

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>  
  ::analyze(const ordinal_type m,
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
      // for small problems, use lapack and no analysis and no permutation
      _perm = Kokkos::create_mirror_view(exec_memory_space(), _h_perm); 
      _peri = Kokkos::create_mirror_view(exec_memory_space(), _h_peri); 

      Kokkos::deep_copy(_perm, _h_perm);
      Kokkos::deep_copy(_peri, _h_peri);
        
      if (_verbose) {
        printf("  Linear system A\n");
        printf("             number of equations:                             %10d\n", _m);
        printf("\n");
        printf("  A is a small problem ( < %d ) and LAPACK is used\n", _small_problem_thres);
        printf("\n");
      }
      return 0;
    } else {
      symbolic_tools_type S(_m, _h_ap, _h_aj, _h_perm, _h_peri);
      S.symbolicFactorize(_verbose);

      _nsupernodes = S.NumSupernodes();
      _stree_level = S.SupernodesTreeLevel();
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

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>  
  ::initialize() {
    if (_verbose) {
      printf("TachoSolver: Initialize\n");
      printf("=======================\n");
    }

    ///
    /// initialize numeric tools
    ///
    if (_m < _small_problem_thres) {
      //_U = value_type_matrix_host("U", _m, _m);
    } else {
      if (_N == nullptr) 
        _N = (numeric_tools_type*) ::operator new (sizeof(numeric_tools_type));
      
      new (_N) numeric_tools_type(_m, _ap, _aj,
                                  _perm, _peri,
                                  _nsupernodes, _supernodes,
                                  _gid_super_panel_ptr, _gid_super_panel_colidx,
                                  _sid_super_panel_ptr, _sid_super_panel_colidx, _blk_super_panel_colidx,
                                  _stree_parent, _stree_ptr, _stree_children, 
                                  _stree_level, _stree_roots);
      
      if (_serial_thres_size < 0) { // set default values
        _serial_thres_size = 64;
      }
      _N->setSerialThresholdSize(_serial_thres_size);
      
      if (_max_num_superblocks < 0) { // set default values
        _max_num_superblocks = 16;
      }
      _N->setMaxNumberOfSuperblocks(_max_num_superblocks);
      
      if (_front_update_mode < 0) { // set default values
        _front_update_mode = 1; // atomic is default
      }
      _N->setFrontUpdateMode(_front_update_mode);
      _N->printMemoryStat(_verbose);

      ///
      /// initialize levelset tools
      ///
      if (_levelset) {
        if        (_variant == 0) {
          if (_L0 == nullptr) 
            _L0 = (levelset_tools_var0_type*) ::operator new (sizeof(levelset_tools_var0_type));
          new (_L0) levelset_tools_var0_type(*_N);
          _L0->initialize(_device_level_cut, _device_factor_thres, _device_solve_thres, _verbose);
          _L0->createStream(_nstreams);
        } else if (_variant == 1) {
          if (_L1 == nullptr) 
            _L1 = (levelset_tools_var1_type*) ::operator new (sizeof(levelset_tools_var1_type));
          new (_L1) levelset_tools_var1_type(*_N);
          _L1->initialize(_device_level_cut, _device_factor_thres, _device_solve_thres, _verbose);
          _L1->createStream(_nstreams);
        } else if (_variant == 2) {
          if (_L2 == nullptr) 
            _L2 = (levelset_tools_var2_type*) ::operator new (sizeof(levelset_tools_var2_type));
          new (_L2) levelset_tools_var1_type(*_N);
          _L2->initialize(_device_level_cut, _device_factor_thres, _device_solve_thres, _verbose);
          _L2->createStream(_nstreams);
        }
      }
    }
    return 0;
  }

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>      
  ::factorize(const value_type_array &ax) {
    if (_verbose) {
      printf("TachoSolver: Factorize\n");
      printf("======================\n");
    }
    if (_m < _small_problem_thres) {
      Kokkos::Impl::Timer timer;

      timer.reset();
      _U = value_type_matrix_host("U", _m, _m);
      auto h_ax = Kokkos::create_mirror_view(host_memory_space(), ax); Kokkos::deep_copy(h_ax, ax);
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

#if !defined (KOKKOS_ENABLE_CUDA)
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      const ordinal_type nthreads = host_space::thread_pool_size(0);
#else
      const ordinal_type nthreads = host_space::impl_thread_pool_size(0);
#endif
#endif
      if (_levelset) {
        if      (_variant == 0) _L0->factorizeCholesky(ax, _verbose);
        else if (_variant == 1) _L1->factorizeCholesky(ax, _verbose);
        else if (_variant == 2) _L2->factorizeCholesky(ax, _verbose);
      } 
#if !defined (KOKKOS_ENABLE_CUDA)
      else if (nthreads == 1) {
        if (_nb < 0) 
          _N->factorizeCholesky_Serial(ax, _verbose);
        else 
          _N->factorizeCholesky_SerialPanel(ax, _nb, _verbose);
      }
#endif
      else {
        const ordinal_type max_dense_size = max(_N->getMaxSupernodeSize(),_N->getMaxSchurSize());
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
            _N->factorizeCholesky_ParallelByBlocks(ax, _mb, _verbose);
          else
            _N->factorizeCholesky_Parallel(ax, _verbose);
        else 
          if (_mb > 0) 
            _N->factorizeCholesky_ParallelByBlocksPanel(ax, _mb, _nb, _verbose);
          else
            _N->factorizeCholesky_ParallelPanel(ax, _nb, _verbose);
      }
    }
    return 0;
  }

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>      
  ::solve(const value_type_matrix &x,
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
      auto h_x = Kokkos::create_mirror_view(host_memory_space(), x); 
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
      if (_levelset) {
        if      (_variant == 0) _L0->solveCholesky(x, b, tt, _verbose);
        else if (_variant == 1) _L1->solveCholesky(x, b, tt, _verbose);
        else if (_variant == 2) _L2->solveCholesky(x, b, tt, _verbose);
      } 
#if !defined (KOKKOS_ENABLE_CUDA)
      else if (nthreads == 1) {
        _N->solveCholesky_Serial(x, b, tt, _verbose);
      }
#endif
      else {
        _N->solveCholesky_Parallel(x, b, tt, _verbose);
      }
    }
    return 0;
  }

  template<typename VT, typename ST>  
  double
  Solver<VT,ST>      
  ::computeRelativeResidual(const value_type_array &ax,
                            const value_type_matrix &x,
                            const value_type_matrix &b) {
    CrsMatrixBase<value_type,device_type> A; 
    A.setExternalMatrix(_m, _m, _nnz, _ap, _aj, ax);

    return Tacho::computeRelativeResidual(A, x, b);
  }

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>          
  ::exportFactorsToCrsMatrix(crs_matrix_type &A) {
    if (_m < _small_problem_thres) {
      typedef ArithTraits<value_type> ats;
      const typename ats::mag_type zero(0);

      /// count nonzero elements in dense U
      const ordinal_type m = _m;
      size_type_array_host h_ap("h_ap", m+1);
      for (ordinal_type i=0;i<m;++i) 
        for (ordinal_type j=0;j<m;++j) 
          h_ap(i+1) += (ats::abs(_U(i,j)) > zero);
        
      /// serial scan; this is a small problem
      h_ap(0) = 0;
      for (ordinal_type i=0;i<m;++i)
        h_ap(i+1) += h_ap(i);

      /// create a host crs matrix
      const ordinal_type nnz = h_ap(m);
      ordinal_type_array_host h_aj(do_not_initialize_tag("h_aj"), nnz);
      value_type_array_host h_ax(do_not_initialize_tag("h_ax"), nnz);

      for (ordinal_type i=0,k=0;i<m;++i) 
        for (ordinal_type j=i;j<m;++j)
          if (ats::abs(_U(i,j)) > zero) {
            h_aj(k) = j;
            h_ax(k) = _U(i,j);
            ++k;
          }

      crs_matrix_type_host h_A;
      h_A.setExternalMatrix(m, m, nnz, h_ap, h_aj, h_ax); 
      ///h_A.showMe(std::cout, true);
      A.clear();
      A.createConfTo(h_A);
      A.copy(h_A);
    } else {
      _N->exportFactorsToCrsMatrix(A, false);
    }
    return 0;
  } 

  template<typename VT, typename ST>  
  int
  Solver<VT,ST>
  ::release() {
    if (_verbose) {
      printf("TachoSolver: Release\n");
      printf("====================\n");
    }

    if (_levelset) {
      if (_variant == 0) {
        if (_L0 != nullptr)
          _L0->release(_verbose);
        delete _L0; _L0 = nullptr;
      } else if (_variant == 1) {
        if (_L1 != nullptr)
          _L1->release(_verbose);
        delete _L1; _L1 = nullptr;
      } else if (_variant == 2) {
        if (_L2 != nullptr)
          _L2->release(_verbose);
        delete _L2; _L2 = nullptr;
      }
    }
    
    {
      if (_N != nullptr) 
        _N->release(_verbose);
      delete _N; _N = nullptr;
    }
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
      
      _U = value_type_matrix_host();
        
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

}

#endif
