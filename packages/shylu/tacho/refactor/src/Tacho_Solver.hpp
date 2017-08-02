#ifndef __TACHO_SOLVER_HPP__
#define __TACHO_SOLVER_HPP__

/// \file Tacho_Solver.hpp
/// \brief solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho.hpp"

namespace Tacho {

  using namespace Experimental;

  template<typename ValueType,
           typename DeviceSpaceType>
  struct Solver {
  public:
    typedef ValueType value_type;
    typedef DeviceSpaceType device_exec_space;

    typedef Graph crs_graph_type; // host only

    typedef SupernodeInfo<value_type,device_exec_space> supernode_info_type;
    typedef typename supernode_info_type::crs_matrix_type crs_matrix_type;

    typedef typename supernode_info_type::size_type_array size_type_array;
    typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
    typedef typename supernode_info_type::value_type_array value_type_array;
    typedef typename supernode_info_type::value_type_matrix value_type_matrix;

#if   defined(HAVE_SHYLUTACHO_METIS)
    typedef GraphTools_Metis graph_tools_type;
#elif defined(HAVE_SHYLUTACHO_SCOTCH)
    typedef GraphTools_Scotch graph_tools_type;
#else
    typedef GraphTools_CAMD graph_tools_type;
#endif

    typedef SymbolicTools symbolic_tools_type;
    typedef NumericTools<value_type,device_exec_space> numeric_tools_type;

  private:
    // problem
    ordinal_type _m;
    size_type _nnz;

    size_type_array _ap;
    ordinal_type_array _aj;

    // permutation
    ordinal_type_array _perm, _peri;

    // tools
    graph_tools_type _G;
    symbolic_tools_type _S;
    numeric_tools_type _N;

    value_type_matrix _U;

    // options
    ordinal_type _verbose;              // print
    ordinal_type _small_problem_thres;  // smaller than this, use lapack
    ordinal_type _serial_thres_size;    // serialization threshold size
    ordinal_type _mb;                   // block size for byblocks algorithms

    // parallelism and memory constraint is made via this parameter
    ordinal_type _max_num_superblocks;  // # of superblocks in the memoyrpool
    
  public:
    Solver()
      : _m(0), _nnz(0),
        _ap(), _aj(),
        _perm(), _peri(),
        _verbose(0),
        _small_problem_thres(1024),
        _serial_thres_size(-1),
        _mb(-1) {}

    Solver(const Solver &b) = default;

    void setVerbose(const ordinal_type verbose = 1) {
      _verbose = verbose;
    }
    void setSerialThresholdSize(const ordinal_type serial_thres_size = -1) {
      _serial_thres_size = serial_thres_size;
    }
    void setBlocksize(const ordinal_type mb = -1) {
      _mb = mb;
    }
    void setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks = -1) {
      _max_num_superblocks = max_num_superblocks;
    }

    /// need to figure out how HTS options are set here

    int analyze(const ordinal_type m,
                const size_type_array &ap,
                const ordinal_type_array &aj) {
      // assign graph structure
      _m = m;
      _nnz = ap(m);
      _ap = ap;
      _aj = aj;

      if (_verbose) {
        printf("TachoSolver: Analyze\n");
        printf("====================\n");
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
        Graph graph(_m, _nnz, _ap, _aj);
        _G = graph_tools_type(graph);
        _G.reorder(_verbose);

        _S = symbolic_tools_type(_m, _ap, _aj, _G.PermVector(), _G.InvPermVector());
        _S.symbolicFactorize(_verbose);
      }
      
      return 0;
    }

    int factorize(const value_type_array &ax) {
      if (_verbose) {
        printf("TachoSolver: Factorize\n");
        printf("======================\n");
      }

      if (_m < _small_problem_thres) {
        Kokkos::Impl::Timer timer;
        _U = value_type_matrix("U", _m, _m);
        timer.reset();
        for (ordinal_type i=0;i<_m;++i) {
          const size_type jbeg = _ap(i), jend = _ap(i+1);
          for (size_type j=jbeg;j<jend;++j) {
            const ordinal_type col = _aj(j);
            if (i <= col) _U(i, col) = ax(j);
          }
        }
        const double t_copy = timer.seconds();

        timer.reset();
        const ordinal_type sched = 0, member = 0;
        Chol<Uplo::Upper,Algo::External>
          ::invoke(sched, member, _U);
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
                                _G.PermVector(), _G.InvPermVector(),
                                _S.NumSupernodes(), _S.Supernodes(),
                                _S.gidSuperPanelPtr(), _S.gidSuperPanelColIdx(),
                                _S.sidSuperPanelPtr(), _S.sidSuperPanelColIdx(), _S.blkSuperPanelColIdx(),
                                _S.SupernodesTreeParent(), _S.SupernodesTreePtr(), _S.SupernodesTreeChildren(),
                                _S.SupernodesTreeRoots());
        
        if (_serial_thres_size < 0) { // set default values
          _serial_thres_size = 64;
        }
        _N.setSerialThresholdSize(_serial_thres_size);

        if (_max_num_superblocks < 0) { // set default values
          _max_num_superblocks = 4;
        }
        _N.setMaxNumberOfSuperblocks(_max_num_superblocks);

        const ordinal_type nthreads = device_exec_space::thread_pool_size(0);
        if (nthreads == 1) {
          _N.factorizeCholesky_Serial(ax, _verbose);
        } else {
          if (_mb < 0) {
            const ordinal_type max_dense_size = max(_N.getMaxSupernodeSize(),_N.getMaxSchurSize());
            if      (max_dense_size < 256)  _mb =  -1;
            else if (max_dense_size < 512)  _mb =  96;
            else if (max_dense_size < 1024) _mb = 128;
            else                            _mb = 256;
          }
          if (_mb > 0)
            _N.factorizeCholesky_ParallelByBlocks(ax, _mb, _verbose);
          else
            _N.factorizeCholesky_Parallel(ax, _verbose);
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
        const ordinal_type sched = 0, member = 0;
        Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
          ::invoke(sched, member, Diag::NonUnit(), 1.0, _U, x);
        Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
          ::invoke(sched, member, Diag::NonUnit(), 1.0, _U, x);
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
        const ordinal_type nthreads = device_exec_space::thread_pool_size(0);
        if (nthreads == 1)
          _N.solveCholesky_Serial(x, b, t, _verbose);
        else
          _N.solveCholesky_Parallel(x, b, t, _verbose);
      }
      return 0;
    }

    double computeRelativeResidual(const value_type_array &ax,
                                   const value_type_matrix &x,
                                   const value_type_matrix &b) {
      crs_matrix_type A;
      A.setExternalMatrix(_m, _m, _nnz, _ap, _aj, ax);

      return _N.computeRelativeResidual(A, x, b);
    }
  };

}

#endif
