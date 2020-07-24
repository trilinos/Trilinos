#ifndef __TACHO_SOLVER_HPP__
#define __TACHO_SOLVER_HPP__

/// \file Tacho_Solver.hpp
/// \brief solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho.hpp"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

namespace Tacho {

  /// forward decl
  class Graph;
#if defined(TACHO_HAVE_METIS)
  class GraphTools_Metis;
#endif
  class GraphTools_CAMD;
  
  class SymbolicTools;
  template<typename ValueType, typename DeviceType> class CrsMatrixBase;
  template<typename ValueType, typename SchedulerType> class NumericTools;
  template<typename ValueType, typename SchedulerType, int Variant> class LevelSetTools;
  
  ///
  /// Tacho Solver interface
  ///
  template<typename ValueType,
           typename SchedulerType>
  struct Solver {
  public:
    typedef ValueType value_type;
    typedef SchedulerType scheduler_type;

    typedef typename UseThisDevice<typename scheduler_type::execution_space>::device_type device_type;
    typedef typename device_type::execution_space exec_space;
    typedef typename device_type::memory_space exec_memory_space;

    typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;
    typedef typename host_device_type::execution_space host_space;
    typedef typename host_device_type::memory_space host_memory_space;

    typedef Kokkos::View<size_type*,device_type> size_type_array;
    typedef Kokkos::View<ordinal_type*,device_type> ordinal_type_array;
    typedef Kokkos::View<value_type*,device_type> value_type_array;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> value_type_matrix; 

    typedef Kokkos::View<size_type*,host_device_type> size_type_array_host;
    typedef Kokkos::View<ordinal_type*,host_device_type> ordinal_type_array_host;
    typedef Kokkos::View<value_type*,host_device_type> value_type_array_host;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,host_device_type> value_type_matrix_host; 

    typedef CrsMatrixBase<value_type,device_type> crs_matrix_type; 
    typedef CrsMatrixBase<value_type,host_device_type> crs_matrix_type_host; 

#if defined(TACHO_HAVE_METIS)
    typedef GraphTools_Metis graph_tools_type;
#else
    typedef GraphTools_CAMD graph_tools_type;
#endif

    typedef SymbolicTools symbolic_tools_type;
    typedef NumericTools<value_type,scheduler_type> numeric_tools_type;
    typedef LevelSetTools<value_type,scheduler_type,0> levelset_tools_var0_type;
    typedef LevelSetTools<value_type,scheduler_type,1> levelset_tools_var1_type;
    typedef LevelSetTools<value_type,scheduler_type,2> levelset_tools_var2_type;

  public:
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
    ordinal_type_array_host _stree_level, _stree_roots;

    // ** numeric factorization output
    numeric_tools_type *_N;

    // ** level set interface
    levelset_tools_var0_type *_L0;
    levelset_tools_var1_type *_L1;
    levelset_tools_var2_type *_L2;

    // small dense matrix
    value_type_matrix_host _U;

    // ** options
    ordinal_type _verbose;              // print
    ordinal_type _small_problem_thres;  // smaller than this, use lapack
    
    // ** tasking options
    ordinal_type _serial_thres_size;    // serialization threshold size    
    ordinal_type _mb;                   // block size for byblocks algorithms
    ordinal_type _nb;                   // panel size for panel algorithms
    ordinal_type _front_update_mode;    // front update mode 0 - lock, 1 - atomic

    // ** levelset options
    bool         _levelset;             // use level set code instead of tasking
    ordinal_type _device_level_cut;     // above this level, matrices are computed on device
    ordinal_type _device_factor_thres;  // bigger than this threshold, device function is used
    ordinal_type _device_solve_thres;   // bigger than this threshold, device function is used
    ordinal_type _variant;              // algorithmic variant in levelset 0: naive, 1: invert diagonals
    ordinal_type _nstreams;             // on cuda, multi streams are used

    // parallelism and memory constraint is made via this parameter
    ordinal_type _max_num_superblocks;  // # of superblocks in the memoyrpool

  public:
    Solver();
    /// delete copy constructor and assignment operator
    /// sharing numeric tools for different inputs does not make sense
    Solver(const Solver &) = delete;
    Solver& operator=(const Solver &) = delete;

    ///
    /// common options
    ///
    void setVerbose(const ordinal_type verbose = 1);
    void setSmallProblemThresholdsize(const ordinal_type small_problem_thres = 1024);
    void setTransposeSolve(const bool transpose);
    void setMatrixType(const int symmetric, // 0 - unsymmetric, 1 - structure sym, 2 - symmetric, 3 - hermitian
                       const bool is_positive_definite);

    ///
    /// tasking options
    ///
    void setSerialThresholdsize(const ordinal_type serial_thres_size = -1);
    void setBlocksize(const ordinal_type mb = -1);
    void setPanelsize(const ordinal_type nb = -1);
    void setFrontUpdateMode(const ordinal_type front_update_mode = 1);
    void setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks = -1);

    ///
    /// Level set tools options
    ///
    void setLevelSetScheduling(const bool levelset);
    void setLevelSetOptionDeviceLevelCut(const ordinal_type device_level_cut);
    void setLevelSetOptionDeviceFunctionThreshold(const ordinal_type device_factor_thres,
                                                  const ordinal_type device_solve_thres);
    void setLevelSetOptionNumStreams(const ordinal_type nstreams);
    void setLevelSetOptionAlgorithmVariant(const ordinal_type variant);

    ///
    /// get interface
    ///
    ordinal_type       getNumSupernodes() const;
    ordinal_type_array getSupernodes() const;
    ordinal_type_array getPermutationVector() const;
    ordinal_type_array getInversePermutationVector() const;

    // internal only
    int analyze();
    int analyze(const ordinal_type m,
                const size_type_array_host &ap,
                const ordinal_type_array_host &aj,
                const ordinal_type_array_host &perm,
                const ordinal_type_array_host &peri);

    template<typename arg_size_type_array,
             typename arg_ordinal_type_array>
    int analyze(const ordinal_type m,
                const arg_size_type_array &ap,
                const arg_ordinal_type_array &aj) {
      _m = m;

      _ap   = Kokkos::create_mirror_view(exec_memory_space(), ap); Kokkos::deep_copy(_ap, ap);
      _aj   = Kokkos::create_mirror_view(exec_memory_space(), aj); Kokkos::deep_copy(_aj, aj);

      _h_ap = Kokkos::create_mirror_view(host_memory_space(), ap); Kokkos::deep_copy(_h_ap, ap);
      _h_aj = Kokkos::create_mirror_view(host_memory_space(), aj); Kokkos::deep_copy(_h_aj, aj);

      _nnz = _h_ap(m);
      
      return analyze();
    }

    int initialize();
    int factorize(const value_type_array &ax);
    int solve(const value_type_matrix &x,
              const value_type_matrix &b,
              const value_type_matrix &t);

    double computeRelativeResidual(const value_type_array &ax,
                                   const value_type_matrix &x,
                                   const value_type_matrix &b);
    
    int exportFactorsToCrsMatrix(crs_matrix_type &A);
    int release();

  };

}

//#include "Tacho_Solver_Impl.hpp"

#endif
