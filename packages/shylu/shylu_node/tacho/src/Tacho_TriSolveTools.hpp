#ifndef __TACHO_TRISOLVE_TOOLS_HPP__
#define __TACHO_TRISOLVE_TOOLS_HPP__

/// \file Tacho_TriSolveTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_External.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_External.hpp"

#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_TeamFunctor_SolveLowerChol.hpp"
#include "Tacho_TeamFunctor_SolveUpperChol.hpp"
#if defined (KOKKOS_ENABLE_CUDA)
#include "cublas_v2.h"
#endif
namespace Tacho {

  ///
  /// Here we do not use a scheduler but all derived types in supernodes 
  /// info are required scheduler  
  ///
  template<typename ValueType, typename SchedulerType>
  class TriSolveTools {
  public:
    typedef ValueType value_type;
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

    ///
    /// device space typedefs
    ///
    typedef SchedulerType scheduler_type;

    typedef typename UseThisDevice<typename scheduler_type::execution_space>::device_type device_type;
    typedef typename device_type::execution_space exec_space;
    typedef typename device_type::memory_space exec_memory_space;

    typedef SupernodeInfo<value_type,scheduler_type> supernode_info_type;

    typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
    typedef typename supernode_info_type::size_type_array size_type_array;
    typedef typename supernode_info_type::value_type_array value_type_array;

    typedef typename supernode_info_type::value_type_matrix value_type_matrix;
    typedef typename supernode_info_type::supernode_type_array supernode_type_array;
    typedef typename supernode_type_array::HostMirror supernode_type_array_host;

    typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;
    typedef typename host_device_type::execution_space host_space;
    typedef typename host_device_type::memory_space host_memory_space;

    typedef Kokkos::View<ordinal_type*,host_device_type> ordinal_type_array_host;
    typedef Kokkos::View<size_type*,host_device_type> size_type_array_host;

  private:

    ///
    /// supernode data structure memory "managed"
    /// this holds all necessary connectivity data
    ///

    // graph ordering input
    ordinal_type_array _perm, _peri;

    ///
    /// supernode info: supernode data structure with "unamanged" view
    /// this is passed into computation algorithm without reference counting
    ///
    supernode_info_type _info;
    ordinal_type_array_host _h_stree_level;

    // supernodes details       
    ordinal_type _nsupernodes;
    supernode_type_array_host _h_supernodes;

    // 0: device level function, 1: team policy, 2: team policy recursive
    ordinal_type _device_function_thres;//, _team_recursive_thres;
    ordinal_type _device_level_cut, _team_recursive_level_cut;
    ordinal_type_array_host _h_compute_mode;
    ordinal_type_array        _compute_mode;

    // level details on host
    ordinal_type _nlevel;
    size_type_array_host _h_level_ptr;
    ordinal_type_array_host _h_level_sids;

    // level sids on device
    ordinal_type_array _level_sids;

    // workspace
    ordinal_type _max_nrhs;
    size_type_array_host _h_buf_ptr;
    size_type_array _buf_ptr;
    value_type_array _buf;

    // cuda stream
#if defined(KOKKOS_ENABLE_CUDA)
    ordinal_type _nstreams;
    cublasHandle_t _cublas_handle;
    //typedef Kokkos::View<cudaStream_t*,Kokkos::HostSpace> cuda_stream_array_host;
    typedef std::vector<cudaStream_t> cuda_stream_array_host;
    cuda_stream_array_host _cuda_streams;
    int _status;
#endif

    ///
    /// statistics
    ///
    struct {
      double t_init, t_prepare, t_solve, t_extra;
      double m_used, m_peak;
    } stat;

  public:

#if defined(KOKKOS_ENABLE_CUDA)
    inline 
    void 
    checkStatus(const char *func, const char *lib) {
      if (_status != 0) {
        printf("Error: %s, %s returns non-zero status %d\n", 
               lib, func, _status);
        std::runtime_error("checkStatus failed");
      }
    }
    inline void checkCuBlasStatus(const char *func) { checkStatus(func, "CuBlas");  }
    inline void checkCudaStatus(const char *func) { checkStatus(func, "Cuda");  }
#endif

    inline
    void
    track_alloc(const double in) {
      stat.m_used += in;
      stat.m_peak  = max(stat.m_peak, stat.m_used);
    }

    inline
    void
    track_free(const double out) {
      stat.m_used -= out;
    }

    inline
    void
    reset_stat() {
      stat.t_solve = 0;
      stat.t_extra = 0;
      stat.m_used = 0;
      stat.m_peak = 0;
    }

    inline
    void
    print_stat_init() {
      printf("  Time\n");
      printf("             time for initialization:                         %10.6f s\n", stat.t_init);
      printf("             time for compute mode classification:            %10.6f s\n", stat.t_prepare);
      printf("             total time spent:                                %10.6f s\n", (stat.t_init+stat.t_prepare));
      printf("\n");
      printf("  Memory\n");
      printf("             memory used in solve:                            %10.2f MB\n", stat.m_used/1024/1024);
      printf("             peak memory used in solve:                       %10.2f MB\n", stat.m_peak/1024/1024);
      printf("\n");
    }
      
    inline
    void
    print_stat_solve() {
      printf("  Time\n");
      printf("             total time spent:                                %10.6f s\n", (stat.t_solve+stat.t_extra));
      printf("\n");
      printf("  Memory\n");
      printf("             memory used in solve:                            %10.2f MB\n", stat.m_used/1024/1024);
      printf("\n");
    }

    inline
    void
    print_stat_memory() {
      printf("  Memory\n"); // better get zero leak
      printf("             leak (or not tracked) memory:                    %10.2f MB\n", stat.m_used/1024/1024);
    }

    inline
    void
    initialize(const ordinal_type device_function_thres,
               // const ordinal_type team_recursive_thres,
               const ordinal_type verbose = 0) {
      Kokkos::Impl::Timer timer;

      timer.reset();

      // # of supernodes
      _nsupernodes = _info.supernodes.extent(0);

      // local host supernodes info 
      _h_supernodes = Kokkos::create_mirror_view(host_memory_space(), _info.supernodes);
      Kokkos::deep_copy(_h_supernodes, _info.supernodes);

      // # of levels
      _nlevel = 0;
      {
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          _nlevel = max(_h_stree_level(sid), _nlevel);
        ++_nlevel;
      }

      // create level ptr
      _h_level_ptr = size_type_array_host("h_level_ptr", _nlevel+1);
      {
        // first count # of supernodes in each level
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          ++_h_level_ptr(_h_stree_level(sid)+1);
        
        // scan 
        for (ordinal_type i=0;i<_nlevel;++i) 
          _h_level_ptr(i+1) += _h_level_ptr(i);
      }

      // fill sids
      _h_level_sids = ordinal_type_array_host(do_not_initialize_tag("h_level_sids"), _nsupernodes);
      {
        size_type_array_host tmp_level_ptr(do_not_initialize_tag("tmp_level_ptr"), _h_level_ptr.extent(0));
        Kokkos::deep_copy(tmp_level_ptr, _h_level_ptr);
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          const ordinal_type lvl = _h_stree_level(sid);
          _h_level_sids(tmp_level_ptr(lvl)++) = sid;
        }
      }
      _level_sids = Kokkos::create_mirror_view(exec_memory_space(), _h_level_sids); 
      Kokkos::deep_copy(_level_sids, _h_level_sids);
      track_alloc(_level_sids.span()*sizeof(ordinal_type));

      // create workspace; this might needs to be created in prepare
      _h_buf_ptr = size_type_array_host(do_not_initialize_tag("h_buf_ptr"), _nsupernodes+1);
      {
        _h_buf_ptr(0) = 0;
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          const auto s = _h_supernodes(sid); 
          // anticipating other algorithms, consider tT and tB space
          _h_buf_ptr(sid+1) = s.n*_max_nrhs; //(s.n-s.m)*_max_nrhs;
        }
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          _h_buf_ptr(sid+1) += _h_buf_ptr(sid);
      }
      _buf_ptr = Kokkos::create_mirror_view(exec_memory_space(), _h_buf_ptr);
      Kokkos::deep_copy(_buf_ptr, _h_buf_ptr);
      track_alloc(_buf_ptr.span()*sizeof(size_type));

      _buf = value_type_array(do_not_initialize_tag("buf"), _h_buf_ptr(_nsupernodes));
      track_alloc(_buf.span()*sizeof(value_type));

      // cuda stream setup
#if defined(KOKKOS_ENABLE_CUDA)
      _status = cublasCreate(&_cublas_handle); checkCuBlasStatus("cublasCreate");
      _nstreams = 0;
#endif
      stat.t_init = timer.seconds();

      timer.reset();

      _h_compute_mode = ordinal_type_array_host(do_not_initialize_tag("h_compute_mode"), _nsupernodes);
      Kokkos::deep_copy(_h_compute_mode, -1);

      _device_function_thres = device_function_thres;
      _device_level_cut = min(4, _nlevel);
      {
        for (ordinal_type lvl=0;lvl<_device_level_cut;++lvl) {
          const ordinal_type 
            pbeg = _h_level_ptr(lvl), 
            pend = _h_level_ptr(lvl+1);
          for (ordinal_type p=pbeg;p<pend;++p) {
            const ordinal_type sid = _h_level_sids(p);
            _h_compute_mode(sid) = 0;
          }
        }
      }
      
      _team_recursive_level_cut = _nlevel;
      //_team_recursive_thres = team_recursive_thres;
      // {
      //   //recursive function on cuda does not work well; need workspace for postorder of children
      //   for (ordinal_type lvl=_device_level_cut;lvl<_nlevel;++lvl) {
      //     const ordinal_type 
      //       pbeg = _h_level_ptr(lvl), 
      //       pend = _h_level_ptr(lvl+1),
      //       pcnt = pend - pbeg;
      //     if (pcnt > _team_recursive_thres) { 
      //       _team_recursive_level_cut = lvl;
      //       const ordinal_type 
      //         rbeg = _h_level_ptr(_team_recursive_level_cut),
      //         rend = _h_level_ptr(_team_recursive_level_cut+1);
      //       for (ordinal_type r=rbeg;r<rend;++r) {
      //         const ordinal_type sid = _h_level_sids(r);
      //         _h_compute_mode(sid) = 2;
      //       }
      //       break;
      //     }
      //   }
      // }  

      {        
        for (ordinal_type lvl=_device_level_cut;lvl<_team_recursive_level_cut;++lvl) {          
          const ordinal_type 
            pbeg = _h_level_ptr(lvl), 
            pend = _h_level_ptr(lvl+1);
          for (ordinal_type p=pbeg;p<pend;++p) {
            const ordinal_type sid = _h_level_sids(p);
            const auto s = _h_supernodes(sid); 
            const ordinal_type m = s.m, n = s.n-s.m;
            if (m > _device_function_thres || 
                n > _device_function_thres) {
              _h_compute_mode(sid) = 0;
            } else {
              _h_compute_mode(sid) = 1;
            }
          }
        }
      }

      _compute_mode = Kokkos::create_mirror_view(exec_memory_space(), _h_compute_mode);
      Kokkos::deep_copy(_compute_mode, _h_compute_mode);
      track_alloc(_compute_mode.span()*sizeof(ordinal_type));

      stat.t_prepare = timer.seconds();

      if (verbose) {
        printf("Summary: TriSolveTools (Initialize)\n");
        printf("===================================\n");
        print_stat_init();
      }
    }

    inline
    void
    release(const ordinal_type verbose = 0) {
      track_free(_buf_ptr.span()*sizeof(size_type));
      track_free(_buf.span()*sizeof(value_type));
      track_free(_compute_mode.span()*sizeof(ordinal_type));
      track_free(_level_sids.span()*sizeof(ordinal_type));
      if (verbose) {
        printf("Summary: TriSolveTools (Release)\n");
        printf("================================\n");
        print_stat_memory();
      }
    }

    TriSolveTools()  {}
    
    TriSolveTools(const TriSolveTools &b) = default;

    TriSolveTools(// input permutation
                  const ordinal_type_array &perm,
                  const ordinal_type_array &peri,
                  // supernodes
                  const supernode_info_type &info,
                  const ordinal_type_array_host &h_stree_level,
                  const ordinal_type max_nrhs)
      : _perm(perm), _peri(peri),
        _info(info),
        _h_stree_level(h_stree_level),
        _max_nrhs(max_nrhs) 
    {}

    virtual~TriSolveTools() {
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previously created streams
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamDestroy(_cuda_streams[i]); checkCudaStatus("cudaStreamDestroy");
      }
      _cuda_streams.clear();
      _status = cublasDestroy(_cublas_handle); checkCuBlasStatus("cublasDestroy");
#endif
    }

    inline
    void
    createCudaStream(const ordinal_type nstreams) {
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previously created streams
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamDestroy(_cuda_streams[i]); checkCudaStatus("cudaStreamDestroy");
      }
      // new streams
      _nstreams = nstreams;
      //_cuda_streams = cuda_stream_array_host(do_not_initialize_tag("cuda streams"), _nstreams);
      _cuda_streams.clear();
      _cuda_streams.resize(_nstreams);
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamCreate(&_cuda_streams[i]); checkCudaStatus("cudaStreamCreate");
      }
#else
      // no-op
#endif
    }
    
    inline
    void
    solveLowerOnDevice(const ordinal_type pbeg, 
                       const ordinal_type pend,
                       const value_type_matrix &t) {
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), zero(0);
      for (ordinal_type p=pbeg,q=0;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_compute_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          _status = cublasSetStream(_cublas_handle, _cuda_streams[q%_nstreams]); checkCuBlasStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          value_type *ptr = s.buf; 
          {
            const ordinal_type m = s.m, n = s.n - s.m;
            if (m > 0) {
              // solve diag
              const ordinal_type offm = s.row_begin;
              UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
              auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
              if (is_host) {
                const ordinal_type member(0); // dummy member for testing
                Trsv<Uplo::Upper,Trans::ConjTranspose,Algo::External>
                  ::invoke(member, Diag::NonUnit(), AL, tT);
              } else {
                /// cublas
#if defined(KOKKOS_ENABLE_CUDA)
                _status = cublasDtrsv(_cublas_handle, CUBLAS_FILL_MODE_UPPER,
                                      CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT,
                                      m, 
                                      AL.data(), AL.stride_1(),
                                      tT.data(), tT.stride_0()); checkCuBlasStatus("cublasDtrsv");
#endif
              }

              if (n > 0) {
                // solve offdiag
                UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
                UnmanagedViewType<value_type_matrix> tB(_buf.data()+_h_buf_ptr(sid), n, nrhs);
                if (is_host) {
                  const ordinal_type member(0); // dummy member for testing
                  Gemv<Trans::ConjTranspose,Algo::External>
                    ::invoke(member, minus_one, AR, tT, zero, tB);
                } else {
                  // cublas
#if defined(KOKKOS_ENABLE_CUDA)
                  _status = cublasDgemv(_cublas_handle, CUBLAS_OP_C,
                                        m, n, 
                                        &minus_one,
                                        AR.data(), AR.stride_1(),
                                        tT.data(), tT.stride_0(),
                                        &zero,
                                        tB.data(), tB.stride_0()); checkCuBlasStatus("cublasDgemv");
#endif

                }
              }
            }
          }
        }
      }
    }
    
    inline
    void
    solveUpperOnDevice(const ordinal_type pbeg,
                       const ordinal_type pend,
                       const value_type_matrix &t) {
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), one(1);
      for (ordinal_type p=pbeg,q=0;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_compute_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          _status = cublasSetStream(_cublas_handle, _cuda_streams[q%_nstreams]); checkCuBlasStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          value_type *ptr = s.buf; 
          {
            const ordinal_type m = s.m, n = s.n - s.m;
            {
              const UnmanagedViewType<value_type_matrix> tB(_buf.data()+_h_buf_ptr(sid), n, nrhs); 
              if (m > 0) {
                // solve
                const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
                const ordinal_type offm = s.row_begin;
                const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
                
                if (n > 0) {
                  const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
                  if (is_host) {
                    const ordinal_type member(0); // dummy member for testing
                    Gemv<Trans::NoTranspose,Algo::External>
                      ::invoke(member, minus_one, AR, tB, one, tT);
                  } else {
                    // cublas
#if defined(KOKKOS_ENABLE_CUDA)
                    _status = cublasDgemv(_cublas_handle, CUBLAS_OP_N,
                                          m, n, 
                                          &minus_one,
                                          AR.data(), AR.stride_1(),
                                          tB.data(), tB.stride_0(),
                                          &one,
                                          tT.data(), tT.stride_0()); checkCuBlasStatus("cublasDgemv");
#endif
                  }
                }
                if (is_host) {
                  const ordinal_type member(0); // dummy member for testing
                  Trsv<Uplo::Upper,Trans::NoTranspose,Algo::External>
                    ::invoke(member, Diag::NonUnit(), AL, tT);
                } else {
                  // cublas
#if defined(KOKKOS_ENABLE_CUDA)
                  _status = cublasDtrsv(_cublas_handle, CUBLAS_FILL_MODE_UPPER,
                                        CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                        m, 
                                        AL.data(), AL.stride_1(),
                                        tT.data(), tT.stride_0()); checkCuBlasStatus("cublasDtrsv");
#endif
                }
              }
            }
          }
        }
      }
    }

    /// 
    /// Level set solve
    /// ---------------

    inline
    void
    solveCholesky(const value_type_matrix &x,   // solution
                  const value_type_matrix &b,   // right hand side
                  const value_type_matrix &t,
                  const ordinal_type verbose = 0) { // temporary workspace (store permuted vectors)
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) ||
                               x.extent(1) != b.extent(1) ||
                               x.extent(0) != t.extent(0) ||
                               x.extent(1) != t.extent(1), std::logic_error,
                               "x, b, t, and w dimensions do not match");

      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() ||
                               x.data() == t.data(), std::logic_error,
                               "x, b, t, and w have the same data pointer");
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;

      // solve U^{H} (U x) = b 
      const ordinal_type nrhs = x.extent(1);
      TACHO_TEST_FOR_EXCEPTION(nrhs > _max_nrhs, std::logic_error,
                               "nrhs is bigger than max nrhs");
      
      Kokkos::Impl::Timer timer;



      // 0. permute and copy b -> t
      timer.reset();
      applyRowPermutationToDenseMatrix(t, b, _perm);      
      stat.t_extra = timer.seconds();

      timer.reset();
      { 
        {
          typedef TeamFunctor_SolveLowerChol<supernode_info_type> functor_type;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveUpdateTag> team_policy_solve_update;

          functor_type functor(_info, 
                               _compute_mode,
                               _level_sids,
                               t,
                               _buf_ptr,
                               _buf);

          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);
          
          //  1. U^{H} w = t
#if 0
          if (_team_recursive_level_cut < _nlevel) {
            const ordinal_type lvl = _team_recursive_level_cut;
            const ordinal_type 
              pbeg = _h_level_ptr(lvl), 
              pend = _h_level_ptr(lvl+1),
              pcnt = pend - pbeg;

            functor.setRange(pbeg, pend);
            team_policy_solve_update policy_solve_update(1,1,1);
            if (is_host) {
              policy_solve_update = team_policy_solve_update(pcnt, 1, 1);
            } else {
              policy_solve_update = team_policy_solve_update(pcnt, 8, 8);
            }
            Kokkos::parallel_for("solve lower and update", policy_solve_update, functor);
          }
#endif          
          {
            for (ordinal_type lvl=(_team_recursive_level_cut-1);lvl>=0;--lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;
              
              solveLowerOnDevice(pbeg, pend, t);
              
              functor.setRange(pbeg, pend);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                policy_solve  = team_policy_solve(pcnt, 32,  8);
                policy_update = team_policy_update(pcnt, 1, 32);
              }
              if (lvl < _device_level_cut) {
                // do nothing
              } else {
                Kokkos::parallel_for("solve lower", policy_solve, functor);
              }
              Kokkos::fence();

              Kokkos::parallel_for("update lower", policy_update, functor); 
              Kokkos::fence();
            }
          }
        } // end of lower tri solve

        {
          typedef TeamFunctor_SolveUpperChol<supernode_info_type> functor_type;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveUpdateTag> team_policy_solve_update;

          functor_type functor(_info, 
                               _compute_mode,
                               _level_sids,
                               t,
                               _buf_ptr,
                               _buf);

          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);
          
          //  2. U t = w;
          {
            for (ordinal_type lvl=0;lvl<_team_recursive_level_cut;++lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;
              
              functor.setRange(pbeg, pend);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                policy_solve  = team_policy_solve(pcnt, 32,  8);
                policy_update = team_policy_update(pcnt, 1, 32);
              }
              Kokkos::parallel_for("update upper", policy_update, functor);              
              Kokkos::fence();

              if (lvl < _device_level_cut) {
                // do nothing
              } else {
                Kokkos::parallel_for("solve lower", policy_solve, functor); 
              }
              solveUpperOnDevice(pbeg, pend, t);
              Kokkos::fence();
            }
          }
#if 0
          if (_team_recursive_level_cut < _nlevel) {
            const ordinal_type lvl = _team_recursive_level_cut;
            const ordinal_type 
              pbeg = _h_level_ptr(lvl), 
              pend = _h_level_ptr(lvl+1),
              pcnt = pend - pbeg;

            functor.setRange(pbeg, pend);
            team_policy_solve_update policy_solve_update(1,1,1);
            if (is_host) {
              policy_solve_update = team_policy_solve_update(pcnt, 1, 1);
            } else {
              policy_solve_update = team_policy_solve_update(pcnt, 8, 8);
            }
            Kokkos::parallel_for("solve lower and update", policy_solve_update, functor);
          }
#endif
        }/// end of upper tri solve

      } // end of solve
      stat.t_solve = timer.seconds();

      // permute and copy t -> x
      timer.reset();
      applyRowPermutationToDenseMatrix(x, t, _peri);
      stat.t_extra += timer.seconds();

      if (verbose) {
        printf("Summary: TriSolveTools (ParallelSolve: %3d)\n", nrhs);
        printf("===========================================\n");
        print_stat_solve();
      }
    }

  };

}
#endif




    // inline
    // void
    // solveCholesky_Serial(const value_type_matrix &x,   // solution
    //                      const value_type_matrix &b,   // right hand side
    //                      const value_type_matrix &t,
    //                      const ordinal_type verbose = 0) { // temporary workspace (store permuted vectors)
    //   TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) ||
    //                            x.extent(1) != b.extent(1) ||
    //                            x.extent(0) != t.extent(0) ||
    //                            x.extent(1) != t.extent(1), std::logic_error,
    //                            "x, b, t, and w dimensions do not match");

    //   TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() ||
    //                            x.data() == t.data(), std::logic_error,
    //                            "x, b, t, and w have the same data pointer");
    //   // solve U^{H} (U x) = b 
    //   const ordinal_type nrhs = 1; //x.extent(1);
    //   const ordinal_type member = 0; // dummy member for testing
    //   const value_type one(1), zero(0);

    //   // 0. permute and copy b -> t
    //   applyRowPermutationToDenseMatrix(t, b, _perm);      

    //   { 
    //     //  1. U^{H} w = t
    //     for (ordinal_type lvl=(_nlevel-1);lvl>=0;--lvl) {
    //       const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl+1);
    //       for (ordinal_type p=pbeg;p<pend;++p) {
    //         const ordinal_type sid = _h_level_sids(p);
    //         const auto &s = _info.supernodes(sid);
    //         value_type *ptr = s.buf; 
    //         {
    //           const ordinal_type m = s.m, n = s.n - s.m;
    //           if (m > 0) {
    //             // solve diag
    //             const ordinal_type offm = s.row_begin;
    //             UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
    //             auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
    //             Trsv<Uplo::Upper,Trans::ConjTranspose,Algo::External>
    //               ::invoke(member, Diag::NonUnit(), AL, tT);
                
    //             if (n > 0) {
    //               // solve offdiag
    //               UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
    //               UnmanagedViewType<value_type_matrix> tB(&_buf(_buf_ptr(sid)), n, nrhs);
    //               Gemv<Trans::ConjTranspose,Algo::External>
    //                 ::invoke(member, -one, AR, tT, zero, tB);

    //               // update
    //               const ordinal_type 
    //                 sbeg = s.sid_col_begin + 1, send = s.sid_col_end - 1;
    //               for (ordinal_type i=sbeg,is=0;i<send;++i) {
    //                 const ordinal_type 
    //                   tbeg = _info.sid_block_colidx(i).second,
    //                   tend = _info.sid_block_colidx(i+1).second;
              
    //                 for (ordinal_type it=tbeg;it<tend;++it,++is) {
    //                   const ordinal_type row = _info.gid_colidx(s.gid_col_begin + it);
    //                   for (ordinal_type j=0;j<nrhs;++j) 
    //                     Kokkos::atomic_add(&t(row,j), tB(is,j));
    //                 }
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }

    //     //  2. U t = w;
    //     for (ordinal_type lvl=0;lvl<_nlevel;++lvl) {
    //       const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl+1);
    //       for (ordinal_type p=pbeg;p<pend;++p) {
    //         const ordinal_type sid = _h_level_sids(p);
    //         const auto &s = _info.supernodes(sid);
    //         value_type *ptr = s.buf; 
    //         {
    //           const ordinal_type m = s.m, n = s.n - s.m;
    //           {
    //             const UnmanagedViewType<value_type_matrix> tB(&_buf(_buf_ptr(sid)), n, nrhs); 
    //             if (n > 0) {
    //               // update
    //               const ordinal_type goffset = s.gid_col_begin + s.m;
    //               for (ordinal_type j=0;j<nrhs;++j) {
    //                 for (ordinal_type i=0;i<n;++i) {
    //                   const ordinal_type row = _info.gid_colidx(i+goffset);
    //                   tB(i,j) = t(row,j);
    //                 }
    //               }
    //             }

    //             if (m > 0) {
    //               // solve
    //               const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
    //               const ordinal_type offm = s.row_begin;
    //               const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());

    //               if (n > 0) {
    //                 const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
    //                 Gemv<Trans::NoTranspose,Algo::External>
    //                   ::invoke(member, -one, AR, tB, one, tT);
    //               }
    //               Trsv<Uplo::Upper,Trans::NoTranspose,Algo::External>
    //                 ::invoke(member, Diag::NonUnit(), AL, tT);
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }

    //   // permute and copy t -> x
    //   applyRowPermutationToDenseMatrix(x, t, _peri);

    //   if (verbose) {
    //     printf("Summary: TriSolveTools (ParallelSolve: %3d)\n", ordinal_type(x.extent(1)));
    //     printf("===========================================\n");
    //   }
    // }
