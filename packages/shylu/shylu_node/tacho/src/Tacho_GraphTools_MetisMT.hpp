#ifndef __TACHO_GRAPH_TOOLS_METIS_MT_HPP__
#define __TACHO_GRAPH_TOOLS_METIS_MT_HPP__

/// \file Tacho_GraphTools_Metis_MT.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(TACHO_HAVE_METIS_MT)
#include "Tacho_Util.hpp"
#include "Tacho_Graph.hpp"

#include "mtmetis.h"

namespace Tacho {

    class GraphTools_MetisMT {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;

      typedef Kokkos::View<mtmetis_vtx_type*,host_exec_space> mtmetis_vtx_type_array;
      typedef Kokkos::View<mtmetis_adj_type*,host_exec_space> mtmetis_adj_type_array;
      typedef Kokkos::View<mtmetis_wgt_type*,host_exec_space> mtmetis_wgt_type_array;
      typedef Kokkos::View<mtmetis_pid_type*,host_exec_space> mtmetis_pid_type_array;

      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array;

    private:
        
      // metis main data structure
      mtmetis_vtx_type _nvts;
      mtmetis_vtx_type_array _xadj;
      mtmetis_adj_type_array _adjncy;
      mtmetis_wgt_type_array _vwgt;
      
      double _options[MTMETIS_NOPTIONS];

      // metis output
      mtmetis_pid_type_array _perm_t, _peri_t;
      ordinal_type_array _perm, _peri;

      // status flag
      bool _is_ordered, _verbose;

    public:
      GraphTools_MetisMT() = default;
      GraphTools_MetisMT(const GraphTools_MetisMT &b) = default;

      ///
      /// construction of scotch graph
      ///
      GraphTools_MetisMT(const Graph &g) {
        _is_ordered = false;
        _verbose = false;
        
        // input 
        _nvts = g.NumRows();

        _xadj   = mtmetis_vtx_type_array("vtx_type_xadj",   g.RowPtr().dimension_0());
        _adjncy = mtmetis_adj_type_array("adj_type_adjncy", g.ColIdx().dimension_0());
        _vwgt   = mtmetis_wgt_type_array();

        const auto &g_row_ptr = g.RowPtr();
        const auto &g_col_idx = g.ColIdx();

        for (ordinal_type i=0;i<static_cast<ordinal_type>(_xadj.dimension_0());++i)
          _xadj(i) = g_row_ptr(i);
        for (ordinal_type i=0;i<static_cast<ordinal_type>(_adjncy.dimension_0());++i)
          _adjncy(i) = g_col_idx(i);

        // default
        for (ordinal_type i=0;i<static_cast<ordinal_type>(MTMETIS_NOPTIONS);++i)
          _options[i] = MTMETIS_VAL_OFF;

        // by default, metis use
        //   # of threads : omp_get_max_threads
        //   seed : (unsigned int)time(NULL)
        // internal verbose options are :
        //   MTMETIS_VERBOSITY_NONE,
        //   MTMETIS_VERBOSITY_LOW,
        //   MTMETIS_VERBOSITY_MEDIUM,
        //   MTMETIS_VERBOSITY_HIGH,
        //   MTMETIS_VERBOSITY_MAXIMUM

        _options[MTMETIS_OPTION_NTHREADS]  = host_exec_space::thread_pool_size(0); // from kokkos
        //_options[MTMETIS_OPTION_SEED]      = 0; // for testing, use the same seed now
        //_options[MTMETIS_OPTION_PTYPE]     = MTMETIS_PTYPE_ND; // when explicit interface is used
        //_options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
        //_options[MTMETIS_OPTION_METIS]     = 1; // flag to use serial metis

        _perm_t = mtmetis_pid_type_array("pid_type_perm", _nvts);
        _peri_t = mtmetis_pid_type_array("pid_type_peri", _nvts);    

        // output
        _perm  = ordinal_type_array("MetisMT::PermutationArray", _nvts);
        _peri  = ordinal_type_array("MetisMT::InvPermutationArray", _nvts);    
      }
      virtual~GraphTools_MetisMT() {}

      ///
      /// setup metis parameters
      ///

      void setVerbose(const bool verbose) { _verbose = verbose; }
      void setOption(const int id, const double value) {
        _options[id] = value;
      }

      ///
      /// reorder by metis
      ///

      void reorder(const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;
        double t_metis = 0; 

        int ierr = 0;

        mtmetis_vtx_type *xadj   = (mtmetis_vtx_type*)_xadj.data();
        mtmetis_adj_type *adjncy = (mtmetis_adj_type*)_adjncy.data();
        mtmetis_wgt_type *vwgt   = (mtmetis_wgt_type*)_vwgt.data();

        mtmetis_pid_type *perm   = (mtmetis_pid_type*)_perm_t.data();
        mtmetis_pid_type *peri   = (mtmetis_pid_type*)_peri_t.data();

        timer.reset();
        ierr = MTMETIS_NodeND(&_nvts, xadj, adjncy, vwgt, _options, 
                              perm, peri);
        t_metis = timer.seconds();

        for (mtmetis_vtx_type i=0;i<_nvts;++i) {
          _perm(i) = _perm_t(i);
          _peri(i) = _peri_t(i);
        }

        TACHO_TEST_FOR_EXCEPTION(ierr != MTMETIS_SUCCESS, 
                                 std::runtime_error,
                                 "Failed in METIS_NodeND");
        _is_ordered = true;

        if (verbose) {
          printf("Summary: GraphTools (MetisMT)\n");
          printf("=============================\n");

          switch (verbose) {
          case 1: {
            printf("  Time\n");
            printf("             time for reordering: %10.6f s\n", t_metis);
            printf("\n");
          }
          }
        }

      }

      ordinal_type_array PermVector()    const { return _perm; }
      ordinal_type_array InvPermVector() const { return _peri; }
        
      std::ostream& showMe(std::ostream &os, const bool detail = false) const {
        std::streamsize prec = os.precision();
        os.precision(4);
        os << std::scientific;

        if (_is_ordered)
          os << " -- MetisMT Ordering -- " << std::endl
             << "  PERM     PERI     " << std::endl;
        else 
          os << " -- Not Ordered -- " << std::endl;

        if (detail) {
          const ordinal_type w = 6, m = _perm.dimension_0();
          for (ordinal_type i=0;i<m;++i)
            os << std::setw(w) << _perm[i] << "   "
               << std::setw(w) << _peri[i] << "   "
               << std::endl;
        }
        os.unsetf(std::ios::scientific);
        os.precision(prec);

        return os;
      }

    };

}
#endif
#endif
