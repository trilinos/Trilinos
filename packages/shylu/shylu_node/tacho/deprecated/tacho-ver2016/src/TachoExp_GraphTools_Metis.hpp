#ifndef __TACHOEXP_GRAPH_TOOLS_METIS_HPP__
#define __TACHOEXP_GRAPH_TOOLS_METIS_HPP__

/// \file TachoExp_GraphTools_Scotch.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(HAVE_SHYLU_NODETACHO_METIS)
#include "TachoExp_Util.hpp"
#include "TachoExp_Graph.hpp"

#include "metis.h"

namespace Tacho {

  namespace Experimental {

    class GraphTools_Metis {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<idx_t*,host_exec_space> ordinal_type_array;

    private:
        
      // metis main data structure
      idx_t _nvts;
      ordinal_type_array _xadj, _adjncy, _vwgt;
      
      idx_t _options[METIS_NOPTIONS];

      // metis output
      ordinal_type_array _perm, _peri;

      // status flag
      bool _is_ordered, _verbose;

    public:
      GraphTools_Metis() = default;
      GraphTools_Metis(const GraphTools_Metis &b) = default;

      ///
      /// construction of scotch graph
      ///
      GraphTools_Metis(const Graph &g) {
        _is_ordered = false;
        _verbose = false;
        
        // input 
        _nvts = g.NumRows();
        _xadj = g.RowPtr();
        _adjncy = g.ColIdx();
        _vwgt = ordinal_type_array();
        
        METIS_SetDefaultOptions(_options);
        _options[METIS_OPTION_NUMBERING] = 0;

        // output
        _perm  = ordinal_type_array("Metis::PermutationArray", _nvts);
        _peri  = ordinal_type_array("Metis::InvPermutationArray", _nvts);    
      }
      virtual~GraphTools_Metis() {}

      ///
      /// setup metis parameters
      ///

      void setVerbose(const bool verbose) { _verbose = verbose; }
      void setOption(const int id, const idx_t value) {
        _options[id] = value;
      }

      ///
      /// reorder by metis
      ///

      void reorder() {
        int ierr = 0;

        idx_t *xadj   = (idx_t*)_xadj.data();
        idx_t *adjncy = (idx_t*)_adjncy.data();
        idx_t *vwgt   = (idx_t*)_vwgt.data();

        idx_t *perm   = (idx_t*)_perm.data();
        idx_t *peri   = (idx_t*)_peri.data();

        ierr = METIS_NodeND(&_nvts, xadj, adjncy, vwgt, _options, 
                            perm, peri);
        TACHO_TEST_FOR_EXCEPTION(ierr != METIS_OK, 
                                 std::runtime_error,
                                 "Failed in METIS_NodeND");
        _is_ordered = true;
      }

      ordinal_type_array PermVector()    const { return _perm; }
      ordinal_type_array InvPermVector() const { return _peri; }
        
      std::ostream& showMe(std::ostream &os, const bool detail = false) const {
        std::streamsize prec = os.precision();
        os.precision(4);
        os << std::scientific;

        if (_is_ordered)
          os << " -- Metis Ordering -- " << std::endl
             << "  PERM     PERI     " << std::endl;
        else 
          os << " -- Not Ordered -- " << std::endl;

        if (detail) {
          const ordinal_type w = 6, m = _perm.extent(0);
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
}
#endif
#endif
