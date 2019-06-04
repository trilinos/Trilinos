#ifndef __TACHOEXP_GRAPH_TOOLS_SCOTCH_HPP__
#define __TACHOEXP_GRAPH_TOOLS_SCOTCH_HPP__

/// \file TachoExp_GraphTools_Scotch.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(HAVE_SHYLU_NODETACHO_SCOTCH)
#include "TachoExp_Util.hpp"
#include "TachoExp_Graph.hpp"

#include "scotch.h"

namespace Tacho {

  namespace Experimental {

    class GraphTools_Scotch {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<SCOTCH_Num*,host_exec_space> ordinal_type_array;

      enum : int { DefaultRandomSeed = -1 };

    private:
        
      // scotch main data structure
      SCOTCH_Graph _graph;
      SCOTCH_Num _strat;
      int _level;
    
      // scotch output
      ordinal_type _cblk;
      ordinal_type_array _perm,_peri,_range,_tree;

      // status flag
      bool _is_ordered, _verbose;

    public:
      GraphTools_Scotch() = default;
      GraphTools_Scotch(const GraphTools_Scotch &b) = default;

      ///
      /// construction of scotch graph
      ///
      GraphTools_Scotch(const Graph &g) {
        _is_ordered = false;
        _verbose = false;
      
        // input 
        const ordinal_type base = 0;
        const ordinal_type m = g.NumRows();
        const size_type nnz = g.NumNonZeros();

        // scotch control parameter
        _strat = 0;
        _level = 0;
      
        // output
        _cblk  = 0;
        _perm  = ordinal_type_array("Scotch::PermutationArray", m);
        _peri  = ordinal_type_array("Scotch::InvPermutationArray", m);
        _range = ordinal_type_array("Scotch::RangeArray", m);
        _tree  = ordinal_type_array("Scotch::TreeArray", m);
      
        // construct scotch graph
        int ierr = 0;
        const SCOTCH_Num *rptr_ptr = reinterpret_cast<const SCOTCH_Num*>(g.RowPtr().data());
        const SCOTCH_Num *cidx_ptr = reinterpret_cast<const SCOTCH_Num*>(g.ColIdx().data());

        ierr = SCOTCH_graphInit(&_graph);
        TACHO_TEST_FOR_EXCEPTION(ierr, std::runtime_error, "Failed in SCOTCH_graphInit");

        ierr = SCOTCH_graphBuild(&_graph,             // scotch graph
                                 base,                // base value
                                 m,                   // # of vertices
                                 rptr_ptr,            // column index array pointer begin
                                 rptr_ptr+1,          // column index array pointer end
                                 NULL,                // weights on vertices (optional)
                                 NULL,                // label array on vertices (optional)
                                 nnz,                 // # of nonzeros
                                 cidx_ptr,            // column index array
                                 NULL);               // edge load array (optional)
        TACHO_TEST_FOR_EXCEPTION(ierr, std::runtime_error, "Failed in SCOTCH_graphBuild");

        ierr = SCOTCH_graphCheck(&_graph);
        TACHO_TEST_FOR_EXCEPTION(ierr, std::runtime_error, "Failed in SCOTCH_graphCheck");
      }
      virtual~GraphTools_Scotch() {
        SCOTCH_graphFree(&_graph);
      }

      ///
      /// setup scotch parameters
      ///

      void setVerbose(const bool verbose) { _verbose = verbose; }
      void setSeed(const int seed = DefaultRandomSeed) {
        if (seed != DefaultRandomSeed) {
          SCOTCH_randomSeed(seed);
          SCOTCH_randomReset();
        }
      }

      void setStrategy(const SCOTCH_Num strat = 0) {
        // a typical choice
        //(SCOTCH_STRATLEVELMAX));//   |
        //SCOTCH_STRATLEVELMIN   |
        //SCOTCH_STRATLEAFSIMPLE |
        //SCOTCH_STRATSEPASIMPLE);
        _strat = strat;
      }

      void setTreeLevel(const unsigned int level = 0) {
        _level = level;
      }

      ///
      /// setup scotch parameters
      ///

      void reorder(const ordinal_type treecut = 0) {
        int ierr = 0;

        // pointers for global graph ordering
        ordinal_type *perm  = _perm.data();
        ordinal_type *peri  = _peri.data();
        ordinal_type *range = _range.data();
        ordinal_type *tree  = _tree.data();

        {
          // set desired tree level
          if (_strat & SCOTCH_STRATLEVELMAX ||
              _strat & SCOTCH_STRATLEVELMIN) {
            TACHO_TEST_FOR_EXCEPTION(_level == 0, 
                                     std::logic_error,
                                     "SCOTCH_STRATLEVEL(MIN/MAX) is used but level is not specified");
          }
          const int level = max(1, _level-treecut);

          SCOTCH_Strat stradat;
          SCOTCH_Num straval = _strat;
        
          ierr = SCOTCH_stratInit(&stradat);
          TACHO_TEST_FOR_EXCEPTION(ierr, 
                                   std::runtime_error,
                                   "Failed in SCOTCH_stratInit");


          // if both are zero, do not build strategy
          if (_strat || _level) {
            if (_verbose)
              std::cout << "GraphTools_Scotch:: User provide a strategy and/or level" << std::endl
                        << "                    strategy = " << _strat << ", level =  " << _level << ", treecut = " << treecut << std::endl
                        << "                    strategy & SCOTCH_STRATLEVELMAX   = " << (_strat & SCOTCH_STRATLEVELMAX) << std::endl
                        << "                    strategy & SCOTCH_STRATLEVELMIN   = " << (_strat & SCOTCH_STRATLEVELMIN) << std::endl
                        << "                    strategy & SCOTCH_STRATLEAFSIMPLE = " << (_strat & SCOTCH_STRATLEAFSIMPLE) << std::endl
                        << "                    strategy & SCOTCH_STRATSEPASIMPLE = " << (_strat & SCOTCH_STRATSEPASIMPLE) << std::endl
                        << std::endl;
            ierr = SCOTCH_stratGraphOrderBuild(&stradat, straval, level, 0.2);
            TACHO_TEST_FOR_EXCEPTION(ierr, 
                                     std::runtime_error,
                                     "Failed in SCOTCH_stratGraphOrderBuild");
          }
          ierr = SCOTCH_graphOrder(&_graph,
                                   &stradat,
                                   perm,
                                   peri,
                                   &_cblk,
                                   range,
                                   tree);
          TACHO_TEST_FOR_EXCEPTION(ierr, 
                                   std::runtime_error,
                                   "Failed in SCOTCH_graphOrder");
          SCOTCH_stratExit(&stradat);
        }
        _is_ordered = true;
      
        if (_verbose)
          std::cout << "GraphTools_Scotch:: # of block partitions of rows/columns = " << _cblk << std::endl;
      }

      ordinal_type_array PermVector()    const { return _perm; }
      ordinal_type_array InvPermVector() const { return _peri; }
    
      ordinal_type_array RangeVector()   const { return _range; }
      ordinal_type_array TreeVector()    const { return _tree; }
    
      ordinal_type NumBlocks() const { return _cblk; }
      ordinal_type TreeLevel() const {
        ordinal_type r_val;
        if (_strat & SCOTCH_STRATLEVELMAX ||
            _strat & SCOTCH_STRATLEVELMIN)
          r_val = _level;
        else
          r_val = 0;
        return r_val;
      }
    
      std::ostream& showMe(std::ostream &os, const bool detail = false) const {
        std::streamsize prec = os.precision();
        os.precision(4);
        os << std::scientific;

        if (_is_ordered)
          os << " -- Scotch Ordering -- " << std::endl
             << "    CBLK   = " << _cblk << std::endl
             << "  PERM     PERI     RANG     TREE" << std::endl;
        else 
          os << " -- Not Ordered -- " << std::endl;

        if (detail) {
          const ordinal_type w = 6, m = _perm.extent(0);
          for (ordinal_type i=0;i<m;++i)
            os << std::setw(w) << _perm[i] << "   "
               << std::setw(w) << _peri[i] << "   "
               << std::setw(w) << _range[i] << "   "
               << std::setw(w) << _tree[i] << std::endl;
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
