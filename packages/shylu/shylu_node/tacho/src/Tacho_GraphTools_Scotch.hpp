#ifndef __TACHO_GRAPH_TOOLS_SCOTCH_HPP__
#define __TACHO_GRAPH_TOOLS_SCOTCH_HPP__

/// \file Tacho_GraphTools_Scotch.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(TACHO_HAVE_SCOTCH)
#include "Tacho_Util.hpp"
#include "Tacho_Graph.hpp"

#include "scotch.h"

namespace Tacho {

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

      void reorder(const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;
        double t_scotch = 0;

        _verbose = verbose;

        const int treecut = 0;
        int ierr = 0;

        // pointers for global graph ordering
        ordinal_type *perm  = _perm.data();
        ordinal_type *peri  = _peri.data();
        ordinal_type *range = _range.data();
        ordinal_type *tree  = _tree.data();

        timer.reset();
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
        t_scotch = timer.seconds();
        _is_ordered = true;
      
        if (_verbose) {
          printf("Summary: GraphTools (Scotch)\n");
          printf("===========================\n");          
          printf("  Time\n");
          printf("             time for reordering: %10.6f s\n", t_scotch);
          printf("\n");
          if (_strat || _level) {
            printf("  User provided strategy ( %d ) and/or level ( %d )\n", _strat, _level);
            printf("             strategy & SCOTCH_STRATLEVELMAX:   %3d\n", (_strat & SCOTCH_STRATLEVELMAX));
            printf("             strategy & SCOTCH_STRATLEVELMIN:   %3d\n", (_strat & SCOTCH_STRATLEVELMIN));
            printf("             strategy & SCOTCH_STRATLEAFSIMPLE: %3d\n", (_strat & SCOTCH_STRATLEAFSIMPLE));
            printf("             strategy & SCOTCH_STRATSEPASIMPLE: %3d\n", (_strat & SCOTCH_STRATSEPASIMPLE));
            printf("\n");
          }
          printf("  Partitions\n");
          printf("             number of block partitions: %3d\n", _cblk);
          printf("\n");
        }
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
#endif
#endif
