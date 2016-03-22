#ifndef __TACHO_GRAPH_TOOLS_SCOTCH_HPP__
#define __TACHO_GRAPH_TOOLS_SCOTCH_HPP__

/// \file Tacho_GraphTools_Scotch.hpp
/// \brief Interface to scotch reordering
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "scotch.h"

namespace Tacho {

  template<typename OrdinalType, typename SizeType = OrdinalType>
  class GraphTools_Scotch {
  public:
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    // force to use host space
    typedef typename
    Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space space_type;

    typedef Kokkos::View<ordinal_type*, space_type> ordinal_type_array;
    typedef Kokkos::View<size_type*,    space_type> size_type_array;

    static const int DefaultRandomSeed = -1;

  private:
    // scotch main data structure
    SCOTCH_Graph _graph;
    SCOTCH_Num _strat;
    int _level;
    
    // scotch input has no diagonal contribution
    ordinal_type _base,_m;
    ordinal_type_array _cidx;

    size_type _nnz;
    size_type_array _rptr;

    // scotch output
    ordinal_type _cblk;
    ordinal_type_array _perm,_peri,_range,_tree;

    // status flag
    bool _is_ordered;

  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (x),
    /// - Callable in KokkosFunctors (x)
    size_type NumNonZeros() const { return _nnz; }
    ordinal_type NumRows() const { return _m; }

    size_type_array RowPtrVector() const { return _rptr; }
    ordinal_type_array ColIndexVector() const { return _cidx; }
    
    ordinal_type_array PermVector()    const { return _perm; }
    ordinal_type_array InvPermVector() const { return _peri; }

    ordinal_type_array RangeVector()   const { return _range; }
    ordinal_type_array TreeVector()    const { return _tree; }

    ordinal_type NumBlocks() const { return _cblk; }
    
    GraphTools_Scotch() = default;
    GraphTools_Scotch(const GraphTools_Scotch &b) = default;
    virtual~GraphTools_Scotch() {
      SCOTCH_graphFree(&_graph);
    }
    
    void setGraph(const ordinal_type m,
                  const size_type_array rptr,
                  const ordinal_type_array cidx) {
      _is_ordered = false;
      _cblk  = 0;

      /// Scotch graph spec
      /// - no diagonals, symmetric

      _base  = 0; 
      _m     = m; 
      _nnz   = rptr[m];

      _rptr  = rptr; 
      _cidx  = cidx; 

      _perm  = ordinal_type_array("Scotch::PermutationArray", _m);
      _peri  = ordinal_type_array("Scotch::InvPermutationArray", _m);
      _range = ordinal_type_array("Scotch::RangeArray", _m);
      _tree  = ordinal_type_array("Scotch::TreeArray", _m);

      _strat = 0;
      _level = 0;

      int ierr = 0;
      ordinal_type *rptr_ptr = reinterpret_cast<ordinal_type*>(_rptr.ptr_on_device());
      ordinal_type *cidx_ptr = reinterpret_cast<ordinal_type*>(_cidx.ptr_on_device());

      ierr = SCOTCH_graphInit(&_graph);TACHO_TEST_FOR_ABORT(ierr, "Failed in SCOTCH_graphInit");
      ierr = SCOTCH_graphBuild(&_graph,             // scotch graph
                               _base,               // base value
                               _m,                  // # of vertices
                               rptr_ptr,            // column index array pointer begin
                               rptr_ptr+1,          // column index array pointer end
                               NULL,                // weights on vertices (optional)
                               NULL,                // label array on vertices (optional)
                               _nnz,                // # of nonzeros
                               cidx_ptr,            // column index array
                               NULL);               // edge load array (optional)
      TACHO_TEST_FOR_ABORT(ierr, "Failed in SCOTCH_graphBuild");  
      ierr = SCOTCH_graphCheck(&_graph);TACHO_TEST_FOR_ABORT(ierr, "Failed in SCOTCH_graphCheck");  
    }

    void setSeed(const int seed = DefaultRandomSeed) {
      if (seed != DefaultRandomSeed) {
        SCOTCH_randomSeed(seed);
        SCOTCH_randomReset();
      }
    }
    
    void setStratGraph(const SCOTCH_Num strat = 0) {
      // a typical choice
      //(SCOTCH_STRATLEVELMAX));//   |
      //SCOTCH_STRATLEVELMIN   |
      //SCOTCH_STRATLEAFSIMPLE |
      //SCOTCH_STRATSEPASIMPLE);
      _strat = strat;
    }

    void setTreeLevel(const int level = 0) {
      _level = level;
    }

    int computeOrdering(const ordinal_type treecut = 0) {
      int ierr = 0;
      
      // pointers for global graph ordering
      const ordinal_type *perm  = _perm.ptr_on_device();
      const ordinal_type *peri  = _peri.ptr_on_device();
      const ordinal_type *range = _range.ptr_on_device();
      const ordinal_type *tree  = _tree.ptr_on_device();

      {
        // set desired tree level
        const int level = (_level ? _level : max(1, int(log2(_m)-treecut))); // level = log2(_nnz)+10;

        SCOTCH_Strat stradat;
        SCOTCH_Num straval = _strat;

        ierr = SCOTCH_stratInit(&stradat);CHKERR(ierr);

        // if both are zero, do not build strategy
        if (_strat) {
          std::cout << "GraphTools_Scotch:: User provide a strategy and/or level" << std::endl
                    << "                    strategy = " << _strat << ", level =  " << _level << std::endl
                    << "                    strategy & SCOTCH_STRATLEVELMAX   = " << (_strat & SCOTCH_STRATLEVELMAX) << std::endl
                    << "                    strategy & SCOTCH_STRATLEVELMIN   = " << (_strat & SCOTCH_STRATLEVELMIN) << std::endl
                    << "                    strategy & SCOTCH_STRATLEAFSIMPLE = " << (_strat & SCOTCH_STRATLEAFSIMPLE) << std::endl
                    << "                    strategy & SCOTCH_STRATSEPASIMPLE = " << (_strat & SCOTCH_STRATSEPASIMPLE) << std::endl
                    << std::endl;
          
          if ((SCOTCH_STRATLEVELMAX & _strat) && _level <= 1) {
            std::cout << "GraphTools_Scotch:: Strategy SCOTCH_STRATLEVELMAX is used with level below 1" << std::endl
                      << "                    Provided strategy is ignored" << std::endl
                      << std::endl;
          } else {
            ierr = SCOTCH_stratGraphOrderBuild(&stradat, straval, level, 0.2);
            TACHO_TEST_FOR_ABORT(ierr, "Failed in SCOTCH_stratGraphOrderBuild");  
          }
        }
        ierr = SCOTCH_graphOrder(&_graph,
                                 &stradat,
                                 perm,
                                 peri,
                                 &_cblk,
                                 range,
                                 tree);
        TACHO_TEST_FOR_ABORT(ierr, "Failed in SCOTCH_graphOrder");  
        SCOTCH_stratExit(&stradat);
      }


      {
        ordinal_type nroot = 0;
        for (ordinal_type i=0;i<_cblk;++i)
          nroot += (_tree[i] == -1);

        if (nroot) {
          _tree [_cblk]   = -1;          // dummy root
          _range[_cblk+1] = _range[_cblk]; // zero range for the dummy root
          
          for (ordinal_type i=0;i<_cblk;++i)
            if (_tree[i] == -1)           // multiple roots becomes children of the dummy root
              _tree[i] = (_cblk+1);
          ++_cblk;                       // include the dummy root
        }
      }
      _is_ordered = true;

      //std::cout << "SCOTCH level = " << level << std::endl;
      //std::cout << "Range   Tree " << std::endl;
      //for (int i=0;i<_cblk;++i)
      //  std::cout << _range[i] << " :: " << i << " " << _tree[i] << std::endl;

      return 0;
    }

    int pruneTree(const ordinal_type cut) {
      if (cut <=0 ) return 0;

      ordinal_type_array work = ordinal_type_array("Scotch::WorkArray", _cblk+1);
      for (ordinal_type iter=0;iter<cut && _cblk > 1;++iter) {
        // horizontal merging
        {
          ordinal_type cnt = 0;
          ordinal_type parent = _tree[0];
          work[0] = cnt;
          for (ordinal_type i=1;i<_cblk;++i) {
            const ordinal_type myparent = _tree[i];
            if (myparent == parent) {
              work[i] = cnt;
            } else {
              parent = _tree[i];
              work[i] = ++cnt;
            }
          }
          work[_cblk] = ++cnt;

          ordinal_type prev = -2;
          const ordinal_type root = _cblk - 1;
          for (ordinal_type i=0;i<root;++i) {
            const ordinal_type myparent = _tree[i];
            const ordinal_type me = work[i];

            _tree[me] = work[myparent];
            if (prev != me) {
              _range[me] = _range[i];
              prev = me;
            }
          }
          {
            const ordinal_type me = work[root];
            _tree[me] = -1;
            _range[me] = _range[root];

            _range[work[root+1]] = _range[root+1];
            _cblk = cnt;
          }
        }

        // vertical merging
        if (_cblk == 2) {
          _tree[0] = -1;
          _range[0] = 0;
          _range[1] = _range[2];
          _cblk = 1;
        } else {
          ordinal_type cnt = 0;
          for (ordinal_type i=0;i<_cblk;++i) {
            const ordinal_type diff = _tree[i+1] - _tree[i];
            work[i] = (diff == 1 ? cnt : cnt++);
          }
          work[_cblk] = cnt;

          ordinal_type prev = -2;
          const ordinal_type root = _cblk - 1;
          for (ordinal_type i=0;i<root;++i) {
            const ordinal_type myparent = _tree[i];
            const ordinal_type me = work[i];

            _tree[me] = work[myparent];
            if (prev != me) {
              _range[me] = _range[i];
              prev = me;
            }
          }
          {
            const ordinal_type me = work[root];
            _tree[me] = -1;
            _range[me] = _range[root];

            _range[work[root+1]] = _range[root+1];
            _cblk = cnt;
          }
        }
      }

      // cleaning
      {
        for (ordinal_type i=(_cblk+1);i<_m;++i) {
          _tree[i] = 0;
          _range[i] = 0;
        }
        _tree[_cblk] = 0;
      }

      return 0;
    }

    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(4);
      os << scientific;

      os << " -- Scotch input -- " << std::endl
         << "    Base Value     = " << _base << std::endl
         << "    # of Rows      = " << _m << std::endl
         << "    # of NonZeros  = " << _nnz << std::endl;

      if (_is_ordered)
        os << " -- Ordering -- " << std::endl
           << "    CBLK   = " << _cblk << std::endl
           << "  PERM     PERI     RANG     TREE" << std::endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i)
        os << setw(w) << _perm[i] << "   "
           << setw(w) << _peri[i] << "   "
           << setw(w) << _range[i] << "   "
           << setw(w) << _tree[i] << std::endl;

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif
