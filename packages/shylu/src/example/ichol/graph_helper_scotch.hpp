#pragma once
#ifndef __GRAPH_HELPER_SCOTCH_HPP__
#define __GRAPH_HELPER_SCOTCH_HPP__

#include "scotch.h"
#include "util.hpp"

namespace Example { 
  
  using namespace std;

  template<class CrsMatrix>
  class GraphHelper_Scotch {
  public:
    typedef typename CrsMatrix::ordinal_type ordinal_type;
    typedef typename CrsMatrix::size_type    size_type;

  private:
    // scotch main data structure
    SCOTCH_Graph _graph;

    // scotch input wihtout diagonal contributions
    ordinal_type _base,_m,*_cidx;
    size_type _nnz,*_rptr;

    // scotch output 
    ordinal_type _cblk,*_perm,*_peri,*_range,*_tree;

    // status flag
    bool _is_ordered;

  public:
    GraphHelper_Scotch(CrsMatrix& A) {

      _is_ordered = false;
      _cblk  = 0;

      // scotch does not allow self-contribution (diagonal term in sparse matrix)
      _base  = A.BaseVal();
      _m     = A.NumRows();
      _nnz   = A.NumNonZeros();       

      _rptr  = new size_type[_m+1]();    
      _cidx  = new ordinal_type[_nnz]();

      _perm  = new ordinal_type[_m]();
      _peri  = new ordinal_type[_m]();
      _range = new ordinal_type[_m+1]();
      _tree  = new ordinal_type[_m]();

      // adjust graph structure
      A.convertGraph(_nnz, _rptr, _cidx);

      int ierr = 0;
      
      ierr = SCOTCH_graphInit(&_graph);CHKERR(ierr);
      ierr = SCOTCH_graphBuild(&_graph,             // scotch graph
                               _base,               // base value
                               _m,                  // # of vertices
                               _rptr,               // column index array pointer begin
                               _rptr+1,             // column index array pointer end
                               NULL,                // weights on vertices (optional)
                               NULL,                // label array on vertices (optional)
                               _nnz,                // # of nonzeros
                               _cidx,               // column index array
                               NULL);CHKERR(ierr);  // edge load array (optional)
      ierr = SCOTCH_graphCheck(&_graph);CHKERR(ierr);
    }
    virtual~GraphHelper_Scotch() {
      SCOTCH_graphFree(&_graph);
      delete _rptr,_cidx,_perm,_peri,_range,_tree;
    }
    ordinal_type getNumBlocks() const {
      return _cblk; 
    }

    // ordinal_type* getPermutationVector() const {
    //   return _perm;
    // }
    // ordinal_type* getInversePermutationVector() const {
    //   return _peri;
    // }

    int computeOrdering() {
      int ierr = 0, level = log2(_nnz)+10;
        
      SCOTCH_Strat stradat;
      SCOTCH_Num straval = (SCOTCH_STRATLEVELMAX   | 
                            SCOTCH_STRATLEVELMIN   | 
                            SCOTCH_STRATLEAFSIMPLE | 
                            SCOTCH_STRATSEPASIMPLE);

      ierr = SCOTCH_stratInit(&stradat);CHKERR(ierr);
      ierr = SCOTCH_stratGraphOrderBuild (&stradat, straval, level, 0.2);CHKERR(ierr);

      ierr = SCOTCH_graphOrder(&_graph, 
                               &stradat, 
                               _perm, 
                               _peri,
                               &_cblk,
                               _range,
                               _tree);CHKERR(ierr);
      SCOTCH_stratExit(&stradat);

      _is_ordered = true;

      return 0;
    }

    int constructTree() {
    }

    int showMe(ostream &os) {
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;
      
      os << " -- Scotch input -- " << endl
         << "    Base Value     = " << _base << endl
         << "    # of Rows      = " << _m << endl
         << "    # of NonZeros  = " << _nnz << endl;
      
      // for (ordinal_type i=0;i<_m;++i) {
      //   size_type jbegin = _rptr[i], jend = _rptr[i+1];
      //   os << endl;
      //   for (size_type j=jbegin;j<jend;++j)
      //     os << i << "  " << _cidx[j] << endl;                                              
      // } 
      // os << endl;
      
      if (_is_ordered) 
        os << " -- Elimination tree -- " << endl
           << "    CBLK   = " << _cblk << endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i) 
        os << setw(w) << _perm[i] << "   " 
           << setw(w) << _peri[i] << "   " 
           << setw(w) << _range[i] << "   "
           << setw(w) << _tree[i] << endl; 
        
      os.unsetf(ios::scientific);
      os.precision(prec);    
    }
    
  };

}

#endif
