#pragma once
#ifndef __GRAPH_HELPER_SCOTCH_HPP__
#define __GRAPH_HELPER_SCOTCH_HPP__

/// \file graph_helper_scotch.hpp
/// \brief Interface to scotch reordering
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "scotch.h"
#include "util.hpp"

namespace Example {

  using namespace std;

  template<class CrsMatrixType>
  class GraphHelper_Scotch : public Disp {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::size_type    size_type;

    typedef typename CrsMatrixType::ordinal_type_array ordinal_type_array;
    typedef typename CrsMatrixType::size_type_array    size_type_array;

  private:
    string _label;

    // scotch main data structure
    SCOTCH_Graph _graph;

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

    void setLabel(string label) { _label = label; }
    string Label() const { return _label; }

    size_type NumNonZeros() const { return _nnz; }
    ordinal_type NumRows() const { return _m; }

    ordinal_type_array PermVector()    const { return _perm; }
    ordinal_type_array InvPermVector() const { return _peri; }

    ordinal_type_array RangeVector()   const { return _range; }
    ordinal_type_array TreeVector()    const { return _tree; }

    ordinal_type NumBlocks() const { return _cblk; }

    GraphHelper_Scotch(const CrsMatrixType& A, 
                       const int seed = GraphHelper::DefaultRandomSeed) {

      _label = "GraphHelper_Scotch::" + A.Label();

      _is_ordered = false;
      _cblk  = 0;

      // scotch does not allow self-contribution (diagonal term in sparse matrix)
      _base  = 0; //A.BaseVal();
      _m     = A.NumRows();
      _nnz   = A.NumNonZeros();

      _rptr  = size_type_array(_label+"::RowPtrArray", _m+1);
      _cidx  = ordinal_type_array(_label+"::ColIndexArray", _nnz);

      _perm  = ordinal_type_array(_label+"::PermutationArray", _m);
      _peri  = ordinal_type_array(_label+"::InvPermutationArray", _m);
      _range = ordinal_type_array(_label+"::RangeArray", _m);
      _tree  = ordinal_type_array(_label+"::TreeArray", _m);

      // create a graph structure without diagonals
      A.convertGraph(_nnz, _rptr, _cidx);

      int ierr = 0;
      ordinal_type *rptr = reinterpret_cast<ordinal_type*>(_rptr.ptr_on_device());
      ordinal_type *cidx = reinterpret_cast<ordinal_type*>(_cidx.ptr_on_device());

      if (seed != GraphHelper::DefaultRandomSeed) {
        SCOTCH_randomSeed(seed);
        SCOTCH_randomReset();
      }

      ierr = SCOTCH_graphInit(&_graph);CHKERR(ierr);
      ierr = SCOTCH_graphBuild(&_graph,             // scotch graph
                               _base,               // base value
                               _m,                  // # of vertices
                               rptr,                // column index array pointer begin
                               rptr+1,              // column index array pointer end
                               NULL,                // weights on vertices (optional)
                               NULL,                // label array on vertices (optional)
                               _nnz,                // # of nonzeros
                               cidx,                // column index array
                               NULL);CHKERR(ierr);  // edge load array (optional)
      ierr = SCOTCH_graphCheck(&_graph);CHKERR(ierr);
    }
    virtual~GraphHelper_Scotch() {
      SCOTCH_graphFree(&_graph);
    }

    int computeOrdering(const int treecut = 15,
                        const int minblksize = 0) {
      int ierr = 0;

      // pointers for global graph ordering
      ordinal_type *perm  = _perm.ptr_on_device();
      ordinal_type *peri  = _peri.ptr_on_device();
      ordinal_type *range = _range.ptr_on_device();
      ordinal_type *tree  = _tree.ptr_on_device();
      
      {
        const int level = max(1, int(log2(_m)-treecut)); // level = log2(_nnz)+10;
        SCOTCH_Strat stradat;
        SCOTCH_Num straval = (SCOTCH_STRATLEVELMAX   |
                              SCOTCH_STRATLEVELMIN   |
                              SCOTCH_STRATLEAFSIMPLE |
                              SCOTCH_STRATSEPASIMPLE);
        
        ierr = SCOTCH_stratInit(&stradat);CHKERR(ierr);
        ierr = SCOTCH_stratGraphOrderBuild (&stradat, straval, level, 0.2);CHKERR(ierr);
        
        ierr = SCOTCH_graphOrder(&_graph,
                                 &stradat,
                                 perm,
                                 peri,
                                 &_cblk,
                                 range,
                                 tree);CHKERR(ierr);
        SCOTCH_stratExit(&stradat);
      }

      // provided blksize is greater than 0, reorder internally
      if (treecut != 0 && minblksize > 0) {
        // graph array
        ordinal_type *rptr = reinterpret_cast<ordinal_type*>(_rptr.ptr_on_device());
        ordinal_type *cidx = reinterpret_cast<ordinal_type*>(_cidx.ptr_on_device());

        // create workspace in
        size_type_array    rptr_work = size_type_array(_label+"::Block::RowPtrArray", _m+1);
        ordinal_type_array cidx_work = ordinal_type_array(_label+"::Block::ColIndexArray", _nnz);

        // create workspace output
        ordinal_type_array perm_work  = ordinal_type_array(_label+"::Block::PermutationArray", _m);
        ordinal_type_array peri_work  = ordinal_type_array(_label+"::Block::InvPermutationArray", _m);
        ordinal_type_array range_work = ordinal_type_array(_label+"::Block::RangeArray", _m);
        ordinal_type_array tree_work  = ordinal_type_array(_label+"::Block::TreeArray", _m);
        
        // scotch input
        ordinal_type *rptr_blk = reinterpret_cast<ordinal_type*>(rptr_work.ptr_on_device());
        ordinal_type *cidx_blk = reinterpret_cast<ordinal_type*>(cidx_work.ptr_on_device());
      
        size_type nnz = 0;
        rptr_blk[0] = nnz;        

        for (ordinal_type iblk=0;iblk<_cblk;++iblk) {
          // allocate graph
          SCOTCH_Graph graph;
          
          ierr = SCOTCH_graphInit(&graph);CHKERR(ierr);
          
          SCOTCH_Strat stradat;
          SCOTCH_Num straval = (/*SCOTCH_STRATLEVELMAX   |
                                  SCOTCH_STRATLEVELMIN   |*/
                                SCOTCH_STRATLEAFSIMPLE |
                                SCOTCH_STRATSEPASIMPLE);
          
          ierr = SCOTCH_stratInit(&stradat);CHKERR(ierr);
          ierr = SCOTCH_stratGraphOrderBuild(&stradat, straval, 0, 0.2);CHKERR(ierr);
          
          const ordinal_type ibegin = range[iblk], iend = range[iblk+1], m = iend - ibegin;

          // scotch output
          ordinal_type cblk_blk = 0;
          
          ordinal_type *perm_blk  = perm_work.ptr_on_device()  + ibegin;
          ordinal_type *peri_blk  = peri_work.ptr_on_device()  + ibegin;
          ordinal_type *range_blk = range_work.ptr_on_device() + ibegin;
          ordinal_type *tree_blk  = tree_work.ptr_on_device()  + ibegin;
          
          // if each blk is greater than the given minblksize, reorder internally
          if (m < minblksize) {
            for (int i=ibegin;i<iend;++i) {
              const ordinal_type ii = peri[i];
              const ordinal_type jbegin = rptr[ii];
              const ordinal_type jend = rptr[ii+1];
              
              for (int j=jbegin;j<jend;++j) {
                const ordinal_type jj = perm[cidx[j]];
                if (ibegin <= jj && jj < iend)  
                  cidx_blk[nnz++] = (jj - ibegin);
              }
              rptr_blk[i+1] = nnz;
            }
            const size_type nnz_blk = nnz - rptr_blk[ibegin];

            ierr = SCOTCH_graphBuild(&graph,             // scotch graph
                                     0,                  // base value
                                     m,                  // # of vertices
                                     &rptr_blk[ibegin],  // column index array pointer begin
                                     &rptr_blk[ibegin]+1,// column index array pointer end
                                     NULL,               // weights on vertices (optional)
                                     NULL,               // label array on vertices (optional)
                                     nnz_blk,            // # of nonzeros
                                     cidx_blk,           // column index array
                                     NULL);CHKERR(ierr); // edge load array (optional)
            ierr = SCOTCH_graphCheck(&graph);CHKERR(ierr);
            ierr = SCOTCH_graphOrder(&graph,
                                     &stradat,
                                     perm_blk,
                                     peri_blk,
                                     &cblk_blk,
                                     range_blk,
                                     tree_blk);CHKERR(ierr);
          } else {
            for (ordinal_type i=0;i<m;++i) {
              perm_blk[i] = i;
              peri_blk[i] = i;
            }
            range_blk[1] = m;
            tree_blk[0] = -1;
          }

          SCOTCH_stratExit(&stradat);
          SCOTCH_graphFree(&graph);

          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type ii = peri_blk[i] + ibegin;
            peri_blk[i] = peri[ii];
          }
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type ii = i + ibegin;
            peri[ii] = peri_blk[i];
          }
          
        }

        for (ordinal_type i=0;i<_m;++i) 
          perm[peri[i]] = i;
      }
    
      _is_ordered = true;

      //cout << "SCOTCH level = " << level << endl;
      //cout << "Range   Tree " << endl;
      //for (int i=0;i<_cblk;++i)
      //  cout << _range[i] << " :: " << i << " " << _tree[i] << endl;

      return 0;
    }

    ostream& showMe(ostream &os) const {
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;

      os << " -- Scotch input -- " << endl
         << "    Base Value     = " << _base << endl
         << "    # of Rows      = " << _m << endl
         << "    # of NonZeros  = " << _nnz << endl;

      if (_is_ordered)
        os << " -- Elimination tree -- " << endl
           << "    CBLK   = " << _cblk << endl
           << "  PERM     PERI     RANG     TREE" << endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i)
        os << setw(w) << _perm[i] << "   "
           << setw(w) << _peri[i] << "   "
           << setw(w) << _range[i] << "   "
           << setw(w) << _tree[i] << endl;

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif
