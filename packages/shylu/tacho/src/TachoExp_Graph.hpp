#ifndef __TACHOEXP_GRAPH_HPP__
#define __TACHOEXP_GRAPH_HPP__

/// \file TachoExp_Graph.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {
    
    ///
    /// graph reordering happens in host space only
    ///
    class Graph {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array;

      // scotch use int and long is misinterpreted in pointers
      typedef Kokkos::View<ordinal_type*,host_exec_space> size_type_array;

    private:
      ordinal_type _m;
      size_type _nnz;

      size_type_array _rptr;
      ordinal_type_array _cidx;

      template<typename SizeTypeArray,
               typename OrdinalTypeArray>
      inline
      void
      init(const ordinal_type m, 
           const size_type nnz, 
           const SizeTypeArray &ap,
           const OrdinalTypeArray &aj) {
        _rptr = size_type_array("Graph::rptr", m + 1);
        _cidx = ordinal_type_array("Graph::cidx", nnz);

        _m = m;
        _nnz = 0;
        for (ordinal_type i=0;i<_m;++i) {
          const size_type jbeg = ap(i), jend = ap(i+1);
          _rptr(i) = _nnz;
          for (size_type j=jbeg;j<jend;++j) {
            // skip diagonal 
            const ordinal_type col = aj(j);
            if (i != col) _cidx(_nnz++) = col;
          }
        }
        _rptr(_m) = _nnz;        
      }

    public:

      Graph() = default;
      Graph(const Graph &b) = default;

      template<typename SizeTypeArray,
               typename OrdinalTypeArray>
      Graph(const ordinal_type m, 
            const size_type nnz, 
            const SizeTypeArray &ap,
            const OrdinalTypeArray &aj) {
        init(m, nnz, ap, aj);
      }
      
      template<typename ValueType, typename SpaceType>
      inline
      Graph(const CrsMatrixBase<ValueType,SpaceType> &A) {
        //
        // host mirroring
        //
        CrsMatrixBase<ValueType,host_exec_space> AA;
        AA.createMirror(A);
        AA.copy(A);
        
        init(AA.NumRows(), AA.NumNonZeros(), AA.RowPtr(), AA.Cols());
      }

      inline
      size_type_array RowPtr() const { return _rptr; }
      
      inline
      ordinal_type_array ColIdx() const { return _cidx; }

      inline
      ordinal_type NumRows() const { return _m; }
      
      inline
      ordinal_type NumNonZeros() const { return _nnz; }

      inline
      void clear() {
        _m = 0; 
        _nnz = 0;
        _rptr = size_type_array();
        _cidx = ordinal_type_array();
      }

      std::ostream& showMe(std::ostream &os, const bool detail = false) const {
        std::streamsize prec = os.precision();
        os.precision(4);
        os << std::scientific;

        os << " -- Graph -- " << std::endl
           << "    # of Rows      = " << _m << std::endl
           << "    # of NonZeros  = " << _nnz << std::endl;

        const int w = 10;
        if (detail) {
          os << std::setw(w) <<  "Row" << "  "
             << std::setw(w) <<  "Col" 
             << std::endl;
          for (ordinal_type i=0;i<_m;++i) {
            const size_type jbeg = _rptr[i], jend = _rptr[i+1];
            for (size_type j=jbeg;j<jend;++j) {
              os << std::setw(w) <<        i << "  "
                 << std::setw(w) << _cidx[j] << "  "
                 << std::endl;
            }
          }
        }

        os.unsetf(std::ios::scientific);
        os.precision(prec);
        
        return os;
      }


    };

  }
}


#endif
