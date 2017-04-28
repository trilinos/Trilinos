#ifndef __TACHOEXP_SYMBOLIC_TOOLS_HPP__
#define __TACHOEXP_SYMBOLIC_TOOLS_HPP__

/// \file Tacho_CrsMatrixTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    class SymbolicTools {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array;
      typedef Kokkos::View<size_type*,host_exec_space> size_type_array;

      // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 42.
      inline
      static void
      computeEliminationTree(const ordinal_type m,
                             const size_type_array &ap,
                             const ordinal_type_array &aj,
                             const ordinal_type_array &perm,
                             const ordinal_type_array &peri,
                             const ordinal_type_array &parent,
                             const ordinal_type_array &ancestor) {
        for (ordinal_type i=0;i<m;++i) {
          parent(i) = -1;
          ancestor(i) = -1;

          const ordinal_type idx = perm(i);
          for (size_type p=ap(idx);p<ap(idx+1);++p) {
            ordinal_type j = peri(aj(p));
            for ( ;j!=-1 && j<i;) {
              const ordinal_type next = ancestor(j);
              ancestor(j) = i;
              if (next == -1) parent(j) = i;
              j = next;
            }
          }
        }
      }

      // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 45.
      inline
      static ordinal_type
      TreeDepthFirstSearch(const ordinal_type j,
                           const ordinal_type c,
                           const ordinal_type_array &head,
                           const ordinal_type_array &next,
                           const ordinal_type_array &post,
                           const ordinal_type_array &stack) {
        ordinal_type top = 0, k = c;
        stack(top) = j;
        while (top >= 0) {
          const ordinal_type p = stack(top);
          const ordinal_type i = head(p);
          if (i == -1) {
            --top;
            post(k++) = p;
          } else {
            head(p) = next(i);
            stack(++top) = i;
          }
        }
        return k;
      }

      // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 45.
      inline
      static void
      computePostOrdering(const ordinal_type m,
                          const ordinal_type_array &parent,
                          const ordinal_type_array &post,
                          const ordinal_type_array &work) {
        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
        auto head  = Kokkos::subview(work, range_type(0*m, 1*m));
        auto next  = Kokkos::subview(work, range_type(1*m, 2*m));
        auto stack = Kokkos::subview(work, range_type(2*m, 3*m));
        
        Kokkos::deep_copy(head, -1);

        for (ordinal_type i=m-1;i>=0;--i) {
          if (parent(i) == -1) continue;
          next(i) = head(parent(i));;
          head(parent(i)) = i;
        }
        ordinal_type k = 0;
        for (ordinal_type i=0;i<m;++i) {
          if (parent(i) != -1) continue;
          k = TreeDepthFirstSearch(i, k, head, next, post, stack);
        }
      }

      // Tim Davis, Algorithm 849: A Concise Sparse Cholesky Factorization Package
      // ACM TOMS Vol 31 No. 4 pp 587--591.
      inline
      static void
      computeFillPatternUpper(const ordinal_type m, 
                              const size_type_array &ap,
                              const ordinal_type_array &aj,
                              const ordinal_type_array &perm,
                              const ordinal_type_array &peri,
                              /* */ size_type_array &up,
                              /* */ ordinal_type_array &uj,
                              const ordinal_type_array &work) {        
        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        auto parent        = Kokkos::subview(work, range_type(0*m, 1*m));
        auto flag          = Kokkos::subview(work, range_type(1*m, 2*m));
        auto upper_row_cnt = Kokkos::subview(work, range_type(2*m, 3*m));

        // count nnz per row
        for (ordinal_type i=0;i<m;++i) {
          parent(i) = -1;
          flag(i) = -1;
          upper_row_cnt(i) = 1; // add diagonal entry

          const ordinal_type ii = perm(i);          
          for (size_type p=ap(ii);p<ap(ii+1);++p) {
            for (ordinal_type j=peri(aj(p));flag(j)!=i && j<i;j=parent(j)) {
              if (parent(j) == -1) parent(j) = i;
              ++upper_row_cnt(j);
              flag(j) = i;
            }
          }
        }

        // prefix scan
        up = size_type_array("up", m+1);
        for (ordinal_type i=0;i<m;++i)
          up(i+1) = up(i) + upper_row_cnt(i);
        
        // fill-in
        uj = ordinal_type_array("uj", up(m));
        Kokkos::deep_copy(uj, -1);
        for (ordinal_type i=0;i<m;++i) {
          parent(i) = -1;
          flag(i) = -1;

          // diagonal entry
          upper_row_cnt(i) = 1; 
          uj(up(i)) = i;

          const ordinal_type ii = perm(i);          
          for (size_type p=ap(ii);p<ap(ii+1);++p) {
            for (ordinal_type j=peri(aj(p));flag(j)!=i && j<i;j=parent(j)) {
              if (parent(j) == -1) parent(j) = i;
              uj(up(j)+upper_row_cnt(j)) = i;
              ++upper_row_cnt(j);
              flag(j) = i;
            }
          }
        }
      }

    private:
      // matrix input
      ordinal_type _m;
      size_type_array _ap;
      ordinal_type_array _aj;

      // graph ordering input
      ordinal_type_array _perm, _peri;

    public:
      SymbolicTools() = default;
      SymbolicTools(const SymbolicTools &b) = default;

      ///
      /// construction
      ///      
      SymbolicTools(const ordinal_type m,
                    const size_type_array &ap,
                    const ordinal_type_array &aj,
                    const ordinal_type_array &perm,
                    const ordinal_type_array &peri) 
        : _m(m),
          _ap(ap),
          _aj(aj),
          _perm(perm),
          _peri(peri) {}

      inline
      void
      SymbolicFactorize(size_type_array &ap, 
                        ordinal_type_array &aj) {
        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        // compute elimination tree
        ordinal_type_array work("work", _m*4);
        ordinal_type_array parent("parent", _m);
        {
          auto post   = Kokkos::subview(work, range_type(0*_m, 1*_m));
          auto w      = Kokkos::subview(work, range_type(1*_m, 4*_m));
          
          // compute elimination tree and its post ordering
          computeEliminationTree(_m, _ap, _aj, _perm, _peri, parent, w);
          computePostOrdering(_m, parent, post, w);
          
          // update permutation vector
          for (ordinal_type i=0;i<_m;++i) w(i) = _perm(i);
          for (ordinal_type i=0;i<_m;++i) _perm(i) = w(post(i));
          for (ordinal_type i=0;i<_m;++i) _peri(_perm(i)) = i;
          
          // update elimination tree according to the updated permutation vector
          computeEliminationTree(_m, _ap, _aj, _perm, _peri, parent, w);
        }

        // compute super node structure
        //computeSuperNodesStructure();

        computeFillPatternUpper(_m, _ap, _aj, _perm, _peri, ap, aj, work);          
      }
        

    };

  }
}
#endif
