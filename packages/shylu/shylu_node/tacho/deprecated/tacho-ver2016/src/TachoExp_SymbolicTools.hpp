#ifndef __TACHOEXP_SYMBOLIC_TOOLS_HPP__
#define __TACHOEXP_SYMBOLIC_TOOLS_HPP__

/// \file TachoExp_SymbolicTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    class SymbolicTools {
    public:
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array;
      typedef Kokkos::View<size_type*,host_exec_space> size_type_array;

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

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
        auto head  = Kokkos::subview(work, range_type(0*m, 1*m));
        auto next  = Kokkos::subview(work, range_type(1*m, 2*m));
        auto stack = Kokkos::subview(work, range_type(2*m, 3*m));

        //Kokkos::deep_copy(head, -1);        
        for (size_type i=0;i<head.extent(0);++i)
          head(i) = -1;

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

      // Joseph, W. H. Liu, Esmond G. Ng, and Barry W. Peyton, 
      // "On Finding Supernodes for Sparse Matrix Computations," 
      // SIAM J. Matrix Anal. Appl., Vol. 14, No. 1, pp. 242-252.
      inline
      static void
      computeSuperNodes(const ordinal_type m, 
                        const size_type_array &ap,
                        const ordinal_type_array &aj,
                        const ordinal_type_array &perm,
                        const ordinal_type_array &peri,
                        const ordinal_type_array &parent,
                        /* */ ordinal_type_array &supernodes,
                        const ordinal_type_array &work) {        
        // workspace
        auto flag  = Kokkos::subview(work, range_type(0*m, 1*m));
        auto count = Kokkos::subview(work, range_type(1*m, 2*m));
        auto prev  = Kokkos::subview(work, range_type(2*m, 3*m));

        // initialize workspace
        Kokkos::deep_copy(work, 0);

        // count # of children of parent
        for (ordinal_type i=0;i<m;++i) if (parent(i) >= 0) ++count(parent(i));

        // parent has more than a child, it becomes a supernode candidate 
        for (ordinal_type i=0;i<m;++i) if (count(i) > 1) flag(i) = true;

        // accumulated subtree sizes are stored in count
        for (ordinal_type i=0;i<m;++i) count(i) = 1;
        for (ordinal_type i=0;i<m;++i) if (parent(i) >= 0) count(parent(i)) += count(i);

        // tree leaves are also supernode candidate (not easy to understand this)
        for (ordinal_type i=0;i<m;++i) {
          const ordinal_type ii = perm(i);
          for (size_type p=ap(ii);p<ap(ii+1);++p) {
            const ordinal_type j = peri(aj(p));
            if (i < j) {
              const ordinal_type k = prev(j);
              if (k < (i-count(i)+1)) flag(i) = true;
              prev(j) = i;
            }
          }
        }

        // count # of supernodes 
        {
          ordinal_type k = 0; 
          flag(k) = true; // supernodes begin

          for (ordinal_type i=0;i<m;++i) k += flag(i);
          supernodes = ordinal_type_array("supernodes", k+1);
          
          // record supernodes 
          k = 0;
          for (ordinal_type i=0;i<m;++i) if (flag(i)) supernodes(k++) = i;
          supernodes(k) = m; // supernodes end
        }
      }

      /// Based on the symbolic factors, allocate pannels
      inline
      static void
      allocateSuperNodes(const ordinal_type m,
                         const size_type_array &ap,
                         const ordinal_type_array &aj,
                         const ordinal_type_array &supernodes,
                         const ordinal_type_array &work,
                         /* */ size_type_array &gid_super_panel_ptr,
                         /* */ ordinal_type_array &gid_super_panel_colidx,
                         /* */ size_type_array &sid_super_panel_ptr,
                         /* */ ordinal_type_array &sid_super_panel_colidx,
                         /* */ ordinal_type_array &blk_super_panel_colidx) {
        const ordinal_type numSuperNodes = supernodes.extent(0) - 1;
          
        // for each supernode
        auto clear_flag = [](const ordinal_type cnt,
                             const ordinal_type_array &id,
                             const ordinal_type_array &flag) {
          for (ordinal_type i=0;i<cnt;++i) flag(id(i)) = 0;
        };

        auto clear_array = [](const ordinal_type cnt,
                              const ordinal_type_array &a) {
          memset(a.data(), 0, cnt*sizeof(typename ordinal_type_array::value_type));
        };

        auto colmap_per_supernode = [](const size_type_array &ap,
                                       const ordinal_type_array &aj,
                                       const ordinal_type sbeg, 
                                       const ordinal_type send,
                                       const ordinal_type_array &cid,
                                       const ordinal_type_array &flag) -> ordinal_type {
          // # of columns accessed by this super node
          ordinal_type cnt = 0;
          
          // loop over super node cols (diagonal block)
          for (ordinal_type col=sbeg;col<send;++col) {
            flag(col) = true; // visitation flag on
            cid(cnt++) = col; // record the column indicies
          }
          
          // visit each super node row (off diagonal block)
          for (ordinal_type i=sbeg;i<send;++i) {
            for (size_type j=ap(i);j<ap(i+1);++j) {
              const ordinal_type col = aj(j);
              if (flag(col)) {
                // already visited; pass on
              } else {
                flag(col) = true;
                cid(cnt++) = col;
              }
            }
          }
          return cnt;
        };

        auto sidmap_per_supernode = [](const ordinal_type ndofs,
                                       const ordinal_type_array &cid,
                                       const ordinal_type_array &sid_colored_in_rows,
                                       const ordinal_type_array &sids_connected_to_this_row,
                                       const ordinal_type_array &count,
                                       const ordinal_type_array &blks_connected_to_this_row) -> ordinal_type {
          ordinal_type cnt = 0;
          for (ordinal_type k=0;k<ndofs;++k) {
            const ordinal_type sid = sid_colored_in_rows(cid(k));
            if (count(sid) == 0) 
              sids_connected_to_this_row(cnt++) = sid;
            ++count(sid);
          }
          if (blks_connected_to_this_row.data() != NULL) {
            for (ordinal_type k=0;k<cnt;++k) {
              const ordinal_type sid = sids_connected_to_this_row(k);
              blks_connected_to_this_row(k) = count(sid);
            }
          }
          return cnt;
        };
        
        auto flag = Kokkos::subview(work, range_type(0*m, 1*m));
        auto cid  = Kokkos::subview(work, range_type(1*m, 2*m));
        auto rid  = Kokkos::subview(work, range_type(2*m, 3*m));
        auto tmp  = Kokkos::subview(work, range_type(3*m, 4*m));
        
        // zeros
        Kokkos::deep_copy(work, 0);

        ///
        /// super color in rows
        ///
        for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
          const ordinal_type sbeg = supernodes(sid), send = supernodes(sid+1);
          for (ordinal_type i=sbeg;i<send;++i) 
            rid(i) = sid;
        }
        
        ///
        /// super_panel_ptr, colidx_super_panel
        ///

        const ordinal_type_array null_ordinal_type_array("null_ordinal_type_array");

        /// count the # of associated columns to all supernodes
        gid_super_panel_ptr = size_type_array("gid_super_panel_ptr", numSuperNodes+1);
        sid_super_panel_ptr= size_type_array("sid_super_panel_ptr", numSuperNodes+1);

        for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
          const ordinal_type sbeg = supernodes(sid), send = supernodes(sid+1);

          const ordinal_type ndofs = colmap_per_supernode(ap, aj, sbeg, send, cid, flag);
          gid_super_panel_ptr(sid+1) = gid_super_panel_ptr(sid) + ndofs;
          clear_flag(ndofs, cid, flag);

          const ordinal_type nsids = sidmap_per_supernode(ndofs, cid, rid, tmp, flag,
                                                          null_ordinal_type_array);
          sid_super_panel_ptr(sid+1) = sid_super_panel_ptr(sid) + nsids + 1;
          clear_flag(nsids, tmp, flag); 

          clear_array(ndofs, cid);
          clear_array(nsids, tmp);
        }
        gid_super_panel_colidx = ordinal_type_array("gid_super_panel_colidx", gid_super_panel_ptr(numSuperNodes));
        sid_super_panel_colidx = ordinal_type_array("sid_super_panel_colidx", sid_super_panel_ptr(numSuperNodes));
        blk_super_panel_colidx = ordinal_type_array("blk_super_panel_colidx", sid_super_panel_ptr(numSuperNodes));
        
        /// sort and column indices per supernode
        for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
          const ordinal_type sbeg = supernodes(sid), send = supernodes(sid+1);

          // *** sort and construct gid map to sparse matrix
          const auto gid_range = range_type(gid_super_panel_ptr(sid), gid_super_panel_ptr(sid+1));
        
          auto colidx = Kokkos::subview(gid_super_panel_colidx, gid_range);
          const ordinal_type ndofs = colmap_per_supernode(ap, aj, sbeg, send, colidx, flag);

          std::sort(colidx.data() + (send-sbeg), colidx.data() + ndofs);
          clear_flag(ndofs, colidx, flag);          
          
          // *** sort associated supernodes to sid; stored in the work array cid
          const auto sid_range = range_type(sid_super_panel_ptr(sid), sid_super_panel_ptr(sid+1));
        
          auto sididx = Kokkos::subview(sid_super_panel_colidx, sid_range);
          auto blkidx = Kokkos::subview(blk_super_panel_colidx, sid_range);

          // cid and tmp are used for temporary storage
          const ordinal_type nsids = sidmap_per_supernode(ndofs, colidx, rid, cid, flag, tmp);

          // flag is now used for permutation vector
          clear_flag(nsids, cid, flag); 
          for (ordinal_type i=0;i<nsids;++i) flag(i) = i;

          std::sort(flag.data(), flag.data()+nsids, 
                    [&](std::size_t i, std::size_t j){ return sididx(i) < sididx(j); });
          std::transform(flag.data(), flag.data()+nsids, 
                         sididx.data(),[&](std::size_t i) { return cid(i); });
          std::transform(flag.data(), flag.data()+nsids, 
                         cid.data(),[&](std::size_t i) { return tmp(i); });

          for (ordinal_type i=0;i<nsids;++i) blkidx(i+1) = blkidx(i) + cid(i);

          clear_array(nsids, cid);
          clear_array(nsids, flag);
        }
      }

      /// construct tree explicitly
      inline
      static void
      computeSuperNodesAssemblyTree(const ordinal_type_array &parent,
                                    const ordinal_type_array &supernodes,
                                    /* */ size_type_array &stree_ptr,
                                    /* */ ordinal_type_array &stree_children,
                                    /* */ ordinal_type_array &stree_roots,
                                    const ordinal_type_array &work) {
        const ordinal_type numSuperNodes = supernodes.extent(0) - 1;
        const ordinal_type m = supernodes(numSuperNodes);

        auto sparent = Kokkos::subview(work, range_type(0*m, 1*m));        
        auto flag    = Kokkos::subview(work, range_type(1*m, 2*m));        

        // color flag with supernodes (for the ease to detect supernode id from dofs)
        for (ordinal_type i=0;i<numSuperNodes;++i) 
          for (ordinal_type j=supernodes(i);j<supernodes(i+1);++j) 
            flag(j) = i;

        // coarse parent into sparent
        for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
          sparent(sid) = -1;
          for (ordinal_type i=supernodes(sid);i<supernodes(sid+1);++i) {
            if (parent(i) >= 0) {
              const ordinal_type sidpar = flag(parent(i));
              if (sidpar != sid) sparent(sid) = sidpar;
            }
          }
        }

        auto clear_array = [](const ordinal_type cnt,
                              const ordinal_type_array &a) {
          memset(a.data(), 0, cnt*sizeof(typename ordinal_type_array::value_type));
        };

        // construct parent - child relations
        {
          clear_array(m, flag);
          ordinal_type cnt = 0;
          for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
            const ordinal_type sidpar = sparent(sid); 
            if (sidpar == -1) 
              ++cnt;
            else
              ++flag(sparent(sid));
          }
          stree_roots = ordinal_type_array("stree_roots", cnt);
        }

        // prefix scan
        {
          stree_ptr = size_type_array("stree_ptr", numSuperNodes + 1);
          for (ordinal_type sid=0;sid<numSuperNodes;++sid)
            stree_ptr(sid+1) = stree_ptr(sid) + flag(sid); 
        }

        {
          clear_array(numSuperNodes, flag);
          ordinal_type cnt = 0;
          stree_children = ordinal_type_array("stree_children", stree_ptr(numSuperNodes));
          for (ordinal_type sid=0;sid<numSuperNodes;++sid) {
            const ordinal_type sidpar = sparent(sid);                       
            if (sidpar == -1) 
              stree_roots(cnt++) = sid;
            else
              stree_children(stree_ptr(sidpar)+flag(sidpar)++) = sid;
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

      // supernodes output
      ordinal_type_array _supernodes;      

      // dof mapping to sparse matrix
      size_type_array _gid_super_panel_ptr;
      ordinal_type_array _gid_super_panel_colidx;

      // supernode map and panel size configuration
      size_type_array _sid_super_panel_ptr;
      ordinal_type_array _sid_super_panel_colidx, _blk_super_panel_colidx;

      // supernode elimination tree (parent - children)
      size_type_array _stree_ptr;
      ordinal_type_array _stree_children, _stree_roots;
      
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
      
      template<typename CrsMatBaseType, typename GraphToolType>
      SymbolicTools(CrsMatBaseType &A,
                    GraphToolType &G) {
        _m = A.NumRows();
        
        _ap   = Kokkos::create_mirror_view(typename host_exec_space::memory_space(), A.RowPtr());
        _aj   = Kokkos::create_mirror_view(typename host_exec_space::memory_space(), A.Cols());
        _perm = Kokkos::create_mirror_view(typename host_exec_space::memory_space(), G.PermVector());
        _peri = Kokkos::create_mirror_view(typename host_exec_space::memory_space(), G.InvPermVector());

        Kokkos::deep_copy(_ap, A.RowPtr());
        Kokkos::deep_copy(_aj, A.Cols());
        Kokkos::deep_copy(_perm, G.PermVector());
        Kokkos::deep_copy(_peri, G.InvPermVector());
      }

      inline
      ordinal_type NumSuperNodes() const { return _supernodes.extent(0) - 1; }
      
      inline
      ordinal_type_array SuperNodes() const { return _supernodes; }
      
      inline 
      size_type_array gidSuperPanelPtr() const { return _gid_super_panel_ptr; }

      inline 
      ordinal_type_array gidSuperPanelColIdx() const { return _gid_super_panel_colidx; }

      inline 
      size_type_array sidSuperPanelPtr() const { return _sid_super_panel_ptr; }

      inline 
      ordinal_type_array sidSuperPanelColIdx() const { return _sid_super_panel_colidx; }

      inline 
      ordinal_type_array blkSuperPanelColIdx() const { return _blk_super_panel_colidx; }
      
      inline 
      size_type_array SuperNodesTreePtr() const { return _stree_ptr; }

      inline 
      ordinal_type_array SuperNodesTreeChildren() const { return _stree_children; }
      
      inline 
      ordinal_type_array SuperNodesTreeRoots() const { return _stree_roots; }
      
      inline
      void
      symbolicFactorize() {
        // compute elimination tree
        ordinal_type_array work("work", _m*4);
        ordinal_type_array parent("parent", _m);
        {
          auto post = Kokkos::subview(work, range_type(0*_m, 1*_m));
          auto w    = Kokkos::subview(work, range_type(1*_m, 4*_m));
          
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
        computeSuperNodes(_m, _ap, _aj, _perm, _peri, parent, _supernodes, work);

        // compute fill pattern
        size_type_array ap;
        ordinal_type_array aj;
        computeFillPatternUpper(_m, _ap, _aj, _perm, _peri, ap, aj, work);          

        // allocate supernodes
        allocateSuperNodes(_m, ap, aj, _supernodes, work, 
                           _gid_super_panel_ptr,
                           _gid_super_panel_colidx,
                           _sid_super_panel_ptr,
                           _sid_super_panel_colidx,
                           _blk_super_panel_colidx);

        // supernode assembly tree
        computeSuperNodesAssemblyTree(parent,                                                 
                                      _supernodes,                                             
                                      _stree_ptr,                                                 
                                      _stree_children,
                                      _stree_roots,
                                      work);              
      }            
    };

  }
}
#endif
