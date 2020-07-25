#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho.hpp"
#include "Tacho_Internal.hpp"
#include "Tacho_CommandLineParser.hpp"

using namespace Tacho;

using host_space = Kokkos::DefaultHostExecutionSpace;
using size_type_1d_view_type = Kokkos::View<size_type*,host_space>;
using ordinal_type_1d_view_type = Kokkos::View<ordinal_type*,host_space>;

int readGraphFile(const std::string &filename,
		  size_type_1d_view_type &ap,
		  ordinal_type_1d_view_type &aj) {
  std::ifstream in;
  in.open(filename);
  if (in.good()) {
    ordinal_type num_nodes(0);
    in >> num_nodes;
    ap = size_type_1d_view_type("rowptr", num_nodes+1);
    for (ordinal_type i=0;i<ordinal_type(num_nodes+1);++i)
      in >> ap(i);

    aj = ordinal_type_1d_view_type("colidx", ap(num_nodes));
    for (ordinal_type i=0;i<ordinal_type(num_nodes);++i)
      for (ordinal_type j=ap(i);j<ap(i+1);++j)
	in >> aj(j);
  } else {
    std::cout <<" Failed to open the file: " << filename << std::endl;
    return -1;
  }
  return 0;
}

int readWeightFile(const std::string &filename,
		   ordinal_type_1d_view_type &aw) {
  std::ifstream in;
  in.open(filename);
  if (in.good()) {
    ordinal_type total_dofs(0), num_nodes(0);
    in >> total_dofs;
    in >> num_nodes;
    aw = ordinal_type_1d_view_type("weight", num_nodes);
    for (ordinal_type i=0;i<num_nodes;++i)
      in >> aw(i);
  } else {
    std::cout <<" Failed to open the file: " << filename << std::endl;
    return -1;
  }
  return 0;
}

int doSymbolicFactorization(const ordinal_type m,
                            const size_type_1d_view_type &ap, 
                            const ordinal_type_1d_view_type &aj,
                            const ordinal_type_1d_view_type &aw,
                            const bool verbose,
                            const bool is_condensed,
                            const bool evaporate) {
  Kokkos::Impl::Timer timer;
  double t = 0.0;
  
  /// Graph construction
  const size_type nnz = ap(m);
  Graph G(m, nnz, ap, aj);
  G.showMe(std::cout, false);

  /// Graph reordering
  GraphTools_Metis T(G);
  timer.reset();
  T.reorder(verbose);
  t = timer.seconds();
  std::cout << "CondensedGraph:: reordering::time = " << t << std::endl;

  /// Symbolic factorization
  const auto perm = T.PermVector();
  const auto peri = T.InvPermVector();
  SymbolicTools S(m, ap, aj, perm, peri);
  S.symbolicFactorize();

  if (is_condensed) {
    timer.reset();
    ///
    /// scan weights
    ///
    size_type_1d_view_type as(do_not_initialize_tag("as"), m+1); /// original one
    size_type_1d_view_type aq(do_not_initialize_tag("aq"), m+1); /// permuted one
    as(0) = 0; aq(0) = 0;
    for (ordinal_type i=0;i<m;++i) {
      as(i+1) = as(i) + aw(i);
      aq(i+1) = aq(i) + aw(peri(i));
    }

    ///
    /// Evaporate perm and peri
    ///
    const size_type nnz_evp = as(m);
    ordinal_type_1d_view_type perm_evp(do_not_initialize_tag("perm_evp"), nnz_evp);
    ordinal_type_1d_view_type peri_evp(do_not_initialize_tag("peri_evp"), nnz_evp);

    for (ordinal_type i=0,iend=perm.extent(0);i<iend;++i) {
      const ordinal_type idx = perm(i);
      const ordinal_type ndof = aw(i);
      ///printf("i %d, idx %d, as %d, aq %d, ndof %d\n", i, idx, as(i), aq(i), ndof);
      for (ordinal_type j=0;j<ndof;++j) {
        const ordinal_type tgt = aq(idx)+j, src = as(i)+j;
        //printf("i %d, j %d, tgt %d, src %d\n", i, j, tgt, src);
        perm_evp(src) = tgt;
        peri_evp(tgt) = src;
      }
    }

    ///
    /// Evaporate supernodes
    ///
    const ordinal_type nsupernodes = S.NumSupernodes();
    const auto supernodes_evp = S.Supernodes();

    if (evaporate) {
      ordinal_type jbeg = supernodes_evp(0);
      for (ordinal_type i=0;i<nsupernodes;++i) {
        const ordinal_type 
          jend = supernodes_evp(i+1);
        ordinal_type ndof(0);
        for (ordinal_type j=jbeg;j<jend;++j) {
          const ordinal_type idx = peri(j);
          ndof += aw(idx);
        }
        supernodes_evp(i+1) = supernodes_evp(i) + ndof;
        jbeg = jend;
      }
    }

    ///
    /// Evaporate for gid colidx
    ///
    const auto gid_spanel_ptr = S.gidSuperPanelPtr();
    const auto gid_spanel_colidx = S.gidSuperPanelColIdx();

    size_type_1d_view_type gid_spanel_ptr_evp(do_not_initialize_tag("gid_spanel_ptr"), nsupernodes+1);
    if (evaporate) {
      gid_spanel_ptr_evp(0) = 0;
      for (ordinal_type i=0;i<nsupernodes;++i) {
        const ordinal_type
          jbeg = gid_spanel_ptr(i),
          jend = gid_spanel_ptr(i+1);
        ordinal_type ndof(0);
        for (ordinal_type j=jbeg;j<jend;++j) { 
          const ordinal_type idx = gid_spanel_colidx(j);
          ndof += (aq(idx+1) - aq(idx));
        }
        gid_spanel_ptr_evp(i+1) = gid_spanel_ptr_evp(i) + ndof;
      }
    }

    ordinal_type_1d_view_type gid_spanel_colidx_evp(do_not_initialize_tag("gid_spanel_colidx"), gid_spanel_ptr_evp(nsupernodes));
    if (evaporate) {
      for (ordinal_type i=0;i<nsupernodes;++i) {
        const ordinal_type
          jbeg = gid_spanel_ptr(i),
          jend = gid_spanel_ptr(i+1),
          offs = gid_spanel_ptr_evp(i);
        ordinal_type cnt(0);
        for (ordinal_type j=jbeg;j<jend;++j) { 
          const ordinal_type idx = gid_spanel_colidx(j);
          const ordinal_type
            kbeg = aq(idx),
            kend = aq(idx+1);
          for (ordinal_type k=kbeg;k<kend;++k,++cnt) 
            gid_spanel_colidx_evp(offs+cnt) = k;
        }
      }
    }

    const auto sid_spanel_ptr_evp = S.sidSuperPanelPtr();
    const auto sid_spanel_colidx_evp = S.sidSuperPanelColIdx();
    const auto blk_spanel_colidx = S.blkSuperPanelColIdx();

    ordinal_type_1d_view_type blk_spanel_colidx_evp(do_not_initialize_tag("blk_spanel_colidx"), sid_spanel_ptr_evp(nsupernodes));
    if (evaporate) {
      for (ordinal_type i=0;i<nsupernodes;++i) {
        const ordinal_type 
          jbeg = sid_spanel_ptr_evp(i), 
          jend = sid_spanel_ptr_evp(i+1)-1;
        const ordinal_type offs = gid_spanel_ptr(i);
        for (ordinal_type j=jbeg;j<jend;++j) {
          const ordinal_type 
            kbeg = blk_spanel_colidx(j),
            kend = blk_spanel_colidx(j+1);
          ordinal_type blk(0);
          for (ordinal_type k=kbeg;k<kend;++k) {
            const ordinal_type idx = gid_spanel_colidx(offs+k);
            const ordinal_type ndof = aq(idx+1) - aq(idx);
            blk += ndof;
          }
          blk_spanel_colidx_evp(j+1) = blk_spanel_colidx_evp(j) + blk;
        }
      }
    }
    
    t = timer.seconds();
    std::cout << "CondensedGraph:: evaporation::time = " << t << std::endl;

    if (verbose) {
      std::cout << "Let's evaporate\n";
      {
        std::cout << "  Permutation Vectors \n";
        for (ordinal_type i=0,iend=perm.extent(0);i<iend;++i) {    
          printf("perm %d -> %d, peri %d -> %d\n", i, perm(i), i, peri(i));
        }
        std::cout << "  Permutation Vectors Evp \n";
        for (ordinal_type i=0,iend=perm_evp.extent(0);i<iend;++i) {
          printf("perm %d -> %d, peri %d -> %d\n", i, perm_evp(i), i, peri_evp(i));
        }
      }

      {
        std::cout << "  Supernodes Evp \n";      
        for (ordinal_type i=0,iend=(nsupernodes+1);i<iend;++i) {    
          printf("i %d, supernode %d\n", i, supernodes_evp(i));
        }
      }

      {
        std::cout << "  Super Panel Sids and Blks Evp \n";
        for (ordinal_type i=0;i<nsupernodes;++i) {
          const ordinal_type 
            jbeg = sid_spanel_ptr_evp(i), 
            jend = sid_spanel_ptr_evp(i+1);
          for (ordinal_type j=jbeg;j<jend;++j) {
            printf("i %d, j %d, sid %d, blk %d\n", i,j,sid_spanel_colidx_evp(j), blk_spanel_colidx_evp(j));
          }
        }
      }
      {
        std::cout << "  Super Panel Gids Evp \n";
        for (ordinal_type i=0;i<nsupernodes;++i) {
          const ordinal_type 
            jbeg = gid_spanel_ptr_evp(i), 
            jend = gid_spanel_ptr_evp(i+1);
          for (ordinal_type j=jbeg;j<jend;++j) {
            printf("i %d, j %d, gid %d\n", i,j,gid_spanel_colidx_evp(j));
          }
        }
      }

    }
  }
  return 0;
}


int main (int argc, char *argv[]) {
  CommandLineParser opts("This example illustrates condensed graph usage");

  bool verbose = true, evaporate = true;
  std::string condensed_file = "condensed_graph.dat";
  std::string weight_file = "condensed_weight.dat";
  std::string graph_file = "graph.dat";

  opts.set_option<std::string>("condensed-file", "Input: condensed graph file", &condensed_file);
  opts.set_option<std::string>("node-weight-file", "Input: weight of nodes file", &weight_file);
  opts.set_option<std::string>("graph-file", "Input: evaporated graph", &graph_file);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("evaporate", "Flag for evaporating", &evaporate);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);

  typedef Kokkos::DefaultHostExecutionSpace host_space;
  int r_val = 0;
  {
    using size_type_1d_view_type = Kokkos::View<size_type*,host_space>;
    using ordinal_type_1d_view_type = Kokkos::View<ordinal_type*,host_space>;

    size_type_1d_view_type A_ap, Ac_ap;
    ordinal_type_1d_view_type A_aj, A_aw, Ac_aj, Ac_aw;

    readGraphFile(graph_file, A_ap, A_aj); 
    readGraphFile(condensed_file, Ac_ap, Ac_aj); 
    readWeightFile(weight_file, Ac_aw); 

    ordinal_type A_m = A_ap.extent(0) - 1, Ac_m = Ac_ap.extent(0) - 1;

    std::cout << "CondensedGraph:: A" << std::endl;
    doSymbolicFactorization(A_m, A_ap, A_aj, A_aw, verbose, false, evaporate);

    std::cout << "CondensedGraph:: Ac" << std::endl;
    doSymbolicFactorization(Ac_m, Ac_ap, Ac_aj, Ac_aw, verbose, true, evaporate);
  }
  Kokkos::finalize();

  return r_val;
}
