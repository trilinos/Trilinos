#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Internal.hpp"
#include "Tacho_CommandLineParser.hpp"

using namespace Tacho;

template<typename SizeType1DViewType,
	 typename OrdinalType1DViewType>
int readGraphFile(const std::string &filename,
		  SizeType1DViewType &ap,
		  OrdinalType1DViewType &aj) {
  std::ifstream in;
  in.open(filename);
  if (in.good()) {
    ordinal_type num_nodes(0);
    in >> num_nodes;
    ap = SizeType1DViewType("rowptr", num_nodes+1);
    for (ordinal_type i=0;i<ordinal_type(num_nodes+1);++i)
      in >> ap(i);

    aj = OrdinalType1DViewType("colidx", ap(num_nodes));
    for (ordinal_type i=0;i<ordinal_type(num_nodes);++i)
      for (ordinal_type j=ap(i);j<ap(i+1);++j)
	in >> aj(j);
  } else {
    std::cout <<" Failed to open the file: " << filename << std::endl;
    return -1;
  }
  return 0;
}

template<typename OrdinalType1DViewType>
int readWeightFile(const std::string &filename,
		   OrdinalType1DViewType &aw) {
  std::ifstream in;
  in.open(filename);
  if (in.good()) {
    ordinal_type total_dofs(0), num_nodes(0);
    in >> total_dofs;
    in >> num_nodes;
    aw = OrdinalType1DViewType("weight", num_nodes);
    for (ordinal_type i=0;i<num_nodes;++i)
      in >> aw(i);
  } else {
    std::cout <<" Failed to open the file: " << filename << std::endl;
    return -1;
  }
  return 0;
}

int main (int argc, char *argv[]) {
  CommandLineParser opts("This example illustrates condensed graph usage");

  bool verbose = true;
  std::string condensed_file = "condensed_graph.dat";
  std::string weight_file = "condensed_weight.dat";
  std::string graph_file = "graph.dat";

  opts.set_option<std::string>("condensed-file", "Input: condensed graph file", &condensed_file);
  opts.set_option<std::string>("node-weight-file", "Input: weight of nodes file", &weight_file);
  opts.set_option<std::string>("graph-file", "Input: evaporated graph", &graph_file);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);

  typedef Kokkos::DefaultHostExecutionSpace host_space;
  int r_val = 0;
  {
    using size_type_1d_view_type = Kokkos::View<size_type*,host_space>;
    using ordinal_type_1d_view_type = Kokkos::View<ordinal_type*,host_space>;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    size_type_1d_view_type A_ap, Ac_ap;
    ordinal_type_1d_view_type A_aj, Ac_aj, Ac_aw;

    readGraphFile(graph_file, A_ap, A_aj);
    readGraphFile(condensed_file, Ac_ap, Ac_aj);
    readWeightFile(weight_file, Ac_aw);

    ordinal_type A_m = A_ap.extent(0) - 1, Ac_m = Ac_ap.extent(0) - 1;
    size_type A_nnz = A_ap(A_m), Ac_nnz = Ac_ap(A_m);
    Graph G (A_m, A_nnz, A_ap, A_aj), Gc(Ac_m, Ac_nnz, Ac_ap, Ac_aj);

    std::cout <<" G " << std::endl;
    G.showMe(std::cout, false);

    std::cout <<" Gc " << std::endl;
    Gc.showMe(std::cout, false);
    
    GraphTools_Metis T(G), Tc(Gc);

    timer.reset();
    Tc.reorder(verbose);
    t = timer.seconds();
    std::cout << "CondensedGraph:: Ac, reordering::time = " << t << std::endl;

    timer.reset();
    T.reorder(verbose);
    t = timer.seconds();
    std::cout << "CondensedGraph:: A, reordering::time = " << t << std::endl;

  }
  Kokkos::finalize();

  return r_val;
}
