// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_CommandLineParser.hpp"
#include "Tacho_Internal.hpp"

using namespace Tacho;

int main(int argc, char *argv[]) {
  CommandLineParser opts("This example program combines data file into a single matrix market file");

  std::string graph_data_file = "graph.dat";
  std::string value_data_file = "value.dat";
  std::string matrix_market_file = "mm-out.mtx";

  opts.set_option<std::string>("graph-file", "Input graph data file ", &graph_data_file);
  opts.set_option<std::string>("value-file", "Input value data file ", &value_data_file);
  opts.set_option<std::string>("matrix-market-file", "Output matrixmarket file", &matrix_market_file);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);

  typedef Kokkos::DefaultHostExecutionSpace host_space;
  {
    typedef double value_type;
    typedef CrsMatrixBase<value_type, host_space> CrsMatrixBaseTypeHost;

    CrsMatrixBaseTypeHost A;
    using ordinal_type_array = typename CrsMatrixBaseTypeHost::ordinal_type_array;
    using size_type_array = typename CrsMatrixBaseTypeHost::size_type_array;
    using value_type_array = typename CrsMatrixBaseTypeHost::value_type_array;

    ordinal_type m(0), nnz(0);
    size_type_array ap;
    ordinal_type_array aj;
    value_type_array ax;
    {
      std::ifstream in;
      in.open(graph_data_file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << graph_data_file << std::endl;
        return -1;
      }

      in >> m;
      ap = size_type_array("ap", m + 1);
      for (ordinal_type i = 0; i < (m + 1); ++i)
        in >> ap(i);

      nnz = ap(m);
      aj = ordinal_type_array("aj", nnz);
      for (ordinal_type k = 0; k < nnz; ++k)
        in >> aj(k);
    }
    {
      std::ifstream in;
      in.open(value_data_file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << value_data_file << std::endl;
        return -1;
      }

      ax = value_type_array("ax", nnz);
      for (ordinal_type k = 0; k < nnz; ++k)
        in >> ax(k);
    }

    A.setExternalMatrix(m, m, nnz, ap, aj, ax);
    {
      std::ofstream out(matrix_market_file);
      MatrixMarket<value_type>::write(out, A);
    }
  }
  Kokkos::finalize();

  return 0;
}
