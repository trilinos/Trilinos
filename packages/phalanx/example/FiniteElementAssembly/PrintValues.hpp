// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_EXAMPLE_PRINT_UTILITIES_HPP
#define PHALANX_EXAMPLE_PRINT_UTILITIES_HPP

#include <string>
#include <ostream>
#include <fstream>

namespace phx_example {

  template<typename ...Properties>
  void printResidual(const Kokkos::View<double*,Properties...>& f,
		     const std::string& description, 
		     const bool& print_to_file,
		     const std::string& filename = "")
  {
    auto host_f = Kokkos::create_mirror_view(f);
    Kokkos::deep_copy(host_f,f);
    typename PHX::exec_space().fence();
    
    std::ostream* os = &std::cout;
    std::ofstream ofs;
    if (print_to_file) {
      ofs.open(filename.c_str());
      os = &ofs;
    }

    std::ios::fmtflags os_flags( os->flags() );
    os->precision(10);

    if (description != "")
      *os << description << std::endl;
    for (int i=0; i < static_cast<int>(host_f.extent(0)); ++i)
      *os << "f(" << i << ") = " << std::fixed << host_f(i) << std::endl;

    os->flags(os_flags);

    if (print_to_file)
      ofs.close();
  }
  template<typename ...Properties, typename CrsMatrix>
  void printResidualAndJacobian(const Kokkos::View<double*,Properties...>& f,
				const CrsMatrix& J,
				const std::string& description, 
				const bool& print_to_file,
				const std::string& filename = "")
  {
    auto host_f = Kokkos::create_mirror_view(f);
    auto host_J_vals = Kokkos::create_mirror_view(J.values);
    auto host_graph = Kokkos::create_mirror(J.graph); // deep_copies automagically
    Kokkos::deep_copy(host_f,f);
    Kokkos::deep_copy(host_J_vals,J.values);
    typename PHX::exec_space().fence();

    std::ostream* os = &std::cout;
    std::ofstream ofs;
    if (print_to_file) {
      ofs.open(filename.c_str());
      os = &ofs;
    }

    if (description != "")
      *os << description << std::endl;
    for (int i=0; i < static_cast<int>(host_f.extent(0)); ++i)
      *os << "f(" << i << ") = " << host_f(i) << std::endl;

    size_t val_index = 0;
    for (typename decltype(host_graph)::size_type row=0; row < host_graph.numRows(); ++row) {
      for (int j=0; j < host_graph.rowConst(row).length; ++j) {
	*os << "J(" << row << "," << host_graph.rowConst(row).colidx(j) << ") = "
	    << host_J_vals(val_index) << std::endl;
	++val_index;
      }
    }
    
    if (print_to_file)
      ofs.close();
  }

}
#endif
