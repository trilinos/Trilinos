// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_Core.hpp"

using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

class MyFunctor {
  Kokkos::View<double***,mem_t> a_;
  Kokkos::View<double***,mem_t> b_;
  Kokkos::View<double***,mem_t> c_;
  
public:  
  MyFunctor(Kokkos::View<double***,mem_t> a,
	    Kokkos::View<double***,mem_t> b,
	    Kokkos::View<double***,mem_t> c)
    : a_(a),b_(b),c_(c) {}

  void evaluate()
  {
    Kokkos::parallel_for(Kokkos::TeamPolicy<exec_t>(a_.extent(0),Kokkos::AUTO()),*this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<exec_t>::member_type& team) const
  {
    const int cell = team.league_rank();
    const int num_pts = a_.extent(1);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_pts), [=] (const int& pt) {
      const int num_eq = a_.extent(2);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,num_eq), [=] (const int& eq) {
  	c_(cell,pt,eq) = a_(cell,pt,eq) + b_(cell,pt,eq);
      });
    });
  }
};

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc,argv);
  {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  Kokkos::View<double***,mem_t> a("a",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> b("b",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> c("c",num_cells,num_pts,num_equations);

  MyFunctor f(a,b,c);
  f.evaluate();

  }
  Kokkos::finalize();
  return 0;
}
