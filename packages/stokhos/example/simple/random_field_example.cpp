// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// random_field_example
//
//  usage:
//     random_field_example
//
//  output:
//     prints out KL expansion for an exponential random field

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_KL_ExponentialRandomField.hpp"

// Functor for evaluating random field on the device
template <typename Value, typename Device = Kokkos::DefaultExecutionSpace>
struct RF {
  typedef Value value_type;
  typedef Device execution_space;

  typedef Stokhos::KL::ExponentialRandomField<value_type,execution_space> rf_type;
  typedef Kokkos::View<value_type*,execution_space> point_type;
  typedef Kokkos::View<value_type*,execution_space> rv_type;
  typedef Kokkos::View<value_type*,execution_space> result_type;
  const rf_type rf;
  const point_type x;
  const rv_type rv;
  const result_type y;

  RF(const rf_type& rf_, const point_type& x_, const rv_type& rv_) :
    rf(rf_), x(x_), rv(rv_), y("y",1)
  {
    Kokkos::parallel_for(1, *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const unsigned i) const {
    y(0) = rf.evaluate(x, rv);
  }
};

int main(int argc, char **argv)
{
  Kokkos::initialize();

  try {

    // Create random field
    int M = 10;
    Teuchos::ParameterList solverParams;
    solverParams.set("Number of KL Terms", M);
    solverParams.set("Mean", 1.0);
    solverParams.set("Standard Deviation", 0.1);
    int ndim = 3;
    Teuchos::Array<double> domain_upper(ndim), domain_lower(ndim),
      correlation_length(ndim);
    for (int i=0; i<ndim; i++) {
      domain_upper[i] = 1.0;
      domain_lower[i] = 0.0;
      correlation_length[i] = 10.0;
    }
    solverParams.set("Domain Upper Bounds", domain_upper);
    solverParams.set("Domain Lower Bounds", domain_lower);
    solverParams.set("Correlation Lengths", correlation_length);
    Stokhos::KL::ExponentialRandomField<double> rf(solverParams);
    rf.print(std::cout);

    // Evaluate random field at a point
    Teuchos::Array<double> x(ndim);
    for (int i=0; i<ndim; i++)
      x[i] = (domain_upper[i] + domain_lower[i])/2.0 +
        0.1*(domain_upper[i] - domain_lower[i])/2.0;
    Teuchos::Array<double> rvar(M);
    for (int i=0; i<M; i++)
      rvar[i] = 1.5;
    double result = rf.evaluate(x, rvar);
    std::cout << "result (host)  = " << result << std::endl;

    // Evaluate random field in a functor on device
    typedef Kokkos::View<double*> view_type;
    typedef view_type::HostMirror host_view_type;
    view_type x_view("x", ndim);
    host_view_type host_x = Kokkos::create_mirror_view(x_view);
    for (int i=0; i<ndim; i++)
      host_x(i) = x[i];
    Kokkos::deep_copy(x_view, host_x);
    view_type rvar_view("rvar", M);
    host_view_type host_rvar = Kokkos::create_mirror_view(rvar_view);
    for (int i=0; i<M; i++)
      host_rvar(i) = rvar[i];
    Kokkos::deep_copy(rvar_view, host_rvar);
    RF<double> rf_func(rf, x_view, rvar_view);
    host_view_type host_y = Kokkos::create_mirror_view(rf_func.y);
    Kokkos::deep_copy(host_y, rf_func.y);
    double result2 = host_y(0);
    std::cout << "result (device)= " << result2 << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  Kokkos::finalize();
}
