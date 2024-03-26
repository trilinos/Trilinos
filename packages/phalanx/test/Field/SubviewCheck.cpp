// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Sacado.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include <limits>
#include <cmath>


TEUCHOS_UNIT_TEST(field, all)
{
  printf("\n");

  using fad = Sacado::Fad::DFad<double>;
  const int num_cell = 5;
  const int num_qp = 4;
  const int num_dim = 3;
  const int num_deriv = 2;

  // Must compile with hierarchic dfad enabled to show the issue. The ContiguousLayout has incorrect values.
  PHX::View<fad***> x("xrk",num_cell,num_qp,num_dim,num_deriv+1);
  auto x_host = Kokkos::create_mirror_view(x);
  
  for (int cell=0; cell < num_cell; ++cell) {
    for (int qp=0; qp < num_qp; ++qp) {
      for (int dim = 0; dim < num_dim; ++dim) {
        x_host(cell,qp,dim).val() = 100.*cell + 10.*qp + 1.*dim;
        for (int deriv = 0; deriv < num_deriv; ++deriv) {
          x_host(cell,qp,dim).fastAccessDx(deriv) = x_host(cell,qp,dim).val() + (1.0*deriv)/10.;
        }
      }
    }
  }
  Kokkos::deep_copy(x,x_host);

  printf("\nPRINTING INITIAL VALUES\n");

  // Check the values on device
  Kokkos::TeamPolicy<PHX::exec_space> tp(num_cell,Kokkos::AUTO(),32);
  Kokkos::parallel_for("print initial values from device",tp,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, num_qp), [&](const int qp) {
      for (int dim=0; dim < num_dim; ++dim) {
        Kokkos::single(Kokkos::PerThread(team),[&] () {
          printf("x(%d,%d,%d).val()=%f,dx(0)=%f,dx(1)=%f\n",cell,qp,dim,x(cell,qp,dim).val(),x(cell,qp,dim).fastAccessDx(0),x(cell,qp,dim).fastAccessDx(1));
        });
      }
    });
  });
  Kokkos::fence();

  printf("\nPRINTING FUNCTOR VALUES\n");

  // Take a subview and assign value to itself
  Kokkos::parallel_for("subview assign",tp,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, num_qp), [&](const int qp) {
        
      // Subview of a DFad with Contiguous layout fails!
      auto x_tmp = Kokkos::subview(x,cell,qp,Kokkos::ALL());

      for (int dim=0; dim < num_dim; ++dim) {

        team.team_barrier();
        Kokkos::single(Kokkos::PerThread(team),[&] () {
          printf("x(%d,%d,%d).val()=%f,dx(0)=%f,dx(1)=%f : x_tmp(%d).val()=%f,dx(0)=%f,dx(1)=%f\n",
                 cell,qp,dim,x(cell,qp,dim).val(),x(cell,qp,dim).fastAccessDx(0),x(cell,qp,dim).fastAccessDx(1),
                 dim,x_tmp(dim).val(),x_tmp(dim).fastAccessDx(0),x_tmp(dim).fastAccessDx(1));
        });


        x(cell,qp,dim) = x_tmp(dim);

        // This works fine, not the assignment operator, but the subview
        // x(cell,qp,dim) = x(cell,qp,dim);
      }
    });
  });
  Kokkos::fence();

  printf("\nPRINTING FINAL VALUES\n");

  // Check the values on device again
  Kokkos::parallel_for("print final values from device",tp,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, num_qp), [&](const int qp) {
      for (int dim=0; dim < num_dim; ++dim) {
        Kokkos::single(Kokkos::PerThread(team),[&] () {
          printf("x(%d,%d,%d).val()=%f,dx(0)=%f,dx(1)=%f\n",cell,qp,dim,x(cell,qp,dim).val(),x(cell,qp,dim).fastAccessDx(0),x(cell,qp,dim).fastAccessDx(1));
        });
      }
    });
  });
  Kokkos::fence();
}
