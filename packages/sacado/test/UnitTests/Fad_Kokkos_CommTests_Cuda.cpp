// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Sacado.hpp"
#include "Fad_CommTests.hpp"

typedef int Ordinal;
typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,10> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,5> Fad_SFadType;

typedef Sacado::CacheFad::DFad<double> CacheFad_DFadType;
typedef Sacado::CacheFad::SLFad<double,10> CacheFad_SLFadType;
typedef Sacado::CacheFad::SFad<double,5> CacheFad_SFadType;

typedef Sacado::ELRFad::DFad<double> ELRFad_DFadType;
typedef Sacado::ELRFad::SLFad<double,10> ELRFad_SLFadType;
typedef Sacado::ELRFad::SFad<double,5> ELRFad_SFadType;

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheFad_DFadType;
typedef Sacado::ELRCacheFad::SLFad<double,10> ELRCacheFad_SLFadType;
typedef Sacado::ELRCacheFad::SFad<double,5> ELRCacheFad_SFadType;
Sacado::Random<double> rnd;

FAD_KOKKOS_COMM_TESTS_CUDA(Fad_DFadType, Fad_DFad)
FAD_KOKKOS_COMM_TESTS_CUDA(Fad_SLFadType, Fad_SLFad)
FAD_KOKKOS_COMM_TESTS_CUDA(Fad_SFadType, Fad_SFad)

FAD_KOKKOS_COMM_TESTS_CUDA(CacheFad_DFadType, CacheFad_DFad)
FAD_KOKKOS_COMM_TESTS_CUDA(CacheFad_SLFadType, CacheFad_SLFad)
FAD_KOKKOS_COMM_TESTS_CUDA(CacheFad_SFadType, CacheFad_SFad)

FAD_KOKKOS_COMM_TESTS_CUDA(ELRFad_DFadType, ELRFad_DFad)
FAD_KOKKOS_COMM_TESTS_CUDA(ELRFad_SLFadType, ELRFad_SLFad)
FAD_KOKKOS_COMM_TESTS_CUDA(ELRFad_SFadType, ELRFad_SFad)

FAD_KOKKOS_COMM_TESTS_CUDA(ELRCacheFad_DFadType, ELRCacheFad_DFad)
FAD_KOKKOS_COMM_TESTS_CUDA(ELRCacheFad_SLFadType, ELRCacheFad_SLFad)
FAD_KOKKOS_COMM_TESTS_CUDA(ELRCacheFad_SFadType, ELRCacheFad_SFad)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize cuda
  Kokkos::InitArguments init_args;
  init_args.device_id = 0;
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize cuda
  Kokkos::finalize();

  return ret;
}
