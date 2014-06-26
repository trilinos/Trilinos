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
#include "Sacado_Fad_SimpleFad.hpp"
#include "Fad_CommTests.hpp"

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;
template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage< Sacado::Fad::DMFad<double> >::defaultPool_ = NULL;

typedef int Ordinal;
typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,10> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,5> Fad_SFadType;
typedef Sacado::Fad::DMFad<double> Fad_DMFadType;
typedef Sacado::Fad::SimpleFad<double> Fad_SimpleFadType;
typedef Sacado::LFad::LogicalSparse<double,bool> Fad_LSType;
Sacado::Random<double> rnd;
FAD_COMM_TESTS(Fad_DFadType, Fad_DFad)
FAD_COMM_TESTS(Fad_SLFadType, Fad_SLFad)
FAD_COMM_TESTS(Fad_SFadType, Fad_SFad)
FAD_COMM_TESTS(Fad_DMFadType, Fad_DMFad)
FAD_COMM_TESTS(Fad_SimpleFadType, Fad_SimpleFad)
//FAD_COMM_TESTS(Fad_LSType, Fad_LogicalSparse)

// DVFad, LFad, Flop

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Sacado::Fad::MemPoolManager<double> poolManager(100);
  Sacado::Fad::MemPool *pool = poolManager.getMemoryPool(5);
  Sacado::Fad::DMFad<double>::setDefaultPool(pool);

  Sacado::Fad::MemPoolManager< Sacado::Fad::DMFad<double> > poolManager2(100);
  Sacado::Fad::MemPool *pool2 = poolManager2.getMemoryPool(5);
  Sacado::Fad::DMFad< Sacado::Fad::DMFad<double> >::setDefaultPool(pool2);

  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
