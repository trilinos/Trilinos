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

#include "TraitsTests.hpp"

#include "Sacado_Fad_SimpleFad.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"

typedef ::testing::Types<
  Sacado::Fad::DFad<double>,
  Sacado::Fad::SFad<double,5>,
  Sacado::Fad::SLFad<double,10>,
  Sacado::Fad::SimpleFad<double>,
  Sacado::Fad::DVFad<double>,
  Sacado::ELRFad::DFad<double>,
  Sacado::ELRFad::SFad<double,5>,
  Sacado::ELRFad::SLFad<double,10>,
  Sacado::CacheFad::DFad<double>,
  Sacado::CacheFad::SFad<double,5>,
  Sacado::CacheFad::SLFad<double,10>,
  Sacado::ELRCacheFad::DFad<double>,
  Sacado::ELRCacheFad::SFad<double,5>,
  Sacado::ELRCacheFad::SLFad<double,10>,
  Sacado::LFad::LogicalSparse<double,bool>,
  Sacado::FlopCounterPack::ScalarFlopCounter<double>,
  Sacado::Tay::Taylor<double>,
  Sacado::Tay::CacheTaylor<double>,
  Sacado::Rad::ADvar<double>,
  Sacado::Rad2::ADvar<double>,
  Sacado::RadVec::ADvar<double>
  > ADTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(Traits, TraitsTests, ADTypes);
