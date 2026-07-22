// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
