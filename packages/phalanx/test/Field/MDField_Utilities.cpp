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


#include <vector>

#include "Phalanx_config.hpp"
#include "Phalanx.hpp"
#include "Phalanx_MDField_Utilities.hpp"
#include "Phalanx_DimTag.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

struct Dim : public PHX::DimTag {
  Dim () {}
  const char* name () const { static const char name[] = "Dim"; return name; }
  static const Dim& tag() { static const Dim tag; return tag; }
};

namespace {
PHX::DataLayout* makeLayout (const int rank, const int* d) {
  switch (rank) {
  case 1: return new PHX::MDALayout<Dim>(
    d[0]);
  case 2: return new PHX::MDALayout<Dim,Dim>(
    d[0], d[1]);
  case 3: return new PHX::MDALayout<Dim,Dim,Dim>(
    d[0], d[1], d[2]);
  case 4: return new PHX::MDALayout<Dim,Dim,Dim,Dim>(
    d[0], d[1], d[2], d[3]);
  case 5: return new PHX::MDALayout<Dim,Dim,Dim,Dim,Dim>(
    d[0], d[1], d[2], d[3], d[4]);
  case 6: return new PHX::MDALayout<Dim,Dim,Dim,Dim,Dim,Dim>(
    d[0], d[1], d[2], d[3], d[4], d[5]);
  case 7: return new PHX::MDALayout<Dim,Dim,Dim,Dim,Dim,Dim,Dim>(
    d[0], d[1], d[2], d[3], d[4], d[5], d[6]);
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "bad rank");
    return 0;
  }
}
} // namespace

void testRank (const int rank) {
  // A random set of dimensions.
  static const int dims[] = {3, 2, 4, 2, 2, 1, 3};

  Teuchos::RCP<PHX::DataLayout> layout = Teuchos::rcp(makeLayout(rank, dims));
  PHX::MDField<PHX::MyTraits::FadType> f("test", layout);
  PHX::MDField<double> d("test", layout);
  const std::vector<PHX::index_size_type> ddims(1, 8);
  d.setFieldData(
    PHX::KokkosViewFactory<double, PHX::Device>::buildView(
      d.fieldTag(), ddims));
  f.setFieldData(
    PHX::KokkosViewFactory<PHX::MyTraits::FadType, PHX::Device>::buildView(
      f.fieldTag(), ddims));
  
  const PHX::MyTraits::FadType
    a(ddims[0], 0, 1),
    b(ddims[0], 0, 2),
    c = a*b;

  PHX::MDFieldIterator<PHX::MyTraits::FadType> fi(f);
  PHX::MDFieldIterator<double> di(d);
  // Test (1) every method in MDFieldIterator and (2) basic FadType read/write.
  unsigned int k;
  for (k = 0; ! fi.done(); ++k) {
    if (k % 2 == 0) *di = static_cast<double>(k);
    else di.ref() = static_cast<double>(k);
    *fi = a;
    *fi *= b;
    TEUCHOS_ASSERT(fi->val() == c.val());
    TEUCHOS_ASSERT(fi->dx(0) == c.dx(0));
    TEUCHOS_ASSERT(*fi == c);
    TEUCHOS_ASSERT(fi.idx() == k);
    ++di; fi++;
  }
  TEUCHOS_ASSERT(k == d.size());
}

TEUCHOS_UNIT_TEST(mdfield, Utilities) {
  Teuchos::RCP<Teuchos::Time>
    total_time = Teuchos::TimeMonitor::getNewTimer("Total Run Time");
  Teuchos::TimeMonitor tm(*total_time);
  
  PHX::InitializeKokkosDevice();

  for (int rank = 1; rank <= 7; ++rank) testRank(rank);

  PHX::FinalizeKokkosDevice();  
  Teuchos::TimeMonitor::summarize();
}
