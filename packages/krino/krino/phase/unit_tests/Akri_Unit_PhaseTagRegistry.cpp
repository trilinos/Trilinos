// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_PhaseTag.hpp>
#include <Akri_PhaseTagRegistry.hpp>
#include <Akri_Surface_Identifier.hpp>

namespace krino {

namespace {
PhaseTag make_phase(const std::vector<std::pair<unsigned, int>> & sides)
{
  PhaseTag phase;
  for (auto && [ls, sign] : sides)
    phase.add(Surface_Identifier(ls), sign);
  return phase;
}
}

TEST(PhaseTagRegistry, EmptyRegistryHasNoPhases)
{
  const PhaseTagRegistry registry;
  EXPECT_EQ(0u, registry.size());
}

TEST(PhaseTagRegistry, GetReturnsEqualPhase)
{
  PhaseTagRegistry registry;
  const PhaseTag phase = make_phase({{0, -1}, {1, +1}});
  EXPECT_EQ(phase, registry.get(phase));
}

TEST(PhaseTagRegistry, EqualPhasesShareOneEntry)
{
  PhaseTagRegistry registry;
  const PhaseTag & first = registry.get(make_phase({{0, -1}, {1, +1}}));
  const PhaseTag & second = registry.get(make_phase({{1, +1}, {0, -1}}));
  EXPECT_EQ(&first, &second);
  EXPECT_EQ(1u, registry.size());
}

TEST(PhaseTagRegistry, DistinctPhasesGetDistinctEntries)
{
  PhaseTagRegistry registry;
  const PhaseTag & inside = registry.get(make_phase({{0, -1}}));
  const PhaseTag & outside = registry.get(make_phase({{0, +1}}));
  EXPECT_NE(&inside, &outside);
  EXPECT_EQ(2u, registry.size());
}

TEST(PhaseTagRegistry, RegisteredPhaseAddressesStayValidAsPhasesAreAdded)
{
  PhaseTagRegistry registry;
  const PhaseTag & first = registry.get(make_phase({{0, -1}}));
  for (unsigned ls=1; ls<100; ++ls)
    registry.get(make_phase({{ls, -1}}));
  EXPECT_EQ(&first, &registry.get(make_phase({{0, -1}})));
  EXPECT_EQ(make_phase({{0, -1}}), first);
}

TEST(PhaseTagRegistry, EmptyPhaseIsEmptyAndShared)
{
  EXPECT_TRUE(PhaseTagRegistry::empty_phase().empty());
  EXPECT_EQ(&PhaseTagRegistry::empty_phase(), &PhaseTagRegistry::empty_phase());
}

} // namespace krino
