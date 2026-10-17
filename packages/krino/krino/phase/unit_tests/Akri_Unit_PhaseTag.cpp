// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <vector>

namespace krino {

namespace {
constexpr int INSIDE = -1;
constexpr int OUTSIDE = +1;
constexpr int INTERFACE = 0;

Surface_Identifier ls(const unsigned id) { return Surface_Identifier(id); }
}

//--------------------------------------------------------------------
// LS_SideTag
//--------------------------------------------------------------------

TEST(LS_SideTag, Accessors)
{
  const LS_SideTag side(ls(3), INSIDE);
  EXPECT_EQ(3u, side.get_ls_identifier().get());
  EXPECT_EQ(INSIDE, side.get_ls_sign());
  EXPECT_FALSE(side.is_interface());
}

TEST(LS_SideTag, IsInterface_WhenSignIsZero)
{
  EXPECT_TRUE(LS_SideTag(ls(0), INTERFACE).is_interface());
  EXPECT_FALSE(LS_SideTag(ls(0), INSIDE).is_interface());
  EXPECT_FALSE(LS_SideTag(ls(0), OUTSIDE).is_interface());
}

TEST(LS_SideTag, SignChangingConstructor_KeepsIdentifier)
{
  const LS_SideTag original(ls(7), INSIDE);
  const LS_SideTag flipped(original, OUTSIDE);
  EXPECT_EQ(7u, flipped.get_ls_identifier().get());
  EXPECT_EQ(OUTSIDE, flipped.get_ls_sign());
}

TEST(LS_SideTag, OppositeSide_NegatesSignKeepsIdentifier)
{
  const LS_SideTag inside(ls(2), INSIDE);
  const LS_SideTag opp = inside.opposite_side();
  EXPECT_EQ(2u, opp.get_ls_identifier().get());
  EXPECT_EQ(OUTSIDE, opp.get_ls_sign());

  // Opposite of an interface side is still the interface.
  EXPECT_EQ(INTERFACE, LS_SideTag(ls(2), INTERFACE).opposite_side().get_ls_sign());
}

TEST(LS_SideTag, EqualityAndInequality)
{
  const LS_SideTag a(ls(1), INSIDE);
  EXPECT_TRUE(a == LS_SideTag(ls(1), INSIDE));
  EXPECT_FALSE(a != LS_SideTag(ls(1), INSIDE));

  EXPECT_TRUE(a != LS_SideTag(ls(1), OUTSIDE)); // different sign
  EXPECT_TRUE(a != LS_SideTag(ls(2), INSIDE));  // different identifier
}

TEST(LS_SideTag, LessThan_OrdersByIdentifierThenSign)
{
  EXPECT_TRUE(LS_SideTag(ls(1), OUTSIDE) < LS_SideTag(ls(2), INSIDE)); // identifier dominates
  EXPECT_TRUE(LS_SideTag(ls(1), INSIDE) < LS_SideTag(ls(1), OUTSIDE)); // tie-break on sign
  EXPECT_FALSE(LS_SideTag(ls(1), INSIDE) < LS_SideTag(ls(1), INSIDE)); // irreflexive
}

TEST(LS_SideTag, Composite_UndeclaredMapsToItself)
{
  // Identifier chosen to be one that no other test declares as a composite constituent.
  const LS_SideTag side(ls(9000), INSIDE);
  const LS_SideTag comp = side.composite();
  EXPECT_EQ(9000u, comp.get_ls_identifier().get());
  EXPECT_EQ(INSIDE, comp.get_ls_sign());
}

TEST(LS_SideTag, DeclareComposite_RemapsConstituentIdentifierPreservingSign)
{
  // the_composite_ls_map is a static global; use identifiers unique to this test.
  LS_SideTag::declare_composite(ls(9101), ls(9100));

  const LS_SideTag comp = LS_SideTag(ls(9101), INSIDE).composite();
  EXPECT_EQ(9100u, comp.get_ls_identifier().get());
  EXPECT_EQ(INSIDE, comp.get_ls_sign()); // sign is preserved through compositing

  // Re-declaring updates the target.
  LS_SideTag::declare_composite(ls(9101), ls(9102));
  EXPECT_EQ(9102u, LS_SideTag(ls(9101), OUTSIDE).composite().get_ls_identifier().get());
}

//--------------------------------------------------------------------
// PhaseTag
//--------------------------------------------------------------------

TEST(PhaseTag, DefaultConstructed_IsEmpty)
{
  PhaseTag phase;
  EXPECT_TRUE(phase.empty());
  EXPECT_TRUE(phase.ls_sides().empty());
}

TEST(PhaseTag, AddAndContain)
{
  PhaseTag phase;
  phase.add(ls(1), INSIDE);

  EXPECT_FALSE(phase.empty());
  EXPECT_TRUE(phase.contain(ls(1), INSIDE));
  EXPECT_FALSE(phase.contain(ls(1), OUTSIDE));
  EXPECT_FALSE(phase.contain(ls(2), INSIDE));
}

TEST(PhaseTag, Clear)
{
  PhaseTag phase;
  phase.add(ls(1), INSIDE);
  phase.clear();
  EXPECT_TRUE(phase.empty());
}

TEST(PhaseTag, AddIsIdempotent)
{
  PhaseTag phase;
  phase.add(ls(1), INSIDE);
  phase.add(ls(1), INSIDE);
  EXPECT_EQ(1u, phase.ls_sides().size());
}

TEST(PhaseTag, AddBothSidesOfNonCompositeLS_KeepsBoth)
{
  // For an ordinary (non-composite) level set, the two sides are independent
  // tags and both are retained; the priority rule does not apply here.
  PhaseTag phase;
  phase.add(ls(1), OUTSIDE);
  phase.add(ls(1), INSIDE);

  EXPECT_EQ(2u, phase.ls_sides().size());
  EXPECT_TRUE(phase.contain(ls(1), INSIDE));
  EXPECT_TRUE(phase.contain(ls(1), OUTSIDE));
}

TEST(PhaseTag, AddPriorityRule_InsideBeatsOutside)
{
  // The priority rule (-1 inside beats +1 outside) only fires when constituent
  // level sets are remapped onto a shared composite. Adding the outside side,
  // then the inside side of the same composite: inside wins.
  // the_composite_ls_map is a static global; use identifiers unique to this test.
  LS_SideTag::declare_composite(ls(9201), ls(9200));
  LS_SideTag::declare_composite(ls(9202), ls(9200));

  PhaseTag phase;
  phase.add(ls(9201), OUTSIDE);
  phase.add(ls(9202), INSIDE);

  EXPECT_EQ(1u, phase.ls_sides().size());
  EXPECT_TRUE(phase.contain(ls(9200), INSIDE));
  EXPECT_FALSE(phase.contain(ls(9200), OUTSIDE));
}

TEST(PhaseTag, AddPriorityRule_OutsideDoesNotDisplaceInside)
{
  // Reverse order: the composite's inside is already present, adding the
  // outside side of another constituent must not override it.
  // the_composite_ls_map is a static global; use identifiers unique to this test.
  LS_SideTag::declare_composite(ls(9211), ls(9210));
  LS_SideTag::declare_composite(ls(9212), ls(9210));

  PhaseTag phase;
  phase.add(ls(9211), INSIDE);
  phase.add(ls(9212), OUTSIDE);

  EXPECT_EQ(1u, phase.ls_sides().size());
  EXPECT_TRUE(phase.contain(ls(9210), INSIDE));
  EXPECT_FALSE(phase.contain(ls(9210), OUTSIDE));
}

TEST(PhaseTag, AddPhase_MergesSides)
{
  PhaseTag a;
  a.add(ls(1), INSIDE);

  PhaseTag b;
  b.add(ls(2), OUTSIDE);
  b.add(ls(3), INSIDE);

  a.add(b);
  EXPECT_EQ(3u, a.ls_sides().size());
  EXPECT_TRUE(a.contain(ls(1), INSIDE));
  EXPECT_TRUE(a.contain(ls(2), OUTSIDE));
  EXPECT_TRUE(a.contain(ls(3), INSIDE));
}

TEST(PhaseTag, AddOppositeSides_FlipsEverySign)
{
  PhaseTag source;
  source.add(ls(1), INSIDE);
  source.add(ls(2), OUTSIDE);

  PhaseTag phase;
  phase.add_opposite_sides(source);

  EXPECT_TRUE(phase.contain(ls(1), OUTSIDE));
  EXPECT_TRUE(phase.contain(ls(2), INSIDE));
}

TEST(PhaseTag, Remove)
{
  PhaseTag phase;
  phase.add(ls(1), INSIDE);
  phase.add(ls(2), OUTSIDE);

  phase.remove(ls(1), INSIDE);
  EXPECT_FALSE(phase.contain(ls(1), INSIDE));
  EXPECT_TRUE(phase.contain(ls(2), OUTSIDE));
}

TEST(PhaseTag, ContainPhase_IsSubsetTest)
{
  PhaseTag superset;
  superset.add(ls(1), INSIDE);
  superset.add(ls(2), INSIDE);

  PhaseTag subset;
  subset.add(ls(1), INSIDE);

  EXPECT_TRUE(superset.contain(subset));
  EXPECT_FALSE(subset.contain(superset));

  // The empty phase is contained by everything.
  EXPECT_TRUE(superset.contain(PhaseTag()));
}

TEST(PhaseTag, IsNonconformal_WhenAnySideIsInterface)
{
  PhaseTag conformal;
  conformal.add(ls(1), INSIDE);
  conformal.add(ls(2), OUTSIDE);
  EXPECT_FALSE(conformal.is_nonconformal());

  PhaseTag nonconformal;
  nonconformal.add(ls(1), INSIDE);
  nonconformal.add(ls(2), INTERFACE);
  EXPECT_TRUE(nonconformal.is_nonconformal());
}

TEST(PhaseTag, Equality_IgnoresInsertionOrder)
{
  PhaseTag a;
  a.add(ls(1), INSIDE);
  a.add(ls(2), OUTSIDE);

  PhaseTag b;
  b.add(ls(2), OUTSIDE);
  b.add(ls(1), INSIDE);

  EXPECT_TRUE(a == b);
  EXPECT_FALSE(a != b);
}

TEST(PhaseTag, Inequality_DifferentSizeOrContents)
{
  PhaseTag a;
  a.add(ls(1), INSIDE);

  PhaseTag bigger;
  bigger.add(ls(1), INSIDE);
  bigger.add(ls(2), INSIDE);
  EXPECT_TRUE(a != bigger); // differing size

  PhaseTag different;
  different.add(ls(1), OUTSIDE);
  EXPECT_TRUE(a != different); // same size, different side
}

TEST(PhaseTag, Intersection_KeepsOnlySharedSides)
{
  PhaseTag a;
  a.add(ls(1), INSIDE);
  a.add(ls(2), INSIDE);
  a.add(ls(3), OUTSIDE);

  PhaseTag b;
  b.add(ls(2), INSIDE);
  b.add(ls(3), INSIDE); // different sign than a's (3, OUTSIDE)
  b.add(ls(4), INSIDE);

  a &= b;

  EXPECT_EQ(1u, a.ls_sides().size());
  EXPECT_TRUE(a.contain(ls(2), INSIDE));
  EXPECT_FALSE(a.contain(ls(1), INSIDE));  // not in b
  EXPECT_FALSE(a.contain(ls(3), OUTSIDE)); // sign mismatch with b
}

TEST(PhaseTag, Intersection_WithEmptyYieldsEmpty)
{
  PhaseTag a;
  a.add(ls(1), INSIDE);
  a &= PhaseTag();
  EXPECT_TRUE(a.empty());
}

}
