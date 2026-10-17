// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_PhaseInput.hpp>

namespace krino {

TEST(SurfaceNames, NewNames_GetIdentifiersInOrder)
{
  SurfaceNames names;
  EXPECT_EQ(Surface_Identifier(0), names.get_identifier("ls1"));
  EXPECT_EQ(Surface_Identifier(1), names.get_identifier("ls2"));
  EXPECT_EQ((std::vector<std::string>{"LS1", "LS2"}), names.get_names());
}

TEST(SurfaceNames, RepeatedName_IsCaseInsensitiveAndKeepsItsIdentifier)
{
  SurfaceNames names;
  const Surface_Identifier ls1 = names.get_identifier("ls1");
  names.get_identifier("ls2");
  EXPECT_EQ(ls1, names.get_identifier("LS1"));
  EXPECT_EQ(2u, names.get_names().size());
}

} // namespace krino
