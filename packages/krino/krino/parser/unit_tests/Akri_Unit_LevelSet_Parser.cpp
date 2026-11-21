// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <Akri_LevelSet_Parser.hpp>
#include <Akri_Parser.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Surface_Manager.hpp>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

class SurfaceNameParserTest : public ::testing::Test
{
protected:
  std::string parse_name(const std::string & yamlText)
  {
    return parse_level_set_or_bounding_surface_name(Parser::Node(YAML::Load(yamlText)));
  }
};

TEST_F(SurfaceNameParserTest, PlainName)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  EXPECT_EQ("LS1", parse_name("name: LS1\n"));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceNameParserTest, MissingNameIsEmptyWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  EXPECT_EQ("", parse_name("distance_variable: D\n"));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceNameParserTest, PrefixIndexAndDigits)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  EXPECT_EQ("LS01", parse_name("name:\n  prefix: LS\n  index: 1\n  digits: 2\n"));
  EXPECT_EQ("LS10", parse_name("name: {prefix: LS, index: 10, digits: 2}\n"));
  EXPECT_EQ("LS007", parse_name("name: {prefix: LS, index: 7, digits: 3}\n"));
  EXPECT_EQ("LS0", parse_name("name: {prefix: LS, index: 0, digits: 1}\n"));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceNameParserTest, PrefixIndexAndMaxIndex)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  EXPECT_EQ("LS01", parse_name("name:\n  prefix: LS\n  index: 1\n  max_index: 39\n"));
  EXPECT_EQ("LS39", parse_name("name: {prefix: LS, index: 39, max_index: 39}\n"));
  EXPECT_EQ("LS000", parse_name("name: {prefix: LS, index: 0, max_index: 100}\n"));
  EXPECT_EQ("LS9", parse_name("name: {prefix: LS, index: 9, max_index: 9}\n"));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceNameParserTest, InvalidComposedNameIsEmptyAndDooms)
{
  for (const std::string yamlText : {
      "name: {index: 1, digits: 2}\n",                          // no prefix
      "name: {prefix: LS, digits: 2}\n",                        // no index
      "name: {prefix: LS, index: 1}\n",                         // neither digits nor max_index
      "name: {prefix: LS, index: 1, digits: 2, max_index: 39}\n", // both
      "name: {prefix: LS, index: -1, digits: 2}\n",             // negative index
      "name: {prefix: LS, index: 1, digits: 0}\n",              // non-positive digits
      "name: {prefix: LS, index: 100, digits: 2}\n",            // index too wide for digits
      "name: {prefix: LS, index: 40, max_index: 39}\n",         // index above max_index
      "name: {prefix: LS, index: one, digits: 2}\n"})           // non-integer index
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ("", parse_name(yamlText)) << yamlText;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

// A level set whose initial surface has a motion is initialized from it again at each step.
class LevelSetMotionParserTest : public ::testing::Test
{
protected:
  LevelSetMotionParserTest() : meta(3)
  {
    auto & coords = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(coords, meta.universal_part(), 3, nullptr);
    meta.set_coordinate_field(&coords);
    AuxMetaData::get_or_create(meta);
  }

  LevelSet & parse_level_set(const std::string & levelSetContents)
  {
    const std::string yamlText = "level_set_interfaces:\n- name: LS\n" + levelSetContents;
    LevelSet_Parser::parse(Parser::Node(YAML::Load(yamlText)), meta, timer);
    return Surface_Manager::get(meta).get_levelset("LS");
  }

  stk::mesh::MetaData meta;
  stk::diag::Timer timer{"LS Parser Timer", sierra::Diag::sierraTimer()};
};

TEST_F(LevelSetMotionParserTest, movingInitialSurface_reinitializedEveryStep)
{
  const LevelSet & ls = parse_level_set(
      "  initial_conditions:\n"
      "  - sphere:\n"
      "    center: [0,0,0]\n"
      "    radius: 0.5\n"
      "    motion:\n"
      "      translational_velocity: [1, 0, 0]\n");
  EXPECT_TRUE(ls.get_reinitialize_every_step());
}

TEST_F(LevelSetMotionParserTest, movingSurfaceInCsgTree_reinitializedEveryStep)
{
  const LevelSet & ls = parse_level_set(
      "  csg_initialization:\n"
      "    union:\n"
      "    - sphere:\n"
      "      center: [0,0,0]\n"
      "      radius: 0.5\n"
      "    - sphere:\n"
      "      center: [1,0,0]\n"
      "      radius: 0.5\n"
      "      motion:\n"
      "        translational_velocity: [1, 0, 0]\n");
  EXPECT_TRUE(ls.get_reinitialize_every_step());
}

TEST_F(LevelSetMotionParserTest, stringFunctionOfTime_reinitializedEveryStep)
{
  const LevelSet & ls = parse_level_set(
      "  initial_conditions:\n"
      "  - string_function:\n"
      "    expression: \"x - t\"\n");
  EXPECT_TRUE(ls.get_reinitialize_every_step());
}

TEST_F(LevelSetMotionParserTest, stringFunctionNotOfTime_notReinitialized)
{
  const LevelSet & ls = parse_level_set(
      "  initial_conditions:\n"
      "  - string_function:\n"
      "    expression: \"x - 0.5\"\n");
  EXPECT_FALSE(ls.get_reinitialize_every_step());
}

TEST_F(LevelSetMotionParserTest, staticInitialSurface_notReinitialized)
{
  const LevelSet & ls = parse_level_set(
      "  initial_conditions:\n"
      "  - sphere:\n"
      "    center: [0,0,0]\n"
      "    radius: 0.5\n");
  EXPECT_FALSE(ls.get_reinitialize_every_step());
}

TEST_F(LevelSetMotionParserTest, motionAndInterfaceVelocity_dooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  parse_level_set(
      "  interface_velocity: [\"1\",\"0\",\"0\"]\n"
      "  initial_conditions:\n"
      "  - sphere:\n"
      "    center: [0,0,0]\n"
      "    radius: 0.5\n"
      "    motion:\n"
      "      translational_velocity: [1, 0, 0]\n");
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

#else

TEST(SurfaceNameParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
