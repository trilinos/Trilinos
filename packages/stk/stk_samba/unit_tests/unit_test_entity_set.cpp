#include <gtest/gtest.h>

#include <samba/entity_block_key.hpp>

#include <iostream>
#include <sstream>

TEST(samba, entity_block_key)
{

  using namespace samba;

  entity_block_key _1={1}, _2={2};

  entity_block_key _3 = entity_block_key::invalid();

  {
    std::ostringstream out;
    out << _3;
    EXPECT_EQ(out.str(),"{entity_block:invalid}");
  }

  EXPECT_EQ( entity_block_key::create(1), _1);
  EXPECT_NE( _1, _2);
  EXPECT_LT( _1, _2);
  EXPECT_LE( _1, _2);
  EXPECT_GT( _2, _1);
  EXPECT_GE( _2, _2);

  {
    std::ostringstream out;
    out << _1;
    EXPECT_EQ(out.str(),"{entity_block:1}");
  }
}



