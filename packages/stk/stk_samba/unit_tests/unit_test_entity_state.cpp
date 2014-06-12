#include <gtest/gtest.h>

#include <samba/entity_state.hpp>

#include <iostream>
#include <sstream>


TEST(samba, entity_state)
{
  using namespace samba;

  {
    entity_state r = entity_state::universe();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:universe}");
  }

  {
    entity_state r = entity_state::modified();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:modified}");
  }

  {
    entity_state r = entity_state::owned();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:owned}");
  }

  {
    entity_state r = entity_state::shared();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:shared}");
  }

  {
    entity_state r = entity_state::ghosted();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:ghosted}");
  }

  {
    entity_state r = entity_state::invalid();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_state:invalid}");
  }

}



