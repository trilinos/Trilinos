#include <gtest/gtest.h>

#include <samba/entity_rank.hpp>

#include <iostream>
#include <sstream>


TEST(samba, entity_rank)
{
  using namespace samba;

  {
    entity_rank r = entity_rank::node();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_rank:node}");
  }

  {
    entity_rank r = entity_rank::edge();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_rank:edge}");
  }

  {
    entity_rank r = entity_rank::face();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_rank:face}");
  }

  {
    entity_rank r = entity_rank::element();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_rank:element}");
  }

  {
    entity_rank r = entity_rank::invalid();

    std::ostringstream out;
    out << r;
    EXPECT_EQ(out.str(),"{entity_rank:invalid}");
  }
}



