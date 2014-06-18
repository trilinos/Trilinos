#include <gtest/gtest.h>

#include <samba/entity_part.hpp>

#include <vector>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>


TEST(samba, entity_part)
{
  using namespace samba;

  {
    entity_part t = entity_block_key::create(0);

    std::ostringstream out;
    out << t;
    EXPECT_EQ(out.str(),"{entity_block:0}");
  }

  {
    entity_part t = entity_rank::node();

    std::ostringstream out;
    out << t;
    EXPECT_EQ(out.str(),"{entity_rank:node}");
  }


  {
    std::vector<entity_part> vec;

    vec.push_back( entity_block_key::create(0));
    vec.push_back( entity_state::ghosted());
    vec.push_back( entity_topology::hex_8());
    vec.push_back( entity_block_key::create(3));
    vec.push_back( entity_rank::element());
    vec.push_back( entity_block_key::create(2));
    vec.push_back( entity_rank::node());
    vec.push_back( entity_topology::quad_4());
    vec.push_back( entity_state::universe());
    vec.push_back( entity_state::modified());
    vec.push_back( entity_block_key::create(1));

    std::sort(vec.begin(), vec.end());


    std::vector<entity_part> gold;

    gold.push_back( entity_rank::node());
    gold.push_back( entity_rank::element());
    gold.push_back( entity_topology::quad_4());
    gold.push_back( entity_topology::hex_8());
    gold.push_back( entity_state::universe());
    gold.push_back( entity_state::modified());
    gold.push_back( entity_state::ghosted());
    gold.push_back( entity_block_key::create(0));
    gold.push_back( entity_block_key::create(1));
    gold.push_back( entity_block_key::create(2));
    gold.push_back( entity_block_key::create(3));

    EXPECT_TRUE(vec == gold);

  }

}

