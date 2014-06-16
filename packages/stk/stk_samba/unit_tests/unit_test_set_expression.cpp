#include <gtest/gtest.h>

#include <samba/set_expression.hpp>

#include <vector>
#include <sstream>

TEST(samba, set_expression)
{
  using namespace samba;

  {
    set_expression expr;

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"{entity_state:invalid}");
  }

  {
    set_expression expr(entity_rank::node());

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"{entity_rank:node}");
  }

  {
    set_expression expr;

    expr = entity_rank::node();

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"{entity_rank:node}");
  }

  {

    entity_rank rank = entity_rank::node();
    set_expression expr;
    expr = rank;

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"{entity_rank:node}");
  }

  {
    set_expression expr =  (entity_rank::node() & entity_topology::hex_8());

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"({entity_rank:node}&{entity_topology:hexahedron_8})");
  }

  {
    entity_rank rank = entity_rank::node();
    entity_topology topo = entity_topology::hex_8();

    set_expression expr =  (rank & topo);

    std::ostringstream out;
    out << expr;
    EXPECT_EQ(out.str(),"({entity_rank:node}&{entity_topology:hexahedron_8})");
  }

  {
    entity_block_key block_1 = entity_block_key::create(1);
    entity_block_key block_2 = entity_block_key::create(2);
    entity_block_key block_3 = entity_block_key::create(3);

    set_expression expr =  ( (entity_rank::node() & entity_topology::hex_8())
                            & (block_1 | block_2 | block_3) )
                          - entity_state::ghosted();

    std::ostringstream out;
    out << expr;
    EXPECT_EQ( out.str()
              ,"((({entity_rank:node}&{entity_topology:hexahedron_8})&(({entity_block:1}|{entity_block:2})|{entity_block:3}))-{entity_state:ghosted})"
             );
  }
}

TEST(samba, set_expression_contains)
{
  using namespace samba;

  entity_block_key block_1 = entity_block_key::create(1);
  entity_block_key block_2 = entity_block_key::create(2);
  entity_block_key block_3 = entity_block_key::create(3);

  set_expression expr =  ( (entity_rank::node() & entity_topology::hex_8())
                          & (block_1 | block_2 | block_3) )
                        - entity_state::ghosted();

  {
    entity_part_vector partition;

    partition.push_back(entity_part(entity_rank::node()));
    partition.push_back(entity_part(entity_topology::hex_8()));
    partition.push_back(entity_part(block_1));

    EXPECT_TRUE(contains(expr,partition));
  }

  {
    entity_part_vector partition;

    partition.push_back(entity_part(entity_rank::node()));
    partition.push_back(entity_part(entity_topology::hex_8()));
    partition.push_back(entity_part(entity_state::ghosted()));
    partition.push_back(entity_part(block_1));

    EXPECT_FALSE(contains(expr,partition));
  }
}
