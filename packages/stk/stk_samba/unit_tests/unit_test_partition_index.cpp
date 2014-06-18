#include <gtest/gtest.h>

#include <samba/utility/debug_message.hpp>
#include <samba/partition_index.hpp>

#include <iostream>
#include <sstream>


TEST(samba, partition_index)
{
  using namespace samba;

  partition_index a;

  a = partition_index::invalid();

  EXPECT_EQ(a.rank(), entity_rank::invalid());
  EXPECT_EQ(a.partition(), partition_id::invalid());
  EXPECT_EQ(a.offset(), partition_offset::invalid());

  entity_rank r = {10};
  partition_id s = {20};
  partition_offset o = {30};

  a = partition_index::create(r,s,o);

  EXPECT_EQ(a.rank(),r);
  EXPECT_EQ(a.partition(),s);
  EXPECT_EQ(a.offset(),o);

}


