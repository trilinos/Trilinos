#include <gtest/gtest.h>

#include <samba/entity_key.hpp>

#include <iostream>
#include <sstream>


TEST(samba, entity_key)
{
  using namespace samba;

  entity_key a;

  a = entity_key::invalid();

  EXPECT_EQ(a.topology(), entity_topology::invalid());
  EXPECT_EQ(a.process(), process_id::invalid());
  EXPECT_EQ(a.local_id(), entity_local_id::invalid());

  entity_topology t = {10};
  process_id s = {20};
  entity_local_id o = {30};

  a = entity_key::create(t,s,o);

  EXPECT_EQ(a.topology(),t);
  EXPECT_EQ(a.process(),s);
  EXPECT_EQ(a.local_id(),o);

}

