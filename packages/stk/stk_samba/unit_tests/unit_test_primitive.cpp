#include <gtest/gtest.h>
#include <boost/mpi.hpp>

#include <samba/utility/macros.hpp>
#include <samba/utility/make_primitive.hpp>
#include <samba/utility/hash_value.hpp>
#include <samba/utility/primitive_output.hpp>
#include <samba/utility/debug_message.hpp>

#include <boost/type_traits.hpp>

#include <iostream>
#include <sstream>

//put inside of samba namespace to test
//argument dependant lookup (ADL)
//for operators
namespace samba {

SAMBA_MAKE_PRIMITIVE(A,A,int);
SAMBA_MAKE_PRIMITIVE_WITH_HASHABLE_TAG(B,B,int,1);
SAMBA_MAKE_PRIMITIVE_WITH_HASHABLE_TAG(C,C,int,2);

} //namespace samba

//the following macros must be called outside any
//namespace without a trailing semi-colon
SAMBA_IS_PRIMITIVE(samba::A)
SAMBA_IS_PRIMITIVE(samba::B)
SAMBA_IS_PRIMITIVE(samba::C)


TEST(samba, primitive_hash_value)
{
  using namespace samba;

  A a = {0};
  B b = {0};
  C c = {0};

  EXPECT_NE( hash_value(a), hash_value(b));
  EXPECT_NE( hash_value(a), hash_value(c));
  EXPECT_NE( hash_value(b), hash_value(c));
}

TEST(samba, primitive_output)
{
  using namespace samba;

  A a = {0};

  std::ostringstream out;
  out << a;
  EXPECT_EQ(out.str(),"{A:0}");
}

TEST(samba, primitive_mpi)
{
  using namespace samba;

  boost::mpi::communicator world;
  if (world.size() > 1) {
    A msg;
    const A gold = A::create(100);

    if (world.rank() == 0) {
      msg = A::create(100);

      world.send(world.size()-1,0,msg);
    }

    if (world.rank()==world.size()-1)
    {
      world.recv(0,0,msg);

      EXPECT_EQ(msg(),gold());
    }
  }
}


