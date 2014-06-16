#include <gtest/gtest.h>

#include <samba/utility.hpp>
#include <iostream>

namespace samba {

SAMBA_MAKE_PRIMITIVE(primitive,tag,int);

} //namespace

//the following macros must be called outside any
//namespace without a trailing semi-colon
SAMBA_IS_PRIMITIVE(samba::primitive)

TEST(samba, pod_support)
{
  using namespace boost;

  EXPECT_TRUE( is_pod<samba::primitive>::value );

}

TEST(samba, nothrow)
{
  using namespace boost;

  EXPECT_TRUE( has_nothrow_copy<samba::primitive>::value );

  EXPECT_TRUE( has_nothrow_constructor<samba::primitive>::value );
}

TEST(samba, is_primitive)
{
  using namespace samba;
  EXPECT_TRUE( is_primitive<samba::primitive>::value );
}

#ifdef SAMBA_ENABLE_PARALLEL
TEST(samba, is_bitwise_serializable)
{
  using namespace boost::serialization;

  EXPECT_TRUE( is_bitwise_serializable<samba::primitive>::value );
}

TEST(samba, is_mpi_datatype)
{
  using namespace boost::mpi;

  EXPECT_TRUE( is_mpi_datatype<samba::primitive>::value );
}
#endif



