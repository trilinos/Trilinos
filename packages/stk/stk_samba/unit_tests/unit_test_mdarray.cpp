#if 0
#include <gtest/gtest.h>

#include <samba/utility/mdarray.hpp>
#include <samba/utility/debug_message.hpp>

#include <iostream>

//DO NOT USE EXPECT_EQ/ASSERT_EQ WHEN TESTING META_FUNCTIONS
//the macro is broken and sometime thinks that a meta_fuction
//is a function pointer and results in hard-to-debug/misleading
//linker errors

TEST(samba, mdarray_reverse_array)
{
  typedef int array[1][2][3][4][5][6][7][8][9];

  typedef samba::detail::reverse_array<array>::type reverse_array;

  EXPECT_TRUE( (boost::extent<reverse_array,0>::value ==  9u) );
  EXPECT_TRUE( (boost::extent<reverse_array,1>::value ==  8u) );
  EXPECT_TRUE( (boost::extent<reverse_array,2>::value ==  7u) );
  EXPECT_TRUE( (boost::extent<reverse_array,3>::value ==  6u) );
  EXPECT_TRUE( (boost::extent<reverse_array,4>::value ==  5u) );
  EXPECT_TRUE( (boost::extent<reverse_array,5>::value ==  4u) );
  EXPECT_TRUE( (boost::extent<reverse_array,6>::value ==  3u) );
  EXPECT_TRUE( (boost::extent<reverse_array,7>::value ==  2u) );
  EXPECT_TRUE( (boost::extent<reverse_array,8>::value ==  1u) );

}

TEST(samba, mdarray_length)
{
  EXPECT_TRUE( (samba::detail::mdarray_length<int[2]>::value == 2u) );
  EXPECT_TRUE( (samba::detail::mdarray_length<int[2][3]>::value == 6u) );
  EXPECT_TRUE( (samba::detail::mdarray_length<int[2][3][4]>::value == 24u) );
  EXPECT_TRUE( (samba::detail::mdarray_length<int[][2][3][4]>::value == 0u) );
}


TEST(samba, mdarray_dimension_helper)
{
  typedef int array[1][2][3][4][5][6][7][8][9];

  EXPECT_EQ( samba::detail::dimension_helper<array>(0), 1u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(1), 2u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(2), 3u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(3), 4u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(4), 5u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(5), 6u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(6), 7u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(7), 8u);
  EXPECT_EQ( samba::detail::dimension_helper<array>(8), 9u);
}

TEST(samba, mdarray_basic)
{
  samba::mdarray<int[2][3][4]> array = {{0}};

  EXPECT_EQ(array.rank(),3u);

  EXPECT_EQ(array.length(),24u);

  EXPECT_EQ(array.stride(),12u);

  EXPECT_EQ(array.dimension(), 2u);
  EXPECT_EQ(array.size(), 2u);

  EXPECT_EQ(array.dimension(0), 2u);
  EXPECT_EQ(array.dimension(1), 3u);
  EXPECT_EQ(array.dimension(2), 4u);

  EXPECT_EQ(array.get_dimension<0>(), 2u);
  EXPECT_EQ(array.get_dimension<1>(), 3u);
  EXPECT_EQ(array.get_dimension<2>(), 4u);

  int count = 0;
  for (size_t i=0; i<array.get_dimension<0>(); ++i) {
  for (size_t j=0; j<array.get_dimension<1>(); ++j) {
  for (size_t k=0; k<array.get_dimension<2>(); ++k) {
    EXPECT_EQ(array[i][j][k],0);
    array[i][j][k] = ++count;
  }
  }
  }

  count = 0;
  for (size_t i=0; i<array.get_dimension<0>(); ++i) {
  for (size_t j=0; j<array.get_dimension<1>(); ++j) {
  for (size_t k=0; k<array.get_dimension<2>(); ++k) {
    EXPECT_EQ(array[i][j][k],++count);
  }
  }
  }


  //TODO after itertor is implemented test iterator
}

#endif

