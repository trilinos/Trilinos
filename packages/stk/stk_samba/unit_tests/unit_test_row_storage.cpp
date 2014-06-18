#include <gtest/gtest.h>

#include <samba/utility/row_storage.hpp>

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

TEST(samba, row_storage)
{
  int num_rows = 5;

  typedef samba::utility::row_storage<int> row_storage_type;

  std::vector<std::vector<int> > vec(num_rows);
  row_storage_type row_store(num_rows);

  std::ostream_iterator<int> out(std::cout, ", ");

  for (int i=0; i<num_rows; ++i) {
    vec[i].resize(2*i);
    for (int j=0; j<2*i; ++j) {
      vec[i][j] = j;
      row_store.insert_column(i,j,j);
    }
  }

  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }


  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }

  num_rows *= 2;

  vec.resize(num_rows);
  row_store.resize_rows(num_rows);

  for (int i=0; i<num_rows; ++i) {
    vec[i].resize(2*i);
    for (int j=0; j<2*i; ++j) {
      vec[i][j] = j;
      row_store.insert_column(i,j,j);
    }
  }

  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }

  row_store.compress();

  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }

  row_store.inflate();

  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }

  num_rows /= 2;

  vec.resize(num_rows);
  row_store.resize_rows(num_rows);

  for (int i=0; i<num_rows; ++i) {
    EXPECT_TRUE( std::equal(vec[i].begin(),vec[i].end(),row_store.begin(i)) );
    //std::copy(row_store.begin(i),row_store.end(i), out);
    //std::cout << "\b\b  " << std::endl;
  }
}


