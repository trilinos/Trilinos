// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_mesh/baseImpl/SlideBucketContents.hpp>
#include <vector>

namespace {

class TestBkt {
public:
  TestBkt(int n) : m_x(n) {}
  size_t size() const { return 1; }
  const int& operator[](size_t /*i*/) const { return m_x; }
private:
  int m_x;
};

TEST(GetBucket, byValue)
{
  std::vector<TestBkt> vec = {TestBkt(9)};
  EXPECT_EQ(9, stk::mesh::impl::get_bucket(vec,0)[0]);
}

TEST(GetBucket, byPtr)
{
  std::vector<TestBkt*> vec;
  TestBkt bkt(9);
  vec.push_back(&bkt);

  EXPECT_EQ(9, stk::mesh::impl::get_bucket(vec,0)[0]);
}

template<typename BucketType>
void call_the_slide_function(unsigned holeBkt, unsigned holeOrd,
                             std::vector<BucketType>& buckets)
{
  auto is_hole = [&](int item){return item == -1;};

  auto overwrite =
    [&](unsigned destBkt, unsigned destOrd, unsigned srcBkt, unsigned srcOrd)
    {buckets[destBkt][destOrd] = buckets[srcBkt][srcOrd];};

  auto remove_last =
    [&]()
    {
      BucketType& lastBucket = buckets.back();
      lastBucket.pop_back();
      if (lastBucket.empty()) {
        buckets.pop_back();
      }
    };

  stk::mesh::impl::slide_contents_to_fill_holes(
    holeBkt, holeOrd, buckets,
    is_hole, overwrite, remove_last
  );
}

TEST(UnitTestSlide, slide_bucket_contents_nominal)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {0, 1, -1, 3}, {6, -1, -1, 10}, {52, -1, 54, 55, 56, 57} };

  std::vector<BucketType> finalBuckets =
  { {0, 1, 3, 6}, {10, 52, 54, 55}, {56, 57} };

  call_the_slide_function(0, 2, buckets);

  EXPECT_EQ(buckets, finalBuckets);
}

TEST(UnitTestSlide, slide_bucket_contents_holeAtEndOfBkt)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {0, -1, 3, -1}, {6, 8, -1, 13}, {52, 53, 54, 55} };

  std::vector<BucketType> finalBuckets =
  { {0, 3, 6, 8}, {13, 52, 53, 54}, {55} };

  call_the_slide_function(0, 1, buckets);

  EXPECT_EQ(buckets, finalBuckets);
}

TEST(UnitTestSlide, slide_bucket_contents_holeAtStartOfBkt)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {-1, -1, 3, -1}, {-1, 8, -1, 13}, {52, 53, 54, 55, 56, 57} };

  std::vector<BucketType> finalBuckets =
  { {3, 8, 13, 52}, {53, 54, 55, 56}, {57} };

  call_the_slide_function(0, 0, buckets);

  EXPECT_EQ(buckets, finalBuckets);
}

TEST(UnitTestSlide, slide_bucket_contents_removeMiddleBucket)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {0, 1}, {-1, -1, -1}, {52, 53} };

  std::vector<BucketType> finalBuckets =
  { {0, 1}, {52, 53} };

  call_the_slide_function(1, 0, buckets);

  EXPECT_EQ(buckets, finalBuckets);
}

TEST(UnitTestSlide, slide_bucket_contents_removeEndBucket)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {0, 1}, {52, 53}, {-1, -1, -1} };

  std::vector<BucketType> finalBuckets =
  { {0, 1}, {52, 53} };

  call_the_slide_function(2, 0, buckets);

  EXPECT_EQ(buckets, finalBuckets);
}

TEST(UnitTestSlide, slide_bucket_contents_removeOnlyBucket)
{
  using BucketType = std::vector<int>;
  std::vector<BucketType> buckets =
  { {-1, -1, -1} };

  call_the_slide_function(0, 0, buckets);

  EXPECT_TRUE(buckets.empty());
}

} // empty namespace
