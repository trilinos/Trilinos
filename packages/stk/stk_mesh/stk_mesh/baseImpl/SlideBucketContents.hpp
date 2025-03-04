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

#ifndef stk_mesh_impl_SlideBucketContents_hpp
#define stk_mesh_impl_SlideBucketContents_hpp

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <vector>

namespace stk {
namespace mesh {
namespace impl {

template<typename T> T& get_bucket(std::vector<T>& buckets, size_t i) { return buckets[i]; }
template<typename T> T& get_bucket(std::vector<T*>& buckets, size_t i) { return *buckets[i]; }

template<typename BucketType,
         class IsHole, class Overwrite, class RemoveLast>
void slide_contents_to_fill_holes(unsigned holeBkt, unsigned holeOrd,
                                  std::vector<BucketType>& buckets,
                                  const IsHole& is_hole,
                                  const Overwrite& overwrite,
                                  const RemoveLast& remove_last)
{
  STK_ThrowRequireMsg(holeBkt < buckets.size(),
                      "holeBkt("<<holeBkt<<") out-of-range, buckets.size()="<<buckets.size());
  STK_ThrowRequireMsg(holeOrd < get_bucket(buckets,holeBkt).size(),
                      "holeOrd("<<holeOrd<<") out-of-range, buckets["<<holeBkt<<"].size()="<<get_bucket(buckets,holeBkt).size());
  STK_ThrowRequire(is_hole(get_bucket(buckets,holeBkt)[holeOrd]));

  unsigned bktIdx = holeBkt;
  unsigned bktOrd = holeOrd+1;
  if (bktOrd >= get_bucket(buckets,bktIdx).size()) {
    ++bktIdx;
    bktOrd = 0;
  }
  if (bktIdx >= buckets.size()) {
    remove_last();
    return;
  }

  unsigned numHolesFilled = 1;
  while(bktIdx < buckets.size()) {
    auto& bkt = get_bucket(buckets,bktIdx);
    while(bktOrd < bkt.size()) {
      if(!is_hole(bkt[bktOrd])) {
        overwrite(holeBkt, holeOrd, bktIdx, bktOrd);
        ++holeOrd;
        if (holeOrd >= get_bucket(buckets,holeBkt).size()) {
          ++holeBkt;
          holeOrd = 0;
        }
      }
      else {
        ++numHolesFilled;
      }
      ++bktOrd;
    }
    ++bktIdx;
    bktOrd = 0;
  }

  for(unsigned i=0; i<numHolesFilled; ++i) {
    remove_last();
  }
}

}}} // end namepsace stk mesh impl

#endif
