 // Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef stk_util_util_RemoveIntersection_hpp
#define stk_util_util_RemoveIntersection_hpp

#include <vector>
#include <algorithm>

namespace stk {
namespace util {

template<typename T1, typename T2, typename COMPARE=std::less<T1> >
void remove_intersection(std::vector<T1>& v1, std::vector<T2>& v2, COMPARE comp=std::less<T1>())
{
    std::vector<T1> intersection(std::min(v1.size(), v2.size()));
    typename std::vector<T1>::iterator it;

    std::sort(v1.begin(), v1.end(), comp);
    std::sort(v2.begin(), v2.end(), comp);

    it=std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), intersection.begin(), comp);
    intersection.resize(it-intersection.begin());

    for (T1 item : intersection) {
        typename std::vector<T1>::iterator found1 = std::lower_bound(v1.begin(), v1.end(), item, comp);
        if (*found1 == item) {
            v1.erase(found1);
        }
        typename std::vector<T2>::iterator found2 = std::lower_bound(v2.begin(), v2.end(), item, comp);
        if (*found2 == item) {
            v2.erase(found2);
        }
    }
}

} //namespace util
} //namespace stk

#endif
