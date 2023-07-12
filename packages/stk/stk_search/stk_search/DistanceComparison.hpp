/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
//

#ifndef STK_SEARCH_UTIL_DISTANCE_COMPARISON_HPP
#define STK_SEARCH_UTIL_DISTANCE_COMPARISON_HPP

#include <cmath>
#include <limits>
#include <vector>

namespace stk {
namespace search {

static constexpr double STK_COMPARISON_EPSILON = .00000001;

inline bool less_than(const double d1, const double d2, const double eps = STK_COMPARISON_EPSILON)
{
  bool returnValue = false;
  if(std::numeric_limits<double>::max() == d2) {
    returnValue = true;
  }
  else if(std::numeric_limits<double>::max() == d1) {
    returnValue = false;
  }
  else if(d1 < d2 - eps) {
    returnValue = true;
  }
  else {
    returnValue = false;
  }

  return returnValue;
}

inline bool eq_to(const double d1, const double d2, const double eps = STK_COMPARISON_EPSILON)
{
  bool returnValue = false;
  if((std::numeric_limits<double>::max() == d2) && (std::numeric_limits<double>::max() == d1)) {
    returnValue = true;
  }
  else if((d1 > d2 - eps) && (d1 < d2 + eps)) {
    returnValue = true;
  }
  return returnValue;
}


inline double distance_sq(const unsigned spatialDim, const double* p1, const double* p2)
{
  double distance = 0.0;
  for(unsigned i = 0; i < spatialDim; ++i) {
    distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }

  return distance;
}

inline double distance(const unsigned spatialDim, const double* p1, const double* p2)
{
  double distanceSquared = distance_sq(spatialDim, p1, p2);

  return std::sqrt(distanceSquared);
}

inline bool centroid_less_than(const std::vector<double>& prevSendingCentroid,
                               const std::vector<double>& currSendingCentroid,
                               const std::vector<double>& recvCentroid, const double eps)
{
  bool returnValue = false;

  const double prevRealDistance =
      distance_sq(recvCentroid.size(), prevSendingCentroid.data(), recvCentroid.data());
  const double currRealDistance =
      distance_sq(recvCentroid.size(), currSendingCentroid.data(), recvCentroid.data());

  if(currRealDistance < prevRealDistance - eps) {
    returnValue = true;
  }
  else if(prevRealDistance < currRealDistance - eps) {
    returnValue = false;
  }
  else {
    returnValue = true;
  }
  return returnValue;
}

template <class SENDMESH, class RECVMESH>
inline bool
less_than_with_centroid_fallback(const double d1, const double d2,
                                 const typename SENDMESH::EntityKey previousSend,
                                 const typename SENDMESH::EntityKey currentSend,
                                 const typename RECVMESH::EntityKey recvObj,
                                 const SENDMESH* mesha, const RECVMESH* meshb,
                                 const double eps = STK_COMPARISON_EPSILON)
{
  bool returnValue = false;
  if(less_than(d1, d2, eps)) {
    returnValue = true;
  }
  else if(d2 < d1 - eps) {
    returnValue = false;
  }
  else if(mesha->is_valid(previousSend)) {
    // Gather the send object's nodal coordinates
    std::vector<double> prevSendingCentroid;
    std::vector<double> currSendingCentroid;
    std::vector<double> recvCentroid;

    mesha->centroid(previousSend, prevSendingCentroid);
    mesha->centroid(currentSend, currSendingCentroid);
    meshb->centroid(recvObj, recvCentroid);
    returnValue = centroid_less_than(prevSendingCentroid, currSendingCentroid, recvCentroid, eps);
  }
  else {
    returnValue = true;
  }

  return returnValue;
}

} // namespace search
} // namespace stk

#endif
