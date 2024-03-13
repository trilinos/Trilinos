// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

namespace MueLu {

size_t LevenshteinDistance(const char* s, size_t len_s, const char* t, size_t len_t) {
  // degenerate cases
  if (len_s == 0) return len_t;
  if (len_t == 0) return len_s;
  if (!strncmp(s, t, std::min(len_s, len_t))) return 0;

  // create two work vectors of integer distances
  size_t len = len_t + 1;
  std::vector<size_t> v0(len);
  std::vector<size_t> v1(len);

  // initialize v0 (the previous row of distances)
  // this row is A[0][i]: edit distance for an empty s
  // the distance is just the number of characters to delete from t
  for (size_t i = 0; i < len; i++)
    v0[i] = i;

  for (size_t i = 0; i < len_s; i++) {
    // calculate v1 (current row distances) from the previous row v0

    // first element of v1 is A[i+1][0]
    //   edit distance is delete (i+1) chars from s to match empty t
    v1[0] = i + 1;

    // use formula to fill in the rest of the row
    for (size_t j = 0; j < len_t; j++) {
      size_t cost = (s[i] == t[j]) ? 0 : 1;
      v1[j + 1]   = std::min(v1[j] + 1,
                             std::min(v0[j + 1] + 1,
                                      v0[j] + cost));
    }

    // copy v1 (current row) to v0 (previous row) for next iteration
    for (size_t j = 0; j < len; j++)
      v0[j] = v1[j];
  }

  return v1[len_t];
}

}  // namespace MueLu
