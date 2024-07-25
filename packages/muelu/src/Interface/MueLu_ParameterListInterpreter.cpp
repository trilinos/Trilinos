// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
