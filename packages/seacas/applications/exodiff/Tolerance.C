// Copyright(C) 2008-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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
#include "Tolerance.h"
#include <cstdlib>     // for abs
#include <sys/types.h> // for int32_t, int64_t

namespace {
  /* See
     http://randomascii.wordpress.com/2012/01/11/tricks-with-the-floating-point-format
     for the potential portability problems with the union and bit-fields below.
  */
  union Float_t {
    explicit Float_t(float num = 0.0f) : f(num) {}
    // Portable extraction of components.
    bool Negative() const { return ((uint32_t)i >> 31) != 0; }

    int32_t i;
    float   f;
  };

  union Double_t {
    explicit Double_t(double num = 0.0) : f(num) {}
    // Portable extraction of components.
    bool Negative() const { return ((uint64_t)i >> 63) != 0; }

    int64_t i;
    double  f;
  };

  bool AlmostEqualUlpsFloat(float A, float B, int maxUlpsDiff)
  {
    Float_t uA(A);
    Float_t uB(B);

    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative()) {
      // Check for equality to make sure +0==-0
      return (A == B);
    }

    // Find the difference in ULPs.
    int ulpsDiff = std::abs(uA.i - uB.i);
    return (ulpsDiff <= maxUlpsDiff);
  }

  bool AlmostEqualUlpsDouble(double A, double B, int maxUlpsDiff)
  {
    Double_t uA(A);
    Double_t uB(B);

    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative()) {
      // Check for equality to make sure +0==-0
      return (A == B);
    }

    // Find the difference in ULPs.
    int ulpsDiff = std::abs(uA.i - uB.i);
    return (ulpsDiff <= maxUlpsDiff);
  }
} // namespace

bool Tolerance::use_old_floor = false;

bool Tolerance::Diff(double v1, double v2) const
{
  if (type == IGNORE_) {
    return false;
  }

  if (use_old_floor) {
    if (fabs(v1 - v2) < floor) {
      return false;
    }
  }
  else {
    if (fabs(v1) <= floor && fabs(v2) <= floor) {
      return false;
    }
  }

  if (type == RELATIVE_) {
    if (v1 == 0.0 && v2 == 0.0) {
      return false;
    }
    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    return fabs(v1 - v2) > value * max;
  }
  if (type == ABSOLUTE_) {
    return fabs(v1 - v2) > value;
  }
  if (type == COMBINED_) {
    // if (Abs(x - y) <= Max(absTol, relTol * Max(Abs(x), Abs(y))))
    // In the current implementation, absTol == relTol;
    // At some point, store both values...
    // In summary, use abs tolerance if both values are less than 1.0;
    // else use relative tolerance.

    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    double tol = 1.0 < max ? max : 1.0;
    return fabs(v1 - v2) >= tol * value;

    // EIGEN type diffs work on the absolute values. They are intended
    // for diffing eigenvectors where one shape may be the inverse of
    // the other and they should compare equal.  Eventually would like
    // to do a better check to ensure that ratio of one shape to other
    // is 1 or -1...
  }
  if (type == ULPS_FLOAT_) {
    return !AlmostEqualUlpsFloat(v1, v2, static_cast<int>(value));
  }
  else if (type == ULPS_DOUBLE_) {
    return !AlmostEqualUlpsDouble(v1, v2, static_cast<int>(value));
  }
  else if (type == EIGEN_REL_) {
    if (v1 == 0.0 && v2 == 0.0) {
      return false;
    }
    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    return fabs(fabs(v1) - fabs(v2)) > value * max;
  }
  else if (type == EIGEN_ABS_) {
    return fabs(fabs(v1) - fabs(v2)) > value;
  }
  else if (type == EIGEN_COM_) {
    // if (Abs(x - y) <= Max(absTol, relTol * Max(Abs(x), Abs(y))))
    // In the current implementation, absTol == relTol;
    // At some point, store both values...
    // In summary, use abs tolerance if both values are less than 1.0;
    // else use relative tolerance.

    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    double tol = 1.0 < max ? max : 1.0;
    return fabs(fabs(v1) - fabs(v2)) >= tol * value;
  }
  else {
    return true;
  }
}

const char *Tolerance::typestr() const
{
  if (type == RELATIVE_) {
    return "relative";
  }
  if (type == ABSOLUTE_) {
    return "absolute";
  }
  else if (type == COMBINED_) {
    return "combined";
  }
  else if (type == ULPS_FLOAT_) {
    return "ulps_float";
  }
  else if (type == ULPS_DOUBLE_) {
    return "ulps_double";
  }
  else if (type == EIGEN_REL_) {
    return "eigenrel";
  }
  else if (type == EIGEN_ABS_) {
    return "eigenabs";
  }
  else if (type == EIGEN_COM_) {
    return "eigencom";
  }
  else {
    return "ignore";
  }
}

const char *Tolerance::abrstr() const
{
  if (type == RELATIVE_) {
    return "rel";
  }
  if (type == ABSOLUTE_) {
    return "abs";
  }
  else if (type == COMBINED_) {
    return "com";
  }
  else if (type == ULPS_FLOAT_) {
    return "upf";
  }
  else if (type == ULPS_DOUBLE_) {
    return "upd";
  }
  else if (type == EIGEN_REL_) {
    return "ere";
  }
  else if (type == EIGEN_ABS_) {
    return "eab";
  }
  else if (type == EIGEN_COM_) {
    return "eco";
  }
  else {
    return "ign";
  }
}

double Tolerance::UlpsDiffFloat(double A, double B) const
{
  Float_t uA(A);
  Float_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative()) {
    // Check for equality to make sure +0==-0
    if (A == B) {
      return 0.0;
    }
    return 2 << 28;
  }

  // Find the difference in ULPs.
  return abs(uA.i - uB.i);
}

double Tolerance::UlpsDiffDouble(double A, double B) const
{
  Double_t uA(A);
  Double_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative()) {
    // Check for equality to make sure +0==-0
    if (A == B) {
      return 0.0;
    }
    return 2 << 28;
  }

  // Find the difference in ULPs.
  return std::abs(uA.i - uB.i);
}
