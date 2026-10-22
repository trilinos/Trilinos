// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_LowerEnvelope_h
#define Akri_LowerEnvelope_h

#include <stk_util/util/ReportHandler.hpp>

#include <array>
#include <vector>

namespace krino {

class LS_Segment
{
public:
  LS_Segment(const double xL, const double xR, int id);

  double length() const { return x_right - x_left; }
  int ls_index() const { return LS_index; }
  double left_endpoint() const { return x_left; }
  double right_endpoint() const { return x_right; }

  friend std::ostream & operator << (std::ostream &os, const LS_Segment & segment);

private:
  double x_left;
  double x_right;
  int LS_index;
};
typedef std::vector<LS_Segment> Segment_Vector;

inline std::ostream & operator << (std::ostream &os, const LS_Segment & segment)
{
  os << "segment LS id = " << segment.LS_index
      << " left point = " << segment.x_left
      << " right point = " << segment.x_right;
  return os;
}

inline std::ostream & operator << (std::ostream &os, const Segment_Vector & segmentVec)
{
  os << "Segments: ";
  for (auto && segment : segmentVec)
    os << "{" << segment << "} ";
  return os;
}

// Returns true and fills crossingPos (absolute position along pos) and crossingSign
// (the sign of the level-set difference just past the crossing) when the level-set pair
// (lsA, lsB) crosses somewhere along the edge. The key difference is phi[lsFirst]-phi[lsSecond]
// where lsFirst/lsSecond are the min/max of (lsA, lsB), so the result is independent of the
// order in which the two indices are passed. Returns false when there is no crossing, in
// which case crossingPos is left unchanged and crossingSign is set to the common sign shared
// by both end nodes.
bool
find_crossing_position_and_sign(const int lsA, const int lsB, const std::vector<double> & pos, const std::vector<std::vector<double>> & phi, double & crossingPos, int & crossingSign);

bool find_lower_envelope_crossing_point(const std::array<std::vector<double>,2> & phi, double & x);

bool find_lower_envelope_crossing_point(const std::array<std::vector<double>,3> & phi, std::array<double,2> & point);

bool find_lower_envelope_crossing_point(const std::array<std::vector<double>,4> & phi, std::array<double,3> & point);

int get_min(const std::vector<double> & phi);

class SegmentLowerEnvelope
{
public:
  // Must be larger than machine epsilon, but puts limit on smallest effective snap tolerance
  static constexpr double MinSize() {return 1.e-12;};

  // linear (2 points)
  static Segment_Vector find_lower_envelope(const std::vector<double> & phi0,
      const std::vector<double> & phi1) { SegmentLowerEnvelope env(0., phi0, 1., phi1); return env.get_segments(); }
  static Segment_Vector find_lower_envelope(const double x0, const std::vector<double> & phi0,
    const double x1, const std::vector<double> & phi1) { SegmentLowerEnvelope env(x0, phi0, x1, phi1); return env.get_segments(); }
  // piecewise linear (>2 points)
  static Segment_Vector find_lower_envelope(const std::vector<double> & pos, const std::vector<std::vector<double>> & phi)
    { SegmentLowerEnvelope env(pos,phi); return env.get_segments(); }

  SegmentLowerEnvelope(const std::vector<double> & phi0, const std::vector<double> & phi1) : SegmentLowerEnvelope(0., phi0, 1., phi1) {}
  SegmentLowerEnvelope(const double x0, const std::vector<double> & phi0,
    const double x1, const std::vector<double> & phi1)
  {
    add_segment(mySegments, x0, phi0, x1, phi1);
    collapse_identical_segments(mySegments);
  }
  SegmentLowerEnvelope(const std::vector<double> & pos, const std::vector<std::vector<double>> & phi)
  {
    STK_ThrowRequire(pos.size() == phi.size() && pos.size() > 1);
    add_segment(mySegments, pos.front(), phi.front(), pos.back(), phi.back());
    collapse_identical_segments(mySegments);
    adjust_piecewise_linear_positions(mySegments, pos, phi);
  }
  const Segment_Vector & get_segments() const { return mySegments; }

private:
  // xEnd is the absolute position of the right end of the whole edge, retained
  // unchanged through the recursion so the infinitesimal-cut handling can tell
  // when the current subinterval reaches the true right end of the edge.
  static void add_segment(Segment_Vector & segments,
    const double x0, const std::vector<double> & phi0,
    const double x1, const std::vector<double> & phi1)
  {
    const int min0 = get_min(phi0);
    const int min1 = get_min(phi1);
    return add_segment(segments, x0, phi0, min0, x1, phi1, min1, x1);
  }

  static void add_segment(Segment_Vector & segments,
    const double x0, const std::vector<double> & phi0, const int min0,
    const double x1, const std::vector<double> & phi1, const int min1,
    const double xEnd);

  static void collapse_identical_segments(Segment_Vector & segments);

  static void adjust_piecewise_linear_positions(Segment_Vector & segments, const std::vector<double> & pos, const std::vector<std::vector<double>> & phi);

  Segment_Vector mySegments;
};

} // namespace krino

#endif // Akri_LowerEnvelope_h
