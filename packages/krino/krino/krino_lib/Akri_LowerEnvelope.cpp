// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <math.h>

#include <Akri_MeshHelpers.hpp>
#include <Akri_LowerEnvelope.hpp>
#include <Akri_Sign.hpp>

#include <limits>

namespace krino{

LS_Segment::LS_Segment(const double xL, const double xR, int id_)
: x_left(xL), x_right(xR), LS_index(id_)
{
}

int
get_min(const std::vector<double> & phi)
{
  double min_value = 0.0;
  int min_index = -1;
  for (unsigned i=0; i<phi.size(); ++i)
  {
    // <= is here to return the max index when multiple phi values are identical
    if (min_index < 0 || phi[i] <= min_value)
    {
      min_index = i;
      min_value = phi[i];
    }
  }
  return min_index;
}


std::pair<double, int>
find_crossing_position_and_sign(const InterfaceID key, const std::vector<double> & pos, const std::vector<std::vector<double>> & phi)
{
  const int num_nodes = pos.size();
  std::vector<double> key_isovar(num_nodes);
  for ( int node = 0; node < num_nodes; ++node )
  {
    key_isovar[node] = phi[node][key.first_ls()]-phi[node][key.second_ls()];
  }

  if ( sign_change(key_isovar[0], key_isovar[num_nodes-1]) )
  {
    for ( int s = 0; s < num_nodes-1; ++s )
    {
      const double ls0 = key_isovar[s];
      const double ls1 = key_isovar[s+1];
      if ( sign_change(ls0, ls1) )
      {
        const double interval_position = ls0 / ( ls0 - ls1 );
        const double abs_position = (1.-interval_position)*pos[s] + interval_position*pos[s+1];
        return std::make_pair(abs_position, sign(ls1));
      }
    }
  }

  const int node0_sign = sign(key_isovar[0]);
  const int node1_sign = sign(key_isovar[num_nodes-1]);
  STK_ThrowRequire(node0_sign == node1_sign);
  return std::make_pair(-1., node0_sign);
}

bool
find_lower_envelope_crossing_point( const std::array<std::vector<double>,2> & phi, double & x)
{
  const int id0 = get_min(phi[0]);
  const int id1 = get_min(phi[1]);

  if ( id0 != id1 )
  {
    const double ls0 = (phi[0][id0]-phi[0][id1]);
    const double ls1 = (phi[1][id0]-phi[1][id1]);

    x = ls0 / ( ls0 - ls1 );
    return true;
  }

  return false;
}

bool
find_lower_envelope_crossing_point(const std::array<std::vector<double>,3> & phi, std::array<double,2> & point)
{
  const int id0 = get_min(phi[0]);
  const int id1 = get_min(phi[1]);
  const int id2 = get_min(phi[2]);

  if (id0 != id1 && id1 != id2 && id0 != id2)
  {
    const double d10 = (phi[0][id1]-phi[0][id0]);
    const double d11 = (phi[1][id1]-phi[1][id0]);
    const double d12 = (phi[2][id1]-phi[2][id0]);
    const double d20 = (phi[0][id2]-phi[0][id0]);
    const double d21 = (phi[1][id2]-phi[1][id0]);
    const double d22 = (phi[2][id2]-phi[2][id0]);

    const double a00 = d11-d10;
    const double a01 = d12-d10;
    const double a10 = d21-d20;
    const double a11 = d22-d20;

    const double b0 = -d10;
    const double b1 = -d20;

    const double det = a00*a11-a01*a10;
    if (det == 0.) return false;

    const double x = (b0*a11-b1*a01)/det;
    const double y = (-b0*a10+b1*a00)/det;
    if (x < 0. || y < 0. || 1.-x-y < 0.) return false;

    point = {{x,y}};
    return true;
  }

  return false;
}

bool
find_lower_envelope_crossing_point(const std::array<std::vector<double>,4> & phi, std::array<double,3> & point)
{
  const int id0 = get_min(phi[0]);
  const int id1 = get_min(phi[1]);
  const int id2 = get_min(phi[2]);
  const int id3 = get_min(phi[3]);

  if (id0 != id1 && id0 != id2 && id0 != id3 && id1 != id2 && id1 != id3 && id2 != id3)
  {
    const double d10 = (phi[0][id1]-phi[0][id0]);
    const double d11 = (phi[1][id1]-phi[1][id0]);
    const double d12 = (phi[2][id1]-phi[2][id0]);
    const double d13 = (phi[3][id1]-phi[3][id0]);
    const double d20 = (phi[0][id2]-phi[0][id0]);
    const double d21 = (phi[1][id2]-phi[1][id0]);
    const double d22 = (phi[2][id2]-phi[2][id0]);
    const double d23 = (phi[3][id2]-phi[3][id0]);
    const double d30 = (phi[0][id3]-phi[0][id0]);
    const double d31 = (phi[1][id3]-phi[1][id0]);
    const double d32 = (phi[2][id3]-phi[2][id0]);
    const double d33 = (phi[3][id3]-phi[3][id0]);

    const double a00 = d11-d10;
    const double a01 = d12-d10;
    const double a02 = d13-d10;
    const double a10 = d21-d20;
    const double a11 = d22-d20;
    const double a12 = d23-d20;
    const double a20 = d31-d30;
    const double a21 = d32-d30;
    const double a22 = d33-d30;
    const double b0 = -d10;
    const double b1 = -d20;
    const double b2 = -d30;
    const double det =  a00*(a22*a11-a21*a12)-a10*(a22*a01-a21*a02)+a20*(a12*a01-a11*a02);
    const double x =( b0*(a22*a11-a21*a12)-b1*(a22*a01-a21*a02)+b2*(a12*a01-a11*a02))/det;
    const double y =(-b0*(a22*a10-a20*a12)+b1*(a22*a00-a20*a02)-b2*(a12*a00-a10*a02))/det;
    const double z =( b0*(a21*a10-a20*a11)-b1*(a21*a00-a20*a01)+b2*(a11*a00-a10*a01))/det;

    point = {{x,y,z}};
    return true;
  }

  return false;
}

void
SegmentLowerEnvelope::adjust_piecewise_linear_positions(Segment_Vector & segments, const std::vector<double> & pos, const std::vector<std::vector<double>> & phi)
{
  if (pos.size() > 2)
  {
    // Even for piecewise edges, the structure of the phases must be no more complicated than that
    // given by the linear version of the edge.  So adjust the linear locations and then remove
    // the degenerate, inverted segements.

    // Use piecewise approximation to find actual transition locations
    for(Segment_Vector::iterator it = segments.begin(); it != segments.end()-1; ++it)
    {
      LS_Segment & cur = *it;
      LS_Segment & next = *(it+1);
      if (cur.ls_index() == next.ls_index()) continue;
      InterfaceID iface(cur.ls_index(), next.ls_index());

      const double crossing_pos = (find_crossing_position_and_sign(iface, pos, phi)).first;
      STK_ThrowRequire(crossing_pos >= 0.);
      cur = LS_Segment(cur.left_endpoint(), crossing_pos, cur.ls_index());
      next = LS_Segment(crossing_pos, next.right_endpoint(), next.ls_index());
    }

    // Now remove inverted segments
    if (segments.size() > 2)
    {
      for(Segment_Vector::iterator it = segments.begin()+1; it != segments.end()-1; ++it)
      {
        LS_Segment & prev = *(it-1);
        LS_Segment & cur = *it;
        LS_Segment & next = *(it+1);

        if (cur.right_endpoint() < cur.left_endpoint())
        {
          InterfaceID iface(prev.ls_index(), next.ls_index());

          const double crossing_pos = (find_crossing_position_and_sign(iface, pos, phi)).first;
          STK_ThrowRequire(crossing_pos >= 0.);
          prev = LS_Segment(prev.left_endpoint(), cur.right_endpoint(), prev.ls_index());
          cur = LS_Segment(cur.right_endpoint(), crossing_pos, prev.ls_index());
          next = LS_Segment(crossing_pos, next.right_endpoint(), next.ls_index());
        }
      }
    }
  }
}

void
SegmentLowerEnvelope::collapse_identical_segments(Segment_Vector & segments)
{
  if (segments.size() > 1)
  {
    Segment_Vector uncollapsed;
    uncollapsed.swap(segments);
    segments.push_back(uncollapsed.front());
    for(Segment_Vector::const_iterator it = uncollapsed.begin()+1; it != uncollapsed.end(); ++it)
    {
      const LS_Segment & prev = *(it-1);
      const LS_Segment & curr = *it;
      if (prev.ls_index() == curr.ls_index())
      {
        segments.back() = LS_Segment(segments.back().left_endpoint(), curr.right_endpoint(), prev.ls_index());
      }
      else
      {
        segments.push_back(curr);
      }
    }
  }
}

void
SegmentLowerEnvelope::add_segment(Segment_Vector & segments,
    const double x0, const std::vector<double> & phi0, const int min0,
    const double x1, const std::vector<double> & phi1, const int min1)
{
  if (min0 == min1)
  {
    const LS_Segment segment(x0, x1, min0);
    segments.push_back(segment);
    return;
  }
  else
  {
    InterfaceID key(min0, min1);
    const double ls0 = (phi0[key.first_ls()]-phi0[key.second_ls()]); // ls >= 0 for phi2 >= phi1
    const double ls1 = (phi1[key.first_ls()]-phi1[key.second_ls()]);
    STK_ThrowRequire(ls0*ls1 <= 0.);
    const double cut = ls0 / (ls0 - ls1);

    {
      const double x_cut = (1.-cut)*x0 + cut*x1;
      STK_ThrowAssert(phi0.size() == phi1.size());
      std::vector<double> phiAtCut(phi0.size());
      for (unsigned i=0; i<phiAtCut.size(); ++i)
        phiAtCut[i] = (1.-cut)*(phi0[i]) + cut*(phi1[i]);
      const int minAtCut = get_min(phiAtCut);

      const bool minCutMatchesEndPoint = minAtCut == min0 || minAtCut == min1;
      const double tol = MinSize()/(x1-x0);
      const bool isInfinitesmalCutAt0 = cut < tol;
      const bool isInfinitesmalCutAt1 = cut > 1.-tol;

      if (isInfinitesmalCutAt0)
      {
        if (segments.empty())
          segments.emplace_back(0., 0., min0);
        if (minCutMatchesEndPoint)
          segments.emplace_back(x0, x1, min1);
        else
          add_segment(segments, x0, phiAtCut, minAtCut, x1, phi1, min1);
      }
      else if (isInfinitesmalCutAt1)
      {
        if (minCutMatchesEndPoint)
          segments.emplace_back(x0, x1, min0);
        else
          add_segment(segments, x0, phi0, min0, x1, phiAtCut, minAtCut);
        if (x1 == 1.)
        {
          if (segments.back().left_endpoint() == 1.)
            segments.pop_back();
          segments.emplace_back(1., 1., min1);
        }
      }
      else
      {
        if (minCutMatchesEndPoint)
          segments.emplace_back(x0, x_cut, min0);
        else
          add_segment(segments, x0, phi0, min0, x_cut, phiAtCut, minAtCut);

        if (minCutMatchesEndPoint)
          segments.emplace_back(x_cut, x1, min1);
        else
          add_segment(segments, x_cut, phiAtCut, minAtCut, x1, phi1, min1);
      }
    }
  }
}

} // namespace krino
