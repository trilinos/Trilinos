// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <math.h>

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_LowerEnvelope.hpp>
#include <Akri_Sign.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <limits>

namespace krino{

CDFEM_Parent_Edge::CDFEM_Parent_Edge(const std::vector<stk::mesh::Entity> & edgeNodes,
    const std::vector<double> & edgeNodePositions)
: my_edge_nodes(edgeNodes),
  my_edge_node_positions(edgeNodePositions)
{
}

double CDFEM_Parent_Edge::MinSize()
{
  return SegmentLowerEnvelope::MinSize();
}

void
CDFEM_Parent_Edge::debug_print_crossings() const
{
  krinolog << " crossings: { ";
  for ( CrossingMap::const_iterator it=my_crossings.begin(); it != my_crossings.end(); ++it )
  {
    krinolog << it->second << " ";
  }
  krinolog << "}" << "\n";
}

void
CDFEM_Parent_Edge::find_crossings(const bool oneLSPerPhase, const std::vector<std::vector<double> > & nodes_isovar)
{
  my_crossings.clear();
  my_crossing_signs.clear();

  const int num_nodes = get_num_nodes();
  STK_ThrowAssert(static_cast<int>(nodes_isovar.size()) == num_nodes);
  const int num_ls = nodes_isovar[0].size();
  if (num_ls > 1 && oneLSPerPhase)
  {
    find_crossings_multiple_levelset(nodes_isovar);
    find_crossings_including_fake_ones(oneLSPerPhase, nodes_isovar);
    return;
  }

  // TODO: respect minimum_internal_edge_size

  for ( int ls_index = 0; ls_index < num_ls; ++ls_index )
  {
    InterfaceID iface(ls_index, ls_index);
      my_crossing_signs[iface] = sign(nodes_isovar[num_nodes-1][ls_index]);
    if( !sign_change(nodes_isovar[0][ls_index], nodes_isovar[num_nodes-1][ls_index]) ) continue;
    for ( int s = 0; s < num_nodes-1; ++s )
    {
      const double ls0 = nodes_isovar[s][ls_index];
      const double ls1 = nodes_isovar[s+1][ls_index];
      if ( sign_change(ls0, ls1) )
      {
        const double interval_position = ls0 / ( ls0 - ls1 );
        const double abs_position = (1.-interval_position)*my_edge_node_positions[s] + interval_position*my_edge_node_positions[s+1];
        my_crossings[iface] = abs_position;
        my_crossing_signs[iface] = sign(ls1);
      }
    }
  }

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    debug_print_crossings();
  }
}

void
CDFEM_Parent_Edge::find_crossings_multiple_levelset(const std::vector<std::vector<double> > & nodes_isovar)
{
  const Segment_Vector lower_envelope = SegmentLowerEnvelope::find_lower_envelope(my_edge_node_positions, nodes_isovar);

  edge_phases.clear();
  for(auto && it : lower_envelope)
  {
    edge_phases.insert(it.ls_index());
  }

  // Create crossings between the segments of the lower_envelope
  for(Segment_Vector::const_iterator it = lower_envelope.begin(); it != lower_envelope.end()-1; ++it)
  {
    const LS_Segment & cur = *it;
    const LS_Segment & next = *(it+1);
    if (cur.ls_index() == next.ls_index()) continue;

    const double cur_right = cur.right_endpoint();
    const double next_left = next.left_endpoint();
    STK_ThrowRequire(cur_right == next_left);
    STK_ThrowRequire(cur.ls_index() != next.ls_index());
    InterfaceID iface(cur.ls_index(), next.ls_index());
    STK_ThrowRequireMsg(my_crossings.find(iface) == my_crossings.end(), "Multiple interface crossing error after pruning.");
    my_crossings[iface] = cur_right;
    my_crossing_signs[iface] = (cur.ls_index() < next.ls_index()) ? 1 : -1;
  }

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    debug_print_crossings();
  }
}

bool
CDFEM_Parent_Edge::have_any_crossings() const
{
  return !my_crossings.empty();
}

std::tuple<double, int, bool>
CDFEM_Parent_Edge::get_crossing_position_and_sign(const bool oneLSPerPhase, const InterfaceID key) const
{
  STK_ThrowRequire(oneLSPerPhase);

  if(have_crossing(key))
  {
    return std::make_tuple(get_crossing_position(key), get_crossing_sign(key), false);
  }

  return std::make_tuple(get_fake_crossing_position(key), get_fake_crossing_sign(key), true);
}

typedef std::array<const std::pair<const InterfaceID,double>*,2> CrossingInterval;

static CrossingInterval get_crossing_interval(const CDFEM_Parent_Edge & edge, const std::pair<const InterfaceID,double> & fakeCrossing)
{
  const double fakeLoc = fakeCrossing.second;
  const std::pair<const InterfaceID,double> * before = nullptr;
  const std::pair<const InterfaceID,double> * after = nullptr;
  for (auto && crossing : edge.get_crossings())
  {
    const double loc = crossing.second;
    if (loc == fakeLoc)
    {
      return {&crossing, &crossing};
    }
    else if (loc < fakeLoc)
    {
      if (before == nullptr || loc > before->second)
        before = &crossing;
    }
    else
    {
      if (after == nullptr || loc < after->second)
        after = &crossing;
    }
  }
  return {before, after};
}

static int get_phase_on_interval(const CDFEM_Parent_Edge & edge, const CrossingInterval & crossingInterval)
{
  const std::pair<const InterfaceID,double> * before = crossingInterval[0];
  const std::pair<const InterfaceID,double> * after = crossingInterval[1];
  STK_ThrowRequire(before != nullptr || after != nullptr);
  if (before != nullptr)
  {
    return (edge.get_crossing_sign(before->first) == -1) ? before->first.first_ls() : before->first.second_ls();
  }
  return (edge.get_crossing_sign(after->first) == -1) ? after->first.second_ls() : after->first.first_ls();
}

static bool fake_crossing_is_actually_real(const CDFEM_Parent_Edge & edge, const CrossingInterval & crossingInterval, const std::pair<const InterfaceID,double> & fakeCrossing)
{
  if (crossingInterval[0] != crossingInterval[1])
  {
    const int phaseOnInterval = get_phase_on_interval(edge, crossingInterval);
    const InterfaceID fakeInterface = fakeCrossing.first;
    const bool crossingIsActuallyReal = (fakeInterface.first_ls() == phaseOnInterval || fakeInterface.second_ls() == phaseOnInterval);
    return crossingIsActuallyReal;
  }
  return false;
}

static void fixup_fake_crossing_location_for_consistency(CDFEM_Parent_Edge & edge, std::pair<const InterfaceID,double> & fakeCrossing)
{
  const CrossingInterval & crossingInterval = get_crossing_interval(edge, fakeCrossing);
  if (fake_crossing_is_actually_real(edge, crossingInterval, fakeCrossing))
  {
    const std::pair<const InterfaceID,double> * before = crossingInterval[0];
    const std::pair<const InterfaceID,double> * after = crossingInterval[1];
    double & loc = fakeCrossing.second;
    const bool useBefore = before ? (after ? (loc-before->second < after->second-loc) : true) : false;
    loc = useBefore ? before->second : after->second;
  }
}

static bool fake_crossing_is_actually_real(const CDFEM_Parent_Edge & edge, const std::pair<const InterfaceID,double> & fakeCrossing)
{
  const CrossingInterval & crossingInterval = get_crossing_interval(edge, fakeCrossing);
  return fake_crossing_is_actually_real(edge, crossingInterval, fakeCrossing);
}

bool CDFEM_Parent_Edge::all_fake_crossings_are_really_fake() const
{
  for (auto && crossing : my_crossings_including_fake)
    if (!have_crossing(crossing.first) && fake_crossing_is_actually_real(*this, crossing))
      return false;
  return true;
}

void CDFEM_Parent_Edge::fixup_fake_crossing_locations_for_consistency()
{
  for (auto && crossing : my_crossings_including_fake)
    if (!have_crossing(crossing.first))
      fixup_fake_crossing_location_for_consistency(*this, crossing);
}

void
CDFEM_Parent_Edge::find_crossings_including_fake_ones(const bool oneLSPerPhase, const std::vector<std::vector<double> > & nodes_isovar)
{
  my_crossings_including_fake.clear();
  my_crossing_signs_including_fake.clear();

  STK_ThrowRequire(oneLSPerPhase);
  const int numLS = nodes_isovar[0].size();
  std::vector<double> lsMins(numLS,std::numeric_limits<double>::max());
  std::vector<double> lsMaxs(numLS,std::numeric_limits<double>::lowest());
  for (auto && nodeIsovar : nodes_isovar)
  {
    for (int i=0; i<numLS; ++i)
    {
      lsMins[i] = std::min(lsMins[i], nodeIsovar[i]);
      lsMaxs[i] = std::max(lsMaxs[i], nodeIsovar[i]);
    }
  }

  for (int ls1=0; ls1<numLS; ++ls1)
  {
    for (int ls2=ls1; ls2<numLS; ++ls2)
    {
      if (!(lsMaxs[ls1] < lsMins[ls2] || lsMaxs[ls2] < lsMins[ls1]))
      {
        InterfaceID iface(ls1, ls2);
        const std::pair<double, int> result = have_crossing(iface) ?
            std::make_pair(get_crossing_position(iface), get_crossing_sign(iface)) :
            find_crossing_position_and_sign(oneLSPerPhase, iface, nodes_isovar);
        if (result.first >= 0.)
        {
          my_crossings_including_fake[iface] = result.first;
          my_crossing_signs_including_fake[iface] = result.second;
        }
      }
    }
  }

  fixup_fake_crossing_locations_for_consistency();
}

std::pair<double, int>
CDFEM_Parent_Edge::find_crossing_position_and_sign(const bool oneLSPerPhase, const InterfaceID key, const std::vector<std::vector<double> > & nodes_isovar) const
{
  STK_ThrowRequire(oneLSPerPhase);
  double crossingPos = -1.;
  int crossingSign = 0;
  krino::find_crossing_position_and_sign(key.first_ls(), key.second_ls(), my_edge_node_positions, nodes_isovar, crossingPos, crossingSign);
  return std::make_pair(crossingPos, crossingSign);
}

} // namespace krino
