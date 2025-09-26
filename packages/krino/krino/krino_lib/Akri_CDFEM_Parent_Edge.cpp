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

double CDFEM_Parent_Edge::get_edge_node_position(stk::mesh::Entity edgeNode) const
{
  for (size_t i=0; i<my_edge_nodes.size(); ++i)
    if (my_edge_nodes[i] == edgeNode)
      return my_edge_node_positions[i];
  STK_ThrowRequireMsg(false, "Could not find edge node.");
  return -1.0;
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

std::vector<std::pair<const InterfaceID,double>*> get_sorted_internal_crossings(CrossingMap & crossings)
{
  std::vector<std::pair<const InterfaceID,double>*> sortedCrossings;
  sortedCrossings.reserve(crossings.size());
  for (auto && crossing : crossings)
    if (crossing.second > 0. && crossing.second < 1.)
      sortedCrossings.push_back(&crossing);

  std::sort(sortedCrossings.begin(), sortedCrossings.end(),
      [](const std::pair<const InterfaceID,double>* crossing0, const std::pair<const InterfaceID,double>* crossing1)
      { return crossing0->second < crossing1->second || (crossing0->second == crossing1->second && crossing0->first < crossing1->first); });
  return sortedCrossings;
}

void collapse_small_internal_segments_while_perserving_topology(CrossingMap &crossings, const double snapTol)
{
  const std::vector<std::pair<const InterfaceID,double>*> sortedCrossings = get_sorted_internal_crossings(crossings);

  bool done = false;
  while (!done)
  {
    size_t minSegment = 0;
    double minSegmentSize = snapTol;
    for (size_t i=1; i<sortedCrossings.size(); ++i)
    {
      const double segmentSize = (sortedCrossings[i]->second - sortedCrossings[i-1]->second);
      if (segmentSize > 0 && segmentSize <= minSegmentSize)
      {
        minSegment = i;
        minSegmentSize = segmentSize;
      }
    }
    done = minSegment == 0;
    if (!done)
    {
      double & loc0 = sortedCrossings[minSegment-1]->second;
      double & loc1 = sortedCrossings[minSegment]->second;
      const double newLoc = 0.5*(loc0+loc1);
      loc0 = newLoc;
      loc1 = newLoc;
    }
  }
}

void CDFEM_Parent_Edge::collapse_small_segments_while_preserving_topology(const double snapTol)
{
  for (auto && crossing : my_crossings)
  {
    if (crossing.second <= snapTol)
      crossing.second = 0.0;
    else if (crossing.second >= 1.-snapTol )
      crossing.second = 1.0;
  }
  collapse_small_internal_segments_while_perserving_topology(my_crossings, snapTol);
}

std::vector<int> get_effective_sorted_parent_node_domains(const int parentNodeIndex, const CrossingMap & crossings, const std::vector<int> & sortedParentNodeDomains)
{
  std::vector<int> effectiveSortedParentNodeDomains = sortedParentNodeDomains;
  double endPoint = (parentNodeIndex == 0) ? 0. : 1.;
  for (auto && crossing : crossings)
  {
    if (std::binary_search(sortedParentNodeDomains.begin(), sortedParentNodeDomains.end(), crossing.first.first_ls()) ||
        std::binary_search(sortedParentNodeDomains.begin(), sortedParentNodeDomains.end(), crossing.first.second_ls()))
    {
      endPoint = (parentNodeIndex == 0) ? std::max(endPoint, crossing.second) : std::min(endPoint, crossing.second);
    }
  }
  for (auto && crossing : crossings)
  {
    const bool inRange = (parentNodeIndex == 0) ? (crossing.second < endPoint) : (crossing.second > endPoint);
    if (inRange)
    {
      effectiveSortedParentNodeDomains.push_back(crossing.first.first_ls());
      effectiveSortedParentNodeDomains.push_back(crossing.first.second_ls());
    }
  }

  stk::util::sort_and_unique(effectiveSortedParentNodeDomains);
  return effectiveSortedParentNodeDomains;
}

void adjust_crossing_locations_based_on_node_captured_domains_for_level_set_per_interface(const std::vector<int> & sortedParentNode0Domains, const std::vector<int> & sortedParentNode1Domains, CrossingMap & crossings)
{
  if (sortedParentNode0Domains.empty() && sortedParentNode1Domains.empty())
      return;

  for (auto && crossing : crossings)
  {
    const InterfaceID & iface = crossing.first;
    STK_ThrowAssert(iface.first_ls() == iface.second_ls());
    const bool in0 = std::binary_search(sortedParentNode0Domains.begin(), sortedParentNode0Domains.end(), iface.first_ls());
    const bool in1 = std::binary_search(sortedParentNode1Domains.begin(), sortedParentNode1Domains.end(), iface.first_ls());

    double & loc = crossing.second;
    if (in0 && (!in1 || loc<0.5))
      loc = 0.;
    else if (in1 && (!in0 || loc>=0.5))
      loc = 1.;
  }
}

void adjust_crossing_locations_based_on_node_captured_domains_for_level_set_per_phase(const int parentNodeIndex, const std::vector<int> & sortedParentNodeDomains, CrossingMap & allCrossings)
{
  if (sortedParentNodeDomains.empty())
      return;

  STK_ThrowAssert(parentNodeIndex == 0 || parentNodeIndex == 1);
  const double nodePos = (parentNodeIndex == 0) ? 0. : 1.;

  double furthestLocationToAdjust = (parentNodeIndex == 0) ? 0. : 1.;
  for (auto && crossing : allCrossings)
  {
    const InterfaceID & iface = crossing.first;
    if (std::binary_search(sortedParentNodeDomains.begin(), sortedParentNodeDomains.end(), iface.first_ls()) &&
        std::binary_search(sortedParentNodeDomains.begin(), sortedParentNodeDomains.end(), iface.second_ls()))
    {
      furthestLocationToAdjust = (parentNodeIndex == 0) ? std::max(furthestLocationToAdjust, crossing.second) : std::min(furthestLocationToAdjust, crossing.second);
      crossing.second = nodePos;
    }
  }

  for (auto && crossing : allCrossings)
  {
    if ((parentNodeIndex == 0 && crossing.second < furthestLocationToAdjust) ||
        (parentNodeIndex == 1 && crossing.second > furthestLocationToAdjust))
    {
      crossing.second = nodePos;
    }
  }
}

void adjust_crossing_locations_based_on_node_captured_domains_for_level_set_per_phase(const std::vector<int> & sortedParentNode0Domains, const std::vector<int> & sortedParentNode1Domains, CrossingMap & allCrossings)
{
  if (sortedParentNode0Domains.empty() && sortedParentNode1Domains.empty())
      return;

  double highestLocToSendTo0 = 0.;
  double lowestLocToSendTo1 = 1.;
  for (auto && crossing : allCrossings)
  {
    const InterfaceID & iface = crossing.first;
    const bool in0 =
        std::binary_search(sortedParentNode0Domains.begin(), sortedParentNode0Domains.end(), iface.first_ls()) &&
        std::binary_search(sortedParentNode0Domains.begin(), sortedParentNode0Domains.end(), iface.second_ls());
    const bool in1 =
        std::binary_search(sortedParentNode1Domains.begin(), sortedParentNode1Domains.end(), iface.first_ls()) &&
        std::binary_search(sortedParentNode1Domains.begin(), sortedParentNode1Domains.end(), iface.second_ls());

    double & loc = crossing.second;
    if (in0 && (!in1 || loc<0.5))
    {
      highestLocToSendTo0 = std::max(loc, highestLocToSendTo0);
      loc = 0.;
    }
    if (in1 && (!in0 || loc>=0.5))
    {
      highestLocToSendTo0 = std::min(loc, lowestLocToSendTo1);
      loc = 1.;
    }
  }

  for (auto && crossing : allCrossings)
  {
    double & loc = crossing.second;
    if (loc < highestLocToSendTo0)
      loc = 0.;
    if (loc > lowestLocToSendTo1)
      loc = 1.;
  }
}

void copy_real_crossing_locations_to_fake_crossing_locations(const CrossingMap & realCrossings, CrossingMap & allCrossings)
{
  for (auto && crossing : allCrossings)
  {
    auto iter = realCrossings.find(crossing.first);
    if (iter != realCrossings.end())
      crossing.second = iter->second;
  }
}

std::pair<int,int> get_begin_and_end_phases(const CrossingMap & realCrossings, const CrossingSignMap & realCrossingSign)
{
  double begLoc = 1.;
  double endLoc = 0.;
  int begPhase = -1;
  int endPhase = -1;
  for (auto && crossing : realCrossings)
  {
    const InterfaceID iface = crossing.first;
    const double loc = crossing.second;
    const int sign = realCrossingSign.at(iface);
    const int fromPhase = (sign == 1) ? iface.first_ls() : iface.second_ls();
    if (loc < begLoc || (loc == begLoc && fromPhase > begPhase))
    {
      begPhase = fromPhase;
      begLoc = loc;
    }
    const int toPhase = (sign == -1) ? iface.first_ls() : iface.second_ls();
    if (loc > endLoc || (loc == endLoc && toPhase > endPhase))
    {
      endPhase = toPhase;
      endLoc = loc;
    }
  }
  return {begPhase,endPhase};
}

static double find_next_crossing_location(const CrossingMap & crossings,
    const CrossingSignMap & crossingSigns,
    const double currentLocation,
    const int currentPhase)
{
  double nextLocation = 1.;
  for (auto && crossing : crossings)
  {
    if (crossing.second >= currentLocation)
    {
      const InterfaceID iface = crossing.first;
      const int fromPhase = (crossingSigns.at(iface) == 1) ? iface.first_ls() : iface.second_ls();
      if (fromPhase == currentPhase)
        nextLocation = std::min(nextLocation, crossing.second);
    }
  }
  return nextLocation;
}

static std::vector<int> find_next_phase_candidates_at_location(const CrossingMap & crossings,
    const CrossingSignMap & crossingSigns,
    const double nextLocation,
    const int currentPhase)
{
  std::vector<int> nextPhaseCandidates;
  for (auto && crossing : crossings)
  {
    if (crossing.second == nextLocation)
    {
      const InterfaceID iface = crossing.first;
      const int sign = crossingSigns.at(iface);
      const int fromPhase = (sign == 1) ? iface.first_ls() : iface.second_ls();
      if (fromPhase == currentPhase)
      {
        const int toPhase = (sign == -1) ? iface.first_ls() : iface.second_ls();
        nextPhaseCandidates.push_back(toPhase);
      }
    }
  }
  return nextPhaseCandidates;
}

std::vector<int> find_next_phase_candidates(const CrossingMap & crossings,
    const CrossingSignMap & crossingSigns,
    const double currentLocation,
    const int currentPhase)
{
  const double nextLocation = find_next_crossing_location(crossings, crossingSigns, currentLocation, currentPhase);
  return find_next_phase_candidates_at_location(crossings, crossingSigns, nextLocation, currentPhase);
}

std::vector<std::pair<InterfaceID,double>> shortest_path_to_end(const std::vector<std::pair<InterfaceID,double>> & pathSoFar,
    const CrossingMap & crossings,
    const CrossingSignMap & crossingSigns,
    const double currentLocation,
    const int currentPhase,
    const int endPhase)
{
  const std::vector<int> nextPhaseCandidates = find_next_phase_candidates(crossings, crossingSigns, currentLocation, currentPhase);
  if (nextPhaseCandidates.empty())
  {
    if (currentPhase == endPhase)
      return pathSoFar;
    else
      return {};
  }

  std::vector<std::pair<InterfaceID,double>> shortestPath;
  size_t shortestPathSize = std::numeric_limits<size_t>::max();
  for (int nextPhase : nextPhaseCandidates)
  {
    std::vector<std::pair<InterfaceID,double>> path = pathSoFar;
    const InterfaceID iface(currentPhase, nextPhase);
    const double nextLocation = crossings.at(iface);
    path.emplace_back(iface, nextLocation);
    const auto fullPath = shortest_path_to_end(path, crossings, crossingSigns, nextLocation, nextPhase, endPhase);
    if (!fullPath.empty() && fullPath.size() < shortestPathSize)
    {
      shortestPath = fullPath;
      shortestPathSize = fullPath.size();
    }
  }
  return shortestPath;
}

std::vector<std::pair<InterfaceID,double>> shortest_path_from_begin_to_end(const CrossingMap & crossings,
    const CrossingSignMap & crossingSigns,
    const int beginPhase,
    const int endPhase)
{
  std::vector<std::pair<InterfaceID,double>> emptyPath;
  double beginLocation = -1.;
  return shortest_path_to_end(emptyPath, crossings, crossingSigns, beginLocation, beginPhase, endPhase);
}

bool determine_real_crossings_from_locations(CrossingMap & realCrossings, CrossingSignMap & realCrossingSigns, std::set<int> & edgePhases, const CrossingMap & allCrossings, const CrossingSignMap & allCrossingSigns)
{
  const auto begEndPhases = get_begin_and_end_phases(realCrossings, realCrossingSigns);

  if (begEndPhases.first == begEndPhases.second)
  {
    return true;
  }

  realCrossings.clear();
  realCrossingSigns.clear();
  edgePhases.clear();

  const auto shortestPath = shortest_path_from_begin_to_end(allCrossings, allCrossingSigns, begEndPhases.first, begEndPhases.second);
  if (shortestPath.empty())
    return false;

  for (auto && crossing : shortestPath)
  {
    realCrossings[crossing.first] = crossing.second;
    realCrossingSigns[crossing.first] = allCrossingSigns.at(crossing.first);
    edgePhases.insert(crossing.first.first_ls());
    edgePhases.insert(crossing.first.second_ls());
  }
  return true;
}

void CDFEM_Parent_Edge::adjust_crossing_locations_based_on_node_captured_domains(const bool oneLSPerPhase, const std::vector<int> & sortedParentNode0Domains, const std::vector<int> & sortedParentNode1Domains)
{
  if (sortedParentNode0Domains.empty() && sortedParentNode1Domains.empty())
    return;
  if (oneLSPerPhase)
  {
    if (my_crossings_including_fake.empty())
      return;

    adjust_crossing_locations_based_on_node_captured_domains_for_level_set_per_phase(sortedParentNode0Domains, sortedParentNode1Domains, my_crossings_including_fake);

    const bool success = determine_real_crossings_from_locations(my_crossings, my_crossing_signs, edge_phases, my_crossings_including_fake, my_crossing_signs_including_fake);
    if (!success)
      krinolog << "Failed to adjust crossings " << *this << stk::diag::dendl;
  }
  else
  {
    if (my_crossings.empty())
      return;
    adjust_crossing_locations_based_on_node_captured_domains_for_level_set_per_interface(sortedParentNode0Domains, sortedParentNode1Domains, my_crossings);
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
  return krino::find_crossing_position_and_sign(key, my_edge_node_positions, nodes_isovar);
}

} // namespace krino
