// Copyright(C) 1999-2017 National Technology & Engineering Solutions
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

#include <Ioss_CodeTypes.h>
#include <algorithm>
#include <cgns/Iocgns_StructuredZoneData.h>

#define OUTPUT std::cerr

namespace {
  struct Range
  {
    Range(int a, int b) : m_beg(a < b ? a : b), m_end(a < b ? b : a), m_reversed(b < a) {}

    int  begin() const { return m_reversed ? m_end : m_beg; }
    int  end() const { return m_reversed ? m_beg : m_end; }
    int  m_beg;
    int  m_end;
    bool m_reversed;
  };

  bool overlaps(const Range &a, const Range &b) { return a.m_beg <= b.m_end && b.m_beg <= a.m_end; }

  bool zgc_overlaps(const Iocgns::StructuredZoneData *zone, const Ioss::ZoneConnectivity &zgc)
  {
    // Note that zone range is nodes and m_ordinal[] is cells, so need to add 1 to range.
    Range z_i(1 + zone->m_offset[0], zone->m_ordinal[0] + zone->m_offset[0] + 1);
    Range z_j(1 + zone->m_offset[1], zone->m_ordinal[1] + zone->m_offset[1] + 1);
    Range z_k(1 + zone->m_offset[2], zone->m_ordinal[2] + zone->m_offset[2] + 1);

    Range gc_i(zgc.m_ownerRangeBeg[0], zgc.m_ownerRangeEnd[0]);
    Range gc_j(zgc.m_ownerRangeBeg[1], zgc.m_ownerRangeEnd[1]);
    Range gc_k(zgc.m_ownerRangeBeg[2], zgc.m_ownerRangeEnd[2]);

    return overlaps(z_i, gc_i) && overlaps(z_j, gc_j) && overlaps(z_k, gc_k);
  }

  bool zgc_donor_overlaps(const Iocgns::StructuredZoneData *zone, const Ioss::ZoneConnectivity &zgc)
  {
    // Note that zone range is nodes and m_ordinal[] is cells, so need to add 1 to range.
    Range z_i(1 + zone->m_offset[0], zone->m_ordinal[0] + zone->m_offset[0] + 1);
    Range z_j(1 + zone->m_offset[1], zone->m_ordinal[1] + zone->m_offset[1] + 1);
    Range z_k(1 + zone->m_offset[2], zone->m_ordinal[2] + zone->m_offset[2] + 1);

    Range gc_i(zgc.m_donorRangeBeg[0], zgc.m_donorRangeEnd[0]);
    Range gc_j(zgc.m_donorRangeBeg[1], zgc.m_donorRangeEnd[1]);
    Range gc_k(zgc.m_donorRangeBeg[2], zgc.m_donorRangeEnd[2]);

    return overlaps(z_i, gc_i) && overlaps(z_j, gc_j) && overlaps(z_k, gc_k);
  }

  Range subset_range(const Range &a, const Range &b)
  {
    Range ret(std::max(a.m_beg, b.m_beg), std::min(a.m_end, b.m_end));
    ret.m_reversed = a.m_reversed || b.m_reversed;
    return ret;
  }

  void zgc_subset_ranges(const Iocgns::StructuredZoneData *zone, Ioss::ZoneConnectivity &zgc)
  {
    if (zgc_overlaps(zone, zgc)) {
      // NOTE: Updates the range and donor_range in zgc

      // Note that zone range is nodes and m_ordinal[] is cells, so need to add 1 to range.
      Range z_i(1 + zone->m_offset[0], zone->m_ordinal[0] + zone->m_offset[0] + 1);
      Range z_j(1 + zone->m_offset[1], zone->m_ordinal[1] + zone->m_offset[1] + 1);
      Range z_k(1 + zone->m_offset[2], zone->m_ordinal[2] + zone->m_offset[2] + 1);

      Range gc_i(zgc.m_ownerRangeBeg[0], zgc.m_ownerRangeEnd[0]);
      Range gc_j(zgc.m_ownerRangeBeg[1], zgc.m_ownerRangeEnd[1]);
      Range gc_k(zgc.m_ownerRangeBeg[2], zgc.m_ownerRangeEnd[2]);

      Range gc_ii = subset_range(z_i, gc_i);
      Range gc_jj = subset_range(z_j, gc_j);
      Range gc_kk = subset_range(z_k, gc_k);

      Ioss::IJK_t range_beg;
      Ioss::IJK_t range_end;
      range_beg[0] = gc_ii.begin();
      range_end[0] = gc_ii.end();
      range_beg[1] = gc_jj.begin();
      range_end[1] = gc_jj.end();
      range_beg[2] = gc_kk.begin();
      range_end[2] = gc_kk.end();

      if (zgc.m_sameRange) {
        zgc.m_ownerRangeBeg = range_beg;
        zgc.m_ownerRangeEnd = range_end;
        zgc.m_donorRangeBeg = zgc.m_ownerRangeBeg;
        zgc.m_donorRangeEnd = zgc.m_ownerRangeEnd;
      }
      else {
        auto d_range_beg    = zgc.transform(range_beg);
        zgc.m_donorRangeEnd = zgc.transform(range_end);
        zgc.m_donorRangeBeg = d_range_beg;
        zgc.m_ownerRangeBeg = range_beg;
        zgc.m_ownerRangeEnd = range_end;
      }
      zgc.m_ownerOffset = {{zone->m_offset[0], zone->m_offset[1], zone->m_offset[2]}};
    }
    else {
      // This zgc does not overlap on this zone, so set all ranges to 0.
      // Still need it in list so can write block out correctly in parallel...
      zgc.m_ownerRangeBeg = {{0, 0, 0}};
      zgc.m_ownerRangeEnd = {{0, 0, 0}};
      zgc.m_donorRangeBeg = {{0, 0, 0}};
      zgc.m_donorRangeEnd = {{0, 0, 0}};
      zgc.m_isActive      = false;
    }
  }

  void zgc_subset_donor_ranges(const Iocgns::StructuredZoneData *don_zone,
                               Ioss::ZoneConnectivity &          zgc)
  {
    // NOTE: Updates the range and donor_range in zgc

    // Note that zone range is nodes and m_ordinal[] is cells, so need to add 1 to range.
    Range z_i(1 + don_zone->m_offset[0], don_zone->m_ordinal[0] + don_zone->m_offset[0] + 1);
    Range z_j(1 + don_zone->m_offset[1], don_zone->m_ordinal[1] + don_zone->m_offset[1] + 1);
    Range z_k(1 + don_zone->m_offset[2], don_zone->m_ordinal[2] + don_zone->m_offset[2] + 1);

    Range gc_i(zgc.m_donorRangeBeg[0], zgc.m_donorRangeEnd[0]);
    Range gc_j(zgc.m_donorRangeBeg[1], zgc.m_donorRangeEnd[1]);
    Range gc_k(zgc.m_donorRangeBeg[2], zgc.m_donorRangeEnd[2]);

    Range gc_ii = subset_range(z_i, gc_i);
    Range gc_jj = subset_range(z_j, gc_j);
    Range gc_kk = subset_range(z_k, gc_k);

    Ioss::IJK_t d_range_beg;
    Ioss::IJK_t d_range_end;
    d_range_beg[0] = gc_ii.begin();
    d_range_end[0] = gc_ii.end();
    d_range_beg[1] = gc_jj.begin();
    d_range_end[1] = gc_jj.end();
    d_range_beg[2] = gc_kk.begin();
    d_range_end[2] = gc_kk.end();

    if (zgc.m_sameRange) {
      zgc.m_donorRangeBeg = d_range_beg;
      zgc.m_donorRangeEnd = d_range_end;
      zgc.m_ownerRangeEnd      = zgc.m_donorRangeEnd;
      zgc.m_ownerRangeBeg      = zgc.m_donorRangeBeg;
    }
    else {
      auto range_beg      = zgc.inverse_transform(d_range_beg);
      zgc.m_ownerRangeEnd      = zgc.inverse_transform(d_range_end);
      zgc.m_ownerRangeBeg      = range_beg;
      zgc.m_donorRangeBeg = d_range_beg;
      zgc.m_donorRangeEnd = d_range_end;
    }
    zgc.m_donorOffset = {{don_zone->m_offset[0], don_zone->m_offset[1], don_zone->m_offset[2]}};
  }

  void propogate_zgc(Iocgns::StructuredZoneData *parent, Iocgns::StructuredZoneData *child,
                     int ordinal, int rank)
  {
    if (rank == 0) {
#if 0 && IOSS_DEBUG_OUTPUT
    OUTPUT << "\t\tPropogating ZGC from " << parent->m_name << " to " << child->m_name << "\n";
#endif
    }
    for (auto zgc : parent->m_zoneConnectivity) {
      if (!zgc.m_intraBlock || zgc_overlaps(child, zgc)) {
        // Modify source and donor range to subset it to new block ranges.
        zgc_subset_ranges(child, zgc);
        child->m_zoneConnectivity.push_back(zgc);
        if (rank == 0) {
#if 0 && IOSS_DEBUG_OUTPUT
        OUTPUT << "\t\tAdding: Donor Zone " << zgc.m_donorName << ":\tConnection Name '" << zgc.m_connectionName
               << "\n";
#endif
        }
      }
      else {
        if (rank == 0) {
#if 0 && IOSS_DEBUG_OUTPUT
        OUTPUT << "\t\tDonor Zone " << zgc.m_donorName << ":\tConnection Name '" << zgc.m_connectionName
               << " does not overlap."
               << "\n";
#endif
        }
      }
    }
  }

  // Add the zgc corresponding to the new communication path between
  // two child zones arising from a parent split along ordinal 'ordinal'
  void add_split_zgc(Iocgns::StructuredZoneData *parent, Iocgns::StructuredZoneData *c1,
                     Iocgns::StructuredZoneData *c2, int ordinal)
  {
    Ioss::IJK_t transform{{1, 2, 3}};

    // Note that range is specified in terms of 'adam' block i,j,k
    // space which is converted to local block i,j,k space
    // via the m_offset[] field on the local block.
    Ioss::IJK_t range_beg{{1 + c1->m_offset[0], 1 + c1->m_offset[1], 1 + c1->m_offset[2]}};
    Ioss::IJK_t range_end{{c1->m_ordinal[0] + c1->m_offset[0] + 1,
                           c1->m_ordinal[1] + c1->m_offset[1] + 1,
                           c1->m_ordinal[2] + c1->m_offset[2] + 1}};

    Ioss::IJK_t donor_range_beg(range_beg);
    Ioss::IJK_t donor_range_end(range_end);

    donor_range_end[ordinal] = donor_range_beg[ordinal] = range_beg[ordinal] = range_end[ordinal];

    auto c1_base = std::to_string(c1->m_adam->m_zone) + "_" + std::to_string(c1->m_zone);
    auto c2_base = std::to_string(c2->m_adam->m_zone) + "_" + std::to_string(c2->m_zone);

    const auto &adam_name = parent->m_adam->m_name;
    c1->m_zoneConnectivity.emplace_back(c1_base + "--" + c2_base, c1->m_zone, adam_name, c2->m_zone,
                                        transform, range_beg, range_end, donor_range_beg,
                                        donor_range_end, c1->m_zone < c2->m_zone, true);
    auto &zgc1 = c1->m_zoneConnectivity.back();
    zgc1.m_sameRange = true;
    zgc1.m_ownerOffset = {{c1->m_offset[0], c1->m_offset[1], c1->m_offset[2]}};
    zgc1.m_donorOffset = {{c2->m_offset[0], c2->m_offset[1], c2->m_offset[2]}};

    // OUTPUT << "Adding c1 " << c1_base << "--" << c2_base << "\n";
    // OUTPUT << c1->m_zoneConnectivity.back() << "\n";

    c2->m_zoneConnectivity.emplace_back(c2_base + "--" + c1_base, c2->m_zone, adam_name, c1->m_zone,
                                        transform, donor_range_beg, donor_range_end, range_beg,
                                        range_end, c2->m_zone < c1->m_zone, true);
    auto &zgc2 = c2->m_zoneConnectivity.back();
    zgc2.m_sameRange = true;
    zgc2.m_ownerOffset = {{c2->m_offset[0], c2->m_offset[1], c2->m_offset[2]}};
    zgc2.m_donorOffset = {{c1->m_offset[0], c1->m_offset[1], c1->m_offset[2]}};
    // OUTPUT << "Adding c2 " << c2_base << "--" << c1_base << "\n";
    // OUTPUT << c2->m_zoneConnectivity.back() << "\n";
  }
} // namespace

namespace Iocgns {

  // ========================================================================
  // Split this StructuredZone along the largest ordinal
  // into two children and return the created zones.
  std::pair<StructuredZoneData *, StructuredZoneData *>
  StructuredZoneData::split(int zone_id, double ratio, int rank)
  {
    assert(is_active());
    if (ratio > 1.0) {
      ratio = 1.0 / ratio;
    }

    // Find ordinal with largest value... Split along that ordinal
    int ordinal = 0;
    if (m_preferentialOrdinal == ordinal) {
      ordinal = 1;
    }
    if (m_ordinal[1] > m_ordinal[ordinal] && m_preferentialOrdinal != 1) {
      ordinal = 1;
    }
    if (m_ordinal[2] > m_ordinal[ordinal] && m_preferentialOrdinal != 2) {
      ordinal = 2;
    }

    assert(ordinal != m_preferentialOrdinal);

    if (m_ordinal[ordinal] <= 1) {
      return std::make_pair(nullptr, nullptr);
    }

    m_child1 = new StructuredZoneData;
    m_child2 = new StructuredZoneData;

    m_child1->m_name             = m_name + "_c1";
    m_child1->m_ordinal          = m_ordinal;
    m_child1->m_ordinal[ordinal] = m_child1->m_ordinal[ordinal] * ratio;
    if (m_child1->m_ordinal[ordinal] == 0) {
      m_child1->m_ordinal[ordinal] = 1;
    }
    m_child1->m_offset = m_offset; // Child1 offsets the same as parent;

    m_child1->m_preferentialOrdinal = m_preferentialOrdinal;
    m_child1->m_zone                = zone_id++;
    m_child1->m_adam                = m_adam;
    m_child1->m_parent              = this;
    m_child1->m_splitOrdinal        = ordinal;
    m_child1->m_sibling             = m_child2;

    m_child2->m_name             = m_name + "_c2";
    m_child2->m_ordinal          = m_ordinal;
    m_child2->m_ordinal[ordinal] = m_ordinal[ordinal] - m_child1->m_ordinal[ordinal];
    assert(m_child2->m_ordinal[ordinal] > 0);
    m_child2->m_offset = m_offset;
    m_child2->m_offset[ordinal] += m_child1->m_ordinal[ordinal];

    m_child2->m_preferentialOrdinal = m_preferentialOrdinal;
    m_child2->m_zone                = zone_id++;
    m_child2->m_adam                = m_adam;
    m_child2->m_parent              = this;
    m_child2->m_splitOrdinal        = ordinal;
    m_child2->m_sibling             = m_child1;

    if (rank == 0) {
#if IOSS_DEBUG_OUTPUT
      OUTPUT << "\nSplit Zone " << m_name << " (" << m_zone << ") Adam " << m_adam->m_name << " ("
             << m_adam->m_zone << ") with intervals " << m_ordinal[0] << " " << m_ordinal[1] << " "
             << m_ordinal[2] << " work = " << work() << " with offset " << m_offset[0] << " "
             << m_offset[1] << " " << m_offset[2] << " along ordinal " << ordinal << " with ratio "
             << ratio << "\n"
             << "\tChild 1: Zone " << m_child1->m_name << " (" << m_child1->m_zone << ") Adam "
             << m_child1->m_adam->m_name << " (" << m_child1->m_adam->m_zone << ") with intervals "
             << m_child1->m_ordinal[0] << " " << m_child1->m_ordinal[1] << " "
             << m_child1->m_ordinal[2] << " work = " << m_child1->work() << " with offset "
             << m_child1->m_offset[0] << " " << m_child1->m_offset[1] << " "
             << m_child1->m_offset[2] << "\n"
             << "\tChild 2: Zone " << m_child2->m_name << " (" << m_child2->m_zone << ") Adam "
             << m_child2->m_adam->m_name << " (" << m_child2->m_adam->m_zone << ") with intervals "
             << m_child2->m_ordinal[0] << " " << m_child2->m_ordinal[1] << " "
             << m_child2->m_ordinal[2] << " work = " << m_child2->work() << " with offset "
             << m_child2->m_offset[0] << " " << m_child2->m_offset[1] << " "
             << m_child2->m_offset[2] << "\n";
#endif
    }

    // Add ZoneGridConnectivity instance to account for split...
    add_split_zgc(this, m_child1, m_child2, ordinal);

    // Propogate parent ZoneGridConnectivities to appropriate children.
    // Split if needed...
    propogate_zgc(this, m_child1, ordinal, rank);
    propogate_zgc(this, m_child2, ordinal, rank);

    return std::make_pair(m_child1, m_child2);
  }

  // If a zgc points to a donor zone which was split (has non-null children),
  // then create two zgc that point to each child.  Update range and donor_range
  void StructuredZoneData::resolve_zgc_split_donor(std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    bool did_split = false;
    do {
      did_split = false;
      std::vector<Ioss::ZoneConnectivity> new_zgc;
      for (auto zgc : m_zoneConnectivity) {
        auto &donor_zone = zones[zgc.m_donorZone - 1];
        if (!donor_zone->is_active()) {
          did_split    = true;
          bool overlap = false;
          auto c1_zgc(zgc);
          c1_zgc.m_donorZone = donor_zone->m_child1->m_zone;
          if (zgc_donor_overlaps(donor_zone->m_child1, c1_zgc)) {
            overlap = true;
            zgc_subset_donor_ranges(donor_zone->m_child1, c1_zgc);
            if (c1_zgc.get_shared_node_count() > 2) {
              new_zgc.push_back(c1_zgc);
            }
          }
          auto c2_zgc(zgc);
          c2_zgc.m_donorZone = donor_zone->m_child2->m_zone;
          if (zgc_donor_overlaps(donor_zone->m_child2, c2_zgc)) {
            overlap = true;
            zgc_subset_donor_ranges(donor_zone->m_child2, c2_zgc);
            if (c2_zgc.get_shared_node_count() > 2) {
              new_zgc.push_back(c2_zgc);
            }
          }
          if (!overlap) {
            // Need to add at least one copy of this zgc even if no overlap
            // so can maintain the original (un parallel decomposed) ranges
            // for use in output...
            c1_zgc.m_ownerRangeBeg = c1_zgc.m_ownerRangeEnd = {{0, 0, 0}};
            c1_zgc.m_donorRangeBeg = c1_zgc.m_donorRangeEnd = {{0, 0, 0}};
            new_zgc.push_back(c1_zgc);
          }
        }
        else {
          new_zgc.push_back(zgc);
        }
      }
      if (did_split) {
        std::swap(m_zoneConnectivity, new_zgc);
      }
    } while (did_split);
  }

  void StructuredZoneData::update_zgc_processor(std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    for (auto &zgc : m_zoneConnectivity) {
      auto &donor_zone     = zones[zgc.m_donorZone - 1];
      zgc.m_donorProcessor = donor_zone->m_proc;
      auto &owner_zone     = zones[zgc.m_ownerZone - 1];
      zgc.m_ownerProcessor = owner_zone->m_proc;
    }
  }
} // namespace Iocgns
