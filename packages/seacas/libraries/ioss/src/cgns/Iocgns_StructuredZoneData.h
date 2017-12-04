/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef IOCGNS_STRUCTUREDZONEDATA_H
#define IOCGNS_STRUCTUREDZONEDATA_H

#include <Ioss_CodeTypes.h>
#include <Ioss_StructuredBlock.h>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace Iocgns {
  class StructuredZoneData
  {
  public:
    StructuredZoneData()
        : m_ordinal{{0, 0, 0}}, m_offset{{0, 0, 0}}, m_preferentialOrdinal(-1), m_zone(0),
          m_adam(nullptr), m_parent(nullptr), m_proc(-1), m_splitOrdinal(0), m_child1(nullptr),
          m_child2(nullptr), m_sibling(nullptr)
    {
    }

    StructuredZoneData(std::string name, int zone, int ni, int nj, int nk)
        : m_name(std::move(name)), m_ordinal{{ni, nj, nk}}, m_offset{{0, 0, 0}},
          m_preferentialOrdinal(-1), m_zone(zone), m_adam(nullptr), m_parent(nullptr), m_proc(-1),
          m_splitOrdinal(0), m_child1(nullptr), m_child2(nullptr), m_sibling(nullptr)
    {
    }

    std::string m_name;
    Ioss::IJK_t m_ordinal;

    // Offset of this block relative to its
    // adam block. ijk_adam = ijk_me + m_offset[ijk];
    Ioss::IJK_t m_offset;

    // If value is 0, 1, or 2, then do not split along that ordinal
    int m_preferentialOrdinal;

    int m_zone;

    // The zone in the undecomposed model that this zone is a
    // descendant of.  If not decomposed, then m_zone == m_adam;
    StructuredZoneData *m_adam;

    // If this zone was the result of splitting another zone, then
    // what is the zone number of that zone.  Zones are kept in a
    // vector and the zone number is its position in that vector+1
    // to make it 1-based and match numbering on file.
    StructuredZoneData *m_parent;

    int m_proc; // The processor this block might be run on...

    // Which ordinal of the parent was split to generate this zone and its sibling.
    int m_splitOrdinal;

    // The two zones that were split off from this zone.
    // Might be reasonable to do a 3-way or n-way split, but for now
    // just do a 2-way.
    StructuredZoneData *m_child1;
    StructuredZoneData *m_child2;

    StructuredZoneData *m_sibling;

    std::vector<Ioss::ZoneConnectivity> m_zoneConnectivity;

    // ========================================================================
    bool is_active() const
    {
      // Zone is active if it hasn't been split.
      return m_child1 == nullptr && m_child2 == nullptr;
    }

    // ========================================================================
    // Assume the "work" or computational effort required for a
    // block is proportional to the number of cells.
    size_t work() const { return m_ordinal[0] * m_ordinal[1] * m_ordinal[2]; }

    std::pair<StructuredZoneData *, StructuredZoneData *> split(int zone_id, double ratio = 0.5,
                                                                int rank = 0);
    void resolve_zgc_split_donor(std::vector<Iocgns::StructuredZoneData *> &zones);
    void update_zgc_processor(std::vector<Iocgns::StructuredZoneData *> &zones);
  };
} // namespace Iocgns

#endif
