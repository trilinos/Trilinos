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

#include <cgns/Iocgns_Defines.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_TerminalColor.h>
#include <Ioss_Utils.h>
#include <cgns/Iocgns_DecompositionData.h>
#include <cgns/Iocgns_Utils.h>
#include <tokenize.h>

#include <cgnsconfig.h>
#include <pcgnslib.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <numeric>

namespace {
  int rank = 0;
#define OUTPUT                                                                                     \
  if (rank == 0)                                                                                   \
  std::cerr

  // ZOLTAN Callback functions...

#if !defined(NO_ZOLTAN_SUPPORT)
  int zoltan_num_dim(void *data, int *ierr)
  {
    // Return dimensionality of coordinate data.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->spatial_dimension();
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    // Return number of objects (element count) on this processor...
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->decomp_elem_count();
  }

  void zoltan_obj_list(void *data, int ngid_ent, int nlid_ent, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR lids, int wdim, float *wgts, int *ierr)
  {
    // Return list of object IDs, both local and global.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    // At the time this is called, we don't have much information
    // These routines are the ones that are developing that
    // information...
    size_t element_count  = zdata->decomp_elem_count();
    size_t element_offset = zdata->decomp_elem_offset();

    *ierr = ZOLTAN_OK;

    if (lids) {
      std::iota(lids, lids + element_count, 0);
    }

    if (wdim) {
      std::fill(wgts, wgts + element_count, 1.0);
    }

    if (ngid_ent == 1) {
      std::iota(gids, gids + element_count, element_offset);
    }
    else if (ngid_ent == 2) {
      int64_t *global_ids = (int64_t *)gids;
      std::iota(global_ids, global_ids + element_count, element_offset);
    }
    else {
      *ierr = ZOLTAN_FATAL;
    }
    return;
  }

  void zoltan_geom(void *data, int ngid_ent, int nlid_ent, int nobj, ZOLTAN_ID_PTR gids,
                   ZOLTAN_ID_PTR lids, int ndim, double *geom, int *ierr)
  {
    // Return coordinates for objects.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    std::copy(zdata->centroids().begin(), zdata->centroids().end(), &geom[0]);

    *ierr = ZOLTAN_OK;
    return;
  }
#endif

  // These are used for structured parallel decomposition...
  void create_zone_data(int cgnsSerFilePtr, int cgnsFilePtr,
                        std::vector<Iocgns::StructuredZoneData *> &zones, MPI_Comm comm)
  {
    Ioss::ParallelUtils par_util(comm);
    int                 myProcessor = par_util.parallel_rank(); // To make error macro work...
    int                 base        = 1;
    int                 num_zones   = 0;

    CGCHECK(cg_nzones(cgnsFilePtr, base, &num_zones));

    std::map<std::string, int> zone_name_map;

    for (cgsize_t zone = 1; zone <= num_zones; zone++) {
      cgsize_t size[9];
      char     zone_name[CGNS_MAX_NAME_LENGTH + 1];
      CGCHECK(cg_zone_read(cgnsFilePtr, base, zone, zone_name, size));
      zone_name_map[zone_name] = zone;

      SMART_ASSERT(size[0] - 1 == size[3])(size[0])(size[3]);
      SMART_ASSERT(size[1] - 1 == size[4])(size[1])(size[4]);
      SMART_ASSERT(size[2] - 1 == size[5])(size[2])(size[5]);

      assert(size[6] == 0);
      assert(size[7] == 0);
      assert(size[8] == 0);

      auto *zone_data = new Iocgns::StructuredZoneData(zone_name, zone, size[3], size[4], size[5]);
      zones.push_back(zone_data);

      // Handle zone-grid-connectivity...
      if (rank == 0) {
        int nconn = 0;
        CGCHECK(cg_n1to1(cgnsSerFilePtr, base, zone, &nconn));
        for (int i = 0; i < nconn; i++) {
          char                    connectname[CGNS_MAX_NAME_LENGTH + 1];
          char                    donorname[CGNS_MAX_NAME_LENGTH + 1];
          std::array<cgsize_t, 6> range;
          std::array<cgsize_t, 6> donor_range;
          Ioss::IJK_t             transform;

          CGCHECK(cg_1to1_read(cgnsSerFilePtr, base, zone, i + 1, connectname, donorname,
                               range.data(), donor_range.data(), transform.data()));

          // Get number of nodes shared with other "previous" zones...
          // A "previous" zone will have a lower zone number this this zone...
          int  donor_zone = -1;
          auto donor_iter = zone_name_map.find(donorname);
          if (donor_iter != zone_name_map.end()) {
            donor_zone = (*donor_iter).second;
          }
          Ioss::IJK_t range_beg{{(int)range[0], (int)range[1], (int)range[2]}};
          Ioss::IJK_t range_end{{(int)range[3], (int)range[4], (int)range[5]}};
          Ioss::IJK_t donor_beg{{(int)donor_range[0], (int)donor_range[1], (int)donor_range[2]}};
          Ioss::IJK_t donor_end{{(int)donor_range[3], (int)donor_range[4], (int)donor_range[5]}};

#if IOSS_DEBUG_OUTPUT
          OUTPUT << "Adding zgc " << connectname << " to " << zone_name << " donor: " << donorname
                 << "\n";
#endif
          zone_data->m_zoneConnectivity.emplace_back(connectname, zone, donorname, donor_zone,
                                                     transform, range_beg, range_end, donor_beg,
                                                     donor_end);
        }
      }
    }

    // If parallel, pack the data on rank 0 and broadcast to all other processors...
#ifdef SEACAS_HAVE_MPI

    if (par_util.parallel_size() > 1) {
      std::vector<int> zgc_size(zones.size());
      // Let each processor know how many zgc each of its zones should have...
      if (rank == 0) {
        for (size_t i = 0; i < zones.size(); i++) {
          zgc_size[i] = (int)zones[i]->m_zoneConnectivity.size();
        }
      }
      MPI_Bcast(zgc_size.data(), (int)zgc_size.size(), MPI_INT, 0, comm);
      int count = std::accumulate(zgc_size.begin(), zgc_size.end(), (int)0);

      // Pack the zgc for all zones on rank=0 and send to all other ranks for unpacking.
      const int         BYTE_PER_NAME = CGNS_MAX_NAME_LENGTH + 1;
      const int         INT_PER_ZGC   = 17;
      std::vector<char> zgc_name(count * 2 * BYTE_PER_NAME);
      std::vector<int>  zgc_data(count * INT_PER_ZGC);

      if (rank == 0) {
        // Pack the data...
        int off_name = 0;
        int off_data = 0;
        int off_cnt  = 0;

        for (auto &zone : zones) {
          for (auto &z : zone->m_zoneConnectivity) {
            strncpy(&zgc_name[off_name], z.m_connectionName.c_str(), BYTE_PER_NAME);
            off_name += BYTE_PER_NAME;
            strncpy(&zgc_name[off_name], z.m_donorName.c_str(), BYTE_PER_NAME);
            off_name += BYTE_PER_NAME;

            off_cnt++;

            zgc_data[off_data++] = z.m_ownerZone;
            zgc_data[off_data++] = z.m_donorZone;

            zgc_data[off_data++] = z.m_ownerRangeBeg[0];
            zgc_data[off_data++] = z.m_ownerRangeBeg[1];
            zgc_data[off_data++] = z.m_ownerRangeBeg[2];
            zgc_data[off_data++] = z.m_ownerRangeEnd[0];
            zgc_data[off_data++] = z.m_ownerRangeEnd[1];
            zgc_data[off_data++] = z.m_ownerRangeEnd[2];

            zgc_data[off_data++] = z.m_donorRangeBeg[0];
            zgc_data[off_data++] = z.m_donorRangeBeg[1];
            zgc_data[off_data++] = z.m_donorRangeBeg[2];
            zgc_data[off_data++] = z.m_donorRangeEnd[0];
            zgc_data[off_data++] = z.m_donorRangeEnd[1];
            zgc_data[off_data++] = z.m_donorRangeEnd[2];

            zgc_data[off_data++] = z.m_transform[0];
            zgc_data[off_data++] = z.m_transform[1];
            zgc_data[off_data++] = z.m_transform[2];
          }
        }
        assert(off_cnt == count);
        assert(count == 0 || (off_data % count == 0));
        assert(count == 0 || (off_data / count == INT_PER_ZGC));
        assert(count == 0 || (off_name % count == 0 && off_name / count / 2 == BYTE_PER_NAME));
      }

      MPI_Bcast(zgc_name.data(), (int)zgc_name.size(), MPI_CHAR, 0, comm);
      MPI_Bcast(zgc_data.data(), (int)zgc_data.size(), MPI_INT, 0, comm);

      if (rank != 0) {
        // Unpack the data...
        int off_name = 0;
        int off_data = 0;
        int off_cnt  = 0;

        for (size_t i = 0; i < zones.size(); i++) {
          auto zgc_cnt = zgc_size[i];
          auto zone    = zones[i];
          for (int j = 0; j < zgc_cnt; j++) {
            off_cnt++;
            std::string name{&zgc_name[off_name]};
            off_name += BYTE_PER_NAME;
            std::string donor_name{&zgc_name[off_name]};
            off_name += BYTE_PER_NAME;

            int         zone_id  = zgc_data[off_data++];
            int         donor_id = zgc_data[off_data++];
            Ioss::IJK_t range_beg{
                {zgc_data[off_data++], zgc_data[off_data++], zgc_data[off_data++]}};
            Ioss::IJK_t range_end{
                {zgc_data[off_data++], zgc_data[off_data++], zgc_data[off_data++]}};
            Ioss::IJK_t donor_beg{
                {zgc_data[off_data++], zgc_data[off_data++], zgc_data[off_data++]}};
            Ioss::IJK_t donor_end{
                {zgc_data[off_data++], zgc_data[off_data++], zgc_data[off_data++]}};
            Ioss::IJK_t transform{
                {zgc_data[off_data++], zgc_data[off_data++], zgc_data[off_data++]}};
            zone->m_zoneConnectivity.emplace_back(name, zone_id, donor_name, donor_id, transform,
                                                  range_beg, range_end, donor_beg, donor_end);
          }
          assert((int)zone->m_zoneConnectivity.size() == zgc_cnt);
        }
        assert(off_cnt == count);
        assert(count == 0 || (off_data % count == 0));
        assert(count == 0 || (off_data / count == INT_PER_ZGC));
        assert(count == 0 || (off_name % count == 0 && off_name / count / 2 == BYTE_PER_NAME));
      }
    }
#endif

    // If there are any Structured blocks, need to iterate them and their 1-to-1 connections
    // and update the donor_zone id for zones that had not yet been processed at the time of
    // definition...
    for (auto &zone : zones) {
      for (auto &conn : zone->m_zoneConnectivity) {
        if (conn.m_donorZone < 0) {
          auto donor_iter = zone_name_map.find(conn.m_donorName);
          assert(donor_iter != zone_name_map.end());
          conn.m_donorZone = (*donor_iter).second;
        }
      }
    }
  }

  void set_line_decomposition(int cgnsFilePtr, const std::string &line_decomposition,
                              std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    // The "line_decomposition" string is a list of 0 or more BC
    // (Family) names.  For all structured zones which this BC
    // touches, the ordinal of the face (i,j,k) will be set such that
    // a parallel decomposition will not split the zone along this
    // ordinal.  For example, if the BC "wall1" has the definition
    // [1->1, 1->5, 1->8], then it is on the constant 'i' face of the
    // zone and therefore, the zone will *not* be split along the 'i'
    // ordinal.

    // Get names of all valid 'bcs' on the mesh
    int base         = 1;
    int num_families = 0;
    CGCHECKNP(cg_nfamilies(cgnsFilePtr, base, &num_families));

    std::vector<std::string> families;
    families.reserve(num_families);
    for (int family = 1; family <= num_families; family++) {
      char name[CGNS_MAX_NAME_LENGTH + 1];
      int  num_bc  = 0;
      int  num_geo = 0;
      CGCHECKNP(cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo));
      if (num_bc > 0) {
        Ioss::Utils::fixup_name(name);
        families.push_back(name);
      }
    }

    // Slit into fields using the commas as delimiters
    auto bcs = Ioss::tokenize(line_decomposition, ",");
    for (auto &bc : bcs) {
      Ioss::Utils::fixup_name(bc);
      if (std::find(families.begin(), families.end(), bc) == families.end()) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: The family/bc name '" << bc
               << "' specified as a line decomposition surface does not exist on this CGNS file.\n";
        errmsg << "             Valid names are: ";
        for (const auto &fam : families) {
          errmsg << "'" << fam << "', ";
        }
        IOSS_ERROR(errmsg);
      }
    }

    for (auto zone : zones) {
      // Read BCs applied to this zone and see if they match any of
      // the BCs in 'bcs' list.  If so, determine the face the BC is
      // applied to and set the m_lineOrdinal to the ordinal
      // perpendicular to this face.
      int izone = zone->m_zone;
      int num_bcs;
      CGCHECKNP(cg_nbocos(cgnsFilePtr, base, izone, &num_bcs));

      for (int ibc = 0; ibc < num_bcs; ibc++) {
        char              boconame[CGNS_MAX_NAME_LENGTH + 1];
        CG_BCType_t       bocotype;
        CG_PointSetType_t ptset_type;
        cgsize_t          npnts;
        cgsize_t          NormalListSize;
        CG_DataType_t     NormalDataType;
        int               ndataset;

        // All we really want from this is 'boconame'
        CGCHECKNP(cg_boco_info(cgnsFilePtr, base, izone, ibc + 1, boconame, &bocotype, &ptset_type,
                               &npnts, nullptr, &NormalListSize, &NormalDataType, &ndataset));

        if (bocotype == CG_FamilySpecified) {
          // Need to get boconame from cg_famname_read
          CGCHECKNP(
              cg_goto(cgnsFilePtr, base, "Zone_t", izone, "ZoneBC_t", 1, "BC_t", ibc + 1, "end"));
          CGCHECKNP(cg_famname_read(boconame));
        }

        Ioss::Utils::fixup_name(boconame);
        if (std::find(bcs.begin(), bcs.end(), boconame) != bcs.end()) {
          cgsize_t range[6];
          CGCHECKNP(cg_boco_read(cgnsFilePtr, base, izone, ibc + 1, range, nullptr));

          // There are some BC that are applied on an edge or a vertex;
          // Don't want those, so filter them out at this time...
          bool i = range[0] == range[3];
          bool j = range[1] == range[4];
          bool k = range[2] == range[5];

          int sum = (i ? 1 : 0) + (j ? 1 : 0) + (k ? 1 : 0);
          // Only set m_lineOrdinal if only a single ordinal selected.
          if (sum == 1) {
            int ordinal = -1;
            if (i) {
              ordinal = 0;
            }
            else if (j) {
              ordinal = 1;
            }
            else if (k) {
              ordinal = 2;
            }
            if (zone->m_lineOrdinal == -1) {
              zone->m_lineOrdinal = ordinal;
#if IOSS_DEBUG_OUTPUT
              OUTPUT << "Setting line ordinal to " << zone->m_lineOrdinal << " on " << zone->m_name
                     << " for surface: " << boconame << "\n";
#endif
            }
            else if (zone->m_lineOrdinal != ordinal && rank == 0) {
              IOSS_WARNING
                  << "CGNS: Zone " << izone << " named " << zone->m_name
                  << " has multiple line decomposition ordinal specifications. Both ordinal "
                  << ordinal << " and " << zone->m_lineOrdinal << " have been specified.  Keeping "
                  << zone->m_lineOrdinal << "\n";
            }
          }
        }
      }
    }
  }

} // namespace

namespace Iocgns {
  template DecompositionData<int>::DecompositionData(const Ioss::PropertyManager &props,
                                                     MPI_Comm                     communicator);
  template DecompositionData<int64_t>::DecompositionData(const Ioss::PropertyManager &props,
                                                         MPI_Comm                     communicator);

  template <typename INT>
  DecompositionData<INT>::DecompositionData(const Ioss::PropertyManager &props,
                                            MPI_Comm                     communicator)
      : DecompositionDataBase(communicator), m_decomposition(props, communicator)
  {
    rank = m_decomposition.m_processor;

    if (props.exists("LOAD_BALANCE_THRESHOLD")) {
      if (props.get("LOAD_BALANCE_THRESHOLD").get_type() == Ioss::Property::STRING) {
        std::string lb_thresh  = props.get("LOAD_BALANCE_THRESHOLD").get_string();
        m_loadBalanceThreshold = std::strtod(lb_thresh.c_str(), nullptr);
      }
      else if (props.get("LOAD_BALANCE_THRESHOLD").get_type() == Ioss::Property::REAL) {
        m_loadBalanceThreshold = props.get("LOAD_BALANCE_THRESHOLD").get_real();
      }
    }
    if (props.exists("LINE_DECOMPOSITION")) {
      m_lineDecomposition = props.get("LINE_DECOMPOSITION").get_string();
    }
  }

  template <typename INT>
  void DecompositionData<INT>::decompose_model(int serFilePtr, int filePtr,
                                               CG_ZoneType_t common_zone_type)
  {
    if (common_zone_type == CG_Unstructured) {
      decompose_unstructured(filePtr);
    }
    else if (common_zone_type == CG_Structured) {
      decompose_structured(serFilePtr, filePtr);
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: The common zone type is not of type Unstructured or Structured "
                "which are the only types currently supported";
      IOSS_ERROR(errmsg);
    }
  }

  template <typename INT>
  void DecompositionData<INT>::decompose_structured(int serFilePtr, int filePtr)
  {
    m_decomposition.show_progress(__func__);
    create_zone_data(serFilePtr, filePtr, m_structuredZones, m_decomposition.m_comm);
    if (m_structuredZones.empty()) {
      return;
    }

    // Determine whether user has specified "line decompositions" for any of the zones.
    // The line decomposition is an ordinal which will not be split during the
    // decomposition.
    if (!m_lineDecomposition.empty()) {
      set_line_decomposition(filePtr, m_lineDecomposition, m_structuredZones);
    }

    size_t work = 0;
    for (const auto &z : m_structuredZones) {
      work += z->work();
      assert(z->is_active());
    }

    size_t px        = 0;
    size_t num_split = 0;
    double avg_work  = (double)work / m_decomposition.m_processorCount;

#if IOSS_DEBUG_OUTPUT
    auto num_active = m_structuredZones.size();
    OUTPUT << "Decomposing structured mesh with " << num_active << " zones for "
           << m_decomposition.m_processorCount << " processors.\nAverage workload is " << avg_work
           << ", Load Balance Threshold is " << m_loadBalanceThreshold << ", Work range "
           << avg_work / m_loadBalanceThreshold << " to " << avg_work * m_loadBalanceThreshold
           << "\n";
#endif

    if (avg_work < 1.0) {
      OUTPUT << "ERROR: Model size too small to distribute over "
             << m_decomposition.m_processorCount << " processors.\n";
      std::exit(EXIT_FAILURE);
    }

#if IOSS_DEBUG_OUTPUT
    OUTPUT << "========================================================================\n";
    OUTPUT << "Pre-Splitting:\n";
#endif
    // Split all blocks where block->work() > avg_work * m_loadBalanceThreshold
    size_t new_zone_id = Utils::pre_split(m_structuredZones, avg_work, m_loadBalanceThreshold, rank,
                                          m_decomposition.m_processorCount);

    // At this point, there should be no zone with block->work() > avg_work * m_loadBalanceThreshold
#if IOSS_DEBUG_OUTPUT
    OUTPUT << "========================================================================\n";
#endif
    do {
      std::vector<size_t> work_vector(m_decomposition.m_processorCount);
      Utils::assign_zones_to_procs(m_structuredZones, work_vector);

      // Calculate workload ratio for each processor...
      px = 0; // Number of processors where workload ratio exceeds threshold.
      std::vector<bool> exceeds(m_decomposition.m_processorCount);
      for (size_t i = 0; i < work_vector.size(); i++) {
        double workload_ratio = double(work_vector[i]) / double(avg_work);
#if IOSS_DEBUG_OUTPUT
        OUTPUT << "\nProcessor " << i << " work: " << work_vector[i]
               << ", workload ratio: " << workload_ratio;
#endif
        if (workload_ratio > m_loadBalanceThreshold) {
          exceeds[i] = true;
          px++;
        }
      }
#if IOSS_DEBUG_OUTPUT
      OUTPUT << "\n\nWorkload threshold exceeded on " << px << " processors.\n";
#endif
      bool single_zone = m_structuredZones.size() == 1;
      if (single_zone) {
        auto active = std::count_if(m_structuredZones.begin(), m_structuredZones.end(),
                                    [](Iocgns::StructuredZoneData *a) { return a->is_active(); });
        if (active >= m_decomposition.m_processorCount) {
          px = 0;
        }
      }
      num_split = 0;
      if (px > 0) {
        auto zone_new(m_structuredZones);
        for (auto zone : m_structuredZones) {
          if (zone->is_active() && exceeds[zone->m_proc]) {
            // Since 'zones' is sorted from most work to least,
            // we just iterate zones and check whether the zone
            // is on a proc where the threshold was exceeded.
            // if so, split the block and set exceeds[proc] to false;
            // Exit the loop when num_split >= px.
            auto children =
                zone->split(new_zone_id, zone->work() / 2.0, m_loadBalanceThreshold, rank);
            if (children.first != nullptr && children.second != nullptr) {
              zone_new.push_back(children.first);
              zone_new.push_back(children.second);

              new_zone_id += 2;
              exceeds[zone->m_proc] = false;
              num_split++;
              if (num_split >= px) {
                break;
              }
            }
          }
        }
        std::swap(zone_new, m_structuredZones);
      }
#if IOSS_DEBUG_OUTPUT
      auto active = std::count_if(m_structuredZones.begin(), m_structuredZones.end(),
                                  [](Iocgns::StructuredZoneData *a) { return a->is_active(); });
      OUTPUT << "Number of active zones = " << active << ", average work = " << avg_work << "\n";
      OUTPUT << "========================================================================\n";
#endif
    } while (px > 0 && num_split > 0);

    std::sort(m_structuredZones.begin(), m_structuredZones.end(),
              [](Iocgns::StructuredZoneData *a, Iocgns::StructuredZoneData *b) {
                return a->m_zone < b->m_zone;
              });

    for (auto zone : m_structuredZones) {
      if (zone->is_active()) {
        zone->resolve_zgc_split_donor(m_structuredZones);
      }
    }

    // Update and Output the processor assignments
    for (auto &zone : m_structuredZones) {
      if (zone->is_active()) {
        zone->update_zgc_processor(m_structuredZones);
#if IOSS_DEBUG_OUTPUT
        auto zone_node_count =
            (zone->m_ordinal[0] + 1) * (zone->m_ordinal[1] + 1) * (zone->m_ordinal[2] + 1);
        OUTPUT << "Zone " << zone->m_name << "(" << zone->m_zone << ") assigned to processor "
               << zone->m_proc << ", Adam zone = " << zone->m_adam->m_zone
               << ", Cells = " << zone->work() << ", Nodes = " << zone_node_count << "\n";
        auto zgcs = zone->m_zoneConnectivity;
        for (auto &zgc : zgcs) {
          OUTPUT << zgc << "\n";
        }
#endif
      }
    }

    // Output the processor assignments in form similar to 'split' file
    if (rank == 0) {
      int z = 1;
      std::cerr
          << "     n    proc  parent    imin    imax    jmin    jmax    kmin     kmax     work\n";
      auto tmp_zone(m_structuredZones);
      std::sort(tmp_zone.begin(), tmp_zone.end(),
                [](Iocgns::StructuredZoneData *a, Iocgns::StructuredZoneData *b) {
                  return a->m_proc < b->m_proc;
                });

      for (auto &zone : tmp_zone) {
        if (zone->is_active()) {
          std::cerr << std::setw(6) << z++ << std::setw(8) << zone->m_proc << std::setw(8)
                    << zone->m_adam->m_zone << std::setw(8) << zone->m_offset[0] + 1 << std::setw(8)
                    << zone->m_ordinal[0] + zone->m_offset[0] + 1 << std::setw(8)
                    << zone->m_offset[1] + 1 << std::setw(8)
                    << zone->m_ordinal[1] + zone->m_offset[1] + 1 << std::setw(8)
                    << zone->m_offset[2] + 1 << std::setw(8)
                    << zone->m_ordinal[2] + zone->m_offset[2] + 1 << std::setw(8) << zone->work()
                    << "\n";
        }
      }
    }

    for (auto &zone : m_structuredZones) {
      if (!zone->is_active()) {
        zone->m_proc = -1;
      }
    }

#if IOSS_DEBUG_OUTPUT
    MPI_Barrier(m_decomposition.m_comm);
    OUTPUT << Ioss::trmclr::green << "Returning from decomposition\n" << Ioss::trmclr::normal;
#endif
  }

  template <typename INT> void DecompositionData<INT>::decompose_unstructured(int filePtr)
  {
    m_decomposition.show_progress(__func__);

    // Initial decomposition is linear where processor #p contains
    // elements from (#p * #element/#proc) to (#p+1 * #element/#proc)

    // ========================================================================
    // Get the number of zones (element blocks) in the mesh...
    int num_zones = 0;
    int base      = 1; // Only single base supported so far.

    {
      int  cell_dimension = 0;
      int  phys_dimension = 0;
      char base_name[CGNS_MAX_NAME_LENGTH + 1];
      CGCHECK2(cg_base_read(filePtr, base, base_name, &cell_dimension, &phys_dimension));
      m_decomposition.m_spatialDimension = phys_dimension;
    }

    CGCHECK2(cg_nzones(filePtr, base, &num_zones));
    m_zones.resize(num_zones + 1); // Use 1-based zones.

    size_t global_cell_node_count = 0;
    size_t global_element_count   = 0;
    for (int zone = 1; zone <= num_zones; zone++) {
      // All zones are "Unstructured" since this was checked prior to
      // calling this function...
      cgsize_t size[3];
      char     zone_name[CGNS_MAX_NAME_LENGTH + 1];
      CGCHECK2(cg_zone_read(filePtr, base, zone, zone_name, size));

      INT total_block_nodes = size[0];
      INT total_block_elem  = size[1];

      m_zones[zone].m_nodeCount     = total_block_nodes;
      m_zones[zone].m_nodeOffset    = global_cell_node_count;
      m_zones[zone].m_name          = zone_name;
      m_zones[zone].m_elementOffset = global_element_count;
      global_cell_node_count += total_block_nodes;
      global_element_count += total_block_elem;
    }

    // Generate element_dist/node_dist --  size m_decomposition.m_processorCount + 1
    // processor p contains all elements/nodes from X_dist[p] .. X_dist[p+1]
    m_decomposition.generate_entity_distributions(global_cell_node_count, global_element_count);

    generate_adjacency_list(filePtr, m_decomposition);

    // Get min and max node used on this processor...
    auto min_max =
        std::minmax_element(m_decomposition.m_adjacency.begin(), m_decomposition.m_adjacency.end());
    INT min_node = *(min_max.first);
    INT max_node = *(min_max.second);
    generate_zone_shared_nodes(filePtr, min_node, max_node);

    // Now iterate adjacency list and update any "zone_shared_node" nodes
    // with their "sharee"
    if (!m_zoneSharedMap.empty()) {
      for (auto &node : m_decomposition.m_adjacency) {
        auto alias = m_zoneSharedMap.find(node);
        if (alias != m_zoneSharedMap.end()) {
          node = (*alias).second;
        }
      }
    }

#if IOSS_DEBUG_OUTPUT
    OUTPUT << "Processor " << m_decomposition.m_processor << " has " << decomp_elem_count()
           << " elements; offset = " << decomp_elem_offset() << "\n";
    OUTPUT << "Processor " << m_decomposition.m_processor << " has " << decomp_node_count()
           << " nodes; offset = " << decomp_node_offset() << ".\n";
#endif

    if (m_decomposition.needs_centroids()) {
      // Get my coordinate data using direct cgns calls
      std::vector<double> x(decomp_node_count());
      std::vector<double> y;
      std::vector<double> z;

      get_file_node_coordinates(filePtr, 0, TOPTR(x));
      if (m_decomposition.m_spatialDimension > 1) {
        y.resize(decomp_node_count());
        get_file_node_coordinates(filePtr, 1, TOPTR(y));
      }
      if (m_decomposition.m_spatialDimension > 2) {
        z.resize(decomp_node_count());
        get_file_node_coordinates(filePtr, 2, TOPTR(z));
      }

      m_decomposition.calculate_element_centroids(x, y, z);
    }

#if !defined(NO_ZOLTAN_SUPPORT)
    float version = 0.0;
    Zoltan_Initialize(0, nullptr, &version);

    Zoltan zz(m_decomposition.m_comm);

    // Register Zoltan Callback functions...
    zz.Set_Num_Obj_Fn(zoltan_num_obj, this);
    zz.Set_Obj_List_Fn(zoltan_obj_list, this);
    zz.Set_Num_Geom_Fn(zoltan_num_dim, this);
    zz.Set_Geom_Multi_Fn(zoltan_geom, this);
#endif

    m_decomposition.decompose_model(
#if !defined(NO_ZOLTAN_SUPPORT)
        zz,
#endif
        m_elementBlocks);

    if (!m_sideSets.empty()) {
      // Create elemGTL map which is used for sidesets (also element sets)
      build_global_to_local_elem_map();
    }

    get_sideset_data(filePtr);

    // Have all the decomposition data needed
    // Can now populate the Ioss metadata...
  }

  template <typename INT>
  void DecompositionData<INT>::generate_zone_shared_nodes(int filePtr, INT min_node, INT max_node)
  {
    // Begin of Zone-Shared node information

    // Modify adjacency list based on shared nodes between zones...
    // Need the map from "global" to "global-shared"
    // * This is not necessarily nodes only on my processor since connectivity can include
    //   nodes other than those I own.
    // * Potentially large number of shared nodes; practically small(?)

    // * Maintain hash map from old id to new (if any)
    // * TODO: Make more scalable

    int base = 1; // Only single base supported so far.

    // Donor zone is always lower numbered, so zone 1 has no donor zone. Start at zone 2.
    for (cgsize_t zone = 2; zone < (cgsize_t)m_zones.size(); zone++) {

      // Determine number of "shared" nodes (shared with other zones)
      int nconn = 0;
      CGCHECK2(cg_nconns(filePtr, base, zone, &nconn));
      for (int i = 0; i < nconn; i++) {
        char                      connectname[CGNS_MAX_NAME_LENGTH + 1];
        CG_GridLocation_t         location;
        CG_GridConnectivityType_t connect_type;
        CG_PointSetType_t         ptset_type;
        cgsize_t                  npnts = 0;
        char                      donorname[CGNS_MAX_NAME_LENGTH + 1];
        CG_ZoneType_t             donor_zonetype;
        CG_PointSetType_t         donor_ptset_type;
        CG_DataType_t             donor_datatype;
        cgsize_t                  ndata_donor;

        CGCHECK2(cg_conn_info(filePtr, base, zone, i + 1, connectname, &location, &connect_type,
                              &ptset_type, &npnts, donorname, &donor_zonetype, &donor_ptset_type,
                              &donor_datatype, &ndata_donor));

        if (connect_type != CG_Abutting1to1 || ptset_type != CG_PointList ||
            donor_ptset_type != CG_PointListDonor) {
          std::ostringstream errmsg;
          errmsg << "ERROR: CGNS: Zone " << zone
                 << " adjacency data is not correct type. Require Abutting1to1 and PointList."
                 << connect_type << "\t" << ptset_type << "\t" << donor_ptset_type;
          IOSS_ERROR(errmsg);
        }

        // Verify data consistency...
        if (npnts != ndata_donor) {
          std::ostringstream errmsg;
          errmsg << "ERROR: CGNS: Zone " << zone << " point count (" << npnts
                 << ") does not match donor point count (" << ndata_donor << ").";
          IOSS_ERROR(errmsg);
        }

        // Get number of nodes shared with other "previous" zones...
        // A "previous" zone will have a lower zone number this this zone...
        std::string dz_name(donorname);
        int         dz = 1;
        for (; dz < zone; dz++) {
          if (m_zones[dz].m_name == dz_name)
            break;
        }

        if (dz != zone) {
#if IOSS_DEBUG_OUTPUT
          if (m_decomposition.m_processor == 0) {
            std::cerr << "Zone " << zone << " shares " << npnts << " nodes with " << donorname
                      << "\n";
          }
#endif
          // The 'ids' in 'points' and 'donors' will be zone-local 1-based.
          std::vector<cgsize_t> points(npnts);
          std::vector<cgsize_t> donors(npnts);

          CGCHECK2(cg_conn_read(filePtr, base, zone, i + 1, TOPTR(points), donor_datatype,
                                TOPTR(donors)));

          for (int j = 0; j < npnts; j++) {
            // Convert to 0-based global id by subtracting 1 and adding zone.m_nodeOffset
            cgsize_t point = points[j] - 1 + m_zones[zone].m_nodeOffset;
            cgsize_t donor = donors[j] - 1 + m_zones[dz].m_nodeOffset;

            // See if 'donor' is mapped to a different node already
            auto donor_map = m_zoneSharedMap.find(donor);
            if (donor_map != m_zoneSharedMap.end()) {
              donor = (*donor_map).second;
            }
            m_zoneSharedMap.insert({point, donor});
#if IOSS_DEBUG_OUTPUT
            if (m_decomposition.m_processor == 0) {
              std::cout << "Inserted " << point << " to " << donor << "\n";
            }
#endif
          }
        }
      }
    }
    // Filter m_zoneSharedMap down to nodes on this processor...
    // This processor contains global zone ids from `min_node` to `max_node`
    // global zone ids are the first entry in m_zoneShardedMap.
    for (auto it = m_zoneSharedMap.cbegin(); it != m_zoneSharedMap.cend(); /* no increment */) {
      if ((*it).first < min_node || (*it).first > max_node) {
        it = m_zoneSharedMap.erase(it);
      }
      else {
        ++it;
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::generate_adjacency_list(int                       filePtr,
                                                       Ioss::Decomposition<INT> &decomposition)
  {
    int base = 1; // Only single base supported so far.

    // Range of elements currently handled by this processor [)
    size_t p_start = decomp_elem_offset();
    size_t p_end   = p_start + decomp_elem_count();

    assert(sizeof(INT) == sizeof(cgsize_t));
    size_t sum    = 0; // Size of adjacency vector.
    size_t offset = 0;

    int num_zones        = 0;
    INT zone_node_offset = 0;

    CGCHECK2(cg_nzones(filePtr, base, &num_zones));
    for (int zone = 1; zone <= num_zones; zone++) {
      cgsize_t size[3];
      char     zone_name[CGNS_MAX_NAME_LENGTH + 1];
      CGCHECK2(cg_zone_read(filePtr, base, zone, zone_name, size));

      INT total_elements = size[1];
      // NOTE: A Zone will have a single set of nodes, but can have
      //       multiple sections each with their own element type...
      //       Keep treating sections as element blocks until we
      //       have handled 'size[1]' number of elements; the remaining
      //       sections are then the boundary faces (?)
      int num_sections = 0;
      CGCHECK2(cg_nsections(filePtr, base, zone, &num_sections));

      size_t last_blk_location = 0;
      for (int is = 1; is <= num_sections; is++) {
        char             section_name[CGNS_MAX_NAME_LENGTH + 1];
        CG_ElementType_t e_type;
        cgsize_t         el_start    = 0;
        cgsize_t         el_end      = 0;
        int              num_bndry   = 0;
        int              parent_flag = 0;

        // Get the type of elements in this section...
        CGCHECK2(cg_section_read(filePtr, base, zone, is, section_name, &e_type, &el_start, &el_end,
                                 &num_bndry, &parent_flag));

        INT num_entity = el_end - el_start + 1;

        if (parent_flag == 0 && total_elements > 0) {
          total_elements -= num_entity;

          // Range of elements in element block b [)
          size_t b_start = offset; // offset is index of first element in this block...
          offset += num_entity;
          size_t b_end = offset;

          int element_nodes;
          CGCHECK2(cg_npe(e_type, &element_nodes));

          if (b_start < p_end && p_start < b_end) {
            // Some of this blocks elements are on this processor...
            size_t overlap = std::min(b_end, p_end) - std::max(b_start, p_start);
            sum += overlap * element_nodes;
          }

          Ioss::BlockDecompositionData block;
          block.zone_          = zone;
          block.section_       = is;
          block.name_          = zone_name;
          block.topologyType   = Utils::map_cgns_to_topology_type(e_type);
          block.nodesPerEntity = element_nodes;
          block.fileCount      = num_entity;
          block.zoneNodeOffset = zone_node_offset;

          last_blk_location = m_elementBlocks.size();
          m_elementBlocks.push_back(block);
        }
        else {
          // This is a boundary-condition -- sideset (?)
          std::string ss_name(section_name);

          Ioss::SetDecompositionData sset;
          sset.zone_            = zone;
          sset.section_         = is;
          sset.name_            = ss_name;
          sset.fileCount        = num_entity;
          sset.topologyType     = Utils::map_cgns_to_topology_type(e_type);
          sset.parentBlockIndex = last_blk_location;
          m_sideSets.push_back(sset);
        }
      }
      zone_node_offset += size[0];
    }
    int block_count = (int)m_elementBlocks.size();

    // Get the global element block index list at this time also.
    // The global element at index 'I' (0-based) is on block B
    // if global_block_index[B] <= I && global_block_index[B+1] < I
    // allocate and TODO: Fill
    m_decomposition.m_fileBlockIndex.reserve(block_count + 1);
    for (auto block : m_elementBlocks) {
      m_decomposition.m_fileBlockIndex.push_back(block.file_count());
    }
    m_decomposition.m_fileBlockIndex.push_back(0);
    Ioss::Utils::generate_index(m_decomposition.m_fileBlockIndex);

    // Make sure 'sum' can fit in INT...
    INT tmp_sum = (INT)sum;
    if ((size_t)tmp_sum != sum) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The decomposition of this mesh requires 64-bit integers, but is being\n"
             << "       run with 32-bit integer code. Please rerun with the property "
                "INTEGER_SIZE_API\n"
             << "       set to 8. The details of how to do this vary with the code that is being "
                "run.\n"
             << "       Contact gdsjaar@sandia.gov for more details.\n";
      OUTPUT << errmsg.str();
      exit(EXIT_FAILURE);
    }

    // Now, populate the vectors...
    decomposition.m_pointer.reserve(decomp_elem_count() + 1);
    decomposition.m_adjacency.reserve(sum);
    offset = 0;
    sum    = 0; // Size of adjacency vector.

    for (auto &block : m_elementBlocks) {
      // Range of elements in element block b [)
      size_t b_start = offset; // offset is index of first element in this block...
      offset += block.file_count();
      size_t b_end = b_start + block.file_count();

      ssize_t overlap      = std::min(b_end, p_end) - std::max(b_start, p_start);
      overlap              = std::max(overlap, (ssize_t)0);
      block.fileCount      = overlap;
      size_t element_nodes = block.nodesPerEntity;
      int    zone          = block.zone_;
      int    section       = block.section_;

      // Get the connectivity (raw) for this portion of elements...
      std::vector<cgsize_t> connectivity(overlap * element_nodes);
      INT                   blk_start = std::max(b_start, p_start) - b_start + 1;
      INT                   blk_end   = blk_start + overlap - 1;
      blk_start                       = blk_start < 0 ? 0 : blk_start;
      blk_end                         = blk_end < 0 ? 0 : blk_end;
#if IOSS_DEBUG_OUTPUT
      OUTPUT << "Processor " << m_decomposition.m_processor << " has " << overlap
             << " elements on element block " << block.name() << "\t(" << blk_start << " to "
             << blk_end << ")\n";
#endif
      block.fileSectionOffset = blk_start;
      CGCHECK2(cgp_elements_read_data(filePtr, base, zone, section, blk_start, blk_end,
                                      TOPTR(connectivity)));
      size_t el          = 0;
      INT    zone_offset = block.zoneNodeOffset;

      for (ssize_t elem = 0; elem < overlap; elem++) {
        decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
        for (size_t k = 0; k < element_nodes; k++) {
          INT node = connectivity[el++] - 1 + zone_offset; // 0-based node
          decomposition.m_adjacency.push_back(node);
        }
      }
      sum += overlap * element_nodes;
    }
    decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
  }

  template <typename INT> void DecompositionData<INT>::get_sideset_data(int filePtr)
  {
    int base = 1; // Only single base supported so far.

    // Get total length of sideset elemlists...
    size_t elemlist_size = 0;
    for (auto &sset : m_sideSets) {
      elemlist_size += sset.file_count();
    }

    // Calculate the max "buffer" size usable for storing sideset
    // elemlists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "decomp_node_count()/2*2" is to
    // equalize the decomp_node_count() among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = std::max(10000, (decomp_node_count() / 2) * 2 * 3 *sizeof(double) /
    //    sizeof(cgsize_t));

    bool subsetting = false; // elemlist_size > max_size;

    if (subsetting) {
      assert(1 == 0);
    }
    else {
      // Can handle reading all sideset elem lists on a single
      // processor simultaneously.
      std::vector<cgsize_t> elemlist(elemlist_size);

      size_t offset = 0;
      for (auto &sset : m_sideSets) {

        // TODO? Possibly rewrite using cgi_read_int_data so can skip reading element connectivity
        int                   nodes_per_face = 4; // FIXME: sb->topology()->number_nodes();
        std::vector<cgsize_t> elements(nodes_per_face *
                                       sset.file_count()); // Not needed, but can't skip

        // We get:
        // *  num_to_get parent elements,
        // *  num_to_get zeros (other parent element for face, but on boundary so 0)
        // *  num_to_get face_on_element
        // *  num_to_get zeros (face on other parent element)
        std::vector<cgsize_t> parent(4 * sset.file_count());

        CGCHECK2(cg_elements_read(filePtr, base, sset.zone(), sset.section(), TOPTR(elements),
                                  TOPTR(parent)));

        // Move from 'parent' to 'elementlist'
        size_t zone_element_id_offset = m_zones[sset.zone()].m_elementOffset;
        for (size_t i = 0; i < sset.file_count(); i++) {
          elemlist[offset++] = parent[i] + zone_element_id_offset;
        }
      }
      SMART_ASSERT(offset == elemlist_size)(offset)(elemlist_size);

      // Each processor now has a complete list of all elems in all
      // sidesets.
      // Determine which of these are owned by the current
      // processor...
      {
        offset = 0;
        for (auto &sset : m_sideSets) {
          size_t ss_beg = offset;
          size_t ss_end = ss_beg + sset.file_count();

          for (size_t n = ss_beg; n < ss_end; n++) {
            cgsize_t elem = elemlist[n];
            // See if elem owned by this processor...
            if (i_own_elem(elem)) {
              // Save elem in this processors elemlist for this set.
              // The saved data is this elems location in the global
              // elemlist for this set.
              sset.entitylist_map.push_back(n - offset);
            }
          }
          offset = ss_end;
        }
      }

      // Each processor knows how many of the sideset elems it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero elem count is
      // the "root" for this sideset.
      {
        std::vector<int> has_elems_local(m_sideSets.size());
        for (size_t i = 0; i < m_sideSets.size(); i++) {
          has_elems_local[i] = m_sideSets[i].entitylist_map.empty() ? 0 : 1;
        }

        std::vector<int> has_elems(m_sideSets.size() * m_decomposition.m_processorCount);
        MPI_Allgather(TOPTR(has_elems_local), has_elems_local.size(), MPI_INT, TOPTR(has_elems),
                      has_elems_local.size(), MPI_INT, m_decomposition.m_comm);

        for (size_t i = 0; i < m_sideSets.size(); i++) {
          m_sideSets[i].hasEntities.resize(m_decomposition.m_processorCount);
          m_sideSets[i].root_ = m_decomposition.m_processorCount;
          for (int p = 0; p < m_decomposition.m_processorCount; p++) {
            if (p < m_sideSets[i].root_ && has_elems[p * m_sideSets.size() + i] != 0) {
              m_sideSets[i].root_ = p;
            }
            m_sideSets[i].hasEntities[p] = has_elems[p * m_sideSets.size() + i];
          }
        }
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_file_node_coordinates(int filePtr, int direction,
                                                         double *data) const
  {
    int      base        = 1; // Only single base supported so far.
    cgsize_t beg         = 0;
    cgsize_t end         = 0;
    cgsize_t offset      = 0;
    cgsize_t node_count  = decomp_node_count();
    cgsize_t node_offset = decomp_node_offset();

    int num_zones = (int)m_zones.size() - 1;
    for (int zone = 1; zone <= num_zones; zone++) {
      end += m_zones[zone].m_nodeCount;

      cgsize_t start  = std::max(node_offset, beg);
      cgsize_t finish = std::min(end, node_offset + node_count);
      cgsize_t count  = (finish > start) ? finish - start : 0;

      // Now adjust start for 1-based node numbering and the start of this zone...
      start  = start - beg + 1;
      finish = finish - beg;
      if (count == 0) {
        start  = 0;
        finish = 0;
      }
#if IOSS_DEBUG_OUTPUT
      OUTPUT << m_decomposition.m_processor << ": reading " << count << " nodes from zone " << zone
             << " starting at " << start << " with an offset of " << offset << " ending at "
             << finish << "\n";
#endif
      double *coords = nullptr;
      if (count > 0) {
        coords = &data[offset];
      }
      CGCHECK2(cgp_coord_read_data(filePtr, base, zone, direction + 1, &start, &finish, coords));
      offset += count;
      beg = end;
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_node_coordinates(int filePtr, double *ioss_data,
                                                    const Ioss::Field &field) const
  {
    std::vector<double> tmp(decomp_node_count());
    if (field.get_name() == "mesh_model_coordinates_x") {
      get_file_node_coordinates(filePtr, 0, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_y") {
      get_file_node_coordinates(filePtr, 1, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_z") {
      get_file_node_coordinates(filePtr, 2, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates") {
      // Data required by upper classes store x0, y0, z0, ... xn,
      // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
      // ..., yn, z0, ..., zn so we have to allocate some scratch
      // memory to read in the data and then map into supplied
      // 'data'

      std::vector<double> ioss_tmp(ioss_node_count());

      // This implementation trades off extra communication for
      // reduced memory overhead.
      // * This method uses 'ioss_node_count' extra memory; 3
      // reads; and 3 communicate_node_data calls.
      //
      // * Other method uses 6*ioss_node_count extra memory; 3 reads;
      // and 1 communicate_node_data call.
      //
      for (int d = 0; d < m_decomposition.m_spatialDimension; d++) {
        get_file_node_coordinates(filePtr, d, TOPTR(tmp));
        communicate_node_data(TOPTR(tmp), TOPTR(ioss_tmp), 1);

        size_t index = d;
        for (size_t i = 0; i < ioss_node_count(); i++) {
          ioss_data[index] = ioss_tmp[i];
          index += m_decomposition.m_spatialDimension;
        }
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_node_field(int filePtr, int step, int field_offset,
                                              double *ioss_data) const
  {
    std::vector<double> tmp(decomp_node_count());

    int      base        = 1; // Only single base supported so far.
    cgsize_t beg         = 0;
    cgsize_t end         = 0;
    cgsize_t offset      = 0;
    cgsize_t node_count  = decomp_node_count();
    cgsize_t node_offset = decomp_node_offset();

    int num_zones = (int)m_zones.size() - 1;
    for (int zone = 1; zone <= num_zones; zone++) {
      end += m_zones[zone].m_nodeCount;

      int solution_index = Utils::find_solution_index(filePtr, base, zone, step, CG_Vertex);

      cgsize_t start  = std::max(node_offset, beg);
      cgsize_t finish = std::min(end, node_offset + node_count);
      cgsize_t count  = (finish > start) ? finish - start : 0;

      // Now adjust start for 1-based node numbering and the start of this zone...
      start  = (count == 0) ? 0 : start - beg + 1;
      finish = (count == 0) ? 0 : finish - beg;

      double * data         = (count > 0) ? &tmp[offset] : nullptr;
      cgsize_t range_min[1] = {start};
      cgsize_t range_max[1] = {finish};

      CGCHECK2(cgp_field_read_data(filePtr, base, zone, solution_index, field_offset, range_min,
                                   range_max, data));

      offset += count;
      beg = end;
    }
    communicate_node_data(TOPTR(tmp), ioss_data, 1);
  }

  template void DecompositionData<int>::get_sideset_element_side(
      int filePtr, const Ioss::SetDecompositionData &sset, int *data) const;
  template void DecompositionData<int64_t>::get_sideset_element_side(
      int filePtr, const Ioss::SetDecompositionData &sset, int64_t *data) const;
  template <typename INT>
  void DecompositionData<INT>::get_sideset_element_side(int                               filePtr,
                                                        const Ioss::SetDecompositionData &sset,
                                                        INT *ioss_data) const
  {
    std::vector<INT> element_side;
    int              base = 1;

    int                   nodes_per_face = 4; // FIXME: sb->topology()->number_nodes();
    std::vector<cgsize_t> nodes(nodes_per_face * sset.file_count());

    // TODO? Possibly rewrite using cgi_read_int_data so can skip reading element connectivity

    // We get:
    // *  num_to_get parent elements,
    // *  num_to_get zeros (other parent element for face, but on boundary so 0)
    // *  num_to_get face_on_element
    // *  num_to_get zeros (face on other parent element)
    std::vector<cgsize_t> parent(4 * sset.file_count());

    CGCHECK2(
        cg_elements_read(filePtr, base, sset.zone(), sset.section(), TOPTR(nodes), TOPTR(parent)));
    // Get rid of 'nodes' list -- not used.
    Ioss::Utils::clear(nodes);

    // Move from 'parent' to 'element_side' and interleave. element, side, element, side, ...
    element_side.reserve(sset.file_count() * 2);
    size_t zone_element_id_offset = m_zones[sset.zone()].m_elementOffset;
    for (size_t i = 0; i < sset.file_count(); i++) {
      element_side.push_back(parent[0 * sset.file_count() + i] + zone_element_id_offset);
      element_side.push_back(parent[2 * sset.file_count() + i]);
    }
    // The above was all on root processor for this side set, now need to send data to other
    // processors that own any of the elements in the sideset.
    communicate_set_data(TOPTR(element_side), ioss_data, sset, 2);
  }

  template void DecompositionData<int>::get_block_connectivity(int filePtr, int *data,
                                                               int blk_seq) const;
  template void DecompositionData<int64_t>::get_block_connectivity(int filePtr, int64_t *data,
                                                                   int blk_seq) const;

  template <typename INT>
  void DecompositionData<INT>::get_block_connectivity(int filePtr, INT *data, int blk_seq) const
  {
    auto                  blk = m_elementBlocks[blk_seq];
    std::vector<cgsize_t> file_conn(blk.file_count() * blk.nodesPerEntity);
    int                   base = 1;
    CGCHECK2(cgp_elements_read_data(filePtr, base, blk.zone(), blk.section(), blk.fileSectionOffset,
                                    blk.fileSectionOffset + blk.file_count() - 1,
                                    TOPTR(file_conn)));

    // Map from zone-local node numbers to global implicit
    for (auto &node : file_conn) {
      node += blk.zoneNodeOffset;
    }

    if (!m_zoneSharedMap.empty()) {
      for (auto &node : file_conn) {
        std::unordered_map<cgsize_t, cgsize_t>::const_iterator alias =
            m_zoneSharedMap.find(node - 1);
        if (alias != m_zoneSharedMap.end()) {
          node = (*alias).second + 1;
        }
      }
    }

    communicate_block_data(TOPTR(file_conn), data, blk, (size_t)blk.nodesPerEntity);
  }

  template void DecompositionData<int>::get_element_field(int filePtr, int solution_index,
                                                          int blk_seq, int field_index,
                                                          double *data) const;
  template void DecompositionData<int64_t>::get_element_field(int filePtr, int solution_index,
                                                              int blk_seq, int field_index,
                                                              double *data) const;

  template <typename INT>
  void DecompositionData<INT>::get_element_field(int filePtr, int solution_index, int blk_seq,
                                                 int field_index, double *data) const
  {
    const auto          blk = m_elementBlocks[blk_seq];
    std::vector<double> cgns_data(blk.file_count());
    int                 base         = 1;
    cgsize_t            range_min[1] = {(cgsize_t)blk.fileSectionOffset};
    cgsize_t            range_max[1] = {(cgsize_t)(blk.fileSectionOffset + blk.file_count() - 1)};

    CGCHECK2(cgp_field_read_data(filePtr, base, blk.zone(), solution_index, field_index, range_min,
                                 range_max, cgns_data.data()));

    communicate_block_data(cgns_data.data(), data, blk, (size_t)1);
  }

  DecompositionDataBase::~DecompositionDataBase()
  {
    for (auto &zone : m_structuredZones) {
      delete zone;
    }
  }

  template void DecompositionDataBase::communicate_node_data(int *file_data, int *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(int64_t *file_data, int64_t *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(double *file_data, double *ioss_data,
                                                             size_t comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_node_data(T *file_data, T *ioss_data,
                                                    size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_node_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_node_data(file_data, ioss_data, comp_count);
    }
  }

  template void DecompositionDataBase::communicate_element_data(int *file_data, int *ioss_data,
                                                                size_t comp_count) const;
  template void DecompositionDataBase::communicate_element_data(int64_t *file_data,
                                                                int64_t *ioss_data,
                                                                size_t   comp_count) const;
  template void DecompositionDataBase::communicate_element_data(double *file_data,
                                                                double *ioss_data,
                                                                size_t  comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_element_data(T *file_data, T *ioss_data,
                                                       size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_element_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_element_data(file_data, ioss_data, comp_count);
    }
  }

  void DecompositionDataBase::get_node_entity_proc_data(void *                    entity_proc,
                                                        const Ioss::MapContainer &node_map,
                                                        bool                      do_map) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->m_decomposition.get_node_entity_proc_data((int *)entity_proc, node_map, do_map);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->m_decomposition.get_node_entity_proc_data((int64_t *)entity_proc, node_map, do_map);
    }
  }

  void DecompositionDataBase::get_block_connectivity(int filePtr, void *data, int blk_seq) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_block_connectivity(filePtr, (int *)data, blk_seq);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_block_connectivity(filePtr, (int64_t *)data, blk_seq);
    }
  }

  void DecompositionDataBase::get_element_field(int filePtr, int solution_index, int blk_seq,
                                                int field_index, double *data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_element_field(filePtr, solution_index, blk_seq, field_index, data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_element_field(filePtr, solution_index, blk_seq, field_index, data);
    }
  }

  void DecompositionDataBase::get_node_field(int filePtr, int step, int field_index,
                                             double *data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_node_field(filePtr, step, field_index, data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_node_field(filePtr, step, field_index, data);
    }
  }

  void DecompositionDataBase::get_sideset_element_side(int                               filePtr,
                                                       const Ioss::SetDecompositionData &sset,
                                                       void *                            data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_sideset_element_side(filePtr, sset, (int *)data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_sideset_element_side(filePtr, sset, (int64_t *)data);
    }
  }
} // namespace Iocgns
