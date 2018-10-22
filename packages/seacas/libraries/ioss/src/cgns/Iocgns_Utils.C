// Copyright(C) 1999-2018 National Technology & Engineering Solutions
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

#include <Ioss_Bar2.h>
#include <Ioss_Bar3.h>
#include <Ioss_Hex20.h>
#include <Ioss_Hex27.h>
#include <Ioss_Hex8.h>
#include <Ioss_Node.h>
#include <Ioss_Pyramid13.h>
#include <Ioss_Pyramid14.h>
#include <Ioss_Pyramid5.h>
#include <Ioss_Quad4.h>
#include <Ioss_Quad8.h>
#include <Ioss_Quad9.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Tet10.h>
#include <Ioss_Tet4.h>
#include <Ioss_Tri3.h>
#include <Ioss_Tri4.h>
#include <Ioss_Tri6.h>
#include <Ioss_Unknown.h>
#include <Ioss_Wedge15.h>
#include <Ioss_Wedge18.h>
#include <Ioss_Wedge6.h>

#include <numeric>
#include <set>

#include <cgns/Iocgns_StructuredZoneData.h>
#include <cgns/Iocgns_Utils.h>

#include <cgnsconfig.h>
#if CG_BUILD_PARALLEL
#include <pcgnslib.h>
#else
#include <cgnslib.h>
#endif

#define CGERR(funcall)                                                                             \
  if ((funcall) != CG_OK) {                                                                        \
    Iocgns::Utils::cgns_error(file_ptr, __FILE__, __func__, __LINE__, -1);                         \
  }

namespace {
  int power_2(int count)
  {
    // Return the maximum power of two which is less than or equal to 'count'
    // count = 15 -> returns 8
    // count = 16 -> returns 16
    // count = 17 -> returns 16

    // Use brute force...
    int pow2 = 1;
    while (pow2 <= count) {
      pow2 *= 2;
    }
    if (pow2 > count) {
      pow2 /= 2;
    }
    return pow2;
  }

  struct Range
  {
    Range(int a, int b) : m_beg(a < b ? a : b), m_end(a < b ? b : a), m_reversed(b < a) {}

    int  m_beg;
    int  m_end;
    bool m_reversed;
  };

  bool overlaps(const Range &a, const Range &b) { return a.m_beg <= b.m_end && b.m_beg <= a.m_end; }

  bool bc_overlaps(const Ioss::StructuredBlock *block, const Ioss::BoundaryCondition &bc)
  {
    Ioss::IJK_t ordinal;
    ordinal[0] = block->get_property("ni").get_int();
    ordinal[1] = block->get_property("nj").get_int();
    ordinal[2] = block->get_property("nk").get_int();

    if (ordinal[0] == 0 && ordinal[1] == 0 && ordinal[2] == 0) {
      return false;
    }

    Ioss::IJK_t offset;
    offset[0] = block->get_property("offset_i").get_int();
    offset[1] = block->get_property("offset_j").get_int();
    offset[2] = block->get_property("offset_k").get_int();

    // Note that block range is nodes and m_ordinal[] is cells, so need to add 1 to range.
    Range z_i(1 + offset[0], ordinal[0] + offset[0] + 1);
    Range z_j(1 + offset[1], ordinal[1] + offset[1] + 1);
    Range z_k(1 + offset[2], ordinal[2] + offset[2] + 1);

    Range gc_i(bc.m_rangeBeg[0], bc.m_rangeEnd[0]);
    Range gc_j(bc.m_rangeBeg[1], bc.m_rangeEnd[1]);
    Range gc_k(bc.m_rangeBeg[2], bc.m_rangeEnd[2]);

    return overlaps(z_i, gc_i) && overlaps(z_j, gc_j) && overlaps(z_k, gc_k);
  }

  Range subset_range(const Range &a, const Range &b)
  {
    Range ret(std::max(a.m_beg, b.m_beg), std::min(a.m_end, b.m_end));
    ret.m_reversed = a.m_reversed || b.m_reversed;
    return ret;
  }

  void bc_subset_range(const Ioss::StructuredBlock *block, Ioss::BoundaryCondition &bc)
  {
    if (bc_overlaps(block, bc)) {
      Ioss::IJK_t ordinal;
      ordinal[0] = block->get_property("ni").get_int();
      ordinal[1] = block->get_property("nj").get_int();
      ordinal[2] = block->get_property("nk").get_int();

      Ioss::IJK_t offset;
      offset[0] = block->get_property("offset_i").get_int();
      offset[1] = block->get_property("offset_j").get_int();
      offset[2] = block->get_property("offset_k").get_int();

      // NOTE: Updates the range in bc

      // Note that block range is nodes and m_ordinal[] is cells, so need to add 1 to range.
      Range z_i(1 + offset[0], ordinal[0] + offset[0] + 1);
      Range z_j(1 + offset[1], ordinal[1] + offset[1] + 1);
      Range z_k(1 + offset[2], ordinal[2] + offset[2] + 1);

      Range gc_i(bc.m_rangeBeg[0], bc.m_rangeEnd[0]);
      Range gc_j(bc.m_rangeBeg[1], bc.m_rangeEnd[1]);
      Range gc_k(bc.m_rangeBeg[2], bc.m_rangeEnd[2]);

      Range gc_ii = subset_range(z_i, gc_i);
      Range gc_jj = subset_range(z_j, gc_j);
      Range gc_kk = subset_range(z_k, gc_k);

      bc.m_rangeBeg[0] = gc_ii.m_reversed ? gc_ii.m_end : gc_ii.m_beg;
      bc.m_rangeEnd[0] = gc_ii.m_reversed ? gc_ii.m_beg : gc_ii.m_end;
      bc.m_rangeBeg[1] = gc_jj.m_reversed ? gc_jj.m_end : gc_jj.m_beg;
      bc.m_rangeEnd[1] = gc_jj.m_reversed ? gc_jj.m_beg : gc_jj.m_end;
      bc.m_rangeBeg[2] = gc_kk.m_reversed ? gc_kk.m_end : gc_kk.m_beg;
      bc.m_rangeEnd[2] = gc_kk.m_reversed ? gc_kk.m_beg : gc_kk.m_end;
    }
    else {
      bc.m_rangeBeg = {{0, 0, 0}};
      bc.m_rangeEnd = {{0, 0, 0}};
    }
  }

  int extract_trailing_int(char *name)
  {
    // 'name' consists of an arbitrary number of characters followed by
    // zero or more digits.  Return the integer value of the contiguous
    // set of trailing digits.
    // Example: Name42 returns 42;  Name_52or_perhaps_3_43 returns 43.

    size_t len   = std::strlen(name);
    int    nstep = 0;
    int    mul   = 1;
    for (size_t d = len; d > 0; d--) {
      if (isdigit(name[d - 1])) {
        nstep += mul * (name[d - 1] - '0');
        mul *= 10;
      }
      else {
        break;
      }
    }
    return nstep;
  }

  ssize_t proc_with_minimum_work(Iocgns::StructuredZoneData *zone, const std::vector<size_t> &work,
                                 std::set<std::pair<int, int>> &proc_adam_map)
  {
    size_t  min_work = std::numeric_limits<size_t>::max();
    ssize_t min_proc = -1;
    for (ssize_t i = 0; i < (ssize_t)work.size(); i++) {
      if (work[i] < min_work &&
          proc_adam_map.find(std::make_pair(zone->m_adam->m_zone, i)) == proc_adam_map.end()) {
        min_work = work[i];
        min_proc = i;
        if (min_work == 0) {
          break;
        }
      }
    }
    return min_proc;
  }
  void validate_blocks(const Ioss::StructuredBlockContainer &structured_blocks) {}

  void add_bc_to_block(Ioss::StructuredBlock *block, const std::string &boco_name,
                       const std::string &fam_name, int ibc, cgsize_t *range, CG_BCType_t bocotype,
                       bool is_parallel_io)
  {
    Ioss::SideSet *sset = block->get_database()->get_region()->get_sideset(fam_name);
    if (sset == nullptr) {
      if (block->get_database()->parallel_rank() == 0) {
        IOSS_WARNING << "On block " << block->name() << ", found the boundary condition named "
                     << boco_name << " in family " << fam_name
                     << ". This family was not previously defined at the top-level of the file"
                     << " which is not normal.  Check your file to make sure this does not "
                        "incdicate a problem "
                     << "with the mesh.\n";
      }

      // Need to create a new sideset since didn't see this earlier.
      auto *db = block->get_database();
      sset     = new Ioss::SideSet(db, fam_name);
      if (sset == nullptr) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: Could not create sideset named '" << fam_name << "' on block '"
               << block->name() << "'.\n";
        IOSS_ERROR(errmsg);
      }
      // Get all previous sidesets to make sure we set a unique id...
      int64_t     max_id   = 0;
      const auto &sidesets = db->get_region()->get_sidesets();
      for (const auto &ss : sidesets) {
        if (ss->property_exists("id")) {
          auto id = ss->get_property("id").get_int();
          max_id  = (id > max_id) ? id : max_id;
        }
      }
      sset->property_add(Ioss::Property("id", max_id + 10));
      sset->property_add(Ioss::Property("guid", db->util().generate_guid(max_id + 10)));
      db->get_region()->add(sset);
    }

    if (sset != nullptr) {
      Ioss::IJK_t range_beg{{(int)std::min(range[0], range[3]), (int)std::min(range[1], range[4]),
                             (int)std::min(range[2], range[5])}};

      Ioss::IJK_t range_end{{(int)std::max(range[0], range[3]), (int)std::max(range[1], range[4]),
                             (int)std::max(range[2], range[5])}};

      // Determine overlap of surface with block (in parallel, a block may
      // be split among multiple processors and the block face this is applied
      // to may not exist on this decomposed block)
      auto        bc   = Ioss::BoundaryCondition(boco_name, fam_name, range_beg, range_end);
      std::string name = std::string(boco_name) + "/" + block->name();

      bc_subset_range(block, bc);
      if (!is_parallel_io && !bc.is_valid()) {
        bc.m_rangeBeg = {{0, 0, 0}};
        bc.m_rangeEnd = {{0, 0, 0}};
      }
      block->m_boundaryConditions.push_back(bc);
      auto sb =
          new Ioss::SideBlock(block->get_database(), name, Ioss::Quad4::name, Ioss::Hex8::name,
                              block->m_boundaryConditions.back().get_face_count());
      sb->set_parent_block(block);
      sset->add(sb);

      int base = block->get_property("base").get_int();
      int zone = block->get_property("zone").get_int();
      sb->property_add(Ioss::Property("base", base));
      sb->property_add(Ioss::Property("zone", zone));
      sb->property_add(Ioss::Property("section", ibc + 1));
      sb->property_add(Ioss::Property("id", sset->get_property("id").get_int()));
      sb->property_add(Ioss::Property(
          "guid", block->get_database()->util().generate_guid(sset->get_property("id").get_int())));

      // Set a property on the sideset specifying the boundary condition type (bocotype)
      // In CGNS, the bocotype is an enum; we store it as the integer value of the enum.
      if (sset->property_exists("bc_type")) {
        // Check that the 'bocotype' value matches the value of the property.
        auto old_bocotype = sset->get_property("bc_type").get_int();
        if (old_bocotype != bocotype && bocotype != CG_FamilySpecified) {
          IOSS_WARNING << "On sideset " << sset->name()
                       << ", the boundary condition type was previously set to " << old_bocotype
                       << " which does not match the current value of " << bocotype
                       << ". It will keep the old value.\n";
        }
      }
      else {
        sset->property_add(Ioss::Property("bc_type", (int)bocotype));
      }
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: StructuredBlock '" << block->name()
             << "' Did not find matching sideset with name '" << boco_name << "'";
      IOSS_ERROR(errmsg);
    }
  }

  void sync_transient_variables_fpp(Ioss::Region *region, int myProcessor)
  {
    // With an fpp read, certain blocks may only be on certain
    // processors -- This consistency is addressed elsewhere; however,
    // if a block is not on a processor, then that block will not have
    // any transient fields.  Need to sync across all processors such
    // that a block has the same fields on all processors.
    //
    // ASSUME: A block will have the same fields in the same order on
    // all processors on which it exists.
    //
    // Do the gather all metadata to proc 0; consolidate and then
    // broadcast back...
    // Need: 'name' and 'VariableType'.  Assume all are double and the
    // size will be processor dependent.
    auto &           sblocks = region->get_structured_blocks();
    std::vector<int> fld_count;
    fld_count.reserve(sblocks.size());
    for (const auto &block : sblocks) {
      fld_count.push_back(block->field_count(Ioss::Field::TRANSIENT));
    }
    auto par = region->get_database()->util();
    par.global_array_minmax(fld_count, Ioss::ParallelUtils::DO_MAX);

    // Determine total number of fields on all blocks...
    int tot_fld = std::accumulate(fld_count.begin(), fld_count.end(), 0);
    // Assuming fields are the same on all processors that have fields...
    std::vector<char> fld_names(tot_fld * 2 * (CGNS_MAX_NAME_LENGTH + 1), 0);

    size_t offset = 0;
    for (size_t i = 0; i < sblocks.size(); i++) {
      const auto &   block = sblocks[i];
      Ioss::NameList fields;
      block->field_describe(Ioss::Field::TRANSIENT, &fields);
      if (!fields.empty()) {
        for (const auto &field_name : fields) {
          const Ioss::Field &field = block->get_fieldref(field_name);
          std::string        type  = field.raw_storage()->name();
          strncpy(&fld_names[offset], field_name.c_str(), CGNS_MAX_NAME_LENGTH);
          offset += CGNS_MAX_NAME_LENGTH + 1;
          strncpy(&fld_names[offset], type.c_str(), CGNS_MAX_NAME_LENGTH);
          offset += CGNS_MAX_NAME_LENGTH + 1;
        }
      }
      else {
        offset += (CGNS_MAX_NAME_LENGTH + 1) * 2 * fld_count[i];
      }
    }

    par.global_array_minmax(fld_names, Ioss::ParallelUtils::DO_MAX);

    // Each processor now should have a consistent list of the field
    // names.  Now need to add the missing fields to the blocks that
    // are not 'native' to this processor...
    //
    for (size_t i = 0; i < sblocks.size(); i++) {
      auto &block = sblocks[i];
      if (block->field_count(Ioss::Field::TRANSIENT) != (size_t)fld_count[i]) {
        // Verify that either has 0 or correct number of fields...
        assert(block->field_count(Ioss::Field::TRANSIENT) == 0);

        // Extract the field name and storage type...
        offset = (CGNS_MAX_NAME_LENGTH + 1) * 2 * i;

        for (int nf = 0; nf < fld_count[i]; nf++) {
          std::string fld_name(&fld_names[offset]);
          offset += CGNS_MAX_NAME_LENGTH + 1;
          std::string fld_type(&fld_names[offset]);
          offset += CGNS_MAX_NAME_LENGTH + 1;

          block->field_add(
              Ioss::Field(fld_name, Ioss::Field::DOUBLE, fld_type, Ioss::Field::TRANSIENT, 0));
        }
      }
      assert(block->field_count(Ioss::Field::TRANSIENT) == (size_t)fld_count[i]);
    }
  }
} // namespace

void Iocgns::Utils::cgns_error(int cgnsid, const char *file, const char *function, int lineno,
                               int processor)
{
  std::ostringstream errmsg;
  errmsg << "CGNS error '" << cg_get_error() << "' at line " << lineno << " in file '" << file
         << "' in function '" << function << "'";
  if (processor >= 0) {
    errmsg << " on processor " << processor;
  }
  errmsg << ". Please report to gdsjaar@sandia.gov if you need help.";
  if (cgnsid > 0) {
#if CG_BUILD_PARALLEL
    // This can cause a hang if not all processors call this routine
    // and then the error is not output...
    //    cgp_close(cgnsid);
#else
    cg_close(cgnsid);
#endif
  }
  IOSS_ERROR(errmsg);
}

CG_ZoneType_t Iocgns::Utils::check_zone_type(int cgnsFilePtr)
{
  // ========================================================================
  // Get the number of zones (element/structured blocks) in the mesh...
  int base      = 1;
  int num_zones = 0;
  CGCHECKNP(cg_nzones(cgnsFilePtr, base, &num_zones));

  CG_ZoneType_t common_zone_type = CG_ZoneTypeNull;

  for (int zone = 1; zone <= num_zones; zone++) {
    CG_ZoneType_t zone_type;
    CGCHECKNP(cg_zone_type(cgnsFilePtr, base, zone, &zone_type));

    if (common_zone_type == CG_ZoneTypeNull) {
      common_zone_type = zone_type;
    }

    if (common_zone_type != zone_type) {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: Zone " << zone << " is not the same zone type as previous zones."
             << " This is currently not allowed or supported (hybrid mesh).";
      IOSS_ERROR(errmsg);
    }
  }
  return common_zone_type;
}

int Iocgns::Utils::get_db_zone(const Ioss::EntityBlock *block)
{
  // Returns the zone of the entity as it appears on the cgns database.
  // Usually, but not always the same as the IOSS zone...
  // Can differ on fpp reads and maybe writes.
  if (block->property_exists("db_zone")) {
    return block->get_property("db_zone").get_int();
  }
  return block->get_property("zone").get_int();
}

size_t Iocgns::Utils::index(const Ioss::Field &field) { return field.get_index() & 0xffffffff; }

void Iocgns::Utils::set_field_index(const Ioss::Field &field, size_t index,
                                    CG_GridLocation_t location)
{
  if (location == CG_CellCenter) {
    index |= CG_CELL_CENTER_FIELD_ID;
  }
  if (location == CG_Vertex) {
    index |= CG_VERTEX_FIELD_ID;
  }
  field.set_index(index);
}

bool Iocgns::Utils::is_cell_field(const Ioss::Field &field)
{
  size_t index = field.get_index();
  if (index & CG_VERTEX_FIELD_ID) {
    return false;
  }
  else if (index & CG_CELL_CENTER_FIELD_ID) {
    return true;
  }
  return !(field.get_name() == "mesh_model_coordinates" ||
           field.get_name() == "mesh_model_coordinates_x" ||
           field.get_name() == "mesh_model_coordinates_y" ||
           field.get_name() == "mesh_model_coordinates_z" ||
           field.get_name() == "cell_node_ids"); // Default to cell field...
}

namespace {
#ifdef SEACAS_HAVE_MPI
  void union_zgc_range(Ioss::ZoneConnectivity &zgc_i, const Ioss::ZoneConnectivity &zgc_j)
  {
    assert(zgc_i.m_transform == zgc_j.m_transform);
    for (int i = 0; i < 3; i++) {
      if (zgc_i.m_transform[i] > 0) {
        zgc_i.m_ownerRangeBeg[i] = std::min(zgc_i.m_ownerRangeBeg[i], zgc_j.m_ownerRangeBeg[i]);
        zgc_i.m_ownerRangeEnd[i] = std::max(zgc_i.m_ownerRangeEnd[i], zgc_j.m_ownerRangeEnd[i]);
      }
      else {
        zgc_i.m_ownerRangeBeg[i] = std::max(zgc_i.m_ownerRangeBeg[i], zgc_j.m_ownerRangeBeg[i]);
        zgc_i.m_ownerRangeEnd[i] = std::min(zgc_i.m_ownerRangeEnd[i], zgc_j.m_ownerRangeEnd[i]);
      }
      zgc_i.m_donorRangeBeg[i] = std::min(zgc_i.m_donorRangeBeg[i], zgc_j.m_donorRangeBeg[i]);
      zgc_i.m_donorRangeEnd[i] = std::max(zgc_i.m_donorRangeEnd[i], zgc_j.m_donorRangeEnd[i]);
    }
  }
#endif

  void consolidate_zgc(const Ioss::Region &region)
  {
    // In parallel, the zgc are not necessarily consistent across processors...
    // and the owner/donor ranges are processor specific.
    // Need to make sure all processors have a consistent list of zgc and the
    // owner/donor ranges contain the union of the ranges on each
    // processor.
    // ...Could do this on a per sb basis, but better to do all at once...
    // Data:
    // CGNS_MAX_NAME_LENGTH - connectionName -- CGNS_MAX_NAME_LENGTH char max
    // 1 - int zone
    // 1 - int donor_zone -- get by mapping donorName to zone
    // 6 cgsize_t[6] ownerRange (can probably use 32-bit int...)
    // 6 cgsize_t[6] donorRange (can probably use 32-bit int...)
    // 3 int[3] transform; (values range from -3 to +3 (could store as single int)
    // CGNS_MAX_NAME_LENGTH characters + 17 ints / connection.

#ifdef SEACAS_HAVE_MPI
    const int BYTE_PER_NAME = CGNS_MAX_NAME_LENGTH;
    const int INT_PER_ZGC   = 17;
    // Gather all to processor 0, consolidate, and then scatter back...
    int         my_count          = 0;
    const auto &structured_blocks = region.get_structured_blocks();
    for (const auto &sb : structured_blocks) {
      my_count += std::count_if(
          sb->m_zoneConnectivity.begin(), sb->m_zoneConnectivity.end(),
          [](const Ioss::ZoneConnectivity &z) { return !z.is_intra_block() && z.is_active(); });
    }

    std::vector<int> rcv_data_cnt;
    region.get_database()->util().all_gather(
        my_count, rcv_data_cnt); // Allgather instead of gather so can bail if count=0
    int count = std::accumulate(rcv_data_cnt.begin(), rcv_data_cnt.end(), 0);
    if (count == 0) {
      for (auto &sb : structured_blocks) {
        sb->m_zoneConnectivity.clear();
      }
      return;
    }

    std::vector<char> snd_zgc_name(my_count * BYTE_PER_NAME);
    std::vector<int>  snd_zgc_data(my_count * INT_PER_ZGC);

    // Pack data for gathering to processor 0...
    int off_name = 0;
    int off_data = 0;
    int off_cnt  = 0;

    // ========================================================================
    auto pack_lambda = [&off_data, &off_name, &off_cnt, &snd_zgc_data,
                        &snd_zgc_name](const std::vector<Ioss::ZoneConnectivity> &zgc) {
      for (const auto &z : zgc) {
        if (!z.is_intra_block() && z.is_active()) {
          strncpy(&snd_zgc_name[off_name], z.m_connectionName.c_str(), BYTE_PER_NAME);
          off_cnt++;
          off_name += BYTE_PER_NAME;

          snd_zgc_data[off_data++] = z.m_ownerZone;
          snd_zgc_data[off_data++] = z.m_donorZone;

          snd_zgc_data[off_data++] = z.m_ownerRangeBeg[0];
          snd_zgc_data[off_data++] = z.m_ownerRangeBeg[1];
          snd_zgc_data[off_data++] = z.m_ownerRangeBeg[2];
          snd_zgc_data[off_data++] = z.m_ownerRangeEnd[0];
          snd_zgc_data[off_data++] = z.m_ownerRangeEnd[1];
          snd_zgc_data[off_data++] = z.m_ownerRangeEnd[2];

          snd_zgc_data[off_data++] = z.m_donorRangeBeg[0];
          snd_zgc_data[off_data++] = z.m_donorRangeBeg[1];
          snd_zgc_data[off_data++] = z.m_donorRangeBeg[2];
          snd_zgc_data[off_data++] = z.m_donorRangeEnd[0];
          snd_zgc_data[off_data++] = z.m_donorRangeEnd[1];
          snd_zgc_data[off_data++] = z.m_donorRangeEnd[2];

          snd_zgc_data[off_data++] = z.m_transform[0];
          snd_zgc_data[off_data++] = z.m_transform[1];
          snd_zgc_data[off_data++] = z.m_transform[2];
        }
      }
    };
    // ========================================================================

    off_data = off_name = off_cnt = 0;
    for (const auto &sb : structured_blocks) {
      pack_lambda(sb->m_zoneConnectivity);
    }
    assert(off_cnt == my_count);
    assert(my_count == 0 || (off_data % my_count == 0));
    assert(my_count == 0 || (off_data / my_count == INT_PER_ZGC));
    assert(my_count == 0 || (off_name % my_count == 0 && off_name / my_count == BYTE_PER_NAME));

    std::vector<char> rcv_zgc_name;
    std::vector<int>  rcv_zgc_data;
    region.get_database()->util().gather(my_count, BYTE_PER_NAME, snd_zgc_name, rcv_zgc_name);
    region.get_database()->util().gather(my_count, INT_PER_ZGC, snd_zgc_data, rcv_zgc_data);

    // Processor 0 now has all the zgc instances from all blocks on all processors.
    std::vector<Ioss::ZoneConnectivity> zgc;
    if (region.get_database()->util().parallel_rank() == 0) {
      zgc.reserve(count);

      // Unpack data...
      off_data = 0;
      off_name = 0;
      for (int i = 0; i < count; i++) {
        std::string name{&rcv_zgc_name[off_name]};
        off_name += BYTE_PER_NAME;
        int         zone  = rcv_zgc_data[off_data++];
        int         donor = rcv_zgc_data[off_data++];
        Ioss::IJK_t range_beg{
            {rcv_zgc_data[off_data++], rcv_zgc_data[off_data++], rcv_zgc_data[off_data++]}};
        Ioss::IJK_t range_end{
            {rcv_zgc_data[off_data++], rcv_zgc_data[off_data++], rcv_zgc_data[off_data++]}};
        Ioss::IJK_t donor_beg{
            {rcv_zgc_data[off_data++], rcv_zgc_data[off_data++], rcv_zgc_data[off_data++]}};
        Ioss::IJK_t donor_end{
            {rcv_zgc_data[off_data++], rcv_zgc_data[off_data++], rcv_zgc_data[off_data++]}};
        Ioss::IJK_t transform{
            {rcv_zgc_data[off_data++], rcv_zgc_data[off_data++], rcv_zgc_data[off_data++]}};
        zgc.emplace_back(name, zone, "", donor, transform, range_beg, range_end, donor_beg,
                         donor_end);
      }
      assert(off_data % count == 0);
      assert(off_data / count == INT_PER_ZGC);
      assert(off_name % count == 0 && off_name / count == BYTE_PER_NAME);

      // Consolidate down to the minimum set that has the union of all ranges.
      for (size_t i = 0; i < zgc.size(); i++) {
        if (zgc[i].m_ownerZone > 0 && zgc[i].m_donorZone > 0) {
          auto owner_zone = zgc[i].m_ownerZone;
          auto donor_zone = zgc[i].m_donorZone;

          for (size_t j = i + 1; j < zgc.size(); j++) {
            if (zgc[j].m_connectionName == zgc[i].m_connectionName &&
                zgc[j].m_ownerZone == owner_zone && zgc[j].m_donorZone == donor_zone) {
              // Found another instance of the "same" zgc...  Union the ranges
              union_zgc_range(zgc[i], zgc[j]);
              assert(zgc[i].is_valid());

              // Flag the 'j' instance so it is processed only this time.
              zgc[j].m_ownerZone = -1;
              zgc[j].m_donorZone = -1;
            }
          }
        }
      }

      // Cull out all 'non-active' zgc instances (owner and donor zone <= 0)
      zgc.erase(std::remove_if(zgc.begin(), zgc.end(),
                               [](Ioss::ZoneConnectivity &z) {
                                 return (z.m_ownerZone == -1 && z.m_donorZone == -1) ||
                                        z.is_intra_block() || !z.is_active();
                               }),
                zgc.end());

      count = (int)zgc.size();
      snd_zgc_name.resize(count * BYTE_PER_NAME);
      snd_zgc_data.resize(count * INT_PER_ZGC);
      // Now have a unique set of zgc over all processors with a union
      // of the ranges on each individual processor.  Pack the data
      // and broadcast back to all processors so all processors can
      // output the same data for Zone Connectivity.
      off_data = off_name = off_cnt = 0;
      pack_lambda(zgc);

      assert(off_cnt == count);
      assert(off_data % count == 0);
      assert(off_data / count == INT_PER_ZGC);
      assert(off_name % count == 0 && off_name / count == BYTE_PER_NAME);
    } // End of processor 0 only processing...

    // Send the list of unique zgc instances to all processors so they can all output.
    MPI_Bcast(&count, 1, MPI_INT, 0, region.get_database()->util().communicator());
    snd_zgc_name.resize(count * BYTE_PER_NAME);
    snd_zgc_data.resize(count * INT_PER_ZGC);
    MPI_Bcast(snd_zgc_name.data(), (int)snd_zgc_name.size(), MPI_BYTE, 0,
              region.get_database()->util().communicator());
    MPI_Bcast(snd_zgc_data.data(), (int)snd_zgc_data.size(), MPI_INT, 0,
              region.get_database()->util().communicator());

    // Now clean out existing ZGC lists for all blocks and add on the consolidated instances.
    // Also create a vector for mapping from zone to sb name.
    std::vector<std::string> sb_names(structured_blocks.size() + 1);
    for (auto &sb : structured_blocks) {
      sb->m_zoneConnectivity.clear();
      auto zone = sb->get_property("zone").get_int();
      assert(zone < (int)sb_names.size());
      sb_names[zone] = sb->name();
    }

    // Unpack data and apply to the correct structured block.
    off_data = 0;
    off_name = 0;
    for (int i = 0; i < count; i++) {
      std::string name{&snd_zgc_name[off_name]};
      off_name += BYTE_PER_NAME;
      int zone = snd_zgc_data[off_data++];
      assert(zone < (int)sb_names.size());
      int donor = snd_zgc_data[off_data++];
      assert(donor < (int)sb_names.size());
      Ioss::IJK_t range_beg{
          {snd_zgc_data[off_data++], snd_zgc_data[off_data++], snd_zgc_data[off_data++]}};
      Ioss::IJK_t range_end{
          {snd_zgc_data[off_data++], snd_zgc_data[off_data++], snd_zgc_data[off_data++]}};
      Ioss::IJK_t donor_beg{
          {snd_zgc_data[off_data++], snd_zgc_data[off_data++], snd_zgc_data[off_data++]}};
      Ioss::IJK_t donor_end{
          {snd_zgc_data[off_data++], snd_zgc_data[off_data++], snd_zgc_data[off_data++]}};
      Ioss::IJK_t transform{
          {snd_zgc_data[off_data++], snd_zgc_data[off_data++], snd_zgc_data[off_data++]}};

      auto sb = structured_blocks[zone - 1];
      assert(sb->get_property("zone").get_int() == zone);
      sb->m_zoneConnectivity.emplace_back(name, zone, sb_names[donor], donor, transform, range_beg,
                                          range_end, donor_beg, donor_end);
    }
#endif
  }
} // namespace

size_t Iocgns::Utils::common_write_meta_data(int file_ptr, const Ioss::Region &region,
                                             std::vector<size_t> &zone_offset, bool is_parallel_io)
{
  // Make sure mesh is not hybrid...
  if (region.mesh_type() == Ioss::MeshType::HYBRID) {
    std::ostringstream errmsg;
    errmsg << "ERROR: CGNS: The mesh on region " << region.name() << " is of type 'hybrid'."
           << " This is currently not allowed or supported.";
    IOSS_ERROR(errmsg);
  }

  int base           = 0;
  int phys_dimension = region.get_property("spatial_dimension").get_int();
  CGERR(cg_base_write(file_ptr, "Base", phys_dimension, phys_dimension, &base));

  CGERR(cg_goto(file_ptr, base, "end"));
  CGERR(cg_descriptor_write("Information", "IOSS: CGNS Writer version -1"));
  CGERR(cg_goto(file_ptr, base, "end"));
  CGERR(cg_dataclass_write(CGNS_ENUMV(Dimensional)));
  CGERR(cg_units_write(CGNS_ENUMV(MassUnitsUserDefined), CGNS_ENUMV(LengthUnitsUserDefined),
                       CGNS_ENUMV(TimeUnitsUserDefined), CGNS_ENUMV(TemperatureUnitsUserDefined),
                       CGNS_ENUMV(AngleUnitsUserDefined)))

  // Output the sidesets as Family_t nodes
  const auto &sidesets = region.get_sidesets();
  for (const auto &ss : sidesets) {
    int fam = 0;
    CGERR(cg_family_write(file_ptr, base, ss->name().c_str(), &fam));

    int         bc_index = 0;
    CG_BCType_t bocotype = CG_BCTypeNull;
    if (ss->property_exists("bc_type")) {
      bocotype = (CG_BCType_t)ss->get_property("bc_type").get_int();
    }

    int64_t id = fam;
    if (ss->property_exists("id")) {
      id = ss->get_property("id").get_int();
    }

    CGERR(cg_fambc_write(file_ptr, base, fam, "FamBC", bocotype, &bc_index));
    CGERR(cg_goto(file_ptr, base, "Family_t", fam, nullptr));
    CGERR(cg_descriptor_write("FamBC_TypeId", std::to_string(bocotype).c_str()));
    CGERR(cg_descriptor_write("FamBC_TypeName", BCTypeName[bocotype]));
    CGERR(cg_descriptor_write("FamBC_UserId", std::to_string(id).c_str()));
    CGERR(cg_descriptor_write("FamBC_UserName", ss->name().c_str()));
  }

  // NOTE: Element Block zone write is deferred to put_field_internal so can
  // generate the node count based on connectivity traversal...
  // Just getting processor element count here...
  const auto &element_blocks = region.get_element_blocks();
  size_t      element_count  = 0;
  for (const auto &eb : element_blocks) {
    int64_t local_count = eb->entity_count();
#ifdef SEACAS_HAVE_MPI
    if (is_parallel_io) {
      int64_t start = 0;
      MPI_Exscan(&local_count, &start, 1, Ioss::mpi_type(start), MPI_SUM,
                 region.get_database()->util().communicator());
      // Of the cells/elements in this zone, this proc handles
      // those starting at 'proc_offset+1' to 'proc_offset+num_entity'
      eb->property_update("proc_offset", start);
    }
#endif
    element_count += (size_t)local_count;
  }

  const auto &structured_blocks = region.get_structured_blocks();
  validate_blocks(structured_blocks);

  // If `is_parallel` and `!is_parallel_io`, then writing file-per-processor
  bool is_parallel = region.get_database()->util().parallel_size() > 1;
  int  rank        = region.get_database()->util().parallel_rank();
  int  zone        = 0;
  for (const auto &sb : structured_blocks) {
    cgsize_t size[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    if (is_parallel_io) {
      size[3] = sb->get_property("ni_global").get_int();
      size[4] = sb->get_property("nj_global").get_int();
      size[5] = sb->get_property("nk_global").get_int();
    }
    else {
      size[3] = sb->get_property("ni").get_int();
      size[4] = sb->get_property("nj").get_int();
      size[5] = sb->get_property("nk").get_int();
    }
    size[0] = size[3] + 1;
    size[1] = size[4] + 1;
    size[2] = size[5] + 1;

    if (is_parallel_io || sb->is_active()) {
      std::string name = sb->name();
      if (is_parallel && !is_parallel_io) {
        name += "_proc-";
        name += std::to_string(rank);
      }
      int db_zone = 0;
      CGERR(cg_zone_write(file_ptr, base, name.c_str(), size, CG_Structured, &db_zone));
      sb->property_update("db_zone", db_zone);

      // Add GridCoordinates Node...
      int grid_idx = 0;
      CGERR(cg_grid_write(file_ptr, base, db_zone, "GridCoordinates", &grid_idx));
    }
    else {
      sb->property_update("db_zone", -1);
    }
    zone++;
    assert(zone > 0);
    zone_offset[zone] = zone_offset[zone - 1] + sb->get_property("cell_count").get_int();
    sb->property_update("zone", zone);
    sb->property_update("base", base);
  }

  if (is_parallel_io || !is_parallel) { // Only for single file output or serial...
    // Create a vector for mapping from sb_name to zone -- used to update zgc instances
    std::map<std::string, int> sb_zone;
    for (const auto &sb : structured_blocks) {
      zone                = sb->get_property("zone").get_int();
      sb_zone[sb->name()] = zone;
    }

    // Update zgc instances to make sure the ownerZone and donorZone are
    // consistent with the zones on the output database (from cg_zone_write call)
    for (const auto &sb : structured_blocks) {
      int owner_zone = sb->get_property("zone").get_int();
      for (auto &zgc : sb->m_zoneConnectivity) {
        int donor_zone  = sb_zone[zgc.m_donorName];
        zgc.m_ownerZone = owner_zone;
        zgc.m_ownerGUID = region.get_database()->util().generate_guid(owner_zone);
        zgc.m_donorZone = donor_zone;
        zgc.m_donorGUID = region.get_database()->util().generate_guid(donor_zone);
      }
    }
  }

  if (is_parallel_io) {
    consolidate_zgc(region);
  }

  for (const auto &sb : structured_blocks) {
    if (!is_parallel_io && !sb->is_active()) {
      continue;
    }

    auto        db_zone = get_db_zone(sb);
    std::string name    = sb->name();
    if (is_parallel && !is_parallel_io) {
      name += "_proc-";
      name += std::to_string(rank);
    }

    // Transfer boundary condition nodes...
    // The bc.m_ownerRange argument needs to be the union of the size on all processors
    // Instead of requiring that of the caller, do the union in this routine.
    // TODO: Calculate it outside of the loop...
    // Need to handle possible range == 0,0,0.  Only affects the beg data...
    std::vector<cgsize_t> bc_range(sb->m_boundaryConditions.size() * 6);
    size_t                idx = 0;
    for (const auto &bc : sb->m_boundaryConditions) {
      for (size_t i = 0; i < 3; i++) {
        if (bc.m_rangeBeg[i] == 0) {
          bc_range[idx++] = std::numeric_limits<int>::min();
        }
        else {
          bc_range[idx++] = -bc.m_rangeBeg[i];
        }
      }
      for (size_t i = 0; i < 3; i++) {
        bc_range[idx++] = bc.m_rangeEnd[i];
      }
    }

    if (is_parallel_io) {
      region.get_database()->util().global_array_minmax(bc_range, Ioss::ParallelUtils::DO_MAX);
    }

    for (idx = 0; idx < bc_range.size(); idx += 6) {
      bc_range[idx + 0] = -bc_range[idx + 0];
      bc_range[idx + 1] = -bc_range[idx + 1];
      bc_range[idx + 2] = -bc_range[idx + 2];
    }

    Ioss::IJK_t offset;
    offset[0] = sb->get_property("offset_i").get_int();
    offset[1] = sb->get_property("offset_j").get_int();
    offset[2] = sb->get_property("offset_k").get_int();

    idx = 0;
    for (const auto &bc : sb->m_boundaryConditions) {
      int bc_idx = 0;
      if (!is_parallel_io) {
        bc_range[idx + 0] -= offset[0];
        bc_range[idx + 1] -= offset[1];
        bc_range[idx + 2] -= offset[2];
        bc_range[idx + 3] -= offset[0];
        bc_range[idx + 4] -= offset[1];
        bc_range[idx + 5] -= offset[2];
      }

      if (is_parallel_io ||
          (bc_range[idx + 3] > 0 && bc_range[idx + 4] > 0 && bc_range[idx + 5] > 0)) {
        CGERR(cg_boco_write(file_ptr, base, db_zone, bc.m_bcName.c_str(), CG_FamilySpecified,
                            CG_PointRange, 2, &bc_range[idx], &bc_idx));
        CGERR(
            cg_goto(file_ptr, base, name.c_str(), 0, "ZoneBC_t", 1, bc.m_bcName.c_str(), 0, "end"));
        CGERR(cg_famname_write(bc.m_famName.c_str()));
        CGERR(cg_boco_gridlocation_write(file_ptr, base, db_zone, bc_idx, CG_Vertex));
      }
      idx += 6;
    }
    // Transfer Zone Grid Connectivity...
    for (const auto &zgc : sb->m_zoneConnectivity) {
      if (zgc.is_valid() && zgc.is_active()) {
        int                zgc_idx = 0;
        std::array<INT, 6> owner_range{{zgc.m_ownerRangeBeg[0], zgc.m_ownerRangeBeg[1],
                                        zgc.m_ownerRangeBeg[2], zgc.m_ownerRangeEnd[0],
                                        zgc.m_ownerRangeEnd[1], zgc.m_ownerRangeEnd[2]}};
        std::array<INT, 6> donor_range{{zgc.m_donorRangeBeg[0], zgc.m_donorRangeBeg[1],
                                        zgc.m_donorRangeBeg[2], zgc.m_donorRangeEnd[0],
                                        zgc.m_donorRangeEnd[1], zgc.m_donorRangeEnd[2]}};

        std::string donor_name   = zgc.m_donorName;
        std::string connect_name = zgc.m_connectionName;
        if (is_parallel && !is_parallel_io) {
          if (zgc.is_intra_block()) {
            connect_name = std::to_string(zgc.m_ownerGUID) + "--" + std::to_string(zgc.m_donorGUID);
          }
          else {
            if (zgc.m_ownerProcessor != zgc.m_donorProcessor) {
              connect_name += "_proc-" + std::to_string(zgc.m_donorProcessor);
            }
          }
          donor_name += "_proc-";
          donor_name += std::to_string(zgc.m_donorProcessor);
          owner_range[0] -= zgc.m_ownerOffset[0];
          owner_range[1] -= zgc.m_ownerOffset[1];
          owner_range[2] -= zgc.m_ownerOffset[2];
          owner_range[3] -= zgc.m_ownerOffset[0];
          owner_range[4] -= zgc.m_ownerOffset[1];
          owner_range[5] -= zgc.m_ownerOffset[2];

          donor_range[0] -= zgc.m_donorOffset[0];
          donor_range[1] -= zgc.m_donorOffset[1];
          donor_range[2] -= zgc.m_donorOffset[2];
          donor_range[3] -= zgc.m_donorOffset[0];
          donor_range[4] -= zgc.m_donorOffset[1];
          donor_range[5] -= zgc.m_donorOffset[2];
        }
        CGERR(cg_1to1_write(file_ptr, base, db_zone, connect_name.c_str(), donor_name.c_str(),
                            owner_range.data(), donor_range.data(), zgc.m_transform.data(),
                            &zgc_idx));
      }
    }
  }
  return element_count;
}

std::string Iocgns::Utils::map_cgns_to_topology_type(CG_ElementType_t type)
{
  std::string topology = "unknown";
  switch (type) {
  case CG_NODE: topology = Ioss::Node::name; break;
  case CG_BAR_2: topology = Ioss::Bar2::name; break;
  case CG_BAR_3: topology = Ioss::Bar3::name; break;
  case CG_TRI_3: topology = Ioss::Tri3::name; break;
  case CG_TRI_6: topology = Ioss::Tri6::name; break;
  case CG_QUAD_4: topology = Ioss::Quad4::name; break;
  case CG_QUAD_8: topology = Ioss::Quad8::name; break;
  case CG_QUAD_9: topology = Ioss::Quad9::name; break;
  case CG_TETRA_4: topology = Ioss::Tet4::name; break;
  case CG_TETRA_10: topology = Ioss::Tet10::name; break;
  case CG_PYRA_5: topology = Ioss::Pyramid5::name; break;
  case CG_PYRA_13: topology = Ioss::Pyramid13::name; break;
  case CG_PYRA_14: topology = Ioss::Pyramid14::name; break;
  case CG_PENTA_6: topology = Ioss::Wedge6::name; break;
  case CG_PENTA_15: topology = Ioss::Wedge15::name; break;
  case CG_PENTA_18: topology = Ioss::Wedge18::name; break;
  case CG_HEXA_8: topology = Ioss::Hex8::name; break;
  case CG_HEXA_20: topology = Ioss::Hex20::name; break;
  case CG_HEXA_27: topology = Ioss::Hex27::name; break;
  default:
    std::cerr << "WARNING: Found topology of type " << cg_ElementTypeName(type)
              << " which is not currently supported.\n";
    topology = Ioss::Unknown::name;
  }
  return topology;
}

CG_ElementType_t Iocgns::Utils::map_topology_to_cgns(const std::string &name)
{
  CG_ElementType_t topo = CG_ElementTypeNull;
  if (name == Ioss::Node::name) {
    topo = CG_NODE;
  }
  else if (name == Ioss::Bar2::name) {
    topo = CG_BAR_2;
  }
  else if (name == Ioss::Bar3::name) {
    topo = CG_BAR_3;
  }
  else if (name == Ioss::Tri3::name) {
    topo = CG_TRI_3;
  }
  else if (name == Ioss::Tri6::name) {
    topo = CG_TRI_6;
  }
  else if (name == Ioss::Quad4::name) {
    topo = CG_QUAD_4;
  }
  else if (name == Ioss::Quad8::name) {
    topo = CG_QUAD_8;
  }
  else if (name == Ioss::Quad9::name) {
    topo = CG_QUAD_9;
  }
  else if (name == Ioss::Tet4::name) {
    topo = CG_TETRA_4;
  }
  else if (name == Ioss::Tet10::name) {
    topo = CG_TETRA_10;
  }
  else if (name == Ioss::Pyramid5::name) {
    topo = CG_PYRA_5;
  }
  else if (name == Ioss::Pyramid13::name) {
    topo = CG_PYRA_13;
  }
  else if (name == Ioss::Pyramid14::name) {
    topo = CG_PYRA_14;
  }
  else if (name == Ioss::Wedge6::name) {
    topo = CG_PENTA_6;
  }
  else if (name == Ioss::Wedge15::name) {
    topo = CG_PENTA_15;
  }
  else if (name == Ioss::Wedge18::name) {
    topo = CG_PENTA_18;
  }
  else if (name == Ioss::Hex8::name) {
    topo = CG_HEXA_8;
  }
  else if (name == Ioss::Hex20::name) {
    topo = CG_HEXA_20;
  }
  else if (name == Ioss::Hex27::name) {
    topo = CG_HEXA_27;
  }
  else {
    std::cerr << "WARNING: Found topology of type " << name
              << " which is not currently supported.\n";
  }
  return topo;
}

void Iocgns::Utils::write_flow_solution_metadata(int file_ptr, Ioss::Region *region, int state,
                                                 int *vertex_solution_index,
                                                 int *cell_center_solution_index,
                                                 bool is_parallel_io)
{
  std::string c_name = "CellCenterSolutionAtStep";
  std::string v_name = "VertexSolutionAtStep";
  std::string step   = std::to_string(state);
  c_name += step;
  v_name += step;

  const auto &nblocks          = region->get_node_blocks();
  const auto &nblock           = nblocks[0];
  bool        has_nodal_fields = nblock->field_count(Ioss::Field::TRANSIENT) > 0;

  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.
  auto sol_lambda = [=](Ioss::EntityBlock *block) {
    int base = block->get_property("base").get_int();
    int zone = get_db_zone(block);
    if (has_nodal_fields) {
      CGERR(cg_sol_write(file_ptr, base, zone, v_name.c_str(), CG_Vertex, vertex_solution_index));
      CGERR(
          cg_goto(file_ptr, base, "Zone_t", zone, "FlowSolution_t", *vertex_solution_index, "end"));
      CGERR(cg_gridlocation_write(CG_Vertex));
      CGERR(cg_descriptor_write("Step", step.c_str()));
    }
    if (block->field_count(Ioss::Field::TRANSIENT) > 0) {
      CGERR(cg_sol_write(file_ptr, base, zone, c_name.c_str(), CG_CellCenter,
                         cell_center_solution_index));
      CGERR(cg_goto(file_ptr, base, "Zone_t", zone, "FlowSolution_t", *cell_center_solution_index,
                    "end"));
      CGERR(cg_descriptor_write("Step", step.c_str()));
    }
  };

  // Use the lambda
  const auto &sblocks = region->get_structured_blocks();
  for (auto &block : sblocks) {
    if (is_parallel_io || block->is_active()) {
      sol_lambda(block);
    }
  }
  // Use the lambda
  const auto &eblocks = region->get_element_blocks();
  for (auto &block : eblocks) {
    sol_lambda(block);
  }
}

int Iocgns::Utils::find_solution_index(int cgnsFilePtr, int base, int zone, int step,
                                       CG_GridLocation_t location)
{
  auto str_step = std::to_string(step);
  int  nsols    = 0;
  CGCHECKNP(cg_nsols(cgnsFilePtr, base, zone, &nsols));
  bool location_matches = false;
  for (int i = 0; i < nsols; i++) {
    CG_GridLocation_t db_location;
    char              db_name[CGNS_MAX_NAME_LENGTH + 1];
    CGCHECKNP(cg_sol_info(cgnsFilePtr, base, zone, i + 1, db_name, &db_location));
    if (location == db_location) {
      location_matches = true;
      // Check if steps match.
      // NOTE: Using non-standard "Descriptor_t" node in FlowSolution_t
      CGCHECKNP(cg_goto(cgnsFilePtr, base, "Zone_t", zone, "FlowSolution_t", i + 1, "end"));
      int descriptor_count = 0;
      CGCHECKNP(cg_ndescriptors(&descriptor_count));

      bool found_step_descriptor = false;
      for (int d = 0; d < descriptor_count; d++) {
        char *db_step = nullptr;
        char  name[CGNS_MAX_NAME_LENGTH + 1];
        CGCHECKNP(cg_descriptor_read(d + 1, name, &db_step));
        if (strcmp(name, "step") == 0) {
          found_step_descriptor = true;
          if (str_step == db_step) {
            cg_free(db_step);
            return i + 1;
          }

          cg_free(db_step);
          break; // Found "step" descriptor, but wasn't correct step...
        }
        cg_free(db_step);
      }
      if (!found_step_descriptor) {
        // There was no Descriptor_t node with the name "step",
        // Try to decode the step from the FlowSolution_t name.
        int nstep = extract_trailing_int(db_name);
        if (nstep == step) {
          return i + 1;
        }
      }
    }
  }

  if (nsols == 1 && location_matches) {
    return 1;
  }
  std::cerr << "WARNING: CGNS: Could not find valid solution index for step " << step << ", zone "
            << zone << ", and location " << GridLocationName[location] << "\n";
  return 0;
}

void Iocgns::Utils::add_sidesets(int cgnsFilePtr, Ioss::DatabaseIO *db)
{
  int base         = 1;
  int num_families = 0;
  CGCHECKNP(cg_nfamilies(cgnsFilePtr, base, &num_families));

  for (int family = 1; family <= num_families; family++) {
    char        name[CGNS_MAX_NAME_LENGTH + 1];
    CG_BCType_t bocotype;
    int         num_bc  = 0;
    int         num_geo = 0;
    CGCHECKNP(cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo));

#if IOSS_DEBUG_OUTPUT
    if (db->parallel_rank() == 0) {
      std::cerr << "Family " << family << " named " << name << " has " << num_bc << " BC, and "
                << num_geo << " geometry references\n";
    }
#endif
    if (num_bc > 0) {
      // Create a sideset...
      std::string ss_name(name); // Use name here before cg_fambc_read call overwrites it...

      CGCHECKNP(cg_fambc_read(cgnsFilePtr, base, family, 1, name, &bocotype));

      CGCHECKNP(cg_goto(cgnsFilePtr, base, "Family_t", family, "end"));
      int ndescriptors = 0;
      int id           = 0;
      CGCHECKNP(cg_ndescriptors(&ndescriptors));
      if (ndescriptors > 0) {
        for (int ndesc = 1; ndesc <= ndescriptors; ndesc++) {
          char  dname[CGNS_MAX_NAME_LENGTH + 1];
          char *dtext;
          CGCHECKNP(cg_descriptor_read(ndesc, dname, &dtext));
          if (strcmp(dname, "FamBC_UserId") == 0) {
            // Convert text in `dtext` to integer...
            id = Ioss::Utils::get_number(dtext);
            cg_free(dtext);
            break;
          }
          cg_free(dtext);
        }
      }
      if (id == 0) {
        id = Ioss::Utils::extract_id(ss_name);
      }
      if (id != 0) {
        auto *ss = new Ioss::SideSet(db, ss_name);
        ss->property_add(Ioss::Property("id", id));
        ss->property_add(Ioss::Property("guid", db->util().generate_guid(id)));
        ss->property_add(Ioss::Property("bc_type", bocotype));
        db->get_region()->add(ss);
      }
      else {
        if (db->parallel_rank() == 0) {
          IOSS_WARNING << "*** WARNING: Skipping BC with name '" << ss_name
                       << "' since FamBC_UserId is equal to 0.\n\n";
        }
      }
    }
  }
}

size_t Iocgns::Utils::resolve_nodes(Ioss::Region &region, int my_processor, bool is_parallel)
{
  // Each structured block has its own set of "cell_nodes"
  // At block boundaries, there are duplicate nodes which need to be resolved for the
  // unstructured mesh output.

  // We need to iterate all of the blocks and then each blocks zgc to determine
  // which nodes are shared between blocks. For all shared nodes, the node in the lowest
  // numbered zone is considered the "owner" and all other nodes are shared.

  // At the end of the routine, each block knows where its nodes fit
  // into the implicit ordering of nodes on this processor. This is
  // given by:
  // implicit_location = block->m_blockLocalNodeIndex[i] (0 <= i < #nodes_in_block)
  // Where 0 <= implicit_location < #nodes_on_processor

  // Create a vector of size which is the sum of the on-processor cell_nodes size for each block
  size_t num_total_cell_nodes = 0;
  auto & blocks               = region.get_structured_blocks();
  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    num_total_cell_nodes += node_count;
  }

  ssize_t              ss_max = std::numeric_limits<ssize_t>::max();
  std::vector<ssize_t> cell_node_map(num_total_cell_nodes, ss_max);

  // Each cell_node location in the cell_node_map is currently initialized to ss_max.
  // Iterate each block and then each blocks non-intra-block (i.e., not
  // due to proc decomps) zgc instances and update cell_node_map
  // such that for each shared node, it points to the owner nodes
  // location.
  for (auto &owner_block : blocks) {
    for (const auto &zgc : owner_block->m_zoneConnectivity) {
      if (!zgc.is_intra_block() &&
          zgc.is_active()) { // Not due to processor decomposition and has faces.
        // NOTE: In parallel, the owner block should exist, but may not have
        // any cells on this processor.  We can access its global i,j,k, but
        // don't store or access any "bulk" data on it.
        auto donor_block = region.get_structured_block(zgc.m_donorName);
        assert(donor_block != nullptr);

        std::vector<int> i_range = zgc.get_range(1);
        std::vector<int> j_range = zgc.get_range(2);
        std::vector<int> k_range = zgc.get_range(3);
        for (auto &k : k_range) {
          for (auto &j : j_range) {
            for (auto &i : i_range) {
              Ioss::IJK_t owner_index{{i, j, k}};
              Ioss::IJK_t donor_index = zgc.transform(owner_index);

              // The nodes as 'index' and 'owner' are contiguous and
              // should refer to the same node. 'owner' should be
              // the owner (unless it is already owned by another
              // block)

              ssize_t owner_global_offset = owner_block->get_global_node_offset(owner_index);
              ssize_t donor_global_offset = donor_block->get_global_node_offset(donor_index);

              if (owner_global_offset > donor_global_offset) {
                if (is_parallel && (zgc.m_donorProcessor != my_processor)) {
                  size_t owner_block_local_offset =
                      owner_block->get_block_local_node_offset(owner_index);
                  owner_block->m_globalIdMap.emplace_back(owner_block_local_offset,
                                                          donor_global_offset + 1);
                }
                else if (!is_parallel || (zgc.m_ownerProcessor != my_processor)) {
                  size_t  owner_local_offset = owner_block->get_local_node_offset(owner_index);
                  ssize_t donor_local_offset = donor_block->get_local_node_offset(donor_index);

                  if (cell_node_map[owner_local_offset] == ss_max) {
                    cell_node_map[owner_local_offset] = donor_local_offset;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Now iterate cell_node_map.  If an entry == ss_max, then it is
  // an owned node and needs to have its index into the unstructed
  // mesh node map set; otherwise, the value points to the owner
  // node, so the index at this location should be set to the owner
  // nodes index.
  size_t index = 0;
  for (auto &node : cell_node_map) {
    if (node == ss_max) {
      node = index++;
    }
    else {
      node = -node;
    }
  }

  for (auto &node : cell_node_map) {
    if (node < 0) {
      node = cell_node_map[-node];
    }
  }

  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    block->m_blockLocalNodeIndex.resize(node_count);

    size_t beg = block->get_node_offset();
    size_t end = beg + node_count;
    for (size_t idx = beg, i = 0; idx < end; idx++) {
      block->m_blockLocalNodeIndex[i++] = cell_node_map[idx];
    }
  }
  return index;
}

std::vector<std::vector<std::pair<size_t, size_t>>>
Iocgns::Utils::resolve_processor_shared_nodes(Ioss::Region &region, int my_processor)
{
  // Determine which nodes are shared across processor boundaries.
  // Only need to check on block boundaries..

  // We need to iterate all of the blocks and then each blocks zgc to determine
  // which nodes are shared between processors. For all shared nodes, the node in the lowest
  // numbered zone is considered the "owner" and all other nodes are shared.

  auto &blocks = region.get_structured_blocks();

  std::vector<std::vector<std::pair<size_t, size_t>>> shared_nodes(blocks.size() + 1);

  for (auto &owner_block : blocks) {
    int  owner_zone = owner_block->get_property("zone").get_int();
    auto owner_ids  = owner_block->get_cell_node_ids(true);
    for (const auto &zgc : owner_block->m_zoneConnectivity) {
      assert(zgc.m_donorProcessor >= 0);
      assert(zgc.m_ownerProcessor >= 0);

      if (zgc.is_active() &&
          (zgc.m_donorProcessor != my_processor || zgc.m_ownerProcessor != my_processor)) {
        // NOTE: In parallel, the donor block should exist, but may not have
        // any cells on this processor.  We can access its global i,j,k, but
        // don't store or access any "bulk" data on it.
        auto donor_block = region.get_structured_block(zgc.m_donorName);
        assert(donor_block != nullptr);
        int donor_zone = donor_block->get_property("zone").get_int();

        auto donor_ids = donor_block->get_cell_node_ids(true);

        std::vector<int> i_range = zgc.get_range(1);
        std::vector<int> j_range = zgc.get_range(2);
        std::vector<int> k_range = zgc.get_range(3);
        for (auto &k : k_range) {
          for (auto &j : j_range) {
            for (auto &i : i_range) {
              Ioss::IJK_t owner_index{{i, j, k}};
              Ioss::IJK_t donor_index = zgc.transform(owner_index);

              // The nodes as 'index' and 'owner' are contiguous and
              // should refer to the same node.

              if (my_processor == zgc.m_ownerProcessor) {
                ssize_t owner_offset = owner_block->get_block_local_node_offset(owner_index);
                shared_nodes[owner_zone].emplace_back(owner_offset, zgc.m_donorProcessor);
              }
              else if (my_processor == zgc.m_donorProcessor) {
                ssize_t donor_offset = donor_block->get_block_local_node_offset(donor_index);
                shared_nodes[donor_zone].emplace_back(donor_offset, zgc.m_ownerProcessor);
              }
            }
          }
        }
      }
    }
#if 1 && IOSS_DEBUG_OUTPUT
    std::cerr << "P" << my_processor << ", Block " << owner_block->name()
              << " Shared Nodes: " << shared_nodes[owner_zone].size() << "\n";
#endif
  }
  return shared_nodes;
}

void Iocgns::Utils::add_structured_boundary_conditions(int                    cgnsFilePtr,
                                                       Ioss::StructuredBlock *block,
                                                       bool                   is_parallel_io)
{
  // `is_parallel_io` is true if all processors reading single file.
  // `is_parallel_io` is false if serial, or each processor reading its own file (fpp)
  if (is_parallel_io) {
    add_structured_boundary_conditions_pio(cgnsFilePtr, block);
  }
  else {
    add_structured_boundary_conditions_fpp(cgnsFilePtr, block);
  }
}

void Iocgns::Utils::add_structured_boundary_conditions_pio(int                    cgnsFilePtr,
                                                           Ioss::StructuredBlock *block)
{
  int base = block->get_property("base").get_int();
  int zone = get_db_zone(block);

  // Called by Parallel run reading single file only.
  // The 'cgnsFilePtr' is for the serial file on processor 0.
  // Read all CGNS data on processor 0 and then broadcast to other processors.
  // Data needed:
  // * boco_name (CGNS_MAX_NAME_LENGTH + 1 chars)
  // * fam_name  (CGNS_MAX_NAME_LENGTH + 1 chars)
  // * data     (cgsize_t * 7) (bocotype + range[6])

  int num_bcs = 0;
  int rank    = block->get_database()->util().parallel_rank();
  if (rank == 0) {
    CGCHECKNP(cg_nbocos(cgnsFilePtr, base, zone, &num_bcs));
  }

#ifdef SEACAS_HAVE_MPI
  int proc = block->get_database()->util().parallel_size();
  if (proc > 1) {
    MPI_Bcast(&num_bcs, 1, MPI_INT, 0, block->get_database()->util().communicator());
  }
#endif

  std::vector<int>  bc_data(7 * num_bcs);
  std::vector<char> bc_names(2 * (CGNS_MAX_NAME_LENGTH + 1) * num_bcs);

  if (rank == 0) {
    int      off_data = 0;
    int      off_name = 0;
    cgsize_t range[6];

    for (int ibc = 0; ibc < num_bcs; ibc++) {
      char              boco_name[CGNS_MAX_NAME_LENGTH + 1];
      char              fam_name[CGNS_MAX_NAME_LENGTH + 1];
      CG_BCType_t       bocotype;
      CG_PointSetType_t ptset_type;
      cgsize_t          npnts;
      cgsize_t          NormalListSize;
      CG_DataType_t     NormalDataType;
      int               ndataset;

      // All we really want from this is 'boco_name'
      CGCHECKNP(cg_boco_info(cgnsFilePtr, base, zone, ibc + 1, boco_name, &bocotype, &ptset_type,
                             &npnts, nullptr, &NormalListSize, &NormalDataType, &ndataset));

      if (bocotype == CG_FamilySpecified) {
        // Get family name associated with this boco_name
        CGCHECKNP(
            cg_goto(cgnsFilePtr, base, "Zone_t", zone, "ZoneBC_t", 1, "BC_t", ibc + 1, "end"));
        CGCHECKNP(cg_famname_read(fam_name));
      }
      else {
        strncpy(fam_name, boco_name, CGNS_MAX_NAME_LENGTH);
      }

      CGCHECKNP(cg_boco_read(cgnsFilePtr, base, zone, ibc + 1, range, nullptr));

      strncpy(&bc_names[off_name], boco_name, CGNS_MAX_NAME_LENGTH + 1);
      off_name += (CGNS_MAX_NAME_LENGTH + 1);
      strncpy(&bc_names[off_name], fam_name, CGNS_MAX_NAME_LENGTH + 1);
      off_name += (CGNS_MAX_NAME_LENGTH + 1);

      bc_data[off_data++] = bocotype;
      bc_data[off_data++] = range[0];
      bc_data[off_data++] = range[1];
      bc_data[off_data++] = range[2];
      bc_data[off_data++] = range[3];
      bc_data[off_data++] = range[4];
      bc_data[off_data++] = range[5];
    }
  }

#ifdef SEACAS_HAVE_MPI
  // Broadcast data to other processors...
  if (proc > 1) {
    MPI_Bcast(bc_names.data(), (int)bc_names.size(), MPI_BYTE, 0,
              block->get_database()->util().communicator());
    MPI_Bcast(bc_data.data(), (int)bc_data.size(), MPI_INT, 0,
              block->get_database()->util().communicator());
  }
#endif

  // Now just unpack the data and run through the same calculations on all processors.
  int off_data = 0;
  int off_name = 0;
  for (int ibc = 0; ibc < num_bcs; ibc++) {
    cgsize_t range[6];

    CG_BCType_t bocotype = (CG_BCType_t)bc_data[off_data++];
    range[0]             = bc_data[off_data++];
    range[1]             = bc_data[off_data++];
    range[2]             = bc_data[off_data++];
    range[3]             = bc_data[off_data++];
    range[4]             = bc_data[off_data++];
    range[5]             = bc_data[off_data++];

    std::string boco_name{&bc_names[off_name]};
    off_name += (CGNS_MAX_NAME_LENGTH + 1);
    std::string fam_name{&bc_names[off_name]};
    off_name += (CGNS_MAX_NAME_LENGTH + 1);

    // There are some BC that are applied on an edge or a vertex;
    // Don't want those (yet?), so filter them out at this time...
    {
      int same_count = (range[0] == range[3] ? 1 : 0) + (range[1] == range[4] ? 1 : 0) +
                       (range[2] == range[5] ? 1 : 0);
      if (same_count != 1) {
        std::cerr << "WARNING: CGNS: Skipping Boundary Condition '" << boco_name << "' on block '"
                  << block->name() << "'. It is applied to "
                  << (same_count == 2 ? "an edge" : "a vertex")
                  << ". This code only supports surfaces.\n";
        continue;
      }
    }

    bool is_parallel_io = true;
    add_bc_to_block(block, boco_name, fam_name, ibc, range, bocotype, is_parallel_io);
  }
}

void Iocgns::Utils::add_structured_boundary_conditions_fpp(int                    cgnsFilePtr,
                                                           Ioss::StructuredBlock *block)
{
  int base = block->get_property("base").get_int();
  int zone = get_db_zone(block);

  // Called by both parallel fpp and serial runs.
  // In parallel, the 'cgnsFilePtr' is specific for each processor

  int num_bcs = 0;
  CGCHECKNP(cg_nbocos(cgnsFilePtr, base, zone, &num_bcs));

  for (int ibc = 0; ibc < num_bcs; ibc++) {
    char              boco_name[CGNS_MAX_NAME_LENGTH + 1];
    char              fam_name[CGNS_MAX_NAME_LENGTH + 1];
    CG_BCType_t       bocotype;
    CG_PointSetType_t ptset_type;
    cgsize_t          npnts;
    cgsize_t          NormalListSize;
    CG_DataType_t     NormalDataType;
    int               ndataset;
    cgsize_t          range[6];

    // All we really want from this is 'boco_name'
    CGCHECKNP(cg_boco_info(cgnsFilePtr, base, zone, ibc + 1, boco_name, &bocotype, &ptset_type,
                           &npnts, nullptr, &NormalListSize, &NormalDataType, &ndataset));

    if (bocotype == CG_FamilySpecified) {
      // Get family name associated with this boco_name
      CGCHECKNP(cg_goto(cgnsFilePtr, base, "Zone_t", zone, "ZoneBC_t", 1, "BC_t", ibc + 1, "end"));
      CGCHECKNP(cg_famname_read(fam_name));
    }
    else {
      strncpy(fam_name, boco_name, CGNS_MAX_NAME_LENGTH);
    }

    CGCHECKNP(cg_boco_read(cgnsFilePtr, base, zone, ibc + 1, range, nullptr));

    // There are some BC that are applied on an edge or a vertex;
    // Don't want those (yet?), so filter them out at this time...
    int same_count = (range[0] == range[3] ? 1 : 0) + (range[1] == range[4] ? 1 : 0) +
                     (range[2] == range[5] ? 1 : 0);
    if (same_count != 1) {
      std::cerr << "WARNING: CGNS: Skipping Boundary Condition '" << boco_name << "' on block '"
                << block->name() << "'. It is applied to "
                << (same_count == 2 ? "an edge" : "a vertex")
                << ". This code only supports surfaces.\n";
      continue;
    }

    int num_proc = block->get_database()->util().parallel_size();
    if (num_proc > 1) {
      // Need to modify range with block offset to put into global space
      Ioss::IJK_t offset;
      offset[0] = block->get_property("offset_i").get_int();
      offset[1] = block->get_property("offset_j").get_int();
      offset[2] = block->get_property("offset_k").get_int();
      range[0] += offset[0];
      range[1] += offset[1];
      range[2] += offset[2];
      range[3] += offset[0];
      range[4] += offset[1];
      range[5] += offset[2];
    }

    bool is_parallel_io = false;
    add_bc_to_block(block, boco_name, fam_name, ibc, range, bocotype, is_parallel_io);
  }
}

void Iocgns::Utils::finalize_database(int cgnsFilePtr, const std::vector<double> &timesteps,
                                      Ioss::Region *region, int myProcessor, bool is_parallel_io)
{
  int base = 1;
  CGCHECK(cg_biter_write(cgnsFilePtr, base, "TimeIterValues", timesteps.size()));

  // Now write the timestep time values...
  CGCHECK(cg_goto(cgnsFilePtr, base, "BaseIterativeData_t", 1, "end"));
  cgsize_t dimtv[1] = {(cgsize_t)timesteps.size()};
  CGCHECK(cg_array_write("TimeValues", CG_RealDouble, 1, dimtv, timesteps.data()));

  // Output the ZoneIterativeData which maps a zones flow solutions to timesteps.
  // One per zone and the number of entries matches the number of timesteps...
  const auto &nblocks = region->get_node_blocks();
  auto &      nblock  = nblocks[0];

  bool has_nodal_fields = nblock->field_count(Ioss::Field::TRANSIENT) > 0;

  cgsize_t          dim[2] = {32, (cgsize_t)timesteps.size()};
  std::vector<char> names(32 * timesteps.size(), ' ');
  for (size_t state = 0; state < timesteps.size(); state++) {
    // This name is the "postfix" or common portion of all FlowSolution names...
    std::string name = "SolutionAtStep" + std::to_string(state + 1);
    std::strncpy(&names[state * 32], name.c_str(), 32);
    for (size_t i = name.size(); i < 32; i++) {
      names[state * 32 + i] = ' ';
    }
  }

  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.
  auto ziter = [=](Ioss::EntityBlock *block) {
    int              zone = get_db_zone(block);
    std::vector<int> indices(timesteps.size());
    bool             has_cell_center_fields = block->field_count(Ioss::Field::TRANSIENT) > 0;
    if (has_cell_center_fields || has_nodal_fields) {
      CGCHECK(cg_ziter_write(cgnsFilePtr, base, zone, "ZoneIterativeData"));
      CGCHECK(cg_goto(cgnsFilePtr, base, "Zone_t", zone, "ZoneIterativeData_t", 1, "end"));
      CGCHECK(cg_array_write("FlowSolutionPointers", CG_Character, 2, dim, names.data()));

      if (has_nodal_fields) {
        int index     = 1;
        int increment = has_cell_center_fields ? 2 : 1;
        for (size_t state = 0; state < timesteps.size(); state++) {
          indices[state] = index;
          index += increment;
        }

        CGCHECK(cg_array_write("VertexSolutionIndices", CG_Integer, 1, &dim[1], indices.data()));
        CGCHECK(cg_descriptor_write("VertexPrefix", "Vertex"));
      }
      if (has_cell_center_fields) {
        int index     = has_nodal_fields ? 2 : 1;
        int increment = has_nodal_fields ? 2 : 1;
        for (size_t state = 0; state < timesteps.size(); state++) {
          indices[state] = index;
          index += increment;
        }

        CGCHECK(cg_array_write("CellCenterIndices", CG_Integer, 1, &dim[1], indices.data()));
        CGCHECK(cg_descriptor_write("CellCenterPrefix", "CellCenter"));
      }
    }
  };

  // Use the lambda...
  const auto &sblocks = region->get_structured_blocks();
  for (auto &block : sblocks) {
    if (is_parallel_io || block->is_active()) {
      ziter(block);
    }
  }

  // Use the lambda...
  const auto &eblocks = region->get_element_blocks();
  for (auto &block : eblocks) {
    ziter(block);
  }
}

void Iocgns::Utils::add_transient_variables(int cgnsFilePtr, const std::vector<double> &timesteps,
                                            Ioss::Region *region, bool enable_field_recognition,
                                            char suffix_separator, int myProcessor,
                                            bool is_parallel_io)
{
  // ==========================================
  // Add transient variables (if any) to all zones...
  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.

  // Assuming that the fields on all steps are the same, but can vary
  // from zone to zone.
  auto sol_iter = [=](Ioss::EntityBlock *block) {
    int b = block->get_property("base").get_int();
    int z = get_db_zone(block);

    int sol_count = 0;
    CGCHECK(cg_nsols(cgnsFilePtr, b, z, &sol_count));
    int sol_per_step = sol_count / (int)timesteps.size();
    assert(sol_count % (int)timesteps.size() == 0);

    for (int sol = 1; sol <= sol_per_step; sol++) {
      char              solution_name[CGNS_MAX_NAME_LENGTH + 1];
      CG_GridLocation_t grid_loc;
      CGCHECK(cg_sol_info(cgnsFilePtr, b, z, sol, solution_name, &grid_loc));

      int field_count = 0;
      CGCHECK(cg_nfields(cgnsFilePtr, b, z, sol, &field_count));

      char **field_names = Ioss::Utils::get_name_array(field_count, CGNS_MAX_NAME_LENGTH + 1);
      for (int field = 1; field <= field_count; field++) {
        CG_DataType_t data_type;
        char          field_name[CGNS_MAX_NAME_LENGTH + 1];
        CGCHECK(cg_field_info(cgnsFilePtr, b, z, sol, field, &data_type, field_name));
        std::strncpy(field_names[field - 1], field_name, CGNS_MAX_NAME_LENGTH);
      }

      // Convert raw field names into composite fields (a_x, a_y, a_z ==> 3D vector 'a')
      std::vector<Ioss::Field> fields;
      if (grid_loc == CG_CellCenter) {
        size_t entity_count = block->entity_count();
        Ioss::Utils::get_fields(entity_count, field_names, field_count, Ioss::Field::TRANSIENT,
                                enable_field_recognition, suffix_separator, nullptr, fields);
        size_t index = 1;
        for (const auto &field : fields) {
          Utils::set_field_index(field, index, grid_loc);
          index += field.raw_storage()->component_count();
          block->field_add(field);
        }
      }
      else {
        assert(grid_loc == CG_Vertex);
        const Ioss::NodeBlock *cnb =
            (block->type() == Ioss::STRUCTUREDBLOCK)
                ? &(dynamic_cast<Ioss::StructuredBlock *>(block)->get_node_block())
                : region->get_node_blocks()[0];
        Ioss::NodeBlock *nb           = const_cast<Ioss::NodeBlock *>(cnb);
        size_t           entity_count = nb->entity_count();
        Ioss::Utils::get_fields(entity_count, field_names, field_count, Ioss::Field::TRANSIENT,
                                enable_field_recognition, suffix_separator, nullptr, fields);
        size_t index = 1;
        for (const auto &field : fields) {
          Utils::set_field_index(field, index, grid_loc);
          index += field.raw_storage()->component_count();
          nb->field_add(field);
        }
      }

      Ioss::Utils::delete_name_array(field_names, field_count);
    }
  };
  // ==========================================

  if (!timesteps.empty()) {
    const auto &sblocks = region->get_structured_blocks();
    for (auto &block : sblocks) {
      if (is_parallel_io || block->is_active()) {
        sol_iter(block);
      }
    }
    const auto &eblocks = region->get_element_blocks();
    for (auto &block : eblocks) {
      sol_iter(block);
    }
    bool is_parallel = region->get_database()->util().parallel_size() > 1;
    if (is_parallel && !is_parallel_io) {
      sync_transient_variables_fpp(region, myProcessor);
    }
  }
}

int Iocgns::Utils::get_step_times(int cgnsFilePtr, std::vector<double> &timesteps,
                                  Ioss::Region *region, double timeScaleFactor, int myProcessor)
{
  int  base          = 1;
  int  num_timesteps = 0;
  char bitername[CGNS_MAX_NAME_LENGTH + 1];
  int  ierr = cg_biter_read(cgnsFilePtr, base, bitername, &num_timesteps);
  if (ierr == CG_NODE_NOT_FOUND) {
    return num_timesteps;
  }
  if (ierr == CG_ERROR) {
    Utils::cgns_error(cgnsFilePtr, __FILE__, __func__, __LINE__, myProcessor);
  }

  if (num_timesteps <= 0) {
    return num_timesteps;
  }

  // Read the timestep time values.
  CGCHECK(cg_goto(cgnsFilePtr, base, "BaseIterativeData_t", 1, "end"));
  std::vector<double> times(num_timesteps);
  CGCHECK(cg_array_read_as(1, CG_RealDouble, times.data()));

  timesteps.reserve(num_timesteps);
  for (int i = 0; i < num_timesteps; i++) {
    region->add_state(times[i] * timeScaleFactor);
    timesteps.push_back(times[i]);
  }
  return num_timesteps;
}

void Iocgns::Utils::assign_zones_to_procs(std::vector<Iocgns::StructuredZoneData *> &all_zones,
                                          std::vector<size_t> &                      work_vector)
{
  for (auto &zone : all_zones) {
    zone->m_proc = -1;
  }

  // Sort zones based on work.  Most work first.. Filtered to active only...
  std::vector<Iocgns::StructuredZoneData *> zones;
  std::copy_if(all_zones.begin(), all_zones.end(), std::back_inserter(zones),
               [](Iocgns::StructuredZoneData *z) { return z->is_active(); });

  std::sort(zones.begin(), zones.end(),
            [](Iocgns::StructuredZoneData *a, Iocgns::StructuredZoneData *b) {
              return a->work() > b->work();
            });

  std::set<std::pair<int, int>> proc_adam_map;

  // On first entry, work_vector will be all zeros.  To avoid any
  // searching, assign the first `nproc` zones to the `nproc` entries
  // in `work_vector`.  Avoids searching...
  assert(zones.size() >= work_vector.size());
  size_t i = 0;
  for (; i < work_vector.size(); i++) {
    auto &zone   = zones[i];
    zone->m_proc = i;
    work_vector[i] += zone->work();
    proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, zone->m_proc));
  }

  for (; i < zones.size(); i++) {
    auto &zone = zones[i];

    // Assign zone to processor with minimum work that does not already have a zone with the same
    // adam zone...
    ssize_t proc = proc_with_minimum_work(zone, work_vector, proc_adam_map);

    // See if any other zone on this processor has the same adam zone...
    if (proc >= 0) {
      auto success = proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, proc));
      if (success.second) {
        zone->m_proc = proc;
        work_vector[proc] += zone->work();
      }
      else {
        std::ostringstream errmsg;
        errmsg << "IOCGNS error: Could not assign zones to processors in " << __func__;
        IOSS_ERROR(errmsg);
      }
    }
    else {
      std::ostringstream errmsg;
      errmsg << "IOCGNS error: Could not assign zones to processors in " << __func__;
      IOSS_ERROR(errmsg);
    }
  }
}

size_t Iocgns::Utils::pre_split(std::vector<Iocgns::StructuredZoneData *> &zones, double avg_work,
                                double load_balance, int proc_rank, int proc_count)
{
  auto   new_zones(zones);
  size_t new_zone_id = zones.size() + 1;

  // See if can split each zone over a set of procs...
  double           total_work = 0.0;
  std::vector<int> splits(zones.size());

  for (size_t i = 0; i < zones.size(); i++) {
    auto   zone = zones[i];
    double work = zone->work();
    total_work += work;
    splits[i] = int(std::round(work / avg_work));
    splits[i] = splits[i] == 0 ? 1 : splits[i];
  }

  int num_splits = std::accumulate(splits.begin(), splits.end(), 0);
  int diff       = proc_count - num_splits;

  while (diff != 0) {
    // Adjust splits so sum is equal to proc_count.
    // Adjust the largest split count(s)
    int    step      = diff < 0 ? -1 : 1;
    size_t min_z     = 0;
    double min_delta = 1.0e27;
    for (size_t i = 0; i < zones.size(); i++) {
      auto   zone = zones[i];
      double work = zone->work();

      if (splits[i] == 0) {
        continue;
      }
      double zone_avg = work / (double)splits[i];
      if ((splits[i] + step) > 0) {
        double delta = std::abs(zone_avg - work / (double)(splits[i] + step));
        if (delta < min_delta) {
          min_delta = delta;
          min_z     = i;
        }
      }
    }
    splits[min_z] += step;
    diff -= step;
  }

  assert(diff == 0);
  assert(std::accumulate(splits.begin(), splits.end(), 0) == proc_count);

  // See if splits result in avg_work for all zones in range...
  double min_avg      = avg_work / load_balance;
  double max_avg      = avg_work * load_balance;
  bool   adaptive_avg = true;
  for (size_t i = 0; i < zones.size(); i++) {
    auto   zone = zones[i];
    double work = zone->work();
    if (splits[i] == 0) {
      adaptive_avg = false;
      break;
    }
    double zone_avg = work / (double)splits[i];
    if (zone_avg < min_avg || zone_avg > max_avg) {
      adaptive_avg = false;
      break;
    }
  }

  if (adaptive_avg) {
    for (size_t i = 0; i < zones.size(); i++) {
      auto zone       = zones[i];
      int  num_active = 0;

      auto work_average = avg_work;
      int  split_cnt    = splits[i];
      if (split_cnt == 1) {
        continue;
      }

      work_average = zone->work() / (double)split_cnt;

      std::vector<std::pair<int, Iocgns::StructuredZoneData *>> active;
      active.push_back(std::make_pair(split_cnt, zone));
      do {
        assert(!active.empty());
        split_cnt = active.back().first;
        zone      = active.back().second;
        active.pop_back();

        if (zone->is_active()) {
          if (split_cnt != 1) {
            int max_power_2 = power_2(split_cnt);
            if (max_power_2 == split_cnt) {
              work_average = zone->work() / 2.0;
              max_power_2 /= 2;
            }
            else {
              work_average = zone->work() / (double(split_cnt) / double(max_power_2));
            }

            auto children = zone->split(new_zone_id, work_average, load_balance, proc_rank);
            if (children.first != nullptr && children.second != nullptr) {
              new_zones.push_back(children.first);
              new_zones.push_back(children.second);
              new_zone_id += 2;
              active.push_back(std::make_pair(split_cnt - max_power_2, children.second));
              active.push_back(std::make_pair(max_power_2, children.first));
              num_active++;
            }
          }
        }
        if (num_active >=
            proc_count) { // Don't split a single zone into more than `proc_count` pieces
          break;
        }
      } while (!active.empty());
    }
  }
  else {
    for (size_t i = 0; i < zones.size(); i++) {
      auto zone       = zones[i];
      int  num_active = 0;
      if (zone->work() <= max_avg) {
        // This zone is already in `new_zones`; just skip doing anything else with it.
      }
      else {
        std::vector<std::pair<int, Iocgns::StructuredZoneData *>> active;

        double work      = zone->work();
        int    split_cnt = int(work / avg_work);

        // Find modulus of work % avg_work and split off that amount
        // which will be < avg_work.
        double mod_work = work - avg_work * split_cnt;
        if (mod_work > max_avg - avg_work) {
          auto children = zone->split(new_zone_id, mod_work, load_balance, proc_rank);
          if (children.first != nullptr && children.second != nullptr) {
            new_zones.push_back(children.first);
            new_zones.push_back(children.second);
            new_zone_id += 2;
            num_active++;
            active.push_back(std::make_pair(split_cnt, children.second));
          }
          else {
            active.push_back(std::make_pair(split_cnt, zone));
          }
        }
        else {
          active.push_back(std::make_pair(split_cnt, zone));
        }

        // The work remaining on this zone should be approximately
        // equally divided among `split_cnt` processors.
        do {
          assert(!active.empty());
          split_cnt = active.back().first;
          zone      = active.back().second;
          active.pop_back();

          if (zone->is_active()) {
            int    max_power_2  = power_2(split_cnt);
            double work_average = 0.0;
            if (max_power_2 == split_cnt) {
              work_average = zone->work() / 2.0;
            }
            else {
              work_average = zone->work() / (double(split_cnt) / double(max_power_2));
            }

            if (max_power_2 != 1) {
              if (max_power_2 == split_cnt) {
                max_power_2 /= 2;
              }
              auto children = zone->split(new_zone_id, work_average, load_balance, proc_rank);
              if (children.first != nullptr && children.second != nullptr) {
                new_zones.push_back(children.first);
                new_zones.push_back(children.second);
                new_zone_id += 2;
                active.push_back(std::make_pair(split_cnt - max_power_2, children.second));
                active.push_back(std::make_pair(max_power_2, children.first));
                num_active++;
              }
            }
          }
          if (num_active >=
              proc_count) { // Don't split a single zone into more than `proc_count` pieces
            break;
          }
        } while (!active.empty());
      }
    }
  }
  std::swap(new_zones, zones);
  if (zones.size() < (size_t)proc_count && load_balance > 1.05) {
    // Tighten up the load_balance factor to get some decomposition going...
    double new_load_balance = (1.0 + load_balance) / 2.0;
    new_zone_id             = pre_split(zones, avg_work, new_load_balance, proc_rank, proc_count);
  }
  return new_zone_id;
}
