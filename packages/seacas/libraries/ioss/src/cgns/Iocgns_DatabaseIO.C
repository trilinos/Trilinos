// CGNS Assumptions:
// * All boundary conditions are listed as Family nodes at the "top" level.
// * Single element block per zone.
// * Single Base.
// * ZoneGridConnectivity is 1to1 with point lists for unstructured

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
#include <Ioss_Utils.h>
#include <bitset>
#include <cassert>
#include <cgns/Iocgns_DatabaseIO.h>
#include <cgns/Iocgns_Utils.h>
#include <cgnslib.h>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <sys/select.h>
#include <tokenize.h>
#include <vector>

#if !defined(CGNSLIB_H)
#error "Could not include cgnslib.h"
#endif

#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_Hex8.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Quad4.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_State.h"
#include "Ioss_StructuredBlock.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"

extern char hdf5_access[64];

#ifdef SEACAS_HAVE_MPI
extern int pcg_mpi_initialized;
#endif

namespace {

  template <typename T> void pack(int &idx, std::vector<int> &pack, T *from, int count)
  {
    for (int i = 0; i < count; i++) {
      pack[idx++] = from[i];
    }
  }

  template <typename T> void unpack(int &idx, T *pack, T *to, int count)
  {
    for (int i = 0; i < count; i++) {
      to[i] = pack[idx++];
    }
  }

  struct SBlock
  {
    SBlock() = default;
    SBlock(char *names, int *data)
    {
      name    = std::string{names};
      int idx = 0;
      proc    = data[idx++];
      unpack(idx, data, range.data(), 3);
      local_zone = data[idx++];
    }
    std::string                      name{};
    int                              proc{-1};
    int                              local_zone{};
    std::vector<std::pair<int, int>> adjacency; // face, proc pairs
    std::array<int, 3>               range{{0, 0, 0}};
    std::array<int, 3>               glob_range{{0, 0, 0}};
    std::array<int, 3>               offset{{0, 0, 0}};
    std::bitset<6>                   face_adj;

    bool split() const { return face_adj.any(); }
  };

  std::pair<std::string, int> decompose_name(const std::string &name, bool is_parallel)
  {
    int         proc = 0;
    std::string zname{name};

    if (is_parallel) {
      // Name should/might be of the form `basename_proc-#`.  Strip
      // off the `_proc-#` portion and return just the basename.
      auto tokens = Ioss::tokenize(zname, "_");
      zname       = tokens[0];
      if (tokens.size() >= 2) {
        size_t idx = tokens.size() - 1;
        if (tokens[idx].substr(0, 5) == "proc-") {
          auto ptoken = Ioss::tokenize(tokens[idx], "-");
          proc        = std::stoi(ptoken[1]);
          idx--;
          zname = tokens[idx];
        }
      }
    }
    return std::make_pair(zname, proc);
  }

#ifdef SEACAS_HAVE_MPI
  void add_zgc_fpp(int cgnsFilePtr, Ioss::StructuredBlock *block,
                   const std::map<std::string, int> &zone_name_map, int myProcessor,
                   bool isParallel)
  {
    int base    = block->get_property("base").get_int();
    int zone    = block->get_property("zone").get_int();
    int db_zone = Iocgns::Utils::get_db_zone(block);
    int nconn   = 0;
    CGCHECK(cg_n1to1(cgnsFilePtr, base, db_zone, &nconn));

    for (int ii = 0; ii < nconn; ii++) {
      char                    connectname[CGNS_MAX_NAME_LENGTH + 1];
      char                    donorname[CGNS_MAX_NAME_LENGTH + 1];
      std::array<cgsize_t, 6> range;
      std::array<cgsize_t, 6> donor_range;
      Ioss::IJK_t             transform;

      CGCHECK(cg_1to1_read(cgnsFilePtr, base, db_zone, ii + 1, connectname, donorname, range.data(),
                           donor_range.data(), transform.data()));

      auto        donorname_proc = decompose_name(donorname, isParallel);
      std::string donor_name     = donorname_proc.first;

      // Get number of nodes shared with other "previous" zones...
      // A "previous" zone will have a lower zone number this this zone...
      int  donor_zone = -1;
      auto donor_iter = zone_name_map.find(donor_name);
      if (donor_iter != zone_name_map.end()) {
        donor_zone = (*donor_iter).second;
      }
      Ioss::IJK_t range_beg{{(int)range[0], (int)range[1], (int)range[2]}};
      Ioss::IJK_t range_end{{(int)range[3], (int)range[4], (int)range[5]}};
      Ioss::IJK_t donor_beg{{(int)donor_range[0], (int)donor_range[1], (int)donor_range[2]}};
      Ioss::IJK_t donor_end{{(int)donor_range[3], (int)donor_range[4], (int)donor_range[5]}};

      Ioss::IJK_t offset;
      offset[0] = block->get_property("offset_i").get_int();
      offset[1] = block->get_property("offset_j").get_int();
      offset[2] = block->get_property("offset_k").get_int();
      range_beg[0] += offset[0];
      range_beg[1] += offset[1];
      range_beg[2] += offset[2];
      range_end[0] += offset[0];
      range_end[1] += offset[1];
      range_end[2] += offset[2];

      auto con_name = decompose_name(connectname, isParallel).first;
      block->m_zoneConnectivity.emplace_back(con_name, zone, donor_name, donor_zone, transform,
                                             range_beg, range_end, donor_beg, donor_end, offset);

      block->m_zoneConnectivity.back().m_ownerProcessor = myProcessor;
      block->m_zoneConnectivity.back().m_donorProcessor = donorname_proc.second;
    }
  }
#endif

#ifdef SEACAS_HAVE_MPI
  int adjacent_block(const SBlock &b, int ijk, std::map<int, int> &proc_block_map)
  {
    // Find a block the the 'left|right|up|down|front|back' (ijk) of blocks[br]
    if (b.face_adj[ijk]) {
      for (auto adj : b.adjacency) {
        if (adj.first == ijk) {
          int proc = adj.second;
          return proc_block_map[proc];
        }
      }
    }
    return -1;
  }

  void set_block_offset(size_t begin, size_t end, std::vector<SBlock> &blocks,
                        std::map<int, int> &proc_block_map)
  {
    for (size_t p = 0; p < (end - begin); p++) {
      for (size_t j = begin; j < end; j++) {
        auto &block = blocks[j];
        // See which blocks are below/left/under this block which means
        // that this blocks offset is affected.
        for (int ijk = 0; ijk < 3; ijk++) {
          int br = adjacent_block(block, ijk, proc_block_map);
          if (br >= 0) {
            block.offset[ijk] = blocks[br].offset[ijk] + blocks[br].range[ijk];
          }
        }
      }
    }
  }

  void set_global_extent(size_t begin, size_t end, std::vector<SBlock> &blocks,
                         std::map<int, int> &proc_block_map)
  {
    // Determine the global ijk extent for the block which is spread over multiple processors
    // and is in the range [begin, end) in blocks.
    Ioss::IJK_t global{0, 0, 0};
    for (int ijk = 0; ijk < 3; ijk++) {
      // Find a block in range [bbeg, bend) with no block to the "left|below|behind
      for (size_t bb = begin; bb < end; bb++) {
        if (blocks[bb].face_adj[ijk] == 0) {
          // No blocks to min 'ijk' direction...
          // Traverse all blocks toward max 'ijk' direction setting offsets and global range.
          size_t iter = 0;
          int    br   = bb;
          do {
            global[ijk] += blocks[br].range[ijk];
            br = adjacent_block(blocks[br], ijk + 3, proc_block_map);
            if (++iter > end - begin) {
              auto               bp = adjacent_block(blocks[br], ijk + 3, proc_block_map);
              std::ostringstream errmsg;
              errmsg << "ERROR: CGNS: Block '" << blocks[bb].name
                     << "' is in infinite loop calculating processor adjacencies for direction "
                     << (ijk == 0 ? 'i' : ijk == 1 ? 'j' : 'k') << " on processors "
                     << blocks[bp].proc << " and " << blocks[br].proc << ".  Check decomposition.";
              IOSS_ERROR(errmsg);
            }
          } while (br >= 0);
          break;
        }
      }
    }
    for (size_t bb = begin; bb < end; bb++) {
      blocks[bb].glob_range = global;
    }
  }

  int find_face(const std::array<cgsize_t, 6> &range)
  {
    // 0,1,2 == min x,y,z; 3,4,5 == Max x,y,z
    bool is_x = range[0] == range[3];
    bool is_y = range[1] == range[4];
#ifndef NDEBUG
    bool is_z = range[2] == range[5];
#endif
    assert(is_x || is_y || is_z);
    assert((is_x ? 1 : 0) + (is_y ? 1 : 0) + (is_z ? 1 : 0) == 1);
    int idx = is_x ? 0 : is_y ? 1 : 2;

    // Which face on this block?
    int face = idx;
    if (range[idx] != 1) {
      face += 3;
    }
    return face;
  }
#endif

#ifdef SEACAS_HAVE_MPI
  bool generate_inter_proc_adjacency(int cgnsFilePtr, int base, int zone, int myProcessor,
                                     const std::string &zone_name, std::vector<int> &adjacency)
  {
    // Handle zone-grid-connectivity... At this point we only want
    // the zgc that are inter-proc between the same "base zone".
    // That is, the zgc which are result of parallel decomp.

    // Stored in format:  "-myproc, -local_zone, face, shared_proc" for each shared face.
    bool zone_added = false;
    int  nconn      = 0;
    CGCHECK(cg_n1to1(cgnsFilePtr, base, zone, &nconn));
    for (int i = 0; i < nconn; i++) {
      char                    connectname[CGNS_MAX_NAME_LENGTH + 1];
      char                    donorname[CGNS_MAX_NAME_LENGTH + 1];
      std::array<cgsize_t, 6> range;
      std::array<cgsize_t, 6> donor_range;
      Ioss::IJK_t             transform;

      CGCHECK(cg_1to1_read(cgnsFilePtr, base, zone, i + 1, connectname, donorname, range.data(),
                           donor_range.data(), transform.data()));

      auto        donorname_proc = decompose_name(donorname, true);
      std::string donor_name     = donorname_proc.first;

      if (donor_name == zone_name) {
        // Determine which face of the zone on this processor is
        // shared with the other processor...
        int face = find_face(range);
        adjacency.push_back(-myProcessor);
        adjacency.push_back(-zone);
        adjacency.push_back(face);
        adjacency.push_back(donorname_proc.second);
        zone_added = true;
      }
    }
    return zone_added;
  }

  void set_adjacency(SBlock &b, std::vector<int> &adjacency)
  {
    // Stored in format:  "-myproc, -local_zone, face, shared_proc" for each shared face.
    for (size_t i = 0; i < adjacency.size(); i += 4) {
      assert(adjacency[i] <= 0); // -proc
      if (adjacency[i] == -b.proc) {
        assert(adjacency[i + 1] < 0);
        if (adjacency[i + 1] == -b.local_zone) {
          b.adjacency.emplace_back(adjacency[i + 2], adjacency[i + 3]);
          b.face_adj.set(adjacency[i + 2]);
        }
      }
      else if (adjacency[i] < -b.proc) {
        return;
      }
    }
  }

#endif

#ifdef SEACAS_HAVE_MPI
  void add_empty_bc(Ioss::SideSet *sset, Ioss::StructuredBlock *block, int base, int zone, int face,
                    const std::string &fam_name, const std::string &boco_name)
  {
    assert(sset != nullptr);

    Ioss::IJK_t empty_range{{0, 0, 0}};

    auto sbc   = Ioss::BoundaryCondition(boco_name, fam_name, empty_range, empty_range);
    sbc.m_face = face;
    block->m_boundaryConditions.push_back(sbc);

    std::string name = boco_name + "/" + block->name();

    auto sb =
        new Ioss::SideBlock(block->get_database(), name, Ioss::Quad4::name, Ioss::Hex8::name, 0);
    sb->set_parent_block(block);
    sset->add(sb);
    sb->property_add(Ioss::Property("base", base));
    sb->property_add(Ioss::Property("zone", zone));
    sb->property_add(Ioss::Property("section", face + 1));
    sb->property_add(Ioss::Property("id", sset->get_property("id").get_int()));
    sb->property_add(Ioss::Property(
        "guid", block->get_database()->util().generate_guid(sset->get_property("id").get_int())));
  }
#endif

} // namespace

namespace Iocgns {

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props)
  {
    dbState = Ioss::STATE_UNKNOWN;

#if IOSS_DEBUG_OUTPUT
    if (myProcessor == 0) {
      std::cout << "CGNS DatabaseIO using " << CG_SIZEOF_SIZE << "-bit integers.\n";
    }
#endif
    openDatabase__();
  }

  DatabaseIO::~DatabaseIO()
  {
    for (auto &gtb : m_globalToBlockLocalNodeMap) {
      delete gtb.second;
    }
    closeDatabase__();
  }

  int DatabaseIO::get_file_pointer() const 
  {
    if (cgnsFilePtr < 0) {
      openDatabase__();
    }
    return cgnsFilePtr;
  }

  void DatabaseIO::openDatabase__() const
  {
    if (cgnsFilePtr < 0) {
      if ((is_input() && properties.exists("MEMORY_READ")) ||
          (!is_input() && properties.exists("MEMORY_WRITE"))) {
        strcpy(hdf5_access, "PARALLEL");
      }
      int mode = is_input() ? CG_MODE_READ : CG_MODE_WRITE;
      CGCHECK(cg_set_file_type(CG_FILE_HDF5));

#ifdef SEACAS_HAVE_MPI
      // Kluge to get fpp and dof CGNS working at same time.
      pcg_mpi_initialized = 0;
#endif
      int ierr = cg_open(decoded_filename().c_str(), mode, &cgnsFilePtr);
      // Will not return if error...
      check_valid_file_open(ierr);
      if ((is_input() && properties.exists("MEMORY_READ")) ||
          (!is_input() && properties.exists("MEMORY_WRITE"))) {
        strcpy(hdf5_access, "NATIVE");
      }

      if (properties.exists("INTEGER_SIZE_API")) {
        int isize = properties.get("INTEGER_SIZE_API").get_int();
        if (isize == 8) {
          set_int_byte_size_api(Ioss::USE_INT64_API);
        }
        if (isize == 4) {
          set_int_byte_size_api(Ioss::USE_INT32_API);
        }
      }
      else if (CG_SIZEOF_SIZE == 64) {
        set_int_byte_size_api(Ioss::USE_INT64_API);
      }

#if 0
      // This isn't currently working since CGNS currently has chunking
      // disabled for HDF5 files and compression requires chunking.
      if (!is_input()) {
	if (properties.exists("COMPRESSION_LEVEL")) {
	  int comp = properties.get("COMPRESSION_LEVEL").get_int();
	  cg_configure(CG_CONFIG_HDF5_COMPRESS, (void*)comp);
	}
      }
#endif
    }
    assert(cgnsFilePtr >= 0);
  }

  void DatabaseIO::closeDatabase__() const
  {
    if (cgnsFilePtr != -1) {
      CGCHECK(cg_close(cgnsFilePtr));
    }
    cgnsFilePtr = -1;
  }

  bool DatabaseIO::check_valid_file_open(int status) const
  {
    int global_status = status;
    if (isParallel) {
      global_status = util().global_minmax(status, Ioss::ParallelUtils::DO_MAX);
    }

    if (global_status != CG_OK) {
      Ioss::IntVector err_status;
      if (isParallel) {
        util().all_gather(status, err_status);
      }
      else {
        err_status.push_back(status);
      }

      // See which processors could not open/create the file...
      std::ostringstream errmsg;
      int                ok_count = 0;
      if (isParallel) {
        ok_count = std::count(err_status.begin(), err_status.end(), CG_OK);
        if (ok_count == 0 && util().parallel_size() > 2) {
          errmsg << "ERROR: Unable to open CGNS decomposed database files:\n\t\t"
                 << Ioss::Utils::decode_filename(get_filename(), 0, util().parallel_size())
                 << "  ...\n\t\t"
                 << Ioss::Utils::decode_filename(get_filename(), util().parallel_size() - 1,
                                                 util().parallel_size())
                 << "\n";
        }
        else {
          errmsg << "ERROR: Unable to open CGNS decomposed database files:\n";
          for (int i = 0; i < util().parallel_size(); i++) {
            if (err_status[i] != CG_OK) {
              errmsg << "\t\t"
                     << Ioss::Utils::decode_filename(get_filename(), i, util().parallel_size())
                     << "\n";
            }
          }
        }
        errmsg << "       for " << (is_input() ? "read" : "write") << " access.\n";
      }
      else {
        errmsg << "ERROR: Unable to open CGNS database '" << get_filename() << "' for "
               << (is_input() ? "read" : "write") << " access.\n";
      }
      if (status != CG_OK) {
        if (ok_count != 0 || util().parallel_size() <= 2) {
          std::ostringstream msg;
          msg << "[" << myProcessor << "] CGNS Error: '" << cg_get_error() << "'\n";
          std::cerr << msg.str();
        }
        else {
          // Since error on all processors, assume the same error on all and only print
          // the error from processor 0.
          if (myProcessor == 0) {
            std::ostringstream msg;
            msg << "CGNS Error: '" << cg_get_error() << "'\n";
            std::cerr << msg.str();
          }
        }
      }

      IOSS_ERROR(errmsg);
      return false;
    }
    return true;
  }

  void DatabaseIO::finalize_database()
  {
    if (is_input()) {
      return;
    }

    if (m_timesteps.empty()) {
      return;
    }

    Utils::finalize_database(get_file_pointer(), m_timesteps, get_region(), myProcessor, false);
  }

  int64_t DatabaseIO::node_global_to_local__(int64_t global, bool /*must_exist*/) const
  {
    return global;
  }

  int64_t DatabaseIO::element_global_to_local__(int64_t global) const { return global; }

  void DatabaseIO::create_structured_block_fpp(int base, int num_zones, size_t &num_node)
  {
    assert(isParallel);
#ifdef SEACAS_HAVE_MPI
    // Each processor may have a different set of zones.  This routine
    // will sync the information such that at return, each procesosr
    // has a consistent set of structuredBlocks defined with the
    // correct local and global, i,j,k ranges and offsets.
    // First each processor sends their zone count to processor 0...

    // name, proc (int) , cell-range (3 int), boundary-with

    // First, get basenames of all zones on all processors so we can
    // work with consistent set...
    int                        id               = 0;
    int                        in               = 0;
    const int                  INT_PER_ZONE     = 5;  // proc, range[3], zone
    const int                  OUT_INT_PER_ZONE = 10; // proc, range[3], glob_range[3], offset[3]
    std::vector<int>           zone_data(num_zones * INT_PER_ZONE);
    std::vector<char>          zone_names(num_zones * (CGNS_MAX_NAME_LENGTH + 1));
    std::map<std::string, int> zone_id_map;
    std::vector<int>           adjacency;

    for (int zone = 1; zone <= num_zones; zone++) {
      cgsize_t size[9];
      char     zname[CGNS_MAX_NAME_LENGTH + 1];
      CGCHECK(cg_zone_read(get_file_pointer(), base, zone, zname, size));

      assert(size[0] - 1 == size[3]);
      assert(size[1] - 1 == size[4]);
      assert(size[2] - 1 == size[5]);

      assert(size[6] == 0);
      assert(size[7] == 0);
      assert(size[8] == 0);

      auto        name_proc = decompose_name(zname, isParallel);
      std::string zone_name = name_proc.first;
      int         proc      = name_proc.second;
      assert(proc == myProcessor);

      zone_data[id++] = proc;
      pack(id, zone_data, &size[3], 3); // Packing 3,4,5
      strncpy(&zone_names[in], zone_name.c_str(), CGNS_MAX_NAME_LENGTH);
      in += CGNS_MAX_NAME_LENGTH + 1;
      zone_id_map[zone_name] = zone;

      // Handle zone-grid-connectivity... At this point we only want
      // the zgc that are inter-proc between the same "base zone".
      // That is, the zgc which are result of parallel decomp.

      // Stored as -P, -Z, f1, p1, -P, -Z, f2, p2, ..., -P, -Z, f1, ...
      generate_inter_proc_adjacency(get_file_pointer(), base, zone, myProcessor, zone_name, adjacency);

      zone_data[id++] = zone;
      assert(id % INT_PER_ZONE == 0);
    }

    // Now gather all information to processor 0
    std::vector<char> all_names;
    std::vector<int>  all_data;
    std::vector<int>  all_adj;
    util().gather(num_zones, CGNS_MAX_NAME_LENGTH + 1, zone_names, all_names);
    int tot_zones = util().gather(num_zones, INT_PER_ZONE, zone_data, all_data);
    util().gather((int)adjacency.size(), 1, adjacency, all_adj);

    if (myProcessor == 0) {
      std::vector<SBlock> blocks;
      int                 off_name = 0;
      int                 off_data = 0;
      for (int i = 0; i < tot_zones; i++) {
        blocks.emplace_back(&all_names[off_name], &all_data[off_data]);
        off_name += CGNS_MAX_NAME_LENGTH + 1;
        off_data += INT_PER_ZONE;

        // Add inter-processor adjacency information to the block
        auto &b = blocks.back();
        set_adjacency(b, all_adj);
      }
      all_adj.clear();
      all_adj.shrink_to_fit();

      // Sort blocks to get similar zones adjacent -- will have same name, but different proc
      std::sort(blocks.begin(), blocks.end(), [](const SBlock &b1, const SBlock &b2) {
        return (b1.name == b2.name ? b1.proc < b2.proc : b1.name < b2.name);
      });

      int                 proc_count = util().parallel_size();
      std::vector<SBlock> resolved_blocks;

      for (size_t i = 0; i < blocks.size(); i++) {
        auto &b = blocks[i];
        if (b.split()) {
          // The blocks it is split with should be adjacent in list.
          // Get range of indices referring to this block and build
          // a map from processor to index, so build that now...
          std::map<int, int> proc_block_map;
          proc_block_map[b.proc] = i;
          size_t j               = i + 1;
          for (; j < blocks.size(); j++) {
            if (blocks[j].name != b.name) {
              break;
            }
            proc_block_map[blocks[j].proc] = j;
          }
          auto bbeg = i;
          auto bend = j;

          // Get global ijk extent in each direction...
          set_global_extent(bbeg, bend, blocks, proc_block_map);

          // Iterate to get correct offset for these blocks on all processors...
          set_block_offset(bbeg, bend, blocks, proc_block_map);

#if IOSS_DEBUG_OUTPUT
          std::cerr << "Range of blocks for " << b.name << " is " << i << " to " << j - 1
                    << " Global I,J,K = " << b.glob_range[0] << " " << b.glob_range[1] << " "
                    << b.glob_range[2] << "\n";
#endif
          // All processors need to know about it...
          for (int p = 0; p < proc_count; p++) {
            auto iter = proc_block_map.find(p);
            if (iter == proc_block_map.end()) {
              SBlock newb;
              newb.name       = b.name;
              newb.proc       = p;
              newb.glob_range = b.glob_range;
              resolved_blocks.push_back(newb);
            }
            else {
              auto   idx  = (*iter).second;
              SBlock newb = blocks[idx];
              resolved_blocks.push_back(newb);
            }
          }
          i = bend - 1;
        }
        else {
          // If not split, then global size = local size and offset = 0
          b.glob_range = b.range;

          // All processors need to know about it...
          for (int p = 0; p < proc_count; p++) {
            SBlock newb;
            newb.name       = b.name;
            newb.proc       = p;
            newb.glob_range = b.glob_range;
            if (p == b.proc) {
              newb.range = b.range;
            }
            resolved_blocks.push_back(newb);
          }
        }
      }

      int num_unique = (int)resolved_blocks.size() / proc_count;

#if IOSS_DEBUG_OUTPUT
      for (const auto &b : resolved_blocks) {
        std::cerr << b.name << " " << b.proc << " " << b.local_zone << " "
                  << " (" << b.range[0] << " " << b.range[1] << " " << b.range[2] << ") ("
                  << b.glob_range[0] << " " << b.glob_range[1] << " " << b.glob_range[2] << ") ("
                  << b.offset[0] << " " << b.offset[1] << " " << b.offset[2] << ") [" << b.face_adj
                  << "]"
                  << "\n";
      }
#endif

      // Data now consistent for all zones.  Send back to their "owning" processor
      tot_zones = num_unique;
      all_names.resize(num_unique * (CGNS_MAX_NAME_LENGTH + 1));
      all_data.resize(resolved_blocks.size() * OUT_INT_PER_ZONE);
      id = in = 0;
      for (int off = 0; off < proc_count; off++) {
        for (int b = 0; b < num_unique; b++) {
          int         idx   = off + b * proc_count;
          const auto &block = resolved_blocks[idx];
          if (off == 0) {
            strncpy(&all_names[in], block.name.c_str(), CGNS_MAX_NAME_LENGTH);
            in += CGNS_MAX_NAME_LENGTH + 1;
          }
          all_data[id++] = block.proc;
          pack(id, all_data, block.range.data(), 3);
          pack(id, all_data, block.glob_range.data(), 3);
          pack(id, all_data, block.offset.data(), 3);
        }
      }
      assert(id % OUT_INT_PER_ZONE == 0);
    }
    MPI_Bcast(&tot_zones, 1, MPI_INT, 0, util().communicator());
    zone_data.resize(tot_zones * OUT_INT_PER_ZONE);
    all_names.resize(tot_zones * (CGNS_MAX_NAME_LENGTH + 1));
    MPI_Bcast(all_names.data(), tot_zones * (CGNS_MAX_NAME_LENGTH + 1), MPI_CHAR, 0,
              util().communicator());
    MPI_Scatter(all_data.data(), tot_zones * OUT_INT_PER_ZONE, MPI_INT, zone_data.data(),
                tot_zones * OUT_INT_PER_ZONE, MPI_INT, 0, util().communicator());

    // Each processor now has a consistent set of structured blocks.
    // Create the Ioss::StructuredBlocks objects and add to region.
    id = in = 0;
    for (int i = 0; i < tot_zones; i++) {
      std::string zone_name(&all_names[in]);
      in += CGNS_MAX_NAME_LENGTH + 1;
      Ioss::IJK_t local_ijk;
      Ioss::IJK_t global_ijk;
      Ioss::IJK_t offset_ijk;

      zone_data[id++]; // proc field. Not currently used.
      unpack(id, zone_data.data(), local_ijk.data(), 3);
      unpack(id, zone_data.data(), global_ijk.data(), 3);
      unpack(id, zone_data.data(), offset_ijk.data(), 3);

      Ioss::StructuredBlock *block =
          new Ioss::StructuredBlock(this, zone_name, 3, local_ijk, offset_ijk, global_ijk);

      // See if this zone exists on this processor's file, or is just for interprocessor
      // consistency.
      int  zone   = tot_zones + i;
      bool native = false;
      auto iter   = zone_id_map.find(zone_name);
      if (iter != zone_id_map.end()) {
        zone   = (*iter).second;
        native = true;
      }

      block->property_add(Ioss::Property("base", base));
      if (native) {
        block->property_add(Ioss::Property("db_zone", zone));
      }
      block->property_add(Ioss::Property("zone", i + 1));
      block->property_add(Ioss::Property("id", i + 1));
      // Note that 'zone' is not consistent among processors
      block->property_add(Ioss::Property("guid", util().generate_guid(i + 1)));
      get_region()->add(block);
      m_zoneNameMap[zone_name] = i + 1;

      if (native) {
        // Handle zone-grid-connectivity...
        add_zgc_fpp(get_file_pointer(), block, m_zoneNameMap, myProcessor, isParallel);

        // Handle boundary conditions...
        Utils::add_structured_boundary_conditions(get_file_pointer(), block, false);
      }

      // Need to get a count of number of unique BC's.
      // Note that possible to assign multiple BC to a single face, so can't do this based on faces
      // Assume that if a BC is on multiple processors, then its name will be the same on all
      // processors.
      // * Gather all names to processor 0;
      // * Get unique ordered set
      // * Broadcast back to each processor
      int               in_bc  = 0;
      size_t            num_bc = block->m_boundaryConditions.size();
      std::vector<char> bc_names(num_bc * (CGNS_MAX_NAME_LENGTH + 1));
      for (size_t ibc = 0; ibc < num_bc; ibc++) {
        std::string name = block->m_boundaryConditions[ibc].m_famName + "/" +
                           block->m_boundaryConditions[ibc].m_bcName;
        strncpy(&bc_names[in_bc], name.c_str(), CGNS_MAX_NAME_LENGTH);
        in_bc += CGNS_MAX_NAME_LENGTH + 1;
      }
      std::vector<char> all_bc_names;
      int tot_names = util().gather(num_bc, CGNS_MAX_NAME_LENGTH + 1, bc_names, all_bc_names);

      if (myProcessor == 0) {
        int                      off_name = 0;
        std::vector<std::string> bc;
        for (int ibc = 0; ibc < tot_names; ibc++) {
          bc.emplace_back(&all_bc_names[off_name]);
          off_name += CGNS_MAX_NAME_LENGTH + 1;
        }
        Ioss::Utils::uniquify(bc);
        tot_names = (int)bc.size();
        all_bc_names.clear();
        all_bc_names.shrink_to_fit();
        bc_names.resize(tot_names * (CGNS_MAX_NAME_LENGTH + 1));
        in_bc = 0;
        for (const auto &name : bc) {
          strncpy(&bc_names[in_bc], name.c_str(), CGNS_MAX_NAME_LENGTH);
          in_bc += CGNS_MAX_NAME_LENGTH + 1;
        }
      }
      MPI_Bcast(&tot_names, 1, MPI_INT, 0, util().communicator());
      bc_names.resize(tot_names * (CGNS_MAX_NAME_LENGTH + 1));
      MPI_Bcast(bc_names.data(), tot_names * (CGNS_MAX_NAME_LENGTH + 1), MPI_CHAR, 0,
                util().communicator());

      std::vector<std::string> bc;
      int                      off_name = 0;
      for (int ibc = 0; ibc < tot_names; ibc++) {
        bc.emplace_back(&bc_names[off_name]);
        off_name += CGNS_MAX_NAME_LENGTH + 1;
      }
      bc_names.clear();
      bc_names.shrink_to_fit();

      // Each processor now has a unique set of BC names for this block.
      // Now create the missing (empty) BC on each processor.
      for (const auto &bc_name : bc) {
        auto split_name = Ioss::tokenize(bc_name, "/");
        assert(split_name.size() == 2);
        bool has_name = false;
        for (auto &sbc : block->m_boundaryConditions) {
          if (sbc.m_bcName == split_name[1]) {
            has_name = true;
            break;
          }
        }
        if (!has_name) {
          // Create an empty BC with that name...
          int            face = -1;
          Ioss::SideSet *sset = get_region()->get_sideset(split_name[0]);
          assert(sset != nullptr);
          add_empty_bc(sset, block, base, zone, face, split_name[0], split_name[1]);
        }
      }

      std::sort(block->m_boundaryConditions.begin(), block->m_boundaryConditions.end(),
                [](const Ioss::BoundaryCondition &b1, const Ioss::BoundaryCondition &b2) {
                  return (b1.m_bcName < b2.m_bcName);
                });
    }
#endif
  }

  void DatabaseIO::create_structured_block(int base, int zone, size_t &num_node)
  {
    assert(!isParallel);

    cgsize_t size[9];
    char     zone_name[CGNS_MAX_NAME_LENGTH + 1];
    CGCHECK(cg_zone_read(get_file_pointer(), base, zone, zone_name, size));

    auto        name_proc = decompose_name(zone_name, isParallel);
    std::string zname     = name_proc.first;
    int         proc      = name_proc.second;
    if (proc != myProcessor) {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: Zone " << zone
             << " has a name that specifies it should be on processor " << proc
             << ", but it is actually on processor " << myProcessor;
      IOSS_ERROR(errmsg);
    }

    m_zoneNameMap[zname] = zone;

    assert(size[0] - 1 == size[3]);
    assert(size[1] - 1 == size[4]);
    assert(size[2] - 1 == size[5]);

    assert(size[6] == 0);
    assert(size[7] == 0);
    assert(size[8] == 0);

    int index_dim = 0;
    CGCHECK(cg_index_dim(get_file_pointer(), base, zone, &index_dim));
    // An Ioss::StructuredBlock corresponds to a CG_Structured zone...
    Ioss::StructuredBlock *block =
        new Ioss::StructuredBlock(this, zname, index_dim, size[3], size[4], size[5]);

    block->property_add(Ioss::Property("base", base));
    block->property_add(Ioss::Property("db_zone", zone));
    block->property_add(Ioss::Property("zone", zone));
    block->property_add(Ioss::Property("id", zone));
    block->property_add(Ioss::Property("guid", zone));
    get_region()->add(block);

    num_node += block->get_property("node_count").get_int();

    // Handle zone-grid-connectivity...
    int nconn = 0;
    CGCHECK(cg_n1to1(get_file_pointer(), base, zone, &nconn));
    for (int i = 0; i < nconn; i++) {
      char                    connectname[CGNS_MAX_NAME_LENGTH + 1];
      char                    donorname[CGNS_MAX_NAME_LENGTH + 1];
      std::array<cgsize_t, 6> range;
      std::array<cgsize_t, 6> donor_range;
      Ioss::IJK_t             transform;

      CGCHECK(cg_1to1_read(get_file_pointer(), base, zone, i + 1, connectname, donorname, range.data(),
                           donor_range.data(), transform.data()));

      auto        donorname_proc = decompose_name(donorname, isParallel);
      std::string donor_name     = donorname_proc.first;

      // Get number of nodes shared with other "previous" zones...
      // A "previous" zone will have a lower zone number this this zone...
      int  donor_zone = -1;
      auto donor_iter = m_zoneNameMap.find(donor_name);
      if (donor_iter != m_zoneNameMap.end()) {
        donor_zone = (*donor_iter).second;
      }
      Ioss::IJK_t range_beg{{(int)range[0], (int)range[1], (int)range[2]}};
      Ioss::IJK_t range_end{{(int)range[3], (int)range[4], (int)range[5]}};
      Ioss::IJK_t donor_beg{{(int)donor_range[0], (int)donor_range[1], (int)donor_range[2]}};
      Ioss::IJK_t donor_end{{(int)donor_range[3], (int)donor_range[4], (int)donor_range[5]}};

      block->m_zoneConnectivity.emplace_back(connectname, zone, donor_name, donor_zone, transform,
                                             range_beg, range_end, donor_beg, donor_end);

      block->m_zoneConnectivity.back().m_ownerProcessor = myProcessor;
      block->m_zoneConnectivity.back().m_donorProcessor = donorname_proc.second;
    }

    // Handle boundary conditions...
    Utils::add_structured_boundary_conditions(get_file_pointer(), block, false);
  }

  size_t DatabaseIO::finalize_structured_blocks()
  {
    const auto &blocks = get_region()->get_structured_blocks();

    int              proc_count = util().parallel_size();
    std::vector<int> my_offsets;
    std::vector<int> all_offsets;

    if (proc_count > 1) {
      my_offsets.reserve(blocks.size() * 3 * proc_count);
#ifndef NDEBUG
      int zone = 1;
#endif
      for (const auto &sb : blocks) {
        assert(sb->get_property("zone").get_int() == zone++);
        my_offsets.push_back(sb->get_property("offset_i").get_int());
        my_offsets.push_back(sb->get_property("offset_j").get_int());
        my_offsets.push_back(sb->get_property("offset_k").get_int());
      }
      util().all_gather(my_offsets, all_offsets);
    }

    // If there are any Structured blocks, need to iterate them and their 1-to-1 connections
    // and update the donor_zone id for zones that had not yet been processed at the time of
    // definition...

    // If parallel, then all need to update the donor offset field since that was not known
    // at time of definition...
    for (auto &block : blocks) {
      for (auto &conn : block->m_zoneConnectivity) {
        if (conn.m_donorZone < 0) {
          auto donor_iter = m_zoneNameMap.find(conn.m_donorName);
          assert(donor_iter != m_zoneNameMap.end());
          conn.m_donorZone = (*donor_iter).second;
        }
        if (proc_count > 1) {
          int         offset = (conn.m_donorProcessor * blocks.size() + (conn.m_donorZone - 1)) * 3;
          Ioss::IJK_t donor_offset{
              {all_offsets[offset + 0], all_offsets[offset + 1], all_offsets[offset + 2]}};

          conn.m_donorOffset = donor_offset;
          conn.m_donorRangeBeg[0] += donor_offset[0];
          conn.m_donorRangeBeg[1] += donor_offset[1];
          conn.m_donorRangeBeg[2] += donor_offset[2];
          conn.m_donorRangeEnd[0] += donor_offset[0];
          conn.m_donorRangeEnd[1] += donor_offset[1];
          conn.m_donorRangeEnd[2] += donor_offset[2];
        }
        conn.m_donorGUID = util().generate_guid(conn.m_donorZone, conn.m_donorProcessor);
        conn.m_ownerGUID = util().generate_guid(conn.m_ownerZone, conn.m_ownerProcessor);
      }
    }

    size_t num_nodes = Utils::resolve_nodes(*get_region(), myProcessor, isParallel);
    return num_nodes;
  }

  void DatabaseIO::create_unstructured_block(int base, int zone, size_t &num_node)
  {
    cgsize_t size[9];
    char     zone_name[CGNS_MAX_NAME_LENGTH + 1];
    CGCHECK(cg_zone_read(get_file_pointer(), base, zone, zone_name, size));
    m_zoneNameMap[zone_name] = zone;

    size_t total_block_nodes = size[0];
    m_blockLocalNodeMap[zone].resize(total_block_nodes, -1);

    // Determine number of "shared" nodes (shared with other zones)
    if (zone > 1) { // Donor zone is always lower numbered, so zone 1 has no donor zone.
      int nconn = 0;
      CGCHECK(cg_nconns(get_file_pointer(), base, zone, &nconn));
      cgsize_t num_shared = 0;
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

        CGCHECK(cg_conn_info(get_file_pointer(), base, zone, i + 1, connectname, &location, &connect_type,
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
        auto donor_iter = m_zoneNameMap.find(donorname);
        if (donor_iter != m_zoneNameMap.end() && (*donor_iter).second < zone) {
          num_shared += npnts;
#if IOSS_DEBUG_OUTPUT
          std::cout << "Zone " << zone << " shares " << npnts << " nodes with " << donorname
                    << "\n";
#endif
          std::vector<cgsize_t> points(npnts);
          std::vector<cgsize_t> donors(npnts);

          CGCHECK(cg_conn_read(get_file_pointer(), base, zone, i + 1, TOPTR(points), donor_datatype,
                               TOPTR(donors)));

          // Fill in entries in m_blockLocalNodeMap for the shared nodes...
          auto &donor_map = m_blockLocalNodeMap[(*donor_iter).second];
          auto &block_map = m_blockLocalNodeMap[zone];
          for (int j = 0; j < npnts; j++) {
            cgsize_t point       = points[j];
            cgsize_t donor       = donors[j];
            block_map[point - 1] = donor_map[donor - 1];
          }
        }
      }
    }

    auto & block_map = m_blockLocalNodeMap[zone];
    size_t offset    = num_node;
    for (size_t i = 0; i < total_block_nodes; i++) {
      if (block_map[i] == -1) {
        block_map[i] = offset++;
      }
    }
    num_node = offset;

    size_t num_elem    = size[1];
    m_zoneOffset[zone] = m_zoneOffset[zone - 1] + num_elem;

    // NOTE: A Zone will have a single set of nodes, but can have
    //       multiple sections each with their own element type...
    //       Keep treating sections as element blocks until we
    //       have handled 'size[1]' number of elements; the remaining
    //       sections are then the boundary faces (?)
    int num_sections = 0;
    CGCHECK(cg_nsections(get_file_pointer(), base, zone, &num_sections));

    // ========================================================================
    // Read the sections and create an element block for the ones that
    // define elements.  Some define boundary conditions...
    Ioss::ElementBlock *eblock = nullptr;

    for (int is = 1; is <= num_sections; is++) {
      char             section_name[CGNS_MAX_NAME_LENGTH + 1];
      CG_ElementType_t e_type;
      cgsize_t         el_start    = 0;
      cgsize_t         el_end      = 0;
      int              num_bndry   = 0;
      int              parent_flag = 0;

      // Get the type of elements in this section...
      CGCHECK(cg_section_read(get_file_pointer(), base, zone, is, section_name, &e_type, &el_start,
                              &el_end, &num_bndry, &parent_flag));

      cgsize_t num_entity = el_end - el_start + 1;

      if (parent_flag == 0 && num_elem > 0) {
        num_elem -= num_entity;
        std::string element_topo = Utils::map_cgns_to_topology_type(e_type);
#if IOSS_DEBUG_OUTPUT
        std::cout << "Added block " << zone_name << ": CGNS topology = '"
                  << cg_ElementTypeName(e_type) << "', IOSS topology = '" << element_topo
                  << "' with " << num_entity << " elements\n";
#endif
        eblock = new Ioss::ElementBlock(this, zone_name, element_topo, num_entity);
        eblock->property_add(Ioss::Property("base", base));
        eblock->property_add(Ioss::Property("zone", zone));
        eblock->property_add(Ioss::Property("db_zone", zone));
        eblock->property_add(Ioss::Property("id", zone));
        eblock->property_add(Ioss::Property("guid", zone));
        eblock->property_add(Ioss::Property("section", is));
        eblock->property_add(Ioss::Property("node_count", (int64_t)total_block_nodes));
        eblock->property_add(Ioss::Property("original_block_order", zone));

        assert(is == 1); // For now, assume each zone has only a single element block.
        bool added = get_region()->add(eblock);
        if (!added) {
          delete eblock;
          eblock = nullptr;
        }
      }
      else {
        // This is a boundary-condition -- sideset (?)
        // See if there is an existing sideset with this name...
        Ioss::SideSet *sset = get_region()->get_sideset(section_name);

        if (sset != nullptr) {
          std::string block_name(zone_name);
          block_name += "/";
          block_name += section_name;
          std::string face_topo = Utils::map_cgns_to_topology_type(e_type);
#if IOSS_DEBUG_OUTPUT
          std::cout << "Added sideset " << block_name << " of topo " << face_topo << " with "
                    << num_entity << " faces\n";
#endif
          std::string parent_topo = eblock == nullptr ? "unknown" : eblock->topology()->name();
          auto sblk = new Ioss::SideBlock(this, block_name, face_topo, parent_topo, num_entity);
          // IF name is of form "surface_" + "#", then extract # and use as id...
          int id = Ioss::Utils::extract_id(block_name);
          if (id != 0) {
            sblk->property_add(Ioss::Property("id", id));
            sblk->property_add(Ioss::Property("guid", id));
          }
          else {
            sblk->property_add(Ioss::Property("id", zone));
            sblk->property_add(Ioss::Property("guid", zone));
          }
          sblk->property_add(Ioss::Property("base", base));
          sblk->property_add(Ioss::Property("zone", zone));
          sblk->property_add(Ioss::Property("section", is));
          if (eblock != nullptr) {
            sblk->set_parent_element_block(eblock);
          }
          sset->add(sblk);
        }
      }
    }
  }

  void DatabaseIO::read_meta_data__()
  {
    // Determine the number of bases in the grid.
    // Currently only handle 1.
    int n_bases = 0;
    CGCHECK(cg_nbases(get_file_pointer(), &n_bases));
    if (n_bases != 1) {
      std::ostringstream errmsg;
      errmsg << "CGNS: Too many bases; only support files with a single bases at this time";
      IOSS_ERROR(errmsg);
    }

    get_step_times__();

    // ========================================================================
    // Get the number of families in the mesh...
    // Will treat these as sidesets if they are of the type "FamilyBC_t"
    Utils::add_sidesets(get_file_pointer(), this);

    // ========================================================================
    // Get the number of zones (element blocks) in the mesh...
    int num_zones = 0;
    int base      = 1;
    CGCHECK(cg_nzones(get_file_pointer(), base, &num_zones));
    m_blockLocalNodeMap.resize(num_zones + 1); // Let's use 1-based zones...
    m_zoneOffset.resize(num_zones + 1);        // Let's use 1-based zones...

    // ========================================================================
    size_t        num_node         = 0;
    CG_ZoneType_t common_zone_type = Utils::check_zone_type(get_file_pointer());

    if (isParallel && common_zone_type == CG_Structured) {
      // Handle the file-per-processor parallel case separately for
      // now. Hopefully can consolidate at some later time.
      create_structured_block_fpp(base, num_zones, num_node);
    }
    else {
      for (int zone = 1; zone <= num_zones; zone++) {
        if (common_zone_type == CG_Structured) {
          create_structured_block(base, zone, num_node);
        }
        else if (common_zone_type == CG_Unstructured) {
          create_unstructured_block(base, zone, num_node);
        }
        else {
          // This should be handled already in check_zone_type...
          std::ostringstream errmsg;
          errmsg << "ERROR: CGNS: Zone " << zone
                 << " is not of type Unstructured or Structured "
                    "which are the only types currently supported";
          IOSS_ERROR(errmsg);
        }
      }
    }

    if (common_zone_type == CG_Structured) {
      num_node = finalize_structured_blocks();
    }

    char basename[CGNS_MAX_NAME_LENGTH + 1];
    int  cell_dimension = 0;
    int  phys_dimension = 0;
    CGCHECK(cg_base_read(get_file_pointer(), base, basename, &cell_dimension, &phys_dimension));
    if (phys_dimension != 3) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The model is " << phys_dimension << "D.  Only 3D models are supported.";
      IOSS_ERROR(errmsg);
    }

    Ioss::NodeBlock *nblock = new Ioss::NodeBlock(this, "nodeblock_1", num_node, phys_dimension);
    nblock->property_add(Ioss::Property("base", base));
    get_region()->add(nblock);
    nodeCount = num_node;

    Utils::add_transient_variables(get_file_pointer(), m_timesteps, get_region(), get_field_recognition(),
                                   get_field_separator(), myProcessor, false);
  }

  void DatabaseIO::write_meta_data()
  {
    int num_zones = get_region()->get_property("element_block_count").get_int() +
                    get_region()->get_property("structured_block_count").get_int();
    m_bcOffset.resize(num_zones + 1);   // use 1-based zones...
    m_zoneOffset.resize(num_zones + 1); // use 1-based zones...

    elementCount = Utils::common_write_meta_data(get_file_pointer(), *get_region(), m_zoneOffset, false);
  }

  void DatabaseIO::get_step_times__()
  {
    Utils::get_step_times(get_file_pointer(), m_timesteps, get_region(), timeScaleFactor, myProcessor);
  }

  void DatabaseIO::write_adjacency_data()
  {
    // Determine adjacency information between unstructured blocks.
    // Could save this information from the input mesh, but then
    // could not read an exodus mesh and write a cgns mesh.
    // However, in long run may still want to read/save input adjacency
    // data if doing cgns -> cgns...  For now, try generating information.

    // If block I is adjacent to block J, then they will share at
    // least 1 "side" (face 3D or edge 2D).
    // Currently, assuming they are adjacent if they share at least one node...

    size_t node_count = get_region()->get_property("node_count").get_int();

    const auto &blocks = get_region()->get_element_blocks();
    for (auto I = blocks.cbegin(); I != blocks.cend(); I++) {
      int base = (*I)->get_property("base").get_int();
      int zone = Iocgns::Utils::get_db_zone(*I);

      const auto &I_map = m_globalToBlockLocalNodeMap[zone];

      // Flag all nodes used by this block...
      std::vector<size_t> I_nodes(node_count);
      for (size_t i = 0; i < I_map->size(); i++) {
        auto global     = I_map->map()[i + 1] - 1;
        I_nodes[global] = i + 1;
      }
      for (auto J = I + 1; J != blocks.end(); J++) {
        int                   dzone = (*J)->get_property("zone").get_int();
        const auto &          J_map = m_globalToBlockLocalNodeMap[dzone];
        std::vector<cgsize_t> point_list;
        std::vector<cgsize_t> point_list_donor;
        for (size_t i = 0; i < J_map->size(); i++) {
          auto global = J_map->map()[i + 1] - 1;
          if (I_nodes[global] > 0) {
            // Have a match between nodes used by two different blocks,
            // They are adjacent...
            point_list.push_back(I_nodes[global]);
            point_list_donor.push_back(i + 1);
          }
        }

        // If point_list non_empty, then output this adjacency node...
        if (!point_list.empty()) {
          int         gc_idx = 0;
          std::string name   = (*I)->name();
          name += "_to_";
          name += (*J)->name();
          const auto &d1_name = (*J)->name();
          CGCHECK(cg_conn_write(get_file_pointer(), base, zone, name.c_str(), CG_Vertex, CG_Abutting1to1,
                                CG_PointList, point_list.size(), TOPTR(point_list), d1_name.c_str(),
                                CG_Unstructured, CG_PointListDonor, CG_DataTypeNull,
                                point_list_donor.size(), TOPTR(point_list_donor), &gc_idx));

          name = (*J)->name();
          name += "_to_";
          name += (*I)->name();
          const auto &d2_name = (*I)->name();

          CGCHECK(cg_conn_write(get_file_pointer(), base, dzone, name.c_str(), CG_Vertex, CG_Abutting1to1,
                                CG_PointList, point_list_donor.size(), TOPTR(point_list_donor),
                                d2_name.c_str(), CG_Unstructured, CG_PointListDonor,
                                CG_DataTypeNull, point_list.size(), TOPTR(point_list), &gc_idx));
        }
      }
    }
  }

  bool DatabaseIO::begin__(Ioss::State state)
  {
    dbState = state;
    return true;
  }

  bool DatabaseIO::end__(Ioss::State state)
  {
    // Transitioning out of state 'state'
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_meta_data();
      }
      break;
    case Ioss::STATE_MODEL:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_adjacency_data();
      }
      break;
    case Ioss::STATE_DEFINE_TRANSIENT:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_results_meta_data();
      }
      break;
    default: // ignore everything else...
      break;
    }

    dbState = Ioss::STATE_UNKNOWN;
    return true;
  }

  bool DatabaseIO::begin_state__(Ioss::Region *region, int state, double time)
  {
    if (is_input()) {
      return true;
    }
    Utils::write_flow_solution_metadata(get_file_pointer(), get_region(), state,
                                        &m_currentVertexSolutionIndex,
                                        &m_currentCellCenterSolutionIndex, false);

    return true;
  }

  bool DatabaseIO::end_state__(Ioss::Region * /* region */, int state, double time)
  {
    if (!is_input()) {
      m_timesteps.push_back(time);
      assert(m_timesteps.size() == (size_t)state);
    }
    return true;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(reg, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role       = field.get_role();
    int                   base       = nb->get_property("base").get_int();
    size_t                num_to_get = field.verify(data_size);
    cgsize_t              first      = 1;

    char basename[CGNS_MAX_NAME_LENGTH + 1];

    // Create a lambda to eliminate lots of duplicate code in coordinate outputs...
    auto coord_lambda = [this, &data, &first, base](const char *ordinate) {
      double *rdata = static_cast<double *>(data);

      for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
        auto &              block_map = m_blockLocalNodeMap[zone];
        cgsize_t            num_coord = block_map.size();
        std::vector<double> coord(num_coord);
        CGCHECK(cg_coord_read(get_file_pointer(), base, zone, ordinate, CG_RealDouble, &first, &num_coord,
                              TOPTR(coord)));

        // Map to global coordinate position...
        for (cgsize_t i = 0; i < num_coord; i++) {
          rdata[block_map[i]] = coord[i];
        }
      }
    };
    // End of lambda...

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
        // Use the lambda...
        coord_lambda("CoordinateX");
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        coord_lambda("CoordinateY");
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        coord_lambda("CoordinateZ");
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        int cell_dimension = 0;
        int phys_dimension = 0;
        CGCHECK(cg_base_read(get_file_pointer(), base, basename, &cell_dimension, &phys_dimension));

        double *rdata = static_cast<double *>(data);

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'
        for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
          auto &              block_map = m_blockLocalNodeMap[zone];
          cgsize_t            num_coord = block_map.size();
          std::vector<double> coord(num_coord);

          // ========================================================================
          // Repetitive code for each coordinate direction; use a lambda to consolidate...
          auto blk_coord_lambda = [this, block_map, base, zone, &coord, first, num_coord,
                                   phys_dimension, &rdata](const char *ord_name, int ordinate) {
            CGCHECK(cg_coord_read(get_file_pointer(), base, zone, ord_name, CG_RealDouble, &first,
                                  &num_coord, coord.data()));

            // Map to global coordinate position...
            for (cgsize_t i = 0; i < num_coord; i++) {
              rdata[phys_dimension * block_map[i] + ordinate] = coord[i];
            }
          };
          // End of lambda...
          // ========================================================================

          blk_coord_lambda("CoordinateX", 0);

          if (phys_dimension >= 2) {
            blk_coord_lambda("CoordinateY", 1);
          }

          if (phys_dimension >= 3) {
            blk_coord_lambda("CoordinateZ", 2);
          }
        }
      }
      else if (field.get_name() == "ids") {
        // Map the local ids in this node block
        // (1...node_count) to global node ids.
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          std::iota(idata, idata + num_to_get, 1);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          std::iota(idata, idata + num_to_get, 1);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      // Locate the FlowSolution node corresponding to the correct state/step/time
      // TODO: do this at read_meta_data() and store...
      int step = get_region()->get_current_state();

      for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
        int   solution_index = Utils::find_solution_index(get_file_pointer(), base, zone, step, CG_Vertex);
        auto &block_map      = m_blockLocalNodeMap[zone];
        cgsize_t num_block_node = block_map.size();

        double *            rdata        = static_cast<double *>(data);
        cgsize_t            range_min[1] = {1};
        cgsize_t            range_max[1] = {num_block_node};
        auto                var_type     = field.transformed_storage();
        int                 comp_count   = var_type->component_count();
        std::vector<double> cgns_data(num_block_node);
        if (comp_count == 1) {
          CGCHECK(cg_field_read(get_file_pointer(), base, zone, solution_index, field.get_name().c_str(),
                                CG_RealDouble, range_min, range_max, cgns_data.data()));

          // Map to global nodal field position...
          for (cgsize_t i = 0; i < num_block_node; i++) {
            rdata[block_map[i]] = cgns_data[i];
          }
        }
        else {
          char field_suffix_separator = get_field_separator();
          for (int i = 0; i < comp_count; i++) {
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_read(get_file_pointer(), base, zone, solution_index, var_name.c_str(),
                                  CG_RealDouble, range_min, range_max, cgns_data.data()));
            for (cgsize_t j = 0; j < num_block_node; j++) {
              auto global                    = block_map[j];
              rdata[comp_count * global + i] = cgns_data[j];
            }
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(nb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(eb, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fb, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      int                   base             = eb->get_property("base").get_int();
      int                   zone             = Iocgns::Utils::get_db_zone(eb);
      int                   sect             = eb->get_property("section").get_int();
      cgsize_t              my_element_count = eb->entity_count();
      Ioss::Field::RoleType role             = field.get_role();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for a CGNS file model.
        // (The 'genesis' portion)

        if (field.get_name() == "connectivity" || field.get_name() == "connectivity_raw") {
          // TODO(gdsjaar): Need to map local to global...
          int element_nodes = eb->topology()->number_nodes();
          assert(field.raw_storage()->component_count() == element_nodes);

          if (my_element_count > 0) {
            int field_byte_size = (field.get_type() == Ioss::Field::INT32) ? 32 : 64;
            if (field_byte_size == CG_SIZEOF_SIZE) {
              cgsize_t *idata = reinterpret_cast<cgsize_t *>(data);
              CGCHECK(cg_elements_read(get_file_pointer(), base, zone, sect, idata, nullptr));
            }
            else {
              std::vector<cgsize_t> connect(element_nodes * num_to_get);
              CGCHECK(cg_elements_read(get_file_pointer(), base, zone, sect, connect.data(), nullptr));
              if (field.get_type() == Ioss::Field::INT32) {
                auto * idata = reinterpret_cast<int *>(data);
                size_t i     = 0;
                for (auto node : connect) {
                  idata[i++] = node;
                }
              }
              else {
                auto * idata = reinterpret_cast<int64_t *>(data);
                size_t i     = 0;
                for (auto node : connect) {
                  idata[i++] = node;
                }
              }
            }
          }

          // Now need to map block-local node connectivity to global nodes...
          // This is done for both connectivity and connectivity_raw
          // since the "global id" is the same as the "local id"
          // The connectivities we currently have are "block local"
          const auto &block_map = m_blockLocalNodeMap[zone];
          if (field.get_type() == Ioss::Field::INT32) {
            int *idata = static_cast<int *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              idata[i] = block_map[idata[i] - 1] + 1;
            }
          }
          else {
            int64_t *idata = static_cast<int64_t *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              idata[i] = block_map[idata[i] - 1] + 1;
            }
          }
        }
        else if (field.get_name() == "ids") {
          // Map the local ids in this element block
          // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
          size_t eb_offset_plus_one = eb->get_offset() + 1;
          if (field.get_type() == Ioss::Field::INT64) {
            int64_t *idata = static_cast<int64_t *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
          else {
            assert(field.get_type() == Ioss::Field::INT32);
            int *idata = static_cast<int *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
        }
        else if (field.get_name() == "implicit_ids") {
          size_t eb_offset_plus_one = eb->get_offset() + 1;
          if (field.get_type() == Ioss::Field::INT64) {
            int64_t *idata = static_cast<int64_t *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
          else {
            assert(field.get_type() == Ioss::Field::INT32);
            int *idata = static_cast<int *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
        }
        else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "input");
        }
      }
      else if (role == Ioss::Field::TRANSIENT) {
        // Locate the FlowSolution node corresponding to the correct state/step/time
        // TODO: do this at read_meta_data() and store...
        int step = get_region()->get_current_state();
        int solution_index =
            Utils::find_solution_index(get_file_pointer(), base, zone, step, CG_CellCenter);

        double * rdata        = static_cast<double *>(data);
        cgsize_t range_min[1] = {1};
        cgsize_t range_max[1] = {my_element_count};

        auto var_type   = field.transformed_storage();
        int  comp_count = var_type->component_count();
        if (comp_count == 1) {
          CGCHECK(cg_field_read(get_file_pointer(), base, zone, solution_index, field.get_name().c_str(),
                                CG_RealDouble, range_min, range_max, rdata));
        }
        else {
          std::vector<double> cgns_data(my_element_count);
          for (int i = 0; i < comp_count; i++) {
            char        field_suffix_separator = get_field_separator();
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_read(get_file_pointer(), base, zone, solution_index, var_name.c_str(),
                                  CG_RealDouble, range_min, range_max, cgns_data.data()));
            for (cgsize_t j = 0; j < my_element_count; j++) {
              rdata[comp_count * j + i] = cgns_data[j];
            }
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "output");
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role = field.get_role();
    int                   base = sb->get_property("base").get_int();
    int                   zone = Iocgns::Utils::get_db_zone(sb);

    cgsize_t num_to_get = field.verify(data_size);

    cgsize_t rmin[3] = {0, 0, 0};
    cgsize_t rmax[3] = {0, 0, 0};

    bool cell_field = Utils::is_cell_field(field);
    if ((cell_field && sb->get_property("cell_count").get_int() == 0) ||
        (!cell_field && sb->get_property("node_count").get_int() == 0)) {
      return 0;
    }

    if (cell_field) {
      assert(num_to_get == sb->get_property("cell_count").get_int());
      if (num_to_get > 0) {
        rmin[0] = 1;
        rmin[1] = 1;
        rmin[2] = 1;

        rmax[0] = rmin[0] + sb->get_property("ni").get_int() - 1;
        rmax[1] = rmin[1] + sb->get_property("nj").get_int() - 1;
        rmax[2] = rmin[2] + sb->get_property("nk").get_int() - 1;
      }
    }
    else {
      // cell nodal field.
      assert(num_to_get == sb->get_property("node_count").get_int());
      if (num_to_get > 0) {
        rmin[0] = 1;
        rmin[1] = 1;
        rmin[2] = 1;

        rmax[0] = rmin[0] + sb->get_property("ni").get_int();
        rmax[1] = rmin[1] + sb->get_property("nj").get_int();
        rmax[2] = rmin[2] + sb->get_property("nk").get_int();
      }
    }

    assert(num_to_get ==
           (rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1));
    double *rdata = static_cast<double *>(data);

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
        CGCHECK(cg_coord_read(get_file_pointer(), base, zone, "CoordinateX", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        CGCHECK(cg_coord_read(get_file_pointer(), base, zone, "CoordinateY", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        CGCHECK(cg_coord_read(get_file_pointer(), base, zone, "CoordinateZ", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        char basename[CGNS_MAX_NAME_LENGTH + 1];
        int  cell_dimension = 0;
        int  phys_dimension = 0;
        CGCHECK(cg_base_read(get_file_pointer(), base, basename, &cell_dimension, &phys_dimension));

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'

        std::vector<double> coord(num_to_get);

        // ========================================================================
        // Repetitive code for each coordinate direction; use a lambda to consolidate...
        auto coord_lambda = [this, base, zone, &coord, rmin, rmax, phys_dimension, num_to_get,
                             &rdata](const char *ord_name, int ordinate) {
          CGCHECK(cg_coord_read(get_file_pointer(), base, zone, ord_name, CG_RealDouble, rmin, rmax,
                                TOPTR(coord)));

          // Map to global coordinate position...
          for (cgsize_t i = 0; i < num_to_get; i++) {
            rdata[phys_dimension * i + ordinate] = coord[i];
          }
        };
        // End of lambda...
        // ========================================================================

        coord_lambda("CoordinateX", 0);

        if (phys_dimension >= 2) {
          coord_lambda("CoordinateY", 1);
        }

        if (phys_dimension == 3) {
          coord_lambda("CoordinateZ", 2);
        }
      }
      else if (field.get_name() == "cell_node_ids") {
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          sb->get_cell_node_ids(idata, true);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          sb->get_cell_node_ids(idata, true);
        }
      }
      else if (field.get_name() == "cell_ids") {
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          sb->get_cell_ids(idata, true);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          sb->get_cell_ids(idata, true);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      auto var_type               = field.transformed_storage();
      int  comp_count             = var_type->component_count();
      char field_suffix_separator = get_field_separator();

      int sol_index = 0;
      int step      = get_region()->get_current_state();
      if (cell_field) {
        sol_index = Utils::find_solution_index(get_file_pointer(), base, zone, step, CG_CellCenter);
      }
      else {
        sol_index = Utils::find_solution_index(get_file_pointer(), base, zone, step, CG_Vertex);
      }

      if (comp_count == 1) {
        CGCHECK(cg_field_read(get_file_pointer(), base, zone, sol_index, field.get_name().c_str(),
                              CG_RealDouble, rmin, rmax, rdata));
      }
      else {
        std::vector<double> cgns_data(num_to_get);
        for (int i = 0; i < comp_count; i++) {
          std::string var_name =
              var_type->label_name(field.get_name(), i + 1, field_suffix_separator);
          CGCHECK(cg_field_read(get_file_pointer(), base, zone, sol_index, var_name.c_str(), CG_RealDouble,
                                rmin, rmax, cgns_data.data()));
          for (cgsize_t j = 0; j < num_to_get; j++) {
            rdata[comp_count * j + i] = cgns_data[j];
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ns, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    int base = sb->get_property("base").get_int();
    int zone = Iocgns::Utils::get_db_zone(sb);
    int sect = sb->get_property("section").get_int();

    ssize_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      int64_t entity_count = sb->entity_count();
      if (num_to_get != entity_count) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Partial field input not yet implemented for side blocks";
        IOSS_ERROR(errmsg);
      }
    }

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "element_side_raw" || field.get_name() == "element_side") {

        // TODO(gdsjaar): ? Possibly rewrite using cgi_read_int_data so can skip reading element
        // connectivity
        int                   nodes_per_face = sb->topology()->number_nodes();
        std::vector<cgsize_t> elements(nodes_per_face * num_to_get); // Not needed, but can't skip

        // We get:
        // *  num_to_get parent elements,
        // *  num_to_get zeros (other parent element for face, but on boundary so 0)
        // *  num_to_get face_on_element
        // *  num_to_get zeros (face on other parent element)
        std::vector<cgsize_t> parent(4 * num_to_get);

        CGCHECK(cg_elements_read(get_file_pointer(), base, zone, sect, TOPTR(elements), TOPTR(parent)));

        size_t offset = m_zoneOffset[zone - 1];
        if (field.get_type() == Ioss::Field::INT32) {
          int *  idata = reinterpret_cast<int *>(data);
          size_t j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            idata[j++] = parent[num_to_get * 0 + i] + offset; // Element
            idata[j++] = parent[num_to_get * 2 + i];
            assert(parent[num_to_get * 1 + i] == 0);
            assert(parent[num_to_get * 3 + i] == 0);
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
        }
        else {
          int64_t *idata = reinterpret_cast<int64_t *>(data);
          size_t   j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            idata[j++] = parent[num_to_get * 0 + i] + offset; // Element
            idata[j++] = parent[num_to_get * 2 + i];
            assert(parent[num_to_get * 1 + i] == 0);
            assert(parent[num_to_get * 3 + i] == 0);
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "input");
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(cs, field, "input");
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region *region, const Ioss::Field &field,
                                         void * /*data*/, size_t /*data_size*/) const
  {
    return Ioss::Utils::field_warning(region, field, "output");
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role = field.get_role();
    int                   base = sb->get_property("base").get_int();
    int                   zone = Iocgns::Utils::get_db_zone(sb);

    cgsize_t num_to_get = field.verify(data_size);

    // In this routine, if isParallel, then writing file-per-processor; not parallel io to single
    // file.
    if (isParallel && num_to_get == 0) {
      return 0;
    }

    if (role == Ioss::Field::MESH) {
      bool cell_field = Utils::is_cell_field(field);

      if (cell_field) {
        assert(num_to_get == sb->get_property("cell_count").get_int());
      }

      double *rdata = static_cast<double *>(data);

      int crd_idx = 0;
      if (field.get_name() == "mesh_model_coordinates_x") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateX", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateY", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateZ", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        assert(!cell_field);
        int phys_dimension = get_region()->get_property("spatial_dimension").get_int();

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'
        std::vector<double> coord(num_to_get);

        // ========================================================================
        // Repetitive code for each coordinate direction; use a lambda to consolidate...
        auto coord_lambda = [this, &coord, num_to_get, phys_dimension, &rdata, base,
                             zone](const char *ord_name, int ordinate) {
          int crd_index = 0;

          // Map to global coordinate position...
          for (cgsize_t i = 0; i < num_to_get; i++) {
            coord[i] = rdata[phys_dimension * i + ordinate];
          }

          CGCHECK(cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, ord_name, TOPTR(coord),
                                 &crd_index));
        };
        // End of lambda...
        // ========================================================================

        coord_lambda("CoordinateX", 0);

        if (phys_dimension >= 2) {
          coord_lambda("CoordinateY", 1);
        }
        if (phys_dimension == 3) {
          coord_lambda("CoordinateZ", 2);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "output");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      double *rdata                  = static_cast<double *>(data);
      int     cgns_field             = 0;
      auto    var_type               = field.transformed_storage();
      int     comp_count             = var_type->component_count();
      char    field_suffix_separator = get_field_separator();
      if (comp_count == 1) {
        CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentCellCenterSolutionIndex,
                               CG_RealDouble, field.get_name().c_str(), rdata, &cgns_field));
        Utils::set_field_index(field, cgns_field, CG_CellCenter);
      }
      else {
        std::vector<double> cgns_data(num_to_get);
        for (int i = 0; i < comp_count; i++) {
          for (cgsize_t j = 0; j < num_to_get; j++) {
            cgns_data[j] = rdata[comp_count * j + i];
          }
          std::string var_name =
              var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

          CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentCellCenterSolutionIndex,
                                 CG_RealDouble, var_name.c_str(), cgns_data.data(), &cgns_field));
          if (i == 0) {
            Utils::set_field_index(field, cgns_field, CG_CellCenter);
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }

    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for a CGNS file model.
        // (The 'genesis' portion)
        if (field.get_name() == "connectivity") {
          // This blocks zone has not been defined.
          // Get the "node block" for this element block...
          int element_nodes = eb->topology()->number_nodes();
          assert(field.raw_storage()->component_count() == element_nodes);

          Ioss::MapContainer nodes;
          nodes.reserve(element_nodes * num_to_get + 1);
          nodes.push_back(0); // Unknown whether one-to-one map.

          if (field.get_type() == Ioss::Field::INT32) {
            auto *idata = reinterpret_cast<int *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              nodes.push_back(idata[i]);
            }
          }
          else {
            auto *idata = reinterpret_cast<int64_t *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              nodes.push_back(idata[i]);
            }
          }
          Ioss::Utils::uniquify(nodes, true);
          assert(nodes[0] == 0);

          // Now, we have the node count and cell count so we can create a zone...
          int      base    = 1;
          int      zone    = 0;
          cgsize_t size[3] = {0, 0, 0};
          size[1]          = eb->entity_count();
          size[0]          = nodes.size() - 1;

          CGCHECK(
              cg_zone_write(get_file_pointer(), base, eb->name().c_str(), size, CG_Unstructured, &zone));
          eb->property_update("db_zone", zone);
          eb->property_update("zone", zone);
          eb->property_update("id", zone);
          eb->property_update("guid", zone);
          eb->property_update("section", 1);
          eb->property_update("base", base);

          // Now we have a valid zone so can update some data structures...
          m_zoneOffset[zone]                = m_zoneOffset[zone - 1] + size[1];
          m_globalToBlockLocalNodeMap[zone] = new Ioss::Map("element", "unknown", myProcessor);
          m_globalToBlockLocalNodeMap[zone]->map().swap(nodes);
          m_globalToBlockLocalNodeMap[zone]->build_reverse_map_no_lock();

          // Need to map global nodes to block-local node connectivity
          const auto &block_map = m_globalToBlockLocalNodeMap[zone];
          block_map->reverse_map_data(data, field, num_to_get * element_nodes);

          if (eb->entity_count() > 0) {
            CG_ElementType_t type            = Utils::map_topology_to_cgns(eb->topology()->name());
            int              sect            = 0;
            int              field_byte_size = (field.get_type() == Ioss::Field::INT32) ? 32 : 64;
            if (field_byte_size == CG_SIZEOF_SIZE) {
              CGCHECK(cg_section_write(get_file_pointer(), base, zone, "HexElements", type, 1, num_to_get,
                                       0, (cgsize_t *)data, &sect));
            }
            else {
              std::vector<cgsize_t> connect;
              connect.reserve(element_nodes * num_to_get);
              if (field.get_type() == Ioss::Field::INT32) {
                auto *idata = reinterpret_cast<int *>(data);
                for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                  connect.push_back(idata[i]);
                }
              }
              else {
                auto *idata = reinterpret_cast<int64_t *>(data);
                for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                  connect.push_back(idata[i]);
                }
              }
              CGCHECK(cg_section_write(get_file_pointer(), base, zone, "HexElements", type, 1, num_to_get,
                                       0, connect.data(), &sect));
            }
            m_bcOffset[zone] += num_to_get;
            eb->property_update("section", sect);
          }
        }
        else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "output");
        }
      }
      else if (role == Ioss::Field::TRANSIENT) {
        int     base                   = eb->get_property("base").get_int();
        int     zone                   = Iocgns::Utils::get_db_zone(eb);
        double *rdata                  = static_cast<double *>(data);
        int     cgns_field             = 0;
        auto    var_type               = field.transformed_storage();
        int     comp_count             = var_type->component_count();
        char    field_suffix_separator = get_field_separator();
        if (comp_count == 1) {
          CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentCellCenterSolutionIndex,
                                 CG_RealDouble, field.get_name().c_str(), rdata, &cgns_field));
          Utils::set_field_index(field, cgns_field, CG_CellCenter);
        }
        else {
          std::vector<double> cgns_data(num_to_get);
          for (int i = 0; i < comp_count; i++) {
            for (size_t j = 0; j < num_to_get; j++) {
              cgns_data[j] = rdata[comp_count * j + i];
            }
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentCellCenterSolutionIndex,
                                   CG_RealDouble, var_name.c_str(), cgns_data.data(), &cgns_field));
            if (i == 0) {
              Utils::set_field_index(field, cgns_field, CG_CellCenter);
            }
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "output");
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fb, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(eb, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    // Instead of outputting a global nodeblock's worth of data,
    // the data is output a "zone" at a time.
    // The m_globalToBlockLocalNodeMap[zone] map is used (Ioss::Map pointer)
    // This map is built during the output of block connectivity,
    // so for cgns unstructured mesh, we need to output ElementBlock connectivity
    // prior to outputting nodal coordinates.
    for (const auto &z : m_globalToBlockLocalNodeMap) {
      if (z.second == nullptr) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: The globalToBlockLocalNodeMap is not defined, so nodal fields "
                  "cannot be output.";
        IOSS_ERROR(errmsg);
      }
    }

    Ioss::Field::RoleType role       = field.get_role();
    int                   base       = 1;
    cgsize_t              num_to_get = field.verify(data_size);

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates" ||
          field.get_name() == "mesh_model_coordinates_x" ||
          field.get_name() == "mesh_model_coordinates_y" ||
          field.get_name() == "mesh_model_coordinates_z") {
        double *rdata = static_cast<double *>(data);

        if (field.get_name() == "mesh_model_coordinates") {
          int spatial_dim = nb->get_property("component_degree").get_int();
          for (const auto &block : m_globalToBlockLocalNodeMap) {
            auto zone = block.first;
            // NOTE: 'block_map' has one more entry than node_count.  First entry is for something
            // else.
            //       'block_map' is 1-based.
            const auto &        block_map = block.second;
            std::vector<double> x(block_map->size());
            std::vector<double> y(block_map->size());
            std::vector<double> z(block_map->size());

            for (size_t i = 0; i < block_map->size(); i++) {
              auto global = block_map->map()[i + 1] - 1;
              x[i]        = rdata[global * spatial_dim + 0];
              if (spatial_dim > 1) {
                y[i] = rdata[global * spatial_dim + 1];
              }
              if (spatial_dim > 2) {
                z[i] = rdata[global * spatial_dim + 2];
              }
            }

            // Create the zone
            // Output this zones coordinates...
            int crd_idx = 0;
            CGCHECK(cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateX", TOPTR(x),
                                   &crd_idx));

            if (spatial_dim > 1) {
              CGCHECK(cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateY",
                                     TOPTR(y), &crd_idx));
            }

            if (spatial_dim > 2) {
              CGCHECK(cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, "CoordinateZ",
                                     TOPTR(z), &crd_idx));
            }
          }
        }
        else {
          // Outputting only a single coordinate value...
          for (const auto &block : m_globalToBlockLocalNodeMap) {
            auto zone = block.first;
            // NOTE: 'block_map' has one more entry than node_count.  First entry is for something
            // else.
            //       'block_map' is 1-based.
            const auto &        block_map = block.second;
            std::vector<double> xyz(block_map->size());

            for (size_t i = 0; i < block_map->size(); i++) {
              auto global = block_map->map()[i + 1] - 1;
              xyz[i]      = rdata[global];
            }

            std::string cgns_name = "Invalid";
            if (field.get_name() == "mesh_model_coordinates_x") {
              cgns_name = "CoordinateX";
            }
            else if (field.get_name() == "mesh_model_coordinates_y") {
              cgns_name = "CoordinateY";
            }
            else if (field.get_name() == "mesh_model_coordinates_z") {
              cgns_name = "CoordinateZ";
            }
            // Create the zone
            // Output this zones coordinates...
            int crd_idx = 0;
            CGCHECK(cg_coord_write(get_file_pointer(), base, zone, CG_RealDouble, cgns_name.c_str(),
                                   TOPTR(xyz), &crd_idx));
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "output");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      double *rdata      = static_cast<double *>(data);
      int     cgns_field = 0;

      for (const auto &block : m_globalToBlockLocalNodeMap) {
        auto zone = block.first;
        // NOTE: 'block_map' has one more entry than node_count.
        // First entry is for something else.  'block_map' is
        // 1-based.
        const auto &        block_map = block.second;
        std::vector<double> blk_data(block_map->size());

        auto var_type   = field.transformed_storage();
        int  comp_count = var_type->component_count();

        if (comp_count == 1) {
          for (size_t j = 0; j < block_map->size(); j++) {
            auto global = block_map->map()[j + 1] - 1;
            blk_data[j] = rdata[global];
          }
          CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentVertexSolutionIndex,
                                 CG_RealDouble, field.get_name().c_str(), blk_data.data(),
                                 &cgns_field));
          Utils::set_field_index(field, cgns_field, CG_Vertex);
        }
        else {
          char field_suffix_separator = get_field_separator();

          for (int i = 0; i < comp_count; i++) {
            for (size_t j = 0; j < block_map->size(); j++) {
              auto global = block_map->map()[j + 1] - 1;
              blk_data[j] = rdata[comp_count * global + i];
            }
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);
            CGCHECK(cg_field_write(get_file_pointer(), base, zone, m_currentVertexSolutionIndex,
                                   CG_RealDouble, var_name.c_str(), blk_data.data(), &cgns_field));
            if (i == 0) {
              Utils::set_field_index(field, cgns_field, CG_Vertex);
            }
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(nb, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ns, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    const Ioss::EntityBlock *parent_block = sb->parent_block();
    if (parent_block == nullptr) {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: SideBlock " << sb->name()
             << " does not have a parent-block specified.  This is required for CGNS output.";
      IOSS_ERROR(errmsg);
    }

    int     base       = parent_block->get_property("base").get_int();
    int     zone       = Iocgns::Utils::get_db_zone(parent_block);
    ssize_t num_to_get = field.verify(data_size);

    Ioss::Field::RoleType role = field.get_role();

    if (role == Ioss::Field::MESH) {
      // Handle the MESH fields required for a CGNS file model.
      // (The 'genesis' portion)
      if (field.get_name() == "element_side") {
        // Get name from parent sideset...
        auto &name = sb->owner()->name();

        CG_ElementType_t type = Utils::map_topology_to_cgns(sb->topology()->name());
        int              sect = 0;

        cgsize_t cg_start = m_bcOffset[zone] + 1;
        cgsize_t cg_end   = m_bcOffset[zone] + num_to_get;
        m_bcOffset[zone] += num_to_get;

        // NOTE: Currently not writing the "ElementConnectivity" data for the
        //       boundary condition.  It isn't used in the read and don't have
        //       the data so would have to generate it.  This may cause problems
        //       with codes that use the downstream data if they base the BC off
        //       of the nodes instead of the element/side info.
        CGCHECK(cg_section_partial_write(get_file_pointer(), base, zone, name.c_str(), type, cg_start,
                                         cg_end, 0, &sect));

        sb->property_update("section", sect);

        size_t                offset = m_zoneOffset[zone - 1];
        std::vector<cgsize_t> parent(4 * num_to_get);

        if (field.get_type() == Ioss::Field::INT32) {
          int *  idata = reinterpret_cast<int *>(data);
          size_t j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            parent[num_to_get * 0 + i] = idata[j++] - offset; // Element
            parent[num_to_get * 2 + i] = idata[j++];
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_ioss_face_to_cgns(sb->parent_element_topology(), num_to_get, parent);
        }
        else {
          int64_t *idata = reinterpret_cast<int64_t *>(data);
          size_t   j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            parent[num_to_get * 0 + i] = idata[j++] - offset; // Element
            parent[num_to_get * 2 + i] = idata[j++];
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_ioss_face_to_cgns(sb->parent_element_topology(), num_to_get, parent);
        }

        CGCHECK(cg_parent_data_write(get_file_pointer(), base, zone, sect, TOPTR(parent)));
        return num_to_get;
      }

      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ss, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(cs, field, "output");
  }

  void DatabaseIO::write_results_meta_data() {}

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::STRUCTUREDBLOCK | Ioss::NODESET |
           Ioss::SIDESET | Ioss::REGION;
  }

} // namespace Iocgns
