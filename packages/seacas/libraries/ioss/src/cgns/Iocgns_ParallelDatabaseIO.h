// Copyright(C) 2015
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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

#ifndef IOSS_Iocgns_ParallelDatabaseIO_h
#define IOSS_Iocgns_ParallelDatabaseIO_h

#include <Ioss_DBUsage.h>               // for DatabaseUsage
#include <Ioss_DatabaseIO.h>            // for DatabaseIO
#include <Ioss_IOFactory.h>             // for IOFactory
#include <Ioss_Map.h>                   // for Map
#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for int64_t
#include <iostream>                     // for ostream
#include <string>                       // for string
#include <memory>
#include <Ioss_State.h>                 // for State
#include <Ioss_CodeTypes.h>

#include <cgns/Iocgns_DecompositionData.h>

#include <cgnslib.h>

namespace Ioss { class CommSet; }
namespace Ioss { class EdgeBlock; }
namespace Ioss { class EdgeSet; }
namespace Ioss { class ElementBlock; }
namespace Ioss { class ElementSet; }
namespace Ioss { class ElementTopology; }
namespace Ioss { class FaceBlock; }
namespace Ioss { class FaceSet; }
namespace Ioss { class Field; }
namespace Ioss { class GroupingEntity; }
namespace Ioss { class NodeBlock; }
namespace Ioss { class NodeSet; }
namespace Ioss { class Region; }
namespace Ioss { class SideBlock; }
namespace Ioss { class SideSet; }
namespace Ioss { class EntityBlock; }

namespace Iocgns {
  
  class ParallelDatabaseIO : public Ioss::DatabaseIO
  {
  public:

    enum class entity_type {NODE, ELEM};

    ParallelDatabaseIO(Ioss::Region *region, const std::string& filename,
		       Ioss::DatabaseUsage db_usage,
		       MPI_Comm communicator,
		       const Ioss::PropertyManager &props);

    ~ParallelDatabaseIO();
    
    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const;

    int64_t node_global_to_local(int64_t global, bool must_exist) const;
    int64_t element_global_to_local(int64_t global) const;

    void release_memory() override;

    void openDatabase() const;
    void closeDatabase() const;

    bool begin(Ioss::State state);
    bool   end(Ioss::State state);

    bool begin_state(Ioss::Region *region, int state, double time);
    bool   end_state(Ioss::Region *region, int state, double time);

    // Metadata-related functions.
    void read_meta_data();

  private:
    int64_t get_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::SideBlock* sb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			       void *data, size_t data_size) const;

    int64_t put_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			       void *data, size_t data_size) const;
    int64_t put_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			       void *data, size_t data_size) const;

    // ID Mapping functions.
    const Ioss::Map& get_map(entity_type type) const;
    const Ioss::Map& get_map(Ioss::Map &entity_map,
			     int64_t entityCount,
			     int64_t file_offset, int64_t file_count,
			     entity_type type) const;

    // Bulk Data

    // MAPS -- Used to convert from local exodusII ids/names to Sierra
    // database global ids/names

    //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
    //               sierra side.   global = nodeMap[local]
    // nodeMap[0] contains: -1 if sequential, 0 if ordering unknown, 1
    // if nonsequential

    mutable Ioss::Map nodeMap;
    mutable Ioss::Map elemMap;

    mutable int cgnsFilePtr;
    size_t nodeCount;
    size_t elementCount;
    
    mutable std::unique_ptr<DecompositionDataBase> decomp;
    
    std::vector<size_t> m_zoneOffset; // Offset for local zone/block element ids to global.
    std::vector<std::vector<cgsize_t>> m_blockLocalNodeMap;
    std::map<std::string, int> m_zoneNameMap;
  };
}
#endif
