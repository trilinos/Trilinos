// Copyright(C) 1999-2010
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

#ifndef SIERRA_Iopg_DatabaseIO_h
#define SIERRA_Iopg_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_DBUsage.h>

#include <string>
#include <assert.h>
#include <iostream>

namespace Ioss {
  class GroupingEntity;
  class Region;
  class EntityBlock;
  class NodeBlock;
  class SideBlock;
  class ElementBlock;
  class NodeSet;
  class SideSet;
  class CommSet;
}

namespace Iopg {
  class IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory* factory();
  private:
    IOFactory();
    Ioss::DatabaseIO* make_IO(const std::string& filename,
			      Ioss::DatabaseUsage db_usage,
			      MPI_Comm communicator) const;
  };

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string& filename,
	       Ioss::DatabaseUsage db_usage, MPI_Comm communicator);
    ~DatabaseIO();

    int node_global_to_local(int /* global */, bool /* must_exist */) const {return 0;}
    int element_global_to_local(int /* global */) const {return 0;}

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const {return 0;}

    void read_meta_data();

    bool begin(Ioss::State state);
    bool   end(Ioss::State state);

    bool begin_state(Ioss::Region *region, int state, double time);
    bool   end_state(Ioss::Region *region, int state, double time);

    std::string title()               const     {return databaseTitle;}
    int    spatial_dimension()   const     {return spatialDimension;}
    int    node_count()          const     {return nodeCount;}
    int    side_count()          const     {return 0;}
    int    element_count()       const     {return elementCount;}
    int    node_block_count()    const     {return nodeBlockCount;}
    int    element_block_count() const     {return elementBlockCount;}
    int    sideset_count()       const     {return sidesetCount;}
    int    nodeset_count()       const     {return nodesetCount;}
    int    maximum_symbol_length() const {return 32;}

    void get_block_adjacencies(const Ioss::ElementBlock *eb,
			       std::vector<std::string> &block_adjacency) const;

    void compute_block_membership(int id, std::vector<std::string> &block_membership) const;
    void compute_block_membership(Ioss::SideBlock *efblock,
				  std::vector<std::string> &block_membership) const;

  private:
    void read_region();
    void read_communication_metadata();
      
    void get_nodeblocks();
    void get_elemblocks();
    void get_nodesets();
    void get_sidesets();
    void get_commsets();
      
    int get_side_connectivity(const Ioss::SideBlock* fb, int id, int side_count,
			      int *fconnect, size_t data_size) const;
    int get_side_distributions(const Ioss::SideBlock* fb, int id,
			       int side_count, double *dist_fact, size_t data_size) const;

    const Ioss::MapContainer& get_node_map() const;
    const Ioss::MapContainer& get_element_map() const;
      
    void compute_block_adjacencies() const;

    int get_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int get_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int get_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int get_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int get_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			   void *data, size_t data_size) const;

    int put_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int put_field_internal(const Ioss::FaceBlock* nb, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int put_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::SideBlock* fb, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::EdgeSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int put_field_internal(const Ioss::FaceSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int put_field_internal(const Ioss::ElementSet* ns, const Ioss::Field& field,
			   void *data, size_t data_size) const {return 0;}
    int put_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
			   void *data, size_t data_size) const;
    int put_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			   void *data, size_t data_size) const;

    // Private member functions
    DatabaseIO(const DatabaseIO& from); // do not implement
    DatabaseIO& operator=(const DatabaseIO& from); // do not implement

    std::string databaseTitle;

    int spatialDimension;
    int nodeCount;
    int elementCount;

    int nodeBlockCount;
    int elementBlockCount;
    int nodesetCount;
    int sidesetCount;

    // Communication Set Data
    int *nodeCmapIds;
    int *nodeCmapNodeCnts;
    int *elemCmapIds;
    int *elemCmapElemCnts;
    int commsetNodeCount;
    int commsetElemCount;

    // MAPS -- Used to convert from local exodusII ids/names to Sierra
    // database global ids/names

    //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
    //               sierra side.   global = nodeMap[local]
    // nodeMap[0] contains: -1 if sequential, 0 if ordering unknown, 1
    // if nonsequential
    mutable Ioss::MapContainer        nodeMap;
    mutable Ioss::MapContainer        reorderNodeMap;
    mutable Ioss::ReverseMapContainer reverseNodeMap;
    // (local==global)

    //---Element Map -- Maps internal (1..NUMEL) ids to global ids used on the
    //               sierra side.   global = elementMap[local]
    // elementMap[0] contains: -1 if sequential, 0 if ordering unknown,
    // 1 if nonsequential
    mutable Ioss::MapContainer        elementMap;
    mutable Ioss::MapContainer        reorderElementMap;
    mutable Ioss::ReverseMapContainer reverseElementMap;

    mutable std::vector<std::vector<bool> > blockAdjacency;

    mutable bool sequentialNG2L; // true if reverse node map is sequential
    mutable bool sequentialEG2L; // true if reverse element map is
    // sequential (local==global)
    mutable bool blockAdjacenciesCalculated; // True if the lazy creation of
    // block adjacencies has been calculated.
  };
}
#endif // SIERRA_Iopg_DatabaseIO_h
