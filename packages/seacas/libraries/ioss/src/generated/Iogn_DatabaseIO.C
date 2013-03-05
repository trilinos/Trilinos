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

#include <generated/Iogn_DatabaseIO.h>
#include <generated/Iogn_GeneratedMesh.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <string>

namespace {
  template <typename INT>
    void map_global_to_local(const Ioss::Map &map, size_t count, size_t stride, INT *data)
  {
    for (size_t i=0; i < count; i+=stride) {
      int64_t local = map.global_to_local(data[i], true);
      data[i] = local;
    }
  }

  template <typename INT>
  void fill_transient_data(const Ioss::GroupingEntity *entity, size_t component_count, double *data, INT /*dummy*/) {
    std::vector<INT> ids;
    entity->get_field_data("ids", ids);

    double *rdata = static_cast<double*>(data);
    if (component_count == 1) {
      for (size_t i=0; i < ids.size(); i++) {
        rdata[i] = sqrt(ids[i]);
      }
    } else {
      for (size_t i=0; i < ids.size(); i++) {
        for (size_t j=0; j < component_count; j++) {
          rdata[i*component_count + j] = j+sqrt(ids[i]);
        }
      }
    }
  }

  void fill_transient_data(const Ioss::GroupingEntity *entity, const Ioss::Field &field, void *data) {
    const Ioss::Field &ids = entity->get_fieldref("ids");
    if (ids.is_type(Ioss::Field::INTEGER)) {
      fill_transient_data(entity, field.raw_storage()->component_count(), (double*)data, (int)0);
    } else {
      fill_transient_data(entity, field.raw_storage()->component_count(), (double*)data, (int64_t)0);
    }
  }

  void fill_constant_data(const Ioss::Field &field, void *data, double value) {
    double *rdata = (double*)data;
    size_t count = field.raw_count();
    size_t component_count = field.raw_storage()->component_count();
    for (size_t i=0; i < count*component_count; i++) {
      rdata[i] = value;
    }
  }
}
namespace Iogn {

  // ========================================================================
  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
  : Ioss::IOFactory("generated")
  {}

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
                                       Ioss::DatabaseUsage db_usage,
                                       MPI_Comm communicator,
                                       const Ioss::PropertyManager &props) const
  { return new DatabaseIO(NULL, filename, db_usage, communicator, props); }

  // ========================================================================
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage,
                         MPI_Comm communicator,
                         const Ioss::PropertyManager &props) :
        Ioss::DatabaseIO(region, filename, db_usage, communicator, props),
        m_generatedMesh(NULL),  spatialDimension(3), nodeCount(0),
        elementCount(0), nodeBlockCount(0),
        elementBlockCount(0), nodesetCount(0), sidesetCount(0),
        nodeMap("node"), elemMap("elem"), m_useVariableDf(true)
  {
    if (is_input()) {
      dbState = Ioss::STATE_UNKNOWN;
    } else {
      std::ostringstream errmsg;
      errmsg << "Generated mesh option is only valid for input mesh.";
      IOSS_ERROR(errmsg);
    }
    if (props.exists("USE_CONSTANT_DF")) {
      m_useVariableDf = false;
    }
  }

  DatabaseIO::~DatabaseIO()
  {
    delete m_generatedMesh;
  }

  void DatabaseIO::release_memory()
  {
    nodeMap.release_memory();
    elemMap.release_memory();
  }

  void DatabaseIO::read_meta_data()
  {
    if(m_generatedMesh == NULL)
    {
      m_generatedMesh = new GeneratedMesh(get_filename(), util().parallel_size(), util().parallel_rank());
    }

    spatialDimension = 3;
    nodeCount = m_generatedMesh->node_count_proc();
    elementCount = m_generatedMesh->element_count_proc();
    elementBlockCount =  m_generatedMesh->block_count();
    nodesetCount = m_generatedMesh->nodeset_count();
    sidesetCount = m_generatedMesh->sideset_count();

    get_step_times();

    get_nodeblocks();
    get_elemblocks();
    get_nodesets();
    get_sidesets();
    get_commsets();

    Ioss::Region *this_region = get_region();
    this_region->property_add(Ioss::Property(std::string("title"), std::string("GeneratedMesh: ") += get_filename()));
    this_region->property_add(Ioss::Property(std::string("spatial_dimension"), 3));
  }

  bool DatabaseIO::begin(Ioss::State /* state */)
  {
    return true;
  }

  bool   DatabaseIO::end(Ioss::State /* state */)
  {
    return true;
  }

  bool DatabaseIO::begin_state(Ioss::Region *region, int /* state */, double time )
  {
    return true;
  }

  bool   DatabaseIO::end_state(Ioss::Region */* region */, int /* state */, double /* time */)
  {
    return true;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb,
                                         const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates") {
        // Cast 'data' to correct size -- double
        double *rdata = static_cast<double*>(data);
        m_generatedMesh->coordinates(rdata);
      }

      else if (field.get_name() == "ids") {
        // Map the local ids in this node block
        // (1...node_count) to global node ids.
        get_node_map().map_implicit_data(data, field, num_to_get, 0);
      }
      else if (field.get_name() == "owning_processor") {
        int *owner = static_cast<int*>(data);
        m_generatedMesh->owning_processor(owner, num_to_get);
      }
      else if (field.get_name() == "connectivity") {
        // Do nothing, just handles an idiosyncracy of the GroupingEntity
      }
      else if (field.get_name() == "connectivity_raw") {
        // Do nothing, just handles an idiosyncracy of the GroupingEntity
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "input");
      }
      return num_to_get;
    } else {
      fill_transient_data(nb, field, data);
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region* /* reg */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    int64_t id = eb->get_property("id").get_int();
    int64_t element_count = eb->get_property("entity_count").get_int();
    Ioss::Field::RoleType role = field.get_role();

    if (role == Ioss::Field::MESH) {
      // Handle the MESH fields required for an ExodusII file model.
      // (The 'genesis' portion)

      if (field.get_name() == "connectivity" || field.get_name() == "connectivity_raw") {
        assert(field.raw_storage()->component_count() == m_generatedMesh->topology_type(id).second);

        // The generated mesh connectivity is returned in a vector.  Ids are global
        if (field.is_type(Ioss::Field::INTEGER)) {
          int *connect = static_cast<int*>(data);
          m_generatedMesh->connectivity(id, connect);
          if (field.get_name() == "connectivity_raw") {
            map_global_to_local(get_node_map(), element_count * field.raw_storage()->component_count(), 1, connect);
          }
        } else {
          int64_t *connect = static_cast<int64_t*>(data);
          m_generatedMesh->connectivity(id, connect);
          if (field.get_name() == "connectivity_raw") {
            map_global_to_local(get_node_map(), element_count * field.raw_storage()->component_count(), 1, connect);
          }
        }
      }
      else if (field.get_name() == "ids") {
        // Map the local ids in this element block
        // (eb_offset+1...eb_offset+1+element_count) to global element ids.
        get_element_map().map_implicit_data(data, field, num_to_get, eb->get_offset());
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "input");
      }
    }

    else if (role == Ioss::Field::ATTRIBUTE) {
      if (element_count > 0) {
        int attribute_count = eb->get_property("attribute_count").get_int();
        if (attribute_count > 0) {
          double *attr = static_cast<double*>(data);
          for (size_t i=0; i < num_to_get; i++) {
            attr[i] = 1.0;
          }
        }
      }
    }

    else if (role == Ioss::Field::TRANSIENT) {
      // Fill the field with arbitrary data...
      fill_transient_data(eb, field, data);
    } else if (role == Ioss::Field::REDUCTION) {
      num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* ef_blk,
                                         const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    int64_t id = ef_blk->get_property("id").get_int();
    size_t entity_count = ef_blk->get_property("entity_count").get_int();
    if (num_to_get != entity_count) {
      std::ostringstream errmsg;
      errmsg << "Partial field input not implemented for side blocks";
      IOSS_ERROR(errmsg);
    }

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {

      if (field.get_name() == "ids") {
        // A sideset' is basically an exodus sideset.  A
        // sideset has a list of elements and a corresponding local
        // element side (1-based) The side id is: side_id =
        // 10*element_id + local_side_number This assumes that all
        // sides in a sideset are boundary sides.
        std::vector<int64_t> elem_side;
        m_generatedMesh->sideset_elem_sides(id, elem_side);
        if (field.is_type(Ioss::Field::INTEGER)) {
          int *ids = static_cast<int*>(data);
          for (size_t i=0; i < num_to_get; i++) {
            ids[i] = 10 * elem_side[2*i+0] + elem_side[2*i+1] + 1;
          }
        } else {
          int64_t *ids = static_cast<int64_t*>(data);
          for (size_t i=0; i < num_to_get; i++) {
            ids[i] = 10 * elem_side[2*i+0] + elem_side[2*i+1] + 1;
          }
        }
      }

      else if (field.get_name() == "element_side" ||
          field.get_name() == "element_side_raw") {
        // Since we only have a single array, we need to allocate an extra
        // array to store all of the data.  Note also that the element_id
        // is the global id but only the local id is stored so we need to
        // map from local_to_global prior to generating the side id...

        std::vector<int64_t> elem_side;
        m_generatedMesh->sideset_elem_sides(id, elem_side);
        if (field.get_name() == "element_side_raw") {
          map_global_to_local(get_element_map(), elem_side.size(), 2, &elem_side[0]);
        }

        if (field.is_type(Ioss::Field::INTEGER)) {
          int *element_side = static_cast<int*>(data);
          for (size_t i=0; i < num_to_get; i++) {
            element_side[2*i+0] = elem_side[2*i+0];
            element_side[2*i+1] = elem_side[2*i+1] + 1;
          }
        } else {
          int64_t *element_side = static_cast<int64_t*>(data);
          for (size_t i=0; i < num_to_get; i++) {
            element_side[2*i+0] = elem_side[2*i+0];
            element_side[2*i+1] = elem_side[2*i+1] + 1;
          }
        }
      }

      else if (field.get_name() == "distribution_factors") {
	if (m_useVariableDf) {
	  fill_transient_data(ef_blk, field, data);
	} else {
	  fill_constant_data(field, data, 1.0);
	}
      }

      else {
        num_to_get = Ioss::Utils::field_warning(ef_blk, field, "input");
      }
    } else if (role == Ioss::Field::TRANSIENT) {
      fill_transient_data(ef_blk, field, data);
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
                                         const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    int64_t id = ns->get_property("id").get_int();
    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {

      if (field.get_name() == "ids" ||
          field.get_name() == "ids_raw") {
        std::vector<int64_t> nodes;
        m_generatedMesh->nodeset_nodes(id, nodes);
        if (field.get_name() == "ids_raw") {
          map_global_to_local(get_node_map(), nodes.size(), 1, &nodes[0]);
        }

        if (field.is_type(Ioss::Field::INTEGER)) {
          int *ids = static_cast<int*>(data);
          std::copy(nodes.begin(), nodes.end(), ids);
        } else {
          int64_t *ids = static_cast<int64_t*>(data);
          std::copy(nodes.begin(), nodes.end(), ids);
        }
      } else if (field.get_name() == "distribution_factors") {
	if (m_useVariableDf) {
	  fill_transient_data(ns, field, data);
	} else {
	  fill_constant_data(field, data, 1.0);
	}
      } else {
        num_to_get = Ioss::Utils::field_warning(ns, field, "input");
      }
    } else if (role == Ioss::Field::TRANSIENT) {
      fill_transient_data(ns, field, data);
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
                                         void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    size_t entity_count = cs->get_property("entity_count").get_int();

    // Return the <entity (node or face), processor> pair
    if (field.get_name() == "entity_processor") {

      // Check type -- node or face
      std::string type = cs->get_property("entity_type").get_string();

      if (type == "node") {
        // Allocate temporary storage space
        Ioss::Int64Vector entities(num_to_get);
        Ioss::IntVector procs(num_to_get);
        m_generatedMesh->node_communication_map(entities, procs);

        // and store in 'data' ...
        if (field.is_type(Ioss::Field::INTEGER)) {
          int* entity_proc = static_cast<int*>(data);

          size_t j=0;
          for (size_t i=0; i < entity_count; i++) {
            entity_proc[j++] = entities[i];
            entity_proc[j++] = procs[i];
          }
        } else {
          int64_t* entity_proc = static_cast<int64_t*>(data);

          size_t j=0;
          for (size_t i=0; i < entity_count; i++) {
            entity_proc[j++] = entities[i];
            entity_proc[j++] = procs[i];
          }
        }
      } else {
        std::ostringstream errmsg;
        errmsg << "Invalid commset type " << type;
        IOSS_ERROR(errmsg);
      }

    } else if (field.get_name() == "ids") {
      // Do nothing, just handles an idiosyncracy of the GroupingEntity
    } else {
      num_to_get = Ioss::Utils::field_warning(cs, field, "input");
    }
    return num_to_get;
  }

  // Input only database -- these will never be called...
  int64_t DatabaseIO::put_field_internal(const Ioss::Region*,      const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock*,const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock*,const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock*,const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock*,   const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet*,     const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet*,     const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet*,     const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet*,     const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet*,     const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock*,   const Ioss::Field&, void*, size_t) const {return -1;}
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet*,     const Ioss::Field&, void*, size_t) const {return -1;}

  const Ioss::Map& DatabaseIO::get_node_map() const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.map.empty()) {
      nodeMap.map.resize(nodeCount+1);

      if (is_input()) {
        std::vector<int64_t> map;
        m_generatedMesh->node_map(map);

        // Map needed for Ioss starts at position 1 since the
        // sequential/non-sequential flag is at position 0...
        std::copy(map.begin(), map.end(), &nodeMap.map[1]);

        // Check for sequential node map.
        // If not, build the reverse G2L node map...
        nodeMap.map[0] = -1;
        for (int64_t i=1; i < nodeCount+1; i++) {
          if (i != nodeMap.map[i]) {
            nodeMap.map[0] = 1;
            break;
          }
        }

        nodeMap.build_reverse_map(myProcessor);

      } else {
        // Output database; nodeMap not set yet... Build a default map.
        for (int64_t i=1; i < nodeCount+1; i++) {
          nodeMap.map[i] = i;
        }
        // Sequential map
        nodeMap.map[0] = -1;
      }
    }
    return nodeMap;
  }

  const Ioss::Map& DatabaseIO::get_element_map() const
  {
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elemMap.map.empty()) {
      elemMap.map.resize(elementCount+1);

      if (is_input()) {
        std::vector<int64_t> map;
        m_generatedMesh->element_map(map);

        // Map needed for Ioss starts at position 1 since the
        // sequential/non-sequential flag is at position 0...
        std::copy(map.begin(), map.end(), &elemMap.map[1]);

        // Check for sequential element map.
        // If not, build the reverse G2L element map...
        elemMap.map[0] = -1;
        for (int64_t i=1; i < elementCount+1; i++) {
          if (i != elemMap.map[i]) {
            elemMap.map[0] = 1;
            break;
          }
        }

        if (elemMap.map[0] == 1) {
          elemMap.build_reverse_map(myProcessor);
        }

      } else {
        // Output database; elementMap not set yet... Build a default map.
        for (int64_t i=1; i < elementCount+1; i++) {
          elemMap.map[i] = i;
        }
        elemMap.map[0] = -1;
      }
    }
    return elemMap;
  }

  void DatabaseIO::get_nodeblocks()
  {
    std::string block_name = "nodeblock_1";
    Ioss::NodeBlock *block = new Ioss::NodeBlock(this, block_name,
                                                 m_generatedMesh->node_count_proc(), 3);
    block->property_add(Ioss::Property("id", 1));
    get_region()->add(block);
    add_transient_fields(block);
  }

  void DatabaseIO::get_step_times()
  {
    int time_step_count = m_generatedMesh->timestep_count();
    for (int i=0; i < time_step_count; i++) {
      get_region()->add_state(i);
    }
  }

  void DatabaseIO::get_elemblocks()
  {
    // Attributes of an element block are:
    // -- id
    // -- name
    // -- element type
    // -- number of elements
    // -- number of attributes per element
    // -- number of nodes per element (derivable from type)
    // -- number of faces per element (derivable from type)
    // -- number of edges per element (derivable from type)

    int block_count = m_generatedMesh->block_count();
    for (int i=0; i < block_count; i++) {
      std::string name = Ioss::Utils::encode_entity_name("block", i+1);
      std::string type = m_generatedMesh->topology_type(i+1).first;
      size_t element_count =  m_generatedMesh->element_count_proc(i+1);
      Ioss::ElementBlock *block = new Ioss::ElementBlock(this, name, type, element_count);

      block->property_add(Ioss::Property("id", i+1));

      // Maintain block order on output database...
      block->property_add(Ioss::Property("original_block_order", i));
      block->property_add(Ioss::Property("global_entity_count", (int64_t)m_generatedMesh->element_count(i+1)));

      if(type == "shell4")
      {
        block->field_add(Ioss::Field("thickness", Ioss::Field::REAL, "scalar", Ioss::Field::ATTRIBUTE,
                                     (int64_t)m_generatedMesh->element_count_proc(i+1)));
      }

      get_region()->add(block);
      add_transient_fields(block);
    }
  }

  void DatabaseIO::get_nodesets()
  {
    // Attributes of a nodeset are:
    // -- id
    // -- name
    // -- number of nodes
    // -- number of distribution factors (see next comment)
    // ----the #distribution factors should equal #nodes or 0, any
    //     other value does not make sense. If it is 0, then a substitute
    //     list will be created returning 1.0 for the factor

    // In a parallel execution, it is possible that a nodeset will have
    // no nodes or distribution factors on a particular processor...

    // Get nodeset metadata
    for (int ins = 0; ins < nodesetCount; ins++) {
      int64_t number_nodes = m_generatedMesh->nodeset_node_count_proc(ins+1);

      std::string name = Ioss::Utils::encode_entity_name("nodelist", ins+1);
      Ioss::NodeSet *nodeset = new Ioss::NodeSet(this, name, number_nodes);
      nodeset->property_add(Ioss::Property("id", ins+1));
      get_region()->add(nodeset);
      add_transient_fields(nodeset);
    }
  }

  void DatabaseIO::get_sidesets()
  {
    m_sideset_names.reserve(sidesetCount);
    for (int ifs = 0; ifs < sidesetCount; ifs++) {
      std::string name = Ioss::Utils::encode_entity_name("surface", ifs+1);
      m_sideset_names.push_back(name);
      Ioss::SideSet *sideset = new Ioss::SideSet(this, name);
      sideset->property_add(Ioss::Property("id", ifs+1));
      get_region()->add(sideset);

      std::string ef_block_name = name + "_quad4";
      std::string side_topo_name = "quad4";
      std::string elem_topo_name = "unknown";
      int64_t number_faces = m_generatedMesh->sideset_side_count_proc(ifs+1);

      Ioss::SideBlock *ef_block = new Ioss::SideBlock(this, ef_block_name,
                                                      side_topo_name,
                                                      elem_topo_name,
                                                      number_faces);
      sideset->add(ef_block);
      ef_block->property_add(Ioss::Property("id", ifs+1));

      std::string storage = "Real[";
      storage += Ioss::Utils::to_string(4);
      storage += "]";
      ef_block->field_add(Ioss::Field("distribution_factors",
                                      Ioss::Field::REAL, storage,
                                      Ioss::Field::MESH, number_faces));
    }
  }

  void DatabaseIO::get_commsets()
  {
    if (util().parallel_size() > 1) {
      // Get size of communication map...
      size_t my_node_count = m_generatedMesh->communication_node_count_proc();

      // Create a single node commset 
      Ioss::CommSet *commset = new Ioss::CommSet(this, "commset_node", "node", my_node_count);
      commset->property_add(Ioss::Property("id", 1));
      get_region()->add(commset);
    }
  }

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::REGION;
  }

  void DatabaseIO::add_transient_fields(Ioss::GroupingEntity *entity)
  {
    Ioss::EntityType type = entity->type();
    size_t entity_count = entity->get_property("entity_count").get_int();
    size_t var_count = m_generatedMesh->get_variable_count(type);
    for (size_t i=0; i < var_count; i++) {
      std::string var_name = entity->type_string() + "_" + Ioss::Utils::to_string(i+1);
      entity->field_add(Ioss::Field(var_name, Ioss::Field::REAL, "scalar", Ioss::Field::TRANSIENT, entity_count));
      std::cerr << "Adding field '" << var_name << "' to '" << entity->name() << "'.\n";
    }
  }
}
