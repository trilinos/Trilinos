// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "faodel/Iofaodel_DatabaseIO.h"
#include "faodel/Iofaodel_FieldSerialization.h"
#include "faodel/Iofaodel_PropertySerialization.h"
#include "faodel/Iofaodel_Utils.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_CommSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_SubSystem.h"
#include "Ioss_Utils.h"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "faodel-services/MPISyncStart.hh"

#if 1
#define PRINT_KEYS
#endif

namespace Iofaodel {
  // namespace {
  // Output a message that the operation is unsupported and die...
  void unsupported(const char *operation)
  {
    std::cerr << "ERROR: Unsupported functionality called: " << operation << '\n';
    std::abort();
  }

  int get_file_pointer() { return 0; }

  const char *Version() { return "Iofaodel_DatabaseIO.C 2010/09/22"; }

  void faodel_error(int exoid, int lineno, int /* processor */)
  {
    std::ostringstream errmsg;

    errmsg << "Faodel error at line " << lineno << " in file '" << Version()
           << "' Please report to gdsjaar@sandia.gov if you need help.";

    IOSS_ERROR(errmsg);
  }
  //} // namespace

  // ========================================================================
  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("faodel") {}

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       Ioss_MPI_Comm                communicator,
                                       const Ioss::PropertyManager &properties) const
  {
    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  /*
   * There are three basic cases to address, manifested in whether any bootstrapping is done and
   * what gets provided to kelpie::Connect()
   *
   *  - Each process wants a local store, potentially backed by an IOM, but isn't interested
   *    in sharing data with any other nodes. This will be satisfied by a localkv: URI. Dirman
   *    isn't strictly necessary.
   *
   *  - Some subset (an improper subset; all nodes might join) of the processes in the IOSS job will
   *    participate in keyspace partitioning (and therefore data sharing).
   *    This will be indicated by dht: or rft: URIs (distributed hash
   *    table or rank-folding table, respectively), and the configuration should specify
   *    which processes in the job are participants. URIs may be annotated to configure IOMs on
   *    each kelpie instance and to configure the pools themselves. Dirman must be started among
   *    the pool processes.
   *
   *  - No Kelpie pool will be hosted by this job. We need the node ID of a remote Dirman so that
   *    this job's calls to Connect() can have their URIs translated appropriately to the node IDs
   *    of the processes comprising the remote Kelpie service. We don't have any input into the
   *    configuration of that service.
   */

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, Ioss_MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props), spatialDimension(3),
        nodeBlockCount(0), elementBlockCount(0), nodesetCount(0), sidesetCount(0),
        commsetNodeCount(0), commsetElemCount(0)
  {
    std::string resource("local:/ioss");

    // Allow an override of the faodel configuration.
    // Two methods: supply the faodel-config-replace property, or set the FAODEL_CONFIG
    // environment variable. FAODEL_CONFIG takes precedence.
    // if FAODEL_CONFIG is set, we will also use (if set) FAODEL_RESOURCE_URI and not attempt
    // to do anything smart to it
    char *faodel_config_env = std::getenv("FAODEL_CONFIG");
    if (faodel_config_env not_eq nullptr) {

      faodel_config.Append(faodel_config_env);
      char *resource_override = std::getenv("FAODEL_RESOURCE_URI");

      if (resource_override not_eq nullptr)
        resource = resource_override;
    }
    else {

      if (props.exists("faodel-resource-uri"))
        resource = props.get("faodel-resource-uri").get_string();

      if (props.exists("faodel-config-replace")) {

        faodel_config.Append(props.get("faodel-config-replace").get_string());
      }
      else {

        // Do some things in common for all setups
        faodel_config.Append(R"(
lunasa.lazy_memory_manager malloc
lunasa.eager_memory_manager malloc
kelpie.ioms iom-posix
kelpie.iom.iom-posix.type posixindividualobjects
x)");

        // if there are other fragments the caller wants to add to the Faodel
        //   config, get it before we potentially stand up dirman
        faodel_config.AppendFromReferences();
        if (props.exists("faodel-config-append")) {
          faodel_config.Append(props.get("faodel-config-append").get_string());
        }

        // Set up IOM if the caller provided a path
        faodel_config.Append("kelpie.iom.iom-posix.path " + [props]() -> std::string {
          if (props.exists("faodel-iom-path")) {
            return props.get("faodel-iom-path").get_string();
          }
          else {
            return "/tmp";
          }
        }());
        resource.append("&iom=iom-posix&behavior=writetoiom");

        //
        // Figure out the setup based on the resource URI
        //
        if (resource.rfind("local:", 0) == 0) {

          // Asking for local kv, no syncstart, no dirman
          faodel_config.Append(R"(
kelpie.type nonet
mpisyncstart.enable false
)");
        }
        else {

          // We aren't local, so we can call syncstart
          //  If we were supplied a kelpie-dirman-nodeid property, we assume we're
          //  being pointed at a dirman that has already been set up. This should be
          //  the job-to-job case, where the dirman root runs in a separate application,
          //  but it's possible for the app to have set up kelpie/dirman itself and wants
          //  to use that.
          //
          //  If that property isn't present, we assume that we need to configure a dirman
          //  for our own use. Right now we assume any such dirman will be rooted at rank 0
          //  and set up a pool that spans all ranks.

          faodel_config.Append(R"(
kelpie.type standard
mpisyncstart.enable true
)");
          if (props.exists("faodel-dirman-nodeid")) {
            faodel_config.Append("dirman.root_node " +
                                 props.get("faodel-dirman-nodeid").get_string());
          }
          else {
            faodel_config.Append("dirman.root_node_mpi 0");
            faodel_config.Append("dirman.type centralized");
            faodel_config.Append("dirman.resources_mpi[] " + resource + " all");
          }

          if (resource.rfind("rft:", 0) == 0) {
            // any necessary RFTn URI decorations
            // Use the rank target if given, else use our rank as a target
            resource.append("&rank=");
            if (props.exists("faodel-rft-target-rank")) {
              resource.append(std::to_string(props.get("faodel-rft-target-rank").get_int()));
            }
            else {
              resource.append(std::to_string(parallel_rank()));
            }
          }
        }
      }
    }

    std::once_flag faodel_bootstrap_flag;
    std::call_once(faodel_bootstrap_flag, [this]() {
      faodel::mpisyncstart::bootstrap();
      try {
        faodel::bootstrap::Start(faodel_config, kelpie::bootstrap);
      }
      catch (const std::exception &xc) {
        std::cerr << "bootstrap exception" << xc.what();
        std::rethrow_exception(std::current_exception());
      }
    });

    // Use the filename to create a bucket in the resource URI
    // This will help to isolate multiple datasets being used at the same URI
    resource.append("&bucket=") += filename;

    try {
      pool = kelpie::Connect(resource);
    }
    catch (const std::exception &xc) {
      std::cerr << "kelpie::Connect exception: " << xc.what();
      std::rethrow_exception(std::current_exception());
    }
    dbState = Ioss::STATE_UNKNOWN;
    std::cerr << "exiting Iofaodel ctor" << std::endl;
  }

  DatabaseIO::~DatabaseIO()
  {
    // Shut down Faodel gracefully
    faodel::bootstrap::Finish();
  }

  bool DatabaseIO::put_properties() const
  {
    // TODO add check to see what's been published before publishing again
    map_properties(*(get_region()), [this](const Ioss::Region &r, const Ioss::GroupingEntity &e,
                                           const Ioss::Property &p) {
      this->pool.Publish(make_key(parallel_rank(), r, e, p), pack_property(r, e, p));
    });
    return true;
  }

  void DatabaseIO::finalize_database() const
  {
    if (this->usage() == Ioss::DatabaseUsage::WRITE_RESTART ||
        this->usage() == Ioss::DatabaseUsage::WRITE_RESULTS ||
        this->usage() == Ioss::DatabaseUsage::WRITE_HISTORY ||
        this->usage() == Ioss::DatabaseUsage::WRITE_HEARTBEAT) {

      // write states to LDO
      pool.Publish(make_states_key(parallel_rank(), *get_region()), pack_states(*get_region()));

      auto sidesets = get_region()->get_sidesets();
      for (auto sideset : sidesets) {
        auto sideblocks = sideset->get_side_blocks();
        for (auto sideblock : sideblocks) {
          auto sideblock_key =
              make_sideblock_key(parallel_rank(), *(get_region()), *sideset, *sideblock);

          /* PUBLISH the SideBlocks that the SideSet references. */
          auto ldo = pack_sideblock(*sideblock);
          pool.Publish(sideblock_key, ldo);
        }
      }

      auto structuredblocks = get_region()->get_structured_blocks();
      for (auto structuredblock : structuredblocks) {
        auto structuredblock_key =
            make_structuredblock_key(parallel_rank(), *(get_region()), *structuredblock);

        /* PUBLISH the attributes of StructuredBlock */
        auto ldo = pack_structuredblock(*structuredblock);
        pool.Publish(structuredblock_key, ldo);
      }

      // write properties to LDOs and publish
      this->put_properties();
    }
  }

  std::string DatabaseIO::get_format() const { return "faodel"; }

  bool DatabaseIO::begin_state_nl(int /* state */, double /* time */) { return false; }

  bool DatabaseIO::end_state_nl(int /* state */, double /* time */) { return false; }

  void DatabaseIO::read_meta_data_nl()
  {
    this->get_step_times_nl();

    this->read_region();

    this->get_edgeblocks();
    this->get_elemblocks();
    this->get_faceblocks();
    this->get_nodeblocks();
    this->get_structuredblocks();

    this->get_edgesets();
    this->get_elemsets();
    this->get_facesets();
    this->get_nodesets();
    this->get_sidesets();
    this->get_commsets();
  }

  void DatabaseIO::get_step_times_nl()
  {
    auto                     search_key = make_states_search_key(parallel_rank(), *get_region());
    kelpie::ObjectCapacities oc;
    pool.List(search_key, &oc);
    if (oc.keys.size() == 1) {
      lunasa::DataObject ldo;
      pool.Need(oc.keys[0], oc.capacities[0], &ldo);

      auto meta = static_cast<Iofaodel::meta_entry_t *>(ldo.GetMetaPtr());

      auto entry = static_cast<Iofaodel::state_entry_t *>(
          static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

      auto data = static_cast<Iofaodel::state_entry_t::basic_type *>(
          static_cast<void *>(entry->data + entry->value.offset));

      for (size_t state(1); state <= entry->count; state++)
        get_region()->add_state(data[state - 1]);
    }
    // TODO
    // else {
    // Report error of not having 1 set of time steps
    // }
  }

  void DatabaseIO::read_region()
  {

    {
      // Region Properties
      kelpie::ObjectCapacities oc;
      auto search_key = property_search_key(parallel_rank(), *(get_region()), *(get_region()));
      pool.List(search_key, &oc);
      this->read_entity_properties(oc, *(this->get_region()));
    }

    {
      // Region Fields
      kelpie::ObjectCapacities oc;
      auto search_key = field_search_key(parallel_rank(), *(get_region()), *(get_region()));
      pool.List(search_key, &oc);
      this->read_entity_fields(oc, *(this->get_region()));
    }

    {
      // Region TRANSIENT Fields
      kelpie::ObjectCapacities oc;
      auto search_key = field_search_key(parallel_rank(), 1, *(get_region()), *(get_region()));
      pool.List(search_key, &oc);
      this->read_entity_fields(oc, *(this->get_region()));
    }
  }

  void DatabaseIO::read_entity_properties(kelpie::ObjectCapacities oc, Ioss::GroupingEntity &entity)
  {
    // TODO do we need to update default properties upon construction?
    // Properties
    for (size_t i = 0; i < oc.keys.size(); i++) {
      lunasa::DataObject ldo;
      pool.Need(oc.keys[i], oc.capacities[i], &ldo);

      auto meta(static_cast<meta_entry_t *>(ldo.GetMetaPtr()));
      auto prop = static_cast<Iofaodel::property_entry_t *>(
          static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

      std::string property_name(prop->data + prop->name.offset, prop->name.size);

      auto value_ptr = static_cast<void *>(prop->data + prop->value.offset);
      if (prop->basic_type == Ioss::Property::BasicType::STRING) {
        std::string value(prop->data + prop->value.offset, prop->value.size);
        entity.property_update(property_name, value);
      }
      else if (prop->basic_type == Ioss::Property::BasicType::INTEGER) {
        entity.property_update(property_name, *(reinterpret_cast<int64_t *>(value_ptr)));
      }
      else if (prop->basic_type == Ioss::Property::BasicType::REAL) {
        entity.property_update(property_name, *(reinterpret_cast<double *>(value_ptr)));
      }
    }
  }

  Ioss::Property DatabaseIO::read_property(lunasa::DataObject &ldo)
  {
    // Properties
    auto meta(static_cast<meta_entry_t *>(ldo.GetMetaPtr()));
    auto prop = static_cast<Iofaodel::property_entry_t *>(
        static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

    std::string     property_name(prop->data + prop->name.offset, prop->name.size);
    Ioss::Property *property = nullptr;

    auto value_ptr = static_cast<void *>(prop->data + prop->value.offset);
    if (prop->basic_type == Ioss::Property::BasicType::STRING) {
      std::string value(prop->data + prop->value.offset, prop->value.size);
      property = new Ioss::Property(property_name, value);
    }
    else if (prop->basic_type == Ioss::Property::BasicType::INTEGER) {
      property = new Ioss::Property(property_name, *(reinterpret_cast<int64_t *>(value_ptr)));
    }
    else if (prop->basic_type == Ioss::Property::BasicType::REAL) {
      property = new Ioss::Property(property_name, *(reinterpret_cast<double *>(value_ptr)));
    }

    return *property;
  }

  void DatabaseIO::read_entity_fields(kelpie::ObjectCapacities oc, Ioss::GroupingEntity &entity)
  {
    // Fields
    for (size_t i = 0; i < oc.keys.size(); i++) {
      lunasa::DataObject ldo;

      pool.Need(oc.keys[i], oc.capacities[i], &ldo);

      auto meta(static_cast<meta_entry_t *>(ldo.GetMetaPtr()));

      auto field = static_cast<field_entry_t *>(
          static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta->value.offset));

      std::string field_name(field->data + field->name.offset, field->name.size);
      std::string field_storage(field->data + field->storage.offset, field->storage.size);

      if (!entity.field_exists(field_name)) {
        entity.field_add(Ioss::Field(field_name, field->basic_type, field_storage, field->role_type,
                                     field->raw_count));
      }
    }
  }

  void DatabaseIO::get_edgeblocks()
  {
    std::string              type_string("EdgeBlock");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_edge_type(false);
      std::string original_edge_type;

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }

          if (oc.keys[i].K2().find("original_edge_type") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            original_edge_type      = property_get_int(ldo);
            have_original_edge_type = true;
          }
        }
      }

      if (have_entity_count && have_original_edge_type) {
        auto block = new Ioss::EdgeBlock(this, entity_name, original_edge_type, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *block);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_elemblocks()
  {
    std::string              type_string("ElementBlock");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_topology_type(false);
      std::string original_topology_type;

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }

          if (oc.keys[i].K2().find("original_topology_type") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            original_topology_type      = property_get_string(ldo);
            have_original_topology_type = true;
          }
        }
      }

      if (have_entity_count && have_original_topology_type) {
        auto block =
            new Ioss::ElementBlock(this, entity_name, original_topology_type, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *block);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *block);

        // Add TRANSIENT fields that aren't created in the CTor
        auto field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
        kelpie::ObjectCapacities field_oc_debug;
        pool.List(field_search_debug, &field_oc_debug);
        this->read_entity_fields(field_oc_debug, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_faceblocks()
  {
    std::string              type_string("FaceBlock");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_topology_type(false);
      std::string original_topology_type;

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }

          if (oc.keys[i].K2().find("original_topology_type") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            original_topology_type      = property_get_string(ldo);
            have_original_topology_type = true;
          }
        }
      }

      if (have_entity_count && have_original_topology_type) {
        auto block = new Ioss::FaceBlock(this, entity_name, original_topology_type, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *block);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_nodeblocks()
  {
    std::string              type_string("NodeBlock");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);
      int64_t component_degree(0);

      for (size_t i = 0; i < oc.keys.size(); i++) {
        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }

          if (oc.keys[i].K2().find("component_degree") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            component_degree      = property_get_int(ldo);
            have_component_degree = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        // Creates a new NodeBlock with default Properties and Fields
        auto block = new Ioss::NodeBlock(this, entity_name, entity_count, component_degree);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *block);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *block);

        // Add TRANSIENT fields that aren't created in the CTor
        auto field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
        kelpie::ObjectCapacities field_oc_debug;
        pool.List(field_search_debug, &field_oc_debug);
        this->read_entity_fields(field_oc_debug, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_structuredblocks()
  {
    std::string              type_string("StructuredBlock");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      std::unordered_map<std::string, int> ctor_properties;
      std::list<std::string>               property_names = {
          "component_degree", "ni",        "nj",       "nk",       "ni_global",
          "nj_global",        "nk_global", "offset_i", "offset_j", "offset_k"};

      for (size_t i = 0; i < oc.keys.size(); i++) {
        for (auto property_name : property_names) {
          // NOTE: ordinary find won't work b/c the names of multiple properties share substrings,
          // e.g., "ni" and "ni_global" (a string find on "ni" will match both)
          if (oc.keys[i].K2().find(entity_name) != std::string::npos) {
            size_t position = oc.keys[i].K2().rfind(property_name);
            if (position != std::string::npos &&
                (position + property_name.size()) == oc.keys[i].K2().size()) {
              lunasa::DataObject ldo;
              pool.Need(oc.keys[i], oc.capacities[i], &ldo);
              ctor_properties[property_name] = property_get_int(ldo);
            }
          }
        }
      }

      if (ctor_properties.size() != property_names.size()) {
        fmt::print(
            stderr,
            "[ERROR] Unable to reconstruct StructuredBlock ({}): FOUND {} of {} properties\n",
            entity_name.c_str(), ctor_properties.size(), property_names.size());
        return;
      }

      auto block = new Ioss::StructuredBlock(
          this, entity_name, ctor_properties["component_degree"], ctor_properties["ni"],
          ctor_properties["nj"], ctor_properties["nk"], ctor_properties["offset_i"],
          ctor_properties["offset_j"], ctor_properties["offset_k"], ctor_properties["ni_global"],
          ctor_properties["nj_global"], ctor_properties["nk_global"]);

      // Add attributes that aren't created in the CTor
      auto attribute_key = make_structuredblock_key(parallel_rank(), *(get_region()), *block);
      lunasa::DataObject ldo;
      pool.Need(attribute_key, &ldo);
      unpack_structuredblock(ldo, *block);

      // Add Properties that aren't created in the CTor
      auto property_search = property_search_key(parallel_rank(), *(get_region()), *block);
      kelpie::ObjectCapacities property_oc;
      pool.List(property_search, &property_oc);
      this->read_entity_properties(property_oc, *block);

      // Add fields that aren't created in the CTor
      auto field_search = field_search_key(parallel_rank(), *(get_region()), *block);
      kelpie::ObjectCapacities field_oc;
      pool.List(field_search, &field_oc);
      this->read_entity_fields(field_oc, *block);

      // Add TRANSIENT fields that aren't created in the CTor
      auto field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
      kelpie::ObjectCapacities field_oc_debug;
      pool.List(field_search_debug, &field_oc_debug);
      this->read_entity_fields(field_oc_debug, *block);

      this->get_region()->add(block);
    }
  }

  void DatabaseIO::get_edgesets()
  {
    std::string              type_string("EdgeSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::EdgeSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *entity);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_elemsets()
  {
    std::string              type_string("ElementSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::ElementSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *entity);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_facesets()
  {
    std::string              type_string("FaceSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < oc.keys.size(); i++) {

        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::FaceSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *entity);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_nodesets()
  {
    std::string              type_string("NodeSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);

      for (size_t i = 0; i < oc.keys.size(); i++) {
        std::string entity_name_search = "/" + entity_name + "/";
        if (oc.keys[i].K2().find(entity_name_search) != std::string::npos) {

          if (oc.keys[i].K2().find("entity_count") != std::string::npos) {
            lunasa::DataObject ldo;
            pool.Need(oc.keys[i], oc.capacities[i], &ldo);
            entity_count      = property_get_int(ldo);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count) {
        auto entity = new Ioss::NodeSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        auto property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities property_oc;
        pool.List(property_search, &property_oc);
        this->read_entity_properties(property_oc, *entity);

        // Add fields that aren't created in the CTor
        auto field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        kelpie::ObjectCapacities field_oc;
        pool.List(field_search, &field_oc);
        this->read_entity_fields(field_oc, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_sidesets()
  {
    std::string              type_string("SideSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      auto entity = new Ioss::SideSet(this, entity_name);

      // Add Properties that aren't created in the CTor
      auto property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
      kelpie::ObjectCapacities property_oc;
      pool.List(property_search, &property_oc);
      this->read_entity_properties(property_oc, *entity);

      // Add fields that aren't created in the CTor
      auto field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
      kelpie::ObjectCapacities field_oc;
      pool.List(field_search, &field_oc);
      this->read_entity_fields(field_oc, *entity);

      // FIND the SideBlocks that this SideSet references
      auto sideblocks_key = sideblocks_search_key(parallel_rank(), *(get_region()), *entity);
      kelpie::ObjectCapacities sideblocks_search_oc;
      pool.List(sideblocks_key, &sideblocks_search_oc);

      for (size_t i = 0; i < sideblocks_search_oc.Size(); i++) {
        lunasa::DataObject ldo;
        pool.Need(sideblocks_search_oc.keys[i], sideblocks_search_oc.capacities[i], &ldo);
        int64_t entity_count = unpack_sideblocks(ldo);

        auto sideblock_name = get_entity_name(sideblocks_search_oc.keys[i], "SideBlock");
        auto property_key   = make_property_key(parallel_rank(), *(get_region()), "SideBlock",
                                                sideblock_name, "STRING", "topology_type");
        lunasa::DataObject property_ldo;
        pool.Need(property_key, &property_ldo);
        Ioss::Property topo_property = this->read_property(property_ldo);

        auto parent_property_key =
            make_property_key(parallel_rank(), *(get_region()), "SideBlock", sideblock_name,
                              "STRING", "parent_topology_type");
        lunasa::DataObject parent_property_ldo;
        pool.Need(parent_property_key, &parent_property_ldo);
        Ioss::Property parent_topo_property = this->read_property(parent_property_ldo);

        auto sideblock = new Ioss::SideBlock(this, sideblock_name, topo_property.get_string(),
                                             parent_topo_property.get_string(), entity_count);

        // Add Properties that aren't created in the CTor
        auto sideblock_property_search =
            property_search_key(parallel_rank(), *(get_region()), *sideblock);
        kelpie::ObjectCapacities sideblock_property_oc;
        pool.List(sideblock_property_search, &sideblock_property_oc);
        this->read_entity_properties(sideblock_property_oc, *sideblock);

        // Add fields that aren't created in the CTor
        auto sideblock_field_search =
            field_search_key(parallel_rank(), *(get_region()), *sideblock);
        kelpie::ObjectCapacities sideblock_field_oc;
        pool.List(sideblock_field_search, &sideblock_field_oc);
        this->read_entity_fields(sideblock_field_oc, *sideblock);

        entity->add(sideblock);
      }

      this->get_region()->add(entity);
    }
  }

  void DatabaseIO::get_commsets()
  {
    std::string              type_string("CommSet");
    kelpie::ObjectCapacities oc;
    auto search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    pool.List(search_key, &oc);

    auto entity_names = get_entity_names(oc.keys, type_string);
    for (auto entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_entity_type(false);
      std::string entity_type("");

      for (size_t i = 0; i < oc.keys.size(); i++) {
        if (oc.keys[i].K2().find(entity_name) != std::string::npos) {
          /* Check if we have all of the necessary attributes */
          if (have_entity_type && have_entity_count)
            break;

          if (!have_entity_count) {
            size_t position = oc.keys[i].K2().rfind("entity_count");
            if (position != std::string::npos) {
              lunasa::DataObject ldo;
              pool.Need(oc.keys[i], oc.capacities[i], &ldo);
              entity_count      = property_get_int(ldo);
              have_entity_count = true;
            }
          }

          if (!have_entity_type) {
            size_t position = oc.keys[i].K2().rfind("entity_type");
            if (position != std::string::npos) {
              lunasa::DataObject ldo;
              pool.Need(oc.keys[i], oc.capacities[i], &ldo);
              entity_type      = property_get_string(ldo);
              have_entity_type = true;
            }
          }
        }
      }

      if (!have_entity_type) {
        fmt::print(stderr, "[ERROR] Unable to reconstruct CommSet (entity_type NOT found)\n");
        return;
      }

      if (!have_entity_count) {
        fmt::print(stderr, "[ERROR] Unable to reconstruct CommSet (entity_count NOT found)\n");
        return;
      }

      auto commset = new Ioss::CommSet(this, entity_name, entity_type, entity_count);

      // Add Properties that aren't created in the CTor
      auto property_search = property_search_key(parallel_rank(), *(get_region()), *commset);
      kelpie::ObjectCapacities property_oc;
      pool.List(property_search, &property_oc);
      this->read_entity_properties(property_oc, *commset);

      // Add fields that aren't created in the CTor
      auto field_search = field_search_key(parallel_rank(), *(get_region()), *commset);
      kelpie::ObjectCapacities field_oc;
      pool.List(field_search, &field_oc);
      this->read_entity_fields(field_oc, *commset);

      // Add TRANSIENT fields that aren't created in the CTor
      auto field_search_transient = field_search_key(parallel_rank(), 1, *(get_region()), *commset);
      kelpie::ObjectCapacities field_oc_transient;
      pool.List(field_search_transient, &field_oc_transient);
      this->read_entity_fields(field_oc_transient, *commset);

      this->get_region()->add(commset);
    }
  }

  void DatabaseIO::read_communication_metadata() {}

  int64_t DatabaseIO::get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*reg, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*nb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*nb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*nb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*eb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*fb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*fs, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*cs, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    // Field::raw_count() should be the same for all mesh_model_coords* fields. It refers to
    // the number of types of type SCALAR or VECTOR_3D etc.
    // The number of double values is the raw_count() * Field::raw_storage().component_count()
    if (field.get_name() == "mesh_model_coordinates" &&
        field.get_role() == Ioss::Field::RoleType::MESH) {

      size_t  num_to_get          = field.verify(data_size);
      auto    component_data_size = field.get_size() / this->spatial_dimension();
      double *data_ptr            = static_cast<double *>(data);
      // The role for all mesh_model_coords should be Ioss::Field::RoleType::MESH
      auto role = field.get_role();

      auto                dim = this->spatial_dimension();
      std::vector<double> data_x(num_to_get);
      std::vector<double> data_y((dim > 1) ? num_to_get : 0);
      std::vector<double> data_z((dim > 2) ? num_to_get : 0);

      {
        auto field_x =
            Ioss::Field("mesh_model_coordinates_x", field.get_type(), "scalar", role, num_to_get);
        this->get_field_internal(*sb, field_x, Data(data_x), component_data_size);
      }

      if (dim > 1) {
        auto field_y =
            Ioss::Field("mesh_model_coordinates_y", field.get_type(), "scalar", role, num_to_get);
        this->get_field_internal(*sb, field_y, Data(data_y), component_data_size);
      }

      if (dim > 2) {
        auto field_z =
            Ioss::Field("mesh_model_coordinates_z", field.get_type(), "scalar", role, num_to_get);
        this->get_field_internal(*sb, field_z, Data(data_z), component_data_size);
      }

      size_t index(0);
      for (auto id(0); id < field.raw_count(); ++id) {
        data_ptr[index++] = data_x[id];
        if (dim > 1)
          data_ptr[index++] = data_y[id];
        if (dim > 2)
          data_ptr[index++] = data_z[id];
      }
    }
    else
      return get_field_internal(*sb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::Assembly *a, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return 0;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::Blob *b, const Ioss::Field &field, void *data,
                                         size_t data_size) const
  {
    return 0;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*reg, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*nb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*eb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*nb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*eb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*fb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ns, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*fs, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*cs, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*sb, field, data, data_size);
    ;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::Assembly *a, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return 0;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::Blob *b, const Ioss::Field &field, void *data,
                                         size_t data_size) const
  {
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f,
                                         void *data, size_t data_size) const
  {
    lunasa::DataObject ldo;
    kelpie::Key        k = make_key(parallel_rank(), *(get_region()), e, f);
    pool.Need(k, &ldo);

    /*
       unpack the LDO to retrieve the field data (user data) and set the output variables
       according to the field type information
       */

    auto meta_ptr(static_cast<meta_entry_t *>(ldo.GetMetaPtr()));

    auto field_ptr(static_cast<field_entry_t *>(
        static_cast<void *>(static_cast<char *>(ldo.GetDataPtr()) + meta_ptr->value.offset)));

    // TODO what other checks do we need here?
    if (data_size != field_ptr->value.size)
      return 1;

    std::memcpy(data, static_cast<void *>(field_ptr->data + field_ptr->value.offset),
                field_ptr->value.size);

    return 0;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f,
                                         void *data, size_t data_size) const
  {
    auto key = make_key(parallel_rank(), *(get_region()), e, f);
    auto ldo = pack_field(*(get_region()), e, f, data, data_size);
    pool.Publish(key, ldo);
    return 0;
  }

} // namespace Iofaodel
