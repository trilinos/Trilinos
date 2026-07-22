// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <tokenize.h>

#include "s3/Ios3_DatabaseIO.h"

#include "s3/Ios3_AwsHelpers.h"
#include "s3/Ios3_FieldSerialization.h"
#include "s3/Ios3_PropertySerialization.h"
#include "s3/Ios3_Utils.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_CommSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_SubSystem.h"
#include "Ioss_Utils.h"

namespace Ios3 {

  // constants copied from exodusII.h
  const std::size_t MAX_STR_LENGTH  = 32; // Maximum length of QA record
  const std::size_t MAX_LINE_LENGTH = 80; // Maximum length of an information record

  // ========================================================================
  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("s3") {}

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       Ioss_MPI_Comm                communicator,
                                       const Ioss::PropertyManager &properties) const
  {
    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, Ioss_MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props),
        bucket_name(Ios3::helpers::cleanBucketName(decoded_filename()))
  {
    Ioss::PropertyManager local_properties(props);
    Ios3::helpers::getPropertiesFromEnvVars(local_properties, util());
    Ios3::helpers::getParamsFromProperties(local_properties, helper_params);

    if (helper_params.enable_aws_tracing) {
      Ios3::helpers::print_params(helper_params);
    }

    helper_context = Ios3::helpers::createContext(helper_params);
    bool success   = Ios3::helpers::createBucket(helper_context, bucket_name);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    dbState = Ioss::STATE_UNKNOWN;
  }

  DatabaseIO::~DatabaseIO() { Ios3::helpers::destroyContext(helper_context); }

  bool DatabaseIO::put_properties() const
  {
    auto property_op = [this](const Ioss::Region &r, const Ioss::GroupingEntity &e,
                              const Ioss::Property &p) {
      key_t key   = make_key(parallel_rank(), r, e, p);
      auto  value = pack_property(r, e, p);
      return Ios3::helpers::putValue(helper_context, bucket_name, key.second, value);
    };

    // TODO add check to see what's been published before publishing again
    int num_failed = map_properties(*(get_region()), property_op);
    if (num_failed > 0) {
      return false;
    }

    return true;
  }

  void DatabaseIO::finalize_database() const
  {
    bool success = true;

    if (this->usage() == Ioss::DatabaseUsage::WRITE_RESTART ||
        this->usage() == Ioss::DatabaseUsage::WRITE_RESULTS ||
        this->usage() == Ioss::DatabaseUsage::WRITE_HISTORY ||
        this->usage() == Ioss::DatabaseUsage::WRITE_HEARTBEAT) {

      {
        key_t key   = make_states_key(parallel_rank(), *get_region());
        auto  value = pack_states(*get_region());
        success     = Ios3::helpers::putValue(helper_context, bucket_name, key.second, value);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
      }

      const auto &sidesets = get_region()->get_sidesets();
      for (const auto &sideset : sidesets) {
        const auto &sideblocks = sideset->get_side_blocks();
        for (const auto &sideblock : sideblocks) {
          key_t sideblock_key =
              make_sideblock_key(parallel_rank(), *(get_region()), *sideset, *sideblock);

          /* PUBLISH the SideBlocks that the SideSet references. */
          auto value = pack_sideblock(*sideblock);
          success =
              Ios3::helpers::putValue(helper_context, bucket_name, sideblock_key.second, value);
          if (!success) {
            IOSS_ERROR("S3 Operation Failed");
          }
        }
      }

      const auto &structuredblocks = get_region()->get_structured_blocks();
      for (const auto &structuredblock : structuredblocks) {
        key_t structuredblock_key =
            make_structuredblock_key(parallel_rank(), *(get_region()), *structuredblock);

        /* PUBLISH the attributes of StructuredBlock */
        auto value = pack_structuredblock(*structuredblock);
        success =
            Ios3::helpers::putValue(helper_context, bucket_name, structuredblock_key.second, value);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
      }

      this->put_properties();
    }
  }

  bool DatabaseIO::begin_nl(Ioss::State state)
  {
    dbState = state;
    return true;
  }

  bool DatabaseIO::end_nl(Ioss::State state)
  {
    // Transitioning out of state 'state'
    assert(state == dbState);
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      if (!is_input()) {
        write_meta_data(open_create_behavior());
      }
      break;
    default: // ignore everything else...
      break;
    }

    dbState = Ioss::STATE_UNKNOWN;

    return true;
  }

  void DatabaseIO::put_qa()
  {
    struct qa_element
    {
      char qa_record[4][MAX_STR_LENGTH + 1];
    };

    size_t num_qa_records = qaRecords.size() / 4;

    int                        num_to_put = sizeof(qa_element) * (num_qa_records + 1);
    std::vector<unsigned char> value(num_to_put);
    auto                       qa = reinterpret_cast<qa_element *>(value.data());
    {
      int j = 0;
      for (size_t i = 0; i < num_qa_records; i++) {
        Ioss::Utils::copy_string(qa[i].qa_record[0], qaRecords[j++], MAX_STR_LENGTH + 1);
        Ioss::Utils::copy_string(qa[i].qa_record[1], qaRecords[j++], MAX_STR_LENGTH + 1);
        Ioss::Utils::copy_string(qa[i].qa_record[2], qaRecords[j++], MAX_STR_LENGTH + 1);
        Ioss::Utils::copy_string(qa[i].qa_record[3], qaRecords[j++], MAX_STR_LENGTH + 1);
      }
    }

    Ioss::Utils::time_and_date(qa[num_qa_records].qa_record[3], qa[num_qa_records].qa_record[2],
                               MAX_STR_LENGTH);

    std::string codename = "unknown";
    std::string version  = "unknown";
    if (get_region()->property_exists("code_name")) {
      codename = get_region()->get_property("code_name").get_string();
    }
    if (get_region()->property_exists("code_version")) {
      version = get_region()->get_property("code_version").get_string();
    }
    Ioss::Utils::copy_string(qa[num_qa_records].qa_record[0], codename, MAX_STR_LENGTH + 1);
    Ioss::Utils::copy_string(qa[num_qa_records].qa_record[1], version, MAX_STR_LENGTH + 1);

    std::string qa_key{"::QA_Records"};
    bool        success = Ios3::helpers::putValue(helper_context, bucket_name, qa_key, value);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }
  }

  void DatabaseIO::put_info()
  {
    using info_line = char[MAX_LINE_LENGTH + 1];
    int total_lines = 0;

    // dump info records, include the product_registry
    // See if the input file was specified as a property on the database...
    std::vector<std::string> input_lines;
    if (get_region()->property_exists("input_file_name")) {
      std::string filename = get_region()->get_property("input_file_name").get_string();
      // Determine size of input file so can embed it in info records...
      Ioss::Utils::input_file(filename, &input_lines, MAX_LINE_LENGTH);
    }

    // Get configuration information for IOSS library.
    // Split into strings and remove empty lines...
    std::string config = Ioss::IOFactory::show_configuration();
    std::replace(std::begin(config), std::end(config), '\t', ' ');
    auto lines = Ioss::tokenize(config, "\n");
    lines.erase(std::remove_if(lines.begin(), lines.end(),
                               [](const std::string &line) { return line.empty(); }),
                lines.end());

    // See if the client added any "information_records"
    size_t info_rec_size = informationRecords.size();
    size_t in_lines      = input_lines.size();
    size_t qa_lines      = 1; // Platform info
    size_t config_lines  = lines.size();

    total_lines = in_lines + qa_lines + info_rec_size + config_lines;

    int                        num_to_put = total_lines * sizeof(info_line);
    std::vector<unsigned char> value(num_to_put);
    auto                       info = reinterpret_cast<info_line *>(value.data());

    int i = 0;
    Ioss::Utils::copy_string(info[i++], Ioss::Utils::platform_information(), MAX_LINE_LENGTH + 1);

    // Copy input file lines into 'info' array...
    for (size_t j = 0; j < input_lines.size(); j++, i++) {
      Ioss::Utils::copy_string(info[i], input_lines[j], MAX_LINE_LENGTH + 1);
    }

    // Copy "information_records" property data ...
    for (size_t j = 0; j < informationRecords.size(); j++, i++) {
      Ioss::Utils::copy_string(info[i], informationRecords[j], MAX_LINE_LENGTH + 1);
    }

    for (size_t j = 0; j < lines.size(); j++, i++) {
      Ioss::Utils::copy_string(info[i], lines[j], MAX_LINE_LENGTH + 1);
    }

    std::string info_key{"::Info_Records"};
    bool        success = Ios3::helpers::putValue(helper_context, bucket_name, info_key, value);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }
  }

  void DatabaseIO::write_meta_data(Ioss::IfDatabaseExistsBehavior behavior)
  {
    if (behavior != Ioss::DB_APPEND && behavior != Ioss::DB_MODIFY) {
      bool omit_qa = false;
      Ioss::Utils::check_set_bool_property(properties, "OMIT_QA_RECORDS", omit_qa);
      if (!omit_qa) {
        put_qa();
      }
      bool omit_info = false;
      Ioss::Utils::check_set_bool_property(properties, "OMIT_INFO_RECORDS", omit_info);
      if (!omit_info) {
        put_info();
      }
    }
  }

  void DatabaseIO::read_meta_data_nl()
  {
    this->get_step_times_nl();

    this->read_region();

    this->get_nodeblocks();
    this->get_edgeblocks();
    this->get_faceblocks();
    this->get_elemblocks();

    this->get_structuredblocks();

    this->get_sidesets();
    this->get_nodesets();
    this->get_edgesets();
    this->get_facesets();
    this->get_elemsets();

    this->get_commsets();

    // We don't yet support Assemblies or Blobs
    // this->get_assemblies();
    // this->get_blobs();
  }

  void DatabaseIO::get_step_times_nl()
  {
    bool success = true;

    key_t                    search_key = make_states_search_key(parallel_rank(), *get_region());
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (success && keys.size() == 1) {
      std::vector<unsigned char> value;
      success = Ios3::helpers::getValue(helper_context, bucket_name, keys[0], value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      auto entry = reinterpret_cast<Ios3::state_entry_t *>(value.data() + sizeof(meta_entry_t));

      auto data = reinterpret_cast<Ios3::state_entry_t::basic_type *>(
          reinterpret_cast<void *>(entry->data + entry->value.offset));

      for (size_t state(1); state <= entry->count; state++) {
        get_region()->add_state(data[state - 1]);
      }
    }
    else {
      IOSS_ERROR("S3 Operation Failed");
    }
  }

  void DatabaseIO::read_region()
  {
    {
      // Region Properties
      key_t search_key = property_search_key(parallel_rank(), *(get_region()), *(get_region()));
      std::vector<std::string> keys;
      bool success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      this->read_entity_properties(keys, *(this->get_region()));
    }

    {
      // Region Fields
      key_t search_key = field_search_key(parallel_rank(), *(get_region()), *(get_region()));
      std::vector<std::string> keys;
      bool success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      this->read_entity_fields(keys, *(this->get_region()));
    }

    {
      // Region TRANSIENT Fields
      key_t search_key = field_search_key(parallel_rank(), 1, *(get_region()), *(get_region()));
      std::vector<std::string> keys;
      bool success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      this->read_entity_fields(keys, *(this->get_region()));
    }

    // Get QA records from database and add to qaRecords...
    {
      struct qa_element
      {
        char qa_record[4][MAX_STR_LENGTH + 1];
      };

      std::vector<unsigned char> value;
      bool success = Ios3::helpers::getValue(helper_context, bucket_name, "::QA_Records", value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      int num_qa = value.size() / sizeof(qa_element);

      if (num_qa > 0) {
        auto qa = reinterpret_cast<qa_element *>(value.data());

        for (int i = 0; i < num_qa; i++) {
          add_qa_record(qa[i].qa_record[0], qa[i].qa_record[1], qa[i].qa_record[2],
                        qa[i].qa_record[3]);
        }
      }
    }

    // Get information records from database and add to informationRecords...
    {
      using info_line = char[MAX_STR_LENGTH + 1];

      std::vector<unsigned char> value;
      bool success = Ios3::helpers::getValue(helper_context, bucket_name, "::Info_Records", value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      int num_info = value.size() / sizeof(info_line);
      if (num_info > 0) {
        auto info = reinterpret_cast<info_line *>(value.data());
        for (int i = 0; i < num_info; i++) {
          add_information_record(info[i]);
        }
      }
    }
  }

  void DatabaseIO::read_entity_properties(std::vector<std::string> &keys,
                                          Ioss::GroupingEntity     &entity)
  {
    // TODO do we need to update default properties upon construction?
    // Properties
    for (size_t i = 0; i < keys.size(); i++) {
      std::vector<unsigned char> value;
      bool success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      auto prop = reinterpret_cast<Ios3::property_entry_t *>(value.data());

      std::string property_name(prop->data + prop->name.offset, prop->name.size);

      auto value_ptr = static_cast<void *>(prop->data + prop->value.offset);
      if (prop->basic_type == Ioss::Property::BasicType::STRING) {
        std::string property_value(prop->data + prop->value.offset, prop->value.size);
        entity.property_update(property_name, property_value);
      }
      else if (prop->basic_type == Ioss::Property::BasicType::INTEGER) {
        entity.property_update(property_name, *(reinterpret_cast<int64_t *>(value_ptr)));
      }
      else if (prop->basic_type == Ioss::Property::BasicType::REAL) {
        entity.property_update(property_name, *(reinterpret_cast<double *>(value_ptr)));
      }
    }
  }

  Ioss::Property DatabaseIO::read_property(std::vector<unsigned char> &value)
  {
    // Properties
    auto prop = reinterpret_cast<Ios3::property_entry_t *>(value.data());

    std::string     property_name(prop->data + prop->name.offset, prop->name.size);
    Ioss::Property *property = nullptr;

    auto value_ptr = reinterpret_cast<void *>(prop->data + prop->value.offset);
    if (prop->basic_type == Ioss::Property::BasicType::STRING) {
      std::string property_value(prop->data + prop->value.offset, prop->value.size);
      property = new Ioss::Property(property_name, property_value);
    }
    else if (prop->basic_type == Ioss::Property::BasicType::INTEGER) {
      property = new Ioss::Property(property_name, *(reinterpret_cast<int64_t *>(value_ptr)));
    }
    else if (prop->basic_type == Ioss::Property::BasicType::REAL) {
      property = new Ioss::Property(property_name, *(reinterpret_cast<double *>(value_ptr)));
    }

    return *property;
  }

  void DatabaseIO::read_entity_fields(std::vector<std::string> &keys, Ioss::GroupingEntity &entity)
  {
    // Fields
    for (size_t i = 0; i < keys.size(); i++) {
      std::vector<unsigned char> value;
      bool success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      auto field = reinterpret_cast<field_entry_t *>(value.data());

      std::string field_name(field->data + field->name.offset, field->name.size);
      std::string field_storage(field->data + field->storage.offset, field->storage.size);

      if (!entity.field_exists(field_name)) {
        entity.field_add(Ioss::Field(field_name, field->basic_type, field_storage, field->role_type,
                                     field->raw_count));
      }
    }
  }

  Ioss::Map &DatabaseIO::get_node_map() const
  {
    if (this->nodeMap.defined()) {
      return this->nodeMap;
    }

    std::string nbone_key_str = "::State::-1::Entity::NodeBlock::Name::nodeblock_1::Field::"
                                "RoleType::MESH::BasicType::INTEGER::Name::ids";

    std::vector<unsigned char> value;
    bool success = Ios3::helpers::getValue(helper_context, bucket_name, nbone_key_str, value);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto field = reinterpret_cast<field_entry_t *>(value.data());

    std::string field_name(field->data + field->name.offset, field->name.size);
    std::string field_storage(field->data + field->storage.offset, field->storage.size);

    this->nodeMap.set_size(field->raw_count);
    this->nodeMap.set_map(reinterpret_cast<int *>(field->data + field->value.offset),
                          field->raw_count, 0);

    this->nodeMap.set_defined(true);
    return this->nodeMap;
  }

  void DatabaseIO::get_edgeblocks()
  {
    bool success = true;

    std::string type_string("EdgeBlock");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_edge_type(false);
      std::string original_edge_type;

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }
          if (keys[i].find("original_edge_type") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            original_edge_type      = property_get_int(value);
            have_original_edge_type = true;
          }
        }
      }

      if (have_entity_count && have_original_edge_type) {
        auto block = new Ioss::EdgeBlock(this, entity_name, original_edge_type, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *block);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_elemblocks()
  {
    bool success = true;

    std::string type_string("ElementBlock");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_topology_type(false);
      bool        have_topology_type(false);
      std::string original_topology_type;
      std::string topology_type;

      for (size_t i = 0; i < keys.size(); i++) {

        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("::entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }

          if (keys[i].find("::original_topology_type") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            original_topology_type      = property_get_string(value);
            have_original_topology_type = true;
          }

          if (keys[i].find("::topology_type") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            topology_type      = property_get_string(value);
            have_topology_type = true;
          }
        }
      }

      if (!have_original_topology_type && !have_topology_type) {
        std::string errmsg;
        fmt::print(errmsg, "Unable to reconstruct ElementBlock ({}): couldn't determine topology\n",
                   entity_name.c_str());
        IOSS_ERROR(errmsg);
      }
      if (have_entity_count) {
        if (!have_original_topology_type) {
          original_topology_type = topology_type;
        }
        Ioss::ElementBlock *block =
            new Ioss::ElementBlock(this, entity_name, original_topology_type, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *block);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *block);

        // Add TRANSIENT fields that aren't created in the CTor
        key_t field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
        std::vector<std::string> field_keys_debug;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, field_search_debug.second,
                                          field_keys_debug);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys_debug, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_faceblocks()
  {
    bool success = true;

    std::string type_string("FaceBlock");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_original_topology_type(false);
      std::string original_topology_type;

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }

          if (keys[i].find("original_topology_type") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            original_topology_type      = property_get_string(value);
            have_original_topology_type = true;
          }
        }
      }

      if (have_entity_count && have_original_topology_type) {
        auto block = new Ioss::FaceBlock(this, entity_name, original_topology_type, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *block);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_nodeblocks()
  {
    bool success = true;

    std::string type_string("NodeBlock");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);
      int64_t component_degree(0);

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }

          if (keys[i].find("component_degree") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            component_degree      = property_get_int(value);
            have_component_degree = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        // Creates a new NodeBlock with default Properties and Fields
        auto block = new Ioss::NodeBlock(this, entity_name, entity_count, component_degree);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *block);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *block);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *block);

        // Add TRANSIENT fields that aren't created in the CTor
        key_t field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
        std::vector<std::string> field_keys_debug;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, field_search_debug.second,
                                          field_keys_debug);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys_debug, *block);

        this->get_region()->add(block);
      }
    }
  }

  void DatabaseIO::get_structuredblocks()
  {
    bool success = true;

    std::string type_string("StructuredBlock");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      std::unordered_map<std::string, int> ctor_properties;
      std::list<std::string>               property_names = {
          "component_degree", "ni",        "nj",       "nk",       "ni_global",
          "nj_global",        "nk_global", "offset_i", "offset_j", "offset_k"};

      for (size_t i = 0; i < keys.size(); i++) {
        for (const auto &property_name : property_names) {
          // NOTE: ordinary find won't work b/c the names of multiple properties share substrings,
          // e.g., "ni" and "ni_global" (a string find on "ni" will match both)
          if (keys[i].find(entity_name) != std::string::npos) {
            size_t position = keys[i].rfind(property_name);
            if (position != std::string::npos &&
                (position + property_name.size()) == keys[i].size()) {
              std::vector<unsigned char> value;
              success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
              if (!success) {
                IOSS_ERROR("S3 Operation Failed");
              }
              ctor_properties[property_name] = property_get_int(value);
            }
          }
        }
      }

      if (ctor_properties.size() != property_names.size()) {
        std::string errmsg;
        fmt::print(errmsg, "Unable to reconstruct StructuredBlock ({}): FOUND {} of {} properties",
                   entity_name.c_str(), ctor_properties.size(), property_names.size());
        IOSS_ERROR(errmsg);
      }

      auto block = new Ioss::StructuredBlock(
          this, entity_name, ctor_properties["component_degree"], ctor_properties["ni"],
          ctor_properties["nj"], ctor_properties["nk"], ctor_properties["offset_i"],
          ctor_properties["offset_j"], ctor_properties["offset_k"], ctor_properties["ni_global"],
          ctor_properties["nj_global"], ctor_properties["nk_global"]);

      // Add attributes that aren't created in the CTor
      key_t attribute_key = make_structuredblock_key(parallel_rank(), *(get_region()), *block);
      std::vector<unsigned char> value;
      success = Ios3::helpers::getValue(helper_context, bucket_name, attribute_key.second, value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      unpack_structuredblock(value, *block);

      // Add Properties that aren't created in the CTor
      key_t property_search = property_search_key(parallel_rank(), *(get_region()), *block);
      std::vector<std::string> property_keys;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                        property_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_properties(property_keys, *block);

      // Add fields that aren't created in the CTor
      key_t field_search = field_search_key(parallel_rank(), *(get_region()), *block);
      std::vector<std::string> field_keys;
      success =
          Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_fields(field_keys, *block);

      // Add TRANSIENT fields that aren't created in the CTor
      key_t field_search_debug = field_search_key(parallel_rank(), 1, *(get_region()), *block);
      std::vector<std::string> field_keys_debug;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, field_search_debug.second,
                                        field_keys_debug);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_fields(field_keys_debug, *block);

      this->get_region()->add(block);
    }
  }

  void DatabaseIO::get_edgesets()
  {
    bool success = true;

    std::string type_string("EdgeSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::EdgeSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *entity);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_elemsets()
  {
    bool success = true;

    std::string type_string("ElementSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::ElementSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *entity);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_facesets()
  {
    bool success = true;

    std::string type_string("FaceSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);
      bool    have_component_degree(false);

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count && have_component_degree) {
        auto entity = new Ioss::FaceSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *entity);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_nodesets()
  {
    bool success = true;

    std::string type_string("NodeSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool    have_entity_count(false);
      int64_t entity_count(0);

      for (size_t i = 0; i < keys.size(); i++) {
        std::string entity_name_search = "::" + entity_name + "::";
        if (keys[i].find(entity_name_search) != std::string::npos) {
          if (keys[i].find("entity_count") != std::string::npos) {
            std::vector<unsigned char> value;
            success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
            if (!success) {
              IOSS_ERROR("S3 Operation Failed");
            }
            entity_count      = property_get_int(value);
            have_entity_count = true;
          }
        }
      }

      if (have_entity_count) {
        auto entity = new Ioss::NodeSet(this, entity_name, entity_count);

        // Add Properties that aren't created in the CTor
        key_t property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> property_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                          property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(property_keys, *entity);

        // Add fields that aren't created in the CTor
        key_t field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
        std::vector<std::string> field_keys;
        success =
            Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys, *entity);

        this->get_region()->add(entity);
      }
    }
  }

  void DatabaseIO::get_sidesets()
  {
    bool success = true;

    std::string type_string("SideSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      auto entity = new Ioss::SideSet(this, entity_name);

      // Add Properties that aren't created in the CTor
      key_t property_search = property_search_key(parallel_rank(), *(get_region()), *entity);
      std::vector<std::string> property_keys;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                        property_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_properties(property_keys, *entity);

      // Add fields that aren't created in the CTor
      key_t field_search = field_search_key(parallel_rank(), *(get_region()), *entity);
      std::vector<std::string> field_keys;
      success =
          Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_fields(field_keys, *entity);

      // FIND the SideBlocks that this SideSet references
      key_t sideblocks_key = sideblocks_search_key(parallel_rank(), *(get_region()), *entity);
      std::vector<std::string> sideblocks_search_keys;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, sideblocks_key.second,
                                        sideblocks_search_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      for (size_t i = 0; i < sideblocks_search_keys.size(); i++) {
        std::vector<unsigned char> value;
        success =
            Ios3::helpers::getValue(helper_context, bucket_name, sideblocks_search_keys[i], value);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        int64_t entity_count = unpack_sideblocks(value);

        auto  sideblock_name = get_entity_name(sideblocks_search_keys[i], "SideBlock");
        key_t property_key   = make_property_key(parallel_rank(), *(get_region()), "SideBlock",
                                                 sideblock_name, "STRING", "topology_type");
        std::vector<unsigned char> property_value;
        success = Ios3::helpers::getValue(helper_context, bucket_name, property_key.second,
                                          property_value);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        Ioss::Property topo_property = this->read_property(property_value);

        key_t parent_property_key =
            make_property_key(parallel_rank(), *(get_region()), "SideBlock", sideblock_name,
                              "STRING", "parent_topology_type");
        std::vector<unsigned char> parent_property_value;
        success = Ios3::helpers::getValue(helper_context, bucket_name, parent_property_key.second,
                                          parent_property_value);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        Ioss::Property parent_topo_property = this->read_property(parent_property_value);

        auto sideblock = new Ioss::SideBlock(this, sideblock_name, topo_property.get_string(),
                                             parent_topo_property.get_string(), entity_count);

        // Add Properties that aren't created in the CTor
        key_t sideblock_property_search =
            property_search_key(parallel_rank(), *(get_region()), *sideblock);
        std::vector<std::string> sideblock_property_keys;
        success = Ios3::helpers::listKeys(
            helper_context, bucket_name, sideblock_property_search.second, sideblock_property_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_properties(sideblock_property_keys, *sideblock);

        // Add fields that aren't created in the CTor
        key_t sideblock_field_search =
            field_search_key(parallel_rank(), *(get_region()), *sideblock);
        std::vector<std::string> sideblock_field_keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name,
                                          sideblock_field_search.second, sideblock_field_keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(sideblock_field_keys, *sideblock);

        // Add TRANSIENT fields that aren't created in the CTor
        key_t field_search_debug =
            field_search_key(parallel_rank(), 1, *(get_region()), *sideblock);
        std::vector<std::string> field_keys_debug;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, field_search_debug.second,
                                          field_keys_debug);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        this->read_entity_fields(field_keys_debug, *sideblock);

        entity->add(sideblock);
      }

      this->get_region()->add(entity);
    }
  }

  void DatabaseIO::get_commsets()
  {
    bool success = true;

    std::string type_string("CommSet");
    key_t       search_key = entity_search_key(parallel_rank(), *(get_region()), type_string);
    std::vector<std::string> keys;
    success = Ios3::helpers::listKeys(helper_context, bucket_name, search_key.second, keys);
    if (!success) {
      IOSS_ERROR("S3 Operation Failed");
    }

    auto entity_names = get_entity_names(keys, type_string);
    for (const auto &entity_name : entity_names) {
      bool        have_entity_count(false);
      int64_t     entity_count(0);
      bool        have_entity_type(false);
      std::string entity_type("");

      for (size_t i = 0; i < keys.size(); i++) {
        if (keys[i].find(entity_name) != std::string::npos) {
          /* Check if we have all of the necessary attributes */
          if (have_entity_type && have_entity_count)
            break;

          if (!have_entity_count) {
            size_t position = keys[i].rfind("entity_count");
            if (position != std::string::npos) {
              std::vector<unsigned char> value;
              success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
              if (!success) {
                IOSS_ERROR("S3 Operation Failed");
              }
              entity_count      = property_get_int(value);
              have_entity_count = true;
            }
          }

          if (!have_entity_type) {
            size_t position = keys[i].rfind("entity_type");
            if (position != std::string::npos) {
              std::vector<unsigned char> value;
              success = Ios3::helpers::getValue(helper_context, bucket_name, keys[i], value);
              if (!success) {
                IOSS_ERROR("S3 Operation Failed");
              }
              entity_type      = property_get_string(value);
              have_entity_type = true;
            }
          }
        }
      }

      if (!have_entity_type) {
        std::string errmsg;
        fmt::print(errmsg, "Unable to reconstruct CommSet (entity_type NOT found)\n");
        IOSS_ERROR(errmsg);
      }

      if (!have_entity_count) {
        std::string errmsg;
        fmt::print(errmsg, "Unable to reconstruct CommSet (entity_count NOT found)\n");
        IOSS_ERROR(errmsg);
      }

      auto commset = new Ioss::CommSet(this, entity_name, entity_type, entity_count);

      // Add Properties that aren't created in the CTor
      key_t property_search = property_search_key(parallel_rank(), *(get_region()), *commset);
      std::vector<std::string> property_keys;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, property_search.second,
                                        property_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_properties(property_keys, *commset);

      // Add fields that aren't created in the CTor
      key_t field_search = field_search_key(parallel_rank(), *(get_region()), *commset);
      std::vector<std::string> field_keys;
      success =
          Ios3::helpers::listKeys(helper_context, bucket_name, field_search.second, field_keys);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_fields(field_keys, *commset);

      // Add TRANSIENT fields that aren't created in the CTor
      key_t field_search_transient =
          field_search_key(parallel_rank(), 1, *(get_region()), *commset);
      std::vector<std::string> field_keys_transient;
      success = Ios3::helpers::listKeys(helper_context, bucket_name, field_search_transient.second,
                                        field_keys_transient);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
      this->read_entity_fields(field_keys_transient, *commset);

      this->get_region()->add(commset);
    }
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*reg, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*nb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*eb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*fb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    bool   success    = true;
    size_t num_to_get = field.verify(data_size);

    int64_t count = get_field_internal(*eb, field, data, data_size);
    if (count > 0) {
      if (field.get_name() == "connectivity_raw") {
        const auto &node_map = get_node_map();
        int        *dup_data = (int *)malloc(data_size);
        memcpy(dup_data, data, data_size);
        node_map.map_data(dup_data, field,
                          field.verify(data_size) * field.raw_storage()->component_count());
        free(dup_data);
      }
      else if (field.get_name() == "connectivity") {
        const auto &node_map = get_node_map();
        int        *dup_data = (int *)malloc(data_size);
        memcpy(dup_data, data, data_size);
        node_map.reverse_map_data(dup_data, field,
                                  field.verify(data_size) * field.raw_storage()->component_count());
        free(dup_data);
      }
    }
    if (count == 0) {
      if (field.get_name() == "connectivity_raw") {
        Ios3::key_t key = make_key(parallel_rank(), *(get_region()), *eb, field, "connectivity");
        std::vector<std::string> keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, key.second, keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        if (keys.size() == 1) {
          std::vector<unsigned char> value;
          success = Ios3::helpers::getValue(helper_context, bucket_name, key.second, value);
          if (!success) {
            IOSS_ERROR("S3 Operation Failed");
          }

          field_entry_t *field_ptr(
              reinterpret_cast<field_entry_t *>(reinterpret_cast<void *>(value.data())));

          // TODO what other checks do we need here?
          if (data_size != field_ptr->value.size) {
            std::string errmsg;
            fmt::print(errmsg, "ERROR: data_size({}) != field_ptr->value.size({})", data_size,
                       field_ptr->value.size);
            IOSS_ERROR(errmsg);
          }

          std::memcpy(data, reinterpret_cast<void *>(field_ptr->data + field_ptr->value.offset),
                      field_ptr->value.size);

          get_node_map().reverse_map_data(
              data, field, field.verify(data_size) * field.raw_storage()->component_count());
        }
      }
      else if (field.get_name() == "connectivity") {
        Ios3::key_t key =
            make_key(parallel_rank(), *(get_region()), *eb, field, "connectivity_raw");
        std::vector<std::string> keys;
        success = Ios3::helpers::listKeys(helper_context, bucket_name, key.second, keys);
        if (!success) {
          IOSS_ERROR("S3 Operation Failed");
        }
        if (keys.size() == 1) {
          std::vector<unsigned char> value;
          success = Ios3::helpers::getValue(helper_context, bucket_name, key.second, value);
          if (!success) {
            IOSS_ERROR("S3 Operation Failed");
          }

          field_entry_t *field_ptr(
              reinterpret_cast<field_entry_t *>(reinterpret_cast<void *>(value.data())));

          // TODO what other checks do we need here?
          if (data_size != field_ptr->value.size) {
            std::string errmsg;
            fmt::print(errmsg, "ERROR: data_size({}) != field_ptr->value.size({})", data_size,
                       field_ptr->value.size);
            IOSS_ERROR(errmsg);
          }

          std::memcpy(data, reinterpret_cast<void *>(field_ptr->data + field_ptr->value.offset),
                      field_ptr->value.size);

          get_node_map().map_data(data, field,
                                  field.verify(data_size) * field.raw_storage()->component_count());
        }
      }
    }

    return num_to_get;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*sb, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ns, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*es, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*fs, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*es, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return get_field_internal(*ss, field, data, data_size);
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);

    if (num_to_get > 0) {
      if (field.get_name() == "entity_processor" || field.get_name() == "entity_processor_raw") {
        get_field_internal(*cs, field, data, data_size);
      }
      else if (field.get_name() == "ids") {
        // Do nothing, just handles an idiosyncrasy of the GroupingEntity
      }
      else {
        num_to_get = Ioss::Utils::field_warning(cs, field, "input");
      }
    }
    return num_to_get;
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
      for (uint64_t id(0); id < field.raw_count(); ++id) {
        data_ptr[index++] = data_x[id];
        if (dim > 1)
          data_ptr[index++] = data_y[id];
        if (dim > 2)
          data_ptr[index++] = data_z[id];
      }
    }
    else {
      return get_field_internal(*sb, field, data, data_size);
    }
    return 0;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*reg, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*nb, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*eb, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*fb, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*eb, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*sb, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ns, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*es, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*fs, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*es, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*ss, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*cs, field, data, data_size);
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    return put_field_internal(*sb, field, data, data_size);
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = f.verify(data_size);

    if (num_to_get > 0) {
      key_t key = make_key(parallel_rank(), *(get_region()), e, f);

      std::vector<unsigned char> value;
      bool success = Ios3::helpers::getValue(helper_context, bucket_name, key.second, value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }

      field_entry_t *field_ptr(
          reinterpret_cast<field_entry_t *>(reinterpret_cast<void *>(value.data())));

      // TODO what other checks do we need here?
      if (data_size != field_ptr->value.size) {
        std::string errmsg;
        fmt::print(errmsg, "ERROR: data_size({}) != field_ptr->value.size({})", data_size,
                   field_ptr->value.size);
        IOSS_ERROR(errmsg);
      }

      std::memcpy(data, reinterpret_cast<void *>(field_ptr->data + field_ptr->value.offset),
                  field_ptr->value.size);
    }

    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::GroupingEntity &e, const Ioss::Field &f,
                                         void *data, size_t data_size) const
  {
    size_t num_to_put = f.verify(data_size);

    if (num_to_put > 0) {
      key_t key     = make_key(parallel_rank(), *(get_region()), e, f);
      auto  value   = pack_field(*(get_region()), e, f, data, data_size);
      bool  success = Ios3::helpers::putValue(helper_context, bucket_name, key.second, value);
      if (!success) {
        IOSS_ERROR("S3 Operation Failed");
      }
    }
    return num_to_put;
  }

} // namespace Ios3
