// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioad_DatabaseIO_h
#define IOSS_Ioad_DatabaseIO_h

#include "Ioss_EntitySet.h"
#include "Ioss_Region.h"  // for Region, SideSetContainer, etc
#include "Ioss_SideSet.h" // for SideBlockContainer, SideSet
#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>

#include "Ioss_Field.h" // for Field, etc
#include <AdiosWrapper.h>

namespace Ioss {
  class Assembly;
  class Blob;
  class GroupingEntity;
  class Region;
  class EntityBlock;
  class NodeBlock;
  class SideBlock;
  class ElementBlock;
  class NodeSet;
  class SideSet;
  class CommSet;
} // namespace Ioss

/** \brief A namespace for the adios database format.
 */
namespace Ioad {

  class DatabaseIO : public Ioss::DatabaseIO
  {
  public:
    DatabaseIO(Ioss::Region *region, const std::string &filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &props);
    ~DatabaseIO();
    DatabaseIO(const DatabaseIO &from) = delete;
    DatabaseIO &operator=(const DatabaseIO &from) = delete;

    const std::string get_format() const override { return "ADIOS2"; }

    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    unsigned entity_field_support() const override;
    int      int_byte_size_db() const override;

  private:
    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
    {
      return -1;
    }
    int64_t get_field_internal(const Ioss::Assembly *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return -1;
    }
    int64_t get_field_internal(const Ioss::Blob *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return -1;
    }
    int64_t get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t get_field_internal_t(const Ioss::GroupingEntity *entity, const Ioss::Field &field,
                                 void *data, size_t data_size) const;
    template <typename T>
    void get_data(void *data, const std::string &encoded_name,
                  bool use_step_selection = false) const;

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override
    {
      return -1;
    }
    int64_t put_field_internal(const Ioss::Assembly *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return -1;
    }
    int64_t put_field_internal(const Ioss::Blob *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override
    {
      return -1;
    }
    template <typename T>
    int64_t put_field_internal_t(T entity, const Ioss::Field &field, void *data,
                                 size_t data_size) const;
    void    define_model(Ioss::Field::RoleType *role = nullptr);
    // Model definition that should not be re-defined when defining transient variables.
    void                    define_global_variables();
    template <typename T> T get_attribute(const std::string &attribute_name);

    template <typename T> void put_data(void *data, const std::string &encoded_name) const;
    template <typename T, typename = typename std::enable_if<
                              !std::is_base_of<Ioss::EntitySet, T>::value, T>::type>
    void put_var_type(const Ioss::Field &field, const std::string &encoded_name,
                      bool transformed_field) const;
    template <typename T>
    void define_model_internal(const Ioss::Field &field, const std::string &encoded_name,
                               const std::string &entity_type, const std::string &field_name);
    template <typename T>
    void define_entity_internal(const T &entity_blocks, Ioss::Field::RoleType *role);

    int get_current_state() const;

    struct BlockInfoType
    {
      std::vector<size_t> steps;
      adios2::Dims        Count;
      size_t              global_size{0};
    };

    struct FieldInfoType
    {
      // Information contained in block infos.
      std::vector<size_t> steps;
      size_t              node_boundaries_size = 0;
      size_t              component_count      = 0;
      // Contained in metavariables
      Ioss::Field::RoleType  role;
      std::string            variable_type;
      Ioss::Field::BasicType basic_type;
      std::string            topology;
      std::string            parent_topology;
    };

    template <typename T> BlockInfoType get_block_infos(const adios2::Variable<T> &var) const;

    template <typename T> FieldInfoType get_variable_infos(const std::string &var_name) const;
    using GlobalMapType = std::map<std::string, std::pair<std::string, std::string>>;
    using EntityMapType = std::map<std::string, GlobalMapType>;
    using FieldsMapType = std::map<std::string, EntityMapType>;

    template <typename T>
    std::string get_property_value(const FieldsMapType &properties_map,
                                   const std::string &entity_type, const std::string &entity_name,
                                   const std::string &property_name) const;

    template <typename T>
    FieldInfoType get_expected_variable_infos_from_map(const EntityMapType &fields_map,
                                                       const std::string &  entity_type,
                                                       const std::string &  entity_name,
                                                       const std::string &  var_name) const;
    FieldInfoType get_variable_infos_from_map(const EntityMapType &fields_map,
                                              const std::string &  entity_type,
                                              const std::string &  entity_name,
                                              const std::string &  var_name) const;

    template <typename T>
    using IsIossEntityBlock =
        typename std::enable_if<std::is_base_of<Ioss::EntityBlock, T>::value>::type;
    template <typename T>
    using IsNotIossEntityBlock =
        typename std::enable_if<!std::is_base_of<Ioss::EntityBlock, T>::value>::type;

    template <typename T, typename = IsIossEntityBlock<T>>
    void define_entity_meta_variables(const std::string &encoded_name);

    template <typename T, typename = IsNotIossEntityBlock<T>, typename = void>
    void define_entity_meta_variables(const std::string &encoded_name);

    void define_field_meta_variables(const std::string &);
    void define_coordinate_frames_internal(const Ioss::CoordinateFrameContainer &coordinate_frames);
    std::string encoded_coordinate_frame_name(Ioss::CoordinateFrame coordinate_frame);

    void put_meta_variables(const std::string &encoded_name, const Ioss::Field &field,
                            const std::string &entity_type, const std::string &field_name) const;
    void write_meta_data();

    template <typename T>
    void add_entity_property(Ioss::GroupingEntity *ge, const std::string &encoded_name,
                             const std::string &var_name);
    void add_entity_properties(Ioss::GroupingEntity *ge, const FieldsMapType &properties_map,
                               std::string name = "");

    void write_properties(const Ioss::GroupingEntity *const entity,
                          const std::string &               encoded_name);

    template <typename T> int64_t write_meta_data_container(const T &entity_blocks);
    std::pair<int64_t, int64_t>
    write_meta_data_sideblockcontainer(const Ioss::SideBlockContainer &entity_blocks);

    template <typename T>
    int64_t     get_entities(const FieldsMapType &fields_map, const FieldsMapType &properties_map);
    std::string get_optional_string_variable(const std::string &field_name,
                                             const std::string &string_variable) const;

    void get_globals(const GlobalMapType &globals_map, const FieldsMapType &properties_map);
    void compute_block_membership__(Ioss::SideBlock *         efblock,
                                    std::vector<std::string> &block_membership) const override;
    void define_properties(const Ioss::GroupingEntity *entity_block,
                           const std::string &         encoded_entity_name);

    void read_meta_data__() override;
    void read_communication_metadata();
    void read_region(const FieldsMapType &fields_map);
    void check_processor_info();
    void check_model();

    int           RankInit();
    bool          begin_state__(int state, double time) override;
    bool          end_state__(int state, double time) override;
    unsigned long rank; // rank needs to be declared first to be initialized before adios_wrapper.
    mutable AdiosWrapper adios_wrapper; // adios_wrapper needs to be declared before bpio
                                        // and bp_engine to be initialized first.
    int           spatialDimension{0};
    int64_t       edgeCount{0};
    int64_t       faceCount{0};
    unsigned long number_proc;
    bool          is_streaming;
    double        previous_time_streaming;
  };
} // namespace Ioad
#endif
