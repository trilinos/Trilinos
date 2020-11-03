// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include "Ioss_CommSet.h"      // for CommSet
#include "Ioss_EdgeBlock.h"    // for EdgeBlock
#include "Ioss_EdgeSet.h"      // for EdgeSet
#include "Ioss_ElementBlock.h" // for ElementBlock
#include "Ioss_ElementSet.h"   // for ElementSet
#include "Ioss_FaceBlock.h"    // for FaceBlock
#include "Ioss_FaceSet.h"      // for FaceSet
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h" // for NodeSet
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_VariableType.h"
#include <Ioss_SubSystem.h>

#include "Ioss_DatabaseIO.h" // for DatabaseIO

#include "Ioss_IOFactory.h" // for IOFactory
#include <init/Ionit_Initializer.h>

#include "adios/Ioad_Constants.h"
#include "adios/Ioad_Helper.h"
#include "adios/Ioad_TemplateToValue.h"

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#else
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 1
#endif
#endif

#include <algorithm>
#include <iostream>
#include <ostream>
#include <stddef.h> // for size_t
#include <stdlib.h> // for rand, srand, RAND_MAX
#include <string>
#include <vector>

//////// Global constants ////////////
std::vector<std::string> ignored_properties = {"database_name"};
std::vector<std::string> ignored_fields     = {"implicit_ids"};
std::vector<std::string> ignore_errors      = {"owning_processor"};

////////////////////////////

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  const int result = Catch::Session().run(argc, argv);
#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif

  return result;
}

#define COMPARE_POINTERS_AND_THROW_WITH_INDEX(obj_type, obj1, obj2, meth, index)                   \
  if (obj1->meth(index) != obj2->meth(index)) {                                                    \
    std::ostringstream error;                                                                      \
    error << #obj_type " " #meth " do not match:\n"                                                \
          << "1:" #index " - " << obj1->meth(index) << "\n"                                        \
          << "2:" #index " - " << obj2->meth(index) << "\n";                                       \
    throw(error.str());                                                                            \
  }

void CompareVariableType(const Ioss::VariableType *var1, const Ioss::VariableType *var2)
{
  COMPARE_POINTERS_AND_THROW_WITH_INDEX(VariableType, var1, var2, name, );
  COMPARE_POINTERS_AND_THROW_WITH_INDEX(VariableType, var1, var2, component_count, );
  COMPARE_POINTERS_AND_THROW_WITH_INDEX(VariableType, var1, var2, suffix_count, );

  for (int i = 1; i <= var1->component_count(); i++) {
    COMPARE_POINTERS_AND_THROW_WITH_INDEX(VariableType, var1, var2, label, i);
    // COMPARE_POINTERS_AND_THROW_WITH_INDEX(VariableType, var1, var2, label_name, i);
  }
}

#define COMPARE_VALUES_AND_THROW(obj_type, obj1, obj2, meth)                                       \
  if (obj1.meth() != obj2.meth()) {                                                                \
    std::ostringstream error;                                                                      \
    error << #obj_type " " #meth " do not match:\n"                                                \
          << "1: " << obj1.meth() << "\n"                                                          \
          << "2: " << obj2.meth() << "\n";                                                         \
    throw(error.str());                                                                            \
  }

template <typename T> bool CompareVectors(const std::vector<T> &v1, const std::vector<T> &v2)
{
  std::vector<T> v3;
  std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v3));
  if (v3.empty()) {
    return true;
  }
  else {
    return false;
  }
}

#define COMPARE_VECTORS_AND_THROW(name, vec1, vec2)                                                \
  if (!CompareVectors(vec1, vec2)) {                                                               \
    std::ostringstream error;                                                                      \
    error << "Vectors " << name << " do not match:\n";                                             \
    throw(error.str());                                                                            \
  }

void CompareFields(const Ioss::Field &field1, const Ioss::Field &field2)
{
  // Check that type is the same for both fields
  COMPARE_VALUES_AND_THROW(Field, field1, field2, get_type);
  CompareVariableType(field1.raw_storage(), field2.raw_storage());
  CompareVariableType(field1.transformed_storage(), field2.transformed_storage());
  COMPARE_VALUES_AND_THROW(Field, field1, field2, raw_count);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, transformed_count);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, get_size);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, get_role);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, has_transform);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, is_valid);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, is_invalid);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, get_name);
  COMPARE_VALUES_AND_THROW(Field, field1, field2, get_type);
}

// Compare fields
// compare properties
void CompareProperties(const Ioss::Property &prop1, const Ioss::Property &prop2)
{
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_name);
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_type);
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, is_valid);
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, is_invalid);
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, is_explicit);
  COMPARE_VALUES_AND_THROW(Property, prop1, prop2, is_implicit);
  switch (prop1.get_type()) {
  case Ioss::Property::BasicType::REAL:
    COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_real);
    break;
  case Ioss::Property::BasicType::INTEGER:
    COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_int);
    break;
  case Ioss::Property::BasicType::POINTER:
    // TODO: Verify that is makes sense.
    COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_pointer);
    break;
  case Ioss::Property::BasicType::STRING:
    COMPARE_VALUES_AND_THROW(Property, prop1, prop2, get_string);
    break;
  case Ioss::Property::BasicType::INVALID:
  default: throw("Unsupported property type " + prop1.get_name());
  }
}

template <typename T, typename TF>
void CompareFieldData(const T entity_block1, const T entity_block2, const std::string &field_name)
{
  std::vector<TF> data1;
  std::vector<TF> data2;
  entity_block1->get_field_data(field_name, data1);
  entity_block2->get_field_data(field_name, data2);
  COMPARE_VECTORS_AND_THROW(field_name, data1, data2);
}

template <typename T>
void CompareFieldData(const T entity_block1, const T entity_block2, const std::string &field_name)
{
  auto field = entity_block1->get_field(field_name);
  switch (field.get_type()) {
  case Ioss::Field::BasicType::DOUBLE:
    CompareFieldData<T, double>(entity_block1, entity_block2, field_name);
    break;
  case Ioss::Field::BasicType::INT32:
    CompareFieldData<T, int32_t>(entity_block1, entity_block2, field_name);
    break;
  case Ioss::Field::BasicType::INT64:
    CompareFieldData<T, int64_t>(entity_block1, entity_block2, field_name);
    break;
  case Ioss::Field::BasicType::COMPLEX:
    // CompareFieldData<T, Complex>(entity_block1, entity_block2, field_name);
    break;
  case Ioss::Field::BasicType::CHARACTER:
    CompareFieldData<T, char>(entity_block1, entity_block2, field_name);
    break;
  default:
    std::ostringstream errmsg;
    throw("INTERNAL ERROR: Invalid field type. Something is wrong in the the input entity block.");
  }
}

template <typename T> void CompareAllProperties(const T &obj1, const T &obj2)
{
  std::vector<std::string> block1_property_list;
  obj1->property_describe(&block1_property_list);
  std::vector<std::string> block2_property_list;
  obj2->property_describe(&block2_property_list);
  for (auto property_name1 : block1_property_list) {
    if (std::find(ignored_properties.begin(), ignored_properties.end(), property_name1) !=
        std::end(ignored_properties)) {
      continue;
    }
    obj1->get_property(property_name1);
    auto property_name2 =
        std::find(std::begin(block2_property_list), std::end(block2_property_list), property_name1);
    if (property_name2 != block2_property_list.end()) {
      // Check that the content is the same.
      auto property1 = obj1->get_property(property_name1);
      auto property2 = obj2->get_property(property_name1);

      CompareProperties(property1, property2);
    }
    else {
      throw("Property name " + property_name1 + " not found in db2");
    }
  }
}

template <typename T>
void CompareContainers(const T &entity_blocks1, const T &entity_blocks2, Ioss::State state)
{
  for (auto entity_block1 : entity_blocks1) {
    // Find corresponding block in vector

    auto entity_block2 =
        std::find_if(entity_blocks2.begin(), entity_blocks2.end(),
                     [=](typename T::value_type e) { return e->name() == entity_block1->name(); });
    if (entity_block2 != entity_blocks2.end()) {
      Ioss::NameList block1_field_names;
      entity_block1->field_describe(&block1_field_names);
      Ioss::NameList block2_field_names;
      (*entity_block2)->field_describe(&block2_field_names);
      // Sorts vectors before comparing them
      std::sort(block1_field_names.begin(), block1_field_names.end());
      std::sort(block2_field_names.begin(), block2_field_names.end());
      COMPARE_VECTORS_AND_THROW(entity_block1->name(), block1_field_names, block2_field_names);

      for (auto name1 : block1_field_names) {
        if (std::find(ignored_fields.begin(), ignored_fields.end(), name1) !=
            std::end(ignored_fields)) {
          continue;
        }
        // Check that the content is the same.
        auto field1 = entity_block1->get_fieldref(name1);
        auto field2 = (*entity_block2)->get_fieldref(name1);
        CompareFields(field1, field2);
        // Negated logical XOR
        if ((state == Ioss::STATE_TRANSIENT) ==
            (field1.get_role() == Ioss::Field::RoleType::TRANSIENT)) {
          CompareFieldData(entity_block1, (*entity_block2), name1);
        }
      }

      CompareAllProperties(entity_block1, (*entity_block2));
    }
  }
}

void CompareRegions(Ioss::Region *region1, Ioss::Region *region2)
{

  CompareAllProperties(region1, region2);
  Ioss::State state = region1->get_state();
  // Compare all containers
  CompareContainers(region1->get_node_blocks(), region2->get_node_blocks(), state);
  CompareContainers(region1->get_edge_blocks(), region2->get_edge_blocks(), state);
  CompareContainers(region1->get_face_blocks(), region2->get_face_blocks(), state);
  CompareContainers(region1->get_element_blocks(), region2->get_element_blocks(), state);
  CompareContainers(region1->get_nodesets(), region2->get_nodesets(), state);
  CompareContainers(region1->get_edgesets(), region2->get_edgesets(), state);
  CompareContainers(region1->get_facesets(), region2->get_facesets(), state);
  CompareContainers(region1->get_elementsets(), region2->get_elementsets(), state);
  CompareContainers(region1->get_commsets(), region2->get_commsets(), state);
  CompareContainers(region1->get_sidesets(), region2->get_sidesets(), state);
  // CompareContainers(region1->get_sideblocks(), region2->get_sideblocks());

  // Is there any global field or map?
}

void CompareDB(Ioss::DatabaseIO *db1, Ioss::DatabaseIO *db2)
{
  std::shared_ptr<Ioss::Region> region1(new Ioss::Region(db1));
  std::shared_ptr<Ioss::Region> region2(new Ioss::Region(db2));
  CompareRegions(region1.get(), region2.get());
  // Region properties have been compared previously.
  int timestep_count = region1->get_property("state_count").get_int();
  for (int step = 1; step <= timestep_count; step++) {
    region1->begin_state(step);
    region2->begin_state(step);

    CompareRegions(region1.get(), region2.get());
    region1->end_state(step);
    region2->end_state(step);
  }
}

template <typename Entity, typename T>
void put_field_data(std::string field_name, int local_size, size_t component_count, Entity *e)
{
  std::vector<T> data;
  size_t         data_size = local_size * component_count;
  data.reserve(data_size);
  for (size_t i = 0; i < data_size; ++i) {
    if (field_name == "owning_processor") {
      int rank = 0;
#ifdef SEACAS_HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      data.push_back(rank);
    }
    else {
      data.push_back(i + 1); // Ids must be >0
    }
  }

  int num_ids_written = e->put_field_data(field_name, data);
  if (std::find(ignore_errors.begin(), ignore_errors.end(), field_name) !=
      std::end(ignore_errors)) {
    return;
  }
  if (local_size != num_ids_written) {
    std::ostringstream msg;
    msg << " FAILED in put_field_data:";
    msg << " field_name = " << field_name;
    msg << " , num_nodes = " << local_size;
    msg << " , num_ids_written = " << num_ids_written;
    throw std::runtime_error(msg.str());
  }
}

template <typename Entity> void write_fields(Entity *e, Ioss::Field::RoleType role)
{
  std::vector<std::string> field_names;
  e->field_describe(&field_names);
  for (auto field_name : field_names) {
    Ioss::Field field = e->get_field(field_name);
    if (field.get_role() != role) {
      continue;
    }
    std::string entity_type = e->type_string();

    if (Ioad::find_field_in_mapset(entity_type, field_name, Ioad::Ignore_fields)) {
      continue;
    }
    size_t component_count;
    size_t local_size;

    if (Ioad::use_transformed_storage(field, entity_type, field_name)) {
      component_count = field.transformed_storage()->component_count();
      local_size      = field.transformed_count();
    }
    else {
      component_count = field.raw_storage()->component_count();
      local_size      = field.raw_count();
    }
    switch (field.get_type()) {
    case Ioss::Field::BasicType::DOUBLE:
      put_field_data<Ioss::NodeBlock, double>(field_name, local_size, component_count, e);
      break;
    case Ioss::Field::BasicType::INT32:
      put_field_data<Ioss::NodeBlock, int32_t>(field_name, local_size, component_count, e);
      break;
    case Ioss::Field::BasicType::INT64:
      put_field_data<Ioss::NodeBlock, int64_t>(field_name, local_size, component_count, e);
      break;
    case Ioss::Field::BasicType::COMPLEX:
      put_field_data<Ioss::NodeBlock, Complex>(field_name, local_size, component_count, e);
      break;
    case Ioss::Field::BasicType::CHARACTER:
      put_field_data<Ioss::NodeBlock, char>(field_name, local_size, component_count, e);
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Invalid field type. "
             << "Something is wrong in the Ioad::DatabaseIO::get_field_internal_t() function. "
             << "Please report.\n";
      IOSS_ERROR(errmsg);
    }
  }
}

void create_phantom(Ioss::DatabaseIO *db)
{
  // Create region.
  Ioss::Region *region = new Ioss::Region(db);
  region->begin_mode(Ioss::STATE_DEFINE_MODEL);
  db->set_region(region);
  // Create a NodeBlock with some fields
  int64_t          node_count = 10;
  Ioss::NodeBlock *node_block = new Ioss::NodeBlock(db, "nodeblock_1", node_count, 3);
  region->add(node_block);
  // Define some transient fields
  std::vector<std::string> field_names{"field1", "field2"};
  for (auto field_name : field_names) {
    // Exodus seems to always save TRANSIENT fields as `REAL`, so for the tests to pass
    // we only use this type.
    Ioss::Field field(field_name, Ioss::Field::BasicType::REAL, "scalar",
                      Ioss::Field::RoleType::TRANSIENT);
    node_block->field_add(field);
  }
  // Add coordinate frames
  std::vector<double>   coords(9, 0);
  Ioss::CoordinateFrame cf1(0, 'a', coords.data());
  region->add(cf1);
  Ioss::CoordinateFrame cf2(1, 'b', coords.data());
  region->add(cf2);
  // Fill up the fields with some values.
  region->end_mode(Ioss::STATE_DEFINE_MODEL);
  region->begin_mode(Ioss::STATE_MODEL);
  write_fields(node_block, Ioss::Field::RoleType::MESH);
  region->end_mode(Ioss::STATE_MODEL);
  region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
  region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  region->begin_mode(Ioss::STATE_TRANSIENT);
  region->add_state(0.1);
  region->begin_state(1);
  write_fields(node_block, Ioss::Field::RoleType::TRANSIENT);
  region->end_state(1);
  region->end_mode(Ioss::STATE_TRANSIENT);
}

void create_database(std::string type, std::string file_name)
{
  std::shared_ptr<Ioss::DatabaseIO> db(
      Ioss::IOFactory::create(type, file_name, Ioss::WRITE_RESULTS, MPI_COMM_WORLD));
  create_phantom(db.get());
  db->closeDatabase();
}

TEST_CASE("Ioad", "[Ioad]")
{

  Ioss::Init::Initializer::initialize_ioss();
  std::string exodus_db_name = "phantom.e";
  std::string adios_db_name  = "phantom.bp";
  create_database("exodus", exodus_db_name);
  create_database("adios", adios_db_name);
  // Database pointers are deleted by their respective region destructors.
  Ioss::DatabaseIO *read_exodus_db =
      Ioss::IOFactory::create("exodus", exodus_db_name, Ioss::READ_MODEL, MPI_COMM_WORLD);
  read_exodus_db->set_int_byte_size_api(Ioss::USE_INT64_API);
  Ioss::DatabaseIO *read_adios_db =
      Ioss::IOFactory::create("adios", adios_db_name, Ioss::READ_MODEL, MPI_COMM_WORLD);
  CompareDB(read_exodus_db, read_adios_db);
}

template <typename T> const std::string get_entity_type_test()
{
  // Use "node" as default entity type to enable factory to create object.
  std::unique_ptr<T> entity(Ioad::NewEntity<T>(nullptr, "", "node", 0));
  return entity->type_string();
}

// NodeBlock has a different constructor...
template <> const std::string get_entity_type_test<Ioss::NodeBlock>()
{
  Ioss::NodeBlock nodeblock(nullptr, "", 0, 1);
  return nodeblock.type_string();
}

// SideSet has a different constructor...
template <> const std::string get_entity_type_test<Ioss::SideSet>()
{
  Ioss::SideSet sideset(nullptr, "");
  return sideset.type_string();
}

template <> const std::string get_entity_type_test<Ioss::SideBlock>()
{
  // default element arbitrarily set to "hex8": We need an element type
  // to be able to construct this sideblock.
  Ioss::SideBlock sideblock(nullptr, "", "node", "hex8", 0);
  return sideblock.type_string();
}

TEST_CASE("Ioad_BlockNames", "[Ioad]")
{
  REQUIRE(get_entity_type_test<Ioss::SideBlock>() == Ioad::get_entity_type<Ioss::SideBlock>());
  REQUIRE(get_entity_type_test<Ioss::SideSet>() == Ioad::get_entity_type<Ioss::SideSet>());
  REQUIRE(get_entity_type_test<Ioss::NodeBlock>() == Ioad::get_entity_type<Ioss::NodeBlock>());
  REQUIRE(get_entity_type_test<Ioss::EdgeBlock>() == Ioad::get_entity_type<Ioss::EdgeBlock>());
  REQUIRE(get_entity_type_test<Ioss::FaceBlock>() == Ioad::get_entity_type<Ioss::FaceBlock>());
  REQUIRE(get_entity_type_test<Ioss::ElementBlock>() ==
          Ioad::get_entity_type<Ioss::ElementBlock>());
  REQUIRE(get_entity_type_test<Ioss::NodeSet>() == Ioad::get_entity_type<Ioss::NodeSet>());
  REQUIRE(get_entity_type_test<Ioss::EdgeSet>() == Ioad::get_entity_type<Ioss::EdgeSet>());
  REQUIRE(get_entity_type_test<Ioss::FaceSet>() == Ioad::get_entity_type<Ioss::FaceSet>());
  REQUIRE(get_entity_type_test<Ioss::ElementSet>() == Ioad::get_entity_type<Ioss::ElementSet>());
  REQUIRE(get_entity_type_test<Ioss::CommSet>() == Ioad::get_entity_type<Ioss::CommSet>());
}
