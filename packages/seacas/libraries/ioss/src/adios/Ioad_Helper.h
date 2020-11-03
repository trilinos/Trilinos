// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioad_Helper_h
#define IOSS_Ioad_Helper_h

#include <map>
#include <set>
#include <string>

#include "Ioss_GroupingEntity.h" // for GroupingEntity
#include "Ioss_SideSet.h"        // for SideBlockContainer, SideSet

namespace Ioad {

  template <typename T>
  using DerivedFromIossGroupingEntity =
      typename std::enable_if<std::is_base_of<Ioss::GroupingEntity, T>::value, bool>::type;

  template <typename T>
  using IossHas3ParametersConstructor =
      decltype(DerivedFromIossGroupingEntity<T>{}, T(nullptr, std::string{}, int64_t{}));

  template <typename T>
  using IossHas4ParametersConstructor = decltype(
      DerivedFromIossGroupingEntity<T>{}, T(nullptr, std::string{}, std::string{}, int64_t{}));

  // Takes an extra unused parameter "entity_type" to match the API of the function
  // "NewEntity" used for objects that are not EntitySets which require that parameter.
  // Having matching APIs allows to call this function from generic templated functions
  // that do not have to be specialized to call this function with a different number of
  // parameters.
  template <typename T>
  auto NewEntity(Ioss::DatabaseIO *io_database, const std::string &my_name,
                 const std::string & /*entity_type*/, size_t       entity_count)
      -> IossHas3ParametersConstructor<T> *
  {
    return new T(io_database, my_name, entity_count);
  }

  template <typename T>
  auto NewEntity(Ioss::DatabaseIO *io_database, const std::string &my_name,
                 const std::string &entity_type, size_t entity_count)
      -> IossHas4ParametersConstructor<T> *
  {
    return new T(io_database, my_name, entity_type, entity_count);
  }

  int find_field_in_mapset(const std::string &entity_type, const std::string &field_name,
                           const std::map<std::string, std::set<std::string>> &mapset);

  std::string              get_property_variable_name(const std::string &property_name);
  std::vector<std::string> properties_to_save(const Ioss::GroupingEntity *const entity_block);

  std::string stringify_side_block_names(const Ioss::SideBlockContainer &sblocks);

  std::string encode_field_name(std::vector<std::string> names);

  std::string encode_sideblock_name(const std::string &type_string, const std::string &name);

  bool is_sideblock_name(const std::string &name);

  bool use_transformed_storage(const Ioss::Field &field, const std::string &entity_type,
                               const std::string &field_name);

} // namespace Ioad

#endif
