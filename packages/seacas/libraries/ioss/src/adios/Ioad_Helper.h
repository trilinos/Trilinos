// Copyright(C) 1999-2010 National Technology & Engineering Solutions
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
