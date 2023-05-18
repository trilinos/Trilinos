// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_ASSEMBLYUTILS_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_ASSEMBLYUTILS_HPP_

#include <memory>                                     // for make_shared
#include <cstddef>                                    // for size_t
#include <cstdint>                                    // for int64_t, uint64_t
#include <unistd.h>                                   // for unlink
#include <iostream>                                   // for operator<<, bas...
#include <limits>                                     // for numeric_limits
#include <string>                                     // for string, basic_s...
#include <typeinfo>                                   // for type_info
#include <utility>                                    // for move, pair
#include <stdexcept>                                  // for logic_error
#include <vector>                                     // for vector, swap
#include <gtest/gtest.h>
#include "stk_util/environment/RuntimeWarning.hpp"    // for RuntimeWarningA...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg
#include "stk_util/util/TreeTraverser.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk
{
namespace unit_test_util
{
struct AssemblyDescription {
  static constexpr unsigned InvalidId = std::numeric_limits<unsigned>::max();

  AssemblyDescription(unsigned id_, const std::vector<std::string>& group)
    : id(id_)
    , partNames(group)
  {
    STK_ThrowRequireMsg(id != InvalidId, "Invalid id: " << id);
    stk::util::sort_and_unique(partNames);
  }

  AssemblyDescription()
    : id(InvalidId)
  {
  }

  void add_part_name(const std::string& partName)
  {
    partNames.push_back(partName);
    stk::util::sort_and_unique(partNames);
  }

  void add_part_names(const std::vector<std::string>& parts)
  {
    partNames.insert(partNames.end(), parts.begin(), parts.end());
    stk::util::sort_and_unique(partNames);
  }

  unsigned id;
  std::vector<std::string> partNames;
  std::vector<unsigned> subAssemblies;
};

template <typename GraphKey, typename GraphData>
class AssemblyGraph : public stk::util::TreeTraverser<GraphKey, typename GraphData::const_iterator> {
 public:
  using Key = GraphKey;
  using Iterator = typename GraphData::const_iterator;

  AssemblyGraph(const GraphData& assemblyData)
    : stk::util::TreeTraverser<GraphKey, Iterator>()
    , m_assemblyData(assemblyData)
  {
  }

  Iterator begin() const override { return m_assemblyData.begin(); }

  Iterator end() const override { return m_assemblyData.end(); }

  const Key& get_node(Iterator iter) const override { return iter->first; }

  size_t size() const override { return m_assemblyData.size(); }

  bool node_exists(const Key& node) const override { return m_assemblyData.find(node) != m_assemblyData.end(); }

  const std::vector<Key>& children(const Key& node) const override
  {
    const auto& iter = m_assemblyData.find(node);
    if(iter != m_assemblyData.end()) {
      return iter->second.subAssemblies;
    }

    return m_emptyChildren;
  }

 private:
  AssemblyGraph() = delete;
  AssemblyGraph(const AssemblyGraph&) = delete;

  const GraphData& m_assemblyData;
  const std::vector<Key> m_emptyChildren;
};

class AssemblyManager {
 public:
  AssemblyManager();
  AssemblyManager(const std::string& baseName);
  AssemblyManager(const std::string& baseName, unsigned startId);

  bool insert_sub_assembly_without_creating_cyclic_graph(AssemblyDescription& parentAssembly, unsigned childId);

  bool add_sub_assembly(unsigned parent, unsigned child);

  std::vector<unsigned> add_alternating_pairwise_assemblies(const std::vector<std::string>& names);

  std::vector<unsigned> add_alternating_part_assemblies(const std::vector<std::string>& names);

  std::vector<unsigned> add_assembly_for_alternating_parts(const std::vector<std::string>& names);

  std::vector<unsigned> add_assembly_per_part(const std::vector<std::string>& names);

  unsigned add_assembly(const std::vector<std::string>& names);

  unsigned add_assembly(const std::string& name);

  const std::unordered_map<unsigned, AssemblyDescription>& get_assembly_descriptions() const;

  std::vector<std::string> get_assembly_names() const;

  std::vector<std::string> get_all_part_names() const;

  const std::string& get_assembly_name(unsigned id);

  const AssemblyDescription& get_assembly(unsigned id) const;

  AssemblyDescription& get_assembly(unsigned id);

  void add_part_names_to_assembly(unsigned id, const std::vector<std::string>& partNames);

  void add_part_name_to_assembly(unsigned id, const std::string& partName);

  std::vector<unsigned> get_assembly_reverse_traversal_list(unsigned id);

  std::vector<std::string> get_assembly_tree_part_names(unsigned id);

  unsigned size(unsigned id) const;

  void print_assembly(const std::string& preamble, unsigned id, std::ostream& out = std::cout);
  void print_assemblies(const std::string& preamble, std::ostream& out = std::cout);

 private:
  std::string create_assembly_name(const unsigned id);

  unsigned assign(const std::vector<std::string>& group);

  bool is_assigned(unsigned id) const;

  unsigned get_unassigned_id() const;

  std::string get_default_base_name() const;

  unsigned get_default_start_id() const;

  std::string m_baseName;
  unsigned m_startId;
  std::unordered_map<unsigned, AssemblyDescription> m_assemblies;
  mutable std::unordered_map<unsigned, std::string> m_parts;

  AssemblyGraph<unsigned, std::unordered_map<unsigned, AssemblyDescription> > m_graph;
};


unsigned create_deep_assembly(AssemblyManager& transferAssembly);

} // namespace unit_test_util
} // namespace stk




#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_ASSEMBLYUTILS_HPP_ */
