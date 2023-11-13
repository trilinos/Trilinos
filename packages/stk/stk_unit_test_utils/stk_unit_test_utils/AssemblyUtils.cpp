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

#include "AssemblyUtils.hpp"

namespace stk
{
namespace unit_test_util
{
  AssemblyManager::AssemblyManager()
    : m_baseName(get_default_base_name())
    , m_startId(get_default_start_id())
    , m_graph(m_assemblies)
  {
  }

  AssemblyManager::AssemblyManager(const std::string& baseName)
    : m_baseName(baseName)
    , m_startId(get_default_start_id())
    , m_graph(m_assemblies)
  {
  }

  AssemblyManager::AssemblyManager(const std::string& baseName, unsigned startId)
    : m_baseName(baseName)
    , m_startId(startId)
    , m_graph(m_assemblies)
  {
  }

  bool AssemblyManager::insert_sub_assembly_without_creating_cyclic_graph(AssemblyDescription& parentAssembly,
                                                                          unsigned childId)
  {
    parentAssembly.subAssemblies.push_back(childId);
    bool isCyclic = m_graph.is_cyclic();

    if(isCyclic) {
      parentAssembly.subAssemblies.pop_back();
    }

    return !isCyclic;
  }

  bool AssemblyManager::add_sub_assembly(unsigned parent, unsigned child)
  {
    auto parentIter = m_assemblies.find(parent);
    if(parentIter == m_assemblies.end()) return false;

    if(m_assemblies.find(child) == m_assemblies.end()) return false;

    AssemblyDescription& parentAssembly = parentIter->second;
    if(std::find(parentAssembly.subAssemblies.begin(), parentAssembly.subAssemblies.end(), child) !=
       parentAssembly.subAssemblies.end()) {
      return true;
    }

    if(!insert_sub_assembly_without_creating_cyclic_graph(parentAssembly, child)) {
      stk::RuntimeWarning() << "Attempting to create cyclic dependency by adding assembly: '"
                            << get_assembly_name(child) << "' as a sub-assembly to: '" << get_assembly_name(parent)
                            << "'";
      return false;
    }

    stk::util::sort_and_unique(parentIter->second.subAssemblies);
    return true;
  }

  std::vector<unsigned>
  AssemblyManager::add_alternating_pairwise_assemblies(const std::vector<std::string>& names)
  {
    std::vector<unsigned> ids;

    if(names.size() < 3) {
      std::cout << "Not enough part names to form alternating pairwise assemblies" << std::endl;
      return ids;
    }

    unsigned numPairs = names.size() - 1;
    ids.reserve(numPairs);

    for(unsigned i = 0; i < numPairs; i++) {
      STK_ThrowRequireMsg(names[i] != names[i + 1],
                      "interp parts {" << names[i] << "," << names[i + 1] << "} must be unique");
      unsigned id = assign(std::vector<std::string>{ names[i], names[i + 1] });
      ids.push_back(id);
    }

    return ids;
  }

  std::vector<unsigned>
  AssemblyManager::add_alternating_part_assemblies(const std::vector<std::string>& names)
  {
    std::vector<unsigned> ids;

    if(names.empty()) {
      std::cout << "Not enough part names to form alternating part assembly" << std::endl;
      return ids;
    }

    for(unsigned i = 0; i < names.size(); i += 2) {
      unsigned id = assign(std::vector<std::string>{ names[i] });
      ids.push_back(id);
    }

    return ids;
  }

  std::vector<unsigned>
  AssemblyManager::add_assembly_for_alternating_parts(const std::vector<std::string>& names)
  {
    std::vector<unsigned> ids;

    if(names.empty()) {
      std::cout << "Not enough part names to form alternating part assembly" << std::endl;
      return ids;
    }

    std::vector<std::string> assemblyPartNames;

    for(unsigned i = 0; i < names.size(); i += 2) {
      assemblyPartNames.push_back(names[i]);
    }

    unsigned id = assign(assemblyPartNames);
    ids.push_back(id);
    return ids;
  }

  std::vector<unsigned>
  AssemblyManager::add_assembly_per_part(const std::vector<std::string>& names)
  {
    std::vector<unsigned> ids;
    for(unsigned i = 0; i < names.size(); i++) {
      unsigned id = assign(std::vector<std::string>{ names[i] });
      ids.push_back(id);
    }
    return ids;
  }

  unsigned AssemblyManager::add_assembly(const std::vector<std::string>& names) { return assign(names); }

  unsigned AssemblyManager::add_assembly(const std::string& name)
  {
    return add_assembly(std::vector<std::string>{ name });
  }

  const std::unordered_map<unsigned, AssemblyDescription>&
  AssemblyManager::get_assembly_descriptions() const
  { return m_assemblies; }

  std::vector<std::string> AssemblyManager::get_assembly_names() const
  {
    std::vector<std::string> names;
    names.reserve(m_parts.size());

    for(const auto& part : m_parts) {
      names.push_back(part.second);
    }

    return names;
  }

  std::vector<std::string> AssemblyManager::get_all_part_names() const
  {
    std::vector<std::string> names;

    for(const auto& assembly : m_assemblies) {
      names.insert(names.end(), assembly.second.partNames.begin(), assembly.second.partNames.end());
    }

    stk::util::sort_and_unique(names);
    return names;
  }

  const std::string& AssemblyManager::get_assembly_name(unsigned id)
  {
    const auto iter = m_parts.find(id);
    STK_ThrowRequireMsg(iter != m_parts.end(), "Could not find assembly with id: " << id);
    return iter->second;
  }

  const AssemblyDescription& AssemblyManager::get_assembly(unsigned id) const
  {
    const auto iter = m_assemblies.find(id);
    STK_ThrowRequireMsg(iter != m_assemblies.end(), "Could not find assembly with id: " << id);
    return iter->second;
  }

  AssemblyDescription& AssemblyManager::get_assembly(unsigned id)
  {
    auto iter = m_assemblies.find(id);
    STK_ThrowRequireMsg(iter != m_assemblies.end(), "Could not find assembly with id: " << id);
    return iter->second;
  }

  void AssemblyManager::add_part_names_to_assembly(unsigned id, const std::vector<std::string>& partNames)
  {
    AssemblyDescription& assembly = get_assembly(id);
    assembly.add_part_names(partNames);
  }

  void AssemblyManager::add_part_name_to_assembly(unsigned id, const std::string& partName)
  {
    AssemblyDescription& assembly = get_assembly(id);
    assembly.add_part_name(partName);
  }

  std::vector<unsigned> AssemblyManager::get_assembly_reverse_traversal_list(unsigned id)
  {
    return m_graph.get_reverse_traversal_list(id);
  }

  std::vector<std::string> AssemblyManager::get_assembly_tree_part_names(unsigned id)
  {
    std::vector<std::string> treePartNames;
    std::vector<unsigned> traversalList = m_graph.get_forward_traversal_list(id);

    for(unsigned leafId : traversalList) {
      AssemblyDescription& assembly = get_assembly(leafId);
      treePartNames.insert(treePartNames.end(), assembly.partNames.begin(), assembly.partNames.end());
    }

    stk::util::sort_and_unique(treePartNames);
    return treePartNames;
  }

  unsigned AssemblyManager::size(unsigned id) const
  {
    const AssemblyDescription& assembly = get_assembly(id);
    return assembly.partNames.size() + assembly.subAssemblies.size();
  }

  void AssemblyManager::print_assembly(const std::string& preamble, unsigned id, std::ostream& out)
  {
    const AssemblyDescription& assembly = get_assembly(id);
    const std::string& assemblyName = get_assembly_name(assembly.id);

    out << preamble << assemblyName << " {";

    size_t numNames = assembly.partNames.size();
    for(size_t i = 0; i < numNames; ++i) {
      out << assembly.partNames[i];

      if(i != numNames - 1) {
        out << ",";
      }
    }

    if(!assembly.subAssemblies.empty()) {
      out << " [";
      size_t numSubAssemblies = assembly.subAssemblies.size();
      for(size_t i = 0; i < numSubAssemblies; ++i) {
        const AssemblyDescription& subAssembly = get_assembly(assembly.subAssemblies[i]);
        const std::string& subAssemblyName = get_assembly_name(subAssembly.id);

        out << subAssemblyName;

        if(i != numSubAssemblies - 1) {
          out << ",";
        }
      }
      out << "] ";
    }

    out << "} " << std::endl;
  }

  void AssemblyManager::print_assemblies(const std::string& preamble, std::ostream& out)
  {
    for(const auto& entry : m_assemblies) {
      print_assembly(preamble, entry.first, out);
    }
  }

  std::string AssemblyManager::create_assembly_name(const unsigned id)
  {
    std::string assemblyName = m_baseName;
    assemblyName += "_";
    assemblyName += std::to_string(id);
    return assemblyName;
  }

  unsigned AssemblyManager::assign(const std::vector<std::string>& group)
  {
    AssemblyDescription assembly(get_unassigned_id(), group);
    std::string assemblyName = create_assembly_name(assembly.id);

    m_assemblies[assembly.id] = assembly;
    m_parts[assembly.id] = assemblyName;

    return assembly.id;
  }

  bool AssemblyManager::is_assigned(unsigned id) const { return m_parts.count(id) > 0; }

  unsigned AssemblyManager::get_unassigned_id() const
  {
    unsigned nextPartId = m_startId;
    while(is_assigned(nextPartId))
      nextPartId++;
    return nextPartId;
  }

  std::string AssemblyManager::get_default_base_name() const { return "assembly"; }

  unsigned AssemblyManager::get_default_start_id() const { return 9000; }


  unsigned create_deep_assembly(AssemblyManager& transferAssembly)
  {
    unsigned parentAssembly = transferAssembly.add_assembly(std::vector<std::string>{});
    unsigned subAssembly1 = transferAssembly.add_assembly(std::vector<std::string>{});
    unsigned subAssembly2 = transferAssembly.add_assembly(std::vector<std::string>{});

    EXPECT_TRUE(transferAssembly.add_sub_assembly(parentAssembly, subAssembly1));
    EXPECT_TRUE(transferAssembly.add_sub_assembly(parentAssembly, subAssembly2));

    unsigned subSubAssembly = transferAssembly.add_assembly(std::vector<std::string>{});
    unsigned subSubSubAssembly1 = transferAssembly.add_assembly(std::vector<std::string>{});
    unsigned subSubSubAssembly2 = transferAssembly.add_assembly(std::vector<std::string>{});

    EXPECT_TRUE(transferAssembly.add_sub_assembly(subSubAssembly, subSubSubAssembly1));
    EXPECT_TRUE(transferAssembly.add_sub_assembly(subSubAssembly, subSubSubAssembly2));

    EXPECT_EQ(2u, transferAssembly.size(parentAssembly));
    EXPECT_TRUE(transferAssembly.add_sub_assembly(subAssembly1, subSubAssembly));
    transferAssembly.add_part_name_to_assembly(subAssembly2, "block_3");

    EXPECT_EQ(2u, transferAssembly.size(subSubAssembly));
    transferAssembly.add_part_names_to_assembly(subSubSubAssembly1, { "block_1", "block_2" });
    transferAssembly.add_part_names_to_assembly(subSubSubAssembly2, { "block_4" });

    return parentAssembly;
  }

} // namespace unit_test_util
} // namespace stk


