// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "ED_SystemInterface.h"
#include "Tolerance.h"
#include "assembly.h"
#include "exoII_read.h"
#include "exo_block.h"
#include "exodusII.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"
#include "node_set.h"
#include "side_set.h"
#include "smart_assert.h"
#include "stringx.h"

namespace {

  template <typename INT>
  bool Check_Nodal(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const std::vector<INT> &node_map,
                   const INT *id_map, bool check_only);
  template <typename INT>
  bool Check_Element_Block(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                           const std::vector<INT> &elmt_map, const std::vector<INT> &node_map);
  template <typename INT>
  bool Check_Nodeset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                     const std::vector<INT> &node_map, bool check_only);
  template <typename INT>
  bool Check_Sideset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                     const std::vector<INT> &elmt_map, bool check_only);

  template <typename INT>
  bool Check_Element_Block_Params(const Exo_Block<INT> *block1, const Exo_Block<INT> *block2);
  template <typename INT>
  bool Check_Element_Block_Connectivity(Exo_Block<INT> *block1, Exo_Block<INT> *block2,
                                        const std::vector<INT> &elmt_map,
                                        const std::vector<INT> &node_map);
  template <typename INT> bool Check_Assembly(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2);
  template <typename INT>
  bool Check_Assembly_Params(const Assembly<INT> *assembly1, const Assembly<INT> *assembly2);
  bool close_compare(const std::string &st1, const std::string &st2);
} // namespace

template <typename INT> bool Check_Global(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
{
  bool is_same = true;
  if (file1.Dimension() != file2.Dimension()) {
    Warning(".. Dimension doesn't agree.\n");
    is_same = false;
  }
  if (file1.Num_Nodes() != file2.Num_Nodes()) {
    if (interFace.map_flag != MapType::PARTIAL) {
      Warning(".. Number of nodes doesn't agree.\n");
      is_same = false;
    }
  }
  if (file1.Num_Elements() != file2.Num_Elements()) {
    if (interFace.map_flag != MapType::PARTIAL) {
      Warning(".. Number of elements doesn't agree.\n");
      is_same = false;
    }
  }
  if (file1.Num_Element_Blocks() != file2.Num_Element_Blocks()) {
    if (interFace.map_flag != MapType::PARTIAL) {
      Warning(".. Number of element blocks doesn't agree.\n");
      is_same = false;
    }
  }
  if (file1.Num_Times() != file2.Num_Times() && !interFace.quiet_flag && !interFace.ignore_steps) {
    Warning(fmt::format(".. First file has {} result times while the second file has {}.\n",
                        file1.Num_Times(), file2.Num_Times()));
  }
  return is_same;
}

template <typename INT>
void Check_Compatible_Meshes(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, bool check_only,
                             const std::vector<INT> &node_map, const std::vector<INT> &elmt_map,
                             const INT *node_id_map)
{
  bool is_diff = false;
  // NOTE: Check_Global is called earlier. Don't repeat call here.
  if (!Check_Nodal(file1, file2, node_map, node_id_map, check_only)) {
    Warning(".. Differences found in mesh nodal coordinates.\n");
    is_diff = true;
  }

  if (!Check_Element_Block(file1, file2, elmt_map, node_map)) {
    Warning(".. Differences found in element block metadata or connectivity.\n");
    is_diff = true;
  }

  if (!Check_Nodeset(file1, file2, node_map, check_only)) {
    Warning(".. Differences found in node set metadata or node lists.\n");
    is_diff = true;
  }

  if (!Check_Sideset(file1, file2, elmt_map, check_only)) {
    Warning(".. Differences found in side set metadata or side lists.\n");
    is_diff = true;
  }

  if (!Check_Assembly(file1, file2)) {
    Warning(fmt::format(".. Differences found in assembly metadata or assembly entity lists. {}\n",
                        interFace.pedantic ? "" : "[ignored]"));
    if (interFace.pedantic) {
      is_diff = true;
    }
  }

  if (is_diff) {
    Warning(".. Differences found in mesh (non-transient) data.  Aborting...\n");
    exit(1);
  }
}

namespace {
  template <typename INT>
  bool Check_Nodal(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const std::vector<INT> &node_map,
                   const INT *id_map, bool check_only)
  {
    bool is_same = true;

    if (interFace.coord_tol.type == ToleranceMode::IGNORE_ || !check_only) {
      return is_same;
    }

    file1.Load_Nodal_Coordinates();
    file2.Load_Nodal_Coordinates();

    const auto   *x1 = file1.X_Coords();
    const double *y1 = x1;
    const double *z1 = x1;
    if (file1.Dimension() > 1) {
      y1 = file1.Y_Coords();
    }
    if (file1.Dimension() > 2) {
      z1 = file1.Z_Coords();
    }

    const auto   *x2 = file2.X_Coords();
    const double *y2 = x2;
    const double *z2 = x2;
    if (file2.Dimension() > 1) {
      y2 = file2.Y_Coords();
    }
    if (file2.Dimension() > 2) {
      z2 = file2.Z_Coords();
    }

    double max = 0.0;
    double norm;
    for (size_t n = 0; n < file1.Num_Nodes() && (is_same || interFace.show_all_diffs); ++n) {
      // Should this node be processed...
      if (node_map.empty() || node_map[n] >= 0) {
        INT    n2 = node_map.empty() ? n : node_map[n];
        double dx = interFace.coord_tol.Delta(x1[n], x2[n2]);
        if (dx > interFace.coord_tol.value) {
          fmt::print("   x coord {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (node {})\n",
                     interFace.coord_tol.abrstr(), x1[n], x2[n2], dx, (size_t)id_map[n]);
          is_same = false;
        }
        norm = (x1[n] - x2[n2]) * (x1[n] - x2[n2]);

        if (file1.Dimension() > 1 && file2.Dimension() > 1) {
          double dy = interFace.coord_tol.Delta(y1[n], y2[n2]);
          if (dy > interFace.coord_tol.value) {
            fmt::print("   y coord {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (node {})\n",
                       interFace.coord_tol.abrstr(), y1[n], y2[n2], dy, (size_t)id_map[n]);
            is_same = false;
          }
          norm += (y1[n] - y2[n2]) * (y1[n] - y2[n2]);
        }

        if (file1.Dimension() > 2 && file2.Dimension() > 2) {
          double dz = interFace.coord_tol.Delta(z1[n], z2[n2]);
          if (dz > interFace.coord_tol.value) {
            fmt::print("   z coord {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (node {})\n",
                       interFace.coord_tol.abrstr(), z1[n], z2[n2], dz, (size_t)id_map[n]);
            is_same = false;
          }
          norm += (z1[n] - z2[n2]) * (z1[n] - z2[n2]);
        }
        max = max < norm ? norm : max;
      } // End of node iteration...
    }

    if (!interFace.quiet_flag && is_same && max > 0.0) {
      max = sqrt(max);
      fmt::print("Maximum difference between nodal coordinates = {}\n", max);
    }

    file1.Free_Nodal_Coordinates();
    file2.Free_Nodal_Coordinates();
    return is_same;
  }

  template <typename INT>
  bool Check_Element_Block(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                           const std::vector<INT> &elmt_map, const std::vector<INT> &node_map)
  {
    bool is_same = true;
    // Verify that element blocks match in the two files...
    for (size_t b = 0; b < file1.Num_Element_Blocks(); ++b) {
      Exo_Block<INT> *block1 = file1.Get_Element_Block_by_Index(b);
      Exo_Block<INT> *block2 = nullptr;
      if (interFace.map_flag != MapType::DISTANCE && interFace.map_flag != MapType::PARTIAL) {
        if (block1 != nullptr) {
          if (interFace.by_name) {
            block2 = file2.Get_Element_Block_by_Name(block1->Name());
          }
          else {
            block2 = file2.Get_Element_Block_by_Id(block1->Id());
          }
          if (block2 == nullptr) {
            Warning(fmt::format(".. Block id {} with name {} exists in first "
                                "file but not the second.\n",
                                block1->Id(), block1->Name()));
            is_same = false;
          }
          else {
            if (!Check_Element_Block_Params(block1, block2)) {
              is_same = false;
            }
            else {
              // Only do this check if Check_Element_Block_Params does not fail.
              // TODO(gdsjaar): Pass in node_map and node_id_map...
              if (!Check_Element_Block_Connectivity(block1, block2, elmt_map, node_map)) {
                is_same = false;
              }
            }
          }
        }
      }
    }
    return is_same;
  }

  template <typename INT> bool Check_Assembly(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
  {
    bool is_same = true;
    if (file1.Num_Assembly() != file2.Num_Assembly()) {
      Warning(fmt::format(".. The number of assemblies ({}) in the first file does not match the "
                          "number ({}) in the second file.\n",
                          file1.Num_Assembly(), file2.Num_Assembly()));
      is_same = false;
    }

    // Verify that assemblies match in the two files...
    size_t matched = 0;
    for (size_t b = 0; b < file1.Num_Assembly(); ++b) {
      Assembly<INT> *assembly1 = file1.Get_Assembly_by_Index(b);
      if (assembly1 != nullptr) {
        Assembly<INT> *assembly2 = nullptr;
        assembly2                = file2.Get_Assembly_by_Name(assembly1->Name());
        if (assembly2 == nullptr) {
          Warning(fmt::format(".. Assembly '{}' with id {} exists in first "
                              "file but not the second.\n",
                              assembly1->Name(), assembly1->Id()));
          is_same = false;
        }
        else {
          matched++;
          if (!Check_Assembly_Params(assembly1, assembly2)) {
            is_same = false;
          }
        }
      }
    }
    if (matched != file2.Num_Assembly()) {
      for (size_t b = 0; b < file2.Num_Assembly(); ++b) {
        Assembly<INT> *assembly2 = file2.Get_Assembly_by_Index(b);
        if (assembly2 != nullptr) {
          Assembly<INT> *assembly1 = nullptr;
          assembly1                = file1.Get_Assembly_by_Name(assembly2->Name());
          if (assembly1 == nullptr) {
            Warning(fmt::format(".. Assembly '{}' with id {} exists in second "
                                "file but not the first.\n",
                                assembly2->Name(), assembly2->Id()));
            is_same = false;
          }
          else {
            if (!Check_Assembly_Params(assembly2, assembly1)) {
              is_same = false;
            }
          }
        }
      }
    }
    return is_same;
  }

  template <typename INT>
  bool Check_Element_Block_Connectivity(Exo_Block<INT> *block1, Exo_Block<INT> *block2,
                                        const std::vector<INT> &elmt_map,
                                        const std::vector<INT> &node_map)
  {

    bool is_same = true;
    SMART_ASSERT(block1 != nullptr && block2 != nullptr);

    block1->Load_Connectivity();
    block2->Load_Connectivity();
    const auto &conn1 = block1->Connectivity();
    const auto &conn2 = block2->Connectivity();

    SMART_ASSERT(block1->Size() == 0 || block1->Num_Nodes_per_Element() == 0 || !conn1.empty());
    SMART_ASSERT(block2->Size() == 0 || block2->Num_Nodes_per_Element() == 0 || !conn2.empty());

    if (interFace.map_flag == MapType::FILE_ORDER || elmt_map.empty()) {
      size_t node_count = block1->Size() * block1->Num_Nodes_per_Element();
      SMART_ASSERT(node_count == block2->Size() * block2->Num_Nodes_per_Element());

      if (interFace.map_flag != MapType::FILE_ORDER && !node_map.empty()) {
        for (size_t e = 0; e < node_count; ++e) {
          if (node_map[conn1[e] - 1] + 1 != conn2[e]) {
            size_t elem = e / block2->Num_Nodes_per_Element();
            size_t node = e % block2->Num_Nodes_per_Element();
            Warning(
                fmt::format(".. Connectivities in block id {} are not the same.\n"
                            "                  First difference is node {} of local element {}\n",
                            block1->Id(), node + 1, elem + 1));
            is_same = false;
            break;
          }
        }
      }
      else {
        for (size_t e = 0; e < node_count; ++e) {
          if (conn1[e] != conn2[e]) {
            size_t elem = e / block2->Num_Nodes_per_Element();
            size_t node = e % block2->Num_Nodes_per_Element();
            Warning(
                fmt::format(".. Connectivities in block id {} are not the same.\n"
                            "                  First difference is node {} of local element {}\n",
                            block1->Id(), node + 1, elem + 1));
            is_same = false;
            break;
          }
        }
      }
    }
    else {
      auto   offset1     = block1->offset();
      auto   offset2     = block2->offset();
      size_t num_element = block1->Size();
      size_t nnpe        = block1->Num_Nodes_per_Element();
      for (size_t e1 = 0; is_same && e1 < num_element; e1++) {
        for (size_t n = 0; is_same && n < nnpe; ++n) {
          size_t off1 = e1 * nnpe + n;
          auto   e2   = elmt_map[offset1 + e1];
          if (e2 >= 0) { // If doing partial map, not all elements have a match
            e2 -= offset2;
            size_t off2   = e2 * nnpe + n;
            auto   n1     = conn1[off1];
            auto   map_n1 = node_map.empty() ? n1 : node_map[n1 - 1] + 1;
            if (map_n1 != conn2[off2]) {
              Warning(
                  fmt::format(".. Connectivities in block id {} are not the same.\n"
                              "                  First difference is node {} of local element {} "
                              "(file1) {} (file2)\n",
                              block1->Id(), n + 1, e1 + 1, e2 + 1));
              is_same = false;
              break;
            }
          }
        }
      }
    }

    block2->Free_Connectivity();
    block1->Free_Connectivity();
    return is_same;
  }

  template <typename INT>
  bool Check_Assembly_Params(const Assembly<INT> *assembly1, const Assembly<INT> *assembly2)
  {
    bool is_same = true;
    SMART_ASSERT(assembly1 && assembly2);

    if (interFace.by_name && assembly1->Id() != assembly2->Id()) {
      Warning(fmt::format(".. Assembly '{}' ids don't agree ({} != {}).\n", assembly1->Name(),
                          assembly1->Id(), assembly2->Id()));
      is_same = false;
    }
    if (!interFace.by_name && assembly1->Name() != assembly2->Name()) {
      if (!assembly1->generatedName_ && !assembly2->generatedName_) {
        Warning(fmt::format(".. Assembly {} names don't agree ('{}' != '{}').\n", assembly1->Id(),
                            assembly1->Name(), assembly2->Name()));
        is_same = false;
      }
    }
    if (assembly1->Type() != assembly2->Type()) {
      Warning(fmt::format(".. Assembly '{}': entity types don't agree ({} != {}).\n",
                          assembly1->Name(), ex_name_of_object(assembly1->Type()),
                          ex_name_of_object(assembly2->Type())));
      is_same = false;
    }
    if (assembly1->Size() != assembly2->Size()) {
      Warning(fmt::format(".. Assembly '{}': number of entities doesn't agree ({} != {}).\n",
                          assembly1->Name(), assembly1->Entities().size(),
                          assembly2->Entities().size()));
      is_same = false;
    }
    if ((assembly1->Type() == assembly2->Type()) &&
        (assembly1->Entities().size() == assembly2->Entities().size())) {
      // Check membership of the entities list...
      if (!std::is_permutation(assembly1->Entities().begin(), assembly1->Entities().end(),
                               assembly2->Entities().begin())) {
        Warning(fmt::format(".. Assembly '{}': entity list on first file ({}) does not match "
                            "entity list on second file ({}).\n",
                            assembly1->Name(), fmt::join(assembly1->Entities(), ", "),
                            fmt::join(assembly2->Entities(), ", ")));
        is_same = false;
      }
    }

    return is_same;
  }

  template <typename INT>
  bool Check_Element_Block_Params(const Exo_Block<INT> *block1, const Exo_Block<INT> *block2)
  {
    bool is_same = true;
    SMART_ASSERT(block1 && block2);

    if (interFace.by_name && block1->Id() != block2->Id()) {
      Warning(fmt::format(".. Block '{}' ids don't agree ({} != {}).\n", block1->Name(),
                          block1->Id(), block2->Id()));
      if (interFace.pedantic) {
        is_same = false;
      }
    }
    if (!interFace.by_name && block1->Name() != block2->Name()) {
      if (!block1->generatedName_ && !block2->generatedName_) {
        Warning(fmt::format(".. Block {} names don't agree ({} != {}).\n", block1->Id(),
                            block1->Name(), block2->Name()));
        if (interFace.pedantic) {
          is_same = false;
        }
      }
    }
    if (!(no_case_equals(block1->Element_Type(), block2->Element_Type()))) {
      if (!interFace.short_block_check ||
          !close_compare(block1->Element_Type(), block2->Element_Type())) {
        Warning(fmt::format(".. Block {}: element types don't agree ({} != {}).\n", block1->Name(),
                            block1->Element_Type(), block2->Element_Type()));
        is_same = false;
      }
    }
    if (block1->Size() != block2->Size()) {
      Warning(fmt::format(".. Block {}: number of elements doesn't agree ({} != {}).\n",
                          block1->Name(), block1->Size(), block2->Size()));
      is_same = false;
    }
    if (block1->Num_Nodes_per_Element() != block2->Num_Nodes_per_Element()) {
      Warning(fmt::format(".. Block {}: number of nodes per element doesn't agree ({} != {}).\n",
                          block1->Name(), block1->Num_Nodes_per_Element(),
                          block2->Num_Nodes_per_Element()));
      is_same = false;
    }
#if 0
    if (block1->Num_Attributes() != block2->Num_Attributes()) {
      Warning(fmt::format(".. Block {}: number of attributes doesn't agree ({} != {}).\n"
                        block1->Name(), block1->Num_Attributes(), block2->Num_Attributes());
      is_same = false;
    }
#endif
    return is_same;
  }

  template <typename INT>
  bool Check_Nodeset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                     const std::vector<INT> &node_map, bool /*unused*/)
  {
    // Currently don't set diff flag for most of these since we
    // can continue (somewhat) with these differences...
    // As get more usage of nodeset/sideset variables, may revisit
    // what is a diff.
    bool is_same = true;
    if (file1.Num_Node_Sets() != file2.Num_Node_Sets()) {
      if (interFace.map_flag != MapType::PARTIAL) {
        Warning(".. Number of nodesets doesn't agree...\n");
        if (interFace.pedantic) {
          is_same = false;
        }
      }
    }
    // Check that the files both contain the same nodesets...
    for (size_t b = 0; b < file1.Num_Node_Sets(); ++b) {
      Node_Set<INT> *set1 = file1.Get_Node_Set_by_Index(b);
      Node_Set<INT> *set2 = nullptr;
      if (interFace.by_name) {
        set2 = file2.Get_Node_Set_by_Name(set1->Name());
      }
      else {
        set2 = file2.Get_Node_Set_by_Id(set1->Id());
      }

      if (set2 == nullptr) {
        Warning(
            fmt::format(".. Nodeset id {} exists in first file but not the second.\n", set1->Id()));
        if (interFace.pedantic) {
          is_same = false;
        }
      }
      else {
        if (set1->Size() != set2->Size()) {
          Warning(fmt::format(
              ".. The node count for nodeset id {} is not the same in the two files ({} != {}).\n",
              set1->Id(), set1->Size(), set2->Size()));
          if (interFace.pedantic) {
            is_same = false;
          }
        }
        if (interFace.by_name && set1->Id() != set2->Id()) {
          Warning(fmt::format(".. Nodeset '{}' ids don't agree ({} != {}).\n", set1->Name(),
                              set1->Id(), set2->Id()));
          if (interFace.pedantic) {
            is_same = false;
          }
        }
        if (!interFace.by_name && set1->Name() != set2->Name()) {
          if (!set1->generatedName_ && !set2->generatedName_) {
            Warning(fmt::format(".. Nodeset {} names don't agree ({} != {}).\n", set1->Id(),
                                set1->Name(), set2->Name()));
            if (interFace.pedantic) {
              is_same = false;
            }
          }
        }
      }
    }

    // Check that can access all nodesets in file2.
    // This should never fail if the above tests pass...
    for (size_t b = 0; b < file2.Num_Node_Sets(); ++b) {
      Node_Set<INT> *set2 = file2.Get_Node_Set_by_Index(b);
      if (set2 == nullptr) {
        Warning(
            fmt::format(".. Could not access the Nodeset with index {} in the second file.\n", b));
        if (interFace.pedantic) {
          is_same = false;
        }
      }
    }

    // Do the following check(s) only if there are nodeset variables...
    // For each nodeset, check that the order of the nodeset nodes is the same.
    // Eventually need to be able to map the order...
    if (!interFace.ns_var_names.empty() || interFace.pedantic) {
      for (size_t b = 0; b < file1.Num_Node_Sets(); ++b) {
        Node_Set<INT> *set1 = file1.Get_Node_Set_by_Index(b);
        Node_Set<INT> *set2 = nullptr;
        if (interFace.by_name) {
          set2 = file2.Get_Node_Set_by_Name(set1->Name());
        }
        else {
          set2 = file2.Get_Node_Set_by_Id(set1->Id());
        }

        if (set2 == nullptr) {
          continue;
        }

        if (!node_map.empty()) {
          set1->apply_map(node_map);
        }

        if ((interFace.pedantic || set1->var_count() > 0) && (set1->Size() == set2->Size())) {
          size_t  node_count = set1->Size();
          int64_t diff       = -1;
          for (size_t i = 0; i < node_count; i++) {
            if (set1->Node_Id(i) != set2->Node_Id(i)) {
              diff = i;
              break;
            }
          }
          if (diff >= 0) {
            Warning(fmt::format(
                ".. The nodelists for nodeset id {} are not the same in the two files.\n"
                "\t\tThe first difference is at position {}: Node {} vs. Node {}.\n",
                set1->Id(), set1->Node_Index(diff) + 1, set1->Node_Id(diff), set2->Node_Id(diff)));
            if (interFace.map_flag != MapType::PARTIAL) {
              is_same = false;
            }
            else {
              Warning(".. The nodelist differences are ignored for the "
                      "partial_map case.\n");
            }
          }
        }
      }
    }
    return is_same;
  }

  template <typename INT>
  bool Check_Sideset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                     const std::vector<INT> &elmt_map, bool /*unused*/)
  {
    // Currently don't set diff flag for most of these since we
    // can continue (somewhat) with these differences...
    // As get more usage of nodeset/sideset variables, may revisit
    // what is a diff.
    bool is_same = true;
    if (file1.Num_Side_Sets() != file2.Num_Side_Sets()) {
      if (interFace.map_flag != MapType::PARTIAL) {
        Warning(".. Number of sidesets doesn't agree...\n");
        if (interFace.pedantic) {
          is_same = false;
        }
      }
    }
    // Check that the files both contain the same sidesets...
    for (size_t b = 0; b < file1.Num_Side_Sets(); ++b) {
      Side_Set<INT> *set1 = file1.Get_Side_Set_by_Index(b);
      Side_Set<INT> *set2 = nullptr;
      if (interFace.by_name) {
        set2 = file2.Get_Side_Set_by_Name(set1->Name());
      }
      else {
        set2 = file2.Get_Side_Set_by_Id(set1->Id());
      }

      if (set2 == nullptr) {
        Warning(
            fmt::format(".. Sideset id {} exists in first file but not the second.\n", set1->Id()));
        if (interFace.pedantic) {
          is_same = false;
        }
      }
      else {
        if (set1->Size() != set2->Size()) {
          Warning(fmt::format(
              ".. The side count for sideset id {} is not the same in the two files ({} != {}).\n",
              set1->Id(), set1->Size(), set2->Size()));
          if (interFace.pedantic) {
            is_same = false;
          }
        }
        if (interFace.by_name && set1->Id() != set2->Id()) {
          Warning(fmt::format(".. Sideset '{}' ids don't agree ({} != {}).\n", set1->Name(),
                              set1->Id(), set2->Id()));
          if (interFace.pedantic) {
            is_same = false;
          }
        }
        if (!interFace.by_name && set1->Name() != set2->Name()) {
          if (!set1->generatedName_ && !set2->generatedName_) {
            Warning(fmt::format(".. Sideset {} names don't agree ({} != {}).\n", set1->Id(),
                                set1->Name(), set2->Name()));
            if (interFace.pedantic) {
              is_same = false;
            }
          }
        }
      }
    }

    for (size_t b = 0; b < file2.Num_Side_Sets(); ++b) {
      Side_Set<INT> *set2 = file2.Get_Side_Set_by_Index(b);
      if (set2 == nullptr) {
        Warning(
            fmt::format(".. Could not access the Sideset with index {} in the second file.\n", b));
        if (interFace.pedantic) {
          is_same = false;
        }
      }
    }

    // Do the following check(s) only if there are sideset variables... (or -pedantic)
    // For each sideset, check that the order of the sideset sides is the same.
    // Eventually need to be able to map the order...
    if (!interFace.ss_var_names.empty() || interFace.pedantic || !interFace.ignore_sideset_df) {
      for (size_t b = 0; b < file1.Num_Side_Sets(); ++b) {
        Side_Set<INT> *set1 = file1.Get_Side_Set_by_Index(b);
        Side_Set<INT> *set2 = nullptr;
        if (interFace.by_name) {
          set2 = file2.Get_Side_Set_by_Name(set1->Name());
        }
        else {
          set2 = file2.Get_Side_Set_by_Id(set1->Id());
        }

        if (set2 == nullptr) {
          continue;
        }

        if (!elmt_map.empty()) {
          set1->apply_map(elmt_map);
        }

        // Don't care if sidesets don't match if there are no variables...
        // If different sizes and pedantic, difference caught above.
        if ((interFace.pedantic || set1->var_count() > 0) && (set1->Size() == set2->Size())) {
          size_t  side_count = set1->Size();
          int64_t diff       = -1;
          for (size_t i = 0; i < side_count; i++) {
            if (set1->Side_Id(i) != set2->Side_Id(i)) {
              diff = i;
              break;
            }
          }
          if (diff >= 0) {
            // If `elmt_map` is not null, then need to unmap the set1 ids to get the local id that
            // appears in the file.  If don't do this, error message is very confusing for the
            // user...
            auto set1_id = set1->Side_Id(diff).first;
            if (!elmt_map.empty()) {
              // Iterate map to find an entry equal to `set1_id`.  Its position is then the file1 id
              // of the element.
              for (size_t i = 0; i < file1.Num_Elements(); i++) {
                if (elmt_map[i] == set1_id - 1) {
                  set1_id = i + 1;
                  break;
                }
              }
            }

            Warning(fmt::format(
                ".. The sidelists for sideset id {} are not the same in the two files.\n"
                "\t\tThe first difference is at position {}: Side {}.{} .vs. Side {}.{}.\n",
                set1->Id(), set1->Side_Index(diff) + 1, set1_id, set1->Side_Id(diff).second,
                set2->Side_Id(diff).first, set2->Side_Id(diff).second));
            if (interFace.map_flag != MapType::PARTIAL) {
              is_same = false;
            }
            else {
              Warning(".. The sidelist differences are ignored for the "
                      "partial_map case.\n");
            }
          }
        }
      }
    }
    return is_same;
  }

  bool close_compare(const std::string &st1, const std::string &st2)
  {
    auto len1 = st1.size();
    auto len2 = st2.size();

    // Check that digits (if any) at end of names match
    while ((isdigit(st1[len1 - 1]) != 0) && (isdigit(st2[len2 - 1]) != 0)) {
      if (st1[len1 - 1] != st2[len2 - 1]) {
        return false;
      }
      len1--;
      len2--;
    }

    // Skip any digits at the end.  It's OK if only one name has
    // digits since we want 'tri' and 'triangle6' to match, but not
    // tri6 and tri3 nor quad9 and tri9
    while (isdigit(st1[len1 - 1]) != 0) {
      len1--;
    }

    while (isdigit(st2[len2 - 1]) != 0) {
      len2--;
    }

    unsigned length = (len1 < len2) ? len1 : len2;
    for (unsigned i = 0; i < length; ++i) {
      if (toupper(st1[i]) != toupper(st2[i])) {
        return false;
      }
    }
    return true;
  }
} // namespace

template bool Check_Global(ExoII_Read<int> &file1, ExoII_Read<int> &file2);
template void Check_Compatible_Meshes(ExoII_Read<int> &file1, ExoII_Read<int> &file2,
                                      bool check_only, const std::vector<int> &node_map,
                                      const std::vector<int> &elmt_map, const int *node_id_map);

template bool Check_Global(ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2);
template void Check_Compatible_Meshes(ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2,
                                      bool check_only, const std::vector<int64_t> &node_map,
                                      const std::vector<int64_t> &elmt_map,
                                      const int64_t              *node_id_map);
