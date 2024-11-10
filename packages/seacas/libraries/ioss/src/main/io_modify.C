// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ionit_Initializer.h"
#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_FileInfo.h"
#include "Ioss_Getline.h"
#include "Ioss_Glob.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StructuredBlock.h"
#include "Ioss_Utils.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <stdint.h>
#include <string>
#include <tokenize.h>
#include <unistd.h>
#include <vector>

#include "Ioss_EntityType.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"
#include "modify_interface.h"

// `io_modify` is only built if Exodus is enabled...
#include "exodus/Ioex_Internals.h"
#include "exodus/Ioex_Utils.h"
#include <exodusII.h>

#if defined(SEACAS_HAVE_CGNS)
#include "cgns/Iocgns_Utils.h"
#endif

#if defined(__IOSS_WINDOWS__)
#include <io.h>

#define isatty _isatty
#endif
// ========================================================================

using real = double;

namespace {
  std::string codename;
  std::string version = "2.07 (2024-04-15)";

  std::vector<Ioss::GroupingEntity *> attributes_modified;

  Ioss::EntityType get_entity_type(const std::string &type)
  {
    if (Ioss::Utils::substr_equal(type, "elementblock")) {
      return Ioss::ELEMENTBLOCK;
    }
    else if (Ioss::Utils::substr_equal(type, "block")) {
      return Ioss::ELEMENTBLOCK;
    }
    else if (Ioss::Utils::substr_equal(type, "structuredblock")) {
      return Ioss::STRUCTUREDBLOCK;
    }
    else if (Ioss::Utils::substr_equal(type, "nodeblock")) {
      return Ioss::NODEBLOCK;
    }
    else if (Ioss::Utils::substr_equal(type, "nodeset")) {
      return Ioss::NODESET;
    }
    else if (Ioss::Utils::substr_equal(type, "nset")) {
      return Ioss::NODESET;
    }
    else if (Ioss::Utils::substr_equal(type, "nodelist")) {
      return Ioss::NODESET;
    }
    else if (Ioss::Utils::substr_equal(type, "sideset")) {
      return Ioss::SIDESET;
    }
    else if (Ioss::Utils::substr_equal(type, "sset")) {
      return Ioss::SIDESET;
    }
    else if (Ioss::Utils::substr_equal(type, "surface")) {
      return Ioss::SIDESET;
    }
    else if (Ioss::Utils::substr_equal(type, "assembly")) {
      return Ioss::ASSEMBLY;
    }
    else if (Ioss::Utils::substr_equal(type, "blob")) {
      return Ioss::BLOB;
    }
    return Ioss::INVALID_TYPE;
  }

  Ioss::NameList get_name_list(const Ioss::Region &region, Ioss::EntityType type);
  void           handle_help(const std::string &topic);
  bool           handle_delete(const Ioss::NameList &tokens, Ioss::Region &region);
  void           handle_list(const Ioss::NameList &tokens, const Ioss::Region &region,
                             bool show_attribute = false);
  void           handle_graph(const Ioss::NameList &tokens, const Ioss::Region &region);
  bool handle_assembly(const Ioss::NameList &tokens, Ioss::Region &region, bool allow_modify);
  bool handle_attribute(const Ioss::NameList &tokens, Ioss::Region &region);
  bool handle_geometry(const Ioss::NameList &tokens, Ioss::Region &region);
  bool handle_time(const Ioss::NameList &tokens, Ioss::Region &region);
  bool handle_rename(const Ioss::NameList &tokens, Ioss::Region &region);
  void update_assembly_info(Ioss::Region &region, const Modify::Interface &interFace);

  void modify_time(Ioss::Region &region, double scale, double offset);

  void offset_filtered_coordinates(Ioss::Region &region, real offset[3],
                                   const std::vector<const Ioss::GroupingEntity *> &blocks);
  void scale_filtered_coordinates(Ioss::Region &region, real scale[3],
                                  const std::vector<const Ioss::GroupingEntity *> &blocks);
  void rotate_filtered_coordinates(Ioss::Region &region, real rotation_matrix[3][3],
                                   const std::vector<const Ioss::GroupingEntity *> &blocks);
  bool update_rotation_matrix(real rotation_matrix[3][3], const std::string &axis, double angle);

  void set_db_properties(const Modify::Interface &interFace, Ioss::DatabaseIO *dbi);

  void info_entity(const Ioss::StructuredBlock *sb, bool show_property = false);
  void info_entity(const Ioss::NodeBlock *nb, bool show_property = false);
  void info_entity(const Ioss::ElementBlock *eb, bool show_property = false);
  void info_entity(const Ioss::NodeSet *ns, bool show_property = false);
  void info_entity(const Ioss::SideSet *ss, bool show_property = false);
  void info_entity(const Ioss::Assembly *as, bool show_property = false);
  void info_entity(const Ioss::Blob *blob, bool show_property = false);
  void info_entity(const Ioss::Region &region, bool show_property = false);
  void info_time(const Ioss::Region &region);

  template <typename T>
  void info_entities(const std::vector<T *> &entities, const Ioss::NameList &tokens,
                     const Ioss::Region &region, const std::string &type,
                     bool show_property = false)
  {
    if (entities.empty()) {
      fmt::print("\n\t*** There are no {} in this model.\n", type);
      return;
    }
    if (tokens.size() == 4 && Ioss::Utils::substr_equal(tokens[2], "matches")) {
      //   0       1           2       3
      // LIST {entity_type} MATCHES {regex}
      auto entity_type = get_entity_type(tokens[1]);
      if (entity_type == Ioss::INVALID_TYPE) {
        fmt::print(stderr, fg(fmt::color::yellow), "WARNING: Unrecognized entity type '{}'.\n",
                   tokens[1]);
      }
      std::regex     reg(tokens[3], std::regex::extended);
      Ioss::NameList names = get_name_list(region, entity_type);

      // Check for match against all names in list...
      bool matched = false;
      for (const auto &name : names) {
        if (std::regex_match(name, reg)) {
          const auto *entity = region.get_entity(name, entity_type);
          const T    *ge     = dynamic_cast<const T *>(entity);
          if (ge != nullptr) {
            info_entity(ge, show_property);
            matched = true;
          }
        }
      }
      if (!matched) {
        fmt::print(stderr, fg(fmt::color::yellow),
                   "WARNING: Regular Expression '{}' did not match any {}\n", tokens[3], type);
      }
    }
    else if (tokens.size() == 4 && Ioss::Utils::substr_equal(tokens[2], "glob")) {
      //   0       1           2       3
      // LIST {entity_type} GLOB {glob}
      auto entity_type = get_entity_type(tokens[1]);
      if (entity_type == Ioss::INVALID_TYPE) {
        fmt::print(stderr, fg(fmt::color::yellow), "WARNING: Unrecognized entity type '{}'.\n",
                   tokens[1]);
      }
      Ioss::glob::glob glob(tokens[3]);
      Ioss::NameList   names = get_name_list(region, entity_type);

      // Check for match against all names in list...
      bool matched = false;
      for (const auto &name : names) {
        if (Ioss::glob::glob_match(name, glob)) {
          const auto *entity = region.get_entity(name, entity_type);
          const T    *ge     = dynamic_cast<const T *>(entity);
          if (ge != nullptr) {
            info_entity(ge, show_property);
            matched = true;
          }
        }
      }
      if (!matched) {
        fmt::print(stderr, fg(fmt::color::yellow),
                   "WARNING: Glob Expression '{}' did not match any {}\n", tokens[3], type);
      }
    }
    else if (tokens.size() == 2 ||
             (tokens.size() == 3 && Ioss::Utils::substr_equal(tokens[2], "list"))) {
      for (const T *ge : entities) {
        info_entity(ge, show_property);
      }
    }
    else {
      for (size_t i = 2; i < tokens.size(); i++) {
        const T *ge = dynamic_cast<T *>(region.get_entity(tokens[i]));
        if (ge != nullptr) {
          info_entity(ge, show_property);
        }
      }
    }
  }

  std::string name(const Ioss::GroupingEntity *entity)
  {
    return entity->type_string() + " '" + entity->name() + "'";
  }

  int64_t id(const Ioss::GroupingEntity *entity)
  {
    int64_t id = entity->get_optional_property("id", -1);
    return id;
  }

  int64_t get_next_assembly_id(const Ioss::Region &region)
  {
    // Determines the current maximum assembly id and then
    // returns ids larger than that value.  The ids it
    // returns are multiples of 100 just because.
    static int64_t next_id = 0;
    if (next_id == 0) {
      const auto &assemblies = region.get_assemblies();
      for (const auto *assembly : assemblies) {
        auto my_id = id(assembly);
        next_id    = std::max(next_id, my_id);
      }
      next_id = (next_id / 100);
    }
    ++next_id;
    return next_id * 100;
  }

  Ioss::PropertyManager set_properties(const Modify::Interface & /* interFace */)
  {
    Ioss::PropertyManager properties{};
    return properties;
  }
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  fmt::print(fg(fmt::color::cyan), "\n *** {}, Version {}\n", codename, version);
  Modify::Interface interFace;
  interFace.parse_options(argc, argv);

  Ioss::Init::Initializer io;

  std::string inpfile    = interFace.filename();
  std::string input_type = interFace.type();

  //========================================================================
  // INPUT ...
  // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
  //========================================================================
  Ioss::PropertyManager properties = set_properties(interFace);
  properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_MODIFY));

  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::WRITE_RESTART,
                                                  Ioss::ParallelUtils::comm_world(), properties);
  if (dbi == nullptr || !dbi->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  set_db_properties(interFace, dbi);

  // NOTE: 'region' owns 'db' pointer at this time...
  Ioss::Region region(dbi, "region");
  region.begin_mode(Ioss::STATE_DEFINE_MODEL);

  region.output_summary(std::cout, true);

  bool from_term    = (isatty(0) != 0 && isatty(1) != 0);
  bool changed      = false;
  bool allow_modify = interFace.modify_existing_assembly();

  while (true) {
    std::string input;
    if (from_term) {
      fmt::print(fg(fmt::terminal_color::magenta), "\n");
      const char *cinput = Ioss::getline_int("COMMAND> ");
      if (cinput && cinput[0] == '\0') {
        break;
      }
      if (cinput) {
        Ioss::gl_histadd(cinput);
      }
      input = cinput;
    }
    else {
      std::getline(std::cin, input);
    }

    // NOTE: getline_int returns the trailing '\n'
    auto tokens = Ioss::tokenize(input, " ,\n");
    if (tokens.empty()) {
      continue;
    }
    if (Ioss::Utils::substr_equal(tokens[0], "exit") ||
        Ioss::Utils::substr_equal(tokens[0], "end")) {
      break;
    }
    if (Ioss::Utils::substr_equal(tokens[0], "quit")) {
      changed = false;
      break;
    }

    if (Ioss::Utils::str_equal(tokens[0], "help")) {
      handle_help(tokens.back());
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "allow") &&
             Ioss::Utils::substr_equal(tokens[1], "modifications")) {
      allow_modify = true;
      fmt::print(fg(fmt::color::cyan), "\t*** Modifications to existing assemblies now allowed.\n");
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "list")) {
      handle_list(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "graph")) {
      handle_graph(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "delete")) {
      handle_delete(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "assembly")) {
      changed |= handle_assembly(tokens, region, allow_modify);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "attribute")) {
      changed |= handle_attribute(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "geometry")) {
      changed |= handle_geometry(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "time")) {
      changed |= handle_time(tokens, region);
    }
    else if (Ioss::Utils::substr_equal(tokens[0], "rename")) {
      changed |= handle_rename(tokens, region);
    }
    else {
      fmt::print(stderr, fg(fmt::color::yellow), "\tWARNING: Unrecognized command: {}\n",
                 tokens[0]);
    }
  }

  if (changed) {
    update_assembly_info(region, interFace);
    fmt::print(fg(fmt::color::cyan),
               "\n\t*** Database assembly structure modified. File update required.\n");
  }
  fmt::print("\n{} execution successful.\n", codename);
  return EXIT_SUCCESS;
}

namespace {
  void info_entity(const Ioss::StructuredBlock *sb, bool show_property)
  {
    fmt::print("\n{} {} [{}, Offset = {}] ", name(sb), fmt::join(sb->get_ijk_global(), "x"),
               fmt::join(sb->get_ijk_local(), "x"), fmt::join(sb->get_ijk_offset(), ", "));

    int64_t num_cell = sb->get_property("cell_count").get_int();
    int64_t num_node = sb->get_property("node_count").get_int();
    fmt::print("{:14} cells, {:14} nodes\n", fmt::group_digits(num_cell),
               fmt::group_digits(num_node));
    if (show_property) {
      Ioss::Utils::info_property(sb, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::Region &region, bool show_property)
  {
    fmt::print("\nRegion (global)\n");
    if (show_property) {
      Ioss::Utils::info_property(&region, Ioss::Property::ATTRIBUTE,
                                 "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::Assembly *as, bool show_property)
  {
    std::string modifier;
    if (as->property_exists("created")) {
      modifier = " [created]";
    }
    else if (as->property_exists("modified")) {
      modifier = " [modified]";
    }

    fmt::print("\n{} id: {:6d}, contains: {} member(s) of type {:>10s}.{}\n\tMembers: ", name(as),
               id(as), as->member_count(), as->contains_string(), modifier);
    for (const auto &mem : as->get_members()) {
      fmt::print("'{}' ", mem->name());
    }
    fmt::print("\n");
    if (show_property) {
      Ioss::Utils::info_property(as, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::Blob *blob, bool show_property)
  {
    fmt::print("\n{} id: {:6d}, contains: {} item(s).\n", name(blob), id(blob),
               blob->entity_count());
    if (show_property) {
      Ioss::Utils::info_property(blob, Ioss::Property::ATTRIBUTE,
                                 "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::ElementBlock *eb, bool show_property)
  {
    int64_t num_elem = eb->entity_count();

    std::string type       = eb->topology()->name();
    int64_t     num_attrib = eb->get_property("attribute_count").get_int();
    fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14} elements, {:3d} attributes.\n", name(eb),
               id(eb), type, fmt::group_digits(num_elem), num_attrib);
    if (show_property) {
      Ioss::Utils::info_property(eb, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::SideSet *ss, bool show_property)
  {
    fmt::print("\n{} id: {:6d}", name(ss), id(ss));
    if (ss->property_exists("bc_type")) {
#if defined(SEACAS_HAVE_CGNS)
      auto bc_type = ss->get_property("bc_type").get_int();
      fmt::print(", boundary condition type: {} ({})\n", BCTypeName[bc_type], bc_type);
#else
      fmt::print(", boundary condition type: {}\n", ss->get_property("bc_type").get_int());
#endif
    }
    const Ioss::SideBlockContainer &fbs = ss->get_side_blocks();
    for (auto &fb : fbs) {
      int64_t count      = fb->entity_count();
      int64_t num_attrib = fb->get_property("attribute_count").get_int();
      int64_t num_dist   = fb->get_property("distribution_factor_count").get_int();
      fmt::print("\t{}, {:8} sides, {:3d} attributes, {:8} distribution factors.\n", name(fb),
                 fmt::group_digits(count), num_attrib, fmt::group_digits(num_dist));
    }
    if (show_property) {
      Ioss::Utils::info_property(ss, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::NodeSet *ns, bool show_property)
  {
    int64_t count      = ns->entity_count();
    int64_t num_attrib = ns->get_property("attribute_count").get_int();
    int64_t num_dist   = ns->get_property("distribution_factor_count").get_int();
    fmt::print("\n{} id: {:6d}, {:8} nodes, {:3d} attributes, {:8} distribution factors.\n",
               name(ns), id(ns), fmt::group_digits(count), num_attrib, fmt::group_digits(num_dist));
    if (show_property) {
      Ioss::Utils::info_property(ns, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::NodeBlock *nb, bool show_property)
  {
    int64_t num_nodes  = nb->entity_count();
    int64_t num_attrib = nb->get_property("attribute_count").get_int();
    fmt::print("\n{} {:14} nodes, {:3d} attributes.\n", name(nb), fmt::group_digits(num_nodes),
               num_attrib);
    if (show_property) {
      Ioss::Utils::info_property(nb, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_time(const Ioss::Region &region)
  {
    fmt::print("\n");
    size_t ts_count = region.get_optional_property("state_count", 0);
    auto   width    = Ioss::Utils::number_width(ts_count);

    for (size_t step = 1; step <= ts_count; step++) {
      fmt::print("\tStep {:{}}, Time = {}\n", step, width, region.get_state_time(step));
    }
  }

  void set_db_properties(const Modify::Interface & /* interFace */, Ioss::DatabaseIO *dbi)
  {
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }
  }

  void handle_help(const std::string &topic)
  {
    bool all = Ioss::Utils::substr_equal(topic, "help");
    if (all) {
      fmt::print(fmt::emphasis::bold, "\n\tHELP");
      fmt::print(" [list | assembly | graph | attribute | regex | glob]\n");
      fmt::print(fmt::emphasis::bold, "\n\tEND | EXIT\n");
      fmt::print("\t\tEnd command input and output changed assembly definitions (if any).\n");
      fmt::print(fmt::emphasis::bold, "\n\tQUIT\n");
      fmt::print("\t\tEnd command input and exit with no rewriting of database. Some changes may "
                 "have already been made.\n");

      fmt::print(fmt::emphasis::bold, "\n\tALLOW MODIFICATIONS\n");
      fmt::print("\t\tBy default, io_modify will only allow creation of new assemblies.\n"
                 "\t\tIf this command is specified, then can modify assemblies that already exist "
                 "in database.\n"
                 "\t\tThis will cause the database to be rewritten. Without this option, it is "
                 "updated in place.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "list")) {
      fmt::print(fmt::emphasis::bold, "\n\tLIST ");
      fmt::print(
          "elementblock|block|structuredblock|assembly|nodeset|sideset|blob|times|summary\n");
      fmt::print(fmt::emphasis::bold, "\tLIST ");
      fmt::print("elementblock|block|structuredblock|assembly|nodeset|sideset|blob "
                 "{{names...}}\n");
      fmt::print(fmt::emphasis::bold, "\tLIST ");
      fmt::print("elementblock|block|structuredblock|assembly|nodeset|sideset|blob "
                 "MATCHES {{regex}}\n");
      fmt::print(fmt::emphasis::bold, "\tLIST ");
      fmt::print("elementblock|block|structuredblock|assembly|nodeset|sideset|blob "
                 "GLOB {{glob}}\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "assembly")) {
      fmt::print("\n\tFor all commands, if an assembly named `name` does not exist, it will be "
                 "created.\n");
      fmt::print(fmt::emphasis::bold, "\tASSEMBLY ");
      fmt::print("{{name}}\n");
      fmt::print("\t\tCreates an empty assembly named `name` if it does not exist.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} ADD {{name1}} {{name2}} ... {{nameL}}\n");
      fmt::print("\t\tAdds the specified entities to the assembly.  All entities must be the same "
                 "type.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} REMOVE {{name1}} {{name2}} ... {{nameL}}\n");
      fmt::print("\t\tRemoves the specified entities from the assembly.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} TYPE {{type}} MATCHES {{regex}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose name matches the {{regex}} will be added.\n");
      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} TYPE {{type}} GLOB {{glob}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose name matches the {{glob}} will be added.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} TYPE {{type}} NAMED {{list of one or more names}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose names are listed will be added.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} TYPE {{type}} RANGE {{id}} TO {{id}} BY {{step}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose id matches the specified range will be added.\n"
                 "\t\tNo message will be output for ids not matching an entity.\n");

      fmt::print(fmt::emphasis::bold, "\n\tASSEMBLY ");
      fmt::print("{{name}} TYPE {{type}} IDS {{id}}, {{id2}}, ..., {{idL}}\n");
      fmt::print(
          "\t\tAdds the entities of the specified type to the assembly.\n"
          "\t\tAll entities whose id matches an id in the list will be added.\n"
          "\t\tA warning message will be output if there is no entity with the requested id.\n");

      fmt::print(fmt::emphasis::bold, "\n\tDELETE ");
      fmt::print("{{name}}\n");
      fmt::print("\t\tRemove the assembly with the specified name.\n"
                 "\t\tCurrently only supported for assemblies created during this execution; not "
                 "for assemblies\n"
                 "\t\texisting on the input database.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "graph")) {
      fmt::print(fmt::emphasis::bold, "\n\tGRAPH OUTPUT ");
      fmt::print("[filename]\n");
      fmt::print(
          "\t\tCreate a 'dot' input file with the structure of the assembly graph.\n"
          "\t\tFile is named 'filename' or defaults to 'assembly.dot' if filename not given.\n");
      fmt::print(fmt::emphasis::bold, "\tGRAPH CHECK\n");
      fmt::print("\t\tCheck validity of assembly graph--are there any cycles.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "attribute")) {
      fmt::print(fmt::emphasis::bold, "\n\tATTRIBUTE ");
      fmt::print("{{ent_name}} ADD {{att_name}} STRING {{values...}}\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE ");
      fmt::print("{{ent_name}} ADD {{att_name}} DOUBLE {{values...}}\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE ");
      fmt::print("{{ent_name}} NAME {{att_name}} INTEGER {{values...}}\n");
      fmt::print("\t\tAdd an attribute to the specified entity 'ent_name'.\n"
                 "\t\tThe attribute will be named 'att_name' with value(s) 'values...'\n"
                 "\t\tCan also modify the values of an existing attribute.'\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE LIST ");
      fmt::print("{{ent_name...}}\n"
                 "\t\tList attributes for the selected entities\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE ");
      fmt::print("{{ent_type}} LIST\n"
                 "\t\tList attributes for all entities in the specified entity type\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE ");
      fmt::print("{{ent_type}} MATCH {{regex}}\n"
                 "\t\tList attributes for all entities in the specified entity type whose name "
                 "matches the regex.\n");
      fmt::print(fmt::emphasis::bold, "\tATTRIBUTE ");
      fmt::print("{{ent_type}} GLOB {{glob}}\n"
                 "\t\tList attributes for all entities in the specified entity type whose name "
                 "matches the glob.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "rename")) {
      fmt::print(fmt::emphasis::bold, "\n\tRENAME ");
      fmt::print("{{name}} TO {{new_name}}\n");
      fmt::print(fmt::emphasis::bold, "\tRENAME ");
      fmt::print("{{type}} {{id}} TO {{new_name}}\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "geometry")) {
      fmt::print(fmt::emphasis::bold, "\n\tGEOMETRY ROTATE ");
      fmt::print("{{X|Y|Z}} {{angle}} ...\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY MIRROR ");
      fmt::print("{{X|Y|Z}} ...\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY SCALE  ");
      fmt::print("{{X|Y|Z}} {{scale_factor}} ...\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY OFFSET ");
      fmt::print("{{X|Y|Z}} {{offset}} ...\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY ROTATE ");
      fmt::print("{{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{angle}}\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY MIRROR ");
      fmt::print("{{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} ...\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY SCALE  ");
      fmt::print("{{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{scale_factor}} ... "
                 "\n");
      fmt::print(fmt::emphasis::bold, "\tGEOMETRY OFFSET ");
      fmt::print("{{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{offset}} ...\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "time")) {
      fmt::print(fmt::emphasis::bold, "\n\tTIME SCALE  ");
      fmt::print("{{scale}}   (T_out = T_in * {{scale}})\n");
      fmt::print(fmt::emphasis::bold, "\tTIME OFFSET ");
      fmt::print("{{offset}}  (T_out = T_in + {{offset}})\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "regex")) {
      fmt::print("\n\tRegular Expression help (used in ASSEMBLY MATCHES and LIST MATCHES and "
                 "ATTRIBUTE LIST MATCHES options)\n"
                 "\t\tSupports \"POSIX Extended Regular Expressions\" syntax\n"
                 "\t\tSee https://www.regular-expressions.info/posix.html\n"
                 "\t\tQuickStart: https://www.regular-expressions.info/quickstart.html\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "glob")) {
      fmt::print(
          "\n\tGlob help (used in ASSEMBLY GLOB and LIST GLOB and ATTRIBUTE LIST GLOB options)\n"
          "\t\t?(pattern-list)   Matches zero or one occurrence of the given patterns\n"
          "\t\t*(pattern-list)   Matches zero or more occurrences of the given patterns\n"
          "\t\t+(pattern-list)   Matches one or more occurrences of the given patterns\n"
          "\t\t@(pattern-list)   Matches one of the given patterns\n"
          "\t\t!(pattern-list)   Matches anything except one of the given patterns\n"
          "\tGlob Examples\n"
          "\t\tblock*    : All names that start with 'block'\n"
          "\t\t[A-Z]*    : All names that start with a capital letter\n");
    }
  }

  void handle_list(const Ioss::NameList &tokens, const Ioss::Region &region, bool show_attribute)
  {
    if (tokens.size() > 1) {
      if (Ioss::Utils::substr_equal(tokens[1], "summary")) {
        region.output_summary(std::cout);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "region")) {
        info_entity(region, show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "elementblock")) {
        const auto &entities = region.get_element_blocks();
        info_entities(entities, tokens, region, "Element Blocks", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "assembly") ||
               Ioss::Utils::substr_equal(tokens[1], "assemblies")) {
        const auto &entities = region.get_assemblies();
        info_entities(entities, tokens, region, "Assemblies", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "nodeset") ||
               Ioss::Utils::substr_equal(tokens[1], "nset")) {
        const auto &entities = region.get_nodesets();
        info_entities(entities, tokens, region, "NodeSets", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "nodeblock")) {
        const auto &entities = region.get_node_blocks();
        info_entities(entities, tokens, region, "Node Block", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "structuredblock")) {
        const auto &entities = region.get_structured_blocks();
        info_entities(entities, tokens, region, "Structured Blocks", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "blocks")) {
        const auto type = region.mesh_type();
        if (type == Ioss::MeshType::UNSTRUCTURED) {
          const auto &entities = region.get_element_blocks();
          info_entities(entities, tokens, region, "Element Blocks", show_attribute);
        }
        else if (type == Ioss::MeshType::STRUCTURED) {
          const auto &entities = region.get_structured_blocks();
          info_entities(entities, tokens, region, "Structured Blocks", show_attribute);
        }
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "sideset") ||
               Ioss::Utils::substr_equal(tokens[1], "sset")) {
        const auto &entities = region.get_sidesets();
        info_entities(entities, tokens, region, "SideSets", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "blobs")) {
        const auto &entities = region.get_blobs();
        info_entities(entities, tokens, region, "Blobs", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "times")) {
        info_time(region);
      }
      else {
        fmt::print(stderr, fg(fmt::color::yellow), "\tWARNING: Unrecognized list option '{}'\n",
                   tokens[1]);
        handle_help("list");
      }
    }
    else {
      handle_help("list");
    }
  }

  class Graph
  {
    std::map<std::string, int>   m_vertices;
    Ioss::NameList               m_vertex;
    std::vector<Ioss::IntVector> m_adj; // Pointer to an array containing adjacency std::lists

    bool is_cyclic_internal(int v, std::vector<bool> &visited, std::vector<bool> &recStack);
    int  vertex(const std::string &node);
    int  size() const { return (int)m_vertices.size(); }

  public:
    void add_edge(const std::string &v, const std::string &w); // to add an edge to graph
    bool is_cyclic(); // returns true if there is a cycle in this graph
  };

  int Graph::vertex(const std::string &node)
  {
    if (m_vertices.find(node) == m_vertices.end()) {
      int index        = size();
      m_vertices[node] = index;
      m_vertex.push_back(node);
      //      fmt::print("Node {} is index {}\n", node, index);
    }
    assert(node == m_vertex[m_vertices[node]]);
    return m_vertices[node];
  }

  void Graph::add_edge(const std::string &v, const std::string &w)
  {
    int vv = vertex(v);
    int ww = vertex(w);
    if ((int)m_adj.size() <= vv) {
      m_adj.resize(vv + 1);
    }
    assert((int)m_adj.size() > vv);
    m_adj[vv].push_back(ww); // Add w to v’s std::list.
  }

  // This function is a variation of DFSUtil() in https://www.geeksforgeeks.org/archives/18212
  bool Graph::is_cyclic_internal(int v, std::vector<bool> &visited, std::vector<bool> &recStack)
  {
    if (visited[v] == false) {
      // Mark the current node as visited and part of recursion stack
      visited[v]  = true;
      recStack[v] = true;

      // Recur for all the m_vertices adjacent to this vertex
      if (v < (int)m_adj.size()) {
        for (auto i = m_adj[v].begin(); i != m_adj[v].end(); ++i) {
          if (!visited[*i] && is_cyclic_internal(*i, visited, recStack)) {
            if (*i != 0 && v != 0) {
              fmt::print(fg(fmt::color::yellow), "\t*** Cycle contains {} -> {}\n", m_vertex[v],
                         m_vertex[*i]);
            }
            return true;
          }
          else if (recStack[*i]) {
            if (*i != 0 && v != 0) {
              fmt::print(fg(fmt::color::yellow), "\t*** Cycle contains {} -> {}\n", m_vertex[v],
                         m_vertex[*i]);
            }
            return true;
          }
        }
      }
    }
    recStack[v] = false; // remove the vertex from recursion stack
    return false;
  }

  // Returns true if the graph contains a cycle, else false.
  // This function is a variation of DFS() in https://www.geeksforgeeks.org/archives/18212
  bool Graph::is_cyclic()
  {
    // Mark all the m_vertices as not visited and not part of recursion
    // stack
    std::vector<bool> visited(size());
    std::vector<bool> recStack(size());
    for (int i = 0; i < size(); i++) {
      visited[i]  = false;
      recStack[i] = false;
    }

    // Call the recursive helper function to detect cycle in different
    // DFS trees
    for (int i = 0; i < size(); i++)
      if (is_cyclic_internal(i, visited, recStack))
        return true;

    return false;
  }

  void handle_graph(const Ioss::NameList &tokens, const Ioss::Region &region)
  {
    if (tokens.size() == 1) {
      handle_help("graph");
      return;
    }

    if (Ioss::Utils::substr_equal(tokens[1], "output")) {
      std::string filename = "assembly.dot";
      if (tokens.size() > 2) {
        filename = tokens[2];
      }
      // Create a file containing a 'dot' representation of the assembly graph.
      std::ofstream dot(filename);
      fmt::print(dot, "digraph assemblies {{\n"
                      "node [shape=rectangle]; ");
      const auto &assemblies = region.get_assemblies();
      for (const auto *assembly : assemblies) {
        fmt::print(dot, "{} ", assembly->name());
      }
      fmt::print(dot, ";\nnode [shape=oval];\n");

      for (const auto *assembly : assemblies) {
        const auto &members = assembly->get_members();
        for (const auto *member : members) {
          fmt::print(dot, "{} -> {}\n", assembly->name(), member->name());
        }
      }
      fmt::print(dot, "}}\n");
      fmt::print(fg(fmt::color::cyan), "\t*** Assembly graph output to '{}'.\n", filename);
    }
    else if (Ioss::Utils::substr_equal(tokens[1], "check")) {
      // Check graph for cycles.
      Graph       assem;
      const auto &assemblies = region.get_assemblies();
      for (const auto *assembly : assemblies) {
        assem.add_edge(region.name(), assembly->name());
      }
      for (const auto *assembly : assemblies) {
        const auto &members = assembly->get_members();
        for (const auto *member : members) {
          assem.add_edge(assembly->name(), member->name());
        }
      }
      if (assem.is_cyclic()) {
        fmt::print(stderr, fg(fmt::color::yellow),
                   "\tWARNING: Assembly graph is invalid--contains at least one cycle.\n",
                   tokens[1]);
      }
      else {
        fmt::print(stderr, fg(fmt::color::cyan), "\t*** Assembly graph is valid--no cycles.\n",
                   tokens[1]);
      }
    }
    else {
      fmt::print(stderr, fg(fmt::color::yellow), "\tWARNING: Unrecognized graph option '{}'\n",
                 tokens[1]);
      handle_help("graph");
    }
  }

  bool handle_delete(const Ioss::NameList &tokens, Ioss::Region &region)
  {
    // Returns true if requested assembly was deleted.
    // False if assembly does not exist, or was not deletable (not created during this run)
    // NOTE: If the assembly is a member of other assemblies, it will be removed from them
    //       Since currently, the removed assembly must have been created during this run,
    //       then it also must have been added to the other assemblies also during this run,
    //       so the modified flag has already been set on them and we shouldn't have to
    //       do anything else here to track whether in modify mode or to mark assemblies
    //       as modified.  Once allow the deletion of existing assemblies, then will have
    //       to do better tracking of what other assemblies were modified as a result of
    //       removing/deleting this assembly.
    if (tokens.size() > 1) {
      Ioss::Assembly *assem = region.get_assembly(tokens[1]);
      if (assem == nullptr) {
        fmt::print(stderr, fg(fmt::color::yellow),
                   "WARNING: Requested Assembly '{}' does not exist.\n", tokens[1]);
        return false;
      }
      else {
        if (assem->property_exists("created")) {
          if (region.remove(assem)) {
            fmt::print(fg(fmt::color::cyan), "\t***Assembly '{}' deleted successfully.\n",
                       tokens[1]);
            return true;
          }
        }
        else {
          fmt::print(stderr, fg(fmt::color::red),
                     "ERROR: Requested Assembly '{}' was not created during this execution.  Not "
                     "deletable.\n",
                     tokens[1]);
          return false;
        }
      }
    }
    else {
      handle_help("delete");
    }
    return false;
  }

  bool handle_attribute(const Ioss::NameList &tokens, Ioss::Region &region)
  {
    //     0          1        2       3         4       5...
    // ATTRIBUTE {{ent_name}} ADD {{att_name}} STRING {{values...}}
    // ATTRIBUTE {{ent_name}} ADD {{att_name}} DOUBLE {{values...}}
    // ATTRIBUTE {{ent_name}} ADD {{att_name}} INTEGER {{values...}}
    // ATTRIBUTE LIST {{ent_name}} ...
    // ATTRIBUTE {{ent_type}} LIST
    // ATTRIBUTE {{ent_type}} MATCH {regex}

    // Get requested entity...
    if (Ioss::Utils::substr_equal(tokens[2], "add")) {
      // Must be at least 6 tokens...
      if (tokens.size() < 6) {
        fmt::print(stderr,
#if !defined __NVCC__
                   fg(fmt::color::red),
#endif
                   "ERROR: ATTRIBUTE Command does not have enough tokens to be valid.\n"
                   "\t\t{}\n",
                   fmt::join(tokens, " "));
        handle_help("attribute");
        return false;
      }

      Ioss::GroupingEntity *ge = nullptr;
      if (Ioss::Utils::substr_equal(tokens[1], "region")) {
        ge = &region;
      }
      else {
        ge = region.get_entity(tokens[1]);
      }
      if (ge == nullptr) {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: Entity '{}' not found.\n", tokens[1]);
        return false;
      }

      // Now get name of attribute/property to create...
      const std::string &att_name = tokens[3];

      // Now, the attribute type and whether vector or scalar...
      size_t value_count = tokens.size() - 5;
      if (Ioss::Utils::substr_equal(tokens[4], "string")) {
        std::string value = tokens[5];
        for (size_t i = 6; i < tokens.size(); i++) {
          value += " " + tokens[i];
        }
        ge->property_add(Ioss::Property(att_name, value, Ioss::Property::ATTRIBUTE));
      }
      else if (Ioss::Utils::substr_equal(tokens[4], "double")) {
        std::vector<double> values(value_count);
        for (size_t i = 0; i < value_count; i++) {
          try {
            double val = std::stod(tokens[i + 5]);
            values[i]  = val;
          }
          catch (const std::exception &x) {
            fmt::print(stderr, fg(fmt::color::red), "ERROR: Invalid double value entered: {}. {}\n",
                       tokens[i + 5], x.what());
            return false;
          }
        }

        ge->property_add(Ioss::Property(att_name, values, Ioss::Property::ATTRIBUTE));
      }
      else if (Ioss::Utils::substr_equal(tokens[4], "integer")) {
        Ioss::IntVector values(value_count);
        for (size_t i = 0; i < value_count; i++) {
          try {
            int val   = std::stoi(tokens[i + 5]);
            values[i] = val;
          }
          catch (const std::exception &x) {
            fmt::print(stderr, fg(fmt::color::red),
                       "ERROR: Invalid integer value entered: {}. {}\n", tokens[i + 5], x.what());
            return false;
          }
        }
        ge->property_add(Ioss::Property(att_name, values, Ioss::Property::ATTRIBUTE));
      }
      else {
        fmt::print(stderr, fg(fmt::color::red),
                   "ERROR: Unrecognized attribute type '{}' for attribute {}.\n", tokens[4],
                   att_name);
        handle_help("attribute");
        return false;
      }
      attributes_modified.push_back(ge);
      return true;
    }
    else if (Ioss::Utils::substr_equal(tokens[1], "list")) {
      // ATTRIBUTE LIST {{ent_name}} ...
      for (size_t i = 2; i < tokens.size(); i++) {
        Ioss::GroupingEntity *ge = nullptr;
        if (Ioss::Utils::substr_equal(tokens[i], "region")) {
          ge = &region;
        }
        else {
          ge = region.get_entity(tokens[i]);
        }
        if (ge != nullptr) {
          std::string prefix =
              fmt::format("\n{} id: {:6d}\n\tAttributes (Reduction): ", name(ge), id(ge));
          Ioss::Utils::info_property(ge, Ioss::Property::ATTRIBUTE, prefix, "\t");
        }
      }
      return false;
    }
    else if (Ioss::Utils::substr_equal(tokens[2], "list")) {
      // ATTRIBUTE {{ent_type}} LIST
      handle_list(tokens, region, true);
      return false;
    }
    else if (Ioss::Utils::substr_equal(tokens[2], "matches")) {
      // ATTRIBUTE {{ent_type}} MATCH {regex}
      handle_list(tokens, region, true);
      return false;
    }
    else if (Ioss::Utils::substr_equal(tokens[2], "glob")) {
      // ATTRIBUTE {{ent_type}} MATCH {regex}
      handle_list(tokens, region, true);
      return false;
    }
    fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized attribute command.\n");
    handle_help("attribute");
    return false;
  }

  bool handle_rename(const Ioss::NameList &tokens, Ioss::Region &region)
  {
    //     0          1        2       3         4
    // RENAME   {{ent_name}}   TO   {{new_name}}
    // RENAME   {{ent_type}} {{id}}    TO    {{new_name}}

    // Must be at least 4 tokens...
    if (tokens.size() < 4) {
      fmt::print(stderr,
#if !defined __NVCC__
                 fg(fmt::color::red),
#endif
                 "ERROR: RENAME Command does not have enough tokens to be valid.\n"
                 "\t\t{}\n",
                 fmt::join(tokens, " "));
      handle_help("rename");
      return false;
    }

    // See if asking for actual entity by name or by type + id
    Ioss::GroupingEntity *ge       = nullptr;
    const std::string    &new_name = tokens[tokens.size() - 1];

    if (tokens.size() == 5 && Ioss::Utils::str_equal(tokens[3], "to")) {
      // Type + ID
      auto entity_type = get_entity_type(tokens[1]);
      if (entity_type == Ioss::INVALID_TYPE) {
        fmt::print(stderr, fg(fmt::color::yellow), "WARNING: Unrecognized entity type '{}'.\n",
                   tokens[1]);
        return false;
      }

      int entity_id = std::stoi(tokens[2]);
      ge            = region.get_entity(entity_id, entity_type);
      if (ge == nullptr) {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: Entity type '{}' with id {} not found.\n",
                   tokens[1], tokens[2]);
        return false;
      }
    }
    else if (tokens.size() == 4 && Ioss::Utils::str_equal(tokens[2], "to")) {
      ge = region.get_entity(tokens[1]);
      if (ge == nullptr) {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: Entity '{}' not found.\n", tokens[1]);
        return false;
      }
    }
    else {
      fmt::print(stderr,
#if !defined __NVCC__
                 fg(fmt::color::yellow),
#endif
                 "\tWARNING: Unrecognized rename syntax '{}'\n", fmt::join(tokens, " "));
      handle_help("rename");
    }

    if (ge != nullptr) {
      int            exoid     = region.get_database()->get_file_pointer();
      auto           ioss_type = ge->type();
      ex_entity_type exo_type  = Ioex::map_exodus_type(ioss_type);
      auto           ierr      = ex_put_name(exoid, exo_type, id(ge), new_name.c_str());
      if (ierr != EX_NOERR) {
        Ioex::exodus_error(exoid, __LINE__, __func__, __FILE__);
      }
      ge->set_name(new_name);
    }
    return false;
  }

  void build_block_list(Ioss::Region &region, const Ioss::GroupingEntity *ge,
                        std::vector<const Ioss::GroupingEntity *> &blocks)
  {
    if (ge) {
      if (ge->type() == Ioss::ELEMENTBLOCK) {
        auto *eb = dynamic_cast<const Ioss::ElementBlock *>(ge);
        if (eb != nullptr) {
          blocks.push_back(eb);
        }
        return;
      }
      else if (ge->type() == Ioss::STRUCTUREDBLOCK) {
        auto *sb = dynamic_cast<const Ioss::StructuredBlock *>(ge);
        if (sb != nullptr) {
          blocks.push_back(sb);
        }
        return;
      }
      else if (ge->type() == Ioss::ASSEMBLY) {
        auto       *as      = dynamic_cast<const Ioss::Assembly *>(ge);
        const auto &members = as->get_members();
        for (const auto *member : members) {
          build_block_list(region, member, blocks);
        }
        return;
      }
    }
  }

  std::vector<int> get_filtered_node_list(Ioss::Region                                    &region,
                                          const std::vector<const Ioss::GroupingEntity *> &blocks)
  {
    const auto type = region.get_database()->get_format();
    if (type != "Exodus") {
      return std::vector<int>();
    }
    auto node_count = region.get_property("node_count").get_int();
    if (blocks.empty() ||
        blocks.size() == (size_t)region.get_property("element_block_count").get_int()) {
      return Ioss::IntVector(node_count, 1);
    }
    else {
      Ioss::IntVector node_filter(node_count);
      // Iterate all element blocks in 'blocks', get connectivity_raw
      // and set `node_filter` to 1 for nodes in connectivity list.
      if (region.get_database()->int_byte_size_api() == 4) {
        Ioss::IntVector connect;
        for (const auto *block : blocks) {
          block->get_field_data("connectivity_raw", connect);
          for (auto node : connect) {
            node_filter[node - 1] = 1;
          }
        }
      }
      else {
        std::vector<int64_t> connect;
        for (const auto *block : blocks) {
          block->get_field_data("connectivity_raw", connect);
          for (auto node : connect) {
            node_filter[node - 1] = 1;
          }
        }
      }
      return node_filter;
    }
  }

  bool handle_time(const Ioss::NameList &tokens, Ioss::Region &region)
  {
    //   0      1       2
    // TIME   SCALE  {{scale}}
    // TIME   OFFSET {{offset}
    if (tokens.size() < 3) {
      fmt::print(stderr,
#if !defined __NVCC__
                 fg(fmt::color::red),
#endif
                 "ERROR: TIME Command does not have enough tokens to be valid.\n"
                 "\t\t{}\n",
                 fmt::join(tokens, " "));
      handle_help("time");
      return false;
    }
    if (Ioss::Utils::substr_equal(tokens[1], "scale")) {
      double scale = std::stod(tokens[2]);
      modify_time(region, scale, 0.0);
      fmt::print(fg(fmt::color::cyan), "\t*** Database time scaled:  T_out = T_in * {}.\n", scale);
      return false;
    }
    if (Ioss::Utils::substr_equal(tokens[1], "offset")) {
      double offset = std::stod(tokens[2]);
      modify_time(region, 1.0, offset);
      fmt::print(fg(fmt::color::cyan), "\t*** Database time offset:  T_out = T_in + {}.\n", offset);
      return false;
    }
    fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized time command.\n");
    handle_help("time");
    return false;
  }

  bool handle_geometry(const Ioss::NameList &tokens, Ioss::Region &region)
  {
    //     0        1        2         3         4       5...
    // GEOMETRY   ROTATE {{X|Y|Z}} {{angle}} ...
    // GEOMETRY   MIRROR {{X|Y|Z}} ...
    // GEOMETRY   SCALE  {{X|Y|Z}} {{scale}} ...
    // GEOMETRY   OFFSET {{X|Y|Z}} {{offset}} ...
    // GEOMETRY   ROTATE {{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{angle}}  ...
    // GEOMETRY   SCALE  {{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{scale}}  ...
    // GEOMETRY   OFFSET {{ELEMENTBLOCKS|BLOCKS|ASSEMBLY}} {{names}} {{X|Y|Z}} {{offset}} ...

    if (tokens.size() < 3) {
      fmt::print(stderr,
#if !defined __NVCC__
                 fg(fmt::color::red),
#endif
                 "ERROR: GEOMETRY Command does not have enough tokens to be valid.\n"
                 "\t\t{}\n",
                 fmt::join(tokens, " "));
      handle_help("geometry");
      return false;
    }

    const auto type = region.mesh_type();

    // See if applying just selected (and connected) blocks.
    size_t                                    idx = 2;
    std::vector<const Ioss::GroupingEntity *> blocks;
    if (Ioss::Utils::substr_equal(tokens[idx], "elementblocks") ||
        Ioss::Utils::substr_equal(tokens[idx], "structuredblocks") ||
        Ioss::Utils::substr_equal(tokens[idx], "blocks") ||
        Ioss::Utils::substr_equal(tokens[idx], "assembly")) {
      // Parse list of block|assembly names...
      idx++;
      while (!(Ioss::Utils::str_equal(tokens[idx], "x") ||
               Ioss::Utils::str_equal(tokens[idx], "y") ||
               Ioss::Utils::str_equal(tokens[idx], "z"))) {
        const auto &name = tokens[idx++];
        auto       *ge =
            region.get_entity(name, type == Ioss::MeshType::UNSTRUCTURED ? Ioss::ELEMENTBLOCK
                                                                         : Ioss::STRUCTUREDBLOCK);
        if (ge == nullptr) {
          ge = region.get_entity(name, Ioss::ASSEMBLY);
        }
        build_block_list(region, ge, blocks);
      }

      // Now filter out any duplicates that might be in the list...
      Ioss::Utils::uniquify(blocks);
    }

    // If blocks is non-empty, then we are applying geometry modification to a subset of the model.
    // In that case, we need to get all blocks that are connected to the user-specified blocks...
    if (!blocks.empty()) {
      if (type == Ioss::MeshType::UNSTRUCTURED) {
        auto tmp(blocks);
        for (const auto *block : blocks) {
          auto *eb = dynamic_cast<const Ioss::ElementBlock *>(block);
          if (eb != nullptr) {
            const auto &connected = eb->get_block_adjacencies();
            for (const auto &connect : connected) {
              auto *elb = region.get_element_block(connect);
              tmp.push_back(elb);
            }
          }
        }
        Ioss::Utils::uniquify(tmp);
        blocks = tmp;
      }
    }

    if (blocks.empty()) {
      fmt::print(fg(fmt::color::cyan), "\n\t*** {} transformation will be applied to ALL blocks.\n",
                 tokens[1]);
    }
    else {
      fmt::print(fg(fmt::color::cyan), "\n\t*** {} transformation will be applied to blocks:\n\t",
                 tokens[1]);
      for (const auto *block : blocks) {
        fmt::print(fg(fmt::color::cyan), "{}, ", block->name());
      }
      fmt::print("\n");
    }

    if (Ioss::Utils::substr_equal(tokens[1], "rotate")) {
      real rotation_matrix[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

      // Get rotation axis...
      do {
        const std::string &axis = tokens[idx++];
        if (axis != "x" && axis != "y" && axis != "z" && axis != "X" && axis != "Y" &&
            axis != "Z") {
          fmt::print(stderr, fg(fmt::color::red),
                     "ERROR: Unrecognized axis {}.  Expected one of xyzXYZ.\n", axis);
          handle_help("geometry");
          return false;
        }
        double angle = std::stod(tokens[idx++]);
        auto   ok    = update_rotation_matrix(rotation_matrix, axis, angle);
        if (!ok) {
          return false;
        }
      } while (idx < tokens.size());

      // Do the rotation...
      rotate_filtered_coordinates(region, rotation_matrix, blocks);
      fmt::print(fg(fmt::color::cyan), "\t*** Database coordinates rotated.\n");
      return false;
    }

    if (Ioss::Utils::substr_equal(tokens[1], "scale")) {
      real scale[3] = {1.0, 1.0, 1.0};

      // Get scale axis and scale factor...
      do {
        const std::string &axis = tokens[idx++];
        if (axis != "x" && axis != "y" && axis != "z" && axis != "X" && axis != "Y" &&
            axis != "Z") {
          fmt::print(stderr, fg(fmt::color::red),
                     "ERROR: Unrecognized axis {}.  Expected one of xyzXYZ.\n", axis);
          handle_help("geometry");
          return false;
        }
        double factor = std::stod(tokens[idx++]);
        if (Ioss::Utils::substr_equal(axis, "x")) {
          scale[0] = factor;
        }
        else if (Ioss::Utils::substr_equal(axis, "y")) {
          scale[1] = factor;
        }
        else if (Ioss::Utils::substr_equal(axis, "z")) {
          scale[2] = factor;
        }
      } while (idx < tokens.size());

      // Do the transformation...
      scale_filtered_coordinates(region, scale, blocks);
      fmt::print(fg(fmt::color::cyan), "\t*** Database coordinate(s) scaled.\n");
      return false;
    }

    if (Ioss::Utils::substr_equal(tokens[1], "mirror")) {
      real scale[3] = {1.0, 1.0, 1.0};

      // Get mirror axis ...
      do {
        const std::string &axis = tokens[idx++];
        if (axis != "x" && axis != "y" && axis != "z" && axis != "X" && axis != "Y" &&
            axis != "Z") {
          fmt::print(stderr, fg(fmt::color::red),
                     "ERROR: Unrecognized axis {}.  Expected one of xyzXYZ.\n", axis);
          handle_help("geometry");
          return false;
        }
        if (Ioss::Utils::substr_equal(axis, "x")) {
          scale[0] = -1.0;
        }
        else if (Ioss::Utils::substr_equal(axis, "y")) {
          scale[1] = -1.0;
        }
        else if (Ioss::Utils::substr_equal(axis, "z")) {
          scale[2] = -1.0;
        }
      } while (idx < tokens.size());

      // Do the transformation...
      scale_filtered_coordinates(region, scale, blocks);
      fmt::print(fg(fmt::color::cyan), "\t*** Database coordinate(s) mirrored.\n");
      return false;
    }

    if (Ioss::Utils::substr_equal(tokens[1], "offset")) {
      real offset[3] = {1.0, 1.0, 1.0};

      // Get offset axis and offset factor...
      do {
        const std::string &axis = tokens[idx++];
        if (axis != "x" && axis != "y" && axis != "z" && axis != "X" && axis != "Y" &&
            axis != "Z") {
          fmt::print(stderr, fg(fmt::color::red),
                     "ERROR: Unrecognized axis {}.  Expected one of xyzXYZ.\n", axis);
          handle_help("geometry");
          return false;
        }
        double factor = std::stod(tokens[idx++]);
        if (Ioss::Utils::substr_equal(axis, "x")) {
          offset[0] = factor;
        }
        else if (Ioss::Utils::substr_equal(axis, "y")) {
          offset[1] = factor;
        }
        else if (Ioss::Utils::substr_equal(axis, "z")) {
          offset[2] = factor;
        }
      } while (idx < tokens.size());

      // Do the transformation...
      offset_filtered_coordinates(region, offset, blocks);
      fmt::print(fg(fmt::color::cyan), "\t*** Database coordinate(s) offset.\n");
      return false;
    }

    fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized geometry command.\n");
    handle_help("geometry");
    return false;
  }

  bool handle_assembly(const Ioss::NameList &tokens, Ioss::Region &region, bool allow_modify)
  {
    bool            changed = false;
    Ioss::Assembly *assem   = nullptr;

    if (tokens.size() > 1) {
      assem = region.get_assembly(tokens[1]);
      if (assem == nullptr) {
        // New assembly...
        assem      = new Ioss::Assembly(region.get_database(), tokens[1]);
        auto my_id = get_next_assembly_id(region);
        assem->property_add(Ioss::Property("id", my_id));
        assem->property_add(Ioss::Property("created", 1));
        region.add(assem);
        fmt::print(fg(fmt::color::cyan), "\t*** Created Assembly '{}' with id {}.\n", tokens[1],
                   my_id);
        // Don't set changed to true; only set if members modified
      }
    }
    else {
      handle_help("assembly");
      return false;
    }

    if (assem == nullptr) {
      fmt::print(stderr, fg(fmt::color::red), "ERROR: Unable to create or access assembly '{}'.\n",
                 tokens[1]);
      return false;
    }

    bool created = false;
    if (assem->property_exists("created")) {
      created = assem->get_property("created").get_int() == 1;
    }
    if (!allow_modify && !created) {
      fmt::print(stderr, fg(fmt::color::red),
                 "ERROR: Unable to modify an existing assembly '{}'.\n\tRestart with "
                 "`--allow_modifications` option or enter 'ALLOW MODIFICATIONS' command\n",
                 tokens[1]);
      return false;
    }

    if (tokens.size() > 2) {
      try {
        if (Ioss::Utils::substr_equal(tokens[2], "add")) {
          // List of names ...
          for (size_t i = 3; i < tokens.size(); i++) {
            auto *member = region.get_entity(tokens[i]);
            if (member != nullptr) {
              if (assem->add(member)) {
                changed = true;
              }
            }
            else {
              fmt::print(fg(fmt::color::yellow), "\t*** Entity '{}' does not exist. Not added.\n",
                         tokens[i]);
            }
          }
        }
        else if (Ioss::Utils::substr_equal(tokens[2], "remove")) {
          // List of names ...
          for (size_t i = 3; i < tokens.size(); i++) {
            auto *member = region.get_entity(tokens[i]);
            if (member != nullptr) {
              if (assem->remove(member)) {
                changed = true;
              }
            }
          }
        }
        else if (Ioss::Utils::substr_equal(tokens[2], "type")) {
          // Determine type of add...
          Ioss::EntityType type = get_entity_type(tokens[3]);
          if (type == Ioss::INVALID_TYPE) {
            fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized entity type: '{}'\n",
                       tokens[3]);
            return changed;
          }
          if (Ioss::Utils::substr_equal(tokens[4], "matches")) {
            // regex match on names
            // Get list of all names for this entity type...
            Ioss::NameList names = get_name_list(region, type);

            std::regex reg(tokens[5], std::regex::extended);

            // Check for match against all names in list...
            for (const auto &name : names) {
              if (std::regex_match(name, reg)) {
                const auto *entity = region.get_entity(name, type);
                if (entity != nullptr) {
                  if (assem->add(entity)) {
                    changed = true;
                  }
                }
              }
            }
          }
          else if (Ioss::Utils::substr_equal(tokens[4], "glob")) {
            // regex match on names
            // Get list of all names for this entity type...
            Ioss::NameList names = get_name_list(region, type);

            Ioss::glob::glob glob(tokens[5]);

            // Check for match against all names in list...
            for (const auto &name : names) {
              if (Ioss::glob::glob_match(name, glob)) {
                const auto *entity = region.get_entity(name, type);
                if (entity != nullptr) {
                  if (assem->add(entity)) {
                    changed = true;
                  }
                }
              }
            }
          }
          else if (Ioss::Utils::substr_equal(tokens[4], "named")) {
            // list of names
            for (size_t i = 5; i < tokens.size(); i++) {
              const auto *entity = region.get_entity(tokens[i], type);
              if (entity != nullptr) {
                if (assem->add(entity)) {
                  changed = true;
                }
              }
            }
          }
          else if (Ioss::Utils::substr_equal(tokens[4], "ids")) {
            // list of ids
            for (size_t i = 5; i < tokens.size(); i++) {
              size_t      id     = std::stoul(tokens[i]);
              const auto *entity = region.get_entity(id, type);
              if (entity != nullptr) {
                if (assem->add(entity)) {
                  changed = true;
                }
              }
            }
          }
          else if (Ioss::Utils::substr_equal(tokens[4], "range")) {
            //     0        1     2     3      4      5    6    7    8     9
            // ASSEMBLY {{name}} TYPE {{type}} RANGE {{id}} TO {{id}} BY {{step}}
            size_t begin = std::stoul(tokens[5]);
            size_t end   = begin;
            size_t step  = 1;
            if (tokens.size() >= 8 && Ioss::Utils::substr_equal(tokens[6], "to")) {
              end = std::stoul(tokens[7]);
            }
            if (tokens.size() >= 10 && Ioss::Utils::substr_equal(tokens[8], "by")) {
              step = std::stoul(tokens[9]);
            }
            for (size_t id = begin; id <= end; id += step) {
              const auto *entity = region.get_entity(id, type);
              if (entity != nullptr) {
                if (assem->add(entity)) {
                  changed = true;
                }
              }
            }
          }
        }
        else {
          fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized assembly option '{}'.\n",
                     tokens[2]);
          return changed;
        }
      }
      catch (const std::exception &x) {
        fmt::print(stderr, fg(fmt::color::red), "{}\n", x.what());
      }
    }

    if (changed) {
      assem->property_add(Ioss::Property("modified", true));
    }
    else {
      fmt::print(stderr, fg(fmt::color::yellow), "WARNING: Command did not modify assembly '{}'\n",
                 assem->name());
    }

    return changed;
  }

  template <typename T> Ioss::NameList get_entity_names(const std::vector<T *> &entity_list)
  {
    Ioss::NameList names;
    names.reserve(entity_list.size());

    for (const auto *entity : entity_list) {
      names.push_back(entity->name());
    }
    return names;
  }

  Ioss::NameList get_name_list(const Ioss::Region &region, Ioss::EntityType type)
  {
    Ioss::NameList names;
    switch (type) {
    case Ioss::ELEMENTBLOCK: {
      const auto &entities = region.get_element_blocks();
      names                = get_entity_names(entities);
    } break;
    case Ioss::NODESET: {
      const auto &entities = region.get_nodesets();
      names                = get_entity_names(entities);
    } break;
    case Ioss::SIDESET: {
      const auto &entities = region.get_sidesets();
      names                = get_entity_names(entities);
    } break;
    case Ioss::ASSEMBLY: {
      const auto &entities = region.get_assemblies();
      names                = get_entity_names(entities);
    } break;
    case Ioss::BLOB: {
      const auto &entities = region.get_blobs();
      names                = get_entity_names(entities);
    } break;
    default: break;
    }
    return names;
  }

#if defined(SEACAS_HAVE_CGNS)
  void update_cgns_assembly_info(Ioss::Region &region, const Modify::Interface & /* interFace */)
  {
    region.end_mode(Ioss::STATE_DEFINE_MODEL);
    int file_ptr = region.get_database()->get_file_pointer();

    fmt::print(fg(fmt::color::cyan), "\n\t*** Database changed. Updating assembly definitions.\n");
    const auto &assemblies = region.get_assemblies();
    for (const auto *assembly : assemblies) {
      if (assembly->property_exists("modified")) {
        if (assembly->property_exists("created")) {
          fmt::print(fg(fmt::color::cyan), "\t*** Creating assembly '{}'\n", assembly->name());
          Iocgns::Utils::output_assembly(file_ptr, assembly, false, true);
        }
        else {
          fmt::print(stderr, fg(fmt::color::yellow),
                     "WARNING: Can not modify existing assembly '{}'  yet.\n", assembly->name());
        }
      }
    }
  }
#endif

  void update_exodus_assembly_info(Ioss::Region &region, const Modify::Interface &interFace)
  {
    std::vector<Ioex::Assembly> ex_assemblies;
    bool                        modify_existing = false;

    region.end_mode(Ioss::STATE_DEFINE_MODEL);
    fmt::print(fg(fmt::color::cyan), "\n\t*** Database changed. Updating assembly definitions.\n");
    const auto &assemblies = region.get_assemblies();
    for (const auto *assembly : assemblies) {
      if (assembly->property_exists("modified")) {
        ex_assemblies.emplace_back(*assembly);
        if (!assembly->property_exists("created")) {
          fmt::print(fg(fmt::color::cyan), "\t*** Modifying assembly '{}'\n", assembly->name());
          modify_existing = true;
        }
        else {
          fmt::print(fg(fmt::color::cyan), "\t*** Creating assembly '{}'\n", assembly->name());
        }
      }
    }

    int exoid = region.get_database()->get_file_pointer();
    if (modify_existing) {
      // Need to create a temporary database to copy the database into.
      // Make sure has same int size as current.
      int                   byte_size = region.get_database()->int_byte_size_db();
      Ioss::PropertyManager properties;
      properties.add(Ioss::Property("INTEGER_SIZE_DB", byte_size));
      std::string       out_file  = interFace.filename() + ".mod";
      std::string       file_type = interFace.type();
      Ioss::DatabaseIO *dbo       = Ioss::IOFactory::create(
          file_type, out_file, Ioss::WRITE_RESTART, Ioss::ParallelUtils::comm_world(), properties);

      if (dbo == nullptr || !dbo->ok(true)) {
        std::exit(EXIT_FAILURE);
      }
      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region reg_out(dbo, "region_tmp");

      int out_exoid = reg_out.get_database()->get_file_pointer();
      Ioex::Internals::update_assembly_data(out_exoid, ex_assemblies, 1);
      Ioex::Internals::copy_database(exoid, out_exoid);
      Ioex::Internals::update_assembly_data(out_exoid, ex_assemblies, 2);
      Ioex::write_reduction_attributes(exoid, attributes_modified);

      // Now, remove old file and replace with new...
      region.get_database()->closeDatabase();
      reg_out.get_database()->closeDatabase();
      Ioss::FileInfo in_file(interFace.filename());
      in_file.remove_file();

      if (std::rename(out_file.c_str(), interFace.filename().c_str()) != 0) {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: Could not update modified file {} to {}.\n",
                   out_file, interFace.filename());
        return;
      }
    }
    else {
      Ioex::Internals::update_assembly_data(exoid, ex_assemblies);
      Ioex::write_reduction_attributes(exoid, attributes_modified);
    }
  }

  void update_assembly_info(Ioss::Region &region, const Modify::Interface &interFace)
  {
    // Determine type of underlying database...
    const auto type = region.get_database()->get_format();
    if (type == "Exodus") {
      update_exodus_assembly_info(region, interFace);
    }
    else if (type == "CGNS") {
#if defined(SEACAS_HAVE_CGNS)
      update_cgns_assembly_info(region, interFace);
#else
      fmt::print(stderr, fg(fmt::color::red), "ERROR: CGNS capability is not enabled.\n");
#endif
    }
    else {
      fmt::print(stderr, fg(fmt::color::red),
                 "ERROR: Can not modify the database '{}' of type '{}'.\n", interFace.filename(),
                 type);
    }
  }

  void rotate_filtered_coordinates(const Ioss::GroupingEntity *nb, real rotation_matrix[3][3],
                                   const std::vector<int> &filter)
  {
    // `filter` is of size number of nodes.  Value = 1 then rotate; value = 0 leave as is.

    // Get original coordinates...
    std::vector<double> coord;
    nb->get_field_data("mesh_model_coordinates", coord);
    size_t node_count = coord.size() / 3;

    // Do the rotation...
    for (size_t i = 0; i < node_count; i++) {
      if (filter.empty() || filter[i] == 1) {
        real x = coord[3 * i + 0];
        real y = coord[3 * i + 1];
        real z = coord[3 * i + 2];

        real xn = x * rotation_matrix[0][0] + y * rotation_matrix[1][0] + z * rotation_matrix[2][0];
        real yn = x * rotation_matrix[0][1] + y * rotation_matrix[1][1] + z * rotation_matrix[2][1];
        real zn = x * rotation_matrix[0][2] + y * rotation_matrix[1][2] + z * rotation_matrix[2][2];

        coord[3 * i + 0] = xn;
        coord[3 * i + 1] = yn;
        coord[3 * i + 2] = zn;
      }
    }

    // Output updated coordinates...
    nb->put_field_data("mesh_model_coordinates", coord);
  }

  void rotate_filtered_coordinates(Ioss::Region &region, real rotation_matrix[3][3],
                                   const std::vector<const Ioss::GroupingEntity *> &blocks)
  {
    const auto type = region.mesh_type();
    if (type == Ioss::MeshType::UNSTRUCTURED) {
      auto             filter = get_filtered_node_list(region, blocks);
      Ioss::NodeBlock *nb     = region.get_node_block("nodeblock_1");
      rotate_filtered_coordinates(nb, rotation_matrix, filter);
    }
    else if (type == Ioss::MeshType::STRUCTURED) {
      if (blocks.empty()) {
        auto sblocks = region.get_structured_blocks();
        for (const auto &sb : sblocks) {
          rotate_filtered_coordinates(sb, rotation_matrix, std::vector<int>());
        }
      }
      else {
        for (const auto &blk : blocks) {
          auto *sb = dynamic_cast<const Ioss::StructuredBlock *>(blk);
          if (sb != nullptr) {
            rotate_filtered_coordinates(sb, rotation_matrix, std::vector<int>());
          }
        }
      }
    }
  }

  void offset_filtered_coordinates(const Ioss::GroupingEntity *nb, real offset[3],
                                   const std::vector<int> &filter)
  {
    // `filter` is of size number of nodes.  Value = 1 transform; value = 0 leave as is.

    // Get original coordinates...
    std::vector<double> coord;
    nb->get_field_data("mesh_model_coordinates", coord);
    size_t node_count = coord.size() / 3;

    // Do the transformation...
    for (size_t i = 0; i < node_count; i++) {
      if (filter.empty() || filter[i] == 1) {
        coord[3 * i + 0] += offset[0];
        coord[3 * i + 1] += offset[1];
        coord[3 * i + 2] += offset[2];
      }
    }

    // Output updated coordinates...
    nb->put_field_data("mesh_model_coordinates", coord);
  }

  void offset_filtered_coordinates(Ioss::Region &region, real offset[3],
                                   const std::vector<const Ioss::GroupingEntity *> &blocks)
  {
    const auto type = region.mesh_type();
    if (type == Ioss::MeshType::UNSTRUCTURED) {
      auto             filter = get_filtered_node_list(region, blocks);
      Ioss::NodeBlock *nb     = region.get_node_block("nodeblock_1");
      offset_filtered_coordinates(nb, offset, filter);
    }
    else if (type == Ioss::MeshType::STRUCTURED) {
      if (blocks.empty()) {
        auto sblocks = region.get_structured_blocks();
        for (const auto &sb : sblocks) {
          offset_filtered_coordinates(sb, offset, std::vector<int>());
        }
      }
      else {
        for (const auto &blk : blocks) {
          auto *sb = dynamic_cast<const Ioss::StructuredBlock *>(blk);
          if (sb != nullptr) {
            offset_filtered_coordinates(sb, offset, std::vector<int>());
          }
        }
      }
    }
  }

  void scale_filtered_coordinates(const Ioss::GroupingEntity *nb, real scale[3],
                                  const std::vector<int> &filter)
  {
    // Get original coordinates...
    std::vector<double> coord;
    nb->get_field_data("mesh_model_coordinates", coord);
    size_t node_count = coord.size() / 3;

    // Do the transformation...
    for (size_t i = 0; i < node_count; i++) {
      if (filter.empty() || filter[i] == 1) {
        coord[3 * i + 0] *= scale[0];
        coord[3 * i + 1] *= scale[1];
        coord[3 * i + 2] *= scale[2];
      }
    }

    // Output updated coordinates...
    nb->put_field_data("mesh_model_coordinates", coord);
  }

  void scale_filtered_coordinates(Ioss::Region &region, real scale[3],
                                  const std::vector<const Ioss::GroupingEntity *> &blocks)
  {
    const auto type = region.mesh_type();
    if (type == Ioss::MeshType::UNSTRUCTURED) {
      auto             filter = get_filtered_node_list(region, blocks);
      Ioss::NodeBlock *nb     = region.get_node_block("nodeblock_1");
      scale_filtered_coordinates(nb, scale, filter);
    }
    else if (type == Ioss::MeshType::STRUCTURED) {
      if (blocks.empty()) {
        auto sblocks = region.get_structured_blocks();
        for (const auto &sb : sblocks) {
          scale_filtered_coordinates(sb, scale, std::vector<int>());
        }
      }
      else {
        for (const auto &blk : blocks) {
          auto *sb = dynamic_cast<const Ioss::StructuredBlock *>(blk);
          if (sb != nullptr) {
            scale_filtered_coordinates(sb, scale, std::vector<int>());
          }
        }
      }
    }
  }

  bool update_rotation_matrix(real rotation_matrix[3][3], const std::string &axis, double angle)
  {
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    if (Ioss::Utils::substr_equal(axis, "x")) {
      n1 = 1;
      n2 = 2;
      n3 = 0;
    }
    else if (Ioss::Utils::substr_equal(axis, "y")) {
      n1 = 2;
      n2 = 0;
      n3 = 1;
    }
    else if (Ioss::Utils::substr_equal(axis, "z")) {
      n1 = 0;
      n2 = 1;
      n3 = 2;
    }
    else {
      fmt::print(stderr, fg(fmt::color::red),
                 "ERROR: Requested rotation axis '{}' is not valid.  Must be 'x', 'y', or 'z'\n",
                 axis);
      return false;
    }

    static real degang  = std::atan2(0.0L, -1.0L) / 180.0;
    real        ang_rad = angle * degang;

    auto cosang = std::cos(ang_rad);
    auto sinang = std::sin(ang_rad);

    real by[3][3];
    real res[3][3];

    by[n1][n1] = cosang;
    by[n2][n1] = -sinang;
    by[n1][n3] = 0.0;
    by[n1][n2] = sinang;
    by[n2][n2] = cosang;
    by[n2][n3] = 0.0;
    by[n3][n1] = 0.0;
    by[n3][n2] = 0.0;
    by[n3][n3] = 1.0;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res[i][j] = rotation_matrix[i][0] * by[0][j] + rotation_matrix[i][1] * by[1][j] +
                    rotation_matrix[i][2] * by[2][j];
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rotation_matrix[i][j] = res[i][j];
      }
    }
    return true;
  }

  void modify_time(Ioss::Region &region, double scale, double offset)
  {
    size_t              ts_count = region.get_optional_property("state_count", 0);
    std::vector<double> times(ts_count);

    int exoid = region.get_database()->get_file_pointer();
    ex_get_all_times(exoid, Data(times));

    for (size_t step = 0; step < ts_count; step++) {
      times[step] = times[step] * scale + offset;
      ex_put_time(exoid, step + 1, &times[step]);
    }
    // Now update the `last_time_attribute`...
    auto max_time = *std::max_element(times.begin(), times.end());
    Ioex::update_last_time_attribute(exoid, max_time);
    (void)region.get_min_time(); // Triggers reloading region stateTimes vector.
  }
} // nameSpace
