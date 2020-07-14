// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
// 
// See packages/seacas/LICENSE for details

#include "modify_interface.h"

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include <Ionit_Initializer.h>
#include <Ioss_Assembly.h>
#include <Ioss_Blob.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_CommSet.h>
#include <Ioss_CoordinateFrame.h>
#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_EdgeBlock.h>
#include <Ioss_EdgeSet.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementSet.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_FaceBlock.h>
#include <Ioss_FaceSet.h>
#include <Ioss_Field.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Getline.h>
#include <Ioss_GroupingEntity.h>
#include <Ioss_IOFactory.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_SideBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Utils.h>
#include <Ioss_VariableType.h>
#include <tokenize.h>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#if defined(SEACAS_HAVE_EXODUS)
#include <exodusII.h>
#include <exodus/Ioex_Internals.h>
#include <exodus/Ioex_Utils.h>
#endif

#if defined(SEACAS_HAVE_CGNS)
#include <Iocgns_Utils.h>
#endif

#if defined(_MSC_VER)
#include <io.h>
#define isatty _isatty
#endif
// ========================================================================

namespace {
  std::string codename;
  std::string version = "0.92 (2020-05-20)";

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
  bool           handle_delete(const std::vector<std::string> &tokens, Ioss::Region &region);
  void           handle_list(const std::vector<std::string> &tokens, const Ioss::Region &region,
                             bool show_attribute = false);
  void           handle_graph(const std::vector<std::string> &tokens, const Ioss::Region &region);
  bool           handle_assembly(const std::vector<std::string> &tokens, Ioss::Region &region,
                                 bool allow_modify);
  bool           handle_attribute(const std::vector<std::string> &tokens, Ioss::Region &region);
  void           update_assembly_info(Ioss::Region &region, const Modify::Interface &interFace);

  void set_db_properties(const Modify::Interface &interFace, Ioss::DatabaseIO *dbi);

  void info_entity(const Ioss::StructuredBlock *sb, bool show_property = false);
  void info_entity(const Ioss::NodeBlock *nb, bool show_property = false);
  void info_entity(const Ioss::ElementBlock *eb, bool show_property = false);
  void info_entity(const Ioss::NodeSet *ns, bool show_property = false);
  void info_entity(const Ioss::SideSet *ss, bool show_property = false);
  void info_entity(const Ioss::Assembly *as, bool show_property = false);
  void info_entity(const Ioss::Blob *blob, bool show_property = false);
  void info_entity(const Ioss::Region &region, bool show_property = false);

  template <typename T>
  void info_entities(const std::vector<T *> &entities, const std::vector<std::string> &tokens,
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
          const T *   ge     = dynamic_cast<const T *>(entity);
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
    int64_t id = -1;
    if (entity->property_exists("id")) {
      id = entity->get_property("id").get_int();
    }
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

  Ioss::PropertyManager set_properties(const Modify::Interface &interFace)
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
  properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND));

  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::WRITE_RESTART,
                                                  (MPI_Comm)MPI_COMM_WORLD, properties);
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

  while (1) {
    std::string input;
    if (from_term) {
      fmt::print(fg(fmt::terminal_color::magenta), "\n");
      const char *cinput = getline_int("COMMAND> ");
      if (cinput[0] == '\0') {
        break;
      }
      if (cinput) {
        gl_histadd(cinput);
      }
      input = cinput;
    }
    else {
      std::getline(std::cin, input);
    }

    // NOTE: getline_int returns the trailing '\n'
    auto tokens = Ioss::tokenize(std::string(input), " ,\n");
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
    else {
      fmt::print(stderr, fg(fmt::color::yellow), "\tWARNING: Unrecognized command: {}\n",
                 tokens[0]);
    }
  }

  if (changed) {
    update_assembly_info(region, interFace);
  }
  else {
    fmt::print(fg(fmt::color::cyan), "\n\t*** Database unchanged. No update required.\n");
  }
  fmt::print("\n{} execution successful.\n", codename);
  return EXIT_SUCCESS;
}

namespace {
  void info_entity(const Ioss::StructuredBlock *sb, bool show_property)
  {
    fmt::print("\n{} {}", name(sb), sb->get_property("ni_global").get_int());

    int64_t num_dim = sb->get_property("component_degree").get_int();
    if (num_dim > 1) {
      fmt::print("x{}", sb->get_property("nj_global").get_int());
    }
    if (num_dim > 2) {
      fmt::print("x{}", sb->get_property("nk_global").get_int());
    }

    fmt::print(" [{}x{}x{}, Offset = {}, {}, {}] ", sb->get_property("ni").get_int(),
               sb->get_property("nj").get_int(), sb->get_property("nk").get_int(),
               sb->get_property("offset_i").get_int(), sb->get_property("offset_j").get_int(),
               sb->get_property("offset_k").get_int());

    int64_t num_cell = sb->get_property("cell_count").get_int();
    int64_t num_node = sb->get_property("node_count").get_int();
    fmt::print("{:14n} cells, {:14n} nodes\n", num_cell, num_node);
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
    for (const auto mem : as->get_members()) {
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

    std::string type       = eb->get_property("topology_type").get_string();
    int64_t     num_attrib = eb->get_property("attribute_count").get_int();
    fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14n} elements, {:3d} attributes.\n", name(eb),
               id(eb), type, num_elem, num_attrib);
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
    for (auto fb : fbs) {
      int64_t count      = fb->entity_count();
      int64_t num_attrib = fb->get_property("attribute_count").get_int();
      int64_t num_dist   = fb->get_property("distribution_factor_count").get_int();
      fmt::print("\t{}, {:8n} sides, {:3d} attributes, {:8n} distribution factors.\n", name(fb),
                 count, num_attrib, num_dist);
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
    fmt::print("\n{} id: {:6d}, {:8n} nodes, {:3d} attributes, {:8n} distribution factors.\n",
               name(ns), id(ns), count, num_attrib, num_dist);
    if (show_property) {
      Ioss::Utils::info_property(ns, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void info_entity(const Ioss::NodeBlock *nb, bool show_property)
  {
    int64_t num_nodes  = nb->entity_count();
    int64_t num_attrib = nb->get_property("attribute_count").get_int();
    fmt::print("\n{} {:14n} nodes, {:3d} attributes.\n", name(nb), num_nodes, num_attrib);
    if (show_property) {
      Ioss::Utils::info_property(nb, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
    }
  }

  void set_db_properties(const Modify::Interface &interFace, Ioss::DatabaseIO *dbi)
  {
    std::string inpfile = interFace.filename();

    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }
  }

  void handle_help(const std::string &topic)
  {
    bool all = Ioss::Utils::substr_equal(topic, "help");
    if (all) {
      fmt::print("\n\tHELP [list | assembly | graph | attribute | regex]\n");
      fmt::print("\n\tEND | EXIT\n");
      fmt::print("\t\tEnd command input and output changed assembly definitions (if any).\n");
      fmt::print("\n\tQUIT\n");
      fmt::print("\t\tEnd command input and exit with no changes to database.\n");

      fmt::print("\n\tALLOW MODIFICATIONS\n");
      fmt::print("\t\tBy default, io_modify will only allow creation of new assemblies.\n"
                 "\t\tIf this command is specified, then can modify assemblies that already exist "
                 "in database.\n"
                 "\t\tThis will cause the database to be rewritten. Without this option, it is "
                 "updated in place.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "regex")) {
      fmt::print("\n\tRegular Expression help (used in ASSEMBLY MATCHES and LIST MATCHES and "
                 "ATTRIBUTE LIST MATCHES options)\n"
                 "\tSupports \"POSIX Extended Regular Expressions\" syntax\n"
                 "\tSee https://www.regular-expressions.info/posix.html\n"
                 "\tQuickStart: https://www.regular-expressions.info/quickstart.html\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "list")) {
      fmt::print(
          "\n\tLIST summary|elementblock|block|structuredblock|assembly|nodeset|sideset|blob\n");
      fmt::print("\n\tLIST elementblock|block|structuredblock|assembly|nodeset|sideset|blob "
                 "{{names...}}\n\n");
      fmt::print("\n\tLIST elementblock|block|structuredblock|assembly|nodeset|sideset|blob "
                 "MATCHES {{regex}}\n\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "assembly")) {
      fmt::print("\n\tFor all commands, if an assembly named `name` does not exist, it will be "
                 "created.\n");
      fmt::print("\tASSEMBLY {{name}}\n");
      fmt::print("\t\tCreates an empty assembly named `name` if it does not exist.\n");

      fmt::print("\n\tASSEMBLY {{name}} ADD {{name1}} {{name2}} ... {{namen}}\n");
      fmt::print("\t\tAdds the specified entities to the assembly.  All entities must be the same "
                 "type.\n");

      fmt::print("\n\tASSEMBLY {{name}} REMOVE {{name1}} {{name2}} ... {{namen}}\n");
      fmt::print("\t\tRemoves the specified entities from the assembly.\n");

      fmt::print("\n\tASSEMBLY {{name}} TYPE {{type}} MATCHES {{regex}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose name matches the {{regex}} will be added.\n");

      fmt::print("\n\tASSEMBLY {{name}} TYPE {{type}} NAMED {{list of one or more names}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose names are listed will be added.\n");

      fmt::print("\n\tASSEMBLY {{name}} TYPE {{type}} RANGE {{id}} TO {{id}} BY {{step}}\n");
      fmt::print("\t\tAdds the entities of the specified type to the assembly.\n"
                 "\t\tAll entities whose id matches the specified range will be added.\n"
                 "\t\tNo message will be output for ids not matching an entity.\n");

      fmt::print("\n\tASSEMBLY {{name}} TYPE {{type}} IDS {{id}}, {{id2}}, ..., {{idn}}\n");
      fmt::print(
          "\t\tAdds the entities of the specified type to the assembly.\n"
          "\t\tAll entities whose id matches an id in the list will be added.\n"
          "\t\tA warning message will be output if there is no entity with the requested id.\n");

      fmt::print("\n\tDELETE {{name}}\n");
      fmt::print("\t\tRemove the assembly with the specified name.\n"
                 "\t\tCurrently only supported for assemblies created during this execution; not "
                 "for assemblies\n"
                 "\t\texisting on the input database.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "graph")) {
      fmt::print("\n\tGRAPH OUTPUT [filename]\n");
      fmt::print(
          "\t\tCreate a 'dot' input file with the structure of the assembly graph.\n"
          "\t\tFile is named 'filename' or defaults to 'assembly.dot' if filename not given.\n");
      fmt::print("\tGRAPH CHECK\n");
      fmt::print("\t\tCheck validity of assembly graph--are there any cycles.\n");
    }
    if (all || Ioss::Utils::substr_equal(topic, "attribute")) {
      fmt::print("\n\tATTRIBUTE {{ent_name}} ADD {{att_name}} STRING {{values...}}\n");
      fmt::print("\tATTRIBUTE {{ent_name}} ADD {{att_name}} DOUBLE {{values...}}\n");
      fmt::print("\tATTRIBUTE {{ent_name}} NAME {{att_name}} INTEGER {{values...}}\n");
      fmt::print("\t\tAdd an attribute to the specified entity ('type' and 'name').\n"
                 "\t\tThe attribute will be named 'att_name' with value(s) 'values...'\n"
                 "\t\tCan also modify the values of an existing attribute.'\n");
      fmt::print("\tATTRIBUTE LIST {{ent_name...}}\n"
                 "\t\tList attributes for the selected entities\n");
      fmt::print("\tATTRIBUTE {{ent_type}} LIST\n"
                 "\t\tList attributes for all entities in the specified entity type\n");
      fmt::print("\tATTRIBUTE {{ent_type}} MATCH {{regex}}\n"
                 "\t\tList attributes for all entities in the specified entity type whose name "
                 "matches the regex.\n");
    }
  }

  void handle_list(const std::vector<std::string> &tokens, const Ioss::Region &region,
                   bool show_attribute)
  {
    if (tokens.size() > 1) {
      if (Ioss::Utils::substr_equal(tokens[1], "summary")) {
        region.output_summary(std::cout);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "region")) {
        info_entity(region, show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "elementblock") ||
               Ioss::Utils::substr_equal(tokens[1], "block")) {
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
      else if (Ioss::Utils::substr_equal(tokens[1], "sideset") ||
               Ioss::Utils::substr_equal(tokens[1], "sset")) {
        const auto &entities = region.get_sidesets();
        info_entities(entities, tokens, region, "SideSets", show_attribute);
      }
      else if (Ioss::Utils::substr_equal(tokens[1], "blobs")) {
        const auto &entities = region.get_blobs();
        info_entities(entities, tokens, region, "Blobs", show_attribute);
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
    std::map<std::string, int>    m_vertices;
    std::vector<std::string>      m_vertex;
    std::vector<std::vector<int>> m_adj; // Pointer to an array containing adjacency std::lists

    bool is_cyclic_internal(int v, std::vector<bool> &visited, std::vector<bool> &recStack);
    int  vertex(const std::string &node);
    int  size() const { return (int)m_vertices.size(); }

  public:
    Graph()  = default; // Constructor
    ~Graph() = default;

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

  void handle_graph(const std::vector<std::string> &tokens, const Ioss::Region &region)
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

  bool handle_delete(const std::vector<std::string> &tokens, Ioss::Region &region)
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
                     "deleteable.\n",
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

  bool handle_attribute(const std::vector<std::string> &tokens, Ioss::Region &region)
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
        fmt::print(stderr, fg(fmt::color::red),
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
      std::string att_name = tokens[3];

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
        std::vector<int> values(value_count);
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
    else {
      fmt::print(stderr, fg(fmt::color::red), "ERROR: Unrecognized attribute command.\n");
      handle_help("attribute");
      return false;
    }
    return false;
  }

  bool handle_assembly(const std::vector<std::string> &tokens, Ioss::Region &region,
                       bool allow_modify)
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
  void           update_cgns_assembly_info(Ioss::Region &region, const Modify::Interface &interFace)
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
	  fmt::print(stderr, fg(fmt::color::yellow), "WARNING: Can not modify existing assembly '{}'  yet.\n",
		     assembly->name());
        }
      }
    }
  }
#endif

#if defined(SEACAS_HAVE_EXODUS)
  void           update_exodus_assembly_info(Ioss::Region &region, const Modify::Interface &interFace)
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
      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(file_type, out_file, Ioss::WRITE_RESTART,
                                                      (MPI_Comm)MPI_COMM_WORLD, properties);

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
#endif

  void update_assembly_info(Ioss::Region &region, const Modify::Interface &interFace)
  {
    // Determine type of underlying database...
    const auto type = region.get_database()->get_format();
    if (type == "Exodus") {
#if defined(SEACAS_HAVE_EXODUS)
      update_exodus_assembly_info(region, interFace);
#else
      fmt::print(stderr, fg(fmt::color::red), "ERROR: Exodus capability is not enabled.\n");
#endif
    }
    else if (type == "CGNS") {
#if defined(SEACAS_HAVE_CGNS)
      update_cgns_assembly_info(region, interFace);
#else
      fmt::print(stderr, fg(fmt::color::red), "ERROR: CGNS capability is not enabled.\n");
#endif    
    }
    else {
      fmt::print(stderr, fg(fmt::color::red), "ERROR: Can not modify the database '{}' of type '{}'.\n",
		 interFace.filename(), type);
    }
  }

} // nameSpace
