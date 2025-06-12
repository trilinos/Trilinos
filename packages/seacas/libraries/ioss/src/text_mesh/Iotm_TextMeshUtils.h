// Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#pragma once

#include "iotm_export.h"

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <cctype>                                    // for toupper, isspace, isdigit
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector
#include <unordered_map>
#include <sstream>                       // for ostringstream
#include <iostream>
#include <functional>
#include <stdexcept>
#include <numeric>

#include "Iotm_TextMeshFuncs.h"
#include "Iotm_TextMeshDataTypes.h"
#include "Iotm_TextMeshEntityGroup.h"
#include "Iotm_TextMeshSideset.h"
#include "Iotm_TextMeshNodeset.h"
#include "Iotm_TextMeshAssembly.h"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace Iotm {
  namespace text_mesh {
    using ErrorHandler = std::function<void(const std::ostringstream &)>;

    class IOTM_EXPORT TextMeshLexer
    {
    public:
      TextMeshLexer() : m_currentIndex(0), m_token(""), m_isNumber(false) {}

      void set_input_string(const std::string &input)
      {
        m_input        = input;
        m_currentIndex = 0;
        read_next_token();
      }

      int get_int()
      {
        read_next_token();
        return std::stoi(m_oldToken);
      }

      unsigned get_unsigned()
      {
        read_next_token();
        return std::stoul(m_oldToken);
      }

      std::string get_string()
      {
        read_next_token();
        return make_uppercase(m_oldToken);
      }

      void get_newline() { read_next_token(); }

      bool has_token() const { return !m_token.empty(); }
      bool has_newline() const { return m_token == "\n"; }
      bool has_number() const { return has_token() && m_isNumber; }
      bool has_string() const { return has_token() && !has_number() && !has_newline(); }

    private:
      void read_next_token()
      {
        m_oldToken = m_token;

        if (char_is_newline()) {
          m_isNumber = false;
          m_token    = "\n";
          m_currentIndex++;
          return;
        }

        m_token    = "";
        m_isNumber = true;

        while (has_more_input()) {
          if (char_is_whitespace()) {
            m_currentIndex++;
            continue;
          }

          if (char_is_comma()) {
            m_currentIndex++;
            break;
          }

          if (char_is_newline()) {
            break;
          }

          m_isNumber &= char_is_digit();
          m_token += current_char();
          m_currentIndex++;
        }
      }

      bool has_more_input() { return m_currentIndex < m_input.size(); }

      bool char_is_whitespace() { return current_char() == ' '; }
      bool char_is_comma() { return current_char() == ','; }
      bool char_is_newline() { return current_char() == '\n'; }
      bool char_is_digit() { return std::isdigit(static_cast<unsigned char>(current_char())); }

      char current_char() { return m_input[m_currentIndex]; }

      std::string make_uppercase(std::string str)
      {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        return str;
      }

      std::string m_input{};
      unsigned    m_currentIndex{0};

      std::string m_oldToken{};
      std::string m_token{};

      bool m_isNumber{false};
    };

    template <typename EntityId, typename Topology> class TextMeshOptionParser
    {
    private:
      static constexpr int INVALID_DIMENSION = -1;
      static constexpr int DEFAULT_DIMENSION = 3;

      enum ParsedOptions {
        PARSED_NONE        = 0,
        PARSED_DIMENSION   = 1L << 0,
        PARSED_COORDINATES = 1L << 1,
        PARSED_SIDESET     = 1L << 2,
        PARSED_NODESET     = 1L << 3,
        PARSED_ASSEMBLY    = 1L << 4
      };

    public:
      TextMeshOptionParser(TextMeshData<EntityId, Topology> &data, unsigned enforcedDimension)
          : m_parsedOptionMask(PARSED_NONE), m_parsedDimension(INVALID_DIMENSION),
            m_constructorEnforcedDimension(enforcedDimension), m_data(data)
      {
      }

      explicit TextMeshOptionParser(TextMeshData<EntityId, Topology> &data)
          : m_parsedOptionMask(PARSED_NONE), m_parsedDimension(INVALID_DIMENSION),
            m_constructorEnforcedDimension(INVALID_DIMENSION), m_data(data)
      {
      }

      void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

      std::string get_mesh_connectivity_description() const
      {
        return m_meshConnectivityDescription;
      }

      void initialize_parse(const std::string &parameters)
      {
        if (!parameters.empty()) {
          std::vector<std::string> optionGroups = get_tokens(parameters, "|");
          parse_options(optionGroups);

          m_meshConnectivityDescription = optionGroups[0];
        }

        validate_dimension();
        set_dimension();
      }

      void finalize_parse()
      {
        set_coordinates();
        m_data.partIds.finalize_parse();
        m_data.sidesets.finalize_parse(m_data);
        m_data.nodesets.finalize_parse();
        m_data.assemblies.finalize_parse();
        validate_sidesets();
        validate_nodesets();
        validate_assemblies();
      }

    private:
      bool parsed_dimension_provided() { return m_parsedOptionMask & PARSED_DIMENSION; }

      bool enforced_dimension_provided()
      {
        return m_constructorEnforcedDimension != INVALID_DIMENSION;
      }

      void validate_dimension()
      {
        if (enforced_dimension_provided()) {
          if (parsed_dimension_provided() && m_constructorEnforcedDimension != m_parsedDimension) {
            std::ostringstream errmsg;
            errmsg << "Error!  An enforced dimension of " << m_constructorEnforcedDimension
                   << " was provided but does not match the parsed value of " << m_parsedDimension
                   << ".";
            m_errorHandler(errmsg);
          }
        }
      }

      void set_dimension()
      {
        if (enforced_dimension_provided()) {
          m_data.spatialDim = m_constructorEnforcedDimension;
        }
        else if (parsed_dimension_provided()) {
          m_data.spatialDim = m_parsedDimension;
        }
        else {
          m_data.spatialDim = DEFAULT_DIMENSION;
        }
      }

      void parse_dimension_option(const std::vector<std::string> &option)
      {
        if (parsed_dimension_provided()) {
          std::ostringstream errmsg;
          errmsg << "Spatial dimension has already been parsed! Check syntax.";
          m_errorHandler(errmsg);
        }

        if (option.size() == 2) {
          m_parsedDimension = std::stoull(option[1]);
          if (m_parsedDimension != 2 && m_parsedDimension != 3) {
            std::ostringstream errmsg;
            errmsg << "Error!  Parsed spatial dimension (" << m_parsedDimension
                   << " not defined to be 2 or 3.";
            m_errorHandler(errmsg);
          }

          m_parsedOptionMask |= PARSED_DIMENSION;
        }
        else {
          std::ostringstream errmsg;
          errmsg << "Error!  Invalid spatial dimension syntax.";
          m_errorHandler(errmsg);
        }
      }

      void deallocate_raw_coordinates()
      {
        std::vector<double> swapVectorForDeallocation;
        m_rawCoordinates.swap(swapVectorForDeallocation);
      }

      void set_coordinates()
      {
        if (parsed_coordinates_provided()) {
          m_data.coords.set_coordinate_data(m_data.spatialDim, m_data.nodeIds, m_rawCoordinates);
          deallocate_raw_coordinates();
        }
      }

      bool parsed_coordinates_provided() { return m_parsedOptionMask & PARSED_COORDINATES; }

      void parse_coordinates_option(const std::vector<std::string> &coordinatesOptionGroup)
      {
        if (parsed_coordinates_provided()) {
          std::ostringstream errmsg;
          errmsg << "Coordinates have already been parsed! Check syntax.";
          m_errorHandler(errmsg);
        }

        if (coordinatesOptionGroup.size() > 1) {
          const std::vector<std::string> &coordinateTokens =
              get_tokens(coordinatesOptionGroup[1], ",");
          m_rawCoordinates.reserve(coordinateTokens.size());
          for (const auto &token : coordinateTokens) {
            double coord = std::stod(token);
            m_rawCoordinates.push_back(coord);
          }

          m_parsedOptionMask |= PARSED_COORDINATES;
        }
      }

      template <typename DataType>
      void check_name_collision_with_entity_sets(const EntityGroupData<DataType> &groupData,
                                                 const std::string               &entityType,
                                                 const std::set<std::string>     &entitySetNames)
      {
        std::string groupName = groupData.name;
        convert_to_uppercase(groupName);

        if (entitySetNames.count(groupName) > 0) {
          std::ostringstream errmsg;
          errmsg << "Error! " << groupData.type << " with id: " << groupData.id
                 << " and name: " << groupData.name << " is referencing " << entityType
                 << " with same name.";
          m_errorHandler(errmsg);
        }
      }

      template <typename SrcDataGroup, typename DestDataGroup>
      void check_name_collision_with_group(const SrcDataGroup  &srcGroup,
                                           const DestDataGroup &destGroup)
      {
        std::set<std::string> groupNames = transform_to_set(destGroup.get_part_names());

        for (const auto &srcGroupData : srcGroup.get_group_data()) {
          check_name_collision_with_entity_sets(srcGroupData, destGroup.get_group_type(),
                                                groupNames);
        }
      }

      void check_sideset_element_reference()
      {
        for (const SidesetData<EntityId, Topology> &sidesetData :
             m_data.sidesets.get_group_data()) {
          for (const std::pair<EntityId, int> &elemSidePair : sidesetData.data) {
            EntityId id = elemSidePair.first;
            if (!std::binary_search(m_data.elementDataVec.begin(), m_data.elementDataVec.end(),
                                    id)) {
              std::ostringstream errmsg;
              errmsg << "Error!  Sideset with id: " << sidesetData.id
                     << " and name: " << sidesetData.name << " has reference to invalid element '"
                     << id << "'.";
              m_errorHandler(errmsg);
            }
          }
        }
      }

      void check_sideset_name_collision()
      {
        check_name_collision_with_group(m_data.sidesets, m_data.partIds);
        check_name_collision_with_group(m_data.sidesets, m_data.nodesets);
        check_name_collision_with_group(m_data.sidesets, m_data.assemblies);
      }

      void validate_sidesets()
      {
        check_sideset_element_reference();
        check_sideset_name_collision();
      }

      void check_nodeset_node_reference()
      {
        for (const NodesetData<EntityId> &nodesetData : m_data.nodesets.get_group_data()) {
          for (const EntityId nodeId : nodesetData.data) {
            if (m_data.nodeIds.count(nodeId) == 0) {
              std::ostringstream errmsg;
              errmsg << "Error!  Nodeset with id: " << nodesetData.id
                     << " and name: " << nodesetData.name << " has reference to invalid node '"
                     << nodeId << "'.";
              m_errorHandler(errmsg);
            }
          }
        }
      }

      void check_nodeset_name_collision()
      {
        check_name_collision_with_group(m_data.nodesets, m_data.partIds);
        check_name_collision_with_group(m_data.nodesets, m_data.sidesets);
        check_name_collision_with_group(m_data.nodesets, m_data.assemblies);
      }

      void validate_nodesets()
      {
        check_nodeset_node_reference();
        check_nodeset_name_collision();
      }

      template <typename T>
      void check_assembly_member_reference_in_group(const AssemblyData &assemblyData,
                                                    const T            &group)
      {
        for (const std::string &entry : assemblyData.data) {
          if (!group.is_registered(entry)) {
            std::ostringstream errmsg;
            errmsg << "Error!  Assembly with id: " << assemblyData.id
                   << " and name: " << assemblyData.name << " has reference to invalid "
                   << group.get_group_type() << " '" << entry << "'.";
            m_errorHandler(errmsg);
          }
        }
      }

      void check_assembly_member_reference()
      {
        for (const AssemblyData &assemblyData : m_data.assemblies.get_group_data()) {
          const AssemblyType assemblyType = assemblyData.get_assembly_type();

          switch (assemblyType) {
          case AssemblyType::BLOCK:
            check_assembly_member_reference_in_group(assemblyData, m_data.partIds);
            break;
          case AssemblyType::SIDESET:
            check_assembly_member_reference_in_group(assemblyData, m_data.sidesets);
            break;
          case AssemblyType::NODESET:
            check_assembly_member_reference_in_group(assemblyData, m_data.nodesets);
            break;
          case AssemblyType::ASSEMBLY:
            check_assembly_member_reference_in_group(assemblyData, m_data.assemblies);
            break;
          default:
            std::ostringstream errmsg;
            errmsg << "Error!  Assembly with id: " << assemblyData.id
                   << " and name: " << assemblyData.name << " does not have a valid assembly type '"
                   << assemblyType << "'.";
            m_errorHandler(errmsg);
          }
        }
      }

      void check_assembly_name_collision()
      {
        check_name_collision_with_group(m_data.assemblies, m_data.partIds);
        check_name_collision_with_group(m_data.assemblies, m_data.sidesets);
        check_name_collision_with_group(m_data.assemblies, m_data.nodesets);
      }

      void check_assembly_cyclic_dependency()
      {
        for (const std::string &assembly : m_data.assemblies.get_part_names()) {
          if (m_data.assemblies.is_cyclic(assembly)) {
            std::ostringstream errmsg;
            errmsg << "Error!  Assembly with name: '" << assembly << "' has a cyclic dependency.";
            m_errorHandler(errmsg);
          }
        }
      }

      void validate_assemblies()
      {
        check_assembly_member_reference();
        check_assembly_name_collision();
        check_assembly_cyclic_dependency();
      }

      void parse_sideset_option(const std::vector<std::string> &sidesetOptionGroup)
      {
        if (sidesetOptionGroup.size() > 1) {
          SidesetParser<EntityId> parser;
          parser.set_error_handler(m_errorHandler);
          parser.parse(sidesetOptionGroup[1]);
          parser.verify_parse();

          SidesetData<EntityId, Topology> *sideset =
              m_data.sidesets.add_group_data(parser.get_name(), parser.get_sideset_data());
          sideset->set_split_type(parser.get_split_type());
          sideset->set_skin_blocks(parser.get_skin_blocks());
          m_parsedOptionMask |= PARSED_SIDESET;
        }
      }

      void parse_nodeset_option(const std::vector<std::string> &nodesetOptionGroup)
      {
        if (nodesetOptionGroup.size() > 1) {
          NodesetParser<EntityId> parser;
          parser.set_error_handler(m_errorHandler);
          parser.parse(nodesetOptionGroup[1]);

          m_data.nodesets.add_group_data(parser.get_name(), parser.get_nodeset_data());
          m_parsedOptionMask |= PARSED_NODESET;
        }
      }

      void parse_assembly_option(const std::vector<std::string> &assemblyOptionGroup)
      {
        if (assemblyOptionGroup.size() > 1) {
          AssemblyParser parser;
          parser.set_error_handler(m_errorHandler);
          parser.parse(assemblyOptionGroup[1]);
          parser.verify_parse();

          AssemblyData *assembly =
              m_data.assemblies.add_group_data(parser.get_name(), parser.get_assembly_data());
          assembly->set_assembly_type(parser.get_assembly_type());
          m_parsedOptionMask |= PARSED_ASSEMBLY;
        }
      }

      void print_help_message(std::ostream &out = std::cout)
      {
        out << "\nValid Options for TextMesh parameter string:\n"
               "\tPROC_ID,ELEM_ID,TOPOLOGY,{NODE CONNECTIVITY LIST}[,PART_NAME[,PART_ID]] "
               "(specifies "
               "element list .. first "
               "argument)\n"
               "\t|coordinates:x_1,y_1[,z_1], x_2,y_2[,z_2], ...., x_n,y_n[,z_n] (specifies "
               "coordinate data)\n"
               "\t|sideset:[name=<name>;] data=elem_1,side_1,elem_2,side_2,....,elem_n,side_n; "
               "[split=<block|topology|none>;] [skin=<block list|all>;]"
               "(specifies sideset data)\n"
               "\t|nodeset:[name=<name>;] data=node_1,node_2,....,node_n (specifies nodeset data)\n"
               "\t|assembly:[name=<name>;] type=<assembly|block|sideset|nodeset>; "
               "member=member_1,...,member_n (specifies assembly hierarchy)\n"
               "\t|dimension:spatialDimension (specifies spatial dimension .. default is 3)\n"
               "\t|help -- show this list\n\n";
      }

      void handle_unrecognized_option(const std::string &optionType)
      {
        std::ostringstream errmsg;
        errmsg << "ERROR: Unrecognized option '" << optionType << "'.  It will be ignored.\n";
        m_errorHandler(errmsg);
      }

      void parse_options(const std::vector<std::string> &optionGroups)
      {
        for (size_t i = 1; i < optionGroups.size(); i++) {
          std::vector<std::string> optionGroup = get_tokens(optionGroups[i], ":");
          std::string              optionType  = optionGroup[0];
          convert_to_lowercase(optionType);

          if (optionType == "coordinates") {
            parse_coordinates_option(optionGroup);
          }
          else if (optionType == "dimension") {
            parse_dimension_option(optionGroup);
          }
          else if (optionType == "sideset") {
            parse_sideset_option(optionGroup);
          }
          else if (optionType == "nodeset") {
            parse_nodeset_option(optionGroup);
          }
          else if (optionType == "assembly") {
            parse_assembly_option(optionGroup);
          }
          else if (optionType == "help") {
            print_help_message();
          }
          else {
            handle_unrecognized_option(optionType);
          }
        }
      }

      unsigned long m_parsedOptionMask{PARSED_NONE};

      int m_parsedDimension{INVALID_DIMENSION};
      int m_constructorEnforcedDimension{INVALID_DIMENSION};

      std::string m_meshConnectivityDescription{};

      std::vector<double> m_rawCoordinates{};
      ErrorHandler        m_errorHandler;

      TextMeshData<EntityId, Topology> &m_data;
    };

    template <typename EntityId, typename TopologyMapping> class TextMeshParser
    {
    private:
      using Topology = typename TopologyMapping::Topology;

    public:
      explicit TextMeshParser(unsigned enforcedDimension)
          : m_optionParser(m_data, enforcedDimension)
      {
        initialize_constructor();
      }

      TextMeshParser() : m_optionParser(m_data) { initialize_constructor(); }

      TextMeshData<EntityId, Topology> parse(const std::string &meshDescription)
      {
        initialize_parse(meshDescription);
        parse_description();
        finalize_parse();
        return m_data;
      }

      void set_error_handler(ErrorHandler errorHandler)
      {
        m_errorHandler = errorHandler;
        m_data.partIds.set_error_handler(errorHandler);
        m_data.coords.set_error_handler(errorHandler);
        m_data.sidesets.set_error_handler(errorHandler);
        m_data.nodesets.set_error_handler(errorHandler);
        m_data.assemblies.set_error_handler(errorHandler);
        m_optionParser.set_error_handler(errorHandler);
      }

    private:
      void initialize_constructor()
      {
        ErrorHandler errorHandler = [](const std::ostringstream &errmsg) {
          default_error_handler(errmsg);
        };
        set_error_handler(errorHandler);
        m_topologyMapping.initialize_topology_map();
      }

      void initialize_connectivity_parse(const std::string &meshDescription)
      {
        m_lexer.set_input_string(meshDescription);
        m_lineNumber = 1;
        validate_required_field(m_lexer.has_token());
      }

      void initialize_parse(const std::string &meshDescription)
      {
        m_optionParser.initialize_parse(meshDescription);
        initialize_connectivity_parse(m_optionParser.get_mesh_connectivity_description());
      }

      void finalize_parse() { m_optionParser.finalize_parse(); }

      void parse_description()
      {
        while (m_lexer.has_token()) {
          ElementData<EntityId, Topology> elemData = parse_element();
          m_data.add_element(elemData);

          validate_no_extra_fields();
          parse_newline();
        }

        std::sort(m_data.elementDataVec.begin(), m_data.elementDataVec.end(),
                  ElementDataLess<EntityId, Topology>());
      }

      ElementData<EntityId, Topology> parse_element()
      {
        ElementData<EntityId, Topology> elem;
        elem.proc       = parse_proc_id();
        elem.identifier = parse_elem_id();
        elem.topology   = parse_topology();
        elem.nodeIds    = parse_node_ids(elem.topology);
        elem.partName   = parse_part(elem.topology);
        return elem;
      }

      int parse_proc_id()
      {
        validate_required_field(m_lexer.has_number());
        return parse_int();
      }

      EntityId parse_elem_id()
      {
        validate_required_field(m_lexer.has_number());
        return parse_unsigned();
      }

      Topology parse_topology()
      {
        validate_required_field(m_lexer.has_string());
        std::string topologyName = parse_string();

        Topology topology = m_topologyMapping.topology(topologyName);
        validate_topology(topology, topologyName);

        return topology;
      }

      std::vector<EntityId> parse_node_ids(const Topology &topology)
      {
        std::vector<EntityId> nodeIds;
        while (m_lexer.has_number()) {
          nodeIds.push_back(parse_unsigned());
        }
        validate_node_count(topology, nodeIds.size());
        return nodeIds;
      }

      std::string parse_part(const Topology &topology)
      {
        std::string partName;

        if (m_lexer.has_string()) {
          partName = parse_string();
        }
        else {
          partName = "block_" + topology.name();
        }

        if (m_lexer.has_number()) {
          unsigned partId = parse_unsigned();
          m_data.partIds.register_part_name_with_id(partName, partId);
        }
        else {
          m_data.partIds.register_part_name(partName);
        }

        return partName;
      }

      int         parse_int() { return m_lexer.get_int(); }
      unsigned    parse_unsigned() { return m_lexer.get_unsigned(); }
      std::string parse_string() { return m_lexer.get_string(); }

      void parse_newline()
      {
        m_lexer.get_newline();
        m_lineNumber++;
      }

      void validate_required_field(bool hasNextRequiredField)
      {
        if (!hasNextRequiredField) {
          std::ostringstream errmsg;
          errmsg
              << "Error!  Each line must contain the following fields (with at least one node):  "
                 "Processor, GlobalId, Element Topology, NodeIds.  Error on line "
              << m_lineNumber << ".";
          m_errorHandler(errmsg);
        }
      }

      void validate_no_extra_fields()
      {
        bool requiredCondition = !m_lexer.has_token() || m_lexer.has_newline();
        if (!requiredCondition) {
          std::ostringstream errmsg;
          errmsg << "Error!  Each line should not contain more than the following fields (with at "
                    "least one node):  "
                    "Processor, GlobalId, Element Topology, NodeIds, Part Name, PartId.  "
                    "Error on line "
                 << m_lineNumber << ".";
          m_errorHandler(errmsg);
        }
      }

      void validate_topology(const Topology &topology, const std::string &providedName)
      {
        if (topology == m_topologyMapping.invalid_topology()) {
          std::ostringstream errmsg;
          errmsg << "Error!  Topology = >>" << providedName << "<< is invalid from line "
                 << m_lineNumber << ".";
          m_errorHandler(errmsg);
        }

        if (!topology.defined_on_spatial_dimension(m_data.spatialDim)) {
          std::ostringstream errmsg;
          errmsg << "Error on input line " << m_lineNumber << ".  Topology = " << topology
                 << " is not defined on spatial dimension = " << m_data.spatialDim
                 << " set in parser.";
          m_errorHandler(errmsg);
        }
      }

      void validate_node_count(const Topology &topology, size_t numNodes)
      {
        size_t numTopologyNodes = topology.num_nodes();
        if (numNodes != numTopologyNodes) {
          std::ostringstream errmsg;
          errmsg << "Error!  The input line appears to contain " << numNodes
                 << " nodes, but the topology " << topology << " needs " << numTopologyNodes
                 << " nodes on line " << m_lineNumber << ".";
          m_errorHandler(errmsg);
        }
      }

      unsigned                         m_lineNumber{0};
      TextMeshData<EntityId, Topology> m_data;
      TextMeshLexer                    m_lexer;
      TopologyMapping                  m_topologyMapping;

      ErrorHandler m_errorHandler;

      TextMeshOptionParser<EntityId, Topology> m_optionParser;
    };

  } // namespace text_mesh
} // namespace Iotm
