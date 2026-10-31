// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CSG_Parser.hpp>

#include <Akri_CSG_Tree.hpp>
#include <Akri_CSG_Tree_Output.hpp>
#include <Akri_Surface.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Akri_Surface_Parser.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>

namespace krino {
namespace CSG_Parser {

static Parser::Node get_sequence_nodes(const std::vector<std::string> & sequenceNames, const Parser::Node& node)
{
  for (auto & seqName : sequenceNames)
  {
    Parser::Node seqNode = node.get_sequence_if_present(seqName);
    if (seqNode)
    {
      return seqNode;
    }
  }
  return Parser::Node();
}

std::vector<CSGTree> parse_operand_sequence(
    const Parser::Node & operandNodes,
    const stk::mesh::MetaData & meta,
    const stk::diag::Timer &lsTimer);

CSGTree parse_CSG_node(const Parser::Node& node, const stk::mesh::MetaData & meta, const stk::diag::Timer &lsTimer)
{
    if (!node)
    {
      stk::RuntimeDoomedAdHoc() << "Invalid csg_structure, missing node.";
      return CSGTree();
    }

    if (!node.is_map())
    {
      stk::RuntimeDoomedAdHoc() << "Invalid csg_structure, expecting map: " << node.info();
      return CSGTree();
    }

    const std::vector<std::string> uniteKeywords{"unite","union"};
    const Parser::Node uniteNodes = get_sequence_nodes(uniteKeywords, node);
    if (uniteNodes)
    {
        const std::vector<CSGTree> operands = parse_operand_sequence(uniteNodes, meta, lsTimer);
        if (operands.empty())
        {
            stk::RuntimeDoomedAdHoc() << "CSG operation 'union|unite' requires at least one operand " << node.info();
            return CSGTree();
        }
        return CSGTree::unite(operands);
    }

    const std::vector<std::string> intersectKeywords{"intersect","intersection"};
    const Parser::Node intersectNodes = get_sequence_nodes(intersectKeywords, node);
    if (intersectNodes)
    {
        const std::vector<CSGTree> operands = parse_operand_sequence(intersectNodes, meta, lsTimer);
        if (operands.empty())
        {
            stk::RuntimeDoomedAdHoc() << "CSG operation 'intersect|intersection' requires at least one operand " << node.info();
            return CSGTree();
        }
        return CSGTree::intersect(operands);
    }

    const std::vector<std::string> subtractKeywords{"subtract","difference","subtraction"};
    const Parser::Node subtractNodes = get_sequence_nodes(subtractKeywords, node);
    if (subtractNodes)
    {
        const std::vector<CSGTree> operands = parse_operand_sequence(subtractNodes, meta, lsTimer);
        if (operands.size() != 2)
        {
            stk::RuntimeDoomedAdHoc() << "CSG operation 'subtract|subtraction|difference' requires exactly two operands " << node.info();
            return CSGTree();
        }
        return CSGTree::subtract(operands[0], operands[1]);
    }

    /*
     * If no CSG operation was found, try parsing this node as a primitive surface.
     * This supports the existing structure:
     *
     *   - cuboid:
     *     center: [...]
     *     dimensions: [...]
     *     rotation: [...]
     */
    std::shared_ptr<Surface> parsedSurface = Surface_Parser::parse(node, meta, lsTimer);
    if (parsedSurface)
    {
      return CSGTree(parsedSurface);
    }

    stk::RuntimeDoomedAdHoc()
        << "Unknown CSG operation or surface definition. Expected one of: "
        << "'union', 'unite', 'intersection', 'intersect', "
        << "'difference', 'subtract', 'subtraction', "
        << "or a valid surface definition "
        << node.info();

    return CSGTree();
}

std::vector<CSGTree> parse_operand_sequence(
    const Parser::Node & operandNodes,
    const stk::mesh::MetaData & meta,
    const stk::diag::Timer &lsTimer)
{
    if (!operandNodes || !operandNodes.is_sequence())
    {
        stk::RuntimeDoomedAdHoc()
            << "Invalid CSG operand list, expected sequence: "
            << operandNodes.info();
        return {};
    }

    std::vector<CSGTree> operands;
    for (auto operandNode : operandNodes)
        operands.push_back(parse_CSG_node(operandNode, meta, lsTimer));

    return operands;
}

CSGTree parse_CSG_tree(const Parser::Node & node, const stk::mesh::MetaData & meta, const stk::diag::Timer &lsTimer)
{
  CSGTree tree;
  const Parser::Node csgTreeNode = node.get_map_if_present("csg_initialization");
  if ( csgTreeNode )
  {
    tree = parse_CSG_node(csgTreeNode, meta, lsTimer);
  }
  return tree;
}

}
}


