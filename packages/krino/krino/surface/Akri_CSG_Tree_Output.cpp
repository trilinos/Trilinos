// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_CSG_Tree.hpp>
#include <Akri_Surface.hpp>

namespace krino {

void output_tree_node(
    std::ostream& os,
    const std::shared_ptr<CSGTreeNode>& node,
    unsigned indent)
{
    os << std::string(indent, ' ');

    if (!node)
    {
        os << "<null node>\n";
        return;
    }

    if (node->node_operator() == CSGOperation::Leaf)
    {
        os << "Surface ";
        if (node->surface())
            os << node->surface()->output();
        else
            os << "ptr=<null>";
        os << "\n";
        return;
    }

    os << CSGTree::operation_name(node->node_operator()) << "\n";

    os << std::string(indent+2, ' ') << "operands:\n";
    for (const auto & operand : node->vector_operands())
      output_tree_node(os, operand, indent + 4);
}

std::string output_CSG_surface_tree(const CSGTree & tree)
{
    if (!tree.root())
        return "CSG_Surface with empty tree";

    std::ostringstream os;
    os << "CSG_Surface ";
    output_tree_node(os, tree.root(), 0);
    return os.str();
}

}



