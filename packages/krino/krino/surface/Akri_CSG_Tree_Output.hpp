// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SURFACE_AKRI_CSG_TREE_OUTPUT_HPP_
#define KRINO_KRINO_SURFACE_AKRI_CSG_TREE_OUTPUT_HPP_
#include <iostream>
#include <memory>

namespace krino {

class CSGTreeNode;

void output_tree_node(
    std::ostream& os,
    const std::shared_ptr<CSGTreeNode>& node,
    unsigned indent);

std::string output_CSG_surface_tree(const CSGTree & tree);

}

#endif /* KRINO_KRINO_SURFACE_AKRI_CSG_TREE_OUTPUT_HPP_ */
