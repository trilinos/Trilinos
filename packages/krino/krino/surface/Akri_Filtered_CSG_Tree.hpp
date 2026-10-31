// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SURFACE_AKRI_FILTERED_CSG_TREE_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FILTERED_CSG_TREE_HPP_
#include <Akri_CSG_Tree.hpp>

namespace krino {

CSGTree filtered_CSG_tree(
    const CSGTree& tree,
    const std::vector<std::shared_ptr<Surface>>& nearbySurfaces);

}

#endif /* KRINO_KRINO_SURFACE_AKRI_FILTERED_CSG_TREE_HPP_ */
