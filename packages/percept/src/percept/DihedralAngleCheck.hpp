// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_DihedralAngleCheck_hpp
#define percept_DihedralAngleCheck_hpp

#include <percept/PerceptMesh.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <unordered_map>

namespace percept {

  class DihedralAngleCheck
  {

  public:

    DihedralAngleCheck(PerceptMesh * eMesh, int print_level=0) : m_eMesh(eMesh), m_needs_delete(false), m_print_level(print_level)  {}

    DihedralAngleCheck(stk::mesh::BulkData *bulk, int print_level=0) : m_eMesh(0), m_needs_delete(true), m_print_level(print_level) {
      m_eMesh = new PerceptMesh(&bulk->mesh_meta_data(), bulk, true);
    }

    ~DihedralAngleCheck() {
      if (m_needs_delete)
        delete m_eMesh;
    }

    // modified from Victor Brunini's Aria code
    using ObtuseAngleSidePair = std::pair<int, int>; // Pair of side ordinals
    using EntityObtuseAngleSidePairs = std::pair<stk::mesh::Entity, std::vector<ObtuseAngleSidePair> >;
    using EntityObtuseAngleSidePairsMap = std::unordered_map<stk::mesh::Entity, std::vector<ObtuseAngleSidePair>, std::hash<stk::mesh::Entity> >;

    void find_simplex_elements_with_obtuse_angles(std::ostream& out = std::cout);

    const EntityObtuseAngleSidePairsMap& get_map() { return m_map; }

    PerceptMesh *get_mesh() { return m_eMesh; }

  private:

    PerceptMesh *m_eMesh;
    bool m_needs_delete;
    int m_print_level;
    EntityObtuseAngleSidePairsMap m_map;

  };



}

#endif
