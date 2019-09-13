// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef RebalanceMesh_hpp
#define RebalanceMesh_hpp

#include <percept/PerceptMesh.hpp>

namespace percept {

  class RebalanceMesh {
  public:
    RebalanceMesh(PerceptMesh& eMesh, ScalarFieldType *weight_field=0, bool debug=false, bool doAdaptedMesh = true );
    // return the current imbalance factor
    double compute_imbalance();
    // return the imbalance factor after rebalancing
    double rebalance();

    static
    void fixup_family_trees(PerceptMesh& eMesh, bool debug);

  protected:
    void build_entities_to_rebalance_list(stk::mesh::EntityVector& vec);
    void set_weights(const std::vector<stk::mesh::Entity>& entities_to_rebalance);

    PerceptMesh& m_eMesh;
    bool m_debug;
    ScalarFieldType *m_weight_field;
    stk::mesh::EntityRank m_rank_to_rebalance;
    bool m_doAdaptedMesh;
  };


}

#endif
