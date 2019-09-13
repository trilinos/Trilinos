// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SpacingFieldUtil_hpp
#define SpacingFieldUtil_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/PerceptMesh.hpp>

  namespace percept {

    class SpacingFieldUtil
    {
    public:
      /// either average or take max of spacing to the nodes
      enum SpacingType { SPACING_AVE, SPACING_MAX };

      SpacingFieldUtil(PerceptMesh& eMesh, SpacingType type=SPACING_AVE) : m_eMesh(eMesh), m_type(type) {}

      /// computes an approximation of the spacing field in x,y,z directions - not rotationally invariant
      void compute_spacing_field();

      /// find spacing in unit vector dir direction (local computation - is rotationally invariant)
      double spacing_at_node_in_direction(const double dir[3], stk::mesh::Entity node, stk::mesh::Selector *element_selector=0);

    private:
      PerceptMesh& m_eMesh;
      const SpacingType m_type;
    };
  }

#endif
#endif
