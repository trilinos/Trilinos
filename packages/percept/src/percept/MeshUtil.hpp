// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_MeshUtil_hpp
#define percept_MeshUtil_hpp


#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/PerceptMesh.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>


  namespace percept {

    class MeshUtil
    {
    public:
      static bool m_debug;

      // high level methods
      static bool facesConsistent(PerceptMesh& eMesh);

      static void checkTopology(percept::PerceptMesh& eMesh);

      // low-level support methods

      static void fillSideNodes(percept::PerceptMesh& eMesh, stk::mesh::Entity element, unsigned iside, std::vector<stk::mesh::EntityId>& side_nodes);

      static void fillSideNodes(percept::PerceptMesh& eMesh, stk::mesh::Entity element, unsigned iside, std::vector<stk::mesh::Entity>& side_nodes);

      static double triFaceArea(percept::PerceptMesh& eMesh, stk::mesh::Entity element, unsigned iside);

      static bool nodesMatch(  std::vector<stk::mesh::EntityId>& side1, std::vector<stk::mesh::EntityId>& side2, bool reverse=false);

      static bool sharesFace(percept::PerceptMesh& eMesh, stk::mesh::Entity element1, stk::mesh::Entity element2, unsigned& iside1, unsigned& iside2);

      static bool facesConsistent1(percept::PerceptMesh& eMesh, stk::mesh::Entity element1, stk::mesh::Entity element2);

    };
  }

#endif

