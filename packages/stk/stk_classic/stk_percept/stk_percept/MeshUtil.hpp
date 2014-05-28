#ifndef stk_percept_MeshUtil_hpp
#define stk_percept_MeshUtil_hpp

/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>


namespace stk {
  namespace percept {

    class MeshUtil
    {
    public:
      static bool m_debug;

      // high level methods
      static bool facesConsistent(PerceptMesh& eMesh);

      static void checkTopology(percept::PerceptMesh& eMesh);

      // low-level support methods
      
      static void fillSideNodes(stk::mesh::Entity& element, unsigned iside, std::vector<stk::mesh::EntityId>& side_nodes);

      static void fillSideNodes(stk::mesh::Entity& element, unsigned iside, std::vector<stk::mesh::Entity *>& side_nodes);

      static double triFaceArea(percept::PerceptMesh& eMesh, stk::mesh::Entity& element, unsigned iside);

      static bool nodesMatch(  std::vector<stk::mesh::EntityId>& side1, std::vector<stk::mesh::EntityId>& side2, bool reverse=false);

      static bool sharesFace(stk::mesh::Entity& element1, stk::mesh::Entity& element2, unsigned& iside1, unsigned& iside2);

      static bool facesConsistent1(percept::PerceptMesh& eMesh, stk::mesh::Entity& element1, stk::mesh::Entity& element2);

    }; 
  }
}

#endif

