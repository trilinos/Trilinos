/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Initializer.h>
#include <Ioss_VariableType.h>
#include <Ioss_Utils.h>

#include <Ioss_ElementTopology.h>

#include <stk_io/IossBridge.hpp>
#include <Shards_CellTopology.hpp>

#include <assert.h>

namespace {
template<typename T, typename U>
int my_assert(T a, T b, U msg) {
  if (a != b) {
    std::cerr << "\tERROR: '" << msg << "' assertion failed: " << a << " is not equal to " << b << "\n";
    return 1;
  }
  return 0;
}
}

namespace {
int testElement(const std::string &name)
{
  int errors = 0;
  Ioss::ElementTopology *element = Ioss::ElementTopology::factory(name);
  if (element == NULL) {
    std::cerr << "\tERROR: Element type '" << name << "' could not be constructed.";
    // Must return since we have a NULL pointer and can't do further tests...
    return 1;
  }
  if (element->name() != name) {
    // This is an alias.  Don't run it through the rest of the tests
    std::cerr << "Element '" << name << "' is an alias for element '" << element->name() << "'\n";
    return 0;
  }

  std::cerr << "Testing element '" << name << "'\n";
  // Currently not supported in shards:
  if (element->name() == "unknown") {
    std::cerr << "\tERROR (EXPECTED): No support for '" << element->name() << "'\n";
    return 0;
  }

  // Get the corresponding shards CellTopologyData* ..
  const CellTopologyData *cell_data = stk::io::map_topology_ioss_to_cell(element);
  if (cell_data == NULL) {
    std::cerr << "\tERROR: Could not find a shards CellTopology corresponding to the Ioss::ElementTopology element '"
              << name << "'.";
    // Must return since we have a NULL pointer and can't do further tests...
    return 1;
  }

  // See if we get the same element back when converting from
  // CellTopologyData to Ioss::ElementToplogy
  std::string new_name = stk::io::map_topology_cell_to_ioss(cell_data, element->spatial_dimension());
  Ioss::ElementTopology *new_element = Ioss::ElementTopology::factory(new_name);
  if (element->name() != new_element->name()) {
    std::cerr << "\tERROR: New name = '" << new_element->name()
              << "' doesn't match old name '" << element->name()
              << "'\n";
    errors++;
  }

  shards::CellTopology cell(cell_data);

  // At this point, 'element' is the Ioss element topology and
  //                'cell' is the corresponding shards CellTopology data pointer.
  // Make sure that they agree on all subcell details...
  // Exceptions:
  // 1. An Ioss Node has 1 node per element; a shards Node has 0 nodes per element...

  errors += my_assert(cell.getNodeCount(),
                      static_cast<unsigned>(element->number_nodes()),
                      "node count");
  errors += my_assert(cell.getVertexCount(),
                      static_cast<unsigned>(element->number_corner_nodes()),
                      "vertex count");

  // NOTE: CellTopology and Ioss disagree on parametric dimension.
  int add_to = element->spatial_dimension() != element->parametric_dimension() && element->is_element() ? 1 : 0;
  errors += my_assert(cell.getDimension(),
                      static_cast<unsigned>(element->parametric_dimension()+add_to),
                      "parametric dimension");
  errors += my_assert(cell.getEdgeCount(),
                      static_cast<unsigned>(element->number_edges()),
                      "edge count");

  // NOTE: Ioss counts edges and faces as boundaries for shell elements
  int add_boundary = 0;
  if (add_to == 1 && element->spatial_dimension() == 3 && element->parametric_dimension() == 2)
    add_boundary = cell.getEdgeCount();

  if (element->name() == "edge2" || element->name() == "edge3")
    add_boundary += 2;
  
  errors += my_assert(cell.getSideCount() + add_boundary,
                      static_cast<unsigned>(element->number_boundaries()),
                      "boundary count");


  // Check face topologies for all elements...
  if (element->is_element()) {
    if (cell.getDimension() == 3) {
      int face_count = element->number_faces();
      for (int i=0; i < face_count; i++) {
        Ioss::ElementTopology *face = element->face_type(i+1);
        const CellTopologyData *cell_face = cell.getCellTopologyData(cell.getDimension()-1,i);
        errors += my_assert(face->name(),
                            stk::io::map_topology_cell_to_ioss(cell_face,face->spatial_dimension()),
                            "face type");

        Ioss::IntVector fcon = element->face_connectivity(i+1);
        size_t node_count = fcon.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "face node connectivity for node " << j << " on face " << i;
          errors += my_assert(fcon[j],
                              static_cast<int>(cell.getNodeMap(cell.getDimension()-1, i, j)),
                              msg.str());
        }
      }

      int edge_count = element->number_edges();
      for (int i=0; i < edge_count; i++) {

        Ioss::IntVector fcon = element->edge_connectivity(i+1);
        size_t node_count = fcon.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "edge node connectivity for node " << j << " on edge " << i;
          errors += my_assert(fcon[j],
                              static_cast<int>(cell.getNodeMap(cell.getDimension()-2, i, j)),
                              msg.str());
        }
      }
    }
    else if (cell.getDimension() == 2) {
      int edge_count = element->number_edges();
      for (int i=0; i < edge_count; i++) {
        Ioss::ElementTopology *edge = element->edge_type(i+1);
        const CellTopologyData *cell_edge = cell.getCellTopologyData(cell.getDimension()-1,i);
        errors += my_assert(edge->name(),
                            stk::io::map_topology_cell_to_ioss(cell_edge, edge->spatial_dimension()),
                            "edge type");

        Ioss::IntVector econ = element->edge_connectivity(i+1);
        size_t node_count = econ.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "edge node connectivity for node " << j << " on edge " << i;
          errors += my_assert(econ[j],
                              static_cast<int>(cell.getNodeMap(cell.getDimension()-1, i, j)),
                              msg.str());
        }
      }

    }
  }
  return errors;
}
}

STKUNIT_UNIT_TEST(UnitTestTopology, testUnit)
{
  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer        initialize_topologies;

  Ioss::NameList elements;
  int element_count = Ioss::ElementTopology::describe(&elements);

  int errors = 0;
  for (int i=0; i < element_count; i++) {
    // FIXME: Need to totally skip tetra7 for now
    if (elements[i] == "tetra7" ||
	elements[i] == "tetra11" ||
	elements[i] == "trishell4" ||
	elements[i] == "trishell7") {
      continue;
    }

    int current_error = testElement(elements[i]);
    if (elements[i] != "node" &&
        elements[i] != "rod3d2" &&
        elements[i] != "rod3d3" &&
        elements[i] != "tri4a" /*FIXME?*/) {
      errors += current_error;
    }
    else {
      if (current_error > 0)
        std::cerr << "\t\tIGNORING " << elements[i] << " ERRORS...\n";
    }
  }
  STKUNIT_ASSERT(errors == 0);
}


