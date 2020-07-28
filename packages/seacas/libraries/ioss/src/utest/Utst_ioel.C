// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_CodeTypes.h>

#include <Ioss_ConcreteVariableType.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Initializer.h>
#include <Ioss_NullEntity.h>
#include <Ioss_Utils.h>
#include <Ioss_VariableType.h>

#include <fmt/ostream.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// ========================================================================
static int  test_all_elements();
static void test_aliases(const Ioss::NameList &elements);
static bool test_element(const std::string &type);
// ========================================================================

int main(int /* argc */, char *argv[])
{
  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer        initialize_topologies;

  int err_count = test_all_elements();

  // Make sure Ioss::NullEntity works.  Not used in IOSS itself,
  // but some clients use it, so need to make sure it compiles
  // correctly.
  std::unique_ptr<Ioss::NullEntity> entity{new Ioss::NullEntity()};
  fmt::print(stderr, "\nThe null entity type is '{}' and it contains '{}'\n", entity->type_string(),
             entity->contains_string());

  fmt::print(stderr, "{}\n", argv[0]);
  if (err_count == 0) {
    fmt::print(stderr, "\nSIERRA execution successful.\n");
    return EXIT_SUCCESS;
  }
  fmt::print(stderr, "\nSIERRA execution failed.\n");
  return EXIT_FAILURE;
}

int test_all_elements()
{
  int err_count = 0;

  Ioss::NameList elements;
  Ioss::ElementTopology::describe(&elements);

  for (const auto &element : elements) {
    fmt::print(stderr, "Testing element: {:<25}", element);
    bool result = test_element(element);
    if (result || element == "unknown" || element == "invalid_topology") {
      fmt::print(stderr, " OK\n");
    }
    else {
      fmt::print(stderr, "\n        element: {:<25}FAIL\n", element);
      err_count++;
    }
  }

  test_aliases(elements);

  // Check that asking for invalid element returns nullptr pointer.
  Ioss::ElementTopology *invalid = Ioss::ElementTopology::factory("Greg", true);
  if (invalid == nullptr) {
    fmt::print(stderr, "Testing request for invalid element: OK\n");
  }
  else {
    fmt::print(stderr, "Testing request for invalid element: FAIL\n");
    err_count++;
  }
  return err_count;
}

bool test_element(const std::string &type)
{
  // NOTE: For true test, should run with purify checking enabled to
  //       ensure we are not running off the end of arrays...

  bool                   result  = true;
  Ioss::ElementTopology *element = Ioss::ElementTopology::factory(type);
  if (element == nullptr) {
    fmt::print(stderr, "ERROR: Element type '{}' could not be constructed.", type);
    // Must return since we have a nullptr pointer and can't do further tests...
    return false;
  }

  // See if the name is an alias for another element (type != name())
  std::string name = element->name();
  fmt::print(stderr, "({})\t\t", name);

  // Check that name is alias for name...
  if (!element->is_alias(type)) {
    fmt::print(stderr, "\n\tName is not valid alias");
    result = false;
  }

  // Check that master element name is an alias...
  if (!element->is_alias(element->master_element_name())) {
    if (element->name() == "edge2d2" || element->name() == "edge2d3" ||
        element->name() == "edge2d4") { // kluge
      fmt::print(stderr, "\n\tMaster element name is not valid alias (ignore) ");
    }
    else {
      fmt::print(stderr, "\n\tMaster element name is not valid alias");
      result = false;
    }
  }

  // Check that the hash id method of selecting the element returns the correct element.
  unsigned int hash_val = Ioss::ElementTopology::get_unique_id(name);
  if (Ioss::ElementTopology::factory(hash_val) != element) {
    fmt::print(stderr, "\n\tElement to hash value conversion is not valid");
    result = false;
  }

  int order = element->order();

  bool homo_edges = element->edges_similar();
  bool homo_faces = element->faces_similar();

  int nn = element->number_nodes();
  if (nn <= 0) {
    fmt::print(stderr, "\n\tInvalid node count");
    result = false;
  }

  int ncn = element->number_corner_nodes();
  if (ncn <= 0 || ncn > nn) {
    fmt::print(stderr, "\n\tInvalid corner node count");
    result = false;
  }

  int ne = element->number_edges();
  if (ne < 0) {
    fmt::print(stderr, "\n\tInvalid edge count");
    result = false;
  }

  int nf = element->number_faces();
  if (nf < 0) {
    fmt::print(stderr, "\n\tInvalid face count");
    result = false;
  }

  // Verify Euler's Formula holds... V-E+F=2
  if (element->parametric_dimension() == 3) {
    int euler = ncn - ne + nf;
    if (euler != 2) {
      fmt::print(stderr, "\n\tEuler's formula violated (V-E+F=2), value = {}\n", euler);
      result = false;
    }
  }

  int nne = element->number_nodes_edge(0);
  if (nne == -1) {
    if (homo_edges) {
      fmt::print(stderr, "\n\tInconsistent edge homogeneity...\n");
      result = false;
    }
    else {
      for (int edge = 1; edge <= ne; edge++) {
        int nnei = element->number_nodes_edge(edge);
        if (nnei < 0 || nnei > nn) {
          fmt::print(stderr, "\n\tInconsistent nodes per edge...\n");
          result = false;
        }
      }
    }
  }
  else {
    if (nne < 0 || nne > nn) {
      fmt::print(stderr, "\n\tInconsistent nodes per edge...\n");
      result = false;
    }
  }

  int nnf = element->number_nodes_face(0);
  if (nnf > nn || nnf < -1) {
    fmt::print(stderr, "\n\tInvalid face node count");
    result = false;
  }

  // Check boundary and other topologies...
  if (nf > 0) {
    for (int i = 0; i <= nf; i++) {
      Ioss::ElementTopology *face = element->face_type(i);
      if (face == nullptr && i > 0) {
        fmt::print(stderr, "\n\tBad face type for face {}", i);
        result = false;
      }
      else if (face == nullptr && i == 0 && homo_faces) {
        fmt::print(stderr, "\n\tHomogenous faces, but null face_type");
        result = false;
      }
      else if (face != nullptr) {
        unsigned int nnfi = element->number_nodes_face(i);
        if (nnfi != static_cast<unsigned int>(face->number_nodes())) {
          fmt::print(stderr, "\n\tNode count mismatch on face {}", i);
          result = false;
        }
        if (i != 0) {
          std::vector<int> conn = element->face_connectivity(i);
          if (nnfi != conn.size()) {
            fmt::print(stderr,
                       "\n\tNode count and face connectivity size "
                       "mismatch on face {}",
                       i);
            result = false;
          }
        }
      }
    }
  }
  // Edges...
  if (ne > 0) {
    for (int i = 0; i <= ne; i++) {
      Ioss::ElementTopology *edge = element->edge_type(i);
      if (edge == nullptr && i > 0) {
        fmt::print(stderr, "\n\tBad edge type for edge {}", i);
        result = false;
      }
      else if (edge == nullptr && i == 0 && homo_edges) {
        fmt::print(stderr, "\n\tHomogenous edges, but null edge_type");
        result = false;
      }
      else if (edge != nullptr) {
        unsigned int nnei = element->number_nodes_edge(i);
        if (nnei != static_cast<unsigned int>(edge->number_nodes())) {
          fmt::print(stderr, "\n\tNode count mismatch on edge {}", i);
          result = false;
        }
        if (i != 0) {
          std::vector<int> conn = element->edge_connectivity(i);
          if (nnei != conn.size()) {
            fmt::print(stderr,
                       "\n\tNode count and edge connectivity size "
                       "mismatch on edge {}",
                       i);
            result = false;
          }
        }
      }
    }
  }

  // Variable types...
  const Ioss::VariableType *vt = Ioss::VariableType::factory(element->name());
  if (vt == nullptr) {
    fmt::print(stderr, "\n\tVariable Type does not exist for this name");
    result = false;
  }
  else {
    // See if component counts match...
    int vt_comp = vt->component_count();
    if (nn != vt_comp) {
      fmt::print(stderr, "\n\tNode count does not match component count");
      result = false;
    }
  }

  // For elements with dimension == 3
  // Get face-edge-order
  // Get face-node-connectivity
  // Foreach edge in face, get nodal connectivity
  //   ensure that node is in face connectivity...
  //
  // For an edge:   1------3------2
  //
  // For a face: Corner nodes are first in connectivity
  //             Center nodes follow.
  // So:
  //
  //                        2           x
  //    3----6----2        / \          x
  //    |         |       /   \         x
  //    7         5      5     4        x
  //    |         |     /       \       x
  //    0----4----1    0----3----1      x
  //
  if (element->parametric_dimension() == 3) {
    for (int i = 1; i <= nf; i++) {
      unsigned int fncn           = element->face_type(i)->number_corner_nodes();
      unsigned int num_edges_face = element->number_edges_face(i);
      if (fncn != num_edges_face) {
        fmt::print(stderr, "\n\tFace corner node count should match edges/face for face {}", i);
        result = false;
      }

      // Nodes defining face...
      std::vector<int> face_conn = element->face_connectivity(i);

      // Edges defining face...
      std::vector<int> face_edge_conn = element->face_edge_connectivity(i);
      if (num_edges_face != face_edge_conn.size()) {
        fmt::print(stderr, "\n\tEdges per face mismatch for face {}", i);
        result = false;
      }
      else {
        for (unsigned int j = 0; j < num_edges_face; j++) {
          // Not implemented in all elements yet...
          std::vector<int> edge_conn = element->edge_connectivity(face_edge_conn[j] + 1);
          // Check that first two nodes in 'edge_conn' match
          // corresponding nodes in 'face_conn'
          if ((edge_conn[0] != face_conn[j] && edge_conn[1] != face_conn[j]) ||
              (edge_conn[0] != face_conn[(j + 1) % fncn] &&
               edge_conn[1] != face_conn[(j + 1) % fncn])) {
            fmt::print(stderr,
                       "\n\tEdge Connectivity does not match face connectivity for edge {} on face "
                       "{}\nEdge: {} {}\nFace: {} {}\n",
                       j + 1, i, edge_conn[0], edge_conn[1], face_conn[j],
                       face_conn[(j + 1) % fncn]);

            result = false;
          }
          if (element->name() == "wedge12" || element->name() == "hex16") {
            continue;
          }

          if (edge_conn.size() != size_t(order) + 1) {
            fmt::print(stderr, "\n\tInvalid edge connectivity count. ({} must equal {})",
                       edge_conn.size(), order + 1);
            result = false;
          }

          if (order > 1 && face_conn.size() < fncn + (order - 1) * num_edges_face) {
            fmt::print(stderr,
                       "\n\tInvalid face connectivity count ({} must be greater than {} + {}*{}).",
                       face_conn.size(), fncn, order - 1, num_edges_face);
            result = false;
          }
          if (order == 2) {
            if (edge_conn[2] != face_conn[fncn + j]) {
              fmt::print(stderr,
                         "\n\tMid-Side Node Edge Connectivity does not match face "
                         "connectivity for edge {} on face {}",
                         j + 1, i);
              result = false;
            }
          }
          else if (order == 3) {
            auto en2 = edge_conn[2];
            auto en3 = edge_conn[3];

            auto fn2 = face_conn[fncn + 2 * j];
            auto fn3 = face_conn[fncn + 2 * j + 1];
            if (en2 != fn2 && en3 != fn2 && en2 != fn3 && en3 != fn3) {
              fmt::print(stderr,
                         "\n\tMid-Side Node Edge Connectivity does not match face "
                         "connectivity for edge {} on face {} ({} {} : {} {})",
                         j + 1, i, en2, en3, fn2, fn3);
              result = false;
            }
          }
        }
      }
    }
  }
  return result;
}

void test_aliases(const Ioss::NameList &elements)
{
  int count = elements.size();
  fmt::print(stderr, "\n\nTesting Element Topology Aliases...\n");

  for (int i = 0; i < count; i++) {
    Ioss::ElementTopology *el_top = Ioss::ElementTopology::factory(elements[i]);
    if (el_top->name() == elements[i]) {
      fmt::print(stderr, "Element: {:<25}({}) has the following aliases:\n\t", elements[i],
                 el_top->master_element_name());
      for (int j = 0; j < count; j++) {
        if (i != j && el_top->is_alias(elements[j])) {
          fmt::print(stderr, "{}, ", elements[j]);
        }
      }
      fmt::print(stderr, "\n");
    }
  }
}
