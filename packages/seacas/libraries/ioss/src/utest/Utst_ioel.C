// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_CodeTypes.h"
#include "Ioss_ConcreteVariableType.h"
#include "Ioss_ElementPermutation.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_Initializer.h"
#include "Ioss_NullEntity.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fmt/core.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// ========================================================================
namespace {
  void test_aliases(const Ioss::NameList &elements);
  void test_element(const std::string &type);
} // namespace
// ========================================================================

int main(IOSS_MAYBE_UNUSED int argc, char **argv)
{
  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer        initialize_topologies;

  Catch::Session session; // There must be exactly one instance
  int            returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) // Indicates a command line error
    return returnCode;

  auto ok = session.run() == 0;

  if (ok) {
    fmt::print(stderr, "\nSIERRA execution successful.\n");
    return EXIT_SUCCESS;
  }
  fmt::print(stderr, "\nSIERRA execution failed.\n");
  return EXIT_FAILURE;
}

TEST_CASE("nullentity")
{
  // Make sure Ioss::NullEntity works.  Not used in IOSS itself,
  // but some clients use it, so need to make sure it compiles
  // correctly.
  auto entity = std::make_unique<Ioss::NullEntity>();
  REQUIRE(entity);
  fmt::print(stderr, "\nThe null entity type is '{}' and it contains '{}'\n", entity->type_string(),
             entity->contains_string());
}

TEST_CASE("invalid element")
{
  // Check that asking for invalid element returns nullptr pointer.
  Ioss::ElementTopology *invalid = Ioss::ElementTopology::factory("Greg", true);
  REQUIRE(invalid == nullptr);
}

TEST_CASE("test_all_elements")
{
  Ioss::NameList elements = Ioss::ElementTopology::describe();

  for (const auto &element : elements) {
    DYNAMIC_SECTION(element) { test_element(element); }
  }
}

TEST_CASE("test aliases")
{
  Ioss::NameList elements = Ioss::ElementTopology::describe();
  test_aliases(elements);
}

namespace {
  void test_element(const std::string &type)
  {
    // NOTE: For true test, should run with purify checking enabled to
    //       ensure we are not running off the end of arrays...

    Ioss::ElementTopology *element = Ioss::ElementTopology::factory(type);
    REQUIRE(element != nullptr);

    // See if the name is an alias for another element (type != name())

    // Check that name is alias for name...
    CHECK(element->is_alias(type));

    // Output the 'shape'
    auto shape = Ioss::Utils::lowercase(Ioss::Utils::shape_to_string(element->shape()));

    std::string name = element->name();
    fmt::print(stderr, "\n{:<25}({})\t[{}]", type, name, shape);

    // Check that master element name is an alias...
    if (element->name() == "edge2d2" || element->name() == "edge2d3" ||
        element->name() == "edge2d4") { // kluge
      CHECK_NOFAIL(element->is_alias(element->master_element_name()));
    }
    else {
      CHECK(element->is_alias(element->master_element_name()));
    }

    if (type == "unknown" || type == "invalid_topology") {
      return;
    }

    // Check that the hash id method of selecting the element returns the correct element.
    unsigned int hash_val = Ioss::ElementTopology::get_unique_id(name);
    CHECK(Ioss::ElementTopology::factory(hash_val) == element);

    int order = element->order();

    bool homo_edges = element->edges_similar();
    bool homo_faces = element->faces_similar();

    int nn = element->number_nodes();
    CHECK(nn > 0);

    int ncn = element->number_corner_nodes();
    CHECK(ncn > 0);
    CHECK(ncn <= nn);

    auto *perm = element->permutation();
    REQUIRE(perm != nullptr);

    CHECKED_IF(element->name() != "node")
    {
      REQUIRE(static_cast<int>(perm->num_permutation_nodes()) == element->number_corner_nodes());
      // Element shape matches permutation type
      REQUIRE(Ioss::Utils::lowercase(Ioss::Utils::shape_to_string(element->shape())) ==
              perm->type());

      CHECKED_IF(perm->num_permutations() > 0)
      {
        const auto &permutation = perm->permutation_indices(0);
        const auto  connect     = element->element_connectivity();
        for (size_t i = 0; i < permutation.size(); i++) {
          // permutation matches element connectivity
          CHECK(permutation[i] == connect[i]);
        }
      }
    }

    int ne = element->number_edges();
    CHECK(ne >= 0);

    int nf = element->number_faces();
    CHECK(nf >= 0);

    // Verify Euler's Formula holds... V-E+F=2
    if (element->parametric_dimension() == 3) {
      int euler = ncn - ne + nf;
      CHECK(euler == 2);
    }

    int nne = element->number_nodes_edge(0);
    if (nne == -1) {
      CHECK(homo_edges == false);
      for (int edge = 1; edge <= ne; edge++) {
        int nnei = element->number_nodes_edge(edge);
        CHECK(nnei >= 0);
        CHECK(nnei <= nn);
      }
    }
    else {
      CHECK(nne >= 0);
      CHECK(nne <= nn);
    }

    int nnf = element->number_nodes_face(0);
    CHECK(nnf <= nn);
    CHECK(nnf >= -1);

    // Check boundary and other topologies...
    if (nf > 0) {
      Ioss::ElementTopology *face0 = element->face_type(0);
      if (homo_faces == false) {
        CHECK(face0 == nullptr);
      }
      else {
        CHECK(face0 != nullptr);
      }

      for (int i = 1; i <= nf; i++) {
        Ioss::ElementTopology *face = element->face_type(i);
        CHECK(face != nullptr);
        CHECKED_IF(face != nullptr)
        {
          unsigned int nnfi = element->number_nodes_face(i);
          CHECK(nnfi == static_cast<unsigned int>(face->number_nodes()));
          std::vector<int> conn = element->face_connectivity(i);
          CHECK(nnfi == conn.size());
        }
      }
    }
    // Edges...
    if (ne > 0) {
      Ioss::ElementTopology *edge0 = element->edge_type(0);
      if (homo_edges == false) {
        CHECK(edge0 == nullptr);
      }
      else {
        CHECK(edge0 != nullptr);
      }

      for (int i = 1; i <= ne; i++) {
        Ioss::ElementTopology *edge = element->edge_type(i);
        CHECK(edge != nullptr);
        CHECKED_IF(edge != nullptr)
        {
          unsigned int nnei = element->number_nodes_edge(i);
          CHECK(nnei == static_cast<unsigned int>(edge->number_nodes()));
          std::vector<int> conn = element->edge_connectivity(i);
          CHECK(nnei == conn.size());
        }
      }
    }

    // Variable types...
    const Ioss::VariableType *vt = Ioss::VariableType::factory(element->name());
    REQUIRE(vt != nullptr);
    CHECKED_IF(vt != nullptr)
    {
      // See if component counts match...
      int vt_comp = vt->component_count();
      CHECK(nn == vt_comp);
    }

    // For elements with parametric dimension == 3
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
        CHECK(fncn == num_edges_face);

        // Nodes defining face...
        std::vector<int> face_conn = element->face_connectivity(i);

        // Edges defining face...
        std::vector<int> face_edge_conn = element->face_edge_connectivity(i);
        REQUIRE(num_edges_face == face_edge_conn.size());
        CHECKED_IF(num_edges_face == face_edge_conn.size())
        {
          for (unsigned int j = 0; j < num_edges_face; j++) {
            // Not implemented in all elements yet...
            std::vector<int> edge_conn = element->edge_connectivity(face_edge_conn[j] + 1);
            // Check that first two nodes in 'edge_conn' match
            // corresponding nodes in 'face_conn'
            bool ok = (edge_conn[0] == face_conn[j] || edge_conn[1] == face_conn[j]) &&
                      (edge_conn[0] == face_conn[(j + 1) % fncn] ||
                       edge_conn[1] == face_conn[(j + 1) % fncn]);
            CHECK(ok);

            if (element->name() == "wedge12" || element->name() == "hex16") {
              continue;
            }

            CHECK(edge_conn.size() == size_t(order) + 1);

            CHECK(order >= 1);
            CHECK(face_conn.size() >= fncn + (order - 1) * num_edges_face);
            CHECKED_IF(order == 2) { CHECK(edge_conn[2] == face_conn[fncn + j]); }
            CHECKED_IF(order == 3)
            {
              auto en2 = edge_conn[2];
              auto en3 = edge_conn[3];

              auto fn2 = face_conn[fncn + 2 * j];
              auto fn3 = face_conn[fncn + 2 * j + 1];
              // Mid-Side Node Edge Connectivity must match face connectivity
              CHECK((en2 == fn2 || en3 == fn2 || en2 == fn3 || en3 == fn3));
            }
          }
        }
      }
    }
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
} // namespace
