#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "gtest/gtest.h"

namespace {

using namespace stk::middle_mesh::mesh;
using namespace stk::middle_mesh::mesh::impl;
using namespace stk::middle_mesh::predicates::impl;

void check_all_entities_same(const std::vector<PointRecord>& records)
{
  assert(records.size() > 0);
  EXPECT_NE(records[0].type, PointClassification::Interior);
  EXPECT_NE(records[0].type, PointClassification::Exterior);

  MeshEntityPtr entity = get_entity(records[0]);

  for (auto& r : records)
  {
    MeshEntityPtr entity2 = get_entity(r);
    EXPECT_EQ(entity->get_id(), entity2->get_id());
    EXPECT_EQ(entity->get_type(), entity2->get_type());
  }
}

void check_records_uniqueness(const std::vector<PointRecord>& records)
{
  EXPECT_GE(records.size(), 1u);

  // test all records types are the same
  PointClassification type = records[0].type;
  for (auto& r : records)
    EXPECT_EQ(r.type, type);

  if (type == PointClassification::Interior)
  {
    EXPECT_EQ(records.size(), 1u);
  } else if (type == PointClassification::Edge)
  {
    MeshEntityPtr edge = get_entity(records[0]);
    int numelExpected  = std::min(edge->count_up(), int(2));
    EXPECT_EQ((int)records.size(), numelExpected);
    check_all_entities_same(records);
  } else if (type == PointClassification::Vert)
  {
    MeshEntityPtr vert = get_entity(records[0]);
    std::vector<MeshEntityPtr> els;
    get_upward(vert, 2, els);
    int numelExpected = std::min(els.size(), size_t(4));
    EXPECT_EQ((int)records.size(), numelExpected);
    check_all_entities_same(records);
  }
}

} // namespace

TEST(PointClassifierNormalWrapper, EigthSphere)
{
  if (stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // when projecting from a sphere to another sphere with larger
  // radius, each point should project onto exactly 1 point on
  // the larger sphere
  std::shared_ptr<Mesh> mesh1 = create_eigth_sphere(5, 5, 0.5, 1.5);

  std::shared_ptr<Mesh> mesh2 = create_eigth_sphere(5, 5, 0.5, 2);

  PointClassifierNormalWrapper c(mesh1);

  std::vector<MeshEntityPtr> elements;
  for (auto& v : mesh1->get_vertices())
    if (v)
    {
      get_upward(v, 2, elements);
      MeshEntityPtr el1 = elements[0];

      std::vector<PointRecord> records;
      for (auto& el2 : mesh2->get_elements())
        if (el2)
        {
          PointRecord r = c.classify(el2, el1, v->get_point_orig(0));
          if (r.type != PointClassification::Exterior)
          {
            records.push_back(r);
          }
        }

      check_records_uniqueness(records);
    }
}
