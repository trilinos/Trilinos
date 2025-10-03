#include <gtest/gtest.h>
#include "Akri_FieldPenaltySensitivities.hpp"
#include "Akri_AuxMetaData.hpp"
#include "Akri_BoundingBoxMesh.hpp"
#include "Akri_MeshHelpers.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/parallel/ParallelVectorConcat.hpp"

namespace
{

void populate_bbox_mesh_with_aura_and_activate(krino::BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, 
  const stk::math::Vector3d & maxCorner, double size)
{
  bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), size);
  bboxMesh.populate_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  krino::activate_all_entities(mesh, krino::AuxMetaData::get(bboxMesh.meta_data()).active_part());
  krino::populate_stk_local_ids(mesh);
}
}

TEST(SignedDistancePenalty, ComputeL2Penalty)
{
  double meshSize = 0.1;
  const double normalization = 4.;
  std::vector<stk::topology> topologies{stk::topology::TRIANGLE_3_2D, stk::topology::TETRAHEDRON_4};
  
  for(auto && topo : topologies)
  {
    const double expVal = topo.dimension() == 3 ? 8. : 8./3.;
    std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(topo, MPI_COMM_WORLD);
    auto & auxMeta = krino::AuxMetaData::get_or_create(bboxMesh->meta_data());
    auto ls1 = auxMeta.register_field("LS1", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());
    auto ls2 = auxMeta.register_field("LS2", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());

    const stk::math::Vector3d minCorner{-1., -1., -1.}; 
    const stk::math::Vector3d maxCorner{1., 1., 1.};
    double vol = 1.;
    for(unsigned d=0; d<topo.dimension(); d++ ) vol *= (maxCorner[d]-minCorner[d]);
    populate_bbox_mesh_with_aura_and_activate(*bboxMesh, minCorner, maxCorner, meshSize);

    for(auto && b : bboxMesh->bulk_data().get_buckets(stk::topology::NODE_RANK, bboxMesh->meta_data().universal_part()))
    {
      for(auto && n : (*b))
      {
        auto * lsVal1 = krino::field_data<double>(ls1, n);
        auto * lsVal2 = krino::field_data<double>(ls2, n);
        lsVal1[0] = 0.;
        lsVal2[0] = 0.;
        auto coords = krino::get_vector_field(bboxMesh->bulk_data(), krino::FieldRef(bboxMesh->meta_data().coordinate_field()), n, topo.dimension());
        for(unsigned d=0; d<topo.dimension(); d++)
        {
          lsVal1[0] += coords[d];
          lsVal2[0] += 2.*coords[d];
        }
      }
    }

    std::map<stk::mesh::EntityId, double> sens;
    double l2pen = krino::compute_l2_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
      normalization, bboxMesh->meta_data().universal_part(), sens);

    EXPECT_NEAR(l2pen, expVal/normalization/normalization/vol, 1e-8);
  }
}

TEST(SignedDistancePenalty, ComputeGradPenalty)
{
  double meshSize = 0.1;
  std::vector<stk::topology> topologies{stk::topology::TRIANGLE_3_2D, stk::topology::TETRAHEDRON_4};
  
  for(auto && topo : topologies)
  {
    const double expVal = topo.dimension() == 3 ? 32. : 32./3.;
    std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(topo, MPI_COMM_WORLD);
    auto & auxMeta = krino::AuxMetaData::get_or_create(bboxMesh->meta_data());
    auto ls1 = auxMeta.register_field("LS1", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());
    auto ls2 = auxMeta.register_field("LS2", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());

    const stk::math::Vector3d minCorner{-1., -1., -1.}; 
    const stk::math::Vector3d maxCorner{1., 1., 1.};
    double vol = 1.;
    for(unsigned d=0; d<topo.dimension(); d++ ) vol *= (maxCorner[d]-minCorner[d]);
    populate_bbox_mesh_with_aura_and_activate(*bboxMesh, minCorner, maxCorner, meshSize);

    for(auto && b : bboxMesh->bulk_data().get_buckets(stk::topology::NODE_RANK, bboxMesh->meta_data().universal_part()))
    {
      for(auto && n : (*b))
      {
        auto * lsVal1 = krino::field_data<double>(ls1, n);
        auto * lsVal2 = krino::field_data<double>(ls2, n);
        lsVal1[0] = 0.;
        lsVal2[0] = 0.;
        auto coords = krino::get_vector_field(bboxMesh->bulk_data(), krino::FieldRef(bboxMesh->meta_data().coordinate_field()), n, topo.dimension());
        for(unsigned d=0; d<topo.dimension(); d++)
        {
          lsVal1[0] += coords[d]*coords[d];
          lsVal2[0] += 2.*coords[d]*coords[d];
        }
      }
    }

    std::map<stk::mesh::EntityId, double> sens;
    double gradPen = krino::compute_gradient_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
      bboxMesh->meta_data().universal_part(), sens);

    EXPECT_NEAR(gradPen, expVal/vol, 1e-2*expVal/vol);
  }
}

TEST(SignedDistancePenalty, L2PenaltySens)
{
  double meshSize = 1.0;
  const double normalization = 4.;
  std::vector<stk::topology> topologies{stk::topology::TRIANGLE_3_2D, stk::topology::TETRAHEDRON_4};
  
  for(auto && topo : topologies)
  {
    std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(topo, MPI_COMM_WORLD);
    auto & auxMeta = krino::AuxMetaData::get_or_create(bboxMesh->meta_data());
    auto ls1 = auxMeta.register_field("LS1", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());
    auto ls2 = auxMeta.register_field("LS2", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());

    const stk::math::Vector3d minCorner{0., 0., 0.}; 
    const stk::math::Vector3d maxCorner{4., 4., 4.};
    populate_bbox_mesh_with_aura_and_activate(*bboxMesh, minCorner, maxCorner, meshSize);

    for(auto && b : bboxMesh->bulk_data().get_buckets(stk::topology::NODE_RANK, bboxMesh->meta_data().universal_part()))
    {
      for(auto && n : (*b))
      {
        auto * lsVal1 = krino::field_data<double>(ls1, n);
        auto * lsVal2 = krino::field_data<double>(ls2, n);
        lsVal1[0] = 0.;
        lsVal2[0] = 0.;
        auto coords = krino::get_vector_field(bboxMesh->bulk_data(), krino::FieldRef(bboxMesh->meta_data().coordinate_field()), n, topo.dimension());
        for(unsigned d=0; d<topo.dimension(); d++)
        {
          lsVal1[0] += coords[d];
          lsVal2[0] += 2.*coords[d];
        }
      }
    }

    std::map<stk::mesh::EntityId, double> sens;
    double l2pen = krino::compute_l2_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
      normalization, bboxMesh->meta_data().universal_part(), sens);

    std::vector<stk::mesh::EntityId> lclSensIds;
    std::vector<double> lclSensVals;
    for(auto && entry : sens)
    {
      lclSensIds.push_back(entry.first);
      lclSensVals.push_back(entry.second);
    }
    std::vector<stk::mesh::EntityId> gblSensIds;
    std::vector<double> gblSensVals;
    stk::parallel_vector_concat(bboxMesh->bulk_data().parallel(), lclSensIds, gblSensIds);
    stk::parallel_vector_concat(bboxMesh->bulk_data().parallel(), lclSensVals, gblSensVals);
    const double fdTol = 1e-4;

    for(unsigned i=0; i<gblSensIds.size(); i++)
    {
      std::map<stk::mesh::EntityId, double> dummySens;
      double origLSVal = 0;
      auto node = bboxMesh->bulk_data().get_entity(stk::topology::NODE_RANK, gblSensIds[i]);
      if(node.is_local_offset_valid())
      {
        auto * lsVal = krino::field_data<double>(ls1, node);
        origLSVal = lsVal[0];
        lsVal[0] += fdTol;
      }
      double l2penFD = krino::compute_l2_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
        normalization, bboxMesh->meta_data().universal_part(), dummySens);
      EXPECT_NEAR(gblSensVals[i], (l2penFD-l2pen)/fdTol, fdTol*std::fabs(gblSensVals[i]) + fdTol*fdTol);
      if(node.is_local_offset_valid()) krino::field_data<double>(ls1, node)[0] = origLSVal;
    }
  }
}

TEST(SignedDistancePenalty, GradPenaltySens)
{
  double meshSize = 1.0;
  std::vector<stk::topology> topologies{stk::topology::TRIANGLE_3_2D, stk::topology::TETRAHEDRON_4};
  
  for(auto && topo : topologies)
  {
    std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(topo, MPI_COMM_WORLD);
    auto & auxMeta = krino::AuxMetaData::get_or_create(bboxMesh->meta_data());
    auto ls1 = auxMeta.register_field("LS1", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());
    auto ls2 = auxMeta.register_field("LS2", krino::FieldType::REAL, stk::topology::NODE_RANK, 1, 1, bboxMesh->meta_data().universal_part());

    const stk::math::Vector3d minCorner{0., 0., 0.}; 
    const stk::math::Vector3d maxCorner{4., 4., 4.};
    populate_bbox_mesh_with_aura_and_activate(*bboxMesh, minCorner, maxCorner, meshSize);

    for(auto && b : bboxMesh->bulk_data().get_buckets(stk::topology::NODE_RANK, bboxMesh->meta_data().universal_part()))
    {
      for(auto && n : (*b))
      {
        auto * lsVal1 = krino::field_data<double>(ls1, n);
        auto * lsVal2 = krino::field_data<double>(ls2, n);
        lsVal1[0] = 0.;
        lsVal2[0] = 0.;
        auto coords = krino::get_vector_field(bboxMesh->bulk_data(), krino::FieldRef(bboxMesh->meta_data().coordinate_field()), n, topo.dimension());
        for(unsigned d=0; d<topo.dimension(); d++)
        {
          lsVal1[0] += coords[d]*coords[d];
          lsVal2[0] += 2.*coords[d]*coords[d];
        }
      }
    }

    std::map<stk::mesh::EntityId, double> sens;
    double gradPen = krino::compute_gradient_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
      bboxMesh->meta_data().universal_part(), sens);

    std::vector<stk::mesh::EntityId> lclSensIds;
    std::vector<double> lclSensVals;
    for(auto && entry : sens)
    {
      lclSensIds.push_back(entry.first);
      lclSensVals.push_back(entry.second);
    }
    std::vector<stk::mesh::EntityId> gblSensIds;
    std::vector<double> gblSensVals;
    stk::parallel_vector_concat(bboxMesh->bulk_data().parallel(), lclSensIds, gblSensIds);
    stk::parallel_vector_concat(bboxMesh->bulk_data().parallel(), lclSensVals, gblSensVals);
    const double fdTol = 1e-5;

    for(unsigned i=0; i<gblSensIds.size(); i++)
    {
      std::map<stk::mesh::EntityId, double> dummySens;
      double origLSVal = 0;
      auto node = bboxMesh->bulk_data().get_entity(stk::topology::NODE_RANK, gblSensIds[i]);
      if(node.is_local_offset_valid())
      {
        auto * lsVal = krino::field_data<double>(ls1, node);
        origLSVal = lsVal[0];
        lsVal[0] += fdTol;
      }
      double gradPenFD = krino::compute_gradient_penalty_between_fields(bboxMesh->bulk_data(), ls1.field(), ls2.field(),
        bboxMesh->meta_data().universal_part(), dummySens);
      const double gradTol = 1e-4;
      EXPECT_NEAR(gblSensVals[i], (gradPenFD-gradPen)/fdTol, gradTol*std::fabs(gblSensVals[i]) + gradTol*gradTol);
      if(node.is_local_offset_valid()) krino::field_data<double>(ls1, node)[0] = origLSVal;
    }
  }
}