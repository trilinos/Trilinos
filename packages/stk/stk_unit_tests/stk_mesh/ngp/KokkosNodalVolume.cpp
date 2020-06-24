/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/NgpAtomics.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include "stk_mesh/base/FieldParallel.hpp"

#include <limits>

double calculate_element_volume(const stk::mesh::Entity* elemNodes,
                                unsigned numElemNodes,
                                const stk::mesh::FieldBase& coords)
{
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  double max[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
  for(unsigned i=0; i<numElemNodes; ++i) {
    const double* val = static_cast<const double*>(stk::mesh::field_data(coords, elemNodes[i]));
    for(int j=0; j<3; ++j) {
      if (val[j] > max[j])
        max[j] = val[j];
      if (val[j] < min[j])
        min[j] = val[j];
    }
  }
  return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

void calculate_nodal_volume_given_elem_nodes_stkmesh(const stk::mesh::Entity* elemNodes,
                                                     const unsigned numNodesPerElem,
                                                     const stk::mesh::FieldBase& coords,
                                                     stk::mesh::Field<double>& nodalVolumeField)
{
  double elemVolumePerNode = calculate_element_volume(elemNodes, numNodesPerElem, coords) / numNodesPerElem;
  for(unsigned j = 0; j < numNodesPerElem; ++j)
  {
    double* nodalVolume = stk::mesh::field_data(nodalVolumeField, elemNodes[j]);
    stk::mesh::atomic_add(nodalVolume, elemVolumePerNode);
  }
}

typedef Kokkos::TeamPolicy<stk::mesh::HostExecSpace, stk::mesh::ScheduleType> HostTeamPolicyType;
typedef HostTeamPolicyType::member_type HostTeamHandleType;

#ifndef KOKKOS_ENABLE_CUDA
void calculate_nodal_volume_stkmesh(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
  const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();

  const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK,
                                                                mesh.mesh_meta_data().locally_owned_part());
  unsigned numBuckets = elemBuckets.size();
  for(int n=0; n<numRepeat; ++n)
  {
    Kokkos::parallel_for(Kokkos::TeamPolicy< stk::mesh::HostExecSpace >(numBuckets, Kokkos::AUTO),
                         [&] (const HostTeamHandleType& team)
    {
      const int elementBucketIndex = team.league_rank();
      const stk::mesh::Bucket &bucket = *elemBuckets[elementBucketIndex];
      const unsigned numNodesPerElem = bucket.topology().num_nodes();
      unsigned numElements = bucket.size();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
      {
        const stk::mesh::Entity* elemNodes = bucket.begin_nodes(i);
        calculate_nodal_volume_given_elem_nodes_stkmesh(elemNodes, numNodesPerElem, coords, nodalVolumeField);
      });
    });
  }

  stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}
#endif

void calculate_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData& mesh,
                                                stk::mesh::Field<double>& nodalVolumeField,
                                                int numRepeat)
{
  const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
  for(int n=0; n<numRepeat; ++n)
  {
    stk::mesh::for_each_entity_run(mesh, stk::topology::ELEM_RANK, [&](const stk::mesh::BulkData& bulk, stk::mesh::MeshIndex elem)
    {
      const stk::mesh::Entity* elemNodes = elem.bucket->begin_nodes(elem.bucket_ordinal);
      const unsigned numNodesPerElem = elem.bucket->topology().num_nodes();
      calculate_nodal_volume_given_elem_nodes_stkmesh(elemNodes, numNodesPerElem, coords, nodalVolumeField);
    });
  }
  stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}





STK_FUNCTION
void calculate_nodal_volume_given_elem_volume(const stk::mesh::NgpMesh &ngpMesh,
                                              stk::mesh::NgpMesh::ConnectedNodes nodes,
                                              double elemVolumePerNode,
                                              const stk::mesh::NgpField<double>& nodalVolume)
{
  const unsigned numNodesPerElem = nodes.size();
  for(unsigned j = 0; j < numNodesPerElem; ++j)
  {
    double* v = &nodalVolume(ngpMesh.fast_mesh_index(nodes[j]), 0);
    stk::mesh::atomic_add(v, elemVolumePerNode);
  }
}

STK_FUNCTION
double calculate_element_volume(const stk::mesh::NgpMesh &ngpMesh,
                                stk::mesh::NgpMesh::ConnectedNodes nodes,
                                const stk::mesh::NgpField<double> &coords)
{
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  double max[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
  const unsigned numElemNodes = nodes.size();
  for(unsigned i=0; i<numElemNodes; ++i) {
    stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);
    for(int j=0; j<3; ++j) {
      double val = coords(nodeIndex, j);
      if (val > max[j])
        max[j] = val;
      if (val < min[j])
        min[j] = val;
    }
  }
  return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

STK_FUNCTION
void calculate_nodal_volume_device(const stk::mesh::NgpMesh &ngpMesh,
                                   const stk::mesh::FastMeshIndex& elem,
                                   const stk::mesh::NgpField<double> &coords,
                                   const stk::mesh::NgpField<double> &nodalVolume)
{
  stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
  const unsigned numNodesPerElem = nodes.size();
  double elemVolumePerNode = calculate_element_volume(ngpMesh, nodes, coords) / numNodesPerElem;
  calculate_nodal_volume_given_elem_volume(ngpMesh, nodes, elemVolumePerNode, nodalVolume);
}

void calculate_nodal_volume(stk::mesh::NgpMesh &ngpMesh, stk::mesh::Selector selector, const stk::mesh::NgpField<double> &coords, stk::mesh::NgpField<double> &nodalVolume)
{
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                 {
                                   calculate_nodal_volume_device(ngpMesh, elem, coords, nodalVolume);
                                 });
  nodalVolume.modify_on_device();
}

void calculate_nodal_volume_entity_loop(stk::mesh::BulkData& mesh,
                                        stk::mesh::Selector selector,
                                        stk::mesh::Field<double>& nodalVolumeField,
                                        int numRepeat)
{
  double start = stk::wall_time();
  const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
  stk::mesh::NgpField<double>& ngpCoords = stk::mesh::get_updated_ngp_field<double>(coords);
  stk::mesh::NgpField<double>& ngpNodalVolume = stk::mesh::get_updated_ngp_field<double>(nodalVolumeField);
  stk::mesh::NgpMesh & ngpMesh = mesh.get_updated_ngp_mesh();
  double middle = stk::wall_time();

  for(int n=0; n<numRepeat; ++n)
    calculate_nodal_volume(ngpMesh, selector, ngpCoords, ngpNodalVolume);

  ngpNodalVolume.sync_to_host();
  stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
  double end = stk::wall_time();
  std::cerr << "Init: " << middle - start << ", Calc: " << end - middle << std::endl;
}

template <typename Algorithm>
void repeat_for_each_entity_loop_for_algorithm(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat,
                                               const Algorithm &algorithm)
{
  stk::mesh::NgpField<double>& staticNodalVolume = stk::mesh::get_updated_ngp_field<double>(nodalVolumeField);
  stk::mesh::NgpField<double>& staticElemVolume = stk::mesh::get_updated_ngp_field<double>(elemVolumePerNodeField);
  stk::mesh::NgpMesh & ngpMesh = mesh.get_updated_ngp_mesh();

  for(int n=0; n<numRepeat; ++n)
    algorithm(ngpMesh, staticElemVolume, staticNodalVolume);

  staticNodalVolume.sync_to_host();
  stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

STK_FUNCTION
void assemble_nodal_volume_from_elem_volume_device(const stk::mesh::NgpMesh &ngpMesh,
                                                   const stk::mesh::FastMeshIndex& elem,
                                                   const double elemVolumePerNode,
                                                   const stk::mesh::NgpField<double> &staticNodalVolume)
{
  stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
  calculate_nodal_volume_given_elem_volume(ngpMesh, nodes, elemVolumePerNode, staticNodalVolume);
}

void assemble_nodal_volume_from_elem_volume_bucket_field_access(stk::mesh::NgpMesh &ngpMesh,
                                                                const stk::mesh::NgpField<double> &elemVolumePerNode,
                                                                stk::mesh::NgpField<double> &staticNodalVolume)
{
  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  unsigned numBuckets = ngpMesh.num_buckets(stk::topology::ELEM_RANK);
  Kokkos::parallel_for(Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA (const TeamHandleType& team)
                       {
                         const int bucketIndex = team.league_rank();
                         const stk::mesh::NgpMesh::BucketType &bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
                         unsigned numElements = bucket.size();
                         double *volumePerNodePtr = &elemVolumePerNode(stk::mesh::FastMeshIndex{bucket.bucket_id(), 0}, 0);
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
                         {
                           stk::mesh::FastMeshIndex elem{bucket.bucket_id(), static_cast<unsigned>(i)};
                           assemble_nodal_volume_from_elem_volume_device(ngpMesh, elem, volumePerNodePtr[i], staticNodalVolume);
                         });
                       });
  staticNodalVolume.modify_on_device();
}

void assemble_nodal_volume_bucket_field_access(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat)
{
  repeat_for_each_entity_loop_for_algorithm(mesh, elemVolumePerNodeField, nodalVolumeField, numRepeat,
                                            [=](stk::mesh::NgpMesh &ngpMesh, const stk::mesh::NgpField<double> &elemVolumePerNode, stk::mesh::NgpField<double> &staticNodalVolume)
  {
    assemble_nodal_volume_from_elem_volume_bucket_field_access(ngpMesh, elemVolumePerNode, staticNodalVolume);
  });
}

void assemble_nodal_volume_from_elem_volume_entity_field_access(stk::mesh::NgpMesh &ngpMesh,
                                                                stk::mesh::Selector selector,
                                                                const stk::mesh::NgpField<double> &elemVolumePerNode,
                                                                stk::mesh::NgpField<double> &staticNodalVolume)
{
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                 {
                                   double volumePerNode = elemVolumePerNode(elem, 0);
                                   assemble_nodal_volume_from_elem_volume_device(ngpMesh, elem, volumePerNode, staticNodalVolume);
                                 });
  staticNodalVolume.modify_on_device();
}

void assemble_nodal_volume_entity_field_access(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat)
{
  stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();
  repeat_for_each_entity_loop_for_algorithm(mesh, elemVolumePerNodeField, nodalVolumeField, numRepeat,
                                            [=](stk::mesh::NgpMesh &ngpMesh, const stk::mesh::NgpField<double> &elemVolumePerNode, stk::mesh::NgpField<double> &staticNodalVolume)
  {
    assemble_nodal_volume_from_elem_volume_entity_field_access(ngpMesh, selector, elemVolumePerNode, staticNodalVolume);
  });
}




void expect_nodal_volume(const stk::mesh::BulkData& mesh, stk::mesh::Selector selector, stk::mesh::Field<double>& nodalVolumeField, stk::mesh::Field<int>& numElemsPerNodeField, int numRepeat)
{
  const stk::mesh::BucketVector& nodeBuckets = mesh.buckets(stk::topology::NODE_RANK);
  for(const stk::mesh::Bucket* bucket : nodeBuckets) {
    const double* nodalVolume = stk::mesh::field_data(nodalVolumeField, *bucket);
    const int *numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *bucket);
    for(size_t i=0; i<bucket->size(); ++i) {
      double expectedVolume = 1.0;
      if(selector(*bucket))
      {
        if (numElemsPerNode[i] == 8) {
          expectedVolume = 1.0;
        }
        else if (numElemsPerNode[i] == 4) {
          expectedVolume = 0.5;
        }
        else if (numElemsPerNode[i] == 2) {
          expectedVolume = 0.25;
        }
        else if (numElemsPerNode[i] == 1) {
          expectedVolume = 0.125;
        }
        else {
          ASSERT_TRUE(false) << "numElemsPerNode = " << numElemsPerNode[i];
        }
        EXPECT_NEAR(numRepeat*expectedVolume, nodalVolume[i], 1.e-9);
      }
    }
  }
}

void count_num_elems_per_node(const stk::mesh::BulkData& mesh, stk::mesh::Selector selector, stk::mesh::Field<int>& numElemsPerNodeField)
{
  stk::mesh::Selector ownedOrShared = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
  const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, ownedOrShared);
  for(const stk::mesh::Bucket* nodeBucket : nodeBuckets) {
    int* numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *nodeBucket);
    for(size_t i=0; i<nodeBucket->size(); ++i) {
      unsigned numElems = nodeBucket->num_elements(i);
      const stk::mesh::Entity *elems = nodeBucket->begin_elements(i);
      for(size_t j=0; j<numElems; ++j)
      {
        const stk::mesh::Bucket &elemBucket = mesh.bucket(elems[j]);
        if(elemBucket.owned() && selector(elemBucket))
        {
          numElemsPerNode[i]++;
        }
      }
    }
  }
  stk::mesh::parallel_sum(mesh, {&numElemsPerNodeField});
}

void fill_elem_volume_per_node(const stk::mesh::BulkData& mesh, stk::mesh::Selector selector, stk::mesh::Field<double>& elemVolumePerNodeField)
{
  const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
  for(const stk::mesh::Bucket *bucket : mesh.get_buckets(stk::topology::ELEM_RANK, selector))
  {
    unsigned numNodesPerElem = bucket->topology().num_nodes();
    double *volume = stk::mesh::field_data(elemVolumePerNodeField, *bucket);
    for(size_t j=0; j<bucket->size(); j++)
    {
      const stk::mesh::Entity *nodes = bucket->begin_nodes(j);
      volume[j] = calculate_element_volume(nodes, bucket->topology().num_nodes(), coords) / numNodesPerElem;
    }
  }
}

class NodalVolumeCalculator : public stk::unit_test_util::MeshFixture
{
protected:
  void create_mesh_and_fields(stk::mesh::BulkData::AutomaticAuraOption auraOption, const stk::mesh::Part &part)
  {
    elemVolumePerNodeField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::ELEM_RANK, "elemVolumePerNodeField");
    stk::mesh::put_field_on_mesh(*elemVolumePerNodeField, part,
                                 (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);
    nodalVolumeField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodal_volume");
    stk::mesh::put_field_on_mesh(*nodalVolumeField, part,
                                 (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);
    numElemsPerNodeField = &get_meta().declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "numElemsPerNode");
    stk::mesh::put_field_on_mesh(*numElemsPerNodeField, part,
                                 (stk::mesh::FieldTraits<stk::mesh::Field<int> >::data_type*) nullptr);
    std::string meshName = stk::unit_test_util::get_mesh_spec("-dim");
    setup_mesh(meshName, auraOption);
    std::cout << "Using mesh: "<<meshName<<", numRepeat: "<<get_num_repeat()<<std::endl;
  }

  void fill_fields(const stk::mesh::Part &part)
  {
    count_num_elems_per_node(get_bulk(), part, *numElemsPerNodeField);
    fill_elem_volume_per_node(get_bulk(), part, *elemVolumePerNodeField);
  }

  void create_mesh_and_fill_fields(stk::mesh::BulkData::AutomaticAuraOption auraOption, const stk::mesh::Part &part)
  {
    create_mesh_and_fields(auraOption, part);
    fill_fields(part);
  }

  void move_half_mesh_to_part(stk::mesh::Part &part)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(get_bulk().mesh_meta_data().locally_owned_part(),
                                     get_bulk().buckets(stk::topology::ELEM_RANK),
                                     elements);
    size_t start = elements.size() / 2;

    stk::mesh::Part *block1 = get_meta().get_part("block_1");
    get_bulk().modification_begin();
    for(size_t i=start; i<elements.size(); i++)
    {
      get_bulk().change_entity_parts(elements[i], stk::mesh::ConstPartVector{&part}, stk::mesh::ConstPartVector{block1});
    }
    get_bulk().modification_end();
  }

  int get_num_repeat()
  {
    return stk::unit_test_util::get_command_line_option<int>("-n", 20);
  }

#ifndef KOKKOS_ENABLE_CUDA
  void test_nodal_volume_stkmesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_mesh_and_fill_fields(auraOption, get_meta().universal_part());
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    calculate_nodal_volume_stkmesh(get_bulk(), *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tstkmesh: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), get_meta().universal_part(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }
#endif

  void test_nodal_volume_bucket_field_access(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_mesh_and_fill_fields(auraOption, get_meta().universal_part());
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    assemble_nodal_volume_bucket_field_access(get_bulk(), *elemVolumePerNodeField, *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tbucket field access: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), get_meta().universal_part(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }

  void test_nodal_volume_entity_field_access(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_mesh_and_fill_fields(auraOption, get_meta().universal_part());
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    assemble_nodal_volume_entity_field_access(get_bulk(), *elemVolumePerNodeField, *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tentity field access: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), get_meta().universal_part(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }

  void test_nodal_volume_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_mesh_and_fill_fields(auraOption, get_meta().universal_part());
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    calculate_nodal_volume_entity_loop(get_bulk(), get_meta().universal_part(), *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tentity_loop: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), get_meta().universal_part(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }

  void test_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_mesh_and_fill_fields(auraOption, get_meta().universal_part());
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    calculate_nodal_volume_stkmesh_entity_loop(get_bulk(), *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tstkmesh_entity_loop: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), get_meta().universal_part(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }

  void test_nodal_volume_subset(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    stk::mesh::Part &block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    create_mesh_and_fields(auraOption, block2);
    move_half_mesh_to_part(block2);
    fill_fields(block2);
    int numRepeat = get_num_repeat();
    double startTime = stk::wall_time();
    calculate_nodal_volume_entity_loop(get_bulk(), block2, *nodalVolumeField, numRepeat);
    double elapsedTime = stk::wall_time() - startTime;
    std::cout<<"\t\tentity_loop: "<<elapsedTime<<" seconds"<<std::endl;
    expect_nodal_volume(get_bulk(), block2, *nodalVolumeField, *numElemsPerNodeField, numRepeat);
  }

  stk::mesh::Field<double> *elemVolumePerNodeField;
  stk::mesh::Field<double> *nodalVolumeField;
  stk::mesh::Field<int> *numElemsPerNodeField;
};
TEST_F(NodalVolumeCalculator, nodalVolumeForEachEntityLoop_primerTest) {
  test_nodal_volume_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
}
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(NodalVolumeCalculator, stkMeshNodalVolume) {
  test_nodal_volume_stkmesh(stk::mesh::BulkData::NO_AUTO_AURA);
}
#endif
TEST_F(NodalVolumeCalculator, nodalVolumeBucketFieldAccess) {
  test_nodal_volume_bucket_field_access(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeEntityFieldAccess) {
  test_nodal_volume_entity_field_access(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeForEachEntityLoop) {
  test_nodal_volume_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, stkMeshNodalVolumeForEachEntityLoop) {
  test_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeSubset) {
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;
  test_nodal_volume_subset(stk::mesh::BulkData::NO_AUTO_AURA);
}


