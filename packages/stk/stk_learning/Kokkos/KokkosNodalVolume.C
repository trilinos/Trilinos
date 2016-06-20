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

#include "mtk_kokkos.h"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
#include "../../stk_mesh/stk_mesh/base/FieldParallel.hpp"

#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda     ExecSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_HAVE_OPENMP
   typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_HAVE_CUDA
   typedef Kokkos::CudaSpace    MemSpace;
#else
   typedef Kokkos::HostSpace    MemSpace;
#endif

double calculate_element_volume()
{
    return 1.0;
}

void calculate_nodal_volume(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField)
{
//    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part());
    for(const stk::mesh::Bucket* bucket : elemBuckets) {
        const unsigned numNodesPerElem = bucket->topology().num_nodes();
        for(size_t i=0; i<bucket->size(); ++i) {
//            stk::mesh::Entity elem = (*bucket)[i];
            const stk::mesh::Entity* elemNodes = bucket->begin_nodes(i);
//            double* nodalCoords = static_cast<double*>(stk::mesh::field_data(coords, elemNodes[j]));
            double elemVolumePerNode = calculate_element_volume()/numNodesPerElem;
            for(unsigned j=0; j<numNodesPerElem; ++j) {
                double* nodalVolume = stk::mesh::field_data(nodalVolumeField, elemNodes[j]);
                *nodalVolume += elemVolumePerNode;
            }
        }
    }
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void expect_nodal_volume(const stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, stk::mesh::Field<int>& numElemsPerNodeField)
{
    stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();
    const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, selector);
    for(const stk::mesh::Bucket* bucket : nodeBuckets) {
        const double* nodalVolume = stk::mesh::field_data(nodalVolumeField, *bucket);
        const int *numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *bucket);
        for(size_t i=0; i<bucket->size(); ++i) {
            double expectedVolume = 1.0;
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

            EXPECT_NEAR(expectedVolume, nodalVolume[i], 1.e-9);
        }
    }
}

void count_num_elems_per_node(const stk::mesh::BulkData& mesh, stk::mesh::Field<int>& numElemsPerNodeField)
{
    stk::mesh::Selector ownedOrShared = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
    const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, ownedOrShared);
    for(const stk::mesh::Bucket* bucket : nodeBuckets) {
        int* numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *bucket);
        for(size_t i=0; i<bucket->size(); ++i) {
            unsigned numElems = bucket->num_elements(i);
            const stk::mesh::Entity *elems = bucket->begin_elements(i);
            for(size_t j=0; j<numElems; ++j)
            {
                if(mesh.bucket(elems[j]).owned())
                {
                    numElemsPerNode[i]++;
                }
            }
        }
    }
    stk::mesh::parallel_sum(mesh, {&numElemsPerNodeField});
}

class NodalVolumeCalculator : public stk::unit_test_util::MeshFixture
{
protected:
    void test_nodal_volume(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        auto& nodalVolumeField = get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodal_volume");
        stk::mesh::put_field(nodalVolumeField, get_meta().universal_part());
        auto& numElemsPerNodeField = get_meta().declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "numElemsPerNode");
        stk::mesh::put_field(numElemsPerNodeField, get_meta().universal_part());
        setup_mesh("generated:2x2x2", auraOption);
        count_num_elems_per_node(get_bulk(), numElemsPerNodeField);

        calculate_nodal_volume(get_bulk(), nodalVolumeField);

        expect_nodal_volume(get_bulk(), nodalVolumeField, numElemsPerNodeField);
    }
};
TEST_F(NodalVolumeCalculator, nodalVolumeWithAura) {
    test_nodal_volume(stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAura) {
    test_nodal_volume(stk::mesh::BulkData::NO_AUTO_AURA);
}
