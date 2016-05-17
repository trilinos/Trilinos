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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
 
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

namespace ngp {

class Bulk {
public:
    STK_FUNCTION
    Bulk(const stk::mesh::BulkData& bulk)
     : mesh_indices(bulk.get_size_of_entity_index_space())
    {
        fill_mesh_indices(stk::topology::NODE_RANK, bulk);
        copy_mesh_indices_to_device();
    }

    STK_FUNCTION
    ~Bulk(){}

    STK_FUNCTION
    const stk::mesh::FastMeshIndex& mesh_index(stk::mesh::Entity entity) const
    {
        return device_mesh_indices(entity.local_offset());
    }

private:
    void fill_mesh_indices(stk::mesh::EntityRank rank, const stk::mesh::BulkData& bulk)
    {
        const stk::mesh::BucketVector& bkts = bulk.buckets(rank);

        for(const stk::mesh::Bucket* bktptr : bkts) {
            const stk::mesh::Bucket& bkt = *bktptr;
            for(size_t i=0; i<bkt.size(); ++i) {
                mesh_indices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex(bkt.bucket_id(), i);
            }
        }
    }

    void copy_mesh_indices_to_device()
    {
        Kokkos::View<stk::mesh::FastMeshIndex*,Kokkos::HostSpace> host_mesh_indices("host_mesh_indices",mesh_indices.size());
        for(size_t i=0; i<mesh_indices.size(); ++i) {
            host_mesh_indices(i) = mesh_indices[i];
        }
        Kokkos::View<stk::mesh::FastMeshIndex*> tmp_device_mesh_indices("tmp_dev_mesh_indices", mesh_indices.size());
        Kokkos::deep_copy(tmp_device_mesh_indices, host_mesh_indices);
        device_mesh_indices = tmp_device_mesh_indices;
    }

    std::vector<stk::mesh::FastMeshIndex> mesh_indices;
    Kokkos::View<const stk::mesh::FastMeshIndex*,Kokkos::MemoryTraits<Kokkos::RandomAccess> > device_mesh_indices;
//    Kokkos::View<const stk::mesh::FastMeshIndex*> device_mesh_indices;
};

template<typename T>
class Field {
public:
    STK_FUNCTION
    Field(stk::mesh::EntityRank rank, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : ngpBulk(bulk), ngpBulkRef(ngpBulk), device_data("device_data",0),
      InvalidFieldValue(-999999.9)
    {
        unsigned num_entities = count_entities(rank, bulk, selector);
        create_device_field_data(num_entities, initialValue);
    }

    STK_FUNCTION
    ~Field(){}

    STK_FUNCTION
    T operator[](stk::mesh::Entity entity) const
    {
        const ngp::Bulk& localRef = ngpBulk;
        unsigned bkt_id = localRef.mesh_index(entity).bucket_id;
        unsigned bkt_ord = localRef.mesh_index(entity).bucket_ord;
        unsigned idx = bkt_id*512 + bkt_ord;
        return device_data(idx);
    }

private:
    unsigned count_entities(stk::mesh::EntityRank rank, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    {
        const stk::mesh::BucketVector& bkts = bulk.get_buckets(rank, selector);
        unsigned num_entities = 0;
        for(const stk::mesh::Bucket* bktptr : bkts) {
            const stk::mesh::Bucket& bkt = *bktptr;
            num_entities += bkt.size();
        }
        return num_entities;
    }

    void create_device_field_data(unsigned num_entities, const T& initialValue)
    {
        Kokkos::View<T*,Kokkos::HostSpace> host_data("host_data", num_entities);
        for(size_t i=0; i<num_entities; ++i) {
            host_data(i) = initialValue;
        }
        Kokkos::View<T*> tmp("tmp", num_entities);
        device_data = tmp;
        Kokkos::deep_copy(device_data, host_data);
    }

    ngp::Bulk ngpBulk;
    const ngp::Bulk& ngpBulkRef;
    Kokkos::View<T*> device_data;
    T InvalidFieldValue;
};

}

void test_field() {

  unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  unsigned dimX = 10;
  unsigned dimY = 10;
  unsigned dimZ = 10;
  std::ostringstream mesh_spec;
  mesh_spec << "generated:"<<dimX<<"x"<<dimY<<"x"<<dimZ;
  stk::unit_test_util::fill_mesh_using_stk_io(mesh_spec.str(), bulk);

  double initialValue = 99.9;
  ngp::Field<double> scalarField(stk::topology::NODE_RANK, initialValue, bulk, meta.locally_owned_part());

  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), nodes);
  Kokkos::View<stk::mesh::Entity*> device_nodes("device_nodes", nodes.size());
  Kokkos::View<stk::mesh::Entity*,Kokkos::HostSpace> host_nodes("host_nodes", nodes.size());

  for(size_t i=0; i<nodes.size(); ++i) {
      host_nodes(i) = nodes[i];
  }
  Kokkos::deep_copy(device_nodes, host_nodes);

  struct timeval begin,end;

  gettimeofday(&begin,NULL);

  int nrepeat = 500;
  double result = 0;

  for(int n=0; n<nrepeat; ++n) {

      result = 0;
      Kokkos::parallel_reduce(nodes.size(), KOKKOS_LAMBDA(int i, double& update) {
        update += scalarField[device_nodes(i)];
      }, result);

  }

  gettimeofday(&end,NULL);

  Kokkos::fence();

  double time = 1.0*(end.tv_sec-begin.tv_sec) +
                1.0e-6*(end.tv_usec-begin.tv_usec);

  std::cout << "Time: " << time << std::endl;

  double goldValue = nodes.size()*initialValue;
  EXPECT_NEAR(goldValue, result, 0.0001);
}

TEST_F(MTK_Kokkos, field_exp) {
  test_field();
}

