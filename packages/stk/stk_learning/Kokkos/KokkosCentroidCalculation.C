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

// EXERCISE 2 Goal:
// Add Kokkos Views and deep copy from host to device
// 1. Define device views
// 2. Define host views using mirro of corresponding device views
// 3. Initialize the host views
// 4. Deep copy host view to device view
// Note: Kokkos::parallel_for() initializations were removed to initialize on host

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "mtk_kokkos.h"
#include <stk_util/stk_config.h>

#if KOKKOS_HAVE_CUDA
typedef double my_double;
#else
typedef long double my_double;
#endif

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


// typedef Kokkos::LayoutLeft   Layout ;

typedef Kokkos::RangePolicy<ExecSpace> range_policy ;

typedef Kokkos::View<my_double*, MemSpace>   DeviceViewVectorType;
typedef Kokkos::View<my_double*, Kokkos::HostSpace>   HostViewVectorType;

typedef Kokkos::TeamPolicy<ExecSpace>               team_policy ;
typedef Kokkos::TeamPolicy<ExecSpace>::member_type  member_type ;

typedef Kokkos::LayoutLeft   Layout ;

typedef Kokkos::View<double**, Layout, Kokkos::HostSpace>   HostViewMatrixType;
typedef Kokkos::View<double**, Layout, MemSpace>   DeviceViewMatrixType;

STK_FUNCTION void get_value(DeviceViewVectorType data, size_t index, my_double& value)
{
    value = data(index);
}

my_double get_average(const std::vector<my_double>& data, int num_repeat)
{
    HostViewVectorType host_data("host_data", data.size());
    DeviceViewVectorType device_data("device_data", data.size());

    for(size_t i=0;i<data.size();++i)
    {
        host_data(i) = data[i];
    }

    Kokkos::deep_copy(device_data, host_data);

    my_double sum = 0;
    for(int j=0;j<num_repeat;++j)
    {
        sum = 0;
        Kokkos::parallel_reduce( data.size(), KOKKOS_LAMBDA ( size_t i, my_double &local_sum )
        {
            my_double value = 0;
            get_value(device_data, i, value);
            local_sum += value;
        },
        sum );
    }

    my_double num = data.size();
    return sum/num;
}

TEST_F(MTK_Kokkos, calculate_average)
{
    size_t size_data = unitTestUtils::get_command_line_option<size_t>("-s", "2");
    std::vector<my_double> data(size_data);
    for(size_t i=0;i<data.size();++i)
        data[i]=i;
    int num_repeat = unitTestUtils::get_command_line_option<size_t>("-n", "1");

    struct timeval begin,end;
    gettimeofday(&begin,NULL);

    my_double average = get_average(data, num_repeat);

    gettimeofday(&end,NULL);

    double time = 1.0*(end.tv_sec-begin.tv_sec) +
                  1.0e-6*(end.tv_usec-begin.tv_usec);
    double Gbytes = 1.0e-9 * double(sizeof(my_double) * ( data.size() )) ;
    printf("bandwidth( %g GB/s )\n", Gbytes * num_repeat / time );

    my_double num_items = data.size();
    my_double sum = (num_items*(num_items-1))/2.0;
    my_double gold_ave = sum/data.size();
    EXPECT_NEAR(gold_ave, average, 0.0000001);
}

typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> CoordFieldType;

void calculate_centroids(const stk::mesh::BulkData& bulkData, const CoordFieldType& coordinates, CoordFieldType& centroid)
{
    const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::ELEM_RANK);
    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();++j)
        {
            stk::mesh::Entity element = bucket[j];

            unsigned num_nodes = bulkData.num_nodes(element);
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
            double xave = 0, yave = 0, zave = 0;
            for(unsigned k=0;k<num_nodes;++k)
            {
                double *node_coords = stk::mesh::field_data(coordinates, nodes[k]);
                xave += node_coords[0];
                yave += node_coords[1];
                zave += node_coords[2];
            }
            xave /= num_nodes;
            yave /= num_nodes;
            zave /= num_nodes;
            double *centroid_values = stk::mesh::field_data(centroid, element);
            centroid_values[0] = xave;
            centroid_values[1] = yave;
            centroid_values[2] = zave;
        }
    }
}

TEST_F(MTK_Kokkos, calculate_centroid_field_on_host)
{
    int num_repeat = unitTestUtils::get_command_line_option<size_t>("-n", "1");
    size_t dim = unitTestUtils::get_command_line_option<size_t>("-d", "10");
    std::ostringstream os;
    os << "generated:" << dim << "x" << dim << "x" << dim << std::endl;

    stk::mesh::MetaData meta(3);

    CoordFieldType& centroid = meta.declare_field<CoordFieldType>(stk::topology::ELEM_RANK, "centroid");
    std::vector<double> init_vec = {0,0,0};
    stk::mesh::put_field(centroid, meta.universal_part(), init_vec.data());

    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

    stk::unit_test_util::fill_mesh_using_stk_io(os.str(), bulk);

    CoordFieldType *coords = meta.get_field<CoordFieldType>(stk::topology::NODE_RANK, "coordinates");

    struct timeval begin,end;
    gettimeofday(&begin,NULL);

    calculate_centroids(bulk, *coords, centroid);

    gettimeofday(&end,NULL);

    double time = 1.0*(end.tv_sec-begin.tv_sec) +
                  1.0e-6*(end.tv_usec-begin.tv_usec);

    size_t num_coords = (dim+1)*(dim+1)*(dim+1)*3;
    double Gbytes = 1.0e-9 * double(sizeof(my_double) * ( num_coords )) ;
    std::cerr << "For mesh with " << num_coords << " coordinates.\n";
    printf("bandwidth( %g GB/s )\n", Gbytes * num_repeat / time );

    stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    double *centroid_values = stk::mesh::field_data(centroid, element1);
    EXPECT_NEAR(0.5, centroid_values[0], 0.000001);
    EXPECT_NEAR(0.5, centroid_values[1], 0.000001);
    EXPECT_NEAR(0.5, centroid_values[2], 0.000001);
}

struct ScratchData
{
    struct my_first_operator{};
    struct my_second_operator{};

    ScratchData()
    {
    }

    void initialize(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, CoordFieldType& centroid, const stk::mesh::Selector& selector)
    {
        stk::mesh::EntityVector nodes;
        stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), nodes);
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::ELEM_RANK), elements);
        unsigned numElements = elements.size();

        elementNodes = DeviceViewMatrixType("DElementsNumNodes", numElements, 8, bulk.mesh_meta_data().spatial_dimension());
        hostElementNodes =  Kokkos::create_mirror_view(elementNodes);

        for (unsigned elemIndex = 0; elemIndex < numElements; ++elemIndex)
        {
            stk::mesh::Entity element = elements[elemIndex];
            const stk::mesh::Entity * nodes = bulk.begin_nodes(element);
            const unsigned numNodesThisElem = bulk.num_nodes(element);
            for(unsigned iNode = 0; iNode < numNodesThisElem; ++iNode)
            {
                stk::mesh::Entity node = nodes[iNode];
                double *node_coords = stk::mesh::field_data(coords, node);
                for(unsigned k=0;k<bulk.mesh_meta_data().spatial_dimension();k++)
                {
                    hostElementNodes(elemIndex, iNode, k) = node_coords[k];
                }
            }
        }

        Kokkos::deep_copy(elementNodes, hostElementNodes);

        elementCentroids = DeviceViewMatrixType("Centroids", numElements, bulk.mesh_meta_data().spatial_dimension());
        hostElementCentroids =  Kokkos::create_mirror_view(elementCentroids);
    }

    void test_centroid_of_element_1()
    {
        for (unsigned k=0 ; k < hostElementCentroids.extent(1) ; ++k)
            EXPECT_NEAR(0.5, hostElementCentroids(0, k), 0.000001);
    }

    void calculate_centroids(int num_repeat, int choice)
    {
        for (int repeat=0 ; repeat<num_repeat ; ++repeat)
        {
            const unsigned num_elements = elementCentroids.extent(0);
            if (choice == 0)
            {
                Kokkos::parallel_for("first op", Kokkos::RangePolicy<my_first_operator>(0,num_elements), *this);
                //Kokkos::parallel_for("second op", Kokkos::RangePolicy<my_second_operator>(0,num_elements), *this);
            }
            else
            {
//                Kokkos::parallel_for("new op", team_policy(num_elements, Kokkos::AUTO), *this);
                //Kokkos::parallel_for("second op", Kokkos::RangePolicy<my_second_operator>(0,num_elements), *this);
            }
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(my_first_operator, const int element) const
    {
        const unsigned dim = elementCentroids.extent(1);
        for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
        {
            double temp = elementCentroids(element, k);
            for(unsigned node=0;node<8;++node) // loop over every node of this element
                temp += elementNodes(element, node, k); // sum the coordinates
            elementCentroids(element, k) = temp * 0.125;
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(my_second_operator, const int i) const
    {
        const unsigned dim = elementCentroids.extent(1);
        for (unsigned k=0 ; k<dim ; ++k)
          elementCentroids(i, k) *= 0.125; // divide by num nodes to get centroid
    }

//    KOKKOS_INLINE_FUNCTION void operator()(const member_type &teamMember) const
//    {
//        unsigned i = teamMember.league_rank();
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 8), [=] (int j)
//        {
//            const unsigned dim = elementCentroids.extent(1);
//            double temp = elementCentroids(i, j);
//            for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
//                temp += elementNodes(i, j, k); // sum the coordinates
//            elementCentroids(i, j) = temp * 0.125;
//        });
//    }

    void copy_centroids_from_host()
    {
        Kokkos::deep_copy(hostElementCentroids, elementCentroids);
    }

    DeviceViewMatrixType elementNodes;
    DeviceViewMatrixType elementCentroids;

    HostViewMatrixType::HostMirror hostElementNodes;
    HostViewMatrixType::HostMirror hostElementCentroids;
};


TEST_F(MTK_Kokkos, calculate_centroid_field_on_device)
{
    int num_repeat = unitTestUtils::get_command_line_option<int>("-n", "1");
    size_t dim = unitTestUtils::get_command_line_option<size_t>("-d", "10");
    int choice = unitTestUtils::get_command_line_option<int>("-c", "0");

    std::ostringstream os;
    os << "generated:" << dim << "x" << dim << "x" << dim << std::endl;

    stk::mesh::MetaData meta(3);

    CoordFieldType& centroid = meta.declare_field<CoordFieldType>(stk::topology::ELEM_RANK, "centroid");
    std::vector<double> init_vec = {0,0,0};
    stk::mesh::put_field(centroid, meta.universal_part(), init_vec.data());

    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

    stk::unit_test_util::fill_mesh_using_stk_io(os.str(), bulk);

    CoordFieldType *coords = meta.get_field<CoordFieldType>(stk::topology::NODE_RANK, "coordinates");

    ScratchData scratch;
    scratch.initialize(bulk, *coords, centroid, meta.locally_owned_part());

    struct timeval begin,end;
    gettimeofday(&begin,NULL);

    scratch.calculate_centroids(num_repeat, choice);

    gettimeofday(&end,NULL);

    double time = 1.0*(end.tv_sec-begin.tv_sec) +
                  1.0e-6*(end.tv_usec-begin.tv_usec);

    size_t num_coords = (dim)*(dim)*(dim)*10*3;
    double Gbytes = 1.0e-9 * double(sizeof(my_double) * ( num_coords )) ;
    std::cerr << "For mesh with " << num_coords << " coordinates.\n";
    printf("bandwidth( %g GB/s )\n", Gbytes * num_repeat / time );

    scratch.copy_centroids_from_host();
    scratch.test_centroid_of_element_1();
}
