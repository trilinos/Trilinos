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


#ifndef _KokkosBulkDataCentroidCalculation_h_
#define _KokkosBulkDataCentroidCalculation_h_

#include "KokkosCentroidCalculation.hpp"

#include <sys/time.h>
#include <gtest/gtest.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/Field.hpp>

typedef Kokkos::View<stk::mesh::Entity*, MemSpace> DeviceViewEntitiesType;

typedef Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> DeviceViewMeshIndicesType;

typedef Kokkos::View<const stk::mesh::FastMeshIndex*, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> > ConstDeviceViewMeshIndicesType;

typedef Kokkos::View<stk::mesh::Entity*, MemSpace> EntityViewType;

typedef Kokkos::View<stk::mesh::Entity**, Layout, MemSpace> BucketConnectivityType;

typedef Kokkos::View<stk::mesh::Entity*, Layout, MemSpace> DeviceViewFlatConnectivityType;

typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;
typedef Kokkos::TeamPolicy<ExecSpace, ScheduleType> TeamPolicyType;
typedef TeamPolicyType::member_type TeamHandleType;

typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> CoordFieldType;

#define STRINGIFY(x) #x

#define DEFINE_OPERATOR(PARALLELSCOPE, TEAMTYPE, EXPANSIONTYPE, ISTEAM)	\
  struct PARALLELSCOPE##_##TEAMTYPE##_##EXPANSIONTYPE##_##operator {\
      const char *name = STRINGIFY( PARALLELSCOPE##_##TEAMTYPE##_##EXPANSIONTYPE##_##operator ); \
      const bool  isTeam = ISTEAM; \
  };\

#define TYPE_OPERATOR(PARALLELSCOPE, TEAMTYPE, EXPANSIONTYPE)	\
  PARALLELSCOPE##_##TEAMTYPE##_##EXPANSIONTYPE##_##operator\

#define NAME_OPERATOR(PARALLELSCOPE, TEAMTYPE, EXPANSIONTYPE)	\
  STRINGIFY( PARALLELSCOPE##_##TEAMTYPE##_##EXPANSIONTYPE##_##operator )\

DEFINE_OPERATOR(bucket , solo, compact, false)
DEFINE_OPERATOR(bucket , solo, unroll , false)
DEFINE_OPERATOR(bucket , team, compact, true)
DEFINE_OPERATOR(bucket , team, unroll , true)
DEFINE_OPERATOR(element, solo, compact, false)
DEFINE_OPERATOR(element, solo, unroll , false)
DEFINE_OPERATOR(element, team, compact, true)
DEFINE_OPERATOR(element, team, unroll , true)

void calculate_centroids_on_host(const stk::mesh::BulkData& bulkData, const CoordFieldType& coordinates, CoordFieldType& centroid, stk::mesh::Selector selector);

template<class FUNCTOR>
class CentroidCalculator {
public:
    CentroidCalculator(FUNCTOR f) : functor(f) {}

    void test_centroid_of_element_1()
    {
        for (unsigned k=0 ; k < functor.hostElementCentroids.extent(1) ; ++k)
            EXPECT_NEAR(0.5, functor.hostElementCentroids(0, k), 0.000001);
    }

    void test_centroid_of_element(CoordFieldType& expectedCentroid, 
				  const stk::mesh::Entity element, 
				  const unsigned elementIndex)
    {
        double *expectedCentroidValues = stk::mesh::field_data(expectedCentroid, element);
        for (unsigned k=0 ; k < functor.hostElementCentroids.extent(1) ; ++k) {
            EXPECT_NEAR(expectedCentroidValues[k], functor.hostElementCentroids(elementIndex, k), 0.000001);
	}
    }

    void test_centroid_of_element(const std::vector<double>& expectedCentroidValues, 
				  const stk::mesh::Entity element)
    {
        std::vector<double> elementCentroid = functor.get_centroid_of_element(element);
        EXPECT_EQ(elementCentroid.size(), elementCentroid);
	
        for (unsigned k=0 ; k < functor.hostElementCentroids.extent(1) ; ++k) {
            EXPECT_NEAR(expectedCentroidValues[k], elementCentroid[k], 0.000001);
	}
    }

    template<class T>
    void execute_solo(T choice, int num_repeat)
    {
        std::cout << "Executing operator: " << choice.name << std::endl;
        size_t N = functor.getNumParallelItems();
        for (int repeat=0 ; repeat<num_repeat ; ++repeat)
	    Kokkos::parallel_for(choice.name, Kokkos::RangePolicy< T >(0,N), functor);
    }
    
    template<class T>
    void execute_team(T choice, int num_repeat, int team_size)
    {
        std::cout << "Executing operator: " << choice.name << std::endl;
        size_t N = functor.getNumParallelItems();
	    for (int repeat=0 ; repeat<num_repeat ; ++repeat)
	        Kokkos::parallel_for(choice.name, Kokkos::TeamPolicy< T >(N, team_size), functor);
    }

    void calculate_centroids(int num_repeat, int choice, int team_size)
    {
        TYPE_OPERATOR(bucket , solo, compact) choice_0;
        TYPE_OPERATOR(bucket , team, compact) choice_1;
        TYPE_OPERATOR(element, solo, compact) choice_2;
        TYPE_OPERATOR(bucket,  team, unroll ) choice_3;
        TYPE_OPERATOR(element, team, unroll ) choice_4;
        TYPE_OPERATOR(element, solo, unroll ) choice_5;
        TYPE_OPERATOR(bucket , solo, unroll ) choice_6;

        switch(choice)
        {
            case 0: execute_solo(choice_0, num_repeat           ); break;
            case 1: execute_team(choice_1, num_repeat, team_size); break;
            case 2: execute_solo(choice_2, num_repeat           ); break;
            case 3: execute_team(choice_3, num_repeat, team_size); break;
            case 4: execute_team(choice_4, num_repeat, team_size); break;
            case 5: execute_solo(choice_5, num_repeat           ); break;
            case 6: execute_solo(choice_6, num_repeat           ); break;
            default: printf("No current implementation available for choice %d\n",choice); break;
        }
    }

    void copy_centroids_to_host()
    {
        Kokkos::deep_copy(functor.hostElementCentroids, functor.elementCentroids);
    }

    FUNCTOR functor;
};

struct MyApp {
    MyApp()
    : meta(3), bulk(nullptr),
      centroid(meta.declare_field<CoordFieldType>(stk::topology::ELEM_RANK, "centroid")),
      hostCentroid(meta.declare_field<CoordFieldType>(stk::topology::ELEM_RANK, "hostCentroid")),
      coords(nullptr)
    {
        num_repeat = stk::unit_test_util::get_command_line_option<int>("-n", 1);
        dim = stk::unit_test_util::get_command_line_option<size_t>("-d", 10);
        choice = stk::unit_test_util::get_command_line_option<int>("-c", 0);
        teamSize = stk::unit_test_util::get_command_line_option<int>("-t", 1);
    
        std::ostringstream os;
        os << "generated:" << dim << "x" << dim << "x" << dim << std::endl;

        std::vector<double> init_vec = {0,0,0};
        stk::mesh::put_field_on_mesh(centroid, meta.universal_part(), init_vec.data());
        stk::mesh::put_field_on_mesh(hostCentroid, meta.universal_part(), init_vec.data());
    
        bulk = new stk::mesh::BulkData(meta, MPI_COMM_WORLD);
    
        stk::io::fill_mesh(os.str(), *bulk);
    
        coords = meta.get_field<CoordFieldType>(stk::topology::NODE_RANK, "coordinates");

	calculate_centroids_on_host(*bulk, *coords, hostCentroid, bulk->mesh_meta_data().locally_owned_part());
    }

    ~MyApp() { delete bulk; }

    void report_bandwidth(double time) const
    {
        size_t num_coords = (dim)*(dim)*(dim)*10*3;
        std::cerr << "For mesh with " << num_coords << " coordinates acceses.\n";
        double Gbytes = 1.0e-9 * double(sizeof(my_double) * ( num_coords )) ;
        std::cerr << "bandwidth(" << Gbytes * num_repeat / time << " GB/s )" << std::endl;
    }

    void report_bandwidth() const
    {
        double time = 1.0*(end.tv_sec-begin.tv_sec) +
                  1.0e-6*(end.tv_usec-begin.tv_usec);

        report_bandwidth(time);
    }
  
    void start_timer()
    {
        gettimeofday(&begin,NULL);
    }
  
    void stop_timer()
    {
        gettimeofday(&end,NULL);
    }

    stk::mesh::MetaData meta;
    stk::mesh::BulkData* bulk;
    CoordFieldType& centroid;
    CoordFieldType& hostCentroid;
    CoordFieldType* coords;
    int num_repeat;
    size_t dim;
    int choice;
    int teamSize;  
    struct timeval begin, end;
};

#endif

