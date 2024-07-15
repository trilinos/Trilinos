// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Core.hpp"

#include<ostream>
#include<string>

extern "C" void kokkosp_begin_fence(const char* name, const uint32_t deviceId,
                                    uint64_t* handle)
{
  std::cout << "Begin fence: " << std::string(name) << std::endl;
}

extern "C" void kokkosp_end_fence(const uint64_t handle)
{
  std::cout << "End fence" << std::endl;
}

extern "C" void kokkosp_begin_deep_copy(Kokkos::Profiling::SpaceHandle dst_handle, const char* dst_name, const void* dst_ptr,
                                        Kokkos::Profiling::SpaceHandle src_handle, const char* src_name, const void* src_ptr,
                                        uint64_t size)
{
  std::cout << "begin deep_copy: src=" << std::string(src_name) << ", dst=" << std::string(dst_name) << std::endl; 
}

extern "C" void kokkosp_end_deep_copy() {
  std::cout << "end deep_copy" << std::endl; 
}

extern "C" void kokkosp_allocate_data(const Kokkos::Profiling::SpaceHandle space, const char* label, const void* const ptr, const uint64_t size)
{
  std::cout << "allocate view: " << std::string(label) << std::endl;
}

extern "C" void kokkosp_deallocate_data(const Kokkos::Profiling::SpaceHandle space, const char* label, const void* const ptr, const uint64_t size)
{
  std::cout << "deallocate view: " << std::string(label) << std::endl;
}

// extern "C" void kokkosp_request_tool_settings(const uint32_t /*num_settings*/,Kokkos_Tools_ToolSettings* settings)
// {
//   std::cout << "setting tool settings to false: " << std::endl;
//   settings->requires_global_fencing = false;
// }

// TEUCHOS_UNIT_TEST(kokkos, PartitioningExecSpaces)
int main (int argc, char** argv)
{
  std::cout << "Starting MAIN\n";

  {
    std::cout << "Begin Setting Kokkos profiling\n";
    using namespace Kokkos::Profiling;
    using namespace Kokkos::Profiling::Experimental;

    // Kokkos::Tools::Experimental::set_request_tool_settings_callback((requestToolSettingsFunction)(kokkosp_request_tool_settings));
    Kokkos::Tools::Experimental::
      set_request_tool_settings_callback([](const uint32_t /*num_settings*/,
                                            Kokkos::Tools::Experimental::ToolSettings* settings)
                                         {
                                           std::cout << "In tools setting, PRE  fencing=" << settings->requires_global_fencing << "\n";
                                           settings->requires_global_fencing = true;
                                           std::cout << "In tools setting, POST fencing=" << settings->requires_global_fencing << "\n";
                                         });

    Kokkos::Tools::Experimental::set_begin_fence_callback((beginFunction)(kokkosp_begin_fence));
    Kokkos::Tools::Experimental::set_end_fence_callback((endFunction)(kokkosp_end_fence));
    Kokkos::Tools::Experimental::set_begin_deep_copy_callback((beginDeepCopyFunction)(kokkosp_begin_deep_copy));
    Kokkos::Tools::Experimental::set_end_deep_copy_callback((endDeepCopyFunction)(kokkosp_end_deep_copy));
    Kokkos::Tools::Experimental::set_allocate_data_callback((allocateDataFunction)(kokkosp_allocate_data));
    Kokkos::Tools::Experimental::set_deallocate_data_callback((deallocateDataFunction)(kokkosp_deallocate_data));
    std::cout << "End Setting Kokkos profiling\n";
  }

  Kokkos::initialize(argc,argv);

  {
    std::vector<Kokkos::DefaultExecutionSpace> streams;
    if (Kokkos::DefaultExecutionSpace().concurrency() >= 3) {
      std::cout << "Using partition_space, concurrency=" << Kokkos::DefaultExecutionSpace().concurrency() << std::endl;
      streams = Kokkos::Experimental::partition_space(Kokkos::DefaultExecutionSpace(),1,1,1);
    }
    else {
      std::cout << "NOT using partition_space, concurrency=" << Kokkos::DefaultExecutionSpace().concurrency() << std::endl;
      for (int i=0; i < 3; ++i)
        streams.push_back(Kokkos::DefaultExecutionSpace());
    }

    std::cout << "Default exec space=" << Kokkos::DefaultExecutionSpace().name() << std::endl;
    Kokkos::DefaultExecutionSpace().print_configuration(std::cout,true);
    std::cout << "streams[0]=" << streams[0].name() << std::endl;
    streams[0].print_configuration(std::cout,true);
    std::cout << "streams[1]=" << streams[1].name() << std::endl;
    streams[1].print_configuration(std::cout,true);
    std::cout << "streams[2]=" << streams[2].name() << std::endl;
    streams[2].print_configuration(std::cout,true);

    const int N = 10000;
    std::cout << "Starting a alloc from main\n";
    Kokkos::View<double*> a(Kokkos::view_alloc(streams[0],"a"),N);
    std::cout << "Starting b alloc from main\n";
    Kokkos::View<double*> b(Kokkos::view_alloc(streams[1],"b"),N);
    std::cout << "Starting c alloc from main\n";
    Kokkos::View<double*> c(Kokkos::view_alloc(streams[2],"c"),N);

    Kokkos::deep_copy(streams[0],a,1.0);
    Kokkos::deep_copy(streams[1],b,2.0);
    Kokkos::deep_copy(streams[2],c,3.0);

    streams[0].fence("fence after deep_copy");
    streams[1].fence("fence after deep_copy");
    streams[2].fence("fence after deep_copy");

    auto policy = Kokkos::RangePolicy<>(streams[2],0,N);
    Kokkos::parallel_for("evaluate c",policy,KOKKOS_LAMBDA(const int i){
      c(i) += a(i) + b(i);
    });

    streams[2].fence("fence after parallel for");
  }

  Kokkos::finalize();

  std::cout << "Exiting MAIN\n";
  return 0;
}
