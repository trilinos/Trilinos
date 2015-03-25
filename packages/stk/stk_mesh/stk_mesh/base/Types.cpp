// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 


#include <stk_mesh/base/Types.hpp>
#include <algorithm>                    // for sort
#include <iomanip>                      // for operator<<, setw
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <stk_util/util/human_bytes.hpp>  // for human_bytes
#include <string>                       // for string, operator<<
namespace stk { namespace parallel { class DistributedIndex; } }


namespace stk { namespace mesh {

namespace {

template <typename Tag>
void assemble_data(std::vector<std::string>& names,
                   std::vector<size_t>& memory,
                   std::string const& label,
                   size_t* peak_sum = NULL,
                   size_t* curr_sum = NULL,
                   size_t* alloc_sum = NULL,
                   size_t* dealloc_sum = NULL)
{
  names.push_back(label);

  memory.push_back(allocator_memory_usage<Tag>::peak_memory());
  if (peak_sum != NULL) {
    *peak_sum += memory.back();
  }

  memory.push_back(allocator_memory_usage<Tag>::current_memory());
  if (curr_sum != NULL) {
    *curr_sum += memory.back();
  }

  memory.push_back(allocator_memory_usage<Tag>::num_allocations());
  if (alloc_sum != NULL) {
    *alloc_sum += memory.back();
  }

  memory.push_back(allocator_memory_usage<Tag>::num_deallocations());
  if (dealloc_sum != NULL) {
    *dealloc_sum += memory.back();
  }
}

}

void print_max_stk_memory_usage( ParallelMachine parallel, int parallel_rank, std::ostream & out)
{
  std::vector<size_t> memory;
  std::vector<std::string> names;

  // Total
  assemble_data<void>(names, memory, "Total");

  assemble_data<SelectorMapTag>(names, memory, "Get Buckets Memoization");

  // Distributed index
  assemble_data<parallel::DistributedIndex>(names, memory, "Distributed Index");

  // FieldData
  assemble_data<FieldDataTag>(names, memory, "Fields");

  // Partition
  assemble_data<PartitionTag>(names, memory, "Partitions");

  // Bucket
  assemble_data<BucketTag>(names, memory, "Buckets");

  // EntityComm
  assemble_data<EntityCommTag>(names, memory, "Entity Comm");

  // BucketRelation
  assemble_data<BucketRelationTag>(names, memory, "Fixed Connectivity");

  // Dynamic Connectivity for specific ranks
  size_t total_dyn_peak   = 0;
  size_t total_dyn_curr   = 0;
  size_t total_dyn_allocs = 0;
  size_t total_dyn_dels   = 0;

  assemble_data<DynamicBucketNodeRelationTag>(names, memory, "Dynamic Node Connectivity", &total_dyn_peak, &total_dyn_curr, &total_dyn_allocs, &total_dyn_dels);

  assemble_data<DynamicBucketEdgeRelationTag>(names, memory, "Dynamic Edge Connectivity", &total_dyn_peak, &total_dyn_curr, &total_dyn_allocs, &total_dyn_dels);

  assemble_data<DynamicBucketFaceRelationTag>(names, memory, "Dynamic Face Connectivity", &total_dyn_peak, &total_dyn_curr, &total_dyn_allocs, &total_dyn_dels);

  assemble_data<DynamicBucketElementRelationTag>(names, memory, "Dynamic Element Connectivity", &total_dyn_peak, &total_dyn_curr, &total_dyn_allocs, &total_dyn_dels);

  assemble_data<DynamicBucketOtherRelationTag>(names, memory, "Dynamic Other Connectivity", &total_dyn_peak, &total_dyn_curr, &total_dyn_allocs, &total_dyn_dels);

  // DynamicBucketRelation
  names.push_back("Dynamic Total Connectivity");
  memory.push_back(total_dyn_peak);
  memory.push_back(total_dyn_curr);
  memory.push_back(total_dyn_allocs);
  memory.push_back(total_dyn_dels);

  // AuxRelation
  assemble_data<AuxRelationTag>(names, memory, "Aux Connectivity");

  // DeletedEntity
  assemble_data<DeletedEntityTag>(names, memory, "Deleted Entities");

  assemble_data<VolatileFastSharedCommMapTag>(names, memory, "Volatile Fast Shared Comm Map");

#if defined ( STK_HAS_MPI )
  std::vector<size_t> max_memory(memory.size()*parallel_machine_size(parallel),0);

  MPI_Gather(static_cast<void*>(&memory[0]), memory.size(), MPI_LONG_LONG,
	           static_cast<void*>(&max_memory[0]), memory.size(), MPI_LONG_LONG,
	           0, parallel);
#else
  std::vector<size_t> max_memory(memory);
#endif
  
  const int nproc = parallel_machine_size(parallel);

  if (parallel_rank == 0) {

    std::ostringstream oss;
    oss << "STK_PROFILE_MEMORY (max, min, median across all processes)\n";

    for (size_t tag=0, num_tags=names.size(), j=0; tag < num_tags; ++tag, j += 4) {
      std::vector<size_t> v0(nproc), v1(nproc), v2(nproc), v3(nproc);
      for (size_t i=0, ii=0; i < nproc*memory.size(); i+=memory.size(), ii++) {
	v0[ii] = max_memory[j+i+0];
	v1[ii] = max_memory[j+i+1];
	v2[ii] = max_memory[j+i+2];
	v3[ii] = max_memory[j+i+3];
      }
      std::sort(v0.begin(), v0.end());
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      std::sort(v3.begin(), v3.end());

      oss << "\n  " << names[tag] << ":\n";
      oss << "             peak = "
	  << std::setw(10) << v0[nproc-1] << "  "
	  << std::setw(10) << v0[0]       << "  "
	  << std::setw(10) << v0[nproc/2] << "\t("
	  << human_bytes(v0[nproc-1]) << "  "
	  << human_bytes(v0[0]) << "  "
	  << human_bytes(v0[nproc/2]) << ")\n";
      oss << "          current = "
	  << std::setw(10) << v1[nproc-1] << "  "
	  << std::setw(10) << v1[0]       << "  "
	  << std::setw(10) << v1[nproc/2] << "\t("
	  << human_bytes(v1[nproc-1]) << "  "
	  << human_bytes(v1[0]) << "  "
	  << human_bytes(v1[nproc/2]) << ")\n";
      oss << "      allocations = " << std::setw(10)
	  << v2[nproc-1] << "  " << std::setw(10)
	  << v2[0]       << "  " << std::setw(10)
	  << v2[nproc/2] << "\n";
      oss << "    deallocations = " << std::setw(10)
	  << v3[nproc-1] << "  " << std::setw(10)
	  << v3[0]       << "  " << std::setw(10)
	  << v3[nproc/2] << "\n";
    }
    out << oss.str() << std::endl;
  }
}

}} // namespace stk::mesh
