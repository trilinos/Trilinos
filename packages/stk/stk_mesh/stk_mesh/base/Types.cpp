#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/parallel/DistributedIndex.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <sstream>
#include <iomanip>

namespace stk { namespace mesh {

void print_max_stk_memory_usage( ParallelMachine parallel, int parallel_rank, std::ostream & out)
{
  std::vector<size_t> memory;
  std::vector<std::string> names;

  // Total
  names.push_back("Total");
  memory.push_back(allocator_memory_usage<void>::peak_memory());
  memory.push_back(allocator_memory_usage<void>::current_memory());
  memory.push_back(allocator_memory_usage<void>::num_allocations());
  memory.push_back(allocator_memory_usage<void>::num_deallocations());

  // Distributed index
  names.push_back("Distributed Index");
  memory.push_back(allocator_memory_usage<parallel::DistributedIndex>::peak_memory());
  memory.push_back(allocator_memory_usage<parallel::DistributedIndex>::current_memory());
  memory.push_back(allocator_memory_usage<parallel::DistributedIndex>::num_allocations());
  memory.push_back(allocator_memory_usage<parallel::DistributedIndex>::num_deallocations());

  // FieldData
  names.push_back("Fields");
  memory.push_back(allocator_memory_usage<FieldDataTag>::peak_memory());
  memory.push_back(allocator_memory_usage<FieldDataTag>::current_memory());
  memory.push_back(allocator_memory_usage<FieldDataTag>::num_allocations());
  memory.push_back(allocator_memory_usage<FieldDataTag>::num_deallocations());

  // Partition
  names.push_back("Partitions");
  memory.push_back(allocator_memory_usage<PartitionTag>::peak_memory());
  memory.push_back(allocator_memory_usage<PartitionTag>::current_memory());
  memory.push_back(allocator_memory_usage<PartitionTag>::num_allocations());
  memory.push_back(allocator_memory_usage<PartitionTag>::num_deallocations());

  // Bucket
  names.push_back("Buckets");
  memory.push_back(allocator_memory_usage<BucketTag>::peak_memory());
  memory.push_back(allocator_memory_usage<BucketTag>::current_memory());
  memory.push_back(allocator_memory_usage<BucketTag>::num_allocations());
  memory.push_back(allocator_memory_usage<BucketTag>::num_deallocations());

  // EntityComm
  names.push_back("Entity Comm");
  memory.push_back(allocator_memory_usage<EntityCommTag>::peak_memory());
  memory.push_back(allocator_memory_usage<EntityCommTag>::current_memory());
  memory.push_back(allocator_memory_usage<EntityCommTag>::num_allocations());
  memory.push_back(allocator_memory_usage<EntityCommTag>::num_deallocations());

  // BucketRelation
  names.push_back("Fixed Connectivity");
  memory.push_back(allocator_memory_usage<BucketRelationTag>::peak_memory());
  memory.push_back(allocator_memory_usage<BucketRelationTag>::current_memory());
  memory.push_back(allocator_memory_usage<BucketRelationTag>::num_allocations());
  memory.push_back(allocator_memory_usage<BucketRelationTag>::num_deallocations());

  // DynamicBucketRelation
  names.push_back("Dynamic Connectivity");
  memory.push_back(allocator_memory_usage<DynamicBucketRelationTag>::peak_memory());
  memory.push_back(allocator_memory_usage<DynamicBucketRelationTag>::current_memory());
  memory.push_back(allocator_memory_usage<DynamicBucketRelationTag>::num_allocations());
  memory.push_back(allocator_memory_usage<DynamicBucketRelationTag>::num_deallocations());

  // AuxRelation
  names.push_back("Aux Connectivity");
  memory.push_back(allocator_memory_usage<AuxRelationTag>::peak_memory());
  memory.push_back(allocator_memory_usage<AuxRelationTag>::current_memory());
  memory.push_back(allocator_memory_usage<AuxRelationTag>::num_allocations());
  memory.push_back(allocator_memory_usage<AuxRelationTag>::num_deallocations());

  // DeletedEntity
  names.push_back("Deleted Entities");
  memory.push_back(allocator_memory_usage<DeletedEntityTag>::peak_memory());
  memory.push_back(allocator_memory_usage<DeletedEntityTag>::current_memory());
  memory.push_back(allocator_memory_usage<DeletedEntityTag>::num_allocations());
  memory.push_back(allocator_memory_usage<DeletedEntityTag>::num_deallocations());

  std::vector<size_t> max_memory(memory.size()*parallel_machine_size(parallel),0);

  MPI_Gather((void*)&memory[0], memory.size(), MPI_LONG_LONG,
	     (void*)&max_memory[0], memory.size(), MPI_LONG_LONG,
	     0, parallel);

  if (parallel_rank == 0) {
    size_t nproc = parallel_machine_size(parallel);
    
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
