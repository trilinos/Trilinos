#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/parallel/DistributedIndex.hpp>

#include <sstream>

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

  std::vector<size_t> max_memory(memory.size(),0);


  all_reduce_max(parallel, &memory[0], &max_memory[0], memory.size());

  if (parallel_rank == 0) {

    std::ostringstream oss;
    oss << "STK_PROFILE_MEMORY (max across all processes)\n";

    for (size_t tag=0, num_tags=names.size(), i=0; tag < num_tags; ++tag, i += 4) {
      oss << "\n  " << names[tag] << ":\n";
      oss << "             peak = " << max_memory[i+0] << " (" << human_bytes(max_memory[i+0]) << ")" << "\n";
      oss << "          current = " << max_memory[i+1] << " (" << human_bytes(max_memory[i+1]) << ")" << "\n";
      oss << "      allocations = " << max_memory[i+2] << "\n";
      oss << "    deallocations = " << max_memory[i+3] << "\n";
    }

    out << oss.str() << std::endl;
  }
}

}} // namespace stk::mesh

