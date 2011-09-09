/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_MemoryUsage_hpp
#define stk_mesh_MemoryUsage_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/BulkData.hpp>

#include <iosfwd>
#include <vector>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

struct MemoryUsage {
  unsigned num_fields;
  unsigned field_bytes;
  unsigned num_parts;
  unsigned part_bytes;
  std::vector<std::string> entity_rank_names;
  std::vector<unsigned> entity_counts;
  unsigned bytes_per_entity;
  std::vector<unsigned> downward_relation_counts;
  std::vector<unsigned> upward_relation_counts;
  unsigned bytes_per_relation;
  std::vector<unsigned> bucket_counts;
  std::vector<unsigned> bucket_bytes;
  size_t total_bytes;
};

void compute_memory_usage(const BulkData& bulk, MemoryUsage& mem_usage);

void print_memory_usage(const MemoryUsage& mem_usage, std::ostream& os);

//----------------------------------------------------------------------


} // namespace mesh
} // namespace stk

#endif
