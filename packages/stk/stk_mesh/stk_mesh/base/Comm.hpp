#ifndef stk_mesh_Comm_hpp
#define stk_mesh_Comm_hpp

#include <vector>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */
//----------------------------------------------------------------------
/** \brief Sort and unique an EntityProc array.
 */
void sort_unique( std::vector<EntityProc> & );

/** \brief Sanity check locally for non-null, same-mesh, off-processor,
 *  and proper ordering.  Return an error string if a problem.
 */
bool verify( const std::vector<EntityProc> & , std::string & );

/** \brief Find the first entry corresponding to the given entity.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , Entity & );

/** \brief Find the first entry corresponding to the given entity type.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , unsigned );

/** \brief Find the first entry corresponding to the given entity and processor.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , const EntityProc & );

/** \brief Find the first entry corresponding to the given entity.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & , Entity & );

/** \brief Find the first entry corresponding to the given entity and processor.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & , const EntityProc & );

//----------------------------------------------------------------------
/** \brief Sanity check on existing or potential parallel relation information.
 *
 *  If the result is invalid then outputs a string with an explanation.
 *  Symmetric version of verification.
 */
bool comm_verify( ParallelMachine ,
                  const std::vector<EntityProc> & ,
                  std::string & );

/** \brief Sanity check on existing or potential parallel relation information.
 *
 *  If the result is invalid then outputs a string with an explanation.
 *  Asymmetric version of verification.
 */
bool comm_verify( ParallelMachine ,
                  const std::vector<EntityProc> & ,
                  const std::vector<EntityProc> & ,
                  std::string & );

//----------------------------------------------------------------------
/** \brief Global counts for a mesh's entities. */
bool comm_mesh_counts( BulkData & ,
                       std::vector<size_t> & counts ,
                       bool = false );

//----------------------------------------------------------------------
/** \brief Verify that the shared entity values are bit-wise identical */

bool comm_verify_shared_entity_values(
  const BulkData & , unsigned , const FieldBase & );

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

