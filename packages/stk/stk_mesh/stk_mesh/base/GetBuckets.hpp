/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_GetBucket_hpp
#define stk_mesh_GetBucket_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Iterators.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------

/* \brief  Select buckets from the input to the output. Buckets in the input
 *         vector will be placed in the output vector if the bucket is
 *         selected by the selector argument.
 *         On entry, the output vector is cleared before being filled with
 *         selected buckets.
 */
void get_buckets( const Selector & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< Bucket * > & output );

/**
 * For all the buckets in the mesh, return a range of selected buckets.
 */
AllSelectedBucketsRange get_buckets( const Selector & selector, const BulkData& mesh );

/**
 * Return a range of all buckets in the mesh
 */
AllBucketsRange get_buckets( const BulkData& mesh );

/**
 * Return a range of buckets for a certain rank.
 */
AllBucketsRange get_buckets( EntityRank entity_rank, const BulkData& mesh );

/**
 * Return a range of selected buckets within a bucket range
 */
AllSelectedBucketsRange get_buckets( const Selector & selector, const AllBucketsRange& range);

/* \brief  Get the parts from the union part vector that the bucket is
 *         contained in. 
 */
void get_involved_parts( const PartVector & union_parts,
                         const Bucket & candidate,
                         PartVector & involved_parts);

//----------------------------------------------------------------------
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

