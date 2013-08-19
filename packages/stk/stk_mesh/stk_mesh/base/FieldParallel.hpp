/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_FieldParallel_hpp
#define stk_mesh_FieldParallel_hpp

//----------------------------------------------------------------------

#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/base/FieldParallel_helpers.hpp>

namespace stk {
namespace mesh {

/**
 * This file contains some helper functions that are part of the Field API.
 */

/** Copy data for the given fields, from owned entities to shared-but-not-owned entities.
 * I.e., shared-but-not-owned entities get an update of the field-data from the owned entity.
*/
void copy_owned_to_shared( const BulkData& mesh,
                           const std::vector< const FieldBase *> & fields );

/** Communicate field data from domain to range.
 *  The fields array must be identical on all processors.
 *  All fields and mesh entities must belong to the same mesh.
 *  If symmetric (domain == range) then from owned to not owned.
 *
 *  Note from ABW: there is no known usage of this function currently (8/19/2013). Possible candidate for deletion.
 *  If you are a developer considering using this function, please contact the STK team.
 */
void communicate_field_data(
  const BulkData& mesh,
  ParallelMachine machine,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector< const FieldBase *> & fields );

/** Send field-data from entities to their ghosts, for a specified 'ghosting'.
 * For entities that are ghosted, this function updates field-data from the
 * original entity to the ghosts.
 */
void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields );

/** Communicate field data among shared entities.
 * This function is a helper function for the parallel_reduce/sum/max/min functions below.
 * This function sends field-data for each shared entity to each sharing processor. When this
 * function is finished, the communicated data is in the 'sparse' CommAll object, not unpacked yet.
 * The data is then unpacked by the calling parallel_* function which also performs the
 * sum/max/min operation on the data before storing it.
 */
void communicate_field_data(
  const BulkData & mesh ,
  const unsigned field_count ,
  const FieldBase * const * fields ,
  CommAll & sparse );

/** debugging function... do we need this?
 */
void communicate_field_data_verify_read( CommAll & );

//----------------------------------------------------------------------
/** Parallel reduction of shared entities' field data.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( field ) );
 *
 *  where the operations are: sum, max, min
 */
template< class OpField >
inline
void parallel_reduce( const BulkData & mesh ,
                      const OpField  & op )
{
  const FieldBase * fields[1] = { & op.field };

  CommAll sparse ;

  communicate_field_data( mesh, 1, fields, sparse );

  op( mesh, sparse );

  // For debugging:
  // communicate_field_data_verify_read( sparse );
}

/** Parallel reduction of shared entities' field data.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( fieldA ) , max( fieldB ) );
 */
template< class OpField1 , class OpField2 >
inline
void parallel_reduce( const BulkData & mesh ,
                      const OpField1 & op1 ,
                      const OpField2 & op2 )
{
  const FieldBase * fields[2] = { & op1.field , & op2.field };

  CommAll sparse ;

  communicate_field_data( mesh, 2, fields, sparse );

  op1( mesh , sparse );
  op2( mesh , sparse );

  // For debugging:
  // communicate_field_data_verify_read( sparse );
}

//----------------------------------------------------------------------

/** Sum (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */
inline
void parallel_sum(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  if (fields.empty()) return;

  const FieldBase *const* fieldsPtr = &fields[0];
  CommAll sparse;
  communicate_field_data(mesh, fields.size(), fieldsPtr, sparse);

  for(size_t i=0; i<fields.size(); ++i) {
      if (fields[i]->type_is<double>()) {
        sum<double>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<float>()) {
        sum<float>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<int>()) {
        sum<int>(*fields[i])(mesh, sparse);
      }
      else {
        ThrowRequireMsg(false, "Error, parallel_sum only operates on fields of type double, float or int.");
      }
  }
}

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
inline
void parallel_max(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  if (fields.empty()) return;

  const FieldBase *const* fieldsPtr = &fields[0];
  CommAll sparse;
  communicate_field_data(mesh, fields.size(), fieldsPtr, sparse);

  for(size_t i=0; i<fields.size(); ++i) {
      if (fields[i]->type_is<double>()) {
        max<double>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<float>()) {
        max<float>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<int>()) {
        max<int>(*fields[i])(mesh, sparse);
      }
      else {
        ThrowRequireMsg(false, "Error, parallel_max only operates on fields of type double, float or int.");
      }
  }
}

/** Communicate and take the minimum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (minimum) field values on each sharing proc.
 */
inline
void parallel_min(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  if (fields.empty()) return;

  const FieldBase *const* fieldsPtr = &fields[0];
  CommAll sparse;
  communicate_field_data(mesh, fields.size(), fieldsPtr, sparse);

  for(size_t i=0; i<fields.size(); ++i) {
      if (fields[i]->type_is<double>()) {
        min<double>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<float>()) {
        min<float>(*fields[i])(mesh, sparse);
      }
      else if (fields[i]->type_is<int>()) {
        min<int>(*fields[i])(mesh, sparse);
      }
      else {
        ThrowRequireMsg(false, "Error, parallel_min only operates on fields of type double, float or int.");
      }
  }
}

} // namespace mesh
} // namespace stk

#endif

