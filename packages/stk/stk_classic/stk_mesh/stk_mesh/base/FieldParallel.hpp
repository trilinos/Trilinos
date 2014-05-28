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

namespace stk {
namespace mesh {

/**
 * This file contains some helper functions that are part of the Field API.
 * These functions are for making certain parallel operations more convenient.
 */

/** Communicate field data from domain to range.
 *  The fields array must be identical on all processors.
 *  All fields and mesh entities must belong to the same mesh.
 *  If symmetric ( & domain == & range) then from owned to not owned.
 */
void communicate_field_data(
  ParallelMachine machine,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector< const FieldBase *> & fields );

void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields );

/** Communicate field data among shared entities */
void communicate_field_data(
  const BulkData & mesh ,
  const unsigned field_count ,
  const FieldBase * fields[] ,
  CommAll & sparse );

void communicate_field_data_verify_read( CommAll & );

//----------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------
/** Parallel reduction of shared entities' field data.
 *  In the anonymous namespace to avoid redundant link symbols.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( field ) );
 *
 *  where the operations are: sum, max, min
 */
template< class OpField >
void parallel_reduce( const BulkData & mesh ,
                      const OpField  & op )
{
  const FieldBase * fields[1] = { & op.field };

  CommAll sparse ;

  communicate_field_data( mesh, 1, fields, sparse );

  op( mesh.entity_comm() , sparse );

  // For debugging:
  // communicate_field_data_verify_read( sparse );
}

/** Parallel reduction of shared entities' field data.
 *  In the anonymous namespace to avoid redundant link symbols.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( fieldA ) , max( fieldB ) );
 */
template< class OpField1 , class OpField2 >
void parallel_reduce( const BulkData & mesh ,
                      const OpField1 & op1 ,
                      const OpField2 & op2 )
{
  const FieldBase * fields[2] = { & op1.field , & op2.field };

  CommAll sparse ;

  communicate_field_data( mesh, 2, fields, sparse );

  op1( mesh.entity_comm() , sparse );
  op2( mesh.entity_comm() , sparse );

  // For debugging:
  // communicate_field_data_verify_read( sparse );
}

//----------------------------------------------------------------------

template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
struct ParallelReduceField {
  typedef Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;

  const field_type & field ;

  ParallelReduceField( const field_type & f ) : field(f) {}
  ParallelReduceField( const ParallelReduceField & p ) : field(p.field) {}

  void operator()( const std::vector<Entity*> & entity_comm ,
                   CommAll & sparse ) const ;

private:
  ParallelReduceField & operator = ( const ParallelReduceField & );
};

template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
void ParallelReduceField< ReduceOp , Type ,  Tag1,  Tag2,  Tag3 ,
                                     Tag4 ,  Tag5,  Tag6,  Tag7 >::
  operator()( const std::vector<Entity*> & entity_comm ,
              CommAll & sparse ) const
{
  typedef EntityArray< field_type > array_type ;

  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin(); i != entity_comm.end() ; ++i ) {
    Entity & entity = **i ;
    array_type array( field , entity );
    Type * const ptr_beg = array.contiguous_data();
    Type * const ptr_end = ptr_beg + array.size();

    if (ptr_beg == NULL || ptr_end == NULL) continue;

    for ( PairIterEntityComm
          ec = entity.comm() ; ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {

      CommBuffer & b = sparse.recv_buffer( ec->proc );

      for ( Type * ptr = ptr_beg ; ptr < ptr_end ; ++ptr ) {
        Type tmp ;
        b.template unpack<unsigned char>( (unsigned char *)(&tmp), sizeof(Type) );
        ReduceOp( ptr , & tmp );
      }
    }
  }
}

}

//----------------------------------------------------------------------

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Sum<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
sum( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<Sum<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Max<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
max( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<Max<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Min<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
min( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<Min<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

} // namespace mesh
} // namespace stk

#endif

