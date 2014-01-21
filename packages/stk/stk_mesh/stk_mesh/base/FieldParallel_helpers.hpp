/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_FieldParallel_helpers_hpp
#define stk_mesh_FieldParallel_helpers_hpp

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

/// with Selector
template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
struct ParallelReduceField {
  typedef Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;

  const field_type & field ;
  const stk::mesh::Selector * selector ;

  ParallelReduceField( const field_type & f, const stk::mesh::Selector * s = 0 ) : field(f), selector(s) {}
  ParallelReduceField( const ParallelReduceField & p ) : field(p.field), selector(p.selector) {}

  void operator()( const BulkData& mesh, CommAll & sparse ) const ;

private:
  ParallelReduceField & operator = ( const ParallelReduceField & );
};

  template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
void ParallelReduceField< ReduceOp ,  Type ,  Tag1,  Tag2,  Tag3 ,
                                     Tag4 ,  Tag5,  Tag6,  Tag7 >::
operator()(const BulkData& mesh, CommAll & sparse ) const
{
  typedef EntityArray< field_type > array_type ;

  const EntityCommListInfoVector& entity_comm = mesh.comm_list();
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin(); i != entity_comm.end() ; ++i ) {
    Entity entity = i->entity;
    const MeshIndex& mi = mesh.mesh_index(entity);
    if (mesh.is_valid(entity) && (0 == selector || (*selector)(mi.bucket) ) ) {
      if(!is_matching_rank(field, *mi.bucket)) continue;

      typename FieldTraits<field_type>::data_type *array = mesh.field_data(field , *mi.bucket, mi.bucket_ordinal);
      Type * const ptr_beg = array;
      Type * const ptr_end = ptr_beg + field_data_size_per_entity(field, *mi.bucket);

      if (ptr_beg == NULL || ptr_end == NULL) continue;

      for ( PairIterEntityComm
              ec = mesh.entity_comm(i->key); ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {

        CommBuffer & b = sparse.recv_buffer( ec->proc );

        for ( Type * ptr = ptr_beg ; ptr < ptr_end ; ++ptr ) {
          Type tmp ;
          b.template unpack<unsigned char>( reinterpret_cast<unsigned char *>(&tmp), sizeof(Type) );
          ReduceOp( ptr , & tmp );
        }
      }
    }
  }
}

  template< class ReduceOp , class Type>
  struct ParallelReduceFieldBase {
    const FieldBase & field ;
    const stk::mesh::Selector * selector ;

    ParallelReduceFieldBase( const FieldBase & f, const stk::mesh::Selector * s = 0 ) : field(f), selector(s) {}
    ParallelReduceFieldBase( const ParallelReduceFieldBase & p ) : field(p.field), selector(p.selector) {}

    void operator()( const BulkData& mesh, CommAll & sparse ) const ;

  private:
    ParallelReduceFieldBase & operator = ( const ParallelReduceFieldBase & );
  };


template< class ReduceOp , class Type >
void ParallelReduceFieldBase< ReduceOp ,  Type >::
operator()(const BulkData& mesh, CommAll & sparse ) const
{
  const EntityCommListInfoVector& entity_comm = mesh.comm_list();
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin(); i != entity_comm.end() ; ++i ) {
    Entity entity = i->entity;
    const Bucket& bucket = mesh.bucket(entity);

    if(!is_matching_rank(field, bucket)) continue;

    if (mesh.is_valid(entity) && (0 == selector || (*selector)(bucket) ) ) {
      Type * const ptr_beg = reinterpret_cast<Type*>(mesh.field_data(field, entity));
      const unsigned num_scalars_per_entity = field_data_size_per_entity(field, bucket)/sizeof(Type);
      Type * const ptr_end = ptr_beg + num_scalars_per_entity;

      if (ptr_beg == NULL || ptr_end == NULL) continue;

      for ( PairIterEntityComm
              ec = mesh.entity_comm(i->key); ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {

        CommBuffer & b = sparse.recv_buffer( ec->proc );

        for ( Type * ptr = ptr_beg ; ptr < ptr_end ; ++ptr ) {
          Type tmp ;
          b.template unpack<unsigned char>( reinterpret_cast<unsigned char *>(&tmp), sizeof(Type) );
          ReduceOp( ptr , & tmp );
        }
      }
    }
  }
}

//----------------------------------------------------------------------


/// with Selector
template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Sum<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
sum( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f, Selector * selector=0 )
{
  return ParallelReduceField<Sum<1>, Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f, selector );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Max<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
max( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f, Selector * selector=0 )
{
  return ParallelReduceField<Max<1>, Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f, selector );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<Min<1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
min( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f, Selector * selector=0 )
{
  return ParallelReduceField<Min<1>, Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f, selector );
}

template<class Type>
ParallelReduceFieldBase<Sum<1>, Type>
inline
sum(const FieldBase& f, Selector * selector=0 )
{
    return ParallelReduceFieldBase<Sum<1>, Type >(f, selector);
}

template<class Type>
ParallelReduceFieldBase<Max<1>, Type>
inline
max(const FieldBase& f, Selector * selector=0 )
{
    return ParallelReduceFieldBase<Max<1>, Type >(f, selector);
}

template<class Type>
ParallelReduceFieldBase<Min<1>, Type>
inline
min(const FieldBase& f, Selector * selector=0 )
{
    return ParallelReduceFieldBase<Min<1>, Type >(f, selector);
}

} // namespace mesh
} // namespace stk

#endif

