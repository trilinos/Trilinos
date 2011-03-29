/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>



#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Transaction.hpp>

namespace stk {
namespace mesh {


void print_bucket_list ( const BucketList &bl , std::ostream &os )
{
  BucketList::const_iterator  cur_bucket = bl.begin();
  while ( cur_bucket != bl.end() )
  {
    os << "Bucket key: ";
//This is bad code, no longer works. Re-visit this if/when the
//Transaction class is ever resurrected.
//    for ( unsigned i = 0 ; i <= (*cur_bucket)->m_key[0] ; i++ )
//      os << (*cur_bucket)->m_key[i] << " ";
//    os << "\n";
//    os << "Entities: ";
//    for ( unsigned i = 0 ; i != (*cur_bucket)->m_size ; i++ )
//      os << (*cur_bucket)->m_entities[i]->identifier() << " ";
//    os << "\n-------------------\n";
    cur_bucket++;
  }
}

/**  Disabled to resolve inconsistently destroyed mesh entities.
 *   Revisit when development is driven by a use case.
 */

#if 0

void Transaction::print_transaction ( unsigned type , std::ostream &os ) const
{
  os << "Transaction details for type = " << type << "\n";
  os << "-=-=- Proc " << m_bulk_data.parallel_rank() << " -=-=-\n";
  print_proc_transaction ( type , os );
  os << std::endl;
}

void Transaction::print_proc_transaction ( unsigned type , std::ostream &os ) const
{
  os << "  Modified has " << m_modified[type].size() << " buckets\n"
     << "  Nowhere has " << m_deleted[type].size() << " buckets\n"
     << "  Inserted has " << m_inserted[type].size() << " buckets\n";
  os << " Modified buckets:\n";
  print_bucket_list ( m_modified[type] , os );
  os << " Nowhere buckets:\n";
  print_bucket_list ( m_deleted[type] , os );
  os << " Inserted buckets:\n";
  print_bucket_list ( m_inserted[type] , os );
}

std::ostream & Transaction::print_stream ( std::ostream &os ) const
{
  for ( unsigned i = 0 ; i != m_modified.size() ; i++ )
    print_transaction ( i , os );
  return os;
}

Transaction::Transaction ( BulkData &bd , TransactionType type ) : m_transaction_type(type) ,
                                                                   m_bulk_data ( bd ) ,
                                                                   m_modified ( ) ,
                                                                   m_deleted ( ) ,
                                                                   m_inserted ( )
{
  allocate_bucket_lists ();
}

Transaction::~Transaction ()
{
  reset ();
}



// This allocation works on the assumption that the entity types are
// packed and enumerated from zero.  This assumption is currently safe
// since bulk data makes the same assumption.  Should this ever
// change in bulk data, the same change must occur here.
void Transaction::allocate_bucket_lists ()
{
  m_modified.resize ( m_bulk_data.mesh_meta_data().entity_rank_count() );
  m_deleted.resize ( m_bulk_data.mesh_meta_data().entity_rank_count() );
  m_inserted.resize ( m_bulk_data.mesh_meta_data().entity_rank_count() );
}


// This method will place entity e in the modified buckets.  Unlike
// modify_entity, this will not add the from-entity relations to the
// modified bucket.  This method is invoked on the to-entity of a new
// or changed relation.  If the entity is already in a transaction
// bucket, this method returns with no change
void Transaction::modify_sole_entity ( Entity &e )
{
  // If no bucket is specified yet, then the entity is inserted and
  // should be ignored
  if ( e.m_bucket == 0 ) return;

  // Ignore spurious calls to internal_change_entity_parts
  if ( e.m_trans_bucket != 0 ) return;

  add_parts_to_partset ( e , m_modified_parts );

  // If this is not an incremental transaction, ignore
  if ( m_transaction_type != INCREMENTAL ) return;

  Bucket *transaction_bucket = get_unfilled_transaction_bucket ( e , m_modified[e.entity_rank()] , MODIFIED );

  add_entity_to_transaction_bucket ( e , transaction_bucket );
}

// This method will add e to the modified bucket and every entity for
// which e is directly related to.  This method is invoked an all
// modification other than the to-entity of a new or modified
// relation.  If the entity is already in a transaction bucket, this
// method returns with no change.
void Transaction::modify_entity ( Entity &e )
{

  // If no bucket is specified yet, then the entity is inserted and
  // should be ignored
  if ( e.m_bucket == 0 ) return;

  // Ignore spurious calls to internal_change_entity_parts
  if ( e.m_trans_bucket != 0 ) return;

  add_parts_to_partset ( e , m_modified_parts );

  // If this is not an incremental transaction, ignore
  if ( m_transaction_type != INCREMENTAL ) return;

  Bucket *transaction_bucket = get_unfilled_transaction_bucket ( e , m_modified[e.entity_rank()] , MODIFIED );

  add_entity_to_transaction_bucket ( e , transaction_bucket );

  PairIterRelation current_relation = e.relations();
  while ( current_relation.first != current_relation.second )
  {
    if ( current_relation->entity_rank() > e.entity_rank() )
      modify_sole_entity ( *(current_relation->entity()) );
    current_relation++;
  }
}

// If an entity is removed from this process, it is placed in the
// deleted bucket.  If the entity is in another transaction bucket, it
// is moved to the deleted bucket and the corresponding parts the
// entity is a member of is moved along with it
void Transaction::delete_entity ( Entity &e )
{
  add_parts_to_partset ( e , m_deleted_parts );

  // Determine if entity has already been deleted after insert
  // If so, return
  if ( m_to_delete.find ( &e ) != m_to_delete.end () )
    return;

  // Mark for deletion if the transaction type is BULK
  if ( m_transaction_type == BULK )
    m_to_delete.insert ( &e );

  // If this is not an incremental transaction, ignore
  if ( m_transaction_type != INCREMENTAL ) return;

  if ( e.m_trans_bucket )
  {
    if ( e.transaction_bucket()->transaction_state() == DELETED ) return;
    if ( e.transaction_bucket()->transaction_state() == MODIFIED )
    {
      swap_entity_between_transaction_buckets ( e , m_modified[e.entity_rank()] , m_deleted[e.entity_rank()] , DELETED );
      return;
    }
    if ( e.transaction_bucket()->transaction_state() == INSERTED )
    {
      remove_entity_from_bucket ( e , m_inserted[e.entity_rank()]  );
      // Need to mark this for deletion at reset since it will not be
      // stored anywhere
      m_to_delete.insert ( &e );
      return;
    }
  }

  Bucket *transaction_bucket = get_unfilled_transaction_bucket ( e , m_deleted[e.entity_rank()] , DELETED );
  add_entity_to_transaction_bucket ( e , transaction_bucket );
}

// If an entity is inserted into the mesh, it is placed in the insert
// bucket.  If it is a member of another bucket, this function does
// nothing.  It should be an error to be in another bucket and this
// should be caught elsewhere.
void Transaction::insert_entity ( Entity &e )
{

  if ( e.m_trans_bucket != 0 ) return;

  add_parts_to_partset ( e , m_inserted_parts );

  // If this is not an incremental transaction, ignore
  if ( m_transaction_type != INCREMENTAL ) return;

  Bucket *transaction_bucket = get_unfilled_transaction_bucket ( e , m_inserted[e.entity_rank()] , INSERTED );
  add_entity_to_transaction_bucket ( e , transaction_bucket );
}

// Upon modification_begin() in bulk data, the transaction purges the
// transaction buckets.  This function simply deallocates the bucket
// and removes entities from the transaction bucket.
void Transaction::purge_map ( BucketListByType &buckets )
{
  BucketListByType::iterator cur_type = buckets.begin();
  while ( cur_type != buckets.end() )
  {
    BucketList::iterator cur_bucket = cur_type->begin();
    while ( cur_bucket != cur_type->end() )
    {
      BucketIterator  cur_entity = (*cur_bucket)->begin();
      while ( cur_entity != (*cur_bucket)->end() )
      {
        cur_entity->m_trans_bucket = NULL;
        cur_entity++;
      }
      Bucket::destroy_bucket ( *cur_bucket );
      cur_bucket++;
    }
    cur_type->clear();
    cur_type++;
  }

}

// This method will purge buckets and delete entities from memory.
// This method uses the internal_destroy_entire_bucket method in
// BulkData which does not, in turn, call the delete_entity
// transaction function.
/*
void Transaction::purge_and_erase_map ( BucketListByType &buckets )
{
  BucketListByType::iterator cur_type = buckets.begin();
  while ( cur_type != buckets.end() )
  {
    BucketList::iterator cur_bucket = cur_type->begin();
    while ( cur_bucket != cur_type->end() )
    {
      m_bulk_data.internal_destroy_entire_bucket ( *cur_bucket );
      cur_bucket++;
    }
    cur_type->clear();
    cur_type++;
  }
}

void Transaction::flush_deletes ()
{
  for ( std::set<Entity *>::iterator  cur_del_entity = m_to_delete.begin() ; cur_del_entity != m_to_delete.end() ; cur_del_entity++ )
    m_bulk_data.internal_expunge_entity ( *cur_del_entity );
  m_to_delete.clear ();
}
*/


void Transaction::flush()
{
  purge_map ( m_modified );
  purge_map ( m_inserted );
  purge_map ( m_deleted );

  m_modified_parts.clear();
  m_deleted_parts.clear();
  m_inserted_parts.clear();
}

// Explicity purge and erase memory as needed.
void Transaction::reset ( TransactionType type )
{
  m_transaction_type = type;
  flush();
}

void Transaction::add_parts_to_partset ( Entity &e , PartSet &pl )
{
  PartVector  parts;
  e.bucket().supersets ( parts );

  for ( PartVector::iterator part_iter = parts.begin(); part_iter != parts.end() ; ++part_iter )
    pl.insert ( *part_iter );


}

void Transaction::translate_partset_to_partvector ( const PartSet &in , PartVector &out ) const
{
  out.resize ( in.size() );
  unsigned i = 0;
  for ( PartSet::const_iterator cur_in = in.begin() ; cur_in != in.end() ; cur_in++ )
  {
    out[i] = *cur_in;
    ++i;
  }
}

// The bucket b is assumed to have enough space to add entity e.  This
// is placed on the end of the array and the size is incremented.
void Transaction::add_entity_to_transaction_bucket ( Entity &e , Bucket *b )
{
  b->m_entities[b->m_size] = &e;
  e.m_trans_bucket = b;
  e.m_trans_bucket_ord = b->m_size;
  b->m_size++;
}

// Entity e is removed from ''from'' and placed in ''to''.
void Transaction::swap_entity_between_transaction_buckets ( Entity &e , BucketList &from , BucketList &to , State s )
{
  Bucket *to_bucket = get_unfilled_transaction_bucket ( e , to , s );
  remove_entity_from_bucket ( e , from );
  add_entity_to_transaction_bucket ( e , to_bucket );
}


// This is a wrapper around the Bucket::declare_bucket.  Each bucket
// has no field data and has a copy of the key from the current bucket
// the entity is in.
Bucket *Transaction::get_unfilled_transaction_bucket ( const unsigned * const key , EntityRank type , BucketList &buckets , State s )
{
  Bucket *new_bucket = Bucket::declare_bucket ( m_bulk_data ,
                                                type ,
                                                key[0] - 1 ,
                                                key+1 ,
                                                m_bulk_data.bucket_capacity() ,
                                                std::vector<FieldBase *> () ,
                                                buckets );

  new_bucket->m_transaction_state = s;

  return new_bucket;
}


// This code was copied from bulk data and used to remove entities
// from buckets.  In order to remove an entity, an appropriate entity
// must be found to copy into the hole.  This logic will find the
// appropriate bucket which has an entity to copy into the bucket.
void Transaction::remove_entity_from_bucket ( Entity &e , BucketList &buckets )
{
  Bucket *k = e.m_trans_bucket;
  unsigned i = e.m_trans_bucket_ord;

  Bucket * const first = k->m_key[ *k->m_key ] ? k->m_bucket : k ;
  Bucket * const last  = first->m_bucket ;

  // Only move if not the last entity being removed

  if ( last != k || k->m_size != i + 1 ) {

    // Not the same bucket or not the last entity

    // Copy last entity in last to ik slot i

    Entity * const entity = last->m_entities[ last->m_size - 1 ];

    k->m_entities[i]     = entity ;
    entity->m_trans_bucket     = k ;
    entity->m_trans_bucket_ord = i ;

  }

  --( last->m_size );

  if ( last->m_size != 0 ) {
    last->m_entities[ last->m_size ] = NULL ;
  }
  else {

    // The current 'last' bucket is to be deleted.
    // The previous 'last' bucket becomes the
    // new 'last' bucket in the family:

    std::vector<Bucket*>::iterator ik = lower_bound(buckets, last->m_key);

    ThrowRequireMsg( ik != buckets.end() && last == *ik,
        "Internal failure during removal of entity " << print_entity_key(e) );

    ik = buckets.erase( ik );

    if ( first != last ) { first->m_bucket = *--ik ; }

    Bucket::destroy_bucket( last );
  }

  e.m_trans_bucket = NULL;
}

#endif


}
}

