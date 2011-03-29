/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Transaction_hpp
#define stk_mesh_Transaction_hpp


#include <set>
#include <vector>
#include <iosfwd>


namespace stk{
namespace mesh {

class Bucket;
class BulkData;
class Entity;
class Part;

    /** \brief Bucket containers used to sort mesh entities by state
     *         so that they can be accessed after modification_end()
     */
    typedef std::vector<Bucket *>                 BucketList;

/** \addtogroup stk_mesh_module
 *  \{
 */


//----------------------------------------------------------------------
/** \brief  Transaction journal of modifications to the bulk data
 *          during a transaction.  Since the modification transaction
 *          guarantees a path independent result of mesh entities when
 *          modification_end() is called, the transaction just notes
 *          the state of altered mesh entities when the transaction
 *          was started.
 */


class Transaction
{
  public:
    /* \brief The following are the variable type and valid values for
     * defining the type of transaction bucket.
     */
    typedef unsigned char                State;
    enum { NOT_IN_TRANSACTION = 0 , MODIFIED = 1, INSERTED = 2, DELETED = 3 };

    /* \brief There are two different types of transactions:
     * incremental and bulk.  This can be set when reseting a
     * transaction.
     */
    enum TransactionType { INCREMENTAL = 1 , BULK = 2 };

  public:

    /** \brief Bucket containers used to sort mesh entities by state
     *         so that they can be accessed after modification_end()
     */
    typedef std::vector<BucketList>               BucketListByType;

    /** \brief Part list for tracking bulk transactions
     */
    typedef std::set<Part *>                      PartSet;


    /** \brief Pretty print the transaction
     */
    std::ostream & print_stream ( std::ostream &os ) const;

    /** \brief Retrieve buckets of a particular type that contain
     * modified entities.
     */
    const BucketList & get_modified_buckets ( unsigned type ) const { return m_modified[type]; }



    /** \brief Retrieve buckets of a particular type that contain
     * deleted entities.
     */
    const BucketList & get_deleted_buckets ( unsigned type ) const { return m_deleted[type]; }

    /** \brief Retrieve buckets of a particular type that contain
     * inserted entities.
     */
    const BucketList & get_inserted_buckets ( unsigned type ) const { return m_inserted[type]; }


    /** \brief Retrieve a part vector of parts whose entities were
     * modified.
     */
    void get_parts_with_modified_entities ( PartVector &pv ) const
    {
      translate_partset_to_partvector ( m_modified_parts , pv );
    }

    /** \brief Retrieve a part vector of parts from which entities
     * were deleted.
     */
    void get_parts_with_deleted_entities ( PartVector &pv ) const
    {
      translate_partset_to_partvector ( m_deleted_parts , pv );
    }

    /** \brief Retrieve a part vector of parts into which entities
     * were inserted.
     */
    void get_parts_with_inserted_entities ( PartVector &pv ) const
    {
      translate_partset_to_partvector ( m_inserted_parts , pv );
    }

  private:

    /** \brief Pretty print helper
     */
    void      print_proc_transaction ( unsigned , std::ostream & ) const;
    /** \brief Pretty print helper
     */
    void  print_transaction ( unsigned , std::ostream & ) const;

    Transaction ();
    Transaction ( BulkData & , TransactionType );
    Transaction & operator = ( const Transaction & );
   ~Transaction ();

    /** \brief Unlike other parts of \ref stk::mesh , Transactions
     *         manage their own memory.  This method is called
     *         during purge() to ensure the bucket containers are
     *         in an appropriate state.  Upon destruction, this
     *         function is not called.
     */
    void allocate_bucket_lists ();

    TransactionType     m_transaction_type;

    BulkData           &m_bulk_data;
    BucketListByType    m_modified;
    BucketListByType    m_deleted;
    BucketListByType    m_inserted;

    PartSet            m_modified_parts;
    PartSet            m_deleted_parts;
    PartSet            m_inserted_parts;

    std::set<Entity *>   m_to_delete;


    void   flush_deletes ();

    /** \brief At modification_begin(), the transaction log is purged.
     * This method will empty the m_modified and m_inserted lists and
     * delete the entities in m_deleted.
     */
    void   reset ( TransactionType type = BULK );
    void   flush ();


    /** \brief Let the transaction log know this entity has been
     * modified in some way.  No relations are followed.
     */
    void   modify_sole_entity ( Entity &e );

    /** \brief Let the transaction log know this entity has been
     * modified in some way.  All entities for which e is in the
     * closure will also be marked modified.
     */
    void   modify_entity ( Entity & e );

    /** \brief Let the transaction log know this entity has been added
     * to the local mesh
     */
    void   insert_entity ( Entity & e );

    /** \brief Let the transaction log know this entity has been
     * removed from the local mesh
     */
    void   delete_entity ( Entity & e );

    /** \brief Helper method to empty m_modified and m_inserted
     */
    void   purge_map ( BucketListByType & );

    /** \brief Helper method to delete entities in m_deleted
     */
    void   purge_and_erase_map ( BucketListByType & );

    /** \brief This method will use the declase bucket method in \ref
     * stk::mesh::Bucket .  This returns the appropriate bucket to
     * store entity information in.
     */
    Bucket *get_unfilled_transaction_bucket ( const unsigned *const , EntityRank , BucketList & , State );

    /** \brief This method will use the declase bucket method in \ref
     * stk::mesh::Bucket .  This returns the appropriate bucket to
     * store entity information in. The function is templated to avoid
     * confusing header inclusions.
     */
    template <class ENTITY>
    Bucket *get_unfilled_transaction_bucket ( const ENTITY &e , BucketList &bm , State s )
     { return get_unfilled_transaction_bucket ( e.bucket().key() , e.entity_rank() , bm , s ); }

    /** \brief This method will add an entity to a transaction bucket
     */
    void   add_entity_to_transaction_bucket ( Entity & , Bucket * );

    /** \brief This method will remove an entity from its current
     * transaction bucket and place it in a new bucket.
     */
    void   swap_entity_between_transaction_buckets ( Entity & , BucketList & , BucketList & , State );

    /** \brief This method will remove an entity from a transaction
     * bucket.
     */
    void   remove_entity_from_bucket ( Entity & , BucketList & );


    /** \brief When a part contains an entity that needs logging, this
     * function will add that part to the particular part set, be it
     * modified, deleted, or inserted
     */
    void   add_parts_to_partset ( Entity & , PartSet & );

    /** \brief This functions creates a PartVector from a PartSet
     */
    void   translate_partset_to_partvector ( const PartSet & , PartVector & ) const;

    // BulkData needs access to modify_entity, insert_entity, and
    // delete_entity.
    friend class BulkData;
    friend std::ostream & operator<< ( std::ostream &os , const Transaction &rhs );

};

inline std::ostream &operator<< ( std::ostream &os , const Transaction &rhs )
{
  return rhs.print_stream ( os );
}

/** \brief Pretty print helper
 */
void  print_bucket_list ( const BucketList & , std::ostream & );
/*
 * \}
 */

}
}



#endif
