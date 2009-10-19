#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>


namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

void print_entry( std::ostream & msg , const EntityProc & e )
{
  msg << "( " ;
  if ( e.first == NULL ) {
    msg << "NULL" ;
  }
  else {
    const MetaData & meta_data = e.first->bucket().mesh().mesh_meta_data();
    print_entity_key( msg , meta_data , e.first->key() );
  }
  msg << " , " << e.second << " )" ;
}

}

//----------------------------------------------------------------------

struct LessEntityProc {
  LessEntityProc() {}

  bool operator()( const EntityProc & lhs , const EntityProc & rhs ) const
  {
    const EntityKey lhs_key = lhs.first->key();
    const EntityKey rhs_key = rhs.first->key();
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
  }

  bool operator()( const EntityProc & lhs , const Entity & rhs ) const
  {
    const EntityKey lhs_key = lhs.first->key();
    const EntityKey rhs_key = rhs.key();
    return lhs_key < rhs_key ;
  }

  bool operator()( const EntityProc & lhs , const unsigned rhs ) const
  {
    const unsigned lhs_type = lhs.first->entity_type();
    return lhs_type < rhs ;
  }
};

struct EqualEntityProc {
  bool operator()( const EntityProc & lhs , const EntityProc & rhs ) const
  { return lhs.first == rhs.first && lhs.second == rhs.second ; }
};

//----------------------------------------------------------------------

void sort_unique( std::vector<EntityProc> & v )
{
  std::vector<EntityProc>::iterator i = v.begin();
  std::vector<EntityProc>::iterator e = v.end();

  std::sort( i , e , LessEntityProc() );
  i = std::unique( i , e , EqualEntityProc() );
  v.erase( i , e );
}

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & v , unsigned t )
{
  const std::vector<EntityProc>::const_iterator i = v.begin();
  const std::vector<EntityProc>::const_iterator e = v.end();
  return lower_bound( i , e , t , LessEntityProc() );
}

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & v , Entity & m )
{
  const std::vector<EntityProc>::const_iterator i = v.begin();
  const std::vector<EntityProc>::const_iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & v , const EntityProc & m )
{
  const std::vector<EntityProc>::const_iterator i = v.begin();
  const std::vector<EntityProc>::const_iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & v , Entity & m )
{
  const std::vector<EntityProc>::iterator i = v.begin();
  const std::vector<EntityProc>::iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & v , const EntityProc & m )
{
  const std::vector<EntityProc>::iterator i = v.begin();
  const std::vector<EntityProc>::iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

//----------------------------------------------------------------------

#if 0

void procs( const std::vector<EntityProc> & v ,
            const std::vector< Entity * > & m ,
            std::vector<unsigned> & p )
{
  typedef std::vector<EntityProc>::const_iterator Iter ;

  p.clear();

  std::vector<Iter> tmp ; tmp.reserve( m.size() );

  for ( unsigned j = 0 ; j < tmp.size() ; ++j ) {
    tmp.push_back( lower_bound( v , *m[j] ) );
  }

  if ( ! tmp.empty() ) {
    const Iter ve = v.end();
    Entity * const e0 = m[0] ;

    for ( ; tmp[0] != ve && tmp[0]->first == e0 ; ++tmp[0] ) {
      const unsigned proc = tmp[0]->second ;
      bool flag = true ;

      // Iterate remaining up to proc

      for ( unsigned j = 1 ; j < tmp.size() ; ++j ) {
        Entity * const ej = m[j] ;
        for ( ; tmp[j] != ve &&
                tmp[j]->first == ej &&
                tmp[j]->second < proc ; ++tmp[j] );
        if ( tmp[j] == ve ||
             tmp[j]->first != ej ||
             tmp[j]->second != proc ) {
          flag = false ;
        }
      }

      if ( flag ) { p.push_back( proc ); }
    }
  }
}

#endif

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool verify( const std::vector<EntityProc> & v , std::string & msg )
{
  const LessEntityProc less_op ;

  bool result = true ;

  if ( ! v.empty() ) {

    const std::vector<EntityProc>::const_iterator e = v.end();
    std::vector<EntityProc>::const_iterator i ;

    std::ostringstream os ;

    for ( i = v.begin() ; result && i != e ; ++i ) {
      if ( ! ( result = i->first != NULL ) ) {
        msg.append( "Contains NULL entries" );
      }
    }

    if ( result ) {
      BulkData   & M = v[0].first->bucket().mesh();
      const unsigned p_size = M.parallel_size();

      for ( i = v.begin() ; result && i != e ; ++i ) {
        if ( ! ( result = & M == & i->first->bucket().mesh() ) ) {
          msg.append( "Contains entries from different meshes" );
        }
      }

      for ( i = v.begin() ; result && i != e ; ++i ) {
        const unsigned p = i->second ;
        if ( ! ( result = p < p_size ) ) {
          os << "Contains entry with bad processor " ;
          print_entry( os , *i );
          msg.append( os.str() );
        }
      }

      EntityProc old( v[0] ) ;

      for ( i = v.begin() ; result && ++i != e ; ) {
        if ( ! ( result = less_op( old , *i ) ) ) {
          os << "Contains out-of-order entries " ;
          print_entry( os , old );
          print_entry( os , *i );
          msg.append( os.str() );
        }
        old = *i ;
      }
    }
  }
  return result ;
}


// REFACTOR: The alignment and preallocation of buffers doesn't mix.
//           The ParallelComm alignment and packing are the problem here.
//           See "x" variable hacks below
//           The Marshall code simply packs tightly.
//
//           The buffer sizing and packing originally matched properly.
//           A recently introduced hacking of this code to make the
//           EntityKey refactor work broke this correctness.
//           This function should have a proper refactoring
//           instead of the hacking that it was given.
//
bool comm_verify( ParallelMachine comm ,
                  const std::vector<EntityProc> & v ,
                  std::string & msg )
{
  static const char method[] = "stk::mesh::comm_verify[symmetric]" ;

  // Verify local ordering, the result flag is parallel inconsistent

  bool result = verify( v , msg );

  // Verify parallel consistency
  // Communicate asymmetric and compare keys and owners

  const int p_rank = parallel_machine_rank( comm );
  const int p_size = parallel_machine_size( comm );

  std::ostringstream os ;

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( p_size , zero );

  CommAll comm_sparse(comm) ;

  const std::vector<EntityProc>::const_iterator i_end = v.end();
  std::vector<EntityProc>::const_iterator i ;

  if ( result ) {
    for ( i = v.begin() ; i != i_end ; ++i ) {
      Entity & e = * i->first ;
      const unsigned  proc  = i->second ;

      CommBuffer & buf = comm_sparse.send_buffer( proc );
      buf.pack(e.key());
      buf.pack(e.owner_rank());
    }
  }

  {
    result = ! comm_sparse.allocate_buffers( p_size / 4 );
  }

  // The result flag is now parallel consistent

  if ( result ) {
    // Fill buffer
    for ( i = v.begin() ; i != i_end ; ++i ) {
      Entity & e = * i->first ;
      const unsigned  proc  = i->second ;

      CommBuffer & buf = comm_sparse.send_buffer( proc );
      buf.pack(e.key());
      buf.pack(e.owner_rank());
    }

    comm_sparse.communicate();

    // Verify symmetry of sizes
//     for ( int j = 0 ; result && j < p_size ; ++j ) {
//       const size_t nrecv = comm_sparse.recv_buffer( j ).remaining();
//       if ( nrecv != send_size[j] ) {
//         os << method ;
//         os << " parallel inconsistency, P" << p_rank ;
//         os << " expected from P" << j ;
//         os << " " << send_size[j] ;
//         os << " but received " ;
//         os << nrecv ;
//         os << " instead" ;
//         result = false ;
//         msg.append( os.str() );
//       }
//     }

    // Verify symmetry of content
    for ( i = v.begin() ; result && i != i_end ; ++i ) {
      Entity & e = * i->first ;

      const EntityKey this_key   = e.key();
      const unsigned        this_owner = e.owner_rank();
      const unsigned        proc  = i->second ;

      CommBuffer & buf = comm_sparse.recv_buffer( proc );
      EntityKey entity_key;
      unsigned owner = ~0u ;
      
      buf.unpack(entity_key);
      buf.unpack( owner );
      
      if ( this_key   != entity_key ||
           this_owner != owner ) {
        const MetaData & meta_data = e.bucket().mesh().mesh_meta_data();
        os << method ;
        os << " parallel inconsistency, P" << p_rank << " has " ;
        print_entity_key( os , meta_data , this_key );
        os << "].owner(P" << this_owner ;
        os << ") versus " ;
        print_entity_key( os , meta_data , entity_key );
        os << "].owner(P" << owner ;
        os << ")" ;
        result = false ;
        msg.append( os.str() );
      }
    }

    // The result flag is now parallel inconsistent

    {
      unsigned flag = result ;
      all_reduce( comm , ReduceMin<1>( & flag ) );
      result = flag ;
    }
  }

  return result ;
}

//----------------------------------------------------------------------

bool comm_verify( ParallelMachine comm ,
                  const std::vector<EntityProc> & send ,
                  const std::vector<EntityProc> & recv ,
                  std::string & msg )
{
  static const char method[] = "stk::mesh::comm_verify[asymmetric]" ;

  // Verify local ordering:

  bool result = verify( send , msg );

  const int p_rank = parallel_machine_rank( comm );
  const int p_size = parallel_machine_size( comm );

  std::ostringstream os ;

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( p_size , zero );
  std::vector<unsigned> recv_size( p_size , zero );

  std::vector<EntityProc>::const_iterator i ;

  if ( result ) {
    for ( i = send.begin() ; i != send.end() ; ++i ) {
      ++( send_size[ i->second ] );
    }
  }

  {
    unsigned msg_max = 0 ;
    unsigned * const p_send = & send_size[0] ;
    unsigned * const p_recv = & recv_size[0] ;
    result = ! comm_sizes( comm , p_size / 4 , msg_max ,
                           p_send , p_recv , ! result );
  }

  // Result flag is now parallel consistent

  if ( result ) {
    send_size.assign( p_size , zero );

    for ( i = recv.begin() ; i != recv.end() ; ++i ) {
      ++( send_size[ i->second ] );
    }

    for ( int j = 0 ; result && j < p_size ; ++j ) {
      const unsigned nrecv = recv_size[j] ;
      if ( nrecv != send_size[j] ) {
        os << method ;
        os << " parallel inconsistency, P" << p_rank ;
        os << " expected from P" << j ;
        os << " " << send_size[j] ;
        os << " but received " ;
        os << nrecv ;
        os << " instead" ;
        result = false ;
        msg.append( os.str() );
      }
    }

    // The result flag is now parallel inconsitent

    {
      unsigned flag = result ;
      all_reduce( comm , ReduceBitAnd<1>( & flag ) );
      result = flag ;
    }
  }
  return result ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_mesh_counts( BulkData & M ,
                       std::vector<size_t> & counts ,
                       bool local_flag )
{
  const size_t zero = 0 ;

  // Count locally owned entities

  const MetaData & S = M.mesh_meta_data();
  const unsigned entity_type_count = S.entity_type_count();
  const size_t   comm_count        = entity_type_count + 1 ;

  std::vector<size_t> local(  comm_count , zero );
  std::vector<size_t> global( comm_count , zero );

  ParallelMachine comm = M.parallel();
  Part & owns = S.locally_owned_part();

  for ( unsigned i = 0 ; i < entity_type_count ; ++i ) {
    const std::vector<Bucket*> & ks = M.buckets( i );

    std::vector<Bucket*>::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( has_superset( **ik , owns ) ) {
        local[i] += (*ik)->size();
      }
    }
  }

  local[ entity_type_count ] = local_flag ;

  all_reduce_sum( comm , & local[0] , & global[0] , comm_count );

  counts.assign( global.begin() , global.begin() + entity_type_count );

  return 0 < global[ entity_type_count ] ;
}

//----------------------------------------------------------------------

namespace {

std::pair< std::vector<EntityProc>::const_iterator ,
           std::vector<EntityProc>::const_iterator >
span( const std::vector<EntityProc> & v , unsigned entity_type )
{
  const unsigned t1 = entity_type ;
  const unsigned t2 = t1 + 1 ;

  std::pair< std::vector<EntityProc>::const_iterator ,
             std::vector<EntityProc>::const_iterator > result ;

  result.first  = lower_bound(v,t1);
  result.second = lower_bound(v,t2);

  return result ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_verify_shared_entity_values(
  const BulkData & M , unsigned t , const FieldBase & f )
{
  const unsigned parallel_size = M.parallel_size();

  const unsigned max_size = f.max_size(t) * f.data_traits().size_of ;

  const std::pair< std::vector<EntityProc>::const_iterator ,
                   std::vector<EntityProc>::const_iterator >
    shares = span( M.shared_entities() , t );

  std::vector<EntityProc>::const_iterator ic ;

  ParallelMachine comm = M.parallel();

  CommAll sparse ;

  {
    const unsigned zero = 0 ;
    std::vector<unsigned> comm_size( parallel_size , zero );

    for ( ic = shares.first ; ic != shares.second ; ++ic ) {
      comm_size[ ic->second ] += max_size ;
    }

    const unsigned * const p_size = & comm_size[0] ;

    sparse.allocate_buffers( comm, parallel_size / 4 , p_size, p_size );
  }

  std::vector<unsigned char> scratch( max_size );

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    const Bucket & k = e.bucket();
    const unsigned this_size = field_data_size( f , k );
    CommBuffer & b = sparse.send_buffer( ic->second );

    unsigned char * ptr = reinterpret_cast<unsigned char *>( field_data( f , e ) );

    b.pack<unsigned char>( ptr , this_size );
    b.skip<unsigned char>( max_size - this_size );
  }

  sparse.communicate();

  unsigned ok = 1 ;

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    const Bucket & k = e.bucket();
    const unsigned this_size = field_data_size( f , k );
    CommBuffer & b = sparse.recv_buffer( ic->second );

    unsigned char * ptr = reinterpret_cast<unsigned char *>( field_data( f , e ) );
    unsigned char * scr = & scratch[0] ;

    b.unpack<unsigned char>( scr , this_size );
    b.skip<unsigned char>( max_size - this_size );

    // Compare data and scratch
    for ( unsigned j = 0 ; ok && j < this_size ; ++j ) {
      ok = ptr[j] == scr[j] ;
    }
  }

  all_reduce( comm , ReduceMin<1>( & ok ) );

  return ok ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

