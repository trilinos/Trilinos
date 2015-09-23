/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>


namespace stk_classic {
namespace mesh {

void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = BulkData::get(ghosts);
  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( std::vector<Entity*>::const_iterator
        i =  mesh.entity_comm().begin() ;
        i != mesh.entity_comm().end() ; ++i ) {
    Entity & e = **i ;
    const bool owned = e.owner_rank() == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;
      e_size += field_data_size( f , e );
    }

    for ( PairIterEntityComm ec = e.comm() ; ! ec.empty() ; ++ec ) {
      if ( ghosts.ordinal() == ec->ghost_id ) {
        if ( owned ) {
          send_size[ ec->proc ] += e_size ;
        }
        else {
          recv_size[ ec->proc ] += e_size ;
        }
      }
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const s_size = & send_size[0] ;
    const unsigned * const r_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), parallel_size / 4 , s_size, r_size);
  }

  // Send packing:

  for ( std::vector<Entity*>::const_iterator
        i =  mesh.entity_comm().begin() ;
        i != mesh.entity_comm().end() ; ++i ) {
    Entity & e = **i ;
    if ( e.owner_rank() == parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(field_data( f , e ));

          for ( PairIterEntityComm ec = e.comm() ; ! ec.empty() ; ++ec ) {

            if ( ghosts.ordinal() == ec->ghost_id ) {
              CommBuffer & b = sparse.send_buffer( ec->proc );
              b.pack<unsigned char>( ptr , size );
            }
          }
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( std::vector<Entity*>::const_iterator
        i =  mesh.entity_comm().begin() ;
        i != mesh.entity_comm().end() ; ++i ) {
    Entity & e = **i ;
    if ( e.owner_rank() != parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(field_data( f , e ));

          for ( PairIterEntityComm ec = e.comm() ; ! ec.empty() ; ++ec ) {

            if ( ghosts.ordinal() == ec->ghost_id ) {
              CommBuffer & b = sparse.recv_buffer( ec->proc );
              b.unpack<unsigned char>( ptr , size );
            }
          }
        }
      }
    }
  }
}

// Heterogeneity?

void communicate_field_data(
  ParallelMachine machine,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector<const FieldBase *> & fields)
{
  if ( fields.empty() ) { return; }

  const unsigned parallel_size = parallel_machine_size( machine );
  const unsigned parallel_rank = parallel_machine_rank( machine );
  const bool     asymmetric    = & domain != & range ;

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  std::vector<EntityProc>::const_iterator i ;

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        e_size += field_data_size( f , e );
      }
      send_size[ p ] += e_size ;
    }
  }

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        e_size += field_data_size( f , e );
      }
      recv_size[ p ] += e_size ;
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const s_size = & send_size[0] ;
    const unsigned * const r_size = & recv_size[0] ;
    sparse.allocate_buffers( machine, parallel_size / 4 , s_size, r_size);
  }

  // Pack for send:

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      CommBuffer & b = sparse.send_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );
        if ( size ) {
          unsigned char * ptr = reinterpret_cast<unsigned char *>(field_data( f , e ));
          b.pack<unsigned char>( ptr , size );
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      CommBuffer & b = sparse.recv_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );
        if ( size ) {
          unsigned char * ptr = reinterpret_cast<unsigned char *>(field_data( f , e ));
          b.unpack<unsigned char>( ptr , size );
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void communicate_field_data(
  const BulkData & mesh ,
  const unsigned field_count ,
  const FieldBase * fields[] ,
  CommAll & sparse )
{
  const std::vector<Entity*> & entity_comm = mesh.entity_comm();

  const unsigned parallel_size = mesh.parallel_size();

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> msg_size( parallel_size , zero );

  size_t j = 0;

  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( std::vector<Entity*>::const_iterator
          i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {
      Entity & e = **i ;
      const unsigned size = field_data_size( f , e );
      if ( size ) {
        for ( PairIterEntityComm
              ec = e.comm() ; ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {
          msg_size[ ec->proc ] += size ;
        }
      }
    }
  }

  // Allocate send and receive buffers:

  {
    const unsigned * const s_size = & msg_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), parallel_size / 4 , s_size, s_size);
  }

  // Pack for send:

  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( std::vector<Entity*>::const_iterator
          i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {
      Entity & e = **i ;
      const unsigned size = field_data_size( f , e );
      if ( size ) {
        unsigned char * ptr =
          reinterpret_cast<unsigned char *>(field_data( f , e ));
        for ( PairIterEntityComm
              ec = e.comm() ; ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {
          CommBuffer & b = sparse.send_buffer( ec->proc );
          b.pack<unsigned char>( ptr , size );
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();
}

void communicate_field_data_verify_read( CommAll & sparse )
{
  std::ostringstream msg ;
  int error = 0 ;
  for ( unsigned p = 0 ; p < sparse.parallel_size() ; ++p ) {
    if ( sparse.recv_buffer( p ).remaining() ) {
      msg << "P" << sparse.parallel_rank()
          << " Unread data from P" << p << std::endl ;
      error = 1 ;
    }
  }
  all_reduce( sparse.parallel() , ReduceSum<1>( & error ) );
  ThrowErrorMsgIf( error, msg.str() );
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk_classic

