
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>


namespace stk {
namespace mesh {

void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = ghosts.mesh();
  const unsigned parallel_size = mesh.parallel_size();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( std::vector<EntityProc>::const_iterator
        i = ghosts.send().begin() ; i != ghosts.send().end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;
      e_size += field_data_size( f , e );
    } 
    send_size[ p ] += e_size ;
  }

  for ( std::vector<Entity*>::const_iterator
        i = ghosts.receive().begin() ; i != ghosts.receive().end() ; ++i ) {
    Entity       & e = **i ;
    const unsigned p = e.owner_rank();

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;
      e_size += field_data_size( f , e );
    } 
    recv_size[ p ] += e_size ;
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const s_size = & send_size[0] ;
    const unsigned * const r_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), parallel_size / 4 , s_size, r_size);
  }

  // Pack for send:

  for ( std::vector<EntityProc>::const_iterator
        i = ghosts.send().begin() ; i != ghosts.send().end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

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

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( std::vector<Entity*>::const_iterator
        i = ghosts.receive().begin() ; i != ghosts.receive().end() ; ++i ) {
    Entity       & e = **i ;
    const unsigned p = e.owner_rank();

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
  ParallelMachine machine ,
  const std::vector<EntityProc> & shared ,
  const unsigned field_count ,
  const FieldBase * fields[] ,
  CommAll & sparse )
{
  const unsigned parallel_size = parallel_machine_size( machine );

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> msg_size( parallel_size , zero );

  size_t j = 0;
  std::vector<EntityProc>::const_iterator i ;

  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( i = shared.begin() ; i != shared.end() ; ++i ) {
      Entity       & e = * i->first ;
      const unsigned p = i->second ;
      msg_size[ p ] += field_data_size( f , e );
    }
  }

  // Allocate send and receive buffers:

  {
    const unsigned * const s_size = & msg_size[0] ;
    sparse.allocate_buffers( machine, parallel_size / 4 , s_size, s_size);
  }

  // Pack for send:
 
  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( i = shared.begin() ; i != shared.end() ; ++i ) {
      Entity       & e = * i->first ;
      const unsigned p = i->second ;

      CommBuffer & b = sparse.send_buffer( p );

      const unsigned size = field_data_size( f , e );

      if ( size ) {
        unsigned char * ptr = reinterpret_cast<unsigned char *>(field_data( f , e ));
        b.pack<unsigned char>( ptr , size );
      }
    }
  }

  // Communicate:

  sparse.communicate();
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

