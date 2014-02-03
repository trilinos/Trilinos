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
#include <stk_mesh/base/FieldParallel.hpp>

namespace stk {
namespace mesh {

void copy_owned_to_shared( const BulkData& mesh, const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  unsigned shared_id = mesh.ghostings()[0]->ordinal();

  const unsigned parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;
    const bool owned = i->owner == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;
     
      if(is_matching_rank(f, e)) {
        e_size += field_bytes_per_entity( f , e );
      }
    }

    if (e_size == 0) {
      continue;
    }

    for ( PairIterEntityComm ec = mesh.entity_comm(i->key) ; ! ec.empty() ; ++ec ) {
      if ( shared_id == ec->ghost_id ) {
        if ( owned ) {
          send_size[ ec->proc ] += e_size ;
        }
        else {
          recv_size[ i->owner ] += e_size ;
          break;
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

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;
    if ( i->owner == parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

        if(!is_matching_rank(f, e)) continue;

        const unsigned size = field_bytes_per_entity( f , e );
	

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));

          for ( PairIterEntityComm ec = mesh.entity_comm(i->key); !ec.empty(); ++ec ) {

            if ( shared_id == ec->ghost_id ) {
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

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;
    if ( i->owner != parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

        if(!is_matching_rank(f, e)) continue;
        const unsigned size = field_bytes_per_entity( f , e );

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));

          for ( PairIterEntityComm ec = mesh.entity_comm(i->key) ; ! ec.empty() ; ++ec ) {

            if ( shared_id == ec->ghost_id ) {
              CommBuffer & b = sparse.recv_buffer( i->owner );
              b.unpack<unsigned char>( ptr , size );
              break;
            }
          }
        }
      }
    }
  }
}

void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = BulkData::get(ghosts);
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;

    const bool owned = i->owner == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      if(is_matching_rank(f, e)) {
        e_size += field_bytes_per_entity( f , e );
      }
    }

    if (e_size == 0) {
      continue;
    }

    for ( PairIterEntityComm ec = mesh.entity_comm(i->key) ; ! ec.empty() ; ++ec ) {
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

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;


    if ( i->owner == parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

        if(!is_matching_rank(f, e)) continue;

        const unsigned size = field_bytes_per_entity( f , e );

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));

          for ( PairIterEntityComm ec = mesh.entity_comm(i->key); !ec.empty(); ++ec ) {

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

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.comm_list().begin() ;
        i != mesh.comm_list().end() ; ++i ) {
    Entity e = i->entity;
    if ( i->owner != parallel_rank ) {

      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
 
        if(!is_matching_rank(f, e)) continue;

        const unsigned size = field_bytes_per_entity( f , e );

        if ( size ) {
          unsigned char * ptr =
            reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));

          for ( PairIterEntityComm ec = mesh.entity_comm(i->key) ; ! ec.empty() ; ++ec ) {

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
  const BulkData& mesh,
  ParallelMachine machine,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector<const FieldBase *> & fields)
{
  if ( fields.empty() ) { return; }

  const int parallel_size = parallel_machine_size( machine );
  const int parallel_rank = parallel_machine_rank( machine );
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
    Entity e = i->first;
    const int p = i->second ;

    if ( asymmetric || parallel_rank == mesh.parallel_owner_rank(e) ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

	if(!is_matching_rank(f, e)) continue;

        e_size += field_bytes_per_entity( f , e );
      }
      send_size[ p ] += e_size ;
    }
  }

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity e = i->first;
    const int p = i->second ;

    if ( asymmetric || p == mesh.parallel_owner_rank(e) ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

	if(!is_matching_rank(f, e)) continue;

        e_size += field_bytes_per_entity( f , e );
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
    Entity e = i->first;
    const int p = i->second ;

    if ( asymmetric || parallel_rank == mesh.parallel_owner_rank(e) ) {
      CommBuffer & b = sparse.send_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

	if(!is_matching_rank(f, e)) continue;

        const unsigned size = field_bytes_per_entity( f , e );
        if ( size ) {
          unsigned char * ptr = reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));
          b.pack<unsigned char>( ptr , size );
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity e = i->first;
    const int p = i->second ;

    if ( asymmetric || p == mesh.parallel_owner_rank(e) ) {
      CommBuffer & b = sparse.recv_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;

	if(!is_matching_rank(f, e)) continue;

        const unsigned size = field_bytes_per_entity( f , e );


        if ( size ) {
          unsigned char * ptr = reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));
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
  const FieldBase * const *fields ,
  CommAll & sparse )
{
  const EntityCommListInfoVector & entity_comm = mesh.comm_list();

  const int parallel_size = mesh.parallel_size();

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> msg_size( parallel_size , zero );

  size_t j = 0;

  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( EntityCommListInfoVector::const_iterator
          i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {
      Entity e = i->entity;

      if(!is_matching_rank(f, e)) continue;

      const unsigned size = field_bytes_per_entity( f , e );
      if ( size ) {
        PairIterEntityComm ec = mesh.entity_comm(i->key);
        for (; ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {
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
    for ( EntityCommListInfoVector::const_iterator
          i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {
      Entity e = i->entity;

      if(!is_matching_rank(f, e)) continue;

      const unsigned size = field_bytes_per_entity( f , e );
      if ( size ) {
        unsigned char * ptr =
          reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , e ));
        PairIterEntityComm ec = mesh.entity_comm(i->key);
        for (; ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {
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
  for ( int p = 0 ; p < sparse.parallel_size() ; ++p ) {
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


/** Sum (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */
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
//----------------------------------------------------------------------

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
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


//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

