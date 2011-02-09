/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityComm.hpp>

namespace stk {
namespace mesh {

void Ghosting::send_list( std::vector< EntityProc > & v ) const
{
  for ( std::vector<Entity*>::const_iterator
        i =  m_mesh.entity_comm().begin() ;
        i != m_mesh.entity_comm().end() ; ++i ){
    Entity * const entity = *i ;
    if ( entity->owner_rank() == m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          v.push_back( EntityProc( entity , ec->proc ) );
        }
      }
    }
  }
}

void Ghosting::receive_list( std::vector< Entity * > & v ) const
{
  for ( std::vector<Entity*>::const_iterator
        i =  m_mesh.entity_comm().begin() ;
        i != m_mesh.entity_comm().end() ; ++i ){
    Entity * const entity = *i ;
    if ( entity->owner_rank() != m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          v.push_back( entity );
        }
      }
    }
  }
}

std::ostream& Ghosting::operator<<(std::ostream& out) const
{
  out << "Ghosting object: name: " << name()
      << ", ordinal: " << ordinal() << "\n";

  out << "  Locally owned entities ghosted on other processors (send list):\n";

  for ( std::vector<Entity*>::const_iterator
        i =  m_mesh.entity_comm().begin() ;
        i != m_mesh.entity_comm().end() ; ++i ){
    Entity * const entity = *i ;
    if ( entity->owner_rank() == m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          out << "    ";
          print_entity_key( out, MetaData::get(m_mesh), entity->key() );
          out << ", sending ghost to " << ec->proc << ", status is: "
              << entity->log_query() << "\n";
        }
      }
    }
  }

  out << "  Entities ghosted on this processor from the owner (recv list):\n";
  for ( std::vector<Entity*>::const_iterator
        i =  m_mesh.entity_comm().begin() ;
        i != m_mesh.entity_comm().end() ; ++i ){
    Entity * const entity = *i ;
    if ( entity->owner_rank() != m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          out << "    ";
          print_entity_key( out, MetaData::get(m_mesh), entity->key() );
          out << ", owner of ghost is " << entity->owner_rank()
              << ", status is: " << entity->log_query() << "\n";
        }
      }
    }
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const Ghosting& rhs)
{
  return rhs.operator<<(out);
}

}
}


