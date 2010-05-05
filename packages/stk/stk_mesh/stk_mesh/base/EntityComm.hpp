/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityComm_hpp
#define stk_mesh_EntityComm_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \brief  Is shared with any other process */
bool in_shared( const Entity & entity );

/** \brief  Is shared with a given process */
bool in_shared( const Entity & entity , unsigned proc );

/** \brief  Is a receive ghost copy of an entity */
bool in_receive_ghost( const Entity & entity );

/** \brief  Is a receive ghost copy of an entity */
bool in_receive_ghost( const Ghosting & ghost , const Entity & entity );

/** \brief  Is sent to a ghost copy of an entity on any process */
bool in_send_ghost( const Entity & entity );

/** \brief  Is sent to a ghost copy of an entity on a given process */
bool in_send_ghost( const Entity & entity , unsigned proc );

/** \brief  Is in ghosting either send to 'p' or receive from 'p' */
bool in_ghost( const Ghosting & ghost , const Entity & entity , unsigned p );

/** \brief  Is in owned closure of the given process,
 *          typically the local process.
 */
bool in_owned_closure( const Entity & entity , unsigned proc );

/** \brief  List of all entity communication processes, sorted */
void comm_procs( const Entity & entity , std::vector<unsigned> & procs );

/** \brief  List of entity communication processes for a given ghost, sorted */
void comm_procs( const Ghosting & ghost ,
                 const Entity & entity , std::vector<unsigned> & procs );


//----------------------------------------------------------------------

void pack_entity_info( CommBuffer & buf , const Entity & entity );
 
void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  unsigned       & owner , 
  PartVector     & parts , 
  std::vector<Relation> & relations );
 

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values( CommBuffer & , Entity & );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values( CommBuffer & , Entity & , std::ostream & error_msg );

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

