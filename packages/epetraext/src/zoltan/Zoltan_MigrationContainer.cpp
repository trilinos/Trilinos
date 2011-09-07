//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationContainer.h$
//
// Purpose        : Static container object allowing C style calls
//		    from Zoltan to access Dynamic MigrationObject's.
//
// Special Notes  : 
//
// Creator        : Robert J. Hoekstra, Parallel Computational Science
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#include <Zoltan_MigrationContainer.h>

#include <Zoltan_MigrationObject.h>


int Zoltan::MigrationContainer::CurrentObject = 0;

std::map< int, Zoltan::MigrationObject * > Zoltan::MigrationContainer::StaticMap;

EPETRAEXT_DEPRECATED
void Zoltan::MigrationContainer::setMigrationID( const int & id )
{
  CurrentObject = id;
}

EPETRAEXT_DEPRECATED
const int & Zoltan::MigrationContainer::getMigrationID()
{
  return CurrentObject;
}

EPETRAEXT_DEPRECATED
bool Zoltan::MigrationContainer::registerMigrationObject( const int & id, 
		Zoltan::MigrationObject * obj_ptr )
{
  if( StaticMap.find( id ) == StaticMap.end() )
  {
    StaticMap[ id ] = obj_ptr;
    return true;
  }
  else
  {
    // Redundant id
    return false;
  }
}

EPETRAEXT_DEPRECATED
Zoltan::MigrationObject * Zoltan::MigrationContainer::getMigrationObject(
	const int & id )
{
  if( StaticMap.find( id ) != StaticMap.end() )
  {
    return StaticMap[ id ];
  }
  else
  {
    // Not found
    return 0;
  }
}

