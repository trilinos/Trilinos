//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryContainer.C$
//
// Purpose        : Static Container object to allow Static (C-style)
//                  functions to access methods of dynamic objects.
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
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

#include <Zoltan_QueryContainer.h>

#include <Zoltan_QueryObject.h>


int Zoltan::QueryContainer::CurrentObject = 0;

std::map< int, Zoltan::QueryObject * > Zoltan::QueryContainer::StaticMap;

EPETRAEXT_DEPRECATED
void Zoltan::QueryContainer::setQueryID( const int & id )
{
  CurrentObject = id;
}

EPETRAEXT_DEPRECATED
const int & Zoltan::QueryContainer::getQueryID()
{
  return CurrentObject;
}

EPETRAEXT_DEPRECATED
bool Zoltan::QueryContainer::registerQueryObject( const int & id, 
		Zoltan::QueryObject * obj_ptr )
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
Zoltan::QueryObject * Zoltan::QueryContainer::getQueryObject(
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

