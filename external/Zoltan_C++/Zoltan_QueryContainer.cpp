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


int Zoltan_QueryContainer::CurrentObject = 0;

std::map< int, Zoltan_QueryObject * > Zoltan_QueryContainer::StaticMap;

void Zoltan_QueryContainer::setQueryID( const int & id )
{
  CurrentObject = id;
}

const int & Zoltan_QueryContainer::getQueryID()
{
  return CurrentObject;
}

bool Zoltan_QueryContainer::registerQueryObject( const int & id, 
	Zoltan_QueryObject * obj_ptr )
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

Zoltan_QueryObject * Zoltan_QueryContainer::getQueryObject(
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

