//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationContainer.h$
//
// Purpose        : Static Container object to allow Static (C-style) 
//		    functions to access methods of dynamic objects.
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

#ifndef ZOLTAN_MIGRATIONCONTAINER_H_
#define ZOLTAN_MIGRATIONCONTAINER_H_

#include <map>

class Zoltan_MigrationObject;

class Zoltan_MigrationContainer
{

public:

  static void setMigrationID( const int & id );

  static const int & getMigrationID();

  static bool registerMigrationObject( const int & id, 
	Zoltan_MigrationObject * obj_ptr );

  static Zoltan_MigrationObject * getMigrationObject( const int & id );

private:

  static int CurrentObject;

  static std::map< int, Zoltan_MigrationObject * > StaticMap;

};

#endif
