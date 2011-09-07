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

#include <EpetraExt_ConfigDefs.h>

#include <map>

namespace Zoltan {

class MigrationObject;

class EPETRAEXT_DEPRECATED MigrationContainer
{

public:

  static void setMigrationID( const int & id );

  static const int & getMigrationID();

  static bool registerMigrationObject( const int & id, 
		  MigrationObject * obj_ptr );

  static MigrationObject * getMigrationObject( const int & id );

private:

  static int CurrentObject;

  static std::map< int, MigrationObject * > StaticMap;

};

} //namespace Zoltan

#endif
