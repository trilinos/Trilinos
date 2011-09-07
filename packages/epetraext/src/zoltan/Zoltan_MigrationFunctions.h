//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationFunctions.h$
//
// Purpose        : Static methods which are directly registered with
//                  Zoltan.  They us the static container to access
//                  the dynamic object methods.
//
// Special Notes  : 
//
// Creator        : Robert J. Hoekstra, Parallel Computational Sciences
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------

#ifndef ZOLTAN_MIGRATIONFUNCTIONS_H_
#define ZOLTAN_MIGRATIONFUNCTIONS_H_

#include <EpetraExt_ConfigDefs.h>

#include <zoltan.h>

namespace Zoltan {

class EPETRAEXT_DEPRECATED MigrationFunctions
{

public:

  static int Object_Size        ( void * data,
				int num_gid_entries,
				int num_lid_entries,
                                ZOLTAN_ID_PTR global_id,
                                ZOLTAN_ID_PTR local_id,
				int * ierr );

  static void Pre_Migrate       ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				ZOLTAN_ID_PTR export_global_ids,
				ZOLTAN_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  static void Mid_Migrate       ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				ZOLTAN_ID_PTR export_global_ids,
				ZOLTAN_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  static void Post_Migrate      ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				ZOLTAN_ID_PTR export_global_ids,
				ZOLTAN_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  static void Pack_Object       ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int destination_processor,
				int size,
				char * buffer,
				int * ierr );

  static void Unpack_Object     ( void * data,
				int num_gid_entries,
				ZOLTAN_ID_PTR global_id,
				int size,
				char * buffer,
				int * ierr );

};

} //namespace Zoltan

#endif
