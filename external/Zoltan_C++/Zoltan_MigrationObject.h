//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationObject.h$
//
// Purpose        : Base Class for dynamic versions of migration
//		    functions to be used by Zoltan.  The application
//		    requires a class derived from this base to
//		    be instantiated and registered with the
//		    Zoltan_LoadBalance object.
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
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------

#ifndef ZOLTAN_MIGRATIONOBJECT_H_
#define ZOLTAN_MIGRATIONOBJECT_H_

#include <lbi_const.h>

class Zoltan_MigrationObject
{

public:

  virtual int Object_Size     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				LB_ID_PTR global_id,
				LB_ID_PTR local_id,
				int * ierr );

  virtual void Pre_Migrate    (	void * data,	
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				LB_ID_PTR import_global_ids,
				LB_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				LB_ID_PTR export_global_ids,
				LB_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  virtual void Mid_Migrate   (	void * data,	
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				LB_ID_PTR import_global_ids,
				LB_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				LB_ID_PTR export_global_ids,
				LB_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  virtual void Post_Migrate   (	void * data,	
				int num_gid_entries,
				int num_lid_entries,
				int num_import,
				LB_ID_PTR import_global_ids,
				LB_ID_PTR import_local_ids,
				int * import_procs,
				int num_export,
				LB_ID_PTR export_global_ids,
				LB_ID_PTR export_local_ids,
				int * export_procs,
				int * ierr );

  virtual void Pack_Object    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				LB_ID_PTR global_id,
				LB_ID_PTR local_id,
				int destination_processor,
				int size,
				char * buffer,
				int * ierr );

  virtual void Unpack_Object  (	void * data,
				int num_gid_entries,
				LB_ID_PTR global_id,
				int size,
				char * buffer,
				int * ierr );

};

#endif
