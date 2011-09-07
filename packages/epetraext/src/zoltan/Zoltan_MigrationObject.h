// @HEADER
// ************************************************************************
// 
//            Zoltan_CPP: An Object-Oriented Interface To Zoltan
//                    Copyright (2001) Sandia Corporation
// 
// Questions? Contact Robert J. Hoekstra (rjhoeks@sandia.gov)
// 
// ************************************************************************
// @HEADER

#ifndef ZOLTAN_MIGRATIONOBJECT_H_
#define ZOLTAN_MIGRATIONOBJECT_H_

#include "EpetraExt_ConfigDefs.h"

#include <zoltan.h>

namespace Zoltan {

//! Zoltan::MigrationObject: A base class from which the user can derive an application specific support class for Zoltan's migration callback functions.
/*! As with Zoltan, the user only need implement those methods used by Zoltan during their
    application executions.  If Zoltan calls an unimplemented method, a fatal error will
    be generated.
*/

class EPETRAEXT_DEPRECATED MigrationObject
{

public:

  //@{ \name Setup Methods

  //! Supports ZOLTAN_OBJ_SIZE_FN_TYPE
  virtual int Object_Size     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  //@}

  //@{ \name Migrate Methods

  //! Supports ZOLTAN_PRE_MIGRATE_FN_TYPE
  virtual void Pre_Migrate    (	void * data,	
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

  //! Supports ZOLTAN_MID_MIGRATE_FN_TYPE
  virtual void Mid_Migrate   (	void * data,	
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

  //! Supports ZOLTAN_POST_MIGRATE_FN_TYPE
  virtual void Post_Migrate   (	void * data,	
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

  //@}

  //@{ \name Pack/Unpack Methods
 
  //! Supports ZOLTAN_PACK_OBJ_FN_TYPE
  virtual void Pack_Object    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int destination_processor,
				int size,
				char * buffer,
				int * ierr );

  //! Supports ZOLTAN_UNPACK_OBJ_FN_TYPE
  virtual void Unpack_Object  (	void * data,
				int num_gid_entries,
				ZOLTAN_ID_PTR global_id,
				int size,
				char * buffer,
				int * ierr );

  //@}

};

} //namespace Zoltan

#endif
