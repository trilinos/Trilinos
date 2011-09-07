//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationObject.C$
//
// Purpose        : Base Class for dynamic versions of migration
//                  functions to be used by Zoltan.  The application
//                  requires a class derived from this base to
//                  be instantiated and registered with the
//                  Zoltan_LoadBalance object.
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

#include <Zoltan_MigrationObject.h>

#include <iostream>

EPETRAEXT_DEPRECATED
int Zoltan::MigrationObject::Object_Size       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Object_Size( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int * )"
	<< " must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationObject::Pre_Migrate      (	void * data,
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
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Pre_Migrate( void *, "
	<< "int, int, int ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int, ZOLTAN_ID_PTR, "
	<< "ZOLTAN_ID_PTR, int *, int ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationObject::Mid_Migrate      (	void * data,
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
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Mid_Migrate( void *, "
	<< "int, int, int ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int, ZOLTAN_ID_PTR, "
	<< "ZOLTAN_ID_PTR, int *, int ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationObject::Post_Migrate     (	void * data,
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
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Post_Migrate( void *, "
	<< "int, int, int ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int, ZOLTAN_ID_PTR, "
	<< "ZOLTAN_ID_PTR, int *, int ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationObject::Pack_Object      (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int destination_processor,
						int size,
						char * buffer,
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Pack_Object( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, char *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationObject::Unpack_Object    ( void * data,
						int num_gid_entries,
						ZOLTAN_ID_PTR global_id,
						int size,
						char * buffer,
						int * ierr )
{
  std::cout << "Error: int Zoltan_MigrationObject::Unpack_Object( void *, "
	<< "int, ZOLTAN_ID_PTR, int, char *, int * )"
	<< " must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}
