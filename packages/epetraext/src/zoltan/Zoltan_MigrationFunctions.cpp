/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationFunctions.C$
//
// Purpose        : Static methods which are directly registered with
//		    Zoltan.  They us the static container to access
//		    the dynamic object methods.
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

#include <Zoltan_MigrationFunctions.h>
#include <Zoltan_MigrationContainer.h>
#include <Zoltan_MigrationObject.h>

EPETRAEXT_DEPRECATED
int Zoltan::MigrationFunctions::Object_Size    (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int * ierr )
{
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  return obj_ptr->Object_Size( data, num_gid_entries, num_lid_entries,
		global_id, local_id, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationFunctions::Pre_Migrate   (	void * data,
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
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  obj_ptr->Pre_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationFunctions::Mid_Migrate   (	void * data,
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
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  obj_ptr->Mid_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationFunctions::Post_Migrate  (	void * data,
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
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  obj_ptr->Post_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationFunctions::Pack_Object   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int destination_processor,
						int size,
						char * buffer,
						int * ierr )
{
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  obj_ptr->Pack_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, destination_processor, size, buffer, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::MigrationFunctions::Unpack_Object (	void * data,
						int num_gid_entries,
						ZOLTAN_ID_PTR global_id,
						int size,
						char * buffer,
						int * ierr )
{
  Zoltan::MigrationObject * obj_ptr =
	Zoltan::MigrationContainer::getMigrationObject(
	Zoltan::MigrationContainer::getMigrationID() );

  obj_ptr->Unpack_Object( data, num_gid_entries, global_id, size,
	buffer, ierr );
}
