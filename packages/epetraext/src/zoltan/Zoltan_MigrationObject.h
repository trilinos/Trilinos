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
