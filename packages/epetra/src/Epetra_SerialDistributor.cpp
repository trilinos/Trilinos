
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Epetra_SerialDistributor.h"


//==============================================================================
// Epetra_SerialDistributor constructor
Epetra_SerialDistributor::Epetra_SerialDistributor(const Epetra_SerialComm & Comm): 
Epetra_Object("Epetra::SerialDistributor")
{}

//==============================================================================
Epetra_SerialDistributor::Epetra_SerialDistributor(const Epetra_SerialDistributor & Plan):
Epetra_Object("Epetra::SerialDistributor")
{}

//==============================================================================
// Epetra_SerialDistributor destructor
Epetra_SerialDistributor::~Epetra_SerialDistributor() {
}


//==============================================================================
//---------------------------------------------------------------------------
//CreateFromSends Method
// - create communication plan given a known list of procs to send to
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::CreateFromSends( const int & NumExportIDs,
			           const int * ExportPIDs,
			           const bool & Deterministic,
			           int & NumRemoteIDs ) {
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::CreateFromRecvs( const int & NumRemoteIDs,
				   const int * RemoteGIDs,
			           const int * RemotePIDs,
			           const bool & Deterministic,
			           int & NumExportIDs,
				   int *& ExportGIDs,
				   int *& ExportPIDs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm Do method
int Epetra_SerialDistributor::Do(char * export_objs, const int & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm DoReverse method
int Epetra_SerialDistributor::DoReverse(char * export_objs,const int & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoPosts(char * export_objs,
				const int & obj_size,
				char * import_objs ) {

  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoWaits(char * export_objs,
			       const int & obj_size,
			       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReversePosts(char * export_objs,
				       const int & obj_size,
				       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReverseWaits(char * export_objs,
		 	           const int & obj_size,
			           char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm Do method
int Epetra_SerialDistributor::Do(char * export_objs, const int * & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm DoReverse method
int Epetra_SerialDistributor::DoReverse(char * export_objs,const int * & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoPosts(char * export_objs,
				const int * & obj_size,
				char * import_objs ) {

  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoWaits(char * export_objs,
			       const int * & obj_size,
			       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReversePosts(char * export_objs,
				       const int * & obj_size,
				       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReverseWaits(char * export_objs,
		 	           const int * & obj_size,
			           char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
void Epetra_SerialDistributor::Print( ostream & os) const
{
  os << "Trivial Distributor" << endl;
  return;
}
