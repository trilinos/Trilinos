
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

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
