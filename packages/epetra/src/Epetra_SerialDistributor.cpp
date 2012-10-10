
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialComm.h"


//==============================================================================
// Epetra_SerialDistributor constructor
Epetra_SerialDistributor::Epetra_SerialDistributor(const Epetra_SerialComm & Comm)
  : Epetra_Object("Epetra::SerialDistributor"),
    nrecvs_(0),
    nsends_(0)
{
  (void)Comm;
}

//==============================================================================
Epetra_SerialDistributor::Epetra_SerialDistributor(const Epetra_SerialDistributor & Plan)
  : Epetra_Object("Epetra::SerialDistributor"),
    nrecvs_(Plan.nrecvs_),
    nsends_(Plan.nsends_)
{
}

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
					       bool Deterministic,
					       int & NumRemoteIDs )
{
  (void)Deterministic;
  NumRemoteIDs = 0;

  //In a SerialDistributor, myproc == 0 by definition.
  int myproc = 0;

  //basically just do a sanity check.
  for(int i=0; i<NumExportIDs; ++i) {
    if (ExportPIDs[i] != myproc) {
      cerr << "Epetra_SerialDistributor::CreateFromSends: ExportPIDs["<<i
          <<"]=="<<ExportPIDs[i]<<", not allowed for serial case."<<endl;
      return(-1);
    }
    ++NumRemoteIDs;
  }

  nrecvs_ = NumRemoteIDs;

  return(0);
}

//==============================================================================
//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::CreateFromRecvs( const int & NumRemoteIDs,
				   const int * RemoteGIDs,
			           const int * RemotePIDs,
				   bool Deterministic,
			           int & NumExportIDs,
				   int *& ExportGIDs,
				   int *& ExportPIDs )
{
  (void)NumRemoteIDs;
  (void)RemoteGIDs;
  (void)RemotePIDs;
  (void)Deterministic;
  (void)NumExportIDs;
  (void)ExportGIDs;
  (void)ExportPIDs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_SerialDistributor::CreateFromRecvs( const int & NumRemoteIDs,
				   const long long * RemoteGIDs,
			           const int * RemotePIDs,
				   bool Deterministic,
			           int & NumExportIDs,
				   long long *& ExportGIDs,
				   int *& ExportPIDs )
{
  (void)NumRemoteIDs;
  (void)RemoteGIDs;
  (void)RemotePIDs;
  (void)Deterministic;
  (void)NumExportIDs;
  (void)ExportGIDs;
  (void)ExportPIDs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
#endif

//==============================================================================
// GSComm_Comm Do method
int Epetra_SerialDistributor::Do(char * export_objs,
                                 int obj_size,
                                 int & len_import_objs,
                                 char *& import_objs )
{
  len_import_objs = obj_size*nrecvs_;
  if (len_import_objs > 0) {
    import_objs = new char[len_import_objs];
  }

  for(int i=0; i<len_import_objs; ++i) import_objs[i] = export_objs[i];

  return(0);
}

//==============================================================================
// GSComm_Comm DoReverse method
int Epetra_SerialDistributor::DoReverse(char * export_objs,
                                        int obj_size,
                                        int & len_import_objs,
                                        char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoPosts(char * export_objs,
                                      int obj_size,
                                      int & len_import_objs,
                                      char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoWaits()
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReversePosts(char * export_objs,
                                             int obj_size,
                                             int & len_import_objs,
                                             char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReverseWaits()
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm Do method
int Epetra_SerialDistributor::Do(char * export_objs,
                                 int obj_size,
                                 int *& sizes,
                                 int & len_import_objs,
                                 char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)sizes;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm DoReverse method
int Epetra_SerialDistributor::DoReverse(char * export_objs,
                                        int obj_size,
                                        int *& sizes,
                                        int & len_import_objs,
                                        char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)sizes;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoPosts(char * export_objs,
                                      int obj_size,
                                      int *& sizes,
                                      int & len_import_objs,
                                      char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)sizes;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_SerialDistributor::DoReversePosts(char * export_objs,
                                             int obj_size,
                                             int *& sizes,
                                             int & len_import_objs,
                                             char *& import_objs )
{
  (void)export_objs;
  (void)obj_size;
  (void)sizes;
  (void)len_import_objs;
  (void)import_objs;
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
void Epetra_SerialDistributor::Print( ostream & os) const
{
  os << "Trivial Distributor" << endl;
  return;
}
