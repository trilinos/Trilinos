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

#include "EpetraExt_ZoltanMpiDistributor.h"
#include "EpetraExt_ZoltanMpiComm.h"


//==============================================================================
// constructor
EPETRAEXT_DEPRECATED
EpetraExt::ZoltanMpiDistributor::ZoltanMpiDistributor(const EpetraExt::ZoltanMpiComm &Comm)
 : Epetra_Object("EpetraExt::ZoltanMpiDistributor"),
   plan_(0),
   comm_(Comm.GetMpiComm()),  
   tag_(Comm.GetMpiTag()),
   epComm_(&Comm)
{
}

//==============================================================================
// constructor
EPETRAEXT_DEPRECATED
EpetraExt::ZoltanMpiDistributor::ZoltanMpiDistributor( const EpetraExt::ZoltanMpiDistributor &Distributor)
 : Epetra_Object("EpetraExt::ZoltanMpiDistributor"),
   plan_(0),
   comm_(Distributor.comm_),  
   tag_(Distributor.tag_),
   epComm_(Distributor.epComm_)
{
  EPETRA_CHK_ERR(Zoltan_Comm_Create_Copy (&plan_, Distributor.plan_));
}

//==============================================================================
// destructor
EPETRAEXT_DEPRECATED
EpetraExt::ZoltanMpiDistributor::~ZoltanMpiDistributor()
{
  Zoltan_Comm_Destroy(&plan_);
}

//==============================================================================
// CreateFromSends Method
// - create communication plan given a known list of procs to send to
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::CreateFromSends (
 const int  &NumExportIDs,
 const int  *ExportPIDs,
 const bool &Deterministic,
 int        &NumRemoteIDs)
{
  EPETRA_CHK_ERR (Zoltan_Comm_Create (&plan_, (int) NumExportIDs,
     (int*) ExportPIDs, comm_, tag_, &NumRemoteIDs));            
  return 0;
}

//==============================================================================
// CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::CreateFromRecvs (
 const int  &NumRemoteIDs,
 const int  *RemoteGIDs,
 const int  *RemotePIDs,
 const bool &Deterministic,
 int  &NumExportIDs,
 int *&ExportGIDs,
 int *&ExportPIDs)
{
    int i, myproc, nprocs;
    MPI_Comm_rank (comm_, &myproc);
    MPI_Comm_size (comm_, &nprocs);
    
    ZoltanMpiDistributor tmpdist (*epComm_);
    
    int *proclist = 0, *imports = 0;
    if (NumRemoteIDs > 0)  {
       proclist = new int [NumRemoteIDs];
       imports  = new int [2 * NumRemoteIDs];
       for (i = 0; i < NumRemoteIDs; i++)  {
          proclist[i]     = RemotePIDs[i];
          imports [2*i]   = RemoteGIDs[i];
          imports [2*i+1] = myproc;
          }
       }
       
    EPETRA_CHK_ERR (Zoltan_Comm_Create (&tmpdist.plan_, (int) NumRemoteIDs,
     proclist, tmpdist.comm_, tmpdist.tag_, &NumExportIDs));
     
    int *exports = 0; 
    if (NumExportIDs > 0)  {
       exports    = new int [2 * NumExportIDs];
       ExportGIDs = new int [NumExportIDs];
       ExportPIDs = new int [NumExportIDs];
       }
    else  {
       ExportGIDs = 0;
       ExportPIDs = 0;
       }       
    
    EPETRA_CHK_ERR (Zoltan_Comm_Do (tmpdist.plan_, tmpdist.tag_,
     (char*) imports, 2*sizeof(int), (char*) exports));
     
    for (i = 0; i < NumExportIDs; i++)  {
       ExportGIDs[i] = exports[2*i];
       ExportPIDs[i] = exports[2*i+1];
       }

    EPETRA_CHK_ERR (Zoltan_Comm_Create (&plan_, (int) NumExportIDs,
     (int*) ExportPIDs, comm_, tag_, (int*) NumRemoteIDs));
     
    if (proclist)  delete [] proclist;
    if (imports)   delete [] imports;
    if (exports)   delete [] exports;       
    return 0; 
}

//==============================================================================
// Do method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::Do (
 char      *exports,
 int       size,
 int       &len_imports,
 char      *imports)
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Post (plan_,tag_,exports,size,imports));
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Wait (plan_,tag_,exports,size,imports));
    return 0;
}

//==============================================================================
// DoReverse method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoReverse (
 char      *exports,
 int       size,
 int       &len_imports,
 char      *imports)
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports,
     size, 0, imports));
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports,
     size, 0, imports));
    return 0;
}
    
//==============================================================================
// Do_Posts Method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoPosts (
 char      *exports,
 int       size,
 int       &len_imports,
 char      *imports)
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Post (plan_, tag_, exports, size, imports));   
    return 0;               
}
    
//==============================================================================
// Do_Waits Method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoWaits ()
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Wait (plan_, tag_, 0, 0, imports));
    return 0;
}

//==============================================================================--------------------------------------------
// DoReverse_Posts Method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoReversePosts (
 char      *exports,
 int       size,
 int       &len_imports,
 char      *imports)
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, size, 0, imports));
    return 0;
}

//==============================================================================
// DoReverse_Waits Method
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoReverseWaits ()
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, 0, 0, 0, 0));
    return 0; 
}

//==============================================================================
// Resize Method 
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::Resize (
 int *sizes)
{
    int * sum_recv_sizes = 0;
    EPETRA_CHK_ERR (Zoltan_Comm_Resize (plan_, sizes, tag_, sum_recv_sizes));
    return 0;
}

//==============================================================================
// Do method with variable size objects
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::Do (
 char       *exports,
 int        obj_size,
 int       *&sizes,
 int        &len_imports,
 char       *imports)
{
    int junk;
    EPETRA_CHK_ERR (Zoltan_Comm_Resize  (plan_, (int*) sizes, tag_, &junk));
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Post (plan_, tag_, exports, 1, imports));
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Wait (plan_, tag_, exports, 1, imports));
    return 0;      
}

//==============================================================================
// DoReverse method with variable size objects
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoReverse (
 char       *exports,
 int        obj_size,
 int       *&sizes,
 int        &len_imports,
 char       *imports)
{
   EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, 1, sizes, imports));
   EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports, 1, sizes, imports));
   return 0;      
}
   
//==============================================================================
// Do_Posts Method with variable size objects
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoPosts (
 char       *exports,
 int        obj_size,
 int       *&sizes,
 int        &len_imports,
 char       *imports)
{
    int junk;
    EPETRA_CHK_ERR (Zoltan_Comm_Resize  (plan_, (int*) sizes, tag_, &junk));
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Post (plan_, tag_, exports, 1, imports));
    return 0;      
}

//==============================================================================
// DoReverse_Posts Method with variable size objects
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanMpiDistributor::DoReversePosts (
 char       *exports,
 int        obj_size,
 int       *&sizes,
 int        &len_imports,
 char       *imports)
{
    EPETRA_CHK_ERR (Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, 1, sizes, imports));
    return 0;      
}

//==============================================================================
// Print method
EPETRAEXT_DEPRECATED
void EpetraExt::ZoltanMpiDistributor::Print( ostream & os) const
{
  int i, j;
  int nsends, *send_procs, *send_lengths, send_nvals, send_max_size, *send_list;
  int nrecvs, *recv_procs, *recv_lengths, recv_nvals, recv_total_size;
  int *recv_list, self_msg;
  
  Zoltan_Comm_Info (plan_, &nsends, 0, 0, &send_nvals, &send_max_size, 0,
   &nrecvs, 0, 0, &recv_nvals, &recv_total_size, 0, &self_msg);

  send_procs   = new int [nsends];
  send_lengths = new int [nsends];
  send_list    = new int [send_nvals];
  recv_procs   = new int [nrecvs];
  recv_lengths = new int [nrecvs];
  recv_list    = new int [recv_nvals];

  Zoltan_Comm_Info (plan_, 0, send_procs, send_lengths, 0, 0, send_list, 0,
  recv_procs, recv_lengths, 0, 0, recv_list, 0);

  os << "nsends: " << nsends << endl;
  os << "procs_to: ";
  for( i = 0; i < nsends; i++ )
    os << " " << send_procs[i];
  os << endl;
  os<< "lengths_to: ";
  for( i = 0; i < nsends; i++ )
    os << " " << send_lengths[i];
  os << endl;
  os << "indices_to: ";
//  int k = 0;
//  for( i = 0; i < nsends; i++ )
//  {
//    for( j = 0; j < send_lengths[i]; j++ )
//      os << " " << plan_->indices_to[j+k];
//    k += send_lengths[i];
//  }
  os << endl;
  os << "nrecvs: " << nrecvs << endl;
  os << "procs_from: ";
  for( i = 0; i < nrecvs; i++ )
    os << " " << recv_procs[i];
  os << endl;
  os << "lengths_from: ";
  for( i = 0; i < nrecvs; i++ )
    os << " " << recv_lengths[i];
  os << endl;
/*
  os << "indices_from: ";
  k = 0;
  for( i = 0; i < nrecvs; i++ )
  {
    for( j = 0; j < recv_lengths[i]; j++ )
      os << " " << plan_->indices_from[j+k];
    k += recv_lengths[i];
  }
*/
  os << "self_msg: " << self_msg << endl;
  os << "max_send_length: " << send_max_size << endl;
  os << "total_recv_length: " << recv_total_size << endl;
  os << endl;

  delete [] send_procs;
  delete [] send_lengths;
  delete [] send_list;
  delete [] recv_procs;
  delete [] recv_lengths;
  delete [] recv_list;
  
  return;
}

