
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
 * distribute copies to the public, perform publicly and display==================================================
// DoReverse method
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_MpiDistributor.h"
#include "Epetra_MpiComm.h"


//==============================================================================
// Epetra_MpiDistributor constructor
Epetra_MpiDistributor::Epetra_MpiDistributor (const Epetra_MpiComm &Comm): 
  Epetra_Object("Epetra::MpiDistributor"),
  plan_(0),
  comm_(Comm.GetMpiComm()),  
  tag_(Comm.GetMpiTag()),
  epComm_(&Comm)
     {
     }

//==============================================================================
// Epetra_MpiDistributor constructor
Epetra_MpiDistributor::Epetra_MpiDistributor (
  const Epetra_MpiDistributor &Distributor):
  Epetra_Object("Epetra::MpiDistributor"),
  plan_(0),
  comm_(Distributor.comm_),  
  tag_(Distributor.tag_),
  epComm_(Distributor.epComm_)
     {
     Zoltan_Comm_Create_Copy (&plan_, Distributor.plan_);
     }

//==============================================================================
// Epetra_MpiDistributor destructor
Epetra_MpiDistributor::~Epetra_MpiDistributor()
   {
   Zoltan_Comm_Destroy (&plan_);
   }

//==============================================================================
// CreateFromSends Method
// - create communication plan given a known list of procs to send to
int Epetra_MpiDistributor::CreateFromSends (
 const int  &NumExportIDs,
 const int  *ExportPIDs,
 const bool &Deterministic,
 int        &NumRemoteIDs)
    {
    Zoltan_Comm_Create (&plan_, (int) NumExportIDs, (int*) ExportPIDs, comm_,
     tag_, &NumRemoteIDs);            
    return 0;
    }

//==============================================================================
// CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
int Epetra_MpiDistributor::CreateFromRecvs (
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
    
    Epetra_MpiDistributor tmpdist (*epComm_);
    
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
       
    Zoltan_Comm_Create (&tmpdist.plan_, (int) NumRemoteIDs, proclist,
     tmpdist.comm_, tmpdist.tag_, &NumExportIDs);
     
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
    
    Zoltan_Comm_Do (tmpdist.plan_, tmpdist.tag_, (char*) imports, 2*sizeof(int),
     (char*) exports);
     
    for (i = 0; i < NumExportIDs; i++)  {
       ExportGIDs[i] = exports[2*i];
       ExportPIDs[i] = exports[2*i+1];
       }

    Zoltan_Comm_Create (&plan_, (int) NumExportIDs, (int*) ExportPIDs, comm_,
     tag_, (int*) NumRemoteIDs);
     
    if (proclist)  delete [] proclist;
    if (imports)   delete [] imports;
    if (exports)   delete [] exports;       
    return 0; 
    }

//==============================================================================
// Do method
int Epetra_MpiDistributor::Do (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Post (plan_, tag_, exports, (int) size, imports);
    Zoltan_Comm_Do_Wait (plan_, tag_, exports, (int) size, imports);
    return 0;
    }

//==============================================================================
// DoReverse method
int Epetra_MpiDistributor::DoReverse (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, (int) size, 0, imports);
    Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports, (int) size, 0, imports);
    return 0;
    }
    
//==============================================================================
// Do_Posts Method
int Epetra_MpiDistributor::DoPosts (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Post (plan_, tag_, exports, (int) size, imports);   
    return 0;               
    }
    
//==============================================================================
// Do_Waits Method
int Epetra_MpiDistributor::DoWaits (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Wait (plan_, tag_, exports, (int) size, imports);
    return 0;
    }

//==============================================================================--------------------------------------------
// DoReverse_Posts Method
int Epetra_MpiDistributor::DoReversePosts (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, (int) size, 0, imports);
    return 0;
    }

//==============================================================================
// DoReverse_Waits Method
int Epetra_MpiDistributor::DoReverseWaits (
 char      *exports,
 const int &size,
 char      *imports)
    {
    Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports, (int) size, 0, imports);
    return 0; 
    }

//==============================================================================
// Resize Method 
int Epetra_MpiDistributor::Resize (
 int *sizes,
 int *sum_recv_sizes)
    {
    Zoltan_Comm_Resize (plan_, sizes, tag_, sum_recv_sizes);
    return 0;
    }

//==============================================================================
// Do method with variable size objects
int Epetra_MpiDistributor::Do (
 char       *exports,
 const int *&size,
 char       *imports)
    {
    int junk;
    Zoltan_Comm_Resize  (plan_, (int*) size, tag_, &junk);
    Zoltan_Comm_Do_Post (plan_, tag_, exports, 1, imports);
    Zoltan_Comm_Do_Wait (plan_, tag_, exports, 1, imports);
    return 0;      
    }

//==============================================================================
// DoReverse method with variable size objects
int Epetra_MpiDistributor::DoReverse (
char       *exports,
const int *&size,
char       *imports)
   {
   Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, 1, (int*) size, imports);
   Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports, 1, (int*) size, imports);
   return 0;      
   }
   
//==============================================================================
// Do_Posts Method with variable size objects
int Epetra_MpiDistributor::DoPosts (
 char       *exports,
 const int *&size,
 char       *imports)
    {
    int junk;
    Zoltan_Comm_Resize  (plan_, (int*) size, tag_, &junk);
    Zoltan_Comm_Do_Post (plan_, tag_, exports, 1, imports);
    return 0;      
    }

//==============================================================================
// Do_Waits Method with variable size objects
int Epetra_MpiDistributor::DoWaits (
 char       *exports,
 const int *&size,
 char       *imports)
    {
    Zoltan_Comm_Do_Wait (plan_, tag_, exports, 1, imports);
    return 0;      
    }

//==============================================================================
// DoReverse_Posts Method with variable size objects
int Epetra_MpiDistributor::DoReversePosts (
 char       *exports,
 const int *&size,
 char       *imports)
    {
    Zoltan_Comm_Do_Reverse_Post (plan_, tag_, exports, 1, (int*)size, imports);
    return 0;      
    }

//==============================================================================
// DoReverse_Waits Method with variable size objects
int Epetra_MpiDistributor::DoReverseWaits (
 char       *exports,
 const int *&size,
 char       *imports)
    {
    Zoltan_Comm_Do_Reverse_Wait (plan_, tag_, exports, 1, (int*)size, imports);
    return 0;      
    }

//==============================================================================
// Print method
void Epetra_MpiDistributor::Print( ostream & os) const
{
  int i, j;
  int nsends, *send_procs, *send_lengths, send_nvals, send_max_size, *send_list;
  int nrecvs, *recv_procs, *recv_lengths, recv_nvals, recv_total_size;
  int *recv_list, self_msg;
  
  Zoltan_Comm_Info (plan_, &nsends, 0, 0, &send_nvals, &send_max_size, 0, &nrecvs,
   0, 0, &recv_nvals, &recv_total_size, 0, &self_msg);

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

