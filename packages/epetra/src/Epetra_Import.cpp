
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Import.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"


//==============================================================================
// Epetra_Import constructor for a Epetra_BlockMap object
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  TargetMap, const Epetra_BlockMap & SourceMap)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(TargetMap),
    SourceMap_(SourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{

  int i;
  
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = SourceMap.NumMyElements();
  int NumTargetIDs = TargetMap.NumMyElements();
  
  int *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int[NumTargetIDs];
    TargetMap.MyGlobalElements(TargetGIDs);
  }
  
  int * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int[NumSourceIDs];
    SourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);
  
  
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  
  // Find count of Target IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (SourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
  
  
  
  // Define remote and permutation lists
  
  int * RemoteGIDs=0;
  RemoteLIDs_ = 0;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (SourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = SourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=TargetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }

  if( NumRemoteIDs_>0 && !SourceMap.DistributedGlobal() )
    ReportError("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)", 1);
  
  // Test for distributed cases
  
  int * RemotePIDs = 0;

  if (SourceMap.DistributedGlobal()) {
    
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];
    int ierr = SourceMap.RemoteIDList(NumRemoteIDs_, RemoteGIDs, RemotePIDs, 0); // Get remote PIDs
    if (ierr) throw ReportError("Error in SourceMap.RemoteIDList call", ierr);

    //Get rid of IDs that don't exist in SourceMap
    if(NumRemoteIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumRemoteIDs_; ++i )
        if( RemotePIDs[i] == -1 ) ++cnt;
      if( cnt ) {
        if( NumRemoteIDs_-cnt ) {
          int * NewRemoteGIDs = new int[NumRemoteIDs_-cnt];
          int * NewRemotePIDs = new int[NumRemoteIDs_-cnt];
          cnt = 0;
          for( i = 0; i < NumRemoteIDs_; ++i )
            if( RemotePIDs[i] != -1 ) {
              NewRemoteGIDs[cnt] = RemoteGIDs[i];
              NewRemotePIDs[cnt] = RemotePIDs[i];
              ++cnt;
            }
          NumRemoteIDs_ = cnt;
          delete [] RemoteGIDs;
          delete [] RemotePIDs;
          RemoteGIDs = NewRemoteGIDs;
          RemotePIDs = NewRemotePIDs;
          ReportError("Warning in Epetra_Import: Target IDs not found in Source Map (Do you want to import to subset of Target Map?)", 1);
        }
        else { //valid RemoteIDs empty
          NumRemoteIDs_ = 0;
          delete [] RemoteGIDs;
          RemoteGIDs = 0;
          delete [] RemotePIDs;
          RemotePIDs = 0;
        }
      }
    }

    //Sort Remote IDs by processor so DoReverses will work
    Epetra_Util util;
    int * tmpPtr[2];
    tmpPtr[0] = RemoteLIDs_, tmpPtr[1] = RemoteGIDs;
    util.Sort(true,NumRemoteIDs_,RemotePIDs,0,0,2,tmpPtr);

    Distor_ = SourceMap.Comm().CreateDistributor();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    bool Deterministic = true;
    ierr = Distor_->CreateFromRecvs( NumRemoteIDs_, RemoteGIDs, RemotePIDs,
                       Deterministic, NumExportIDs_, ExportLIDs_, ExportPIDs_ );
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);

    // Export IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumExportIDs_; i++) {
      if (ExportPIDs_[i] < 0) throw ReportError("TargetMap requested a GID that is not in the SourceMap.", -1);
      ExportLIDs_[i] = SourceMap.LID(ExportLIDs_[i]);
    }
  }

  if( NumRemoteIDs_>0 ) delete [] RemoteGIDs;
  if( NumRemoteIDs_>0 ) delete [] RemotePIDs;

  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Epetra_Import copy constructor 
Epetra_Import::Epetra_Import(const Epetra_Import & Importer)
  : Epetra_Object(Importer),
    TargetMap_(Importer.TargetMap_),
    SourceMap_(Importer.SourceMap_),
    NumSameIDs_(Importer.NumSameIDs_),
    NumPermuteIDs_(Importer.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Importer.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Importer.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Importer.NumSend_),
    NumRecv_(Importer.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Importer.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Importer.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Importer.RemoteLIDs_[i];
  }

  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Importer.ExportLIDs_[i];
      ExportPIDs_[i] = Importer.ExportPIDs_[i];
    }
  }

  if (Importer.Distor_!=0) Distor_ = Importer.Distor_->Clone();

}

//==============================================================================
// Epetra_Import destructor 
Epetra_Import::~Epetra_Import()
{
  if( Distor_ != 0 ) delete Distor_;
  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by Distor_
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}
//=============================================================================
void Epetra_Import::Print(ostream & os) const
{

  os << endl << endl << "Source Map:" << endl << endl;
  SourceMap_.Print(os);
  
  os << endl << endl << "Target Map:" << endl << endl;
  TargetMap_.Print(os);
  
  os << endl << endl << "Distributor:" << endl << endl;
  if (Distor_==0) os << "  Is empty...." << endl;
  else Distor_->Print(os);
  
  os << "Number of Same IDs = " << NumSameIDs_ << endl;

  os << "Number of Permute IDs = " << NumPermuteIDs_ << endl;

  os << "Number of Export IDs = " << NumExportIDs_ << endl;

  os << "Number of Remote IDs = " << NumRemoteIDs_ << endl;
  
  os << "Epetra_Import Print Needs attention!!!!" << endl;
  return;
}

