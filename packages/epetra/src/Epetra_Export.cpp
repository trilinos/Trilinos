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

#include "Epetra_ConfigDefs.h"
#include "Epetra_Export.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"
#include "Epetra_Import.h"
#include <vector>

#ifdef HAVE_MPI
#include "Epetra_MpiDistributor.h"
#endif


//==============================================================================
// Epetra_Export constructor function for a Epetra_BlockMap object
template<typename int_type>
void Epetra_Export::Construct( const Epetra_BlockMap &  sourceMap, const Epetra_BlockMap & targetMap)
{

  int i;

  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumExportIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be exported.

  int NumSourceIDs = sourceMap.NumMyElements();
  int NumTargetIDs = targetMap.NumMyElements();

  int_type *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int_type[NumTargetIDs];
    targetMap.MyGlobalElements(TargetGIDs);
  }

  int_type * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int_type[NumSourceIDs];
    sourceMap.MyGlobalElements(SourceGIDs);
  }

  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);

  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;

  // Find count of Source IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumExportIDs_ = 0;
  for (i=NumSameIDs_; i< NumSourceIDs; i++)
    if (targetMap.MyGID(SourceGIDs[i])) NumPermuteIDs_++; // Check if Source GID is a local Target GID
    else NumExportIDs_++; // If not, then it is remote

  // Define remote and permutation lists

  int_type * ExportGIDs = 0;
  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportGIDs = new int_type[NumExportIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }

  NumPermuteIDs_ = 0;
  NumExportIDs_ = 0;
  for (i=NumSameIDs_; i< NumSourceIDs; i++) {
    if (targetMap.MyGID(SourceGIDs[i])) {
      PermuteFromLIDs_[NumPermuteIDs_] = i;
      PermuteToLIDs_[NumPermuteIDs_++] = targetMap.LID(SourceGIDs[i]);
    }
    else {
      //NumSend_ +=sourceMap.ElementSize(i); // Count total number of entries to send
      NumSend_ +=sourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
      ExportGIDs[NumExportIDs_] = SourceGIDs[i];
      ExportLIDs_[NumExportIDs_++] = i;
    }
  }

  if ( NumExportIDs_>0 && !sourceMap.DistributedGlobal())
    ReportError("Warning in Epetra_Export: Serial Export has remote IDs. (Exporting from Subset of Source Map)", 1);

  // Test for distributed cases
  int ierr = 0;

  if (sourceMap.DistributedGlobal()) {

    if (NumExportIDs_>0) ExportPIDs_ = new int[NumExportIDs_];
    ierr = targetMap.RemoteIDList(NumExportIDs_, ExportGIDs, ExportPIDs_, 0); // Get remote PIDs
    if( ierr ) throw ReportError("Error in Epetra_BlockMap::RemoteIDList", ierr);

    //Get rid of IDs not in Target Map
    if(NumExportIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumExportIDs_; ++i )
  if( ExportPIDs_[i] == -1 ) ++cnt;
      if( cnt ) {
  int_type * NewExportGIDs = 0;
  int * NewExportPIDs = 0;
  int * NewExportLIDs = 0;
  int cnt1 = NumExportIDs_-cnt;
  if (cnt1) {
    NewExportGIDs = new int_type[cnt1];
    NewExportPIDs = new int[cnt1];
    NewExportLIDs = new int[cnt1];
  }
  cnt = 0;
  for( i = 0; i < NumExportIDs_; ++i )
    if( ExportPIDs_[i] != -1 ) {
      NewExportGIDs[cnt] = ExportGIDs[i];
      NewExportPIDs[cnt] = ExportPIDs_[i];
      NewExportLIDs[cnt] = ExportLIDs_[i];
      ++cnt;
          }
  assert(cnt==cnt1); // Sanity test
  NumExportIDs_ = cnt;
  delete [] ExportGIDs;
  delete [] ExportPIDs_;
  delete [] ExportLIDs_;
  ExportGIDs = NewExportGIDs;
  ExportPIDs_ = NewExportPIDs;
  ExportLIDs_ = NewExportLIDs;
  ReportError("Warning in Epetra_Export: Source IDs not found in Target Map (Do you want to export from subset of Source Map?)", 1 );
      }
    }

    //Make sure Export IDs are ordered by processor
    Epetra_Util util;

    if(targetMap.GlobalIndicesLongLong()) {
      // FIXME (mfh 11 Jul 2013) This breaks ANSI aliasing rules, if
      // int_type != long long.  On some compilers, it results in
      // warnings such as this: "dereferencing type-punned pointer
      // will break strict-aliasing rules".
      util.Sort(true,NumExportIDs_,ExportPIDs_,0,0,1,&ExportLIDs_, 1, (long long **)&ExportGIDs);
    }
    else if(targetMap.GlobalIndicesInt()) {
      int* ptrs[2] = {ExportLIDs_, (int*) ExportGIDs};
      util.Sort(true,NumExportIDs_,ExportPIDs_,0,0, 2,&ptrs[0], 0, 0);
    }
    else {
      throw ReportError("Epetra_Import::Epetra_Import: GlobalIndices Internal Error", -1);
    }

    Distor_ = sourceMap.Comm().CreateDistributor();

    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.

    ierr = Distor_->CreateFromSends( NumExportIDs_, ExportPIDs_, true, NumRemoteIDs_);
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromSends()", ierr);

    // Use comm plan with ExportGIDs to find out who is sending to us and
    // get proper ordering of GIDs for remote entries
    // (that we will convert to LIDs when done).

    if (NumRemoteIDs_>0) RemoteLIDs_ = new int[NumRemoteIDs_]; // Allocate space for LIDs in target that are
    // going to get something from off-processor.
    char * cRemoteGIDs = 0; //Do will alloc memory for this object
    int LenCRemoteGIDs = 0;
    ierr = Distor_->Do(reinterpret_cast<char *> (ExportGIDs),
    sizeof( int_type ),
    LenCRemoteGIDs,
    cRemoteGIDs);
    if (ierr) throw ReportError("Error in Epetra_Distributor.Do()", ierr);
    int_type * RemoteGIDs = reinterpret_cast<int_type*>(cRemoteGIDs);

    // Remote IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumRemoteIDs_; i++) {
      RemoteLIDs_[i] = targetMap.LID(RemoteGIDs[i]);
      //NumRecv_ += targetMap.ElementSize(RemoteLIDs_[i]); // Count total number of entries to receive
      NumRecv_ += targetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
    }

    if (LenCRemoteGIDs>0) delete [] cRemoteGIDs;
  }
  if (NumExportIDs_>0) delete [] ExportGIDs;
  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;

  return;
}

// Epetra_Export constructor for a Epetra_BlockMap object
Epetra_Export::Epetra_Export( const Epetra_BlockMap &  sourceMap, const Epetra_BlockMap & targetMap)
  : Epetra_Object("Epetra::Export"),
    TargetMap_(targetMap),
    SourceMap_(sourceMap),
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
  if(!targetMap.GlobalIndicesTypeMatch(sourceMap))
    throw ReportError("Epetra_Export::Epetra_Export: GlobalIndicesTypeMatch failed", -1);

  if(targetMap.GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    Construct<int>(sourceMap, targetMap);
#else
    throw ReportError("Epetra_Export::Epetra_Export: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  else if(targetMap.GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    Construct<long long>(sourceMap, targetMap);
#else
    throw ReportError("Epetra_Export::Epetra_Export: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  else
    throw ReportError("Epetra_Export::Epetra_Export: Bad global indices type", -1);
}

//==============================================================================
// Epetra_Export copy constructor
Epetra_Export::Epetra_Export(const Epetra_Export & Exporter)
  : Epetra_Object(Exporter),
     TargetMap_(Exporter.TargetMap_),
    SourceMap_(Exporter.SourceMap_),
    NumSameIDs_(Exporter.NumSameIDs_),
    NumPermuteIDs_(Exporter.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Exporter.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Exporter.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Exporter.NumSend_),
    NumRecv_(Exporter.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Exporter.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Exporter.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Exporter.RemoteLIDs_[i];
  }

  TargetMap().Comm().Barrier();
  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Exporter.ExportLIDs_[i];
      ExportPIDs_[i] = Exporter.ExportPIDs_[i];
    }
  }

  if (Exporter.Distor_!=0) Distor_ = Exporter.Distor_->Clone();

}

//==============================================================================
// Epetra_Export destructor
Epetra_Export::~Epetra_Export()
{
  if( Distor_ != 0 ) delete Distor_;

  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by GSPlan
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}

//==============================================================================
// Epetra_Export pseudo-copy constructor.
Epetra_Export::Epetra_Export(const Epetra_Import& Importer):
  TargetMap_(Importer.SourceMap_), //reverse
  SourceMap_(Importer.TargetMap_),//reverse
  NumSameIDs_(Importer.NumSameIDs_),
  NumPermuteIDs_(Importer.NumPermuteIDs_),
  PermuteToLIDs_(0),
  PermuteFromLIDs_(0),
  NumRemoteIDs_(Importer.NumExportIDs_),//reverse
  RemoteLIDs_(0),
  NumExportIDs_(Importer.NumRemoteIDs_),//reverse
  ExportLIDs_(0),
  ExportPIDs_(0),
  NumSend_(Importer.NumRecv_),//reverse
  NumRecv_(Importer.NumSend_),//revsese
  Distor_(0)
{
  // Reverse the permutes
  if (NumPermuteIDs_ > 0) {
    PermuteToLIDs_   = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (int i = 0; i < NumPermuteIDs_; ++i) {
      PermuteFromLIDs_[i] = Importer.PermuteToLIDs_[i];
      PermuteToLIDs_[i]   = Importer.PermuteFromLIDs_[i];
    }
  }

  // Copy the exports to the remotes
  if (NumRemoteIDs_ > 0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (int i = 0; i < NumRemoteIDs_; ++i) RemoteLIDs_[i] = Importer.ExportLIDs_[i];
  }

  // Copy the remotes to the exports
  if (NumExportIDs_ > 0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (int i = 0; i < NumExportIDs_; ++i) ExportLIDs_[i] = Importer.RemoteLIDs_[i];

    // Extract the RemotePIDs from the Distributor
#ifdef HAVE_MPI
    Epetra_MpiDistributor *D=dynamic_cast<Epetra_MpiDistributor*>(&Importer.Distributor());
    if(!D) throw ReportError("Epetra_Export: Can't have ExportPIDs w/o an Epetra::MpiDistributor.",-1);

    // Get the distributor's data
    const int NumReceives  = D->NumReceives();
    const int *ProcsFrom   = D->ProcsFrom();
    const int *LengthsFrom = D->LengthsFrom();

    // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
    // MpiDistributor so it ought to duplicate that effect.
    int i =0, j = 0;
    for (i = 0, j = 0; i < NumReceives; ++i) {
      const int pid = ProcsFrom[i];
      for (int k = 0; k < LengthsFrom[i]; ++k) {
        ExportPIDs_[j] = pid;
        j++;
      }
    }
#else
    throw ReportError("Epetra_Export: Can't have ExportPIDs w/o an Epetra::MpiDistributor.",-2);
#endif
  }//end NumExportIDs>0

  if (Importer.Distor_!=0) Distor_ = Importer.Distor_->ReverseClone();

}

//=============================================================================
void Epetra_Export::Print(std::ostream & os) const
{
  // mfh 05 Jan 2012: The implementation of Print() I found here
  // previously didn't print much at all, and it included a message
  // saying that it wasn't finished ("Epetra_Export Print needs
  // attention!!!!").  What you see below is a port of
  // Tpetra::Export::print, which does have a full implementation.
  // This should allow a side-by-side comparison of Epetra_Export with
  // Tpetra::Export.

  // If true, then copy the array data and sort it before printing.
  // Otherwise, leave the data in its original order.
  //
  // NOTE: Do NOT sort the arrays in place!  Only sort in the copy.
  // Epetra depends on the order being preserved, and some arrays'
  // orders are coupled.
  const bool sortIDs = true;

  const Epetra_Comm& comm = SourceMap_.Comm();
  const int myRank = comm.MyPID();
  const int numProcs = comm.NumProc();

  if (myRank == 0) {
    os << "Export Data Members:" << std::endl;
  }
  // We don't need a barrier before this for loop, because Proc 0 is
  // the first one to do anything in the for loop anyway.
  for (int p = 0; p < numProcs; ++p) {
    if (myRank == p) {
      os << "Image ID       : " << myRank << std::endl;

      os << "permuteFromLIDs:";
      if (PermuteFromLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> permuteFromLIDs (NumPermuteIDs_);
  std::copy (PermuteFromLIDs_, PermuteFromLIDs_ + NumPermuteIDs_,
       permuteFromLIDs.begin());
  if (sortIDs) {
    std::sort (permuteFromLIDs.begin(), permuteFromLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumPermuteIDs_; ++i) {
    os << permuteFromLIDs[i];
    if (i < NumPermuteIDs_ - 1) {
      os << " ";
    }
  }
  os << "}";
      }
      os << std::endl;

      os << "permuteToLIDs  :";
      if (PermuteToLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> permuteToLIDs (NumPermuteIDs_);
  std::copy (PermuteToLIDs_, PermuteToLIDs_ + NumPermuteIDs_,
       permuteToLIDs.begin());
  if (sortIDs) {
    std::sort (permuteToLIDs.begin(), permuteToLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumPermuteIDs_; ++i) {
    os << permuteToLIDs[i];
    if (i < NumPermuteIDs_ - 1) {
      os << " ";
    }
  }
  os << "}";
      }
      os << std::endl;

      os << "remoteLIDs     :";
      if (RemoteLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> remoteLIDs (NumRemoteIDs_);
  std::copy (RemoteLIDs_, RemoteLIDs_ + NumRemoteIDs_,
       remoteLIDs.begin());
  if (sortIDs) {
    std::sort (remoteLIDs.begin(), remoteLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumRemoteIDs_; ++i) {
    os << remoteLIDs[i];
    if (i < NumRemoteIDs_ - 1) {
      os << " ";
    }
  }
  os << "}";
      }
      os << std::endl;

      // If sorting for output, the export LIDs and export PIDs have
      // to be sorted together.  We can use Epetra_Util::Sort, using
      // the PIDs as the keys to match Tpetra::Export.
      std::vector<int> exportLIDs (NumExportIDs_);
      std::vector<int> exportPIDs (NumExportIDs_);
      if (ExportLIDs_ != NULL) {
  std::copy (ExportLIDs_, ExportLIDs_ + NumExportIDs_, exportLIDs.begin());
  std::copy (ExportPIDs_, ExportPIDs_ + NumExportIDs_, exportPIDs.begin());

  if (sortIDs && NumExportIDs_ > 0) {
    int* intCompanions[1]; // Input for Epetra_Util::Sort().
    intCompanions[0] = Epetra_Util_data_ptr(exportLIDs);
    Epetra_Util::Sort (true, NumExportIDs_, Epetra_Util_data_ptr(exportPIDs),
           0, (double**) NULL, 1, intCompanions, 0, 0);
  }
      }

      os << "exportLIDs     :";
      if (ExportLIDs_ == NULL) {
  os << " NULL";
      } else {
  os << " {";
  for (int i = 0; i < NumExportIDs_; ++i) {
    os << exportLIDs[i];
    if (i < NumExportIDs_ - 1) {
      os << " ";
    }
  }
  os << "}";
      }
      os << std::endl;

      os << "exportImageIDs :";
      if (ExportPIDs_ == NULL) {
  os << " NULL";
      } else {
  os << " {";
  for (int i = 0; i < NumExportIDs_; ++i) {
    os << exportPIDs[i];
    if (i < NumExportIDs_ - 1) {
      os << " ";
    }
  }
  os << "}";
      }
      os << std::endl;

      os << "numSameIDs     : " << NumSameIDs_ << std::endl;
      os << "numPermuteIDs  : " << NumPermuteIDs_ << std::endl;
      os << "numRemoteIDs   : " << NumRemoteIDs_ << std::endl;
      os << "numExportIDs   : " << NumExportIDs_ << std::endl;

      // Epetra keeps NumSend_ and NumRecv_, whereas in Tpetra, these
      // are stored in the Distributor object.  This is why we print
      // them here.
      os << "Number of sends: " << NumSend_ << std::endl;
      os << "Number of recvs: " << NumRecv_ << std::endl;
    } // if my rank is p

    // A few global barriers give I/O a chance to complete.
    comm.Barrier();
    comm.Barrier();
    comm.Barrier();
  } // for each rank p

  // The original implementation printed the Maps first.  We moved
  // printing the Maps to the end, for easy comparison with the output
  // of Tpetra::Export::print().
  if (myRank == 0) {
    os << std::endl << std::endl << "Source Map:" << std::endl << std::flush;
  }
  comm.Barrier();
  SourceMap_.Print(os);
  comm.Barrier();

  if (myRank == 0) {
    os << std::endl << std::endl << "Target Map:" << std::endl << std::flush;
  }
  comm.Barrier();
  TargetMap_.Print(os);
  comm.Barrier();

  if (myRank == 0) {
    os << std::endl << std::endl << "Distributor:" << std::endl << std::flush;
  }
  comm.Barrier();
  if (Distor_ == NULL) {
    if (myRank == 0) {
      os << " is NULL." << std::endl;
    }
  } else {
    Distor_->Print(os); // Printing the Distributor is itself distributed.
  }
  comm.Barrier();
}

//----------------------------------------------------------------------------
Epetra_Export& Epetra_Export::operator=(const Epetra_Export& src)
{
  (void)src;
  //not currently supported
  bool throw_err = true;
  if (throw_err) {
    throw ReportError("Epetra_Export::operator= not supported.",-1);
  }
  return(*this);
}
