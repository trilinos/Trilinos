// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_IMPORT_UTIL_HPP
#define TPETRA_IMPORT_UTIL_HPP

/// \file Tpetra_Import_Util.hpp
/// \brief Internal functions and macros designed for use with
///   Tpetra::Import and Tpetra::Export objects.
/// \warning The functions in this file are implementation details of
///   Tpetra.  We make no promises of backwards compatibility.

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

namespace Tpetra {
  namespace Import_Util {
    /// \brief For each GID in the TargetMap, find who owns the GID in the SourceMap.
    ///
    /// This only uses the Distributor and does not communicate.  It
    /// returns (as an output argument) an array of (PID,GID) pairs.
    /// If use_minus_one_for_local is true, any GIDs owned by this
    /// processor get -1 instead of their PID.
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void
    getPidGidPairs (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
                    Teuchos::Array< std::pair<int,GlobalOrdinal> >& gpids,
                    bool use_minus_one_for_local);

    //! Like getPidGidPairs, but just gets the PIDs, ordered by the column Map.
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void
    getPids (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
             Teuchos::Array<int>& pids,
             bool use_minus_one_for_local);

    //! Like getPidGidPairs, but just gets the PIDs, ordered by the column Map.
    // Like the above, but without the resize
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void
    getPids (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
             Teuchos::ArrayView<int>& pids,
             bool use_minus_one_for_local);


    /// \brief Get a list of remote PIDs from an importer in the order
    ///   corresponding to the remote LIDs.
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void
    getRemotePIDs (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
                   Teuchos::Array<int>& RemotePIDs);

  } // namespace Import_Util
} // namespace Tpetra

namespace Tpetra {
namespace Import_Util {

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
getPidGidPairs (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
                Teuchos::Array< std::pair<int,GlobalOrdinal> >& gpids,
                bool use_minus_one_for_local)
{
  // Put the (PID,GID) pair in member of Importer.TargetMap() in
  // gpids.  If use_minus_one_for_local==true, put in -1 instead of
  // MyPID.
  const Tpetra::Distributor& D = Importer.getDistributor();

  LocalOrdinal ii;
  size_t  i,j,k;
  int mypid = Importer.getTargetMap()->getComm()->getRank();
  size_t N  = Importer.getTargetMap()->getLocalNumElements();

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getProcsFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Resize the outgoing data structure
  gpids.resize(N);

  // Start by claiming that I own all the data
  LocalOrdinal lzero = Teuchos::ScalarTraits<LocalOrdinal>::zero();
  if(use_minus_one_for_local)
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) gpids[ii]=std::make_pair(-1,Importer.getTargetMap()->getGlobalElement(ii));
  else
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) gpids[ii]=std::make_pair(mypid,Importer.getTargetMap()->getGlobalElement(ii));

  // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
  // MpiDistributor so it ought to duplicate that effect.
  for(i=0,j=0; i<NumReceives; i++){
    int pid=ProcsFrom[i];
    for(k=0; k<LengthsFrom[i]; k++){
      if(pid!=mypid) gpids[RemoteLIDs[j]].first=pid;
      j++;
    }
  }
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
getPids (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
         Teuchos::Array<int>& pids,
         bool use_minus_one_for_local)
{
  // Resize the outgoing data structure
  pids.resize(Importer.getTargetMap()->getLocalNumElements());
  Teuchos::ArrayView<int> v_pids = pids();
  getPids(Importer,v_pids,use_minus_one_for_local);
}


template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
getPids (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
         Teuchos::ArrayView<int>& pids,
         bool use_minus_one_for_local)
{
  const Tpetra::Distributor & D=Importer.getDistributor();

  LocalOrdinal ii;
  size_t  i,j,k;
  int mypid = Importer.getTargetMap()->getComm()->getRank();
  size_t N  = Importer.getTargetMap()->getLocalNumElements();
  if(N!=(size_t)pids.size()) throw std::runtime_error("Tpetra::Import_Util::getPids(): Incorrect size for output array");

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getProcsFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Start by claiming that I own all the data
  LocalOrdinal lzero = Teuchos::ScalarTraits<LocalOrdinal>::zero();
  if(use_minus_one_for_local)
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) pids[ii]=-1;
  else
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) pids[ii]=mypid;

  // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
  // MpiDistributor so it ought to duplicate that effect.
  for(i=0,j=0; i<NumReceives; i++){
    int pid=ProcsFrom[i];
    for(k=0; k<LengthsFrom[i]; k++){
      if(pid!=mypid) pids[RemoteLIDs[j]]=pid;
      j++;
    }
  }
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
getRemotePIDs (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer,
               Teuchos::Array<int>& RemotePIDs)
{
  const Tpetra::Distributor& D = Importer.getDistributor();

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getProcsFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Resize the outgoing data structure
  RemotePIDs.resize(Importer.getNumRemoteIDs());

  // Now, for each remote ID, record who actually owns it.  This loop
  // follows the operation order in the MpiDistributor so it ought to
  // duplicate that effect.
  size_t i,j,k;
  for (i = 0, j = 0; i < NumReceives; ++i) {
    const int pid = ProcsFrom[i];
    for (k = 0; k < LengthsFrom[i]; ++k) {
      RemotePIDs[j] = pid;
      j++;
    }
  }
}

 
/* Check some of the validity of an Import object
   WARNING: This is a debugging routine only. */
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool
checkImportValidity (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer)
{
  using Teuchos::RCP;
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > source = Importer.getSourceMap();
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > target = Importer.getTargetMap();
  RCP<const Teuchos::Comm<int> > comm = source->getComm();

  // For now, do not check validity of a locally replicated source map (just return true)
  if (!source->isDistributed()) return true;

  int global_is_valid=0;
  bool is_valid=true;
 
  // We check validity by going through each ID in the source map one by one, broadcasting the sender's PID and then all
  // receivers check it.
  LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  const int MyPID    = comm->getRank();
  const int NumProcs = comm->getSize();

  GlobalOrdinal minSourceGID = source->getMinAllGlobalIndex();
  GlobalOrdinal maxSourceGID = source->getMaxAllGlobalIndex();
  GlobalOrdinal minTargetGID = target->getMinAllGlobalIndex();
  GlobalOrdinal maxTargetGID = target->getMaxAllGlobalIndex();

  std::ostringstream os;

  /***********************************************/
  /*              Check recv side                */
  /***********************************************/

  Teuchos::ArrayView<const LocalOrdinal> permuteTarget = Importer.getPermuteToLIDs();
  Teuchos::ArrayView<const LocalOrdinal> remoteLIDs    = Importer.getRemoteLIDs();
  Teuchos::ArrayView<const LocalOrdinal> exportLIDs    = Importer.getExportLIDs();
  Teuchos::ArrayView<const LocalOrdinal> exportPIDs    = Importer.getExportPIDs();
  Teuchos::Array<int> remotePIDs; getRemotePIDs(Importer,remotePIDs);

  // Generate remoteGIDs
  Teuchos::Array<GlobalOrdinal> remoteGIDs(remoteLIDs.size());
  for(size_t i=0; i<(size_t)remoteLIDs.size(); i++) {
    remoteGIDs[i] = target->getGlobalElement(remoteLIDs[i]);
    if(remoteGIDs[i]<0) {
      os<<MyPID<<"ERROR3: source->getGlobalElement(remoteLIDs[l]) is invalid GID="<<remoteGIDs[i]<<" LID= "<<remoteLIDs[i]<<std::endl;
      is_valid=false;
    }
}
  // Generate exportGIDs
  Teuchos::Array<GlobalOrdinal> exportGIDs(exportLIDs.size(),-1);
  for(size_t i=0; i<(size_t)exportLIDs.size(); i++) {
    exportGIDs[i] = source->getGlobalElement(exportLIDs[i]);
    exportGIDs[i]=source->getGlobalElement(exportLIDs[i]);
    if(exportGIDs[i]<0) {
      os<<MyPID<<"ERROR3: source->getGlobalElement(exportLIDs[l]) is invalid GID="<<exportGIDs[i]<<" LID= "<<exportLIDs[i]<<std::endl;
      is_valid=false;
    }
  }
  
 // Zeroth order test: Remote *** GID *** and Export **GID**'s should be disjoint.  
  for( auto &&rgid : remoteGIDs) {
    if(std::find(exportGIDs.begin(),exportGIDs.end(),rgid) != exportGIDs.end()) {
        is_valid = false;
        os<<MyPID<<"ERROR0: Overlap between remoteGIDs and exportGIDs "<<rgid<<std::endl;
      }
  }

  int TempPID , OwningPID;
  for(GlobalOrdinal i=minSourceGID; i<maxSourceGID; i++) {
    // Get the (source) owner.
    // Abuse reductions to make up for the fact we don't know the owner is.
    // NOTE: If nobody owns this guy, it we'll get -1.
    LocalOrdinal slid = source->getLocalElement(i);    
    if(slid == LO_INVALID) TempPID = -1;
    else TempPID = MyPID;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX,TempPID, Teuchos::outArg(OwningPID));

    // Check to see if I have this guy in the target.  If so, make sure I am receiving him from the owner
    LocalOrdinal tlid = target->getLocalElement(i);    

    if(tlid != LO_INVALID) {
      // This guy is in my target map, now to check if I'm receiving him from the owner (which I now know)
      bool is_ok = false;
      
      // This guy is not in the SourceMap at all.  Weird, but acceptable.
      if(OwningPID == -1) continue;

      if (OwningPID == MyPID) {
	// I own this guy
	if((size_t) tlid < Importer.getNumSameIDs()) {
	  // Check sames
	  is_ok = true;
	}
	else {
	  // Check permutes
	  for (size_t j=0; j<(size_t)permuteTarget.size(); j++) {
	    if(tlid == permuteTarget[j]) {
	      is_ok=true; 
	      break;
	    }
	  }
	}
      }
      else {
	// Check remotes
	bool already_hit = false;
	for(size_t j=0; j<(size_t)remoteGIDs.size(); j++) {
	  if(i == remoteGIDs[j]) {
	    // No double hits please
	    if(already_hit) {
	      is_ok=false; 
	      break;
	    }
	    // GID's match:  Do procs?
	    if(OwningPID == remotePIDs[j]) {
	      is_ok = true;
	      already_hit = true;
	    }
	  }
	}
      }
      if(!is_ok) {
        os<<MyPID<<"  ERROR1: GID "<<i<<" should be remoted from PID "<<OwningPID<<" but isn't."<<std::endl;
	is_valid=false;
      }
    }

  }//end for loop

  /***********************************************/
  /*              Check send side                */
  /***********************************************/
  Teuchos::Array<int> local_proc_mask(NumProcs,0), global_proc_mask(NumProcs,0);


  for(GlobalOrdinal i=minTargetGID; i<maxTargetGID; i++) {

    // If I have the target guy, set the proc mask
    LocalOrdinal tlid = target->getLocalElement(i);    
    LocalOrdinal slid = source->getLocalElement(i);    

    if(tlid==LO_INVALID) local_proc_mask[MyPID] = 0;
    else local_proc_mask[MyPID] = 1;

    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,NumProcs, &local_proc_mask[0],&global_proc_mask[0]);


    if(slid !=LO_INVALID) {
      // If I own this unknown on the src I should check to make sure I'm exporting it to each guy in the global_proc_mask who wants it
      for(int j=0; j<NumProcs; j++) {
	if(j==MyPID) continue; // skip self unknowns
	if(global_proc_mask[j]==1) {
	  bool is_ok = false;
	  // This guy needs the unknown
	  bool already_hit = false;
	  for(size_t k=0; k<(size_t)exportPIDs.size(); k++) {
	    if (exportPIDs[k] == j && source->getGlobalElement(exportLIDs[k]) == i) {
	    // No double hits please
	    if(already_hit) {
	      is_ok=false; 
	      break;
	    }
	    else {
	      is_ok=true;
	      already_hit=true;
	    }
	    }
	  }
	  if(!is_ok) {
	    os<<MyPID<<" ERROR2: GID "<<i<<" should be sent to PID "<<j<<" but isn't"<<std::endl;
	    is_valid=false;
	  }
	}
      }
    }
  }
  
  // cbl check that for each of my remote GIDs I receive a corresponding export id. 

  Teuchos::Array<int> proc_num_exports_recv(NumProcs,0);

  Teuchos::Array<int> remoteGIDcount(remoteGIDs.size(),0);

  int allexpsiz=0;
  Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,exportGIDs.size(),  Teuchos::outArg(allexpsiz));
  
  for(int i=0;i<allexpsiz;++i) {
    Teuchos::Array<GlobalOrdinal> myexpgid(NumProcs,-2), yourexpgid(NumProcs,-2);
    Teuchos::Array<int> myexppid(NumProcs,-2), yourexppid(NumProcs,-2);
    if(i<exportGIDs.size()) {
      myexpgid[MyPID] = exportGIDs[i];
      myexppid[MyPID] = exportPIDs[i];
    }
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_MAX,NumProcs, &myexpgid[0],&yourexpgid[0]);
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,NumProcs, &myexppid[0],&yourexppid[0]);
    for(int p=0;p<NumProcs;++p) { // check one to one and onto
      GlobalOrdinal cgid = yourexpgid[p];
      // ignore -2's. 
      if(cgid == -2) continue;
      if(cgid < 0) {
        os<<MyPID<<" ERROR4: received exportGID is invalid "<<cgid<<std::endl;
        is_valid=false;
      }
      bool foundit=false;
      for(size_t k=0;k<(size_t)remoteGIDs.size();++k) {
        if(cgid == remoteGIDs[k] && yourexppid[p] == MyPID ) {
          if(p != remotePIDs[k]) {
            os<<MyPID<<" ERROR5: receive export from wrong pid: got "<<p<<" expected: "<<remotePIDs[k]<<std::endl;
            is_valid = false;
          }
          remoteGIDcount[k]++;
          if(foundit) {
            os<<MyPID<<" ERROR6: found multiple GIDs from correct pid: GID  "<<remoteGIDs[k]<<std::endl;
            is_valid = false;
          }
          foundit = true;
        }
      }
      if(!foundit &&  yourexppid[p] == MyPID ) {
        os<<MyPID<<" ERROR7: receive gid  "<<cgid<<" that is not in my remote gid list, from pid  "<<p<<std::endl;
        is_valid = false;
      }

    }
  }
  // now check that remoteGIDcount is only 1's.
  for(size_t i = 0; i< (size_t) remoteGIDcount.size(); ++i) {
    int rc = remoteGIDcount[i];
    if(rc == 1) continue;
    os<<MyPID<<" ERROR8: my remote at "<<i<<" gid "<<remoteGIDs[i]<<" has count "<<rc<<std::endl;
    is_valid = false;
  }


  // Do a reduction on the final bool status
  Teuchos::reduceAll<int,int> (*comm, Teuchos::REDUCE_MIN,(int)is_valid, Teuchos::outArg(global_is_valid));

  if(!global_is_valid) {
    std::cerr<<os.str()<<std::flush;
    Importer.print(std::cout);
  }

  return global_is_valid>0;
}


/* Check to see that the import object's target map respects AztecOO/ML ordering */
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool
checkAztecOOMLOrdering (const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>& Importer) 
{

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > source = Importer.getSourceMap();
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > target = Importer.getTargetMap();
	
  // Get the (pid, gid) pairs (with locals flagged as -1)
  Teuchos::Array<std::pair<int, GlobalOrdinal> > gpids;
  getPidGidPairs (Importer, gpids, true);

  bool is_ok=true;
  bool is_local = true;
  int last_PID = -1;
  GlobalOrdinal INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  GlobalOrdinal last_GID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
	
  for(size_t i=0; i<(size_t) gpids.size(); i++) {
    int pid = gpids[i].first;
    GlobalOrdinal gid = gpids[i].second;
     
    if(is_local == false && pid == -1) {
      // Out-of-order local
      is_ok = false;
      break;
    }
    else if(pid == -1) {
      // Local, matching PID
      if(source->getGlobalElement(i) != target->getGlobalElement(i)) {
        // GIDs don't match, though the PIDs do
        is_ok = false;
        break;
      }
    }
    else {
      // Off-rank section
      is_local = false;
      if(last_PID == -1) {
        last_PID = pid;
        last_GID = gid;
      }
      else if(pid < last_PID) {
        // PIDs out of order
        is_ok=false;
        break;
      }
      else if(pid == last_PID) {
        if(gid < last_GID) {
          // Same rank, but gids out of order
                is_ok=false;
                break;
        }
      }
      else {
        // New rank
        last_PID = pid;
        last_GID  = gid;
      }
    }
  }
 
  return is_ok;
}


/* Check to see that the import object's target map respects AztecOO/ML ordering */
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool 
checkBlockConsistency(const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map, size_t block_size) {
  bool consistent_block = true;

  for (size_t lid = 0; lid < map.getLocalNumElements (); lid += block_size) {
    auto lbid = floor(double(lid)/block_size);
    auto gbid = floor(double(map.getGlobalElement(lid))/block_size);

    for (size_t lid2 = 1; lid2 < block_size; ++lid2) {
      auto lbid_n = floor(double(lid+lid2)/block_size);
      auto gbid_n = floor(double(map.getGlobalElement(lid+lid2))/block_size);
      if (consistent_block)
        consistent_block = (lbid == lbid_n);
      if (consistent_block)
        consistent_block = (gbid == gbid_n);
    }
  }

  return consistent_block;
}


} // namespace Import_Util
} // namespace Tpetra

#endif // TPETRA_IMPORT_UTIL_HPP
