/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "TpetraNew_Map.hpp"
#include "TpetraNew_Import.hpp"
#include "TpetraNew_Export.hpp"
#include "Tpetra_Distributor.hpp"
#include <iterator>
#include <sstream>

// #include "Tpetra_CrsGraph.hpp"
// #include "Tpetra_CrsMatrix.hpp"
// #include "Tpetra_MultiVector.hpp"
// #include "Tpetra_Vector.hpp"
// #include "Tpetra_Import_Util.hpp"

namespace {

  void
  getRemotePIDs (const TpetraNew::Import& Importer,
		 Teuchos::Array<int>& RemotePIDs)
  {
    using LO = TpetraNew::Import::local_ordinal_type;
    
    const auto& D = Importer.getDistributor();
    Teuchos::ArrayView<const LO> RemoteLIDs = Importer.getRemoteLIDs();

    // Get the distributor's data
    const size_t NumReceives = D.getNumReceives();
    Teuchos::ArrayView<const int> ProcsFrom = D.getProcsFrom();
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
  
  bool
  checkImportValidity (const TpetraNew::Import& Importer)
  {
    using LO = TpetraNew::Import::local_ordinal_type;
    using GO = TpetraNew::Import::local_ordinal_type;    
    
    auto source = Importer.getSourceMap();
    auto target = Importer.getTargetMap();
    auto comm = source->getComm();

    // For now, do not check validity of a locally replicated source map (just return true)
    if (!source->isDistributed()) return true;

    int global_is_valid=0;
    bool is_valid=true;
 
    // We check validity by going through each ID in the source map one by one, broadcasting the sender's PID and then all
    // receivers check it.
    LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    const int MyPID    = comm->getRank();
    const int NumProcs = comm->getSize();

    GO minSourceGID = source->getMinAllGlobalIndex();
    GO maxSourceGID = source->getMaxAllGlobalIndex();
    GO minTargetGID = target->getMinAllGlobalIndex();
    GO maxTargetGID = target->getMaxAllGlobalIndex();

    std::ostringstream os;

    /***********************************************/
    /*              Check recv side                */
    /***********************************************/

    Teuchos::ArrayView<const LO> permuteTarget = Importer.getPermuteToLIDs();
    Teuchos::ArrayView<const LO> remoteLIDs    = Importer.getRemoteLIDs();
    Teuchos::ArrayView<const LO> exportLIDs    = Importer.getExportLIDs();
    Teuchos::ArrayView<const LO> exportPIDs    = Importer.getExportPIDs();
    Teuchos::Array<int> remotePIDs; getRemotePIDs(Importer,remotePIDs);

    // Generate remoteGIDs
    Teuchos::Array<GO> remoteGIDs(remoteLIDs.size());
    for(size_t i=0; i<(size_t)remoteLIDs.size(); i++) {
      remoteGIDs[i] = target->getGlobalIndex(remoteLIDs[i]);
      if(remoteGIDs[i]<0) {
	os<<MyPID<<"ERROR3: source->getGlobalIndex(remoteLIDs[l]) is invalid GID="<<remoteGIDs[i]<<" LID= "<<remoteLIDs[i]<<std::endl;
	is_valid=false;
      }
    }
    // Generate exportGIDs
    Teuchos::Array<GO> exportGIDs(exportLIDs.size(),-1);
    for(size_t i=0; i<(size_t)exportLIDs.size(); i++) {
      exportGIDs[i] = source->getGlobalIndex(exportLIDs[i]);
      exportGIDs[i]=source->getGlobalIndex(exportLIDs[i]);
      if(exportGIDs[i]<0) {
	os<<MyPID<<"ERROR3: source->getGlobalIndex(exportLIDs[l]) is invalid GID="<<exportGIDs[i]<<" LID= "<<exportLIDs[i]<<std::endl;
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
    for(GO i=minSourceGID; i<maxSourceGID; i++) {
      // Get the (source) owner.
      // Abuse reductions to make up for the fact we don't know the owner is.
      // NOTE: If nobody owns this guy, it we'll get -1.
      LO slid = source->getLocalIndex(i);    
      if(slid == LO_INVALID) TempPID = -1;
      else TempPID = MyPID;
      Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX,TempPID, Teuchos::outArg(OwningPID));

      // Check to see if I have this guy in the target.  If so, make sure I am receiving him from the owner
      LO tlid = target->getLocalIndex(i);    

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


    for(GO i=minTargetGID; i<maxTargetGID; i++) {

      // If I have the target guy, set the proc mask
      LO tlid = target->getLocalIndex(i);    
      LO slid = source->getLocalIndex(i);    

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
	      if (exportPIDs[k] == j && source->getGlobalIndex(exportLIDs[k]) == i) {
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
      Teuchos::Array<GO> myexpgid(NumProcs,-2), yourexpgid(NumProcs,-2);
      Teuchos::Array<int> myexppid(NumProcs,-2), yourexppid(NumProcs,-2);
      if(i<exportGIDs.size()) {
	myexpgid[MyPID] = exportGIDs[i];
	myexppid[MyPID] = exportPIDs[i];
      }
      Teuchos::reduceAll<int,GO>(*comm,Teuchos::REDUCE_MAX,NumProcs, &myexpgid[0],&yourexpgid[0]);
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,NumProcs, &myexppid[0],&yourexppid[0]);
      for(int p=0;p<NumProcs;++p) { // check one to one and onto
	GO cgid = yourexpgid[p];
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

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( ImportExport_new, basic ) {
    using map_type = TpetraNew::Map;    
    using GO = map_type::global_ordinal_type;
    
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid ();
    auto comm = Tpetra::getDefaultComm();
    // create Maps
    auto source = TpetraNew::createContigMap (INVALID, 10, comm);
    auto target = TpetraNew::createContigMap (INVALID,  5, comm);
    auto importer = TpetraNew::createImport (source, target);

    auto same = importer->getNumSameIDs();
    auto permute = importer->getNumPermuteIDs();
    auto remote = importer->getNumRemoteIDs();
    auto sum = same + permute + remote;
    auto expectedSum = target->getMyNumIndices();
    TEST_EQUALITY( sum, expectedSum );

    const bool isvalid = checkImportValidity (*importer);
    TEST_ASSERT( isvalid );
  }

#if 0
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, GetNeighborsForward, Scalar )
  {
    // import with the importer to duplicate
    // export with the exporter to add and reduce
    
    using ST = Teuchos::ScalarTraits<Scalar>;
    using MV = Tpetra::MultiVector<Scalar>;
    using LO = TpetraNew::Map::local_ordinal_type;
    using GO = TpetraNew::Map::global_ordinal_type;    

    auto comm = Tpetra::getDefaultComm();    
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;

    // create a Map    
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();    
    const LO numLocal = 1;
    const LO numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<GO> neighbors;
    if (myImageID != 0) {
      neighbors.push_back(myImageID-1);
    }
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) {
      neighbors.push_back(myImageID+1);
    }
    // two maps: one has one entries per process, the other is the 1-D neighbors
    auto smap = TpetraNew::createContigMap (INVALID, numLocal, comm);
    auto tmap = rcp (new TpetraNew::Map (INVALID, neighbors (), 0, comm));
    
    for (size_t tnum=0; tnum < 2; ++tnum) {
      RCP<MV> mvMine, mvWithNeighbors;
      // for tnum=0, these are contiguously allocated multivectors
      // for tnum=1, these are non-contiguous views of multivectors
      if (tnum == 0) {
        mvMine = rcp(new MV(smap,numVectors));
        mvWithNeighbors = rcp(new MV(tmap,numVectors));
      }
      else {
        MV mineParent(smap,2+numVectors);
	MV neigParent(tmap,2+numVectors);
        TEUCHOS_TEST_FOR_EXCEPTION(numVectors != 5, std::logic_error, "Test assumption broken.");
        mvMine = mineParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
        mvWithNeighbors = neigParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
      }
      // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
      for (LO j=0; j<numVectors; ++j) {
        mvMine->replaceLocalValue(0,j,static_cast<Scalar>(myImageID + j*numImages));
      }
      // create Import from smap to tmap, Export from tmap to smap, test them
      RCP<const Import<LO,GO,Node> > importer =
        Tpetra::createImport<LO,GO,Node>(smap,tmap);
      RCP<const Export<LO,GO,Node> > exporter =
        Tpetra::createExport<LO,GO,Node>(tmap,smap);
      bool local_success = true;
      // importer testing
      TEST_ASSERT( importer->getSourceMap() == smap );
      TEST_ASSERT( importer->getTargetMap() == tmap );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      // exporter testing
      TEST_ASSERT( exporter->getSourceMap() == tmap );
      TEST_ASSERT( exporter->getTargetMap() == smap );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      // import neighbors, test their proper arrival
      //                   [ 0    n     2n    3n    4n ]
      // mvWithNeighbors = [...  ....  ....  ....  ....]
      //                   [n-1  2n-1  3n-1  4n-1  5n-1]
      mvWithNeighbors->doImport(*mvMine,*importer,REPLACE);
      if (myImageID == 0) {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)); // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(j*numImages)+ST::one()); // neighbor
        }
      }
      else if (myImageID == numImages-1) {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
        }
      }
      else {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),2,static_cast<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
        }
      }

      // export values, test
      mvMine->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      mvMine->doExport(*mvWithNeighbors,*exporter,ADD);
      if (myImageID == 0 || myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and one neighbor: double original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(2.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and two neighbors: triple original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(3.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      success &= local_success;
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, GetNeighborsBackward, Scalar )
  {
    // import with the exporter to duplicate
    // export with the importer to add and reduce

    using ST = Teuchos::ScalarTraits<Scalar>;
    using MV = Tpetra::MultiVector<Scalar>;
    using LO = TpetraNew::Map::local_ordinal_type;
    using GO = TpetraNew::Map::global_ordinal_type;    

    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    
    // create a Map    
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    const LO numLocal = 1;
    const LO numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<GO> neighbors;
    if (myImageID != 0) {
      neighbors.push_back(myImageID-1);
    }
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) {
      neighbors.push_back(myImageID+1);
    }
    // two maps: one has one entries per process, the other is the 1-D neighbors
    auto smap = TpetraNew::createContigMap (INVALID, numLocal, comm);
    auto tmap = rcp (new TpetraNew::Map (INVALID, neighbors (), 0, comm));
    for (size_t tnum=0; tnum < 2; ++tnum) {
      RCP<MV> mvMine, mvWithNeighbors;
      // for tnum=0, these are contiguously allocated multivectors
      // for tnum=1, these are non-contiguous views of multivectors
      if (tnum == 0) {
        mvMine = rcp(new MV(smap,numVectors));
        mvWithNeighbors = rcp(new MV(tmap,numVectors));
      }
      else {
        MV mineParent(smap,2+numVectors);
	MV neigParent(tmap,2+numVectors);
        TEUCHOS_TEST_FOR_EXCEPTION(numVectors != 5, std::logic_error, "Test assumption broken.");
        mvMine = mineParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
        mvWithNeighbors = neigParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
      }
      // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
      for (LO j=0; j<numVectors; ++j) {
        mvMine->replaceLocalValue(0,j,static_cast<Scalar>(myImageID + j*numImages));
      }
      // create Import from smap to tmap, Export from tmap to smap, test them
      auto importer = Tpetra::createImport<LO, GO, Node> (smap, tmap);
      auto exporter = Tpetra::createExport<LO, GO, Node> (tmap, smap);
      bool local_success = true;
      // importer testing
      TEST_EQUALITY_CONST( importer->getSourceMap() == smap, true );
      TEST_EQUALITY_CONST( importer->getTargetMap() == tmap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      // exporter testing
      TEST_EQUALITY_CONST( exporter->getSourceMap() == tmap, true );
      TEST_EQUALITY_CONST( exporter->getTargetMap() == smap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      // import neighbors, test their proper arrival
      //                   [ 0    n     2n    3n    4n ]
      // mvWithNeighbors = [...  ....  ....  ....  ....]
      //                   [n-1  2n-1  3n-1  4n-1  5n-1]
      mvWithNeighbors->doImport(*mvMine,*exporter,REPLACE);
      if (myImageID == 0) {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)); // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(j*numImages)+ST::one()); // neighbor
        }
      }
      else if (myImageID == numImages-1) {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
        }
      }
      else {
        for (LO j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),2,static_cast<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
        }
      }
      // export values, test
      mvMine->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      mvMine->doExport(*mvWithNeighbors,*importer,ADD);
      if (myImageID == 0 || myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and one neighbor: double original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(2.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and two neighbors: triple original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(3.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      success &= local_success;
    }
    //
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST( ImportExport, AbsMax )
  {
    // test ABSMAX CombineMode
    // test with local and remote entries, as copyAndPermute() and unpackAndCombine() both need to be tested
    
    using TpetraNew::createContigMap;
    using TpetraNew::createImport;    
    using TpetraNew::createNonContigMap;
    using Teuchos::RCP;
    using Teuchos::rcp;    
    using Teuchos::tuple;
    using SC = Tpetra::Vector<>::scalar_type;
    using Vec = Tpetra::Vector<SC>;
    using GO = TpetraNew::Map::global_ordinal_type;    
    
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;

    auto smap = createContigMap (INVALID, 1, comm);
    const GO myOnlyGID = smap->getGlobalIndex (0);
    auto dmap = createNonContigMap (tuple<GO> (myOnlyGID, (myOnlyGID+1) % numImages), comm);
    RCP<Vec> srcVec = rcp (new Vec (smap));
    srcVec->putScalar (-1.0);
    RCP<Vec> dstVec = rcp (new Vec (dmap));
    dstVec->putScalar (-3.0);
    // first item of dstVec is local (w.r.t. srcVec), while the second is remote
    // ergo, during the import:
    // - the first will be over-written (by 1.0) from the source, while
    // - the second will be "combined", i.e., abs(max(1.0,3.0)) = 3.0 from the dest
    auto importer = createImport (smap, dmap);
    dstVec->doImport (*srcVec,*importer,Tpetra::ABSMAX);
    TEST_COMPARE_ARRAYS( tuple<SC>(-1.0,3.0), dstVec->get1dView() )
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, ExportReverse )
  {
    // This test reproduces an issue seen in Github Issue #114.
    // As of time of checkin, this test will fail on CUDA but pass on other platforms.
    // This is intentional
    using GO = TpetraNew::Map::global_ordinal_type;    
    typedef TpetraNew::Map Tpetra_Map;
    typedef TpetraNew::Import Tpetra_Import;
    typedef Tpetra::Vector<int> IntVector;
    
    auto comm = Tpetra::getDefaultComm();
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();

    int NumProcs = comm->getSize();
    int MyPID    = comm->getRank();

    // This problem only works on 4 procs
    if(NumProcs!=4) {TEST_EQUALITY(true,true);return;}

    // Problem setup
    int num_per_proc;
    if(MyPID==0) num_per_proc=7;
    else num_per_proc=6;

    GO from_gids_p0[7] = {0,1,2,3,4,5,6};
    GO to_gids_p0[7]   = {0,4,8,12,16,20,24};

    GO from_gids_p1[6] = {7,8,9,10,11,12};
    GO to_gids_p1[6]   = {1,5,9,13,17,21};

    GO from_gids_p2[6] = {13,14,15,16,17,18};
    GO to_gids_p2[6]   = {2,6,10,14,18,22};

    GO from_gids_p3[6] = {19,20,21,22,23,24};
    GO to_gids_p3[6]   = {3,7,11,15,19,23};

    // Correctness check array
    int who_owns[25];
    for(int i=0; i<7; i++)
      who_owns[to_gids_p0[i]] = 0;
    for(int i=0; i<6; i++) {
      who_owns[to_gids_p1[i]] = 1;
      who_owns[to_gids_p2[i]] = 2;
      who_owns[to_gids_p3[i]] = 3;
    }

    GO *from_ptr, *to_ptr;
    if(MyPID==0)      {from_ptr=&from_gids_p0[0]; to_ptr=&to_gids_p0[0];}
    else if(MyPID==1) {from_ptr=&from_gids_p1[0]; to_ptr=&to_gids_p1[0];}
    else if(MyPID==2) {from_ptr=&from_gids_p2[0]; to_ptr=&to_gids_p2[0];}
    else if(MyPID==3) {from_ptr=&from_gids_p3[0]; to_ptr=&to_gids_p3[0];}
    else exit(-1);

    Teuchos::ArrayView<GO> myfromgids(from_ptr,num_per_proc);
    Teuchos::ArrayView<GO> mytogids(to_ptr,num_per_proc);

    // FromMap (from.getRowMap() from Zoltan2)
    RCP<Tpetra_Map> FromMap = rcp(new Tpetra_Map(INVALID,myfromgids,0,comm));

    // ToMap (tmap from Zoltan2)
    RCP<Tpetra_Map> ToMap = rcp(new Tpetra_Map(INVALID,mytogids,0,comm));

    // Importer
    Tpetra_Import Importer(FromMap,ToMap);

    // Duplicating what Zoltan2/Tpetra Does
    IntVector FromVector(FromMap);
    IntVector ToVector(ToMap);
    ToVector.putScalar(MyPID);
    FromVector.putScalar(-666);

    FromVector.doExport(ToVector,Importer,Tpetra::REPLACE);

    Teuchos::ArrayRCP<const int> f_rcp = FromVector.getData();
    Teuchos::ArrayView<const int> f_view = f_rcp();
    Teuchos::ArrayRCP<const int> t_rcp = ToVector.getData();
    Teuchos::ArrayView<const int> t_view = t_rcp();

    // Check the "FromAnswer" answer against who_owns
    bool all_is_well=true;
    for(size_t i=0; i<FromMap->getNodeNumElements(); i++) {
      if(f_view[i] != who_owns[FromMap->getGlobalIndex(i)]) {
        std::cerr<<"["<<MyPID<<"] ERROR: Ownership of GID"<<FromMap->getGlobalIndex(i)<<" is incorrect!"<<std::endl;
        all_is_well=false;
      }
    }
    TEST_EQUALITY(all_is_well,true);

    
    const bool isvalid = checkImportValidity(Importer);
    if (!isvalid) {
      std::ostringstream oss;
      Importer.print(oss);
      std::cout<<oss.str()<<std::endl;
    }

    TEST_EQUALITY(isvalid,true);
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_S( SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ImportExport, GetNeighborsForward,  SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ImportExport, GetNeighborsBackward, SC )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_S )

#endif // 0

} // namespace (anonymous)


