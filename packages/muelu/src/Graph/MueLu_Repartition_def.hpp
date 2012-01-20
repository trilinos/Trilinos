#ifndef MUELU_REPARTITION_DEF_HPP
#define MUELU_REPARTITION_DEF_HPP

#include "MueLu_Repartition_decl.hpp"


#include "MueLu_Level.hpp"
#include "Teuchos_Hashtable.hpp"
#include <Teuchos_DefaultMpiComm.hpp>
// #include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Repartition<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Repartition()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Repartition<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Repartition() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Repartition<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    // TODO: declare input for factory
    //currentLevel.DeclareInput(varName_,factory_,this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Repartition<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {

    using Teuchos::Array;
    using Teuchos::ArrayRCP;

    // *********************************************************************
    // Step 1: Determine the global size of each partition.
    // *********************************************************************

    // Length of vector "decomposition" is local number of DOFs.  Its entries are partition numbers each DOF belongs to.
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = currentLevel.Get<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition");
    RCP<const Teuchos::Comm<int> > comm = decomposition->getMap()->getComm();
    GO numPartitions = currentLevel.Get<GO>("number of partitions");

    int mypid = comm->getRank();

    // Use a hashtable to record how many local rows belong to each partition.
    RCP<Teuchos::Hashtable<GO,GO> > hashTable;
    hashTable = rcp(new Teuchos::Hashtable<GO,GO>(numPartitions + numPartitions/2));
    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    for (int i=0; i<decompEntries.size(); ++i) {
      if (hashTable->containsKey(decompEntries[i])) {
        GO count = hashTable->get(decompEntries[i]);
        ++count;
        hashTable->put(decompEntries[i],count);
        //std::cout << "pid " << mypid << " found " << decompEntries[i] << ", count now " << count << std::endl;
      } else {
        hashTable->put(decompEntries[i],1);
        //std::cout << "pid " << mypid << " did not find " << decompEntries[i] << ", count now 1" << std::endl;
      }
    }
    decompEntries = Teuchos::null;

    Teuchos::Array<GO> allPartitionsIContributeTo;
    Teuchos::Array<GO> allLocalPartSize;
    hashTable->arrayify(allPartitionsIContributeTo,allLocalPartSize);

/*
//FIXME debugging output
    sleep(1); comm->barrier();
    if (mypid == 0) std::cout << "hash table entries\n====================" << std::endl;
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        for (int j=0; j<allPartitionsIContributeTo.size(); ++j)
          std::cout << "pid " << mypid << " has " << allLocalPartSize[j] << " entries for partition "
                    << allPartitionsIContributeTo[j] << std::endl;;
      }
      comm->barrier();
    }
    sleep(1); comm->barrier();
*/

    GO indexBase = decomposition->getMap()->getIndexBase();

    // Source map is overlapping.  GIDs owned by this pid are the partition numbers found above.
    RCP<Map> sourceMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                           allPartitionsIContributeTo(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

    // Store # of local DOFs in each partition in a vector based on above map.
    RCP<Xpetra::Vector<GO,LO,GO,NO> > localPartSizeVec = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(sourceMap,false);
    ArrayRCP<GO> data;
    if (localPartSizeVec->getLocalLength() > 0)
      data = localPartSizeVec->getDataNonConst(0);
    for (int i=0; i<hashTable->size(); ++i)
      data[i] = allLocalPartSize[i];
    data = Teuchos::null;

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(-1);

    // Target map is nonoverlapping.  Pid k has GID N if and only if k owns partition N.
    GO myPartitionNumber;
    Array<int> partitionOwners;
    DeterminePartitionPlacement(currentLevel,myPartitionNumber,partitionOwners);

    GO numDofsThatStayWithMe=0;
    Teuchos::Array<GO> partitionsIContributeTo;
    Teuchos::Array<GO> localPartSize;
    for (int i=0; i<allPartitionsIContributeTo.size(); ++i)
    {
      if (allPartitionsIContributeTo[i] != myPartitionNumber) {
        partitionsIContributeTo.push_back(allPartitionsIContributeTo[i]);
        localPartSize.push_back(allLocalPartSize[i]);
      }
      else 
        numDofsThatStayWithMe = allLocalPartSize[i];
    }

    GO numPartitionsISendTo = hashTable->size();
    // I'm a partition owner, so don't count my own.
    if (myPartitionNumber >= 0 && numPartitionsISendTo>0) numPartitionsISendTo--;
    //sanity check
    std::cout << "pid " << mypid << " : numPartitionsISendTo = " << numPartitionsISendTo 
                                 << ", partitionsIContributeTo.size() = " << partitionsIContributeTo.size()
                                 << std::endl;
    assert(numPartitionsISendTo == partitionsIContributeTo.size());

    std::cout << mypid << " here 1" << std::endl; sleep(1); comm->barrier();

    Array<int> partitionOwnersISendTo;
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      //if (partitionsIContributeTo[i] != myPartitionNumber)
      partitionOwnersISendTo.push_back(partitionOwners[partitionsIContributeTo[i]]);
    }

    /*
       for (int j=0; j<comm->getSize(); ++j) {
       if (j==mypid) {
       std::cout << "pid " << mypid << " ";
       if (myPartitionNumber >= 0) std::cout << "owns partition " << myPartitionNumber << "." << std::endl;
       else                        std::cout << "does not own a partition." << std::endl;
       std::cout << "     " << "global partition list (partition #, owner)" << std::endl;
       for (int i=0; i<partitionOwners.size(); ++i)
       std::cout << "(" << i << ", " << partitionOwners[i] << ")" << std::endl;
       std::cout << std::endl;
       }
       comm->barrier();
       }
     */

    std::cout << mypid << " here 2" << std::endl; sleep(1); comm->barrier();

    Array<GO> localMapElement;
    if (myPartitionNumber >= 0)
            localMapElement.push_back(myPartitionNumber);
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                    numPartitions,
                    localMapElement(),
                    decomposition->getMap()->getIndexBase(),
                    comm);

    std::cout << mypid << " here 3" << std::endl; sleep(1); comm->barrier();

    RCP<const Export> exporter = ExportFactory::Build( sourceMap,targetMap);

    // If this pid owns a partition, globalPartSizeVec has one local entry that is the global size of said partition.
    // If this pid doesn't own a partition, globalPartSizeVec is locally empty.
    RCP<Xpetra::Vector<GO,LO,GO,NO> > globalPartSizeVec = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(targetMap);
    globalPartSizeVec->doExport(*localPartSizeVec,*exporter,Xpetra::ADD);
    int myPartitionSize = 0;
    if (globalPartSizeVec->getLocalLength() > 0) {
            data = globalPartSizeVec->getDataNonConst(0);
            myPartitionSize = data[0];
    }
    sleep(1); comm->barrier();
    if (comm->getRank() == 0) std::cout << "=========\npartition size\n==========" << std::endl;
    comm->barrier();
    globalPartSizeVec->describe(*fos,Teuchos::VERB_EXTREME);
    sleep(1);
    comm->barrier();

    // Calculate how many PIDs (other than myself) contribute to my partition.
    RCP<Xpetra::Vector<GO,LO,GO,NO> > howManyPidsSendToThisPartitionVec = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(targetMap);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > partitionsISendTo = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(sourceMap,false);
    if (partitionsISendTo->getLocalLength() > 0)
      data = partitionsISendTo->getDataNonConst(0);
    for (int i=0; i<data.size(); ++i) {
      // don't count myself as someone I send to... (sourceMap is based on allPartitionsIContributeTo)
      // TODO last thing I changed!!
      if (sourceMap->getGlobalElement(i) != myPartitionNumber)
        data[i] = 1;
      else 
        data[i] = 0;
    }
    data = Teuchos::null;
    int howManyPidsSendToThisPartition = 0;
    howManyPidsSendToThisPartitionVec->doExport(*partitionsISendTo,*exporter,Xpetra::ADD);
    if (howManyPidsSendToThisPartitionVec->getLocalLength() > 0) {
      data = howManyPidsSendToThisPartitionVec->getDataNonConst(0);
      howManyPidsSendToThisPartition = data[0];
    }
    sleep(1); comm->barrier();
    if (comm->getRank() == 0) std::cout << "=========\n#pids that contribute to each partition\n==========" << std::endl;
    comm->barrier();
    //howManyPidsSendToThisPartitionVec->describe(*fos,Teuchos::VERB_EXTREME);

    //FIXME debugging output
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "partition owner " << mypid << " receiving from " << howManyPidsSendToThisPartition
                  << " pids" << std::endl;;
      }
      comm->barrier();
    }
    sleep(1);
    comm->barrier();

    //
    // Calculate which PIDs contribute to my partition, and how much each contributes,
    // with ireceive/send/wait cycle.
    // FIXME Jan.12.2012 Teuchos::Comm methods ireceive and wait don't work (bugs 5483 and 5484), so for
    // now use the raw MPI methods.

    // First post non-blocking receives.
    ArrayRCP<GO> pidsIReceiveFrom(howManyPidsSendToThisPartition);
    ArrayRCP<GO> numDofsIReceiveFromOnePid(howManyPidsSendToThisPartition);
    for (int j=0; j<numDofsIReceiveFromOnePid.size(); ++j)
      numDofsIReceiveFromOnePid[j] = -99;
    //ArrayRCP<char> recvBuffer = Teuchos::arcp_reinterpret_cast<char>(numDofsIReceiveFromOnePid);

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    if (tmpic == Teuchos::null)
      throw(Exceptions::RuntimeError("Cannot cast base Teuchos::Comm to Teuchos::MpiComm object."));
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
    ArrayRCP<MPI_Request> requests(howManyPidsSendToThisPartition);
    //std::cout << "pid " << mypid << " posting irecv" << std::endl; comm->barrier();
    //for (int j=0; j<comm->getSize(); ++j) {
    //  if (j==mypid) {
    //  std::cout << "pid " << mypid << " posting " << howManyPidsSendToThisPartition << " receives" << std::endl;
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      MPI_Irecv((void*)&(numDofsIReceiveFromOnePid[i]),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,*rawMpiComm,&requests[i]);
      /*
      //this does NOT work
      ArrayView<int> ttt(numDofsIReceiveFromOnePid.view(i,1));
      ArrayView<char> recvBuffer = Teuchos::av_reinterpret_cast<char>(ttt);
      comm->ireceive(recvBuffer,-1);
      */
      /*
      //this does NOT work either
      ArrayView<int> ttt(&numDofsIReceiveFromOnePid[i],1);
      ArrayView<char> recvBuffer = Teuchos::av_reinterpret_cast<char>(ttt);
      comm->ireceive(recvBuffer,-1);
      */
      //comm->ireceive(recvBuffer.view(i,sizeof(int)/sizeof(char)),-1);
      // comm->ireceive(recvBuffer(),-1);
    }
    //} //if (j==mypid)
    //comm->barrier();
    //} //for (int j=0; j<comm->getSize(); ++j)
    //sleep(1);  comm->barrier();

    // Next post sends.

    //for (int j=0; j<comm->getSize(); ++j) {
    //  if (j==mypid)
    //  {
    //    std::cout << "pid " << mypid << " posting sends" << std::endl; comm->barrier();
        for (int i=0; i< partitionsIContributeTo.size(); ++i) {
          //std::cout << "     " << " sending size " << localPartSize[i] << " to " <<  partitionOwnersISendTo[i] << std::endl;
          //if (partitionsIContributeTo[i] != myPartitionNumber) {
          comm->send(sizeof(GO), (char*)&localPartSize[i], partitionOwnersISendTo[i]);
          //}
        }
    //  }
    //  comm->barrier();
    //}
    //sleep(1);  comm->barrier();

    // Finally do waits.

    ArrayRCP<MPI_Status> status(howManyPidsSendToThisPartition);
    //std::cout << "pid " << mypid << " posting waits" << std::endl; comm->barrier();
    for (int i=0; i<howManyPidsSendToThisPartition; ++i)
      MPI_Wait(&requests[i],&status[i]);

    //std::cout << "pid " << mypid << " done waiting" << std::endl; comm->barrier();
    for (int i=0; i<pidsIReceiveFrom.size(); ++i)
      pidsIReceiveFrom[i] = status[i].MPI_SOURCE;

    // =================================================================================================
    // Calculate partial offsets for permutation row map, via MPI_Scan based on global partition sizes.
    // Communicate those offsets back to respective PIDS of unpermuted matrix using ireceive/send/wait cycle.
    // =================================================================================================
    // Partition numbers are used as tags to ensure messages go to correct PIDs.
    // Use MPI_DOUBLE type to avoid any overflow problems.
    // partitionSizeOffset is the first GID in this partition.
    // gidOffsets is an array of GID offsets.  gidOffsets[i] is the first GID for the dofs received from PID partitionOwnersISendTo[i].
    // Note: Quantities "numPartitionsISendTo" and "howManyPidsSendToThisPartition" do not include me
    double partitionSizeOffset;
    double ttt = myPartitionSize;
    MPI_Scan(&ttt, &partitionSizeOffset, 1, MPI_DOUBLE, MPI_SUM, *rawMpiComm);
    partitionSizeOffset -= myPartitionSize;

    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << " : partition GID offset = " << partitionSizeOffset << std::endl;
      }
      comm->barrier();
    }

    ArrayRCP<double> gidOffsets;
    if (howManyPidsSendToThisPartition > 0) {
      gidOffsets = ArrayRCP<double>(howManyPidsSendToThisPartition);
      gidOffsets[0] = partitionSizeOffset;
    }
    for (int i=1; i<howManyPidsSendToThisPartition; ++i) {
      gidOffsets[i] = gidOffsets[i-1] + numDofsIReceiveFromOnePid[i-1];
    }

    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << " : offsets I send to contributors " << std::endl;
        for (int j=0; j<howManyPidsSendToThisPartition; ++j) {
          std::cout << "     to pid " << pidsIReceiveFrom[j] << ", offset " << gidOffsets[j] << std::endl;
        }
      }
      comm->barrier();
    }

    // Send GID offsets back to their respective PIDs.

    requests = ArrayRCP<MPI_Request>(numPartitionsISendTo);
    ArrayRCP<double> gidOffsetsForPartitionsIContributeTo(numPartitionsISendTo);

    /*
       Quantities "numPartitionsISendTo" and "howManyPidsSendToThisPartition" do not include me
     */

    // First post receives on contributing PIDs
    for (int i=0; i< numPartitionsISendTo; ++i) {
      //if (partitionsIContributeTo[i] != myPartitionNumber) {
      int msgTag = partitionsIContributeTo[i];
      MPI_Irecv((void*)&(gidOffsetsForPartitionsIContributeTo[i]),1,MPI_DOUBLE,partitionOwnersISendTo[i],msgTag,*rawMpiComm,&requests[i]);
      //}
    }

    sleep(1); comm->barrier();

    // Now do sends by partition owners
    for (int j=0; j<comm->getSize(); ++j) {
    if (mypid == j) {
    std::cout << "pid " << mypid << " now posting sends" << std::endl;
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      int msgTag = myPartitionNumber;
      MPI_Send((void*)&gidOffsets[i] , 1 , MPI_DOUBLE , pidsIReceiveFrom[i], msgTag, *rawMpiComm);
      std::cout << "    sending offset " << gidOffsets[i] << " to " << pidsIReceiveFrom[i] << std::endl;
    }
    } //if (mypid == j)
    comm->barrier();
    } //for (int j=0; j<comm->getSize())

    sleep(1); comm->barrier();

    status = ArrayRCP<MPI_Status>(numPartitionsISendTo);
    for (int i=0; i<numPartitionsISendTo; ++i)
      MPI_Wait(&requests[i],&status[i]);

    //!FIXME debugging output
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << " received the following GID offsets" << std::endl;
        for (int j=0; j<numPartitionsISendTo; ++j) {
          std::cout << "   " << gidOffsetsForPartitionsIContributeTo[j] << " (from partition " << status[j].MPI_TAG << ")";
          std::cout << "   (partitionsIContributeTo[" << j << "] = " << partitionsIContributeTo[j];
          std::cout << std::endl;
        }
      }
      comm->barrier();
    }
    //!FIXME end of debugging output

    /*
      Now set up a GID scheme that is the same for both the original unpermuted system and the permuted system.
      This will allow us to create the importer we need to transfer column IDs.
    */

    // original unpermuted system

    // store offsets for easy random access
    hashTable = rcp(new Teuchos::Hashtable<GO,GO>(partitionsIContributeTo.size() + partitionsIContributeTo.size()/2));
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      hashTable->put(partitionsIContributeTo[i],(GO)gidOffsetsForPartitionsIContributeTo[i]);
    }
    //store gid offset for those dofs that will remain with me
    if (myPartitionNumber > -1)
      hashTable->put(myPartitionNumber,((GO)partitionSizeOffset) + myPartitionSize - numDofsThatStayWithMe);
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    Array<GO> uniqueGIDsBeforePermute;
    for (int i=0; i<decompEntries.size(); ++i) {
      GO gid = hashTable->get(decompEntries[i]);
      uniqueGIDsBeforePermute.push_back(gid);
      gid++;
      hashTable->put(decompEntries[i],gid);
    }

    if (mypid == 0) {std::cout << "===============\nunique GIDs before permute\n===============" << std::endl;} sleep(1);comm->barrier();
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << "(#dofs that stay with me = " << numDofsThatStayWithMe << ")" << std::endl;
        for (int j=0; j<uniqueGIDsBeforePermute.size(); ++j) {
          std::cout << "   " << uniqueGIDsBeforePermute[j] << "(" << decompEntries[j] << ")" << std::endl;
        }
        std::cout << std::endl;
      }
      comm->barrier();
    }

    sleep(1); comm->barrier();

    // permuted system

    Array<GO> uniqueGIDsAfterPermute;
    for (int i=0; i<myPartitionSize; ++i) {
        uniqueGIDsAfterPermute.push_back((GO)partitionSizeOffset+i);
    }

    if (mypid == 0) {std::cout << "===============\nunique GIDs after permute\n===============" << std::endl;} sleep(1);comm->barrier();
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << "(partition #" << myPartitionNumber << ")" << std::endl;
        for (int j=0; j<uniqueGIDsAfterPermute.size(); ++j) {
          std::cout << "   " << uniqueGIDsAfterPermute[j] << std::endl;
        }
        std::cout << std::endl;
      }
      comm->barrier();
    }

    sleep(1); comm->barrier();

    // Create importer and necessary vectors to communicate column GIDs for the permutation matrix.

    //TODO we should really supply the global size as another sanity check
    sourceMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsBeforePermute(),
                    indexBase,
                    comm);

    targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsAfterPermute(),
                    indexBase,
                    comm);

    RCP<const Import> importer = ImportFactory::Build( sourceMap,targetMap);

    RCP<Xpetra::Vector<GO,LO,GO,NO> > sourceVec = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(sourceMap);
    ArrayRCP<GO> vectorData;
    if (sourceVec->getLocalLength() > 0) {
      vectorData = sourceVec->getDataNonConst(0);
    }

    //load source vector with unpermuted GIDs.
    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A");
    RCP<const Map> originalRowMap = A->getRowMap();

    //sanity check
    assert(vectorData.size() == (GO) originalRowMap->getNodeNumElements());


    for (LO i=0; i<vectorData.size(); ++i) {
      GO gid = originalRowMap->getGlobalElement(i);
      if (gid == (GO) Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid())
        throw(Exceptions::RuntimeError("Encountered an unexpected error with A's rowmap."));
      vectorData[i] = gid;
    }
    vectorData = Teuchos::null;

    RCP<Xpetra::Vector<GO,LO,GO,NO> > targetVec = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(targetMap);
    targetVec->doImport(*sourceVec,*importer,Xpetra::INSERT);

    sleep(1);comm->barrier();
    if (mypid == 0) std::cout << "global column ids for permutation matrix" << std::endl;
    comm->barrier();
    targetVec->describe(*fos,Teuchos::VERB_EXTREME);
    sleep(1);comm->barrier();

    //Build permutation matrix itself.
    RCP<Operator> permutationMatrix = OperatorFactory::Build(targetMap,1);
    Array<SC> matrixEntry(1);
    matrixEntry[0] = 1.0;
    if (targetVec->getLocalLength() > 0) {
      vectorData = targetVec->getDataNonConst(0);
    }
    for (int i=0; i< uniqueGIDsAfterPermute.size(); ++i) {
      permutationMatrix->insertGlobalValues(uniqueGIDsAfterPermute[i], vectorData(i,1),matrixEntry());
    }
    vectorData = Teuchos::null;

    permutationMatrix->fillComplete(A->getDomainMap(),targetMap);

    sleep(1);comm->barrier();
    if (mypid == 0) std::cout << "~~~~~~ permutation matrix ~~~~~~" << std::endl;
    comm->barrier();
    permutationMatrix->describe(*fos,Teuchos::VERB_EXTREME);
    sleep(1);comm->barrier();

  } //Build

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Repartition<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::DeterminePartitionPlacement(Level & currentLevel, GO &myPartitionNumber,
  Array<int> &partitionOwners) const
  {
    //FIXME This currently makes pid i the owner of partition i.  We must have better logic to minimize data movement.
    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A");
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
    GO numPartitions = currentLevel.Get<GO>("number of partitions");
    if (comm->getRank() < numPartitions) {
      myPartitionNumber = comm->getRank();
    } else {
      myPartitionNumber = -1;
    }
    for (int i=0; i<numPartitions; ++i)
      partitionOwners.push_back(i);
  } //CalculatePartitionOwners

} // namespace MueLu

#define MUELU_REPARTITION_SHORT
#endif // MUELU_REPARTITION_DEF_HPP
