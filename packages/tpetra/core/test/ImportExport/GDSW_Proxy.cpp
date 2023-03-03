#include "GDSW_Proxy.hpp"
#include "Kokkos_ArithTraits.hpp"

using Teuchos::RCP;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrix(RCP<const CrsMatrix> inputMatrix, 
                   RCP<const Map> outputRowMap, 
                   RCP<CrsMatrix> & outputMatrix)
{
    RCP<Tpetra::Distributor> distributor;
    std::vector<GO> ownedRowGIDs;
    std::vector<LO> localRowsSend, localRowsSendBegin, localRowsRecv, localRowsRecvBegin;
    // ownedRowGIDs = rows of inputMatrix contained in outputRowMap
    // localRowsSend[indicesSend] = localRows of sending process i corresponding to
    //                              localRowsRecv[indicesRecv], where
    //               indicesSend = localRowsSendBegin[i]:localRowsSendBegin[i+1]-1
    //               indicesRecv = localRowsRecvBegin[i]:localRowsRecvBegin[i+1]-1
    constructDistributor(inputMatrix, outputRowMap, distributor, ownedRowGIDs,
                         localRowsSend, localRowsSendBegin,
                         localRowsRecv, localRowsRecvBegin);
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap(outputRowMap, distributor, targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData(inputMatrix, outputRowMap, distributor, targetMapGIDs, 
                          targetMapGIDsBegin, ownedRowGIDs,
                          localRowsSend, localRowsSendBegin, localRowsRecv,
                          localRowsRecvBegin, outputMatrix);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrixFromImporter(RCP<const CrsMatrix> inputMatrix, 
                               RCP<const Import> importer,
                               RCP<CrsMatrix> & outputMatrix)
{
  RCP<const Map> outputRowMap = importer->getTargetMap();
  RCP<const Map> sourceMap = importer->getSourceMap();
  RCP<Tpetra::Distributor> distributor = Teuchos::rcpFromRef(importer->getDistributor());
  Teuchos::ArrayView<const LO> localRowsSend = importer->getExportLIDs();
  Teuchos::ArrayView<const LO> localRowsRecv = importer->getRemoteLIDs();


  // Scansum to get the begin arrays
  size_t numSends = distributor->getNumSends();
  size_t numRecvs = distributor->getNumReceives();
  Teuchos::ArrayView<const size_t> lengthsTo = distributor->getLengthsTo();
  Teuchos::ArrayView<const size_t> lengthsFrom = distributor->getLengthsFrom();
  std::vector<LO> localRowsSendBegin(numSends+1);
  std::vector<LO> localRowsRecvBegin(numRecvs+1);
  for (size_t i=0; i<numSends; i++)
    localRowsSendBegin[i+1] = localRowsSendBegin[i]  + lengthsTo[i];
  for (size_t i=0; i<numRecvs; i++)
    localRowsRecvBegin[i+1] = localRowsRecvBegin[i]  + lengthsFrom[i];


  // Combine sames and permutes for ownedRowGIDs
  // QUESTION: Does iSM handle permutes correctly?
  size_t numSames = importer->getNumSameIDs();
  size_t numPermutes = importer->getNumPermuteIDs();
  Teuchos::ArrayView<const LO> permutes = importer->getPermuteFromLIDs();
  std::vector<GO> ownedRowGIDs(numSames + numPermutes);
  for(size_t i=0; i <numSames; i++)
    ownedRowGIDs[i] = sourceMap->getGlobalElement(i);
  for(size_t i=0; i <numPermutes; i++)
    ownedRowGIDs[i] = numSames + sourceMap->getGlobalElement(permutes[i]);
  


  //    std::vector<LO> localRowsSend, localRowsSendBegin, localRowsRecv, localRowsRecvBegin;
    // ownedRowGIDs = rows of inputMatrix contained in outputRowMap
    // localRowsSend[indicesSend] = localRows of sending process i corresponding to
    //                              localRowsRecv[indicesRecv], where
    //               indicesSend = localRowsSendBegin[i]:localRowsSendBegin[i+1]-1
    //               indicesRecv = localRowsRecvBegin[i]:localRowsRecvBegin[i+1]-1
    //    constructDistributor(inputMatrix, outputRowMap, distributor, ownedRowGIDs,
    //                         localRowsSend, localRowsSendBegin,
    //                         localRowsRecv, localRowsRecvBegin);
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap(outputRowMap, distributor, targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData(inputMatrix, outputRowMap, distributor, targetMapGIDs, 
                          targetMapGIDsBegin, ownedRowGIDs,
                          localRowsSend, localRowsSendBegin, localRowsRecv,
                          localRowsRecvBegin, outputMatrix);
}



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class LOVector>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateMatrixData(RCP<const CrsMatrix> inputMatrix, 
                      RCP<const Map> rowMap,
                      RCP<Tpetra::Distributor> distributor,
                      const std::vector<GO> & targetMapGIDs, 
                      const std::vector<LO> & targetMapGIDsBegin,
                      const std::vector<GO> & ownedRowGIDs,
                      const LOVector & localRowsSend,
                      const std::vector<LO> & localRowsSendBegin,
                      const LOVector & localRowsRecv,
                      const std::vector<LO> & localRowsRecvBegin,
                      RCP<CrsMatrix> & outputMatrix)
{
    const size_t numSends = distributor->getNumSends();
    const size_t numRecvs = distributor->getNumReceives();
    TEUCHOS_TEST_FOR_EXCEPTION(numSends != localRowsSendBegin.size()-1,std::runtime_error,
                    "invalid size of localRowsSendBegin");
    TEUCHOS_TEST_FOR_EXCEPTION(numRecvs != localRowsRecvBegin.size()-1,std::runtime_error,
                    "invalid size of localRowsRecvBegin");
    const size_t numRowsSendTotal = localRowsSend.size();
    std::vector<size_t> count(numRowsSendTotal, 0);
    IndicesViewT indices;
    ValuesViewT values;
    RCP<const Map> columnMap = inputMatrix->getColMap();
    auto IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    auto ILO = Teuchos::OrdinalTraits<LO>::invalid();
    RCP<const Teuchos::Comm<int> > serialComm = RCP( new Teuchos::SerialComm<int>() );
    for (size_t i=0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            const LO localRow = localRowsSend[j];
            inputMatrix->getLocalRowView(localRow, indices, values);
            const int sz = indices.size();
            for (int k=0; k<sz; k++) {
                const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
                LO targetIndex = targetMap->getLocalElement(colGlobalID);
                if (targetIndex != ILO) count[j]++;
            }
        }
    }
    /*
      std::cout << "local rows and send counts for myPID = " << myPID << std::endl;
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "--- for sending to proc " << procsTo[i] << " ---" << std::endl;
      for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << localRowsSend[j] << " " << count[j] << std::endl;
      }
      }
    */
    size_t sumCount = 0;
    for (size_t i=0; i<count.size(); i++) sumCount += count[i];
    std::vector<LO> localColIDs(sumCount);
    std::vector<SC> matrixValues(sumCount);
    sumCount = 0;
    for (size_t i=0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            const LO localRow = localRowsSend[j];
            inputMatrix->getLocalRowView(localRow, indices, values);
            const int sz = indices.size();
            for (int k=0; k<sz; k++) {
                const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
                LO targetIndex = targetMap->getLocalElement(colGlobalID);
                if (targetIndex != ILO) {
                    localColIDs[sumCount] = targetIndex;
                    matrixValues[sumCount++] = values[k];
                }
            }
        }
    }
    const size_t numRowsRecv = localRowsRecvBegin[numRecvs];
    std::vector<size_t> recvSize(numRecvs), recvRowNumNonzeros(numRowsRecv);
    for (size_t i=0; i<numRecvs; i++)  {
        recvSize[i] = localRowsRecvBegin[i+1] - localRowsRecvBegin[i];
    }
    std::vector<size_t> sendSize(numSends);
    for (size_t i=0; i<numSends; i++) {
        sendSize[i] = localRowsSendBegin[i+1] - localRowsSendBegin[i];
    }
    distributor->doPostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                 Teuchos::ArrayView<const size_t>(sendSize),
                                 Kokkos::View<size_t*, Kokkos::HostSpace>(recvRowNumNonzeros.data(), recvRowNumNonzeros.size()),
                                 Teuchos::ArrayView<const size_t>(recvSize));
    /*
      const int myPID = rowMap->getComm()->getRank();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "row globalIDs and counts for proc " << myPID << " received from proc "
      << procsFrom[i] << std::endl;
      for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << rowMap->getGlobalElement(localRowsRecv[j]) << " "
      << recvRowNumNonzeros[j] << std::endl;
      }
      }
    */ 
    std::vector<size_t> sourceSize(numSends);
    for (size_t i=0; i<numSends; i++) {
        size_t procCount(0);
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            procCount += count[j];
        }
        sourceSize[i] = procCount;
    }
    std::vector<size_t> targetSize(numRecvs);
    size_t numTerms(0);
    for (size_t i=0; i<numRecvs; i++) {
        size_t procCount(0);
        for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
            procCount += recvRowNumNonzeros[j];
        }
        targetSize[i] = procCount;
        numTerms += procCount;
    }
    std::vector<LO> columnsRecv(numTerms);
    std::vector<SC> valuesRecv(numTerms);
    distributor->doPostsAndWaits(Kokkos::View<const LO*, Kokkos::HostSpace>(localColIDs.data(), localColIDs.size()),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<LO*, Kokkos::HostSpace>(columnsRecv.data(), columnsRecv.size()),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    
    using KSX = typename Kokkos::ArithTraits<SC>::val_type;
    const KSX* matrixValues_K = reinterpret_cast<const KSX*>(matrixValues.data());
    KSX* valuesRecv_K = reinterpret_cast<KSX*>(valuesRecv.data());
    const size_t sizeSend = matrixValues.size();
    const size_t sizeRecv = valuesRecv.size();
    distributor->doPostsAndWaits(Kokkos::View<const KSX*, Kokkos::HostSpace>(matrixValues_K, sizeSend),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<KSX*, Kokkos::HostSpace>(valuesRecv_K, sizeRecv),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    const size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> rowCount(numRows, 0);
    // rowCounts for owned rows
    RCP<const Map> colMapSource = inputMatrix->getColMap();
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        TEUCHOS_TEST_FOR_EXCEPTION(localRowSource == ILO, std::runtime_error,"globalID not found in rowMap1to1");
        TEUCHOS_TEST_FOR_EXCEPTION(localRowTarget == ILO, std::runtime_error,"globalID not found in rowMap");
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) rowCount[localRowTarget]++;
        }
    }
    // rowCounts for received rows
    for (LO i=0; i<localRowsRecvBegin[numRecvs]; i++) {
        const LO localRow = localRowsRecv[i];
        rowCount[localRow] = recvRowNumNonzeros[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "rowGIDs and numEntries for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": " << rowCount[i] << std::endl;
      }
    */
    outputMatrix = 
        RCP( new CrsMatrix(rowMap, rowMap, Teuchos::ArrayView<const size_t>(rowCount)) );
    std::vector<LO> indicesVec(numRows);
    std::vector<SC> valuesVec(numRows);
    // insert owned constributions
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        indicesVec.resize(0);
        valuesVec.resize(0);
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) {
                indicesVec.push_back(localColTarget);
                valuesVec.push_back(values[j]);
            }
        }
        outputMatrix->insertLocalValues(localRowTarget, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    // insert received contributions
    numTerms = 0;
    for (LO j=0; j<localRowsRecvBegin[numRecvs]; j++) {
        const LO localRow = localRowsRecv[j];
        indicesVec.resize(recvRowNumNonzeros[j]);
        valuesVec.resize(recvRowNumNonzeros[j]);
        for (size_t k=0; k<recvRowNumNonzeros[j]; k++) {
            indicesVec[k] = columnsRecv[numTerms];
            valuesVec[k] = valuesRecv[numTerms++];
        }
        outputMatrix->insertLocalValues(localRow, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    outputMatrix->fillComplete(rowMap1to1, rowMap1to1);
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "matrix entries for proc (global indices) = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": ";
      outputMatrix->getLocalRowView(i, indices, values);
      for (int j=0; j<indices.size(); j++) {
      std::cout << rowMap->getGlobalElement(indices[j]) << "("
      << values[j] << ") ";
      }
      std::cout << std::endl;
      }
    */
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateRowMap(RCP<const Map> rowMap, 
                  RCP<Tpetra::Distributor> distributor, 
                  std::vector<GO> & rowMapGIDs, 
                  std::vector<LO> & rowMapGIDsBegin)
{
    size_t numSends = distributor->getNumSends();
    size_t numRecvs = distributor->getNumReceives();
    size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> targetValues(numSends);
    std::vector<size_t> targetSizes(numSends, 1);
    std::vector<size_t> sourceValues(numRecvs, numRows);
    std::vector<size_t> sourceSizes(numRecvs, 1); 

    distributor->doReversePostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(sourceValues.data(), sourceValues.size()),
                                        1,
                                        Kokkos::View<size_t*, Kokkos::HostSpace>(targetValues.data(), targetValues.size()));
    int numTerms(0);
    for (size_t i=0; i<targetValues.size(); i++) numTerms += targetValues[i];
    rowMapGIDs.resize(numTerms);
    std::vector<GO> globalIDsSource(numRecvs*numRows);
    numTerms = 0;
    for (size_t i=0; i<numRecvs; i++) {
        for (size_t j=0; j<numRows; j++) {
            globalIDsSource[numTerms++] = rowMap->getGlobalElement(j);
        }
    }
    distributor->doReversePostsAndWaits(Kokkos::View<const GO*, Kokkos::HostSpace>(globalIDsSource.data(), globalIDsSource.size()),
                                        Teuchos::ArrayView<const size_t>(sourceValues),
                                        Kokkos::View<GO*, Kokkos::HostSpace>(rowMapGIDs.data(), rowMapGIDs.size()),
                                        Teuchos::ArrayView<const size_t>(targetValues));
    rowMapGIDsBegin.resize(numSends+1, 0);
    for (size_t i=0; i<numSends; i++) {
        rowMapGIDsBegin[i+1] = rowMapGIDsBegin[i] + targetValues[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "map globalIDs for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numSends; i++) {
      for (LO j=rowMapGIDsBegin[i]; j<rowMapGIDsBegin[i+1]; j++) {
      std::cout << rowMapGIDs[j] << " ";
      }
      std::cout << std::endl;
      }
    */
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
constructDistributor(RCP<const CrsMatrix> inputMatrix, 
                     RCP<const Map> rowMap, 
                     RCP<Tpetra::Distributor> & distributor,
                     std::vector<GO> & ownedRowGIDs,
                     std::vector<LO> & localRowsSend,
                     std::vector<LO> & localRowsSendBegin,
                     std::vector<LO> & localRowsRecv,
                     std::vector<LO> & localRowsRecvBegin)
{
    const LO numRows = rowMap->getLocalNumElements();
    std::vector<GO> globalIDs(numRows);
    for (LO i=0; i<numRows; i++) globalIDs[i] = rowMap->getGlobalElement(i);
    std::vector<int> remotePIDs(numRows);
    std::vector<LO> remoteLocalRows(numRows);
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    rowMap1to1->getRemoteIndexList(Teuchos::ArrayView<GO>(globalIDs),
                                   Teuchos::ArrayView<int>(remotePIDs),
                                   Teuchos::ArrayView<LO>(remoteLocalRows));
    std::vector<int> offProcessorPIDs(numRows);
    const int myPID = rowMap->getComm()->getRank();
    size_t numOffProcessorRows(0), numOnProcessorRows(0);
    ownedRowGIDs.resize(numRows);
    for (LO i=0; i<numRows; i++) {
        if (remotePIDs[i] != myPID) {
            globalIDs[numOffProcessorRows] = globalIDs[i];
            remotePIDs[numOffProcessorRows] = remotePIDs[i];
            remoteLocalRows[numOffProcessorRows++] = remoteLocalRows[i];
        }
        else {
            ownedRowGIDs[numOnProcessorRows++] = globalIDs[i];
        }
    }
    remotePIDs.resize(numOffProcessorRows);
    remoteLocalRows.resize(numOffProcessorRows);
    ownedRowGIDs.resize(numOnProcessorRows);
    std::vector<int> recvPIDs;
    getUniqueEntries(remotePIDs, recvPIDs);
    size_t numRecvs = recvPIDs.size();
    std::map<int,int> offProcessorMap;
    for (size_t i=0; i<numRecvs; i++) {
        offProcessorMap.emplace(recvPIDs[i], i);
    }
    std::vector<GO> recvGIDs(numRecvs);
    std::vector<size_t> count(numRecvs, 0);
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        recvGIDs[iter->second] = globalIDs[i];
        count[iter->second]++;
    }
    localRowsRecvBegin.resize(numRecvs+1, 0);
    for (size_t i=0; i<numRecvs; i++) {
        localRowsRecvBegin[i+1] = localRowsRecvBegin[i] + count[i];
        count[i] = 0;
    }
    localRowsRecv.resize(numOffProcessorRows);
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        const int index = localRowsRecvBegin[iter->second] + count[iter->second];
        localRowsRecv[index] = remoteLocalRows[i];
        count[iter->second]++;
    }
    /*
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "proc " << myPID << " local rows to be received from proc " << recvPIDs[i]
      << std::endl;
      for (int j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << localRowsRecv[j] << " ";
      }
      std::cout << std::endl;
      }
    */
    Teuchos::Array<GO> sendGIDs;
    Teuchos::Array<int> sendPIDs;
    distributor = RCP( new Tpetra::Distributor(rowMap->getComm()) );
    distributor->createFromRecvs(Teuchos::ArrayView<const GO>(recvGIDs),
                                 Teuchos::ArrayView<const int>(recvPIDs),
                                 sendGIDs, sendPIDs);
    TEUCHOS_TEST_FOR_EXCEPTION(distributor->hasSelfMessage() == true,std::runtime_error,
                               "distributor hasSelfMessage error");
    TEUCHOS_TEST_FOR_EXCEPTION(distributor->getNumReceives() != numRecvs,std::runtime_error,
                               "inconsistent numRecvs");
    /*
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      Teuchos::ArrayView<const size_t> lengthsFrom = distributor->getLengthsFrom();
      Teuchos::ArrayView<const size_t> lengthsTo = distributor->getLengthsTo();
      std::cout << "procsTo for myPID = " << myPID << ": ";
      for (int i=0; i<procsTo.size(); i++) std::cout << procsTo[i] << " ";
      std::cout << std::endl;
      std::cout << "procsFrom for myPID = " << myPID << ": ";
      for (int i=0; i<procsFrom.size(); i++) std::cout << procsFrom[i] << " ";
      std::cout << std::endl;
      std::cout << "lengthsFrom for myPID = " << myPID << ": ";
      for (int i=0; i<lengthsFrom.size(); i++) std::cout << lengthsFrom[i] << " ";
      std::cout << std::endl;
      std::cout << "lengthsTo for myPID = " << myPID << ": ";
      for (int i=0; i<lengthsTo.size(); i++) std::cout << lengthsTo[i] << " ";
      std::cout << std::endl;
    */
    size_t numSends = distributor->getNumSends();
    std::vector<size_t> targetValues(numSends);
    std::vector<size_t> targetSizes(numSends, 1);
    std::vector<size_t> sourceSizes(numRecvs, 1); 
    distributor->doReversePostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                        1,
                                        Kokkos::View<size_t*, Kokkos::HostSpace>(targetValues.data(), targetValues.size()));
    localRowsSendBegin.resize(numSends+1, 0);
    for (size_t i=0; i<numSends; i++) {
        localRowsSendBegin[i+1] = localRowsSendBegin[i] + targetValues[i];
    }
    int numTerms = localRowsSendBegin[numSends];
    localRowsSend.resize(numTerms);
    distributor->doReversePostsAndWaits(Kokkos::View<const int*, Kokkos::HostSpace>(localRowsRecv.data(), localRowsRecv.size()),
                                        Teuchos::ArrayView<const size_t>(count),
                                        Kokkos::View<int*, Kokkos::HostSpace>(localRowsSend.data(), localRowsSend.size()),
                                        Teuchos::ArrayView<const size_t>(targetValues));
    /*
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "proc " << myPID << " globalIDs for rows of matrix to be sent to proc " 
      << procsTo[i] << std::endl;
      for (int j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << rowMap1to1->getGlobalElement(localRowsSend[j]) << " ";
      }
      std::cout << std::endl;
      }
    */
    // switch localRowsRecv to on-processor rather than off-processor localRows
    for (size_t i=0; i<numRecvs; i++) count[i] = 0;
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        const int index = localRowsRecvBegin[iter->second] + count[iter->second];
        localRowsRecv[index] = rowMap->getLocalElement(globalIDs[i]);
        count[iter->second]++;
    }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getUniqueEntries(const std::vector<int> & vector, 
                 std::vector<int> & vectorUnique)
{
    vectorUnique = vector;
    std::sort(vectorUnique.begin(), vectorUnique.end());
    auto iter = std::unique(vectorUnique.begin(), vectorUnique.end());
    vectorUnique.erase(iter, vectorUnique.end());
}


