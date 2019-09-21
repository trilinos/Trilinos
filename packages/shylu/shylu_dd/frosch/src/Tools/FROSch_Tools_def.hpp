//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_TOOLS_DEF_HPP
#define _FROSCH_TOOLS_DEF_HPP

#include <FROSch_Tools_decl.hpp>

#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <typename LO,typename GO>
    OverlappingData<LO,GO>::OverlappingData(GO gid,
                                            int pid,
                                            LO lid) :
    GID_ (gid),
    PIDs_ (1,pid),
    LIDs_ (1,lid)
    {

    }

    template <typename LO,typename GO>
    int OverlappingData<LO,GO>::Merge(const RCP<OverlappingData<LO,GO> > od) const
    {
        FROSCH_ASSERT(GID_ == od->GID_,"FROSch::OverlappingData : ERROR: GID_ != od->GID_");
        for (typename IntVec::iterator it = od->PIDs_.begin(); it != od->PIDs_.end(); it++) {
            PIDs_.push_back(*it);
        }
        for (typename LOVec::iterator it = od->LIDs_.begin(); it != od->LIDs_.end(); it++) {
            LIDs_.push_back(*it);
        }
        return 0;
    }

    template <typename LO,typename GO>
    int MergeList(Array<RCP<OverlappingData<LO,GO> > > &odList)
    {
        int numPackages = 0;
        std::sort(odList.begin(),
                  odList.end(),
                  [] (const RCP<OverlappingData<LO,GO> >& lhs, const RCP<OverlappingData<LO,GO> >& rhs) {
                      return lhs->GID_ < rhs->GID_;
                  }
                  );
        odList.erase(std::unique(odList.begin(),
                                 odList.end(),
                                 [] (const RCP<OverlappingData<LO,GO> > lhs, const RCP<OverlappingData<LO,GO> > rhs) {
                                     if (lhs->GID_ == rhs->GID_) {
                                         lhs->Merge(rhs);
                                         return true;
                                     } else {
                                         return false;
                                     }
                                 }
                                 ),odList.end());

        typename Array<RCP<OverlappingData<LO,GO> > >::iterator tmpODPtrVecIt;
        for (tmpODPtrVecIt=odList.begin(); tmpODPtrVecIt!=odList.end(); tmpODPtrVecIt++) {
            numPackages += (*tmpODPtrVecIt)->PIDs_.size()*(*tmpODPtrVecIt)->PIDs_.size();
        }
        return numPackages;
    }

    template <typename LO,typename GO,typename NO>
    LowerPIDTieBreak<LO,GO,NO>::LowerPIDTieBreak(CommPtr comm,
                                                 ConstXMapPtr originalMap,
                                                 UN levelID) :
    MpiComm_ (comm),
    OriginalMap_ (originalMap),
    ElementCounter_ (),
    OverlappingDataList_ (2*(OriginalMap_->getGlobalNumElements()/MpiComm_->getSize())),
    ComponentsSubdomains_ (OriginalMap_->getNodeNumElements()),
    LevelID_ (levelID)
    {
        FROSCH_TIMER_START_LEVELID(lowerPIDTieBreakTime,"LowerPIDTieBreak::LowerPIDTieBreak");
    }

    template <typename LO,typename GO,typename NO>
    int LowerPIDTieBreak<LO,GO,NO>::sendDataToOriginalMap()
    {
        FROSCH_TIMER_START_LEVELID(sendDataToOriginalMapTime,"LowerPIDTieBreak::sendDataToOriginalMap");
        // This is done analogously to  DistributedNoncontiguousDirectory<LO, GO, NT>::DistributedNoncontiguousDirectory(const map_type& map,const tie_break_type& tie_break)
        typedef typename ArrayView<const GO>::size_type size_type;

        OverlappingDataList_.resize(ElementCounter_);

        int numPackages = MergeList(OverlappingDataList_);

        Tpetra::Distributor distor (MpiComm_);

        const int packetSize = 2;
        IntVec sendImageIDs(numPackages);
        GOVec exportEntries(packetSize*numPackages); // data to send out
        {
            typename OverlappingDataPtrVec::iterator tmpODPtrVecIt;
            typename IntVec::iterator tmpIntVecIt;
            typename IntVec::iterator tmpIntVecIt2;
            typename LOVec::iterator tmpLOVecIt;

            size_type exportIndex = 0;
            size_type exportIndex2 = 0;
            for (tmpODPtrVecIt = OverlappingDataList_.begin(); tmpODPtrVecIt!=OverlappingDataList_.end(); tmpODPtrVecIt++) {
                tmpLOVecIt = (*tmpODPtrVecIt)->LIDs_.begin();
                for (tmpIntVecIt = (*tmpODPtrVecIt)->PIDs_.begin(); tmpIntVecIt != (*tmpODPtrVecIt)->PIDs_.end(); tmpIntVecIt++) {
                    for (tmpIntVecIt2 = (*tmpODPtrVecIt)->PIDs_.begin(); tmpIntVecIt2 != (*tmpODPtrVecIt)->PIDs_.end(); tmpIntVecIt2++) {
                        exportEntries[exportIndex++] = (*tmpODPtrVecIt)->GID_;
                        exportEntries[exportIndex++] = as<GO>(*tmpIntVecIt2);
                        sendImageIDs[exportIndex2++] = *tmpIntVecIt;
                    }
                    tmpLOVecIt++;
                }
            }
        }
        distor.createFromSends(sendImageIDs);

        GOVec importElements(packetSize*distor.getTotalReceiveLength());

        distor.doPostsAndWaits(exportEntries().getConst(),packetSize,importElements());

        LO length = importElements.size()/packetSize;
        for (LO i=0; i<length; i++) {
            ComponentsSubdomains_[OriginalMap_->getLocalElement(importElements[2*i])].push_back(importElements[2*i+1]);
        }

        return 0;
    }

    template <typename LO,typename GO,typename NO>
    std::size_t LowerPIDTieBreak<LO,GO,NO>::selectedIndex(GO GID,
                                                          const std::vector<std::pair<int,LO> > & pid_and_lid) const
    {
        // Always choose index of pair with smallest pid
        const std::size_t numLids = pid_and_lid.size();
        std::size_t idx = 0;
        int minpid = pid_and_lid[0].first;
        std::size_t minidx = 0;

        for (idx = 0; idx < numLids; ++idx) {
            // Organize the overlapping data to send it back to the original processes
            OverlappingDataList_[ElementCounter_].reset(new OverlappingData<LO,GO>(GID,
                                                                                   pid_and_lid[idx].first,
                                                                                   pid_and_lid[idx].second));
            ElementCounter_++;
            if (pid_and_lid[idx].first < minpid) {
                minpid = pid_and_lid[idx].first;
                minidx = idx;
            }
        }
        return minidx;
    }

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildUniqueMap(const RCP<const Map<LO,GO,NO> > map,
                                             bool useCreateOneToOneMap,
                                             RCP<Tpetra::Details::TieBreak<LO,GO> > tieBreak)
    {
        FROSCH_TIMER_START(buildUniqueMapTime,"BuildUniqueMap");
        if (useCreateOneToOneMap && map->lib()==UseTpetra) {
            // Obtain the underlying Tpetra Map
            const RCP<const TpetraMap<LO,GO,NO> >& xTpetraMap = rcp_dynamic_cast<const TpetraMap<LO,GO,NO> >(map);
            RCP<const Tpetra::Map<LO,GO,NO> > tpetraMap = xTpetraMap->getTpetra_Map();

            RCP<const Tpetra::Map<LO,GO,NO> > tpetraMapUnique;
            if (tieBreak.is_null()) {
                tpetraMapUnique = createOneToOne(tpetraMap);
            } else {
                tpetraMapUnique = createOneToOne(tpetraMap,*tieBreak);
            }

            RCP<const TpetraMap<LO,GO,NO> > xTpetraMapUnique(new const TpetraMap<LO,GO,NO>(tpetraMapUnique));
            return rcp_dynamic_cast<const Map<LO,GO,NO> >(xTpetraMapUnique);
        } else { // This is an alternative implementation to createOneToOneMap()
            if (map->lib()==UseEpetra && map->getComm()->getRank()==0) std::cout << "FROSch::BuildUniqueMap : WARNING: createOneToOneMap() does not exist for Epetra => Using a different implementation" << std::endl;

            RCP<Vector<GO,LO,GO,NO> > myIndices = VectorFactory<GO,LO,GO,NO>::Build(map);
            myIndices->putScalar(map->getComm()->getRank()+1);

            RCP<Map<LO,GO,NO> > linearMap = MapFactory<LO,GO,NO>::Build(map->lib(),map->getMaxAllGlobalIndex()+1,0,map->getComm());
            RCP<Vector<GO,LO,GO,NO> > globalIndices = VectorFactory<GO,LO,GO,NO>::Build(linearMap);

            RCP<Import<LO,GO,NO> > importer = ImportFactory<LO,GO,NO>::Build(map,linearMap);
            RCP<Import<LO,GO,NO> > importer2 = ImportFactory<LO,GO,NO>::Build(linearMap,map); // AH 10/16/2017: Ist der notwendig??? Mit Epetra ging es auch ohne einen zweiten XImport und stattdessen mit einem Export
            globalIndices->doImport(*myIndices,*importer,INSERT);
            myIndices->putScalar(0);
            myIndices->doImport(*globalIndices,*importer2,ADD);

            Array<GO> uniqueVector;
            for (unsigned i=0; i<myIndices->getLocalLength(); i++) {
                if (myIndices->getData(0)[i] == map->getComm()->getRank()+1) {
                    uniqueVector.push_back(map->getGlobalElement(i));
                }
            }
            return MapFactory<LO,GO,NO>::Build(map->lib(),-1,uniqueVector(),0,map->getComm()); // We need this setup for maps with offset (with MaxGID+1 not being the number of global elements), or we need an allreduce to determine the number of global elements from uniqueVector
            //        return MapFactory<LO,GO,NO>::Build(map->lib(),map->getMaxAllGlobalIndex()+1,uniqueVector(),0,map->getComm());
        }
    }

    template <class SC,class LO,class GO,class NO>
    ArrayRCP<RCP<const Map<LO,GO,NO> > > BuildRepeatedSubMaps(RCP<const Matrix<SC,LO,GO,NO> > matrix,
                                                              ArrayRCP<RCP<const Map<LO,GO,NO> > > subMaps)
    {
        FROSCH_TIMER_START(buildRepeatedSubMapsTime,"BuildRepeatedSubMaps");
        ArrayRCP<RCP<Map<LO,GO,NO> > > repeatedSubMaps(subMaps.size());
        for (unsigned i = 0; i < subMaps.size(); i++) {
            RCP<Matrix<SC,LO,GO,NO> > subMatrixII;
            ArrayView<GO> indI = av_const_cast<GO> ( subMaps[i]->getNodeElementList() );

            BuildSubmatrix(matrix,indI,subMatrixII);

            repeatedSubMaps[i] = BuildRepeatedMap(subMatrixII);
        }

        return repeatedSubMaps;
    }

    template <class LO,class GO,class NO>
    ArrayRCP<RCP<const Map<LO,GO,NO> > > BuildRepeatedSubMaps(RCP<const CrsGraph<LO,GO,NO> > graph,
                                                              ArrayRCP<RCP<const Map<LO,GO,NO> > > subMaps)
    {
        FROSCH_TIMER_START(buildRepeatedSubMapsTime,"BuildRepeatedSubMaps");
        ArrayRCP<RCP<Map<LO,GO,NO> > > repeatedSubMaps(subMaps.size());
        for (unsigned i = 0; i < subMaps.size(); i++) {
            RCP<CrsGraph<LO,GO,NO> > subGraphII;
            ArrayView<GO> indI = av_const_cast<GO> ( subMaps[i]->getNodeElementList() );

            BuildSubgraph(graph,indI,subGraphII);

            repeatedSubMaps[i] = BuildRepeatedMap(subGraphII);
        }

        return repeatedSubMaps;
    }

    template <class SC,class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConst(RCP<const Matrix<SC,LO,GO,NO> > matrix)
    {
        FROSCH_TIMER_START(buildRepeatedMapNonConstTime,"BuildRepeatedMapNonConst");
        RCP<Map<LO,GO,NO> > uniqueMap = MapFactory<LO,GO,NO>::Build(matrix->getRowMap(),1);
        RCP<const Map<LO,GO,NO> > overlappingMap = uniqueMap.getConst();
        ExtendOverlapByOneLayer<SC,LO,GO,NO>(matrix,overlappingMap,matrix,overlappingMap);

        RCP<Matrix<SC,LO,GO,NO> > tmpMatrix = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,matrix->getGlobalMaxNumRowEntries());

        RCP<Import<LO,GO,NO> > scatter;
        RCP<Export<LO,GO,NO> > gather = ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);

        if (tmpMatrix->getRowMap()->lib()==UseEpetra) {
            scatter = ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);
            tmpMatrix->doImport(*matrix,*scatter,ADD);
        } else {
            tmpMatrix->doImport(*matrix,*gather,ADD);
        }

        Array<SC> one(1,ScalarTraits<SC>::one());
        Array<GO> myPID(1,uniqueMap->getComm()->getRank());

        RCP<Matrix<SC,LO,GO,NO> > commMat = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        RCP<Matrix<SC,LO,GO,NO> > commMatTmp = MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            if (uniqueMap->getLocalElement(globalRow)<0) {
                ArrayView<const GO> indices;
                ArrayView<const SC> values;
                tmpMatrix->getGlobalRowView(globalRow,indices,values);

                LO j=0;
                while (j<indices.size() && overlappingMap->getLocalElement(indices[j])>=0) {
                    j++;
                }
                if (j!=indices.size()) {
                    commMat->insertGlobalValues(overlappingMap->getGlobalElement(i),myPID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                }
            }
        }
        commMat->fillComplete();
        commMatTmp->doExport(*commMat,*gather,INSERT);

        RCP<Matrix<SC,LO,GO,NO> > commMatTmp2 = MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            GO globalRow = uniqueMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            commMatTmp->getGlobalRowView(globalRow,indices,values);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    Array<GO> pID(1,indices[j]);
                    if (pID<myPID) {
                        commMatTmp2->insertGlobalValues(globalRow,pID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                    }
                }
            }
        }
        commMatTmp2->fillComplete();
        commMatTmp.reset();
        commMat = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);

        commMat->doImport(*commMatTmp2,*gather,ADD);

        ArrayView<const GO> myGlobalElements = uniqueMap->getNodeElementList();
        Array<GO> repeatedIndices(uniqueMap->getNodeNumElements());
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            repeatedIndices.at(i) = myGlobalElements[i];
        }

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            commMat->getGlobalRowView(globalRow,indices,values);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    GO pID = indices[j];
                    if (pID==myPID[0]) {
                        repeatedIndices.push_back(globalRow);
                    }
                }
            }
        }
        sortunique(repeatedIndices);
        return MapFactory<LO,GO,NO>::Build(tmpMatrix->getRowMap()->lib(),-1,repeatedIndices(),0,matrix->getRowMap()->getComm());
    }

    template <class SC,class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMap(RCP<const Matrix<SC,LO,GO,NO> > matrix)
    {
        return BuildRepeatedMapNonConst(matrix).getConst();
    }

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConst(RCP<const CrsGraph<LO,GO,NO> > graph)
    {
        FROSCH_TIMER_START(buildRepeatedMapNonConstTime,"BuildRepeatedMapNonConst");
        RCP<Map<LO,GO,NO> > uniqueMap = MapFactory<LO,GO,NO>::Build(graph->getRowMap(),1);
        RCP<const Map<LO,GO,NO> > overlappingMap = uniqueMap.getConst();
        ExtendOverlapByOneLayer<LO,GO,NO>(graph,overlappingMap,graph,overlappingMap);

        RCP<CrsGraph<LO,GO,NO> > tmpGraph = CrsGraphFactory<LO,GO,NO>::Build(overlappingMap,graph->getGlobalMaxNumRowEntries());

        RCP<Import<LO,GO,NO> > scatter;
        RCP<Export<LO,GO,NO> > gather = ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);

        if (tmpGraph->getRowMap()->lib()==UseEpetra) {
            scatter = ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);
            tmpGraph->doImport(*graph,*scatter,ADD);
        } else {
            tmpGraph->doImport(*graph,*gather,ADD);
        }

        Array<GO> myPID(1,uniqueMap->getComm()->getRank());

        RCP<CrsGraph<LO,GO,NO> > commGraph = CrsGraphFactory<LO,GO,NO>::Build(overlappingMap,10);
        RCP<CrsGraph<LO,GO,NO> > commGraphTmp = CrsGraphFactory<LO,GO,NO>::Build(uniqueMap,10);

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            if (uniqueMap->getLocalElement(globalRow)<0) {
                ArrayView<const GO> indices;
                tmpGraph->getGlobalRowView(globalRow,indices);

                LO j=0;
                while (j<indices.size() && overlappingMap->getLocalElement(indices[j])>=0) {
                    j++;
                }
                if (j!=indices.size()) {
                    commGraph->insertGlobalIndices(overlappingMap->getGlobalElement(i),myPID());
                }
            }
        }
        commGraph->fillComplete();
        commGraphTmp->doExport(*commGraph,*gather,INSERT);

        RCP<CrsGraph<LO,GO,NO> > commGraphTmp2 = CrsGraphFactory<LO,GO,NO>::Build(uniqueMap,10);
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            GO globalRow = uniqueMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            commGraphTmp->getGlobalRowView(globalRow,indices);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    Array<GO> pID(1,indices[j]);
                    if (pID<myPID) {
                        commGraphTmp2->insertGlobalIndices(globalRow,pID());
                    }
                }
            }
        }
        commGraphTmp2->fillComplete();
        commGraphTmp.reset();
        commGraph = CrsGraphFactory<LO,GO,NO>::Build(overlappingMap,10);

        commGraph->doImport(*commGraphTmp2,*gather,ADD);

        ArrayView<const GO> myGlobalElements = uniqueMap->getNodeElementList();
        Array<GO> repeatedIndices(uniqueMap->getNodeNumElements());
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            repeatedIndices.at(i) = myGlobalElements[i];
        }

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            commGraph->getGlobalRowView(globalRow,indices);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    GO pID = indices[j];
                    if (pID==myPID[0]) {
                        repeatedIndices.push_back(globalRow);
                    }
                }
            }
        }
        sortunique(repeatedIndices);
        return MapFactory<LO,GO,NO>::Build(tmpGraph->getRowMap()->lib(),-1,repeatedIndices(),0,graph->getRowMap()->getComm());
    }

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMap(RCP<const CrsGraph<LO,GO,NO> > graph)
    {
        return BuildRepeatedMapNonConst(graph).getConst();
    }

    /*
    template <class SC,class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMap(RCP<Matrix<SC,LO,GO,NO> > matrix)
    {
        RCP<Map<LO,GO,NO> > uniqueMap = MapFactory<LO,GO,NO>::Build(matrix->getRowMap(),1);
        RCP<Map<LO,GO,NO> > overlappingMap = uniqueMap;
        ExtendOverlapByOneLayer<SC,LO,GO,NO>(matrix,overlappingMap);

        RCP<Matrix<SC,LO,GO,NO> > tmpMatrix = matrix;
        matrix = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,2*tmpMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);

        matrix->doImport(*tmpMatrix,*scatter,ADD);

        Array<SC> one(1,ScalarTraits<SC>::one());
        Array<GO> myPID(1,uniqueMap->getComm()->getRank());

        RCP<Matrix<SC,LO,GO,NO> > commMat = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        RCP<Matrix<SC,LO,GO,NO> > commMatTmp = MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        RCP<Export<LO,GO,NO> > commExporter = ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            if (uniqueMap->getLocalElement(globalRow)<0) {
                ArrayView<const GO> indices;
                ArrayView<const SC> values;
                matrix->getGlobalRowView(globalRow,indices,values);

                LO j=0;
                while (j<indices.size() && overlappingMap->getLocalElement(indices[j])>=0) {
                    j++;
                }
                if (j!=indices.size()) {
                    commMat->insertGlobalValues(overlappingMap->getGlobalElement(i),myPID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                }
            }
        }
        commMat->fillComplete();
        commMatTmp->doExport(*commMat,*commExporter,INSERT);

        RCP<Matrix<SC,LO,GO,NO> > commMatTmp2 = MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            GO globalRow = uniqueMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            commMatTmp->getGlobalRowView(globalRow,indices,values);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    Array<GO> pID(1,indices[j]);
                    if (pID<myPID) {
                        commMatTmp2->insertGlobalValues(globalRow,pID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                    }
                }
            }
        }
        commMatTmp2->fillComplete();
        commMatTmp.reset();
        commMat = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        commMat->doImport(*commMatTmp2,*commExporter,ADD);

        ArrayView<const GO> myGlobalElements = uniqueMap->getNodeElementList();
        Array<GO> repeatedIndices(uniqueMap->getNodeNumElements());
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            repeatedIndices.at(i) = myGlobalElements[i];
        }

        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            commMat->getGlobalRowView(globalRow,indices,values);

            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    GO pID = indices[j];
                    if (pID==myPID[0]) {
                        repeatedIndices.push_back(globalRow);
                    }
                }
            }
        }
        sortunique(repeatedIndices);
        return MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,repeatedIndices(),0,matrix->getRowMap()->getComm());
    }
    */

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer_Old(RCP<const Matrix<SC,LO,GO,NO> > inputMatrix,
                                    RCP<const Map<LO,GO,NO> > inputMap,
                                    RCP<const Matrix<SC,LO,GO,NO> > &outputMatrix,
                                    RCP<const Map<LO,GO,NO> > &outputMap)
    {
        FROSCH_TIMER_START(extendOverlapByOneLayer_OldTime,"ExtendOverlapByOneLayer_Old");
        RCP<Matrix<SC,LO,GO,NO> > tmpMatrix = MatrixFactory<SC,LO,GO,NO>::Build(inputMap,inputMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(inputMatrix->getRowMap(),inputMap);
        tmpMatrix->doImport(*inputMatrix,*scatter,ADD);
        outputMatrix = tmpMatrix.getConst();

        Array<GO> indicesOverlappingSubdomain(0);
        for (unsigned i=0; i<inputMap->getNodeNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            outputMatrix->getGlobalRowView(inputMap->getGlobalElement(i),indices,values);

            for (LO j=0; j<indices.size(); j++) {
                indicesOverlappingSubdomain.push_back(indices[j]);
            }
        }
        sortunique(indicesOverlappingSubdomain);
        outputMap = MapFactory<LO,GO,NO>::Build(inputMap->lib(),-1,indicesOverlappingSubdomain(),0,inputMap->getComm());
        tmpMatrix->fillComplete(inputMatrix->getDomainMap(),inputMatrix->getRangeMap());

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(RCP<const Matrix<SC,LO,GO,NO> > inputMatrix,
                                RCP<const Map<LO,GO,NO> > inputMap,
                                RCP<const Matrix<SC,LO,GO,NO> > &outputMatrix,
                                RCP<const Map<LO,GO,NO> > &outputMap)
    {
        FROSCH_TIMER_START(extendOverlapByOneLayerTime,"ExtendOverlapByOneLayer");
        RCP<Matrix<SC,LO,GO,NO> > tmpMatrix = MatrixFactory<SC,LO,GO,NO>::Build(inputMap,inputMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(inputMatrix->getRowMap(),inputMap);
        tmpMatrix->doImport(*inputMatrix,*scatter,ADD);
        tmpMatrix->fillComplete(inputMatrix->getDomainMap(),inputMatrix->getRangeMap());

        outputMatrix = tmpMatrix.getConst();
        outputMap = outputMatrix->getColMap();

        return 0;
    }

    template <class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(RCP<const CrsGraph<LO,GO,NO> > inputGraph,
                                RCP<const Map<LO,GO,NO> > inputMap,
                                RCP<const CrsGraph<LO,GO,NO> > &outputGraph,
                                RCP<const Map<LO,GO,NO> > &outputMap)
    {
        FROSCH_TIMER_START(extendOverlapByOneLayerTime,"ExtendOverlapByOneLayer");
        RCP<CrsGraph<LO,GO,NO> > tmpGraph = CrsGraphFactory<LO,GO,NO>::Build(inputMap,inputGraph->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(inputGraph->getRowMap(),inputMap);
        tmpGraph->doImport(*inputGraph,*scatter,ADD);
        tmpGraph->fillComplete(inputGraph->getDomainMap(),inputGraph->getRangeMap());

        outputGraph = tmpGraph.getConst();
        outputMap = outputGraph->getColMap();

        return 0;
    }

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > AssembleMaps(ArrayView<RCP<Map<LO,GO,NO> > > mapVector,
                                     ArrayRCP<ArrayRCP<LO> > &partMappings)
    {
        FROSCH_TIMER_START(assembleMapsTime,"AssembleMaps");
        FROSCH_ASSERT(mapVector.size()>0,"Length of mapVector is == 0!");
        LO i = 0;
        LO localstart = 0;
        LO sizetmp = 0;
        LO size = 0;
        GO globalstart = 0;

        partMappings = ArrayRCP<ArrayRCP<LO> >(mapVector.size());

        ArrayRCP<GO> assembledMapTmp(0);
        for (unsigned j=0; j<mapVector.size(); j++) {
            sizetmp = mapVector[j]->getNodeNumElements();
            partMappings[j] = ArrayRCP<LO>(sizetmp);

            size += sizetmp;
            assembledMapTmp.resize(size);

            localstart = i;
            while (i<localstart+sizetmp) {
                partMappings[j][i-localstart] = i;
                assembledMapTmp[i] = globalstart + mapVector[j]->getGlobalElement(i-localstart);
                i++;
            }
            //std::cout << mapVector[j]->getMaxAllGlobalIndex() << std::endl;
            /*
            globalstart += mapVector[j]->getMaxAllGlobalIndex();

            if (mapVector[0]->lib()==UseEpetra || mapVector[j]->getGlobalNumElements()>0) {
                globalstart += 1;
            }
             */

            globalstart += std::max(mapVector[j]->getMaxAllGlobalIndex(),(GO)-1)+1; // AH 04/05/2018: mapVector[j]->getMaxAllGlobalIndex() can result in -2147483648 if the map is empty on the process => introducing max(,)

            //if (mapVector[j]->getComm()->getRank() == 0) std::cout << std::endl << globalstart << std::endl;
        }
        return MapFactory<LO,GO,NO>::Build(mapVector[0]->lib(),-1,assembledMapTmp(),0,mapVector[0]->getComm());
    }

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > MergeMapsNonConst(ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVector)
    {
        FROSCH_TIMER_START(mergeMapsNonConstTime,"MergeMapsNonConst");
        FROSCH_ASSERT(!mapVector.is_null(),"mapVector is null!");
        FROSCH_ASSERT(mapVector.size()>0,"Length of mapVector is == 0!");

        Array<GO> elementList(mapVector[0]->getNodeElementList());
        GO tmpOffset = 0;
        for (unsigned i=1; i<mapVector.size(); i++) {
            LO nodeNumElements = mapVector[i]->getNodeNumElements();
            tmpOffset += mapVector[i-1]->getMaxAllGlobalIndex()+1;

            Array<GO> subElementList(nodeNumElements);
            for (LO j=0; j<nodeNumElements; j++) {
                subElementList.at(j) = mapVector[i]->getGlobalElement(j)+tmpOffset;
            }

            elementList.insert(elementList.end(),subElementList.begin(),subElementList.end());
        }
        return MapFactory<LO,GO,NO>::Build(mapVector[0]->lib(),-1,elementList(),0,mapVector[0]->getComm());
    }

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > MergeMaps(ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVector)
    {
        return MergeMapsNonConst(mapVector).getConst();
    }

    template <class LO,class GO,class NO>
    int BuildDofMapsVec(const ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVec,
                        ArrayRCP<unsigned> dofsPerNodeVec,
                        ArrayRCP<FROSch::DofOrdering> dofOrderingVec,
                        ArrayRCP<RCP<const Map<LO,GO,NO> > > &nodesMapVec,
                        ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > >&dofMapsVec)
    {
        FROSCH_TIMER_START(buildDofMapsVecTime,"BuildDofMapsVec");
        unsigned numberBlocks = mapVec.size();
        nodesMapVec = ArrayRCP<RCP<const Map<LO,GO,NO> > > (numberBlocks);
        dofMapsVec = ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > > (numberBlocks);

        GO tmpOffset = 0;
        for (unsigned i = 0 ; i < numberBlocks; i++) {
            BuildDofMaps(mapVec[i],dofsPerNodeVec[i],dofOrderingVec[i],nodesMapVec[i],dofMapsVec[i],tmpOffset);
            tmpOffset += mapVec[i]->getMaxAllGlobalIndex()+1;
        }

        return 0;
    }


    template <class LO,class GO,class NO>
    int BuildDofMaps(const RCP<const Map<LO,GO,NO> > map,
                     unsigned dofsPerNode,
                     unsigned dofOrdering,
                     RCP<const Map<LO,GO,NO> > &nodesMap,
                     ArrayRCP<RCP<const Map<LO,GO,NO> > > &dofMaps,
                     GO offset)
    {
        FROSCH_TIMER_START(buildDofMapsTime,"BuildDofMaps");
        //if (map->getComm()->getRank()==0) std::cout << "WARNING: BuildDofMaps is yet to be tested...\n";
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(map->getGlobalNumElements()%dofsPerNode==0 && map->getNodeNumElements()%dofsPerNode==0,"ERROR: The number of dofsPerNode does not divide the number of global dofs in the map!");

        Array<GO> nodes(map->getNodeNumElements()/dofsPerNode);
        Array<ArrayRCP<GO> > dofs(dofsPerNode);
        for (unsigned j=0; j<dofsPerNode; j++) {
            dofs[j] = ArrayRCP<GO>(map->getNodeNumElements()/dofsPerNode);
        }
        if (dofOrdering==0) {
            for (unsigned i=0; i<nodes.size(); i++) {
                nodes[i] = map->getGlobalElement(dofsPerNode*i)/dofsPerNode;
                for (unsigned j=0; j<dofsPerNode; j++) {
                    dofs[j][i] = dofsPerNode*nodes[i]+j+offset;
                }
            }
        } else if (dofOrdering == 1) {
            GO numGlobalIDs = map->getMaxAllGlobalIndex()+1;
            for (unsigned i=0; i<nodes.size(); i++) {
                nodes[i] = map->getGlobalElement(i);
                for (unsigned j=0; j<dofsPerNode; j++) {
                    dofs[j][i] = nodes[i]+j*numGlobalIDs/dofsPerNode+offset;
                }
            }
        } else {
            FROSCH_ASSERT(false,"dofOrdering unknown.");
        }
        nodesMap = MapFactory<LO,GO,NO>::Build(map->lib(),-1,nodes(),0,map->getComm());

        dofMaps = ArrayRCP<RCP<const Map<LO,GO,NO> > >(dofsPerNode);
        for (unsigned j=0; j<dofsPerNode; j++) {
            dofMaps[j] = MapFactory<LO,GO,NO>::Build(map->lib(),-1,dofs[j](),0,map->getComm());
        }
        return 0;
    }

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildMapFromDofMaps(const ArrayRCP<RCP<const Map<LO,GO,NO> > > &dofMaps,
                                                  unsigned dofsPerNode,
                                                  unsigned dofOrdering)
    {
        FROSCH_TIMER_START(buildMapFromDofMapsTime,"BuildMapFromDofMaps");
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(!dofMaps.is_null(),"dofMaps.is_null().");
        FROSCH_ASSERT(dofMaps.size()==dofsPerNode,"dofMaps.size!=dofsPerNode.");
        for (unsigned i=0; i<dofMaps.size(); i++) {
            FROSCH_ASSERT(dofMaps[i]->getGlobalNumElements()%dofsPerNode==0 && dofMaps[i]->getNodeNumElements()%dofsPerNode==0,"ERROR: The number of dofsPerNode does not divide the number of global dofs in the dofMaps!");
        }

        unsigned numNodes = dofMaps[0]->getNodeNumElements();
        Array<GO> globalIDs(numNodes);
        if (dofOrdering==0) {
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<numNodes; j++) {
                    globalIDs[dofsPerNode*j+i] = dofMaps[i]->getGlobalElement(j);
                }
            }
        } else if (dofOrdering == 1) {
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<numNodes; j++) {
                    globalIDs[j+i*numNodes] = dofMaps[i]->getGlobalElement(j);
                }
            }
        } else {
            FROSCH_ASSERT(false,"dofOrdering unknown.");
        }
        return MapFactory<LO,GO,NO>::Build(dofMaps[0]->lib(),-1,globalIDs(),0,dofMaps[0]->getComm());
    }

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildMapFromNodeMap(RCP<Map<LO,GO,NO> > &nodesMap,
                                            unsigned dofsPerNode,
                                            unsigned dofOrdering)
    {
        FROSCH_TIMER_START(buildMapFromNodeMapTime,"BuildMapFromNodeMap");
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(!nodesMap.is_null(),"nodesMap.is_null().");

        unsigned numNodes = nodesMap->getNodeNumElements();
        Array<GO> globalIDs(dofsPerNode*numNodes);
        if (dofOrdering==0) {
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<numNodes; j++) {
                    globalIDs[dofsPerNode*j+i] = dofsPerNode*nodesMap->getGlobalElement(j)+i;
                }
            }
        } else if (dofOrdering == 1) {
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<numNodes; j++) {
                    globalIDs[j+i*numNodes] = nodesMap->getGlobalElement(j)+i*nodesMap->getGlobalNumElements();
                }
            }
        } else {
            FROSCH_ASSERT(false,"dofOrdering unknown.");
        }
        return MapFactory<LO,GO,NO>::Build(nodesMap->lib(),-1,globalIDs(),0,nodesMap->getComm());
    }

//    template <class LO,class GO,class NO>
//    ArrayRCP<RCP<Map<LO,GO,NO> > > BuildSubMaps(RCP<const Map<LO,GO,NO> > &fullMap,
//                                                                          ArrayRCP<GO> maxSubGIDVec){
//
//          RCP<Map<LO,GO,NO> > fullMapNonConst = rcp_dynamic_cast<Map<LO,GO,NO> >(fullMap);
//    }

    template <class LO,class GO,class NO>
    ArrayRCP<RCP<Map<LO,GO,NO> > > BuildSubMaps(RCP<const Map<LO,GO,NO> > &fullMap,
                                                ArrayRCP<GO> maxSubGIDVec)
    {
        FROSCH_TIMER_START(buildSubMapsTime,"BuildSubMaps");
        ArrayRCP<RCP<Map<LO,GO,NO> > > subMaps(maxSubGIDVec.size());

        Array<Array<GO> > indicesSubMaps(maxSubGIDVec.size());
        ArrayView<const GO> nodeElementList = fullMap->getNodeElementList();
        for (unsigned i = 0; i<fullMap->getNodeNumElements(); i++) {
            LO subMapNumber = -1;
            for (unsigned j = (maxSubGIDVec.size()); j > 0; j--) {
                if (nodeElementList[i] <= maxSubGIDVec[j-1]) {
                    subMapNumber = j-1;
                }
            }
            FROSCH_ASSERT(subMapNumber>-1,"Index not resolved on sub map.");
            indicesSubMaps[subMapNumber].push_back(nodeElementList[i]);
        }
        for (unsigned j = 0 ; j < maxSubGIDVec.size(); j++) {
            subMaps[j] = MapFactory<LO,GO,NO>::Build(fullMap->lib(),-1,indicesSubMaps[j](),0,fullMap->getComm());
        }
        return subMaps;
    }

    template <class SC,class LO,class GO,class NO>
    ArrayRCP<GO> FindOneEntryOnlyRowsGlobal(RCP<const Matrix<SC,LO,GO,NO> > matrix,
                                            RCP<const Map<LO,GO,NO> > repeatedMap)
    {
        FROSCH_TIMER_START(findOneEntryOnlyRowsGlobalTime,"FindOneEntryOnlyRowsGlobal");
        RCP<Matrix<SC,LO,GO,NO> > repeatedMatrix = MatrixFactory<SC,LO,GO,NO>::Build(repeatedMap,2*matrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(matrix->getRowMap(),repeatedMap);
        repeatedMatrix->doImport(*matrix,*scatter,ADD);

        ArrayRCP<GO> oneEntryOnlyRows(repeatedMatrix->getNodeNumRows());
        LO tmp = 0;
        LO nnz;
        GO row;
        for (unsigned i=0; i<repeatedMatrix->getNodeNumRows(); i++) {
            row = repeatedMap->getGlobalElement(i);
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            repeatedMatrix->getGlobalRowView(row,indices,values);
            nnz = indices.size();
            if (indices.size()==1) {
                oneEntryOnlyRows[tmp] = row;
                tmp++;
            } else {
                for (LO j=0; j<values.size(); j++) {
                    if (fabs(values[j])<1.0e-12) {
                        nnz--;
                    }
                }
                if (nnz == 1) {
                    oneEntryOnlyRows[tmp] = row;
                    tmp++;
                }
            }
        }
        oneEntryOnlyRows.resize(tmp);
        return oneEntryOnlyRows;
    }


    template <class SC,class LO>
    bool ismultiple(ArrayView<SC> A,
                    ArrayView<SC> B)
    {
        Array<LO> zeros;
        Array<LO> nonzeros;

        FROSCH_ASSERT(A.size()==B.size(),"Cannot be multiple (size)");

        // Search for non-zeros
        for (unsigned i=0; i<A.size(); i++) {
            if (fabs(A[i])<1.0e-12) {
                zeros.push_back(i);
            } else {
                nonzeros.push_back(i);
            }
        }

        // Search for non-zeros
        for (unsigned i=0; i<zeros.size(); i++) {
            if (fabs(B[zeros.at(i)])>=1.0e-12) {
                return false;
            }
        }

        // Check nonzeros for multiple
        double mult = B[nonzeros.at(0)]/A[nonzeros.at(0)], mult2;
        for (unsigned i=1; i<nonzeros.size(); i++) {
            mult2 = B[nonzeros.at(i)]/A[nonzeros.at(i)];
            if (fabs(mult2-mult)>=1.0e-12) {
                return false;
            }
        }
        return true;
    }

    template<class T>
    inline void sortunique(T &v)
    {
        std::sort(v.begin(),v.end());
        v.erase(std::unique(v.begin(),v.end()),v.end());
    }

    template <class SC, class LO,class GO,class NO>
    RCP<MultiVector<SC,LO,GO,NO> > ModifiedGramSchmidt(RCP<const MultiVector<SC,LO,GO,NO> > multiVector,
                                                       ArrayView<unsigned> zero)
    {
        FROSCH_TIMER_START(modifiedGramSchmidtTime,"ModifiedGramSchmidt");
        /*
         n = size(V,1);
         k = size(V,2);
         U = zeros(n,k);
         U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
         for i = 2:k
         U(:,i) = V(:,i);
         for j = 1:i-1
         U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
         end
         U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
         end
         */
        unsigned numVec = multiVector->getNumVectors();
        Array<unsigned> arrayZero(0);
        RCP<MultiVector<SC,LO,GO,NO> > resultMultiVector;
        if (numVec>0) {
            unsigned itmp = 0;
            SC en = ScalarTraits<SC>::zero();
            SC de = ScalarTraits<SC>::zero();
            SC norm = ScalarTraits<SC>::zero();
            RCP<MultiVector<SC,LO,GO,NO> > tmpMultiVector = MultiVectorFactory<SC,LO,GO,NO>::Build(multiVector->getMap(),numVec);
            for (unsigned i=0; i<numVec; i++) {
                tmpMultiVector->getVectorNonConst(i-itmp)->update(ScalarTraits<SC>::one(),*multiVector->getVector(i),ScalarTraits<SC>::zero());
                for (unsigned j=0; j<i-itmp; j++) {
                    en = tmpMultiVector->getVector(i-itmp)->dot(*tmpMultiVector->getVector(j));
                    de = tmpMultiVector->getVector(j)->dot(*tmpMultiVector->getVector(j));
                    tmpMultiVector->getVectorNonConst(i-itmp)->update(-en/de,*tmpMultiVector->getVector(j),ScalarTraits<SC>::one());
                }
                norm = tmpMultiVector->getVector(i-itmp)->norm2();
                if (norm<1.0e-10) {
                    arrayZero.push_back(i);
                    itmp++;
                } else {
                    //tmpMultiVector->getVectorNonConst(i-itmp)->scale(1.0/norm);
                }
            }
            resultMultiVector = MultiVectorFactory<SC,LO,GO,NO>::Build(multiVector->getMap(),numVec);
            for (unsigned i=0; i<numVec-itmp; i++) {
                resultMultiVector->getVectorNonConst(i)->update(ScalarTraits<SC>::one(),*tmpMultiVector->getVector(i),ScalarTraits<SC>::zero());
            }
        }
        zero = arrayZero();
        return resultMultiVector;
    }

    template <class SC, class LO,class GO,class NO>
    RCP<const MultiVector<SC,LO,GO,NO> > BuildNullSpace(unsigned dimension,
                                                        unsigned nullSpaceType,
                                                        RCP<const Map<LO,GO,NO> > repeatedMap,
                                                        unsigned dofsPerNode,
                                                        ArrayRCP<RCP<const Map<LO,GO,NO> > > dofsMaps,
                                                        RCP<const MultiVector<SC,LO,GO,NO> > nodeList)
    {
        FROSCH_TIMER_START(buildNullSpaceTime,"BuildNullSpace");
        /*
         Here, the nodeList has to be ordered in accordence to the dofsMaps.
         */
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode.");

        RCP<MultiVector<SC,LO,GO,NO> > nullSpaceBasis;
        if (nullSpaceType==0) { // n-dimensional Laplace
            nullSpaceBasis = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,dofsPerNode);
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = ScalarTraits<SC>::one();
                }
            }
        } else if (nullSpaceType==1) { // linear elasticity
            FROSCH_ASSERT(!nodeList.is_null(),"nodeList.is_null()==true. Cannot build the null space for linear elasticity.");
            FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"nodeList->getNumVectors()!=dimension.");
            FROSCH_ASSERT(dofsPerNode==dimension,"dofsPerNode==dimension.");

            if (dimension==2) {
                nullSpaceBasis = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,3);
                // translations
                for (unsigned i=0; i<2; i++) {
                    for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                        nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = ScalarTraits<SC>::one();
                    }
                }
                // rotation
                for (unsigned j=0; j<dofsMaps[0]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(2)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = -nodeList->getData(1)[j];
                    nullSpaceBasis->getDataNonConst(2)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = nodeList->getData(0)[j];
                }
            } else if (dimension==3) {
                nullSpaceBasis = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,6);
                // translations
                for (unsigned i=0; i<3; i++) {
                    for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                        nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = ScalarTraits<SC>::one();
                    }
                }
                // rotations
                for (unsigned j=0; j<dofsMaps[0]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = nodeList->getData(1)[j];
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = -nodeList->getData(0)[j];
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = ScalarTraits<SC>::zero();

                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = -nodeList->getData(2)[j];
                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = ScalarTraits<SC>::zero();
                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = nodeList->getData(0)[j];

                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = ScalarTraits<SC>::zero();
                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = nodeList->getData(2)[j];
                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = -nodeList->getData(1)[j];
                }
            }
        } else {
            FROSCH_ASSERT(false,"NullSpaceType unknown.");
        }
        return nullSpaceBasis;
    }

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class SC,class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > ConvertToXpetra<SC,LO,GO,NO>::ConvertMap(UnderlyingLib lib,
                                                                 const Epetra_BlockMap &map,
                                                                 RCP<const Comm<int> > comm)
    {
        FROSCH_ASSERT(false,"FROSch::ConvertToXpetra : ERROR: Needs specialization.");
        return null;
    }

    template <class SC,class LO,class GO,class NO>
    RCP<Matrix<SC,LO,GO,NO> > ConvertToXpetra<SC,LO,GO,NO>::ConvertMatrix(UnderlyingLib lib,
                                                                          Epetra_CrsMatrix &matrix,
                                                                          RCP<const Comm<int> > comm)
    {
        FROSCH_ASSERT(false,"FROSch::ConvertToXpetra : ERROR: Needs specialization.");
        return null;
    }

    template <class SC,class LO,class GO,class NO>
    RCP<MultiVector<SC,LO,GO,NO> > ConvertToXpetra<SC,LO,GO,NO>::ConvertMultiVector(UnderlyingLib lib,
                                                                                    Epetra_MultiVector &vector,
                                                                                    RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMultiVectorTime,"ConvertToXpetra::ConvertMultiVector");
        RCP<Map<LO,GO,NO> > map = ConvertToXpetra<SC,LO,GO,NO>::ConvertMap(lib,vector.Map(),comm);
        RCP<MultiVector<SC,LO,GO,NO> > xMultiVector = MultiVectorFactory<SC,LO,GO,NO>::Build(map,vector.NumVectors());
        for (LO i=0; i<vector.NumVectors(); i++) {
            for (LO j=0; j<vector.MyLength(); j++) {
                xMultiVector->getDataNonConst(i)[j] = vector[i][j];
            }
        }
        return xMultiVector;
    }

    template <class SC,class LO,class NO>
    RCP<Map<LO,int,NO> > ConvertToXpetra<SC,LO,int,NO>::ConvertMap(UnderlyingLib lib,
                                                                   const Epetra_BlockMap &map,
                                                                   RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMapTime,"ConvertToXpetra::ConvertMap");
        ArrayView<int> mapArrayView(map.MyGlobalElements(),map.NumMyElements());
        return MapFactory<LO,int,NO>::Build(lib,-1,mapArrayView,0,comm);
    }

    template <class SC,class LO,class NO>
    RCP<Matrix<SC,LO,int,NO> > ConvertToXpetra<SC,LO,int,NO>::ConvertMatrix(UnderlyingLib lib,
                                                                            Epetra_CrsMatrix &matrix,
                                                                            RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMatrixTime,"ConvertToXpetra::ConvertMatrix");
        RCP<Map<LO,int,NO> > rowMap = ConvertToXpetra<SC,LO,int,NO>::ConvertMap(lib,matrix.RowMap(),comm);
        RCP<Matrix<SC,LO,int,NO> > xmatrix = MatrixFactory<SC,LO,int,NO>::Build(rowMap,matrix.MaxNumEntries());
        for (unsigned i=0; i<xmatrix->getNodeNumRows(); i++) {
            LO numEntries;
            LO* indices;
            SC* values;
            matrix.ExtractMyRowView(i,numEntries,values,indices);

            Array<int> indicesArray(numEntries);
            ArrayView<SC> valuesArrayView(values,numEntries);
            for (LO j=0; j<numEntries; j++) {
                indicesArray[j] = matrix.ColMap().GID(indices[j]);
            }
            xmatrix->insertGlobalValues(matrix.RowMap().GID(i),indicesArray(),valuesArrayView);
        }
        xmatrix->fillComplete();
        return xmatrix;
    }

    template <class SC,class LO,class NO>
    RCP<MultiVector<SC,LO,int,NO> > ConvertToXpetra<SC,LO,int,NO>::ConvertMultiVector(UnderlyingLib lib,
                                                                                      Epetra_MultiVector &vector,
                                                                                      RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMultiVectorTime,"ConvertToXpetra::ConvertMultiVector");
        RCP<Map<LO,int,NO> > map = ConvertToXpetra<SC,LO,int,NO>::ConvertMap(lib,vector.Map(),comm);
        RCP<MultiVector<SC,LO,int,NO> > xMultiVector = MultiVectorFactory<SC,LO,int,NO>::Build(map,vector.NumVectors());
        for (LO i=0; i<vector.NumVectors(); i++) {
            for (LO j=0; j<vector.MyLength(); j++) {
                xMultiVector->getDataNonConst(i)[j] = vector[i][j];
            }
        }
        return xMultiVector;
    }

    template <class SC,class LO,class NO>
    RCP<Map<LO,long long,NO> > ConvertToXpetra<SC,LO,long long,NO>::ConvertMap(UnderlyingLib lib,
                                                                               const Epetra_BlockMap &map,
                                                                               RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMapTime,"ConvertToXpetra::ConvertMap");
        ArrayView<long long> mapArrayView(map.MyGlobalElements64(),map.NumMyElements());
        return MapFactory<LO,long long,NO>::Build(lib,-1,mapArrayView,0,comm);
    }

    template <class SC,class LO,class NO>
    RCP<Matrix<SC,LO,long long,NO> > ConvertToXpetra<SC,LO,long long,NO>::ConvertMatrix(UnderlyingLib lib,
                                                                                        Epetra_CrsMatrix &matrix,
                                                                                        RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMatrixTime,"ConvertToXpetra::ConvertMatrix");
        RCP<Map<LO,long long,NO> > rowMap = ConvertToXpetra<SC,LO,long long,NO>::ConvertMap(lib,matrix.RowMap(),comm);
        RCP<Matrix<SC,LO,long long,NO> > xmatrix = MatrixFactory<SC,LO,long long,NO>::Build(rowMap,matrix.MaxNumEntries());
        for (unsigned i=0; i<xmatrix->getNodeNumRows(); i++) {
            LO numEntries;
            LO* indices;
            SC* values;
            matrix.ExtractMyRowView(i,numEntries,values,indices);

            Array<long long> indicesArray(numEntries);
            ArrayView<SC> valuesArrayView(values,numEntries);
            for (LO j=0; j<numEntries; j++) {
                indicesArray[j] = matrix.ColMap().GID64(indices[j]);
            }
            xmatrix->insertGlobalValues(matrix.RowMap().GID64(i),indicesArray(),valuesArrayView);
        }
        xmatrix->fillComplete();
        return xmatrix;
    }

    template <class SC,class LO,class NO>
    RCP<MultiVector<SC,LO,long long,NO> > ConvertToXpetra<SC,LO,long long,NO>::ConvertMultiVector(UnderlyingLib lib,
                                                                                                  Epetra_MultiVector &vector,
                                                                                                  RCP<const Comm<int> > comm)
    {
        FROSCH_TIMER_START(convertMultiVectorTime,"ConvertToXpetra::ConvertMultiVector");
        RCP<Map<LO,long long,NO> > map = ConvertToXpetra<SC,LO,long long,NO>::ConvertMap(lib,vector.Map(),comm);
        RCP<MultiVector<SC,LO,long long,NO> > xMultiVector = MultiVectorFactory<SC,LO,long long,NO>::Build(map,vector.NumVectors());
        for (LO i=0; i<vector.NumVectors(); i++) {
            for (LO j=0; j<vector.MyLength(); j++) {
                xMultiVector->getDataNonConst(i)[j] = vector[i][j];
            }
        }
        return xMultiVector;
    }
#endif

    template <class Type>
    RCP<Type> ExtractPtrFromParameterList(ParameterList& paramList,
                                          std::string namePtr)
    {
        FROSCH_TIMER_START(extractPtrFromParameterListTime,"ExtractPtrFromParameterList");
        RCP<Type> pointer = null;

        if (paramList.isParameter(namePtr) == false)
            return pointer;

        if (paramList.isType<decltype(pointer)>(namePtr)) {
            pointer = paramList.get<decltype(pointer)>(namePtr);
        }

        return pointer;
    }

    template <class Type>
    ArrayRCP<Type> ExtractVectorFromParameterList(ParameterList& paramList,
                                                  std::string nameVector)
    {
        FROSCH_TIMER_START(extractVectorFromParameterListTime,"ExtractVectorFromParameterList");
        ArrayRCP<Type> vector = null;

        if (paramList.isParameter(nameVector) == false)
            return vector;

        if (paramList.isType<decltype(vector)>(nameVector)) {
            vector = paramList.get<decltype(vector)>(nameVector);
        }

        return vector;
    }

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class LO,class GO,class NO>
    RCP<Epetra_Map> ConvertToEpetra(const Map<LO,GO,NO> &map,
                                    RCP<Epetra_Comm> epetraComm)
    {
        FROSCH_TIMER_START(convertToEpetraTime,"ConvertToEpetra");
        ArrayView<const GO> elementList = map.getNodeElementList();

        GO numGlobalElements = map.getGlobalNumElements();

        if(elementList.size()>0)
            return rcp(new Epetra_Map(numGlobalElements,elementList.size(),elementList.getRawPtr(),0,*epetraComm));
        else
            return rcp(new Epetra_Map(numGlobalElements,0,NULL,0,*epetraComm));
    }

    template <class SC, class LO, class GO,class NO>
    RCP<Epetra_MultiVector> ConvertToEpetra(const MultiVector<SC,LO,GO,NO> &vector,
                                            RCP<Epetra_Comm> epetraComm)
    {
        FROSCH_TIMER_START(convertToEpetraTime,"ConvertToEpetra");
        RCP<Epetra_Map> map = ConvertToEpetra<LO,GO,NO>(*vector.getMap(),epetraComm);
        RCP<Epetra_MultiVector > multiVector(new Epetra_MultiVector(*map,vector.getNumVectors()));
        for (LO i=0; i<vector.getNumVectors(); i++) {
            for (LO j=0; j<vector.getLocalLength(); j++) {
                (*multiVector)[i][j] = vector.getData(i)[j];
            }
        }
        return multiVector;
    }

    template <class SC, class LO,class GO, class NO>
    RCP<Epetra_CrsMatrix> ConvertToEpetra(const Matrix<SC,LO,GO,NO> &matrix,
                                          RCP<Epetra_Comm> epetraComm)
    {
        FROSCH_TIMER_START(convertToEpetraTime,"ConvertToEpetra");
        RCP<Epetra_Map> map = ConvertToEpetra<LO,GO,NO>(*matrix.getMap(),epetraComm);
        RCP<Epetra_CrsMatrix> matrixEpetra(new Epetra_CrsMatrix(::Copy,*map,matrix.getGlobalMaxNumRowEntries()));
        ArrayView<const SC> valuesArrayView;
        ArrayView<const LO> indicesArrayView;

        for (LO i=0; i<(LO) matrix.getRowMap()->getNodeNumElements(); i++) {
            matrix.getLocalRowView(i, indicesArrayView, valuesArrayView);
            Array<GO> indicesGlobal(indicesArrayView.size());
            for (LO j=0; j<indicesArrayView.size(); j++) {
                indicesGlobal[j] = matrix.getColMap()->getGlobalElement(indicesArrayView[j]);
            }
            if (indicesArrayView.size()>0) {
                matrixEpetra->InsertGlobalValues(matrix.getRowMap()->getGlobalElement(i), indicesArrayView.size(), &(valuesArrayView[0]), &(indicesGlobal[0]));
            }
        }
        matrixEpetra->FillComplete();
        return matrixEpetra;
    }
#endif

    template <class LO>
    Array<LO> GetIndicesFromString(std::string string, LO dummy)
    {
        Array<LO> indices(0);
        for (unsigned i=0; i<string.length(); i++) {
            indices.push_back((LO) stoi(string.substr(i,i+1)));
        }
        return indices;
    }

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
    template <class SC, class LO,class GO,class NO>
    int RepartionMatrixZoltan2(RCP<Matrix<SC,LO,GO,NO> > &crsMatrix,
                               RCP<ParameterList> parameterList)
    {
        FROSCH_TIMER_START(repartionMatrixZoltan2Time,"RepartionMatrixZoltan2");
        RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));

        using inputAdapter    = Zoltan2::XpetraCrsMatrixAdapter<CrsMatrix<SC,LO,GO,NO> >;

        RCP<CrsMatrixWrap<SC,LO,GO,NO> > tmpCrsWrap = rcp_dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO> >(crsMatrix);
        RCP<CrsMatrix<SC,LO,GO,NO> > tmpCrsMatrix = tmpCrsWrap->getCrsMatrix();
        inputAdapter adaptedMatrix(tmpCrsMatrix);

        RCP<Zoltan2::PartitioningProblem<inputAdapter> > problem =
             rcp(new Zoltan2::PartitioningProblem<inputAdapter>(&adaptedMatrix, parameterList.get()));

        problem->solve();

        RCP<CrsMatrix<SC,LO,GO,NO> > matrixRepartition;
        adaptedMatrix.applyPartitioningSolution(*tmpCrsMatrix,matrixRepartition,problem->getSolution());

        RCP<CrsMatrixWrap<SC,LO,GO,NO> > tmpCrsWrap2 = rcp(new CrsMatrixWrap<SC,LO,GO,NO>(matrixRepartition));
        crsMatrix = rcp_dynamic_cast<Matrix<SC,LO,GO,NO> >(tmpCrsWrap2);
        return 0;
    }
#endif
} // namespace FROSch

#endif
