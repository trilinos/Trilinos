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

namespace FROSch {
    
    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildUniqueMap(const Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > map)
    {
        Teuchos::RCP<Xpetra::Vector<GO,LO,GO,NO> > myIndices = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map);
        myIndices->putScalar(map->getComm()->getRank()+1);
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > linearMap = Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),map->getMaxAllGlobalIndex()+1,0,map->getComm());
        Teuchos::RCP<Xpetra::Vector<GO,LO,GO,NO> > globalIndices = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(linearMap);
        
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > importer = Xpetra::ImportFactory<LO,GO,NO>::Build(map,linearMap);
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > importer2 = Xpetra::ImportFactory<LO,GO,NO>::Build(linearMap,map); // AH 10/16/2017: Ist der notwendig??? Mit Epetra ging es auch ohne einen zweiten Importer und stattdessen mit einem Export
        globalIndices->doImport(*myIndices,*importer,Xpetra::INSERT);
        myIndices->putScalar(0);
        myIndices->doImport(*globalIndices,*importer2,Xpetra::ADD);
        
        Teuchos::Array<GO> uniqueVector;
        for (unsigned i=0; i<myIndices->getLocalLength(); i++) {
            if (myIndices->getData(0)[i] == map->getComm()->getRank()+1) {
                uniqueVector.push_back(map->getGlobalElement(i));
            }
        }
        
        return Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),-1,uniqueVector(),0,map->getComm());
    }
    
    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildRepeatedMap(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > matrix)
    {
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > uniqueMap = Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap(),1);
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > overlappingMap = uniqueMap;
        ExtendOverlapByOneLayer<SC,LO,GO,NO>(matrix,overlappingMap);
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > tmpMatrix = matrix;
        matrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,2*tmpMatrix->getGlobalMaxNumRowEntries());
#ifdef Tpetra_issue_1752
        // AH 10/10/2017: Can we get away with using just one importer/exporter after the Tpetra issue is fixed?
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);
        Teuchos::RCP<Xpetra::Export<LO,GO,NO> > gather = Xpetra::ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);
        
        matrix->doImport(*tmpMatrix,*scatter,Xpetra::ADD);
#else
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter;
        Teuchos::RCP<Xpetra::Export<LO,GO,NO> > gather = Xpetra::ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);
        
        if (matrix->getRowMap()->lib()==Xpetra::UseEpetra) {
            scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);
            matrix->doImport(*tmpMatrix,*scatter,Xpetra::ADD);
        } else {
            matrix->doImport(*tmpMatrix,*gather,Xpetra::ADD);
        }
#endif
        
        Teuchos::Array<SC> one(1,Teuchos::ScalarTraits<SC>::one());
        Teuchos::Array<GO> myPID(1,uniqueMap->getComm()->getRank());
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        
        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            if (uniqueMap->getLocalElement(globalRow)<0) {
                Teuchos::ArrayView<const GO> indices;
                Teuchos::ArrayView<const SC> values;
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
        commMatTmp->doExport(*commMat,*gather,Xpetra::INSERT);
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp2 = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            GO globalRow = uniqueMap->getGlobalElement(i);
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            commMatTmp->getGlobalRowView(globalRow,indices,values);
            
            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    Teuchos::Array<GO> pID(1,indices[j]);
                    if (pID<myPID) {
                        commMatTmp2->insertGlobalValues(globalRow,pID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                    }
                }
            }
        }
        commMatTmp2->fillComplete();
        commMatTmp.reset();
        commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
#ifdef Tpetra_issue_1752
        commMat->doImport(*commMatTmp2,*scatter,Xpetra::ADD);
#else
        commMat->doImport(*commMatTmp2,*gather,Xpetra::ADD);
#endif
        
        Teuchos::ArrayView<const GO> myGlobalElements = uniqueMap->getNodeElementList();
        Teuchos::Array<GO> repeatedIndices(uniqueMap->getNodeNumElements());
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            repeatedIndices.at(i) = myGlobalElements[i];
        }
        
        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
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
        return Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,repeatedIndices(),0,matrix->getRowMap()->getComm());
    }
    
    /*
    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildRepeatedMap(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > matrix)
    {
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > uniqueMap = Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap(),1);
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > overlappingMap = uniqueMap;
        ExtendOverlapByOneLayer<SC,LO,GO,NO>(matrix,overlappingMap);
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > tmpMatrix = matrix;
        matrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,2*tmpMatrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(uniqueMap,overlappingMap);
        
        matrix->doImport(*tmpMatrix,*scatter,Xpetra::ADD);
        
        Teuchos::Array<SC> one(1,Teuchos::ScalarTraits<SC>::one());
        Teuchos::Array<GO> myPID(1,uniqueMap->getComm()->getRank());
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);
        
        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            if (uniqueMap->getLocalElement(globalRow)<0) {
                Teuchos::ArrayView<const GO> indices;
                Teuchos::ArrayView<const SC> values;
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
        commMatTmp->doExport(*commMat,*commExporter,Xpetra::INSERT);
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp2 = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,10);
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            GO globalRow = uniqueMap->getGlobalElement(i);
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            commMatTmp->getGlobalRowView(globalRow,indices,values);
            
            if (indices.size()>0) {
                for (LO j=0; j<indices.size(); j++) {
                    Teuchos::Array<GO> pID(1,indices[j]);
                    if (pID<myPID) {
                        commMatTmp2->insertGlobalValues(globalRow,pID(),one()); // Hier werden immer nur einzelne Werte eingefügt -> Das geht schneller, wenn es Zeilenwise gemacht wird
                    }
                }
            }
        }
        commMatTmp2->fillComplete();
        commMatTmp.reset();
        commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,10);
        commMat->doImport(*commMatTmp2,*commExporter,Xpetra::ADD);
        
        Teuchos::ArrayView<const GO> myGlobalElements = uniqueMap->getNodeElementList();
        Teuchos::Array<GO> repeatedIndices(uniqueMap->getNodeNumElements());
        for (unsigned i=0; i<uniqueMap->getNodeNumElements(); i++) {
            repeatedIndices.at(i) = myGlobalElements[i];
        }
        
        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            GO globalRow = overlappingMap->getGlobalElement(i);
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
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
        return Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,repeatedIndices(),0,matrix->getRowMap()->getComm());
    }
    */
    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &overlappingMatrix,
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &overlappingMap)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > tmpMatrix = overlappingMatrix;
        overlappingMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,2*tmpMatrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(tmpMatrix->getRowMap(),overlappingMap);
        overlappingMatrix->doImport(*tmpMatrix,*scatter,Xpetra::ADD);
        
        Teuchos::Array<GO> indicesOverlappingSubdomain(0);
        for (unsigned i=0; i<overlappingMap->getNodeNumElements(); i++) {
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            overlappingMatrix->getGlobalRowView(overlappingMap->getGlobalElement(i),indices,values);
            
            for (LO j=0; j<indices.size(); j++) {
                indicesOverlappingSubdomain.push_back(indices[j]);
            }
        }
        sortunique(indicesOverlappingSubdomain);
        overlappingMap = Xpetra::MapFactory<LO,GO,NO>::Build(overlappingMap->lib(),-1,indicesOverlappingSubdomain(),0,overlappingMap->getComm());
        overlappingMatrix->fillComplete(tmpMatrix->getDomainMap(),tmpMatrix->getRangeMap());
        
        return 0;
    }
    
    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > AssembleMaps(Teuchos::ArrayView<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > mapVector,
                                                      Teuchos::ArrayRCP<Teuchos::ArrayRCP<LO> > &partMappings)
    {
        FROSCH_ASSERT(mapVector.size()>0,"Length of mapVector is == 0!");
        LO i = 0;
        LO localstart = 0;
        LO sizetmp = 0;
        LO size = 0;
        GO globalstart = 0;
        
        partMappings = Teuchos::ArrayRCP<Teuchos::ArrayRCP<LO> >(mapVector.size());

        Teuchos::ArrayRCP<GO> assembledMapTmp(0);
        for (unsigned j=0; j<mapVector.size(); j++) {
            sizetmp = mapVector[j]->getNodeNumElements();
            partMappings[j] = Teuchos::ArrayRCP<LO>(sizetmp);
            
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
            
            if (mapVector[0]->lib()==Xpetra::UseEpetra || mapVector[j]->getGlobalNumElements()>0) {
                globalstart += 1;
            }
             */
            
            globalstart += std::max(mapVector[j]->getMaxAllGlobalIndex(),-1)+1; // AH 04/05/2018: mapVector[j]->getMaxAllGlobalIndex() can result in -2147483648 if the map is empty on the process => introducing max(,)
            
            //if (mapVector[j]->getComm()->getRank() == 0) std::cout << std::endl << globalstart << std::endl;
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(mapVector[0]->lib(),-1,assembledMapTmp(),0,mapVector[0]->getComm());
    }
    
    template <class LO,class GO,class NO>
    int BuildDofMaps(const Teuchos::RCP<Xpetra::Map<LO,GO,NO> > map,
                     unsigned dofsPerNode,
                     unsigned dofOrdering,
                     Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &nodesMap,
                     Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > &dofMaps)
    {
        //if (map->getComm()->getRank()==0) std::cout << "WARNING: BuildDofMaps is yet to be tested...\n";
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(map->getGlobalNumElements()%dofsPerNode==0 && map->getNodeNumElements()%dofsPerNode==0,"ERROR: The number of dofsPerNode does not divide the number of global dofs in the map!");
        
        Teuchos::Array<GO> nodes(map->getNodeNumElements()/dofsPerNode);
        Teuchos::Array<Teuchos::ArrayRCP<GO> > dofs(dofsPerNode);
        for (unsigned j=0; j<dofsPerNode; j++) {
            dofs[j] = Teuchos::ArrayRCP<GO>(map->getNodeNumElements()/dofsPerNode);
        }
        if (dofOrdering==0) {
            for (unsigned i=0; i<nodes.size(); i++) {
                nodes[i] = map->getGlobalElement(dofsPerNode*i)/dofsPerNode;
                for (unsigned j=0; j<dofsPerNode; j++) {
                    dofs[j][i] = dofsPerNode*nodes[i]+j;
                }
            }
        } else if (dofOrdering == 1) {
            for (unsigned i=0; i<nodes.size(); i++) {
                nodes[i] = map->getGlobalElement(i);
                for (unsigned j=0; j<dofsPerNode; j++) {
                    dofs[j][i] = nodes[i]+j*(map->getMaxAllGlobalIndex()+1)/dofsPerNode;
                }
            }
        } else {
            FROSCH_ASSERT(0!=0,"dofOrdering unknown.");
        }
        nodesMap = Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),-1,nodes(),0,map->getComm());
        
        dofMaps = Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > >(dofsPerNode);
        for (unsigned j=0; j<dofsPerNode; j++) {
            dofMaps[j] = Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),-1,dofs[j](),0,map->getComm());
        }
        return 0;
    }
    
    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildMapFromDofMaps(const Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > &dofMaps,
                                                             unsigned dofsPerNode,
                                                             unsigned dofOrdering)
    {
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(!dofMaps.is_null(),"dofMaps.is_null().");
        FROSCH_ASSERT(dofMaps.size()==dofsPerNode,"dofMaps.size!=dofsPerNode.");
        for (unsigned i=0; i<dofMaps.size(); i++) {
            FROSCH_ASSERT(dofMaps[i]->getGlobalNumElements()%dofsPerNode==0 && dofMaps[i]->getNodeNumElements()%dofsPerNode==0,"ERROR: The number of dofsPerNode does not divide the number of global dofs in the dofMaps!");
        }
        
        unsigned numNodes = dofMaps[0]->getNodeNumElements();
        Teuchos::Array<GO> globalIDs(numNodes);
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
            FROSCH_ASSERT(0!=0,"dofOrdering unknown.");
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(dofMaps[0]->lib(),-1,globalIDs(),0,dofMaps[0]->getComm());
    }
    
    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildMapFromNodeMap(Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &nodesMap,
                                                             unsigned dofsPerNode,
                                                             unsigned dofOrdering)
    {
        FROSCH_ASSERT(dofOrdering==0 || dofOrdering==1,"ERROR: Specify a valid DofOrdering.");
        FROSCH_ASSERT(!nodesMap.is_null(),"nodesMap.is_null().");
        
        unsigned numNodes = nodesMap->getNodeNumElements();
        Teuchos::Array<GO> globalIDs(dofsPerNode*numNodes);
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
            FROSCH_ASSERT(0!=0,"dofOrdering unknown.");
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(nodesMap->lib(),-1,globalIDs(),0,nodesMap->getComm());
    }
    
    template <class SC,class LO,class GO,class NO>
    Teuchos::ArrayRCP<GO> FindOneEntryOnlyRowsGlobal(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &matrix,
                                                     Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &repeatedMap)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > repeatedMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(repeatedMap,2*matrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(matrix->getRowMap(),repeatedMap);
        repeatedMatrix->doImport(*matrix,*scatter,Xpetra::ADD);
        
        Teuchos::ArrayRCP<GO> oneEntryOnlyRows(repeatedMatrix->getNodeNumRows());
        LO tmp = 0;
        LO nnz;
        GO row;
        for (unsigned i=0; i<repeatedMatrix->getNodeNumRows(); i++) {
            row = repeatedMap->getGlobalElement(i);
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
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
    bool ismultiple(Teuchos::ArrayView<SC> A,
                    Teuchos::ArrayView<SC> B)
    {
        Teuchos::Array<LO> zeros;
        Teuchos::Array<LO> nonzeros;
        
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
    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > ModifiedGramSchmidt(Teuchos::RCP<const Xpetra::MultiVector<SC,LO,GO,NO> > multiVector,
                                                                        Teuchos::ArrayView<unsigned> zero)
    {
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
        Teuchos::Array<unsigned> arrayZero(0);
        Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > resultMultiVector;
        if (numVec>0) {
            unsigned itmp = 0;
            SC en = 0.0;
            SC de = 0.0;
            SC norm = 0.0;
            Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > tmpMultiVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(multiVector->getMap(),numVec);
            for (unsigned i=0; i<numVec; i++) {
                tmpMultiVector->getVectorNonConst(i-itmp)->update(1.0,*multiVector->getVector(i),0.0);
                for (unsigned j=0; j<i-itmp; j++) {
                    en = tmpMultiVector->getVector(i-itmp)->dot(*tmpMultiVector->getVector(j));
                    de = tmpMultiVector->getVector(j)->dot(*tmpMultiVector->getVector(j));
                    tmpMultiVector->getVectorNonConst(i-itmp)->update(-en/de,*tmpMultiVector->getVector(j),1.0);
                }
                norm = tmpMultiVector->getVector(i-itmp)->norm2();
                if (norm<1.0e-10) {
                    arrayZero.push_back(i);
                    itmp++;
                } else {
                    //tmpMultiVector->getVectorNonConst(i-itmp)->scale(1.0/norm);
                }
            }
            resultMultiVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(multiVector->getMap(),numVec);
            for (unsigned i=0; i<numVec-itmp; i++) {
                resultMultiVector->getVectorNonConst(i)->update(1.0,*tmpMultiVector->getVector(i),0.0);
            }
        }
        zero = arrayZero();
        return resultMultiVector;
    }
    
    template <class SC, class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > BuildNullSpace(unsigned dimension,
                                                                   unsigned nullSpaceType,
                                                                   Teuchos::RCP<Xpetra::Map<LO,GO,NO> > repeatedMap,
                                                                   unsigned dofsPerNode,
                                                                   Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > dofsMaps,
                                                                   Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > nodeList)
    {
        /*
         Here, the nodeList has to be ordered in accordence to the dofsMaps.
         */
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode.");
        
        Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > nullSpaceBasis;
        if (nullSpaceType==0) { // n-dimensional Laplace
            nullSpaceBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,dofsPerNode);
            for (unsigned i=0; i<dofsPerNode; i++) {
                for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = 1.0;
                }
            }
        } else if (nullSpaceType==1) { // linear elasticity
            FROSCH_ASSERT(!nodeList.is_null(),"nodeList.is_null()==true. Cannot build the null space for linear elasticity.");
            FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"nodeList->getNumVectors()!=dimension.");
            FROSCH_ASSERT(dofsPerNode==dimension,"dofsPerNode==dimension.");
            
            if (dimension==2) {
                nullSpaceBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,3);
                // translations
                for (unsigned i=0; i<2; i++) {
                    for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                        nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = 1.0;
                    }
                }
                // rotation
                for (unsigned j=0; j<dofsMaps[0]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(2)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = -nodeList->getData(1)[j];
                    nullSpaceBasis->getDataNonConst(2)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = nodeList->getData(0)[j];
                }
            } else if (dimension==3) {
                nullSpaceBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,6);
                // translations
                for (unsigned i=0; i<3; i++) {
                    for (unsigned j=0; j<dofsMaps[i]->getNodeNumElements(); j++) {
                        nullSpaceBasis->getDataNonConst(i)[repeatedMap->getLocalElement(dofsMaps[i]->getGlobalElement(j))] = 1.0;
                    }
                }
                // rotations
                for (unsigned j=0; j<dofsMaps[0]->getNodeNumElements(); j++) {
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = nodeList->getData(1)[j];
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = -nodeList->getData(0)[j];
                    nullSpaceBasis->getDataNonConst(3)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = 0.0;
                    
                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = -nodeList->getData(2)[j];
                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = 0.0;
                    nullSpaceBasis->getDataNonConst(4)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = nodeList->getData(0)[j];
                    
                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[0]->getGlobalElement(j))] = 0.0;
                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[1]->getGlobalElement(j))] = nodeList->getData(2)[j];
                    nullSpaceBasis->getDataNonConst(5)[repeatedMap->getLocalElement(dofsMaps[2]->getGlobalElement(j))] = -nodeList->getData(1)[j];
                }
            }
        } else {
            FROSCH_ASSERT(0!=0,"NullSpaceType unknown.");
        }
        return nullSpaceBasis;
    }
    
    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > ConvertToXpetra(Xpetra::UnderlyingLib lib,
                                                         const Epetra_BlockMap &map,
                                                         Teuchos::RCP<const Teuchos::Comm<int> > comm)
    {
        Teuchos::ArrayView<GO> mapArrayView(map.MyGlobalElements(),map.NumMyElements());
        return Xpetra::MapFactory<LO,GO,NO>::Build(lib,-1,mapArrayView,0,comm);
    }
    
    template <class SC, class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > ConvertToXpetra(Xpetra::UnderlyingLib lib,
                                                               Epetra_CrsMatrix &matrix,
                                                               Teuchos::RCP<const Teuchos::Comm<int> > comm)
    {
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > rowMap = ConvertToXpetra<LO,GO,NO>(lib,matrix.RowMap(),comm);
        
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > xmatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(rowMap,matrix.MaxNumEntries());
        for (unsigned i=0; i<xmatrix->getNodeNumRows(); i++) {
            LO numEntries;
            GO* indices;
            SC* values;
            matrix.ExtractMyRowView(i,numEntries,values,indices);
            
            Teuchos::Array<GO> indicesArray(numEntries);
            Teuchos::ArrayView<SC> valuesArrayView(values,numEntries);
            for (LO j=0; j<numEntries; j++) {
                indicesArray[j] = matrix.ColMap().GID(indices[j]);
            }
            xmatrix->insertGlobalValues(matrix.RowMap().GID(i),indicesArray(),valuesArrayView);
        }
        xmatrix->fillComplete();
        return xmatrix;
    }
    
    template <class SC, class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > ConvertToXpetra(Xpetra::UnderlyingLib lib,
                                                                    Epetra_MultiVector &vector,
                                                                    Teuchos::RCP<const Teuchos::Comm<int> > comm)
    {
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > map = ConvertToXpetra<LO,GO,NO>(lib,vector.Map(),comm);
        Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xMultiVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(map,vector.NumVectors());
        for (LO i=0; i<vector.NumVectors(); i++) {
            for (LO j=0; j<vector.MyLength(); j++) {
                xMultiVector->getDataNonConst(i)[j] = vector[i][j];
            }
        }
        return xMultiVector;
    }
}

#endif
