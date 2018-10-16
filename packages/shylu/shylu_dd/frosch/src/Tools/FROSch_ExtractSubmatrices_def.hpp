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

#ifndef _FROSCH_EXTRACTSUBMATRICES_DEF_HPP
#define _FROSCH_EXTRACTSUBMATRICES_DEF_HPP

#include <FROSch_ExtractSubmatrices_decl.hpp>
namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                           Teuchos::RCP<Xpetra::Map<LO,GO,NO> > map)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > subdomainMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,Xpetra::ADD);
        //cout << *subdomainMatrix << std::endl;
        Teuchos::RCP<const Teuchos::Comm<LO> > SerialComm = rcp(new Teuchos::MpiComm<LO>(MPI_COMM_SELF));
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > localSubdomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),map->getNodeNumElements(),0,SerialComm);
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > localSubdomainMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getNodeMaxNumRowEntries());
        
        for (unsigned i=0; i<localSubdomainMap->getNodeNumElements(); i++) {
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);
            
            LO size = indices.size();
            if (size>0) {
                Teuchos::Array<GO> indicesLocal;
                Teuchos::Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesLocal.push_back(localIndex);
                        valuesLocal.push_back(values[j]);
                    }
                }
                localSubdomainMatrix->insertGlobalValues(i,indicesLocal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();        
        return localSubdomainMatrix;
    }
    
    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                           Teuchos::RCP<Xpetra::Map<LO,GO,NO> > map,
                                                                           SC value)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > subdomainMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,Xpetra::ADD);
        //cout << *subdomainMatrix << std::endl;
        Teuchos::RCP<const Teuchos::Comm<LO> > SerialComm = rcp(new Teuchos::MpiComm<LO>(MPI_COMM_SELF));
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > localSubdomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(map->lib(),map->getNodeNumElements(),0,SerialComm);
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > localSubdomainMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getNodeMaxNumRowEntries());
        
        for (unsigned i=0; i<localSubdomainMap->getNodeNumElements(); i++) {
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);
            
            LO size = indices.size();
            if (size>0) {
                Teuchos::Array<GO> indicesGlobal;
                Teuchos::Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesGlobal.push_back(localIndex);
                        valuesLocal.push_back(value);
                    }
                }
                localSubdomainMatrix->insertGlobalValues(i,indicesGlobal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();
        return localSubdomainMatrix;
    }
    
    template <class SC,class LO,class GO,class NO>
    int UpdateLocalSubdomainMatrix(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > globalMatrix,
                                   Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &map,
                                   Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &localSubdomainMatrix)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > subdomainMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,Xpetra::ADD);
        
        localSubdomainMatrix->resumeFill();
        for (unsigned i=0; i<map->getNodeNumElements(); i++) {
            Teuchos::ArrayView<const GO> indices;
            Teuchos::ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);
            
            LO size = indices.size();
            if (size>0) {
                Teuchos::Array<LO> indicesLocal;
                Teuchos::Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesLocal.push_back(localIndex);
                        valuesLocal.push_back(values[j]);
                    }
                }
                localSubdomainMatrix->replaceLocalValues(i,indicesLocal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrices(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > k,
                         Teuchos::ArrayView<GO> indI,
                         Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &kII,
                         Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &kIJ,
                         Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &kJI,
                         Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &kJJ)
    {
        // We need four Maps
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapI = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI(),0,k->getRowMap()->getComm());
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapILocal = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI.size(),0,k->getRowMap()->getComm());
        
        Teuchos::Array<GO> indJ;
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            if (mapI->getLocalElement(i)<0) {
                indJ.push_back(i);
            }
        }
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapJ = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indJ(),0,k->getRowMap()->getComm());
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapJLocal = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indJ.size(),0,k->getRowMap()->getComm());
        Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > colMap = k->getColMap();
        kII = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,std::min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        kIJ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,std::min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));
        kJI = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,std::min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        kJJ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,std::min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));
        
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            Teuchos::ArrayView<const LO> indices;
            Teuchos::ArrayView<const SC> values;
            
            k->getLocalRowView(i,indices,values);
            //cout << numEntries << std::endl;
            Teuchos::Array<GO> indicesI;
            Teuchos::Array<SC> valuesI;
            Teuchos::Array<GO> indicesJ;
            Teuchos::Array<SC> valuesJ;
            LO tmp1=mapI->getLocalElement(i);
            LO tmp2=0;
            if (tmp1>=0) {
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = mapI->getLocalElement(colMap->getGlobalElement(indices[j]));
                    if (tmp2>=0) {
                        indicesI.push_back(tmp2);
                        valuesI.push_back(values[j]);
                    } else {
                        indicesJ.push_back(mapJ->getLocalElement(colMap->getGlobalElement(indices[j])));
                        valuesJ.push_back(values[j]);
                    }
                }
                //cout << k->getRowMap().Comm().getRank() << " " << tmp1 << " numEntries " << numEntries << " indicesI.size() " << indicesI.size() << " indicesJ.size() " << indicesJ.size() << std::endl;
                kII->insertGlobalValues(tmp1,indicesI(),valuesI());
                kIJ->insertGlobalValues(tmp1,indicesJ(),valuesJ());
            } else  {
                tmp1=mapJ->getLocalElement((GO) i);
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = mapI->getLocalElement(colMap->getGlobalElement(indices[j]));
                    if (tmp2>=0) {
                        indicesI.push_back(tmp2);
                        valuesI.push_back(values[j]);
                    } else {
                        indicesJ.push_back(mapJ->getLocalElement(colMap->getGlobalElement(indices[j])));
                        valuesJ.push_back(values[j]);
                    }
                }
                kJI->insertGlobalValues(tmp1,indicesI(),valuesI());
                kJJ->insertGlobalValues(tmp1,indicesJ(),valuesJ());
            }
        }
        
        kII->fillComplete(mapILocal,mapILocal);
        kIJ->fillComplete(mapJLocal,mapILocal);
        kJI->fillComplete(mapILocal,mapJLocal);
        kJJ->fillComplete(mapJLocal,mapJLocal);
        
        return 0;
    }
    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrix(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > k,
                       Teuchos::ArrayView<GO> indI,
                       Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &kII)
    {
                Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));
     
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapI = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI(),0,k->getRowMap()->getComm());
        
//        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapILocal = Xpetra::MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI.size(),0,k->getRowMap()->getComm());

        kII = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(mapI,std::min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        GO maxGID = mapI->getMaxAllGlobalIndex();
        GO minGID = mapI->getMinAllGlobalIndex();
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            Teuchos::ArrayView<const LO> indices;
            Teuchos::ArrayView<const SC> values;
            
            k->getLocalRowView(i,indices,values);
            //cout << numEntries << std::endl;
            Teuchos::Array<GO> indicesI;
            Teuchos::Array<SC> valuesI;
            
            LO tmp1=mapI->getLocalElement(k->getRowMap()->getGlobalElement(i));
            GO tmp2=0;
            if (tmp1>=0) {
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = k->getColMap()->getGlobalElement(indices[j]);
                    if (minGID<=tmp2 && tmp2<=maxGID) {
                        indicesI.push_back(tmp2);
                        valuesI.push_back(values[j]);
                    }
                }
                //cout << k->getRowMap().Comm().getRank() << " " << tmp1 << " numEntries " << numEntries << " indicesI.size() " << indicesI.size() << " indicesJ.size() " << indicesJ.size() << std::endl;
                kII->insertGlobalValues(mapI->getGlobalElement(tmp1),indicesI(),valuesI());
            }
        }
        kII->fillComplete(mapI,mapI);
        
        return 0;
    }
    

}

#endif
