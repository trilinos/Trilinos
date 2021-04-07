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

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map)
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime,"ExtractLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);
        //cout << *subdomainMatrix << endl;
        RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
        RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(map->lib(),map->getNodeNumElements(),0,SerialComm);
        RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getGlobalMaxNumRowEntries());

        for (unsigned i=0; i<localSubdomainMap->getNodeNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<GO> indicesLocal;
                Array<SC> valuesLocal;
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
        return localSubdomainMatrix.getConst();
    }

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map,
                                                                SC value)
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime,"ExtractLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);
        //cout << *subdomainMatrix << endl;
        RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
        RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(map->lib(),map->getNodeNumElements(),0,SerialComm);
        RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getGlobalMaxNumRowEntries());

        for (unsigned i=0; i<localSubdomainMap->getNodeNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<GO> indicesGlobal;
                Array<SC> valuesLocal;
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
        return localSubdomainMatrix.getConst();
    }

    template <class SC,class LO,class GO,class NO>
    int UpdateLocalSubdomainMatrix(RCP<Matrix<SC,LO,GO,NO> > globalMatrix,
                                   RCP<Map<LO,GO,NO> > &map,
                                   RCP<Matrix<SC,LO,GO,NO> > &localSubdomainMatrix)
    {
        FROSCH_DETAILTIMER_START(updateLocalSubdomainMatrixTime,"UpdateLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);

        localSubdomainMatrix->resumeFill();
        for (unsigned i=0; i<map->getNodeNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<LO> indicesLocal;
                Array<SC> valuesLocal;
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
    int BuildSubmatrices(RCP<const Matrix<SC,LO,GO,NO> > k,
                         ArrayView<GO> indI,
                         RCP<Matrix<SC,LO,GO,NO> > &kII,
                         RCP<Matrix<SC,LO,GO,NO> > &kIJ,
                         RCP<Matrix<SC,LO,GO,NO> > &kJI,
                         RCP<Matrix<SC,LO,GO,NO> > &kJJ)
    {
        FROSCH_DETAILTIMER_START(buildSubmatricesTime,"BuildSubmatrices");
        // We need four Maps
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI(),0,k->getRowMap()->getComm());
        RCP<Map<LO,GO,NO> > mapILocal = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI.size(),0,k->getRowMap()->getComm());

        Array<GO> indJ;
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            if (mapI->getLocalElement(i)<0) {
                indJ.push_back(i);
            }
        }

        RCP<Map<LO,GO,NO> > mapJ = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indJ(),0,k->getRowMap()->getComm());
        RCP<Map<LO,GO,NO> > mapJLocal = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indJ.size(),0,k->getRowMap()->getComm());
        RCP<const Map<LO,GO,NO> > colMap = k->getColMap();
        kII = MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        kIJ = MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));
        kJI = MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        kJJ = MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));

        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            ArrayView<const LO> indices;
            ArrayView<const SC> values;

            k->getLocalRowView(i,indices,values);
            //cout << numEntries << endl;
            Array<GO> indicesI;
            Array<SC> valuesI;
            Array<GO> indicesJ;
            Array<SC> valuesJ;
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
                //cout << k->getRowMap().Comm().getRank() << " " << tmp1 << " numEntries " << numEntries << " indicesI.size() " << indicesI.size() << " indicesJ.size() " << indicesJ.size() << endl;
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
    int BuildSubmatrix(RCP<Matrix<SC,LO,GO,NO> > k,
                       ArrayView<GO> indI,
                       RCP<Matrix<SC,LO,GO,NO> > &kII)
    {
        FROSCH_DETAILTIMER_START(buildSubmatrixTime,"BuildSubmatrix");
        //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI(),0,k->getRowMap()->getComm());

        kII = MatrixFactory<SC,LO,GO,NO>::Build(mapI,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        GO maxGID = mapI->getMaxAllGlobalIndex();
        GO minGID = mapI->getMinAllGlobalIndex();
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            ArrayView<const LO> indices;
            ArrayView<const SC> values;

            k->getLocalRowView(i,indices,values);

            Array<GO> indicesI;
            Array<SC> valuesI;

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
                //cout << k->getRowMap().Comm().getRank() << " " << tmp1 << " numEntries " << numEntries << " indicesI.size() " << indicesI.size() << " indicesJ.size() " << indicesJ.size() << endl;
                kII->insertGlobalValues(mapI->getGlobalElement(tmp1),indicesI(),valuesI());
            }
        }
        kII->fillComplete(mapI,mapI);

        return 0;
    }

    template <class LO,class GO,class NO>
    int BuildSubgraph(RCP<CrsGraph<LO,GO,NO> > k,
                      ArrayView<GO> indI,
                      RCP<CrsGraph<LO,GO,NO> > &kII)
    {
        FROSCH_DETAILTIMER_START(buildSubgraphTime,"BuildSubgraph");
        //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),-1,indI(),0,k->getRowMap()->getComm());

        kII = CrsGraphFactory<LO,GO,NO>::Build(mapI,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        GO maxGID = mapI->getMaxAllGlobalIndex();
        GO minGID = mapI->getMinAllGlobalIndex();
        for (unsigned i=0; i<k->getNodeNumRows(); i++) {
            ArrayView<const LO> indices;

            k->getLocalRowView(i,indices);

            Array<GO> indicesI;

            LO tmp1=mapI->getLocalElement(k->getRowMap()->getGlobalElement(i));
            GO tmp2=0;
            if (tmp1>=0) {
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = k->getColMap()->getGlobalElement(indices[j]);
                    if (minGID<=tmp2 && tmp2<=maxGID) {
                        indicesI.push_back(tmp2);
                    }
                }
                //cout << k->getRowMap().Comm().getRank() << " " << tmp1 << " numEntries " << numEntries << " indicesI.size() " << indicesI.size() << " indicesJ.size() " << indicesJ.size() << endl;
                kII->insertGlobalValues(mapI->getGlobalElement(tmp1),indicesI());
            }
        }
        kII->fillComplete(mapI,mapI);

        return 0;
    }
}

#endif
