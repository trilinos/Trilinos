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

//#ifdef HAVE_MPI
#include "mpi.h"

#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "FROSch_LocalPartitionOfUnityBasis_def.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
//using namespace Belos;

int main(int argc, char *argv[])
{
    
    MPI_Init(&argc,&argv);
    
    {
        
        RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
        
        Array<GO> nodes(9);
        Array<GO> dofs(27);
        Array<GO> dofs1(9);
        Array<GO> dofs2(9);
        Array<GO> dofs3(9);
        
        for (UN i=0; i<9; i++) {
            nodes[i] = i;
            dofs[3*i] = 3*i;
            dofs[3*i+1] = 3*i+1;
            dofs[3*i+2] = 3*i+2;
            dofs1[i] = 3*i;
            dofs2[i] = 3*i+1;
            dofs3[i] = 3*i+2;
        }
        
        RCP<Map<LO,GO,NO> > nodesMap = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,nodes(),0,SerialComm);
        RCP<Map<LO,GO,NO> > dofsMap = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,dofs(),0,SerialComm);
        ArrayRCP<RCP<Map<LO,GO,NO> > > dofsMaps (3);
        dofsMaps[0] = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,dofs1(),0,SerialComm);
        dofsMaps[1] = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,dofs2(),0,SerialComm);
        dofsMaps[2] = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,dofs3(),0,SerialComm);
        
        ArrayRCP<RCP<MultiVector<SC,LO,GO,NO> > > partitionOfUnity(2);
        
        RCP<MultiVector<SC,LO,GO,NO> > tmpVec = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,2);
        tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[0]),0,1.0);
        tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[0]),0,1.0);
        tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[0]),0,1.0);
        
        tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[8]),1,1.0);
        tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[8]),1,1.0);
        tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[8]),1,1.0);
        partitionOfUnity[0] = tmpVec;
        
        tmpVec = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,1);
        for (UN i=1; i<8; i++) {
            tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),0,1.0);
            tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),0,1.0);
            tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),0,1.0);
        }
        partitionOfUnity[1] = tmpVec;
        
        RCP<MultiVector<SC,LO,GO,NO> > nullspace = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,6);
        for (UN i=0; i<9; i++) {
            nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),0,1.0);
            nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),1,1.0);
            nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),2,1.0);
            
            nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),3,double(i));
            nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),3,-double(i));
            nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),3,0.0);
            
            nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),4,-double(i));
            nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),4,0.0);
            nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),4,double(i));
            
            nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),5,0.0);
            nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),5,double(i));
            nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),5,-double(i));
        }
        

        Array<GO> partitionOfUnityMapVec1(2);
        partitionOfUnityMapVec1[0] = 0;
        partitionOfUnityMapVec1[1] = 1;
        Array<GO> partitionOfUnityMapVec2(1);
        partitionOfUnityMapVec1[0] = 0;
        
        ArrayRCP<RCP<Map<LO,GO,NO> > > partitionOfUnityMaps(2);
        partitionOfUnityMaps[0] = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,partitionOfUnityMapVec1(),0,SerialComm);
        partitionOfUnityMaps[1] = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,partitionOfUnityMapVec2(),0,SerialComm);
        
        RCP<ParameterList> parameterList = getParametersFromXmlFile("Parameters.xml");
        
        LocalPartitionOfUnityBasis<SC,LO,GO,NO> TestBasis(SerialComm,SerialComm,3,parameterList,nullspace,partitionOfUnity,partitionOfUnityMaps);
        TestBasis.buildLocalPartitionOfUnityBasis();
        
        RCP<CoarseSpace<SC,LO,GO,NO> > coarseSpace = TestBasis.getLocalPartitionOfUnitySpace();
        
        for (UN i=0; i<coarseSpace->getLocalBasis()->getLocalLength(); i++) {
            for (UN j=0; j<coarseSpace->getLocalBasis()->getNumVectors(); j++) {
                std::cout << coarseSpace->getLocalBasis()->getData(j)[i] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";        
        //        RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));
        //        partitionOfUnity->describe(*fancy,Teuchos::VERB_EXTREME);
        //        globalBasis->describe(*fancy,Teuchos::VERB_EXTREME);

    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);
}
