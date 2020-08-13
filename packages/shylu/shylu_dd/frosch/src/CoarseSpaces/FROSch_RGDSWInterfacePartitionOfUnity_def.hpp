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

#ifndef _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DEF_HPP

#include <FROSch_RGDSWInterfacePartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RGDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::RGDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                                                                CommPtr serialComm,
                                                                                UN dimension,
                                                                                UN dofsPerNode,
                                                                                ConstXMapPtr nodesMap,
                                                                                ConstXMapPtrVecPtr dofsMaps,
                                                                                ParameterListPtr parameterList,
                                                                                Verbosity verbosity,
                                                                                UN levelID) :
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dimension,dofsPerNode,nodesMap,dofsMaps,parameterList,verbosity,levelID)
    {
        FROSCH_DETAILTIMER_START_LEVELID(rGDSWInterfacePartitionOfUnityTime,"RGDSWInterfacePartitionOfUnity::RGDSWInterfacePartitionOfUnity");
        this->UseVertices_ = false;
        this->UseShortEdges_ = false;
        this->UseStraightEdges_ = false;
        this->UseEdges_ = false;
        this->UseFaces_ = false;

        if (!this->ParameterList_->get("Type","Full").compare("Full")) {
            UseRoots_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Roots")) {
            UseRoots_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Custom")) {
            UseRoots_ = this->ParameterList_->sublist("Custom").get("Roots",false);
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWInterfacePartitionOfUnity : ERROR: Specify a valid Type.");
        }

        if (!this->ParameterList_->get("Distance Function","Constant").compare("Constant")) {
            DistanceFunction_ = ConstantDistanceFunction;
        } else if (!this->ParameterList_->get("Distance Function","Constant").compare("Inverse Euclidean")) {
            DistanceFunction_ = InverseEuclideanDistanceFunction;
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWInterfacePartitionOfUnity : ERROR: Specify a valid Distance Function.");
        }
        this->LocalPartitionOfUnity_ = ConstXMultiVectorPtrVecPtr(1);
        this->PartitionOfUnityMaps_ = ConstXMapPtrVecPtr(1);
    }

    template <class SC,class LO,class GO,class NO>
    int RGDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::computePartitionOfUnity(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computePartitionOfUnityTime,"RGDSWInterfacePartitionOfUnity::computePartitionOfUnity");
        // Interface
        UN dofsPerNode = this->DDInterface_->getInterface()->getEntity(0)->getDofsPerNode();
        UN numInterfaceDofs = dofsPerNode*this->DDInterface_->getInterface()->getEntity(0)->getNumNodes();

        this->DDInterface_->buildEntityHierarchy();

        this->DDInterface_->computeDistancesToRoots(this->DDInterface_->getDimension(),nodeList,DistanceFunction_);

        this->DDInterface_->buildEntityMaps(false,
                                            false,
                                            false,
                                            false,
                                            false,
                                            UseRoots_,
                                            false);

        EntitySetVector_ = this->DDInterface_->getEntitySetVector();

        // Map
        if (UseRoots_) {
            Roots_ = this->DDInterface_->getRoots();
            this->PartitionOfUnityMaps_[0] = Roots_->getEntityMap();
        }

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "> RGDSW Interface Partition Of Unity " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")" << right
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "========================================================================================="
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Roots" << right
            << " | " << setw(41) << boolalpha << UseRoots_ << noboolalpha
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << endl;
        }

        // Build Partition Of Unity Vectors
        XMapPtr serialInterfaceMap = MapFactory<LO,GO,NO>::Build(this->DDInterface_->getNodesMap()->lib(),numInterfaceDofs,0,this->SerialComm_);

        if (UseRoots_ && Roots_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Roots_->getNumEntities());

            // Loop over EntitySetVector_
            for (UN i=0; i<EntitySetVector_.size(); i++) {
                // Loop over entities
                for (UN j=0; j<EntitySetVector_[i]->getNumEntities(); j++) {
                    InterfaceEntityPtr tmpEntity = EntitySetVector_[i]->getEntity(j);
                    LO rootID = tmpEntity->getRootID();
                    UN numRoots = tmpEntity->getRoots()->getNumEntities();
                    if (rootID==-1) {
                        FROSCH_ASSERT(numRoots!=0,"rootID==-1 but numRoots==0!");
                        for (UN m=0; m<numRoots; m++) {
                            InterfaceEntityPtr tmpRoot = tmpEntity->getRoots()->getEntity(m);
                            LO index = tmpRoot->getRootID();
                            // Offspring: loop over nodes
                            for (UN l=0; l<tmpEntity->getNumNodes(); l++) {
                                SC value = tmpEntity->getDistanceToRoot(l,m)/tmpEntity->getDistanceToRoot(l,numRoots);
                                for (UN k=0; k<dofsPerNode; k++) {
                                    tmpVector->replaceLocalValue(tmpEntity->getGammaDofID(l,k),index,value*ScalarTraits<SC>::one());
                                }
                            }
                        }
                    } else {
                        // Coarse node: loop over nodes
                        for (UN l=0; l<EntitySetVector_[i]->getEntity(j)->getNumNodes(); l++) {
                            for (UN k=0; k<dofsPerNode; k++) {
                                tmpVector->replaceLocalValue(tmpEntity->getGammaDofID(l,k),rootID,ScalarTraits<SC>::one());
                            }
                        }
                    }
                }
            }
            this->LocalPartitionOfUnity_[0] = tmpVector;
        }

        return 0;
    }

}

#endif
