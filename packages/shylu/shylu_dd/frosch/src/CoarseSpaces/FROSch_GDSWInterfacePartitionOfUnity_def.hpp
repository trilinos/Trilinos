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

#ifndef _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DEF_HPP

#include <FROSch_GDSWInterfacePartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::GDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                                                              CommPtr serialComm,
                                                                              UN dimension,
                                                                              UN dofsPerNode,
                                                                              ConstXMapPtr nodesMap,
                                                                              ConstXMapPtrVecPtr dofsMaps,
                                                                              ParameterListPtr parameterList,
                                                                              Verbosity verbosity,
                                                                              UN levelID) :
    InterfacePartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dimension,dofsPerNode,nodesMap,dofsMaps,parameterList,verbosity,levelID),
    UseVertices_ (false),
    UseShortEdges_ (false),
    UseStraightEdges_ (false),
    UseEdges_ (false),
    UseFaces_ (false),
    Vertices_ (),
    ShortEdges_ (),
    StraightEdges_ (),
    Edges_ (),
    Faces_ ()
    {
        FROSCH_TIMER_START_LEVELID(gDSWInterfacePartitionOfUnityTime,"GDSWInterfacePartitionOfUnity::GDSWInterfacePartitionOfUnity");
        if (!this->ParameterList_->get("Type","Full").compare("Full")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Vertices")) {
            UseVertices_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("ShortEdges")) {
            UseShortEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("StraightEdges")) {
            UseStraightEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Edges")) {
            UseEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("AllEdges")) {
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoVertices")) {
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoShortEdges")) {
            UseVertices_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoStraightEdges")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoEdges")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoAllEdges")) {
            UseVertices_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoFaces")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Custom")) {
            UseVertices_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseShortEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseStraightEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseFaces_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
        } else {
            FROSCH_ASSERT(false,"FROSch::GDSWInterfacePartitionOfUnity : ERROR: Specify a valid Type.");
        }
        this->LocalPartitionOfUnity_ = XMultiVectorPtrVecPtr(5);
        this->PartitionOfUnityMaps_ = XMapPtrVecPtr(5);
    }

    template <class SC,class LO,class GO,class NO>
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::~GDSWInterfacePartitionOfUnity()
    {

    }

    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                                                         ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(removeDirichletNodesTime,"GDSWInterfacePartitionOfUnity::removeDirichletNodes");
        if (!dirichletBoundaryDofs.is_null()) {
            GOVec tmpDirichletBoundaryDofs(dirichletBoundaryDofs());
            sortunique(tmpDirichletBoundaryDofs);
            this->DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());
            this->DDInterface_->sortVerticesEdgesFaces(nodeList);
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::sortInterface(ConstXMatrixPtr matrix,
                                                                  ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(sortInterfaceTime,"GDSWInterfacePartitionOfUnity::sortInterface");
        if (this->ParameterList_->get("Test Unconnected Interface",true)) {
            if (matrix.is_null()) {
                if (this->Verbose_) std::cout << "FROSch::GDSWInterfacePartitionOfUnity : WARNING: divideUnconnectedEntities() cannot be performed without the matrix." << std::endl;
            } else this->DDInterface_->divideUnconnectedEntities(matrix);
        }
        this->DDInterface_->sortVerticesEdgesFaces(nodeList);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::computePartitionOfUnity(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(computePartitionOfUnityTime,"GDSWInterfacePartitionOfUnity::computePartitionOfUnity");
        // Interface
        UN dofsPerNode = this->DDInterface_->getInterface()->getEntity(0)->getDofsPerNode();
        UN numInterfaceDofs = dofsPerNode*this->DDInterface_->getInterface()->getEntity(0)->getNumNodes();

        this->DDInterface_->buildEntityMaps(UseVertices_,
                                            UseShortEdges_,
                                            UseStraightEdges_,
                                            UseEdges_,
                                            UseFaces_,
                                            false,
                                            false);

        // Maps
        if (UseVertices_) {
            Vertices_ = this->DDInterface_->getVertices();
            this->PartitionOfUnityMaps_[0] = Vertices_->getEntityMap();
        }

        if (UseShortEdges_) {
            ShortEdges_ = this->DDInterface_->getShortEdges();
            this->PartitionOfUnityMaps_[1] = ShortEdges_->getEntityMap();
        }

        if (UseStraightEdges_) {
            StraightEdges_ = this->DDInterface_->getStraightEdges();
            this->PartitionOfUnityMaps_[2] = StraightEdges_->getEntityMap();
        }

        if (UseEdges_) {
            Edges_ = this->DDInterface_->getEdges();
            this->PartitionOfUnityMaps_[3] = Edges_->getEntityMap();
        }

        if (UseFaces_) {
            Faces_ = this->DDInterface_->getFaces();
            this->PartitionOfUnityMaps_[4] = Faces_->getEntityMap();
        }

        if (this->Verbose_) {
            std::cout << std::boolalpha << "\n\
    ------------------------------------------------------------------------------\n\
     GDSW Interface Partition Of Unity (GDSW IPOU)\n\
    ------------------------------------------------------------------------------\n\
      Vertices                                    --- " << UseVertices_ << "\n\
      ShortEdges                                  --- " << UseShortEdges_ << "\n\
      StraightEdges                               --- " << UseStraightEdges_ << "\n\
      Edges                                       --- " << UseEdges_ << "\n\
      Faces                                       --- " << UseFaces_ << "\n\
    ------------------------------------------------------------------------------\n" << std::noboolalpha;
        }

        // Build Partition Of Unity Vectors
        XMapPtr serialInterfaceMap = MapFactory<LO,GO,NO>::Build(this->DDInterface_->getNodesMap()->lib(),numInterfaceDofs,0,this->SerialComm_);

        if (UseVertices_ && Vertices_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Vertices_->getNumEntities());

            for (UN i=0; i<Vertices_->getNumEntities(); i++) {
                for (UN j=0; j<Vertices_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Vertices_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }

            this->LocalPartitionOfUnity_[0] = tmpVector;
        }

        if (UseShortEdges_ && ShortEdges_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,ShortEdges_->getNumEntities());

            for (UN i=0; i<ShortEdges_->getNumEntities(); i++) {
                for (UN j=0; j<ShortEdges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(ShortEdges_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }

            this->LocalPartitionOfUnity_[1] = tmpVector;
        }

        if (UseStraightEdges_ && StraightEdges_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,StraightEdges_->getNumEntities());

            for (UN i=0; i<StraightEdges_->getNumEntities(); i++) {
                for (UN j=0; j<StraightEdges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(StraightEdges_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }
            this->LocalPartitionOfUnity_[2] = tmpVector;
        }

        if (UseEdges_ && Edges_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Edges_->getNumEntities());

            for (UN i=0; i<Edges_->getNumEntities(); i++) {
                for (UN j=0; j<Edges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Edges_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }

            this->LocalPartitionOfUnity_[3] = tmpVector;
        }

        if (UseFaces_ && Faces_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Faces_->getNumEntities());

            for (UN i=0; i<Faces_->getNumEntities(); i++) {
                for (UN j=0; j<Faces_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Faces_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }

            this->LocalPartitionOfUnity_[4] = tmpVector;
        }

        return 0;
    }

}

#endif
