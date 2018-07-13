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

#ifndef _FROSCH_GDSWCOARSEOPERATOR_DECL_HPP
#define _FROSCH_GDSWCOARSEOPERATOR_DECL_HPP

#include <FROSch_HarmonicCoarseOperator_def.hpp>


namespace FROSch {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class GDSWCoarseOperator : public HarmonicCoarseOperator<SC,LO,GO,NO> {
        
    public:
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::CommPtr CommPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::MapPtr MapPtr;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::MapPtrVecPtr MapPtrVecPtr;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::MapPtrVecPtr2D MapPtrVecPtr2D;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::CrsMatrixPtr CrsMatrixPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr MultiVectorPtr;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr MultiVectorPtrVecPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr ParameterListPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::DDInterfacePtr DDInterfacePtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr EntitySetPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::SubdomainSolverPtr SubdomainSolverPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::UN UN;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr UNVecPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::LOVec LOVec;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr LOVecPtr;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D LOVecPtr2D;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::GOVec GOVec;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr GOVecPtr;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::GOVecView GOVecView;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr2D GOVecPtr2D;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::SCVec SCVec;
        typedef typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr SCVecPtr;
        
        typedef typename SchwarzOperator<SC,LO,GO,NO>::BoolVecPtr BoolVecPtr;
        
        
        GDSWCoarseOperator(CrsMatrixPtr k,
                           ParameterListPtr parameterList);
        
        virtual int initialize()
        {
            FROSCH_ASSERT(0!=0,"GDSWCoarseOperator cannot be built without a repeated Map");
            return 0;
        };
        
        int initialize(UN dimension,
                       MapPtr repeatedMap);
        
        int initialize(UN dimension,
                       MapPtr repeatedMap,
                       GOVecPtr dirichletBoundaryDofs);
        
        int initialize(UN dimension,
                       UN dofsPerNode,
                       MapPtr repeatedNodesMap,
                       MapPtrVecPtr RepeatedDofMaps);
        
        int initialize(UN dimension,
                       UN dofsPerNode,
                       MapPtr repeatedNodesMap,
                       MapPtrVecPtr RepeatedDofMaps,
                       GOVecPtr dirichletBoundaryDofs);
        
        int initialize(UN dimension,
                       UN dofsPerNode,
                       MapPtr repeatedNodesMap,
                       MapPtrVecPtr RepeatedDofMaps,
                       MultiVectorPtr nodeList);
        
        int initialize(UN dimension,
                       UN dofsPerNode,
                       MapPtr repeatedNodesMap,
                       MapPtrVecPtr RepeatedDofMaps,
                       GOVecPtr dirichletBoundaryDofs,
                       MultiVectorPtr nodeList);
        
        void describe(Teuchos::FancyOStream &out,
                      const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
        
        std::string description() const;
        
    protected:
        
        int buildCoarseSpace(UN dimension,
                             MapPtr nodesMap);
        
        int buildCoarseSpace(UN dimension,
                             MapPtr nodesMap,
                             GOVecPtr dirichletBoundaryDofs); // Das kann man auch mit in den Fall davor reinnehmen ?!
        
        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             MapPtr nodesMap,
                             MapPtrVecPtr dofsMaps);
        
        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             MapPtr nodesMap,
                             MapPtrVecPtr dofsMaps,
                             GOVecPtr dirichletBoundaryDofs);
        
        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             MapPtr nodesMap,
                             MapPtrVecPtr dofsMaps,
                             MultiVectorPtr nodeList);
        
        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             MapPtr nodesMap,
                             MapPtrVecPtr dofsMaps,
                             GOVecPtr dirichletBoundaryDofs,
                             MultiVectorPtr nodeList);
        
        virtual int resetCoarseSpaceBlock(UN blockId,
                                          UN dimension,
                                          UN dofsPerNode,
                                          MapPtr nodesMap,
                                          MapPtrVecPtr dofsMaps,
                                          GOVecPtr dirichletBoundaryDofs,
                                          MultiVectorPtr nodeList);
        
        int phiGammaGDSW(UN blockId,
                         bool buildRotations,
                         UN dimension,
                         UN dofsPerNode,
                         MultiVectorPtr nodeList,
                         LOVecPtr2D partMappings,
                         EntitySetPtr vertices,
                         EntitySetPtr shortEdges,
                         EntitySetPtr straightEdges,
                         EntitySetPtr edges,
                         EntitySetPtr faces,
                         BoolVecPtr coarseSpaceFunctions);
        
        
        DDInterfacePtr DDInterface_;
        
    };
    
}

#endif
