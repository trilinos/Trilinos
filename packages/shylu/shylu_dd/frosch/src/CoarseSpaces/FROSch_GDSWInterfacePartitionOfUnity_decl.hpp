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

#ifndef _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <FROSch_InterfacePartitionOfUnity_def.hpp>

namespace FROSch {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC, LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC, LO, GO>::node_type>
    class GDSWInterfacePartitionOfUnity : public InterfacePartitionOfUnity<SC,LO,GO,NO> {
        
    public:
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::CommPtr CommPtr;

        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::Map Map;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::MapPtr MapPtr;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::MapPtrVecPtr MapPtrVecPtr;
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::CrsMatrix CrsMatrix;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::CrsMatrixPtr CrsMatrixPtr;
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVector MultiVector;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMultiVectorPtr ConstMultiVectorPtr;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVectorPtr MultiVectorPtr;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVectorPtrVecPtr MultiVectorPtrVecPtr;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMultiVectorPtrVecPtr ConstMultiVectorPtrVecPtr;

        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::ParameterListPtr ParameterListPtr;
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr DDInterfacePtr;
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr EntitySetPtr;
        
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::UN UN;

        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::GOVec GOVec;
        typedef typename GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::GOVecView GOVecView;
        
        
        GDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                      CommPtr serialComm,
                                      UN dimension,
                                      UN dofsPerNode,
                                      MapPtr nodesMap,
                                      MapPtrVecPtr dofsMaps,
                                      ParameterListPtr parameterList);
        
        virtual ~GDSWInterfacePartitionOfUnity();
        
        virtual int removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                         MultiVectorPtr nodeList);
        
        virtual int sortInterface(CrsMatrixPtr matrix,
                                  MultiVectorPtr nodeList);
        
        virtual int computePartitionOfUnity();
        
    protected:
        
        bool UseVertices_;
        bool UseShortEdges_;
        bool UseStraightEdges_;
        bool UseEdges_;
        bool UseFaces_;
        
        EntitySetPtr Vertices_;
        EntitySetPtr ShortEdges_;
        EntitySetPtr StraightEdges_;
        EntitySetPtr Edges_;
        EntitySetPtr Faces_;
        
    };
    
}

#endif
