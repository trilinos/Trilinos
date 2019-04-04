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

#ifndef _FROSCH_PARTITIONOFUNITYBASIS_DECL_hpp
#define _FROSCH_PARTITIONOFUNITYBASIS_DECL_hpp

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <Xpetra_Operator.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include <FROSch_CoarseSpace_def.hpp>

#include "FROSch_Tools_def.hpp"

namespace FROSch {

    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC, LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC, LO, GO>::node_type>
    class LocalPartitionOfUnityBasis {
        
    public:

        typedef Teuchos::RCP<const Teuchos::Comm<int> > CommPtr;
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        typedef Teuchos::ArrayRCP<MultiVectorPtr> MultiVectorPtrVecPtr;
        typedef Teuchos::ArrayRCP<MultiVectorPtrVecPtr> MultiVectorPtrVecPtr2D;                
        
        typedef Teuchos::RCP<Teuchos::ParameterList> ParameterListPtr;
        
        typedef Teuchos::RCP<CoarseSpace<SC,LO,GO,NO> > CoarseSpacePtr;
        
        typedef unsigned UN;
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        typedef Teuchos::ArrayRCP<LOVecPtr> LOVecPtr2D;
        
        typedef Teuchos::ArrayRCP<bool> BoolVecPtr;
        typedef Teuchos::ArrayRCP<BoolVecPtr> BoolVecPtr2D;

        
        LocalPartitionOfUnityBasis(CommPtr mpiComm,
                                   CommPtr serialComm,
                                   UN dofsPerNode,
                                   ParameterListPtr parameterList,
                                   MultiVectorPtr nullSpaceBasis = MultiVectorPtr(),
                                   MultiVectorPtrVecPtr partitionOfUnity = MultiVectorPtrVecPtr(),
                                   MapPtrVecPtr partitionOfUnityMaps = MapPtrVecPtr());
        
//        virtual ~LocalPartitionOfUnityBasis();
        
        int addPartitionOfUnity(MultiVectorPtrVecPtr partitionOfUnity,
                                MapPtrVecPtr partitionOfUnityMaps);
        
        int addGlobalBasis(MultiVectorPtr nullSpaceBasis);
        
        int buildLocalPartitionOfUnityBasis();
        
        MultiVectorPtrVecPtr getPartitionOfUnity() const;
        
        MultiVectorPtr getNullspaceBasis() const;
        
        CoarseSpacePtr getLocalPartitionOfUnitySpace() const;
        
    protected:
        
        CommPtr MpiComm_;
        CommPtr SerialComm_;
        
        UN DofsPerNode_;
        
        ParameterListPtr ParameterList_;
        
        CoarseSpacePtr LocalPartitionOfUnitySpace_;
        
        MultiVectorPtrVecPtr PartitionOfUnity_;
        MultiVectorPtr NullspaceBasis_;
                
        MapPtrVecPtr PartitionOfUnityMaps_;
        
    };
    
}

#endif
