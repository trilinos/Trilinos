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

#ifndef _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <FROSch_GDSWInterfacePartitionOfUnity_def.hpp>


namespace FROSch {

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class RGDSWInterfacePartitionOfUnity : public GDSWInterfacePartitionOfUnity<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::CommPtr;

        using Map                           = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::Map ;
        using MapPtr                        = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::MapPtr;
        using ConstMapPtr                   = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMapPtr;
        using MapPtrVecPtr                  = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::MapPtrVecPtr;
        using ConstMapPtrVecPtr             = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMapPtrVecPtr;

        using CrsMatrix                     = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::CrsMatrix;
        using CrsMatrixPtr                  = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::CrsMatrixPtr;
        using ConstCrsMatrixPtr             = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstCrsMatrixPtr;

        using MultiVector                   = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVector;
        using ConstMultiVectorPtr           = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMultiVectorPtr;
        using MultiVectorPtr                = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVectorPtr;
        using MultiVectorPtrVecPtr          = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::MultiVectorPtrVecPtr;
        using ConstMultiVectorPtrVecPtr     = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstMultiVectorPtrVecPtr;

        using ParameterListPtr              = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr;
        
        using EntitySetPtr                  = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetPtrVecPtr            = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::EntitySetPtrVecPtr;
        
        using InterfaceEntityPtr            = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::InterfaceEntityPtr;

        using UN                            = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::UN;

        using GOVec                         = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::GOVec;
        using GOVecView                     = typename InterfacePartitionOfUnity<SC,LO,GO,NO>::GOVecView;

    public:

        RGDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                       CommPtr serialComm,
                                       UN dimension,
                                       UN dofsPerNode,
                                       ConstMapPtr nodesMap,
                                       ConstMapPtrVecPtr dofsMaps,
                                       ParameterListPtr parameterList,
                                       Verbosity verbosity = All);

        virtual int computePartitionOfUnity();

        virtual int computePartitionOfUnity(ConstMultiVectorPtr nodeList);

    protected:

        bool UseCoarseNodes_;

        EntitySetPtr CoarseNodes_;

        EntitySetPtrVecPtr EntitySetVector_;

        DistanceFunction DistanceFunction_;
    };

}

#endif
