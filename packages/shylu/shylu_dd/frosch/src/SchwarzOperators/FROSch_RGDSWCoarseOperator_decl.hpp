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

#ifndef _FROSCH_RGDSWCOARSEOPERATOR_DECL_HPP
#define _FROSCH_RGDSWCOARSEOPERATOR_DECL_HPP

#include <FROSch_GDSWCoarseOperator_def.hpp>

namespace FROSch {

    template <class SC = Xpetra::Operator<>::scalar_type,
              class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
              class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
              class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class RGDSWCoarseOperator : public GDSWCoarseOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using MapPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::MapPtr;
        using ConstMapPtr             = typename SchwarzOperator<SC,LO,GO,NO>::ConstMapPtr;
        using MapPtrVecPtr            = typename SchwarzOperator<SC,LO,GO,NO>::MapPtrVecPtr;
        using ConstMapPtrVecPtr       = typename SchwarzOperator<SC,LO,GO,NO>::ConstMapPtrVecPtr;

        using CrsMatrixPtr            = typename SchwarzOperator<SC,LO,GO,NO>::CrsMatrixPtr;
        using ConstCrsMatrixPtr       = typename SchwarzOperator<SC,LO,GO,NO>::ConstCrsMatrixPtr;

        using MultiVectorPtr          = typename SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr;
        using ConstMultiVectorPtr     = typename SchwarzOperator<SC,LO,GO,NO>::ConstMultiVectorPtr;
        using MultiVectorPtrVecPtr    = typename SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr;

        using ParameterListPtr        = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr          = typename SchwarzOperator<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr            = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetPtrVecPtr      = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrVecPtr;

        using InterfaceEntityPtr      = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtr;

        using UN                      = typename SchwarzOperator<SC,LO,GO,NO>::UN;

        using LOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using LOVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr;
        using GOVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr2D;

        using SCVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;

    public:

        RGDSWCoarseOperator(ConstCrsMatrixPtr k,
                            ParameterListPtr parameterList);

        virtual int resetCoarseSpaceBlock(UN blockId,
                                          UN dimension,
                                          UN dofsPerNode,
                                          ConstMapPtr nodesMap,
                                          ConstMapPtrVecPtr dofsMaps,
                                          GOVecPtr dirichletBoundaryDofs,
                                          ConstMultiVectorPtr nodeList);


    protected:

        virtual MultiVectorPtrVecPtr computeTranslations(UN blockId,
                                                         EntitySetPtr coarseNodes,
                                                         EntitySetPtrVecPtr entitySetVector,
                                                         DistanceFunction distanceFunction = ConstantDistanceFunction);

        virtual MultiVectorPtrVecPtr computeRotations(UN blockId,
                                                      UN dimension,
                                                      ConstMultiVectorPtr nodeList,
                                                      EntitySetPtr coarseNodes,
                                                      EntitySetPtrVecPtr entitySetVector,
                                                      DistanceFunction distanceFunction = ConstantDistanceFunction);

    };

}

#endif
