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

#ifndef _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP
#define _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP

#include <FROSch_CoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class HarmonicCoarseOperator : public CoarseOperator<SC,LO,GO,NO> {

    protected:

        using XMapPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr            = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr           = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using XMapPtrVecPtr2D         = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D    = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr              = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr         = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr    = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr   = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;

        using ParameterListPtr        = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using CoarseSpacePtr          = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtr;
        using CoarseSpacePtrVecPtr    = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtrVecPtr;

        using EntitySetPtr            = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetConstPtr       = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetConstPtr;

        using SubdomainSolverPtr      = typename SchwarzOperator<SC,LO,GO,NO>::SubdomainSolverPtr;

        using UN                      = typename SchwarzOperator<SC,LO,GO,NO>::UN;
        using UNVec                   = typename SchwarzOperator<SC,LO,GO,NO>::UNVec;
        using UNVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr;
        using ConstUNVecView          = typename SchwarzOperator<SC,LO,GO,NO>::ConstUNVecView;

        using LOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using LOVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecView               = typename SchwarzOperator<SC,LO,GO,NO>::GOVecView;
        using GOVec2D                 = typename SchwarzOperator<SC,LO,GO,NO>::GOVec2D;
        
        using SCVec                   = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;
        using SCVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;

    public:

        HarmonicCoarseOperator(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        virtual int initialize() = 0;

        XMapPtr computeCoarseSpace(CoarseSpacePtr coarseSpace);

    protected:

        int intializeCoarseMap();
        
        int assembleInterfaceCoarseSpace();

        int addZeroCoarseSpaceBlock(ConstXMapPtr dofsMap);

        int computeVolumeFunctions(UN blockId,
                                   UN dimension,
                                   ConstXMapPtr nodesMap,
                                   ConstXMultiVectorPtr nodeList,
                                   EntitySetConstPtr interior);

        virtual XMultiVectorPtrVecPtr computeTranslations(UN blockId,
                                                          EntitySetConstPtr entitySet);

        virtual XMultiVectorPtrVecPtr computeRotations(UN blockId,
                                                       UN dimension,
                                                       ConstXMultiVectorPtr nodeList,
                                                       EntitySetConstPtr entitySet,
                                                       UN discardRotations = 0);

        virtual XMultiVectorPtr computeExtensions(ConstXMapPtr localMap,
                                                  ConstXMapPtr coarseMap,
                                                  GOVecView indicesGammaDofsAll,
                                                  GOVecView indicesIDofsAll,
                                                  XMatrixPtr kII,
                                                  XMatrixPtr kIGamma);


        SubdomainSolverPtr ExtensionSolver_;

        CoarseSpacePtrVecPtr InterfaceCoarseSpaces_ = CoarseSpacePtrVecPtr(0);
        CoarseSpacePtr AssembledInterfaceCoarseSpace_;

        UNVecPtr Dimensions_ = UNVecPtr(0);
        UNVecPtr DofsPerNode_ = UNVecPtr(0);

        LOVecPtr2D GammaDofs_ = GammaDofs_(0);
        LOVecPtr2D IDofs_ = IDofs_(0);

        ConstXMapPtrVecPtr2D DofsMaps_ = DofsMaps_(0); // notwendig??

        UN NumberOfBlocks_ = 0;
    };

}

#endif
