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

#ifndef _FROSCH_COARSESPACE_DECL_HPP
#define _FROSCH_COARSESPACE_DECL_HPP

//#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include <FROSch_Tools_def.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class CoarseSpace {

    protected:

        using XMap                  = Map<LO,GO,NO>;
        using XMapPtr               = RCP<XMap>;
        using ConstXMapPtr          = RCP<const XMap>;
        using XMapPtrVec            = Array<XMapPtr>;

        using XMatrix               = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr            = RCP<XMatrix>;

        using XMultiVector          = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr       = RCP<XMultiVector>;
        using XMultiVectorPtrVec    = Array<XMultiVectorPtr>;

        using ParameterListPtr      = RCP<ParameterList>;

        using UN                    = unsigned;

        using LOVec                 = Array<LO>;
        using GOVec                 = Array<GO>;
        using LOVecPtr              = ArrayRCP<LO>;
        using LOVecPtr2D            = ArrayRCP<LOVecPtr>;

        using SCVec                 = Array<SC>;

    public:

        CoarseSpace();

        int addSubspace(XMapPtr subspaceBasisMap,
                        XMultiVectorPtr subspaceBasis = null);

        int assembleCoarseSpace();

        int buildGlobalBasisMatrix(ConstXMapPtr rowMap,
                                   ConstXMapPtr repeatedMap,
                                   SC treshold);

        int clearCoarseSpace();

        int checkForLinearDependencies();

        bool hasUnassembledMaps() const;
        
        bool hasBasisMap() const;

        XMapPtr getBasisMap() const;

        bool hasAssembledBasis() const;

        XMultiVectorPtr getAssembledBasis() const;

        bool hasGlobalBasisMatrix() const;

        XMatrixPtr getGlobalBasisMatrix() const;

    protected:

        ConstXMapPtr SerialRowMap_;

        XMapPtrVec UnassembledBasesMaps_;

        XMultiVectorPtrVec UnassembledSubspaceBases_;

        XMapPtr AssembledBasisMap_;

        XMultiVectorPtr AssembledBasis_;

        XMatrixPtr GlobalBasisMatrix_;

    };

}

#endif
