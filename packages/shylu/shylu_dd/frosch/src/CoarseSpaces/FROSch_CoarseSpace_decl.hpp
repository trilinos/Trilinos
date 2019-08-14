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

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

//#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include <FROSch_Tools_def.hpp>


namespace FROSch {

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class CoarseSpace {

    protected:

        using Map                   = Xpetra::Map<LO,GO,NO>;
        using MapPtr                = Teuchos::RCP<Map>;
        using ConstMapPtr           = Teuchos::RCP<const Map>;
        using MapPtrVec             = Teuchos::Array<MapPtr>;

        using CrsMatrix             = Xpetra::Matrix<SC,LO,GO,NO>;
        using CrsMatrixPtr          = Teuchos::RCP<CrsMatrix>;

        using MultiVector           = Xpetra::MultiVector<SC,LO,GO,NO>;
        using MultiVectorPtr        = Teuchos::RCP<MultiVector>;
        using MultiVectorPtrVec     = Teuchos::Array<MultiVectorPtr>;

        using ParameterListPtr      = Teuchos::RCP<Teuchos::ParameterList>;

        using UN                    = unsigned;

        using LOVec                 = Teuchos::Array<LO>;
        using GOVec                 = Teuchos::Array<GO>;
        using LOVecPtr              = Teuchos::ArrayRCP<LO>;
        using LOVecPtr2D            = Teuchos::ArrayRCP<LOVecPtr>;

        using SCVec                 = Teuchos::Array<SC>;

    public:

        CoarseSpace();

        int addSubspace(MapPtr subspaceBasisMap,
                        MultiVectorPtr subspaceBasis = Teuchos::null);

        int assembleCoarseSpace();

        int buildGlobalBasisMatrix(ConstMapPtr rowMap,
                                   ConstMapPtr repeatedMap,
                                   SC treshold);

        int clearCoarseSpace();

        int checkForLinearDependencies();

        bool hasBasisMap() const;

        MapPtr getBasisMap() const;

        bool hasAssembledBasis() const;

        MultiVectorPtr getAssembledBasis() const;

        bool hasGlobalBasisMatrix() const;

        CrsMatrixPtr getGlobalBasisMatrix() const;

    protected:

        ConstMapPtr SerialRowMap_;

        MapPtrVec UnassembledBasesMaps_;

        MultiVectorPtrVec UnassembledSubspaceBases_;

        MapPtr AssembledBasisMap_;

        MultiVectorPtr AssembledBasis_;

        CrsMatrixPtr GlobalBasisMatrix_;

    };

}

#endif
