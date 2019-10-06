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

#ifndef _FROSCH_RGDSWPRECONDITIONER_DECL_HPP
#define _FROSCH_RGDSWPRECONDITIONER_DECL_HPP

#include <FROSch_OneLevelPreconditioner_def.hpp>
#include <FROSch_AlgebraicOverlappingPreconditioner_def.hpp>
#include <FROSch_RGDSWCoarseOperator_def.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class RGDSWPreconditioner : public AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr             = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr        = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMatrixPtr                = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMultiVectorPtr;

        using ParameterListPtr          = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using RGDSWCoarseOperatorPtr    = typename SchwarzPreconditioner<SC,LO,GO,NO>::RGDSWCoarseOperatorPtr;

        using UN                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::UN;

        using GOVecPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVecPtr;

    public:

        RGDSWPreconditioner(ConstXMatrixPtr k,
                            ParameterListPtr parameterList);

        int initialize(bool useDefaultParameters = true);

        int initialize(ConstXMapPtr repeatedMap,
                       bool useDefaultParameters = true);

        int initialize(GOVecPtr &dirichletBoundaryDofs,
                       bool useDefaultParameters = true);

        int initialize(ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs,
                       bool useDefaultParameters = true);

        int initialize(UN dimension,
                       int overlap);

        int initialize(UN dimension,
                       int overlap,
                       ConstXMapPtr repeatedMap);

        int initialize(UN dimension,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       ConstXMultiVectorPtr &nodeList);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs,
                       ConstXMultiVectorPtr &nodeList);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        std::string description() const; // @suppress("Type cannot be resolved")
        
        virtual int resetMatrix(ConstXMatrixPtr &k);

    protected:

        RGDSWCoarseOperatorPtr CoarseOperator_;
    };

}

#endif
