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

#ifndef _FROSCH_OVERLAPPINGOPERATOR_DECL_HPP
#define _FROSCH_OVERLAPPINGOPERATOR_DECL_HPP

#include <FROSch_SchwarzOperator_def.hpp>


namespace FROSch {

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class OverlappingOperator : public SchwarzOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr               = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using MapPtr                = typename SchwarzOperator<SC,LO,GO,NO>::MapPtr;
        using ConstMapPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstMapPtr;

        using CrsMatrixPtr          = typename SchwarzOperator<SC,LO,GO,NO>::CrsMatrixPtr;
        using ConstCrsMatrixPtr     = typename SchwarzOperator<SC,LO,GO,NO>::ConstCrsMatrixPtr;

        using MultiVector           = typename SchwarzOperator<SC,LO,GO,NO>::MultiVector;
        using MultiVectorPtr        = typename SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr;

        using ImporterPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ImporterPtr;
        using ExporterPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ExporterPtr;

        using ParameterListPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using SubdomainSolverPtr    = typename SchwarzOperator<SC,LO,GO,NO>::SubdomainSolverPtr;

        using SCVecPtr              = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;

        using UN                    = typename SchwarzOperator<SC,LO,GO,NO>::UN;

    public:

        OverlappingOperator(ConstCrsMatrixPtr k,
                            ParameterListPtr parameterList);

        ~OverlappingOperator();

        virtual int initialize() = 0;

        virtual int compute() = 0;

        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          bool usePreconditionerOnly,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const;

    protected:

        enum CombinationType {Averaging,Full,Restricted};

        virtual int initializeOverlappingOperator();

        virtual int computeOverlappingOperator();

        ConstCrsMatrixPtr OverlappingMatrix_;

        ConstMapPtr OverlappingMap_;

        ImporterPtr Scatter_;

        SubdomainSolverPtr SubdomainSolver_;

        MultiVectorPtr Multiplicity_;

        CombinationType Combine_;

        int LevelID_;

    };

}

#endif
