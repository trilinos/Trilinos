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

#ifndef _FROSCH_ONELEVELPRECONDITIONER__DECL_HPP
#define _FROSCH_ONELEVELPRECONDITIONER__DECL_HPP

#include <FROSch_SchwarzPreconditioner_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;
    
    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class OneLevelPreconditioner : public SchwarzPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;

        using XMatrixPtr                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVector;

        using ParameterListPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using SumOperatorPtr                    = typename SchwarzPreconditioner<SC,LO,GO,NO>::SumOperatorPtr;
        using MultiplicativeOperatorPtr         = typename SchwarzPreconditioner<SC,LO,GO,NO>::MultiplicativeOperatorPtr;
        using OverlappingOperatorPtr            = typename SchwarzPreconditioner<SC,LO,GO,NO>::OverlappingOperatorPtr;
        using AlgebraicOverlappingOperatorPtr   = typename SchwarzPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingOperatorPtr;

    public:

        OneLevelPreconditioner(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        virtual int initialize(bool useDefaultParameters = true);

        virtual int initialize(int overlap = -1,
                               bool buildRepeatedMap = false);

        virtual int initialize(int overlap,
                               ConstXMapPtr repeatedMap);

        virtual int compute();

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

        virtual ConstXMapPtr getDomainMap() const;

        virtual ConstXMapPtr getRangeMap() const;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        virtual std::string description() const;

        virtual int resetMatrix(ConstXMatrixPtr &k);

    protected:

        ConstXMatrixPtr K_;

        SumOperatorPtr SumOperator_;
        MultiplicativeOperatorPtr MultiplicativeOperator_;
        OverlappingOperatorPtr OverlappingOperator_;
        bool UseMultiplicative_ = false;
    };

}

#endif
