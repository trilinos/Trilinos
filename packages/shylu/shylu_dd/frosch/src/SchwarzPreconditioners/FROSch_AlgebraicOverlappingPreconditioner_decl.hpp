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

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DECL_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DECL_HPP

#include <FROSch_SchwarzPreconditioner_def.hpp>

namespace FROSch {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class AlgebraicOverlappingPreconditioner : public SchwarzPreconditioner<SC,LO,GO,NO> {
        
    public:
        
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::MapPtr MapPtr;
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstMapPtr ConstMapPtr;
        
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::CrsMatrixPtr CrsMatrixPtr;
        
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::MultiVector MultiVector;
        
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr ParameterListPtr;
        
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::SumOperatorPtr SumOperatorPtr;
        typedef typename SchwarzPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingOperatorPtr AlgebraicOverlappingOperatorPtr;
        
        
        AlgebraicOverlappingPreconditioner(CrsMatrixPtr k,
                                           ParameterListPtr parameterList);
        
        virtual int initialize(bool useDefaultParameters = true);
        
        virtual int initialize(int overlap,
                               MapPtr repeatedMap);
        
        virtual int compute();
        
        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const;
        
        virtual ConstMapPtr getDomainMap() const;
        
        virtual ConstMapPtr getRangeMap() const;
        
        virtual void describe(Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
        
        virtual std::string description() const;
        
        
    protected:
        
        CrsMatrixPtr K_;
        
        SumOperatorPtr SumOperator_;
        AlgebraicOverlappingOperatorPtr FirstLevelOperator_;
        
    };
    
}

#endif
