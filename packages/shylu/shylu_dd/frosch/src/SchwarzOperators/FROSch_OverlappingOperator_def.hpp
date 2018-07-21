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

#ifndef _FROSCH_OVERLAPPINGOPERATOR_DEF_HPP
#define _FROSCH_OVERLAPPINGOPERATOR_DEF_HPP

#include <FROSch_OverlappingOperator_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    OverlappingOperator<SC,LO,GO,NO>::OverlappingOperator(CrsMatrixPtr k,
                                                          ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k,parameterList),
    OverlappingMatrix_ (),
    OverlappingMap_ (),
    Scatter_ (),
    SubdomainSolver_ ()
    {
        
    }
    
    template <class SC,class LO,class GO,class NO>
    OverlappingOperator<SC,LO,GO,NO>::~OverlappingOperator()
    {
        SubdomainSolver_.reset();
    }
    
    // Y = alpha * A^mode * X + beta * Y
    template <class SC,class LO,class GO,class NO>
    void OverlappingOperator<SC,LO,GO,NO>::apply(const MultiVector &x,
                                                 MultiVector &y,
                                                 bool usePreconditionerOnly,
                                                 Teuchos::ETransp mode,
                                                 SC alpha,
                                                 SC beta) const
    {
        FROSCH_ASSERT(this->IsComputed_,"ERROR: OverlappingOperator has to be computed before calling apply()");
        
        MultiVectorPtr xTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
        *xTmp = x;
        
        if (!usePreconditionerOnly && mode == Teuchos::NO_TRANS) {
            this->K_->apply(x,*xTmp,mode,1.0,0.0);
        }
        
        MultiVectorPtr xOverlap = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(OverlappingMap_,x.getNumVectors());
        MultiVectorPtr yOverlap = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(OverlappingMatrix_->getDomainMap(),x.getNumVectors());
        
        xOverlap->doImport(*xTmp,*Scatter_,Xpetra::INSERT);
        
        xOverlap->replaceMap(OverlappingMatrix_->getRangeMap());
        
        SubdomainSolver_->apply(*xOverlap,*yOverlap,mode,1.0,0.0);
        
        yOverlap->replaceMap(OverlappingMap_);
        
        xTmp->putScalar(0.0);
        xTmp->doExport(*yOverlap,*Scatter_,Xpetra::ADD);
        
        if (!usePreconditionerOnly && mode != Teuchos::NO_TRANS) {
            this->K_->apply(*xTmp,*xTmp,mode,1.0,0.0);
        }
        y.update(alpha,*xTmp,beta);
    }
    
    template <class SC,class LO,class GO,class NO>
    int OverlappingOperator<SC,LO,GO,NO>::initializeOverlappingOperator()
    {
        OverlappingMatrix_ = ExtractLocalSubdomainMatrix(OverlappingMatrix_,OverlappingMap_);

        Scatter_ = Xpetra::ImportFactory<LO,GO,NO>::Build(this->getDomainMap(),OverlappingMap_);

        SubdomainSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(OverlappingMatrix_,sublist(this->ParameterList_,"Solver")));
        SubdomainSolver_->initialize();

        return 0; // RETURN VALUE
    }
    
    template <class SC,class LO,class GO,class NO>
    int OverlappingOperator<SC,LO,GO,NO>::computeOverlappingOperator()
    {
        int ret = FROSch::UpdateLocalSubdomainMatrix(this->K_,OverlappingMap_,OverlappingMatrix_);
        FROSCH_ASSERT(ret==0,"UpdateLocalSubdomainMatrix failed");
        
        ret = SubdomainSolver_->compute();
        
        return ret; // RETURN VALUE
    }
    
}

#endif
