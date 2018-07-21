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
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    OverlappingOperator<SC,LO,GO,NO>::OverlappingOperator(CrsMatrixPtr k,
                                                          ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k,parameterList),
    OverlappingMatrix_ (),
    OverlappingMap_ (),
    RepeatedMap_(),
    Scatter_(),
    ScatterRestricted_(),
    SubdomainSolver_ (),
    Multiplicity_(),
    Restricted_(this->ParameterList_->get("Restricted",false))
#ifdef OVERLAPPING_TIMER
    ,OverlappingOperator_Init_Timer(Teuchos::TimeMonitor::getNewCounter("Overlapping Operator: Initialize")),
    OverlappingOperator_Compute_Timer(Teuchos::TimeMonitor::getNewCounter("Overlapping Operator: Compute")),
    OverlappingOperator_Apply_Timer(Teuchos::TimeMonitor::getNewCounter("Overlapping Operator: Apply"))
#endif
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

        Teuchos::TimeMonitor OverlappingOperator_Apply_TimeMonitor(*OverlappingOperator_Apply_Timer);
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

        xTmp->putScalar(0.0);
        yOverlap->replaceMap(OverlappingMap_);

        if (Restricted_){
            MultiVectorPtr yOverlapRestricted = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(RepeatedMap_,x.getNumVectors());
            yOverlapRestricted->doExport(*yOverlap,*ScatterRestricted_,Xpetra::INSERT);
            xTmp->doExport(*yOverlapRestricted,*Scatter_,Xpetra::ADD);
        }
        else{
            xTmp->doExport(*yOverlap,*Scatter_,Xpetra::ADD);
        }
        if (this->ParameterList_->get("Averaging",false)) {
            ConstSCVecPtr scaling = Multiplicity_->getData(0);
            for (unsigned j=0; j<xTmp->getNumVectors(); j++) {
                SCVecPtr values = xTmp->getDataNonConst(j);
                for (unsigned i=0; i<values.size(); i++) {
                    values[i] = values[i] / scaling[i];
                }
            }
        }
        
        if (!usePreconditionerOnly && mode != Teuchos::NO_TRANS) {
            this->K_->apply(*xTmp,*xTmp,mode,1.0,0.0);
        }
        y.update(alpha,*xTmp,beta);
    }
    
    template <class SC,class LO,class GO,class NO>
    int OverlappingOperator<SC,LO,GO,NO>::initializeOverlappingOperator()
    {
        Teuchos::TimeMonitor OverlappingOperator_Init_TimeMonitor(*OverlappingOperator_Init_Timer);
        
        if (Restricted_) {
            ScatterRestricted_ = Xpetra::ImportFactory<LO,GO,NO>::Build(RepeatedMap_,OverlappingMap_);
            Scatter_ = Xpetra::ImportFactory<LO,GO,NO>::Build(this->getDomainMap(),RepeatedMap_);
        }
        else{
            Scatter_ = Xpetra::ImportFactory<LO,GO,NO>::Build(this->getDomainMap(),OverlappingMap_);
        }
        if (this->ParameterList_->get("Averaging",false)) {
            Multiplicity_ = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(this->getRangeMap(),1);
            MultiVectorPtr multiplicityRepeated;
            if (Restricted_) {
                multiplicityRepeated = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(RepeatedMap_,1);
            }
            else{
                multiplicityRepeated = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(OverlappingMap_,1);
            }

            multiplicityRepeated->putScalar(1.);

            ExporterPtr multiplicityExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(multiplicityRepeated->getMap(),this->getRangeMap());
            Multiplicity_->doExport(*multiplicityRepeated,*multiplicityExporter,Xpetra::ADD);
        }
        
//        SubdomainSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(OverlappingMatrix_,sublist(this->ParameterList_,"Solver")));
//        SubdomainSolver_->initialize();
        return 0; // RETURN VALUE
    }
    
    template <class SC,class LO,class GO,class NO>
    int OverlappingOperator<SC,LO,GO,NO>::computeOverlappingOperator()
    {
        Teuchos::TimeMonitor OverlappingOperator_Compute_TimeMonitor(*OverlappingOperator_Compute_Timer);
        if (this->IsComputed_) {// already computed once and we want to recycle the information. That is why we reset OverlappingMatrix_ to K_, because K_ has been reset at this point
            OverlappingMatrix_ = this->K_;
        }
        
        OverlappingMatrix_ = ExtractLocalSubdomainMatrix(OverlappingMatrix_,OverlappingMap_);
        
        SubdomainSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(OverlappingMatrix_,sublist(this->ParameterList_,"Solver")));
        SubdomainSolver_->initialize();

        int ret = SubdomainSolver_->compute();

        return ret; // RETURN VALUE
    }
    
}

#endif
