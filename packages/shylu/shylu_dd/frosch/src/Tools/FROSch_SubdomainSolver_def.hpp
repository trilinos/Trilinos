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

#ifndef _FROSCH_SUBDOMAINSOLVER_DEF_hpp
#define _FROSCH_SUBDOMAINSOLVER_DEF_hpp

#include <FROSch_SubdomainSolver_decl.hpp>

namespace FROSch {
    
    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::SubdomainSolver(CrsMatrixPtr k,
                                                  ParameterListPtr parameterList) :
    K_ (k),
    ParameterList_ (parameterList),
    EpetraLinearProblem_ (),
    AmesosSolver_ (),
    IsInitialized_ (false),
    IsComputed_ (false)
    {
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            FROSCH_ASSERT(K_->getRowMap()->lib()==Xpetra::UseEpetra,"UnderlyingLib!=Xpetra::UseEpetra");
            // AH 10/18/2017: Dies könnten wir nach initialize() verschieben, oder?
            Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            Xpetra::EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<Xpetra::EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
            EpetraCrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
            
            EpetraMultiVectorPtr xTmp;
            EpetraMultiVectorPtr bTmp;
            
            EpetraLinearProblem_.reset(new Epetra_LinearProblem(epetraMat.get(),xTmp.get(),bTmp.get()));
            
            Amesos amesosFactory;
            AmesosSolver_.reset(amesosFactory.Create(ParameterList_->get("Solver","Mumps"),*EpetraLinearProblem_));
            AmesosSolver_->SetParameters(ParameterList_->sublist("Amesos")); 
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==Xpetra::UseEpetra) {
                Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                Xpetra::EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<Xpetra::EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
                EpetraCrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
                
                EpetraMultiVectorPtr xTmp;
                EpetraMultiVectorPtr bTmp;
                
                Amesos2SolverEpetra_ = Amesos2::create<EpetraCrsMatrix,EpetraMultiVector>(ParameterList_->get("Solver","Mumps"),epetraMat,xTmp,bTmp);
                Amesos2SolverEpetra_->setParameters(sublist(ParameterList_,"Amesos2"));
            } else if (K_->getRowMap()->lib()==Xpetra::UseTpetra) {
                Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                TpetraCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrixNonConst();
                
                TpetraMultiVectorPtr xTmp;
                TpetraMultiVectorPtr bTmp;
                
                Amesos2SolverTpetra_ = Amesos2::create<Tpetra::CrsMatrix<SC,LO,GO,NO>,Tpetra::MultiVector<SC,LO,GO,NO> >(ParameterList_->get("Solver","Mumps"),tpetraMat,xTmp,bTmp);
                Amesos2SolverTpetra_->setParameters(sublist(ParameterList_,"Amesos2"));
            } else {
                FROSCH_ASSERT(0!=0,"This can't happen...");
            }
        } else {
            FROSCH_ASSERT(0!=0,"SolverType nicht bekannt...");
        }
    }
    
    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::~SubdomainSolver()
    {
        AmesosSolver_.reset();
        EpetraLinearProblem_.reset();
    }
    
    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::initialize()
    {
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
            AMESOS_CHK_ERR(AmesosSolver_->SymbolicFactorization());
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==Xpetra::UseEpetra) {
                IsInitialized_ = true;
                IsComputed_ = false;
                Amesos2SolverEpetra_->symbolicFactorization();
            } else {
                IsInitialized_ = true;
                IsComputed_ = false;
                Amesos2SolverTpetra_->symbolicFactorization();
            }
        } else {
            FROSCH_ASSERT(0!=0,"SolverType nicht bekannt...");
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::compute()
    {
        FROSCH_ASSERT(IsInitialized_,"!IsInitialized_.");
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            IsComputed_ = true;
            AMESOS_CHK_ERR(AmesosSolver_->NumericFactorization());
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==Xpetra::UseEpetra) {
                IsComputed_ = true;
                Amesos2SolverEpetra_->numericFactorization();
            } else {
                IsComputed_ = true;
                Amesos2SolverTpetra_->numericFactorization();
            }
        } else {
            FROSCH_ASSERT(0!=0,"SolverType nicht bekannt...");
        }
        return 0;
    }
    
    // Y = alpha * A^mode * X + beta * Y
    template<class SC,class LO,class GO,class NO>
    void SubdomainSolver<SC,LO,GO,NO>::apply(const MultiVector &x,
                                             MultiVector &y,
                                             Teuchos::ETransp mode,
                                             SC alpha,
                                             SC beta) const
    {
        FROSCH_ASSERT(IsComputed_,"!IsComputed_.");
        
        MultiVectorPtr yTmp;
        
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            const Xpetra::EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const Xpetra::EpetraMultiVectorT<GO,NO> *>(&x);
            Teuchos::RCP<Epetra_MultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();
            
            yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
            *yTmp = y;
            Xpetra::EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<Xpetra::EpetraMultiVectorT<GO,NO> *>(yTmp.get());
            Teuchos::RCP<Epetra_MultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();
            
            EpetraLinearProblem_->SetLHS(epetraMultiVectorY.get());
            EpetraLinearProblem_->SetRHS(epetraMultiVectorX.get());
            
            EpetraLinearProblem_->GetMatrix()->SetUseTranspose(mode==Teuchos::TRANS);
            AmesosSolver_->Solve();
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==Xpetra::UseEpetra) {
                const Xpetra::EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const Xpetra::EpetraMultiVectorT<GO,NO> *>(&x);
                Teuchos::RCP<Epetra_MultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();
                
                yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
                *yTmp = y;
                Xpetra::EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<Xpetra::EpetraMultiVectorT<GO,NO> *>(yTmp.get());
                Teuchos::RCP<Epetra_MultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();
                
                Amesos2SolverEpetra_->setX(epetraMultiVectorY);
                Amesos2SolverEpetra_->setB(epetraMultiVectorX);
                
                Amesos2SolverEpetra_->solve(); // Was ist, wenn man mit der transponierten Matrix lösen will
            } else {
                const Xpetra::TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const Xpetra::TpetraMultiVector<SC,LO,GO,NO> *>(&x);
                TpetraMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();
                
                yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
                *yTmp = y;
                const Xpetra::TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const Xpetra::TpetraMultiVector<SC,LO,GO,NO> *>(yTmp.get());
                TpetraMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();
                
                Amesos2SolverTpetra_->setX(tpetraMultiVectorY);
                Amesos2SolverTpetra_->setB(tpetraMultiVectorX);
                
                Amesos2SolverTpetra_->solve(); // Was ist, wenn man mit der transponierten Matrix lösen will
            }
        } else {
            FROSCH_ASSERT(0!=0,"SolverType nicht bekannt...");
        }
        y.update(alpha,*yTmp,beta);
    }
    
    template<class SC,class LO,class GO,class NO>
    typename SubdomainSolver<SC,LO,GO,NO>::ConstMapPtr SubdomainSolver<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }
    
    template<class SC,class LO,class GO,class NO>
    typename SubdomainSolver<SC,LO,GO,NO>::ConstMapPtr SubdomainSolver<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }
    
    template<class SC,class LO,class GO,class NO>
    void SubdomainSolver<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(0!=0,"describe() has be implemented properly...");
    }
    
    template<class SC,class LO,class GO,class NO>
    std::string SubdomainSolver<SC,LO,GO,NO>::description() const
    {
        return "Subdomain Solver";
    }
    
    template<class SC,class LO,class GO,class NO>
    bool SubdomainSolver<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_;
    }
    
    template<class SC,class LO,class GO,class NO>
    bool SubdomainSolver<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_;
    }
    
}

#endif
