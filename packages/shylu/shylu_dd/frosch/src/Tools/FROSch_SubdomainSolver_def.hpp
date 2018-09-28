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
                                                  ParameterListPtr parameterList,
                                                  GOVecPtr blockCoarseSize) :
    K_ (k),
    ParameterList_ (parameterList),
    EpetraLinearProblem_ (),
    AmesosSolver_ (),
    MueLuFactory_ (),
    MueLuHierarchy_ (),
    BelosLinearProblem_(),
    BelosSolverManager_(),
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
                ParameterListPtr parameterList = sublist(ParameterList_,"Amesos2");
                parameterList->setName("Amesos2");
                Amesos2SolverEpetra_->setParameters(parameterList);
            } else if (K_->getRowMap()->lib()==Xpetra::UseTpetra) {
                Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                TpetraCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrixNonConst();
                
                TpetraMultiVectorPtr xTmp;
                TpetraMultiVectorPtr bTmp;
                
                Amesos2SolverTpetra_ = Amesos2::create<Tpetra::CrsMatrix<SC,LO,GO,NO>,Tpetra::MultiVector<SC,LO,GO,NO> >(ParameterList_->get("Solver","Mumps"),tpetraMat,xTmp,bTmp);
                ParameterListPtr parameterList = sublist(ParameterList_,"Amesos2");
                parameterList->setName("Amesos2");
                Amesos2SolverTpetra_->setParameters(parameterList);
            } else {
                FROSCH_ASSERT(0!=0,"This can't happen...");
            }
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            
            MueLuFactory_ = Teuchos::rcp(new MueLu::ParameterListInterpreter<SC,LO,GO,NO>(parameterList->sublist("MueLu").sublist("MueLu Parameter")));
            Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > nullspace;

            if (!ParameterList_->sublist("MueLu").get("NullSpace","Laplace").compare("Laplace")) {
                nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(K_->getRowMap(), 1);
                nullspace->putScalar(1.);
            }
            else if (!ParameterList_->sublist("MueLu").get("NullSpace","Laplace").compare("SPP")) {
                FROSCH_ASSERT(blockCoarseSize.size()==2,"Wrong size of blockCoarseSize for MueLu nullspace...");
                unsigned dofs = (unsigned) ParameterList_->sublist("MueLu").get("Dimension",2);
                nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(K_->getRowMap(), dofs+1);
                //nullspace of upper part
                for (unsigned j=0; j<nullspace->getLocalLength(); j++) {
                    GO globIndex = nullspace->getMap()->getGlobalElement(j);
                    if (globIndex<=(GO)(dofs*blockCoarseSize[0]-1)) {
                        unsigned vecIndex = (globIndex)%dofs;
                        nullspace->getDataNonConst(vecIndex)[j] = 1.;
                    }
                    else{
                        nullspace->getDataNonConst(dofs)[j] = 1.;
                    }
                }
            }
            MueLuHierarchy_ = MueLuFactory_->CreateHierarchy();
            MueLuHierarchy_->GetLevel(0)->Set("A",K_);
            MueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace);
            
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xSolution;// = FROSch::ConvertToXpetra<SC, LO, GO, NO>(Xpetra::UseTpetra,*this->solution_,TeuchosComm);
            Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xRightHandSide;// = FROSch::ConvertToXpetra<SC, LO, GO, NO>(Xpetra::UseTpetra,*residualVec_,TeuchosComm);//hier residualVec. Bei linProb rhs_
            
            Teuchos::RCP<Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > > OpK = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(K_));
            
            
            BelosLinearProblem_.reset(new Belos::LinearProblem<SC,Xpetra::MultiVector<SC,LO,GO,NO>,Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,xRightHandSide));
            
            Belos::SolverFactory<SC,Xpetra::MultiVector<SC,LO,GO,NO>,Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > > belosFactory;
            ParameterListPtr solverParameterList = sublist(ParameterList_,"Belos");
            
            BelosSolverManager_ = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));
            
            BelosSolverManager_->setProblem(BelosLinearProblem_);
            
            
        } else {
            FROSCH_ASSERT(0!=0,"SolverType unknown...");
        }
    }
    
    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::~SubdomainSolver()
    {
        AmesosSolver_.reset();
        EpetraLinearProblem_.reset();
        
        Amesos2SolverEpetra_.reset();
        Amesos2SolverTpetra_.reset();
        
        MueLuFactory_.reset();
        MueLuHierarchy_.reset();
        
        BelosLinearProblem_.reset();
        BelosSolverManager_.reset();
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
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            IsInitialized_ = true;
            IsComputed_ = false;
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
        } else {
            FROSCH_ASSERT(0!=0,"SolverType unknown...");
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
            
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            MueLuFactory_->SetupHierarchy(*MueLuHierarchy_);
            MueLuHierarchy_->IsPreconditioner(false);
            IsComputed_ = true;
            
            
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            ParameterListPtr solverParameterList = sublist(ParameterList_,"Belos");
            if (solverParameterList->get("OneLevelPreconditioner",false)) {

                Teuchos::RCP<FROSch::OneLevelPreconditioner<SC,LO,GO,NO> > prec = Teuchos::rcp(new OneLevelPreconditioner<SC,LO,GO,NO> (K_, solverParameterList));

                bool buildRepeatedMap = solverParameterList->get("Build Repeated Map",false);
                prec->initialize(solverParameterList->get("Overlap",1),buildRepeatedMap);
                prec->compute();
                Teuchos::RCP<Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > > OpP = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(prec));

                if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
                    BelosLinearProblem_->setLeftPrec(OpP);
                    
                } else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
                    BelosLinearProblem_->setRightPrec(OpP);
                    
                } else {
                    FROSCH_ASSERT(0!=0,"PreconditionerPosition unknown...");
                }
            }
            IsComputed_ = true;
            
        } else {
            FROSCH_ASSERT(0!=0,"SolverType unknown...");
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
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());

            int mgridSweeps = ParameterList_->sublist("MueLu").get("mgridSweeps",-1);
            if (mgridSweeps>0) {
                MueLuHierarchy_->Iterate(x,*yTmp,mgridSweeps);
            }
            else{
                typename Teuchos::ScalarTraits<SC>::magnitudeType tol = ParameterList_->sublist("MueLu").get("tol",1.e-6);
                MueLuHierarchy_->Iterate(x,*yTmp,tol);
            }
            y = *yTmp;
            
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            
            ConstMultiVectorPtr xPtr = Teuchos::rcpFromRef(x);
            yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
            BelosLinearProblem_->setProblem(yTmp,xPtr);
            BelosSolverManager_->solve();
            y = *yTmp;
            
        } else {
            FROSCH_ASSERT(0!=0,"SolverType unknown...");
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
