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

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::SubdomainSolver(ConstXMatrixPtr k,
                                                  ParameterListPtr parameterList,
                                                  GOVecPtr blockCoarseSize) :
    K_ (k),
    ParameterList_ (parameterList),
    YTmp_ (),
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    EpetraLinearProblem_ (),
#endif
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
    AmesosSolver_ (),
#endif
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
    MueLuFactory_ (),
    MueLuHierarchy_ (),
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
    BelosLinearProblem_(),
    BelosSolverManager_(),
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
    Ifpack2Preconditioner_ (),
#endif
    IsInitialized_ (false),
    IsComputed_ (false)
    {
        FROSCH_TIMER_START(subdomainSolverTime,"SubdomainSolver::SubdomainSolver");
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS

          FROSCH_ASSERT(K_->getRowMap()->lib()==UseEpetra,"UnderlyingLib!=UseEpetra");
            // AH 10/18/2017: Dies könnten wir nach initialize() verschieben, oder?
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
            ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
            TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

            EMultiVectorPtr xTmp;
            EMultiVectorPtr bTmp;

            EpetraLinearProblem_.reset(new ELinearProblem(epetraMat.get(),xTmp.get(),bTmp.get()));

            Amesos amesosFactory;

            AmesosSolver_.reset(amesosFactory.Create(ParameterList_->get("Solver","Mumps"),*EpetraLinearProblem_));

            AmesosSolver_->SetParameters(ParameterList_->sublist("Amesos"));
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Amesos");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
                ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
                TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

                EMultiVectorPtr xTmp;
                EMultiVectorPtr bTmp;

                Amesos2SolverEpetra_ = Amesos2::create<ECrsMatrix,EMultiVector>(ParameterList_->get("Solver","Mumps"),epetraMat,xTmp,bTmp);
                ParameterListPtr parameterList = sublist(ParameterList_,"Amesos2");
                parameterList->setName("Amesos2");
                Amesos2SolverEpetra_->setParameters(parameterList);
#else
                ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Epetra");
#endif
            } else if (K_->getRowMap()->lib()==UseTpetra) {
                const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
                TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

                TMultiVectorPtr xTmp;
                TMultiVectorPtr bTmp;

                Amesos2SolverTpetra_ = Amesos2::create<TCrsMatrix,TMultiVector>(ParameterList_->get("Solver","Mumps"),tpetraMat,xTmp,bTmp);
                ParameterListPtr parameterList = sublist(ParameterList_,"Amesos2");
                parameterList->setName("Amesos2");
                Amesos2SolverTpetra_->setParameters(parameterList);
            } else {
                FROSCH_ASSERT(false, "This can't happen. Either use Epetra or Tetra linear algebra stack.");
            }
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
            MueLuFactory_ = rcp(new MueLu::ParameterListInterpreter<SC,LO,GO,NO>(parameterList->sublist("MueLu").sublist("MueLu Parameter")));
            RCP<XMultiVector> nullspace;

            if (!ParameterList_->sublist("MueLu").get("NullSpace","Laplace").compare("Laplace")) {
                nullspace = XMultiVectorFactory::Build(K_->getRowMap(), 1);
                nullspace->putScalar(1.);
            }
            else if (!ParameterList_->sublist("MueLu").get("NullSpace","Laplace").compare("SPP")) { // Hier matrix zu block matrix konvertieren
                FROSCH_ASSERT(blockCoarseSize.size()==2,"Wrong size of blockCoarseSize for MueLu nullspace...");
                unsigned dofs = (unsigned) ParameterList_->sublist("MueLu").get("Dimension",2);
                nullspace = XMultiVectorFactory::Build(K_->getRowMap(), dofs+1);
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
            MueLuHierarchy_ = MueLuFactory_->CreateHierarchy(); // Das vor den if block
            MueLuHierarchy_->GetLevel(0)->Set("A",K_); // Das in den if block
            MueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace);
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "MueLu");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Ifpack2")) {
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
            FROSCH_ASSERT(K_->getRowMap()->lib()==UseTpetra,"FROSch::SubdomainSolver : ERROR: Ifpack2 is not compatible with Epetra.")

            // Convert matrix to Tpetra
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
            TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

            Ifpack2::Details::OneLevelFactory<TRowMatrix> ifpack2Factory;
            Ifpack2Preconditioner_ = ifpack2Factory.create(ParameterList_->get("Solver","FILU"),tpetraMat);

            ParameterListPtr parameterList = sublist(ParameterList_,"Ifpack2");
            parameterList->setName("Ifpack2");
            Ifpack2Preconditioner_->setParameters(*parameterList);
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Ifpack2");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
            RCP<XMultiVector> xSolution;// = FROSch::ConvertToXpetra<SC, LO, GO, NO>(UseTpetra,*this->solution_,TeuchosComm);
            RCP<XMultiVector> xRightHandSide;// = FROSch::ConvertToXpetra<SC, LO, GO, NO>(UseTpetra,*residualVec_,TeuchosComm);//hier residualVec. Bei linProb rhs_

            RCP<const Belos::OperatorT<XMultiVector> > OpK = rcp_dynamic_cast<const Belos::XpetraOp<SC,LO,GO,NO> >(K_);
            TEUCHOS_TEST_FOR_EXCEPT(OpK.is_null());

            BelosLinearProblem_.reset(new Belos::LinearProblem<SC,XMultiVector,Belos::OperatorT<XMultiVector> >(OpK,xSolution,xRightHandSide));

            Belos::SolverFactory<SC,XMultiVector,Belos::OperatorT<XMultiVector> > belosFactory;
            ParameterListPtr solverParameterList = sublist(ParameterList_,"Belos");

            BelosSolverManager_ = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));

            BelosSolverManager_->setProblem(BelosLinearProblem_);
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Belos");
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
    }

    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::~SubdomainSolver()
    {
        FROSCH_TIMER_START(subdomainSolverTime,"SubdomainSolver::~SubdomainSolver");
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        AmesosSolver_.reset();
#endif
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
        EpetraLinearProblem_.reset();
#endif

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
        Amesos2SolverEpetra_.reset();
#endif
        Amesos2SolverTpetra_.reset();

#ifdef HAVE_SHYLU_DDFROSCH_MUELU
        MueLuFactory_.reset();
        MueLuHierarchy_.reset();
#endif

#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        Ifpack2Preconditioner_.reset();
#endif

#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        BelosLinearProblem_.reset();
        BelosSolverManager_.reset();
#endif
    }

    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START(initializeTime,"SubdomainSolver::initialize");
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
            AMESOS_CHK_ERR(AmesosSolver_->SymbolicFactorization());
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                IsInitialized_ = true;
                IsComputed_ = false;
                Amesos2SolverEpetra_->symbolicFactorization();
#endif
            } else {
                IsInitialized_ = true;
                IsComputed_ = false;
                Amesos2SolverTpetra_->symbolicFactorization();
            }
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            IsInitialized_ = true;
            IsComputed_ = false;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Ifpack2")) {
            IsInitialized_ = true;
            IsComputed_ = false;
            Ifpack2Preconditioner_->initialize();
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START(computeTime,"SubdomainSolver::compute");
        FROSCH_ASSERT(this->IsInitialized_,"ERROR: SubdomainSolver has to be initialized before calling compute()");
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            IsComputed_ = true;
            AMESOS_CHK_ERR(AmesosSolver_->NumericFactorization());
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                IsComputed_ = true;
                Amesos2SolverEpetra_->numericFactorization();
#endif
            } else {
                IsComputed_ = true;
                Amesos2SolverTpetra_->numericFactorization();
            }
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            MueLuFactory_->SetupHierarchy(*MueLuHierarchy_);
            MueLuHierarchy_->IsPreconditioner(false);
            IsComputed_ = true;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Ifpack2")) {
            Ifpack2Preconditioner_->compute();
            IsComputed_ = true;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {
            ParameterListPtr solverParameterList = sublist(ParameterList_,"Belos");
            if (solverParameterList->get("OneLevelPreconditioner",false)) {

                RCP<FROSch::OneLevelPreconditioner<SC,LO,GO,NO> > prec = rcp(new OneLevelPreconditioner<SC,LO,GO,NO> (K_, solverParameterList));

                bool buildRepeatedMap = solverParameterList->get("Build Repeated Map",false);
                prec->initialize(solverParameterList->get("Overlap",1),buildRepeatedMap);
                prec->compute();
                RCP<Belos::OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(prec));

                if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
                    BelosLinearProblem_->setLeftPrec(OpP);

                } else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
                    BelosLinearProblem_->setRightPrec(OpP);

                } else {
                    FROSCH_ASSERT(false,"PreconditionerPosition unknown...");
                }
            }
            IsComputed_ = true;
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        return 0;


    }

    // Y = alpha * A^mode * X + beta * Y
    template<class SC,class LO,class GO,class NO>
    void SubdomainSolver<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                             XMultiVector &y,
                                             ETransp mode,
                                             SC alpha,
                                             SC beta) const
    {
        FROSCH_TIMER_START(applyTime,"SubdomainSolver::apply");
        FROSCH_ASSERT(IsComputed_,"!IsComputed_.");

#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos")) {
            const EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const EpetraMultiVectorT<GO,NO> *>(&x);
            RCP<EMultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();

            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
            *YTmp_ = y;
            EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<EpetraMultiVectorT<GO,NO> *>(YTmp_.get());
            RCP<EMultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();

            EpetraLinearProblem_->SetLHS(epetraMultiVectorY.get());
            EpetraLinearProblem_->SetRHS(epetraMultiVectorX.get());

            EpetraLinearProblem_->GetMatrix()->SetUseTranspose(mode==TRANS);
            AmesosSolver_->Solve();
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                const EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const EpetraMultiVectorT<GO,NO> *>(&x);
                RCP<EMultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();

                if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
                *YTmp_ = y;
                EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<EpetraMultiVectorT<GO,NO> *>(YTmp_.get());
                RCP<EMultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();

                Amesos2SolverEpetra_->setX(epetraMultiVectorY);
                Amesos2SolverEpetra_->setB(epetraMultiVectorX);

                Amesos2SolverEpetra_->solve(); // Was ist, wenn man mit der transponierten Matrix lösen will
#endif
            } else {
                const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
                TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

                if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
                *YTmp_ = y;
                const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(YTmp_.get());
                TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

                Amesos2SolverTpetra_->setX(tpetraMultiVectorY);
                Amesos2SolverTpetra_->setB(tpetraMultiVectorX);

                Amesos2SolverTpetra_->solve(); // Was ist, wenn man mit der transponierten Matrix lösen will
            }
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
        } else if (!ParameterList_->get("SolverType","Amesos").compare("MueLu")) {
            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());

            int mgridSweeps = ParameterList_->sublist("MueLu").get("mgridSweeps",-1);
            if (mgridSweeps>0) {
                MueLuHierarchy_->Iterate(x,*YTmp_,mgridSweeps);
            }
            else{
                typename ScalarTraits<SC>::magnitudeType tol = ParameterList_->sublist("MueLu").get("tol",1.e-6);
                MueLuHierarchy_->Iterate(x,*YTmp_,tol);
            }
            y = *YTmp_;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Ifpack2")) {
            const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
            TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
            *YTmp_ = y;
            const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(YTmp_.get());
            TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

            Ifpack2Preconditioner_->apply(*tpetraMultiVectorX,*tpetraMultiVectorY,mode,alpha,beta);
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos").compare("Belos")) {

            ConstXMultiVectorPtr xPtr = rcpFromRef(x);
            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
            BelosLinearProblem_->setProblem(YTmp_,xPtr);
            BelosSolverManager_->solve();
            y = *YTmp_;
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        y.update(alpha,*YTmp_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    typename SubdomainSolver<SC,LO,GO,NO>::ConstXMapPtr SubdomainSolver<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template<class SC,class LO,class GO,class NO>
    typename SubdomainSolver<SC,LO,GO,NO>::ConstXMapPtr SubdomainSolver<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template<class SC,class LO,class GO,class NO>
    void SubdomainSolver<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has be implemented properly...");
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
