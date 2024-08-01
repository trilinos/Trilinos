// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SUBDOMAINSOLVER_DEF_HPP
#define _FROSCH_SUBDOMAINSOLVER_DEF_HPP

#include <FROSch_SubdomainSolver_decl.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::SubdomainSolver(ConstXMatrixPtr k,
                                                  ParameterListPtr parameterList,
                                                  string description,
                                                  GOVecPtr blockCoarseSize) :
    K_ (k),
    ParameterList_ (parameterList),
    Description_ (description),
    IsInitialized_ (false),
    IsComputed_ (false)
    {
        FROSCH_TIMER_START_SUBDOMAINSOLVER(subdomainSolverTime,"SubdomainSolver::SubdomainSolver");
        FROSCH_ASSERT(!K_.is_null(),"FROSch::SubdomainSolver: K_ is null.");
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos")) {
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
          FROSCH_ASSERT(K_->getRowMap()->lib()==UseEpetra,"UnderlyingLib!=UseEpetra");
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
            // AH 10/18/2017: Dies könnten wir nach initialize() verschieben, oder?
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
            ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
            TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

            EMultiVectorPtr xTmp;
            EMultiVectorPtr bTmp;

            EpetraLinearProblem_.reset(new ELinearProblem(epetraMat.get(),xTmp.get(),bTmp.get()));

            Amesos amesosFactory;

            AmesosSolver_.reset(amesosFactory.Create(ParameterList_->get("Solver","Klu"),*EpetraLinearProblem_));

            AmesosSolver_->SetParameters(ParameterList_->sublist("Amesos"));
#endif
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Amesos");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
                ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
                TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

                EMultiVectorPtr xTmp;
                EMultiVectorPtr bTmp;

                Amesos2SolverEpetra_ = Amesos2::create<ECrsMatrix,EMultiVector>(ParameterList_->get("Solver","Klu"),epetraMat,xTmp,bTmp);
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

                Amesos2SolverTpetra_ = Amesos2::create<TCrsMatrix,TMultiVector>(ParameterList_->get("Solver","Klu"),tpetraMat,xTmp,bTmp);
                ParameterListPtr parameterList = sublist(ParameterList_,"Amesos2");
                parameterList->setName("Amesos2");
                Amesos2SolverTpetra_->setParameters(parameterList);
            } else {
                FROSCH_ASSERT(false, "This can't happen. Either use Epetra or Tetra linear algebra stack.");
            }
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("MueLu")) {
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
                ConstXMapPtr nullspaceMap = nullspace->getMap();
                //nullspace of upper part
                for (unsigned j=0; j<nullspace->getLocalLength(); j++) {
                    GO globIndex = nullspaceMap->getGlobalElement(j);
                    if (globIndex<=(GO)(dofs*blockCoarseSize[0]-1)) {
                        unsigned vecIndex = (globIndex)%dofs;
                        nullspace->getDataNonConst(vecIndex)[j] = 1.;
                    } else {
                        nullspace->getDataNonConst(dofs)[j] = 1.;
                    }
                }
            }
            MueLuHierarchy_ = MueLuFactory_->CreateHierarchy(); // Das vor den if block
            MueLuHierarchy_->GetLevel(0)->Set("A", Teuchos::rcp_const_cast<XMatrix>(K_)); // Das in den if block
            MueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace);
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "MueLu");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Ifpack2")) {
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
            FROSCH_ASSERT(K_->getRowMap()->lib()==UseTpetra,"FROSch::SubdomainSolver: Ifpack2 is not compatible with Epetra.")

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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {
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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Thyra")) {
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
            FROSCH_ASSERT(K_->getRowMap()->lib()==UseTpetra,"SubdomainSolver cannot use Epetra for Thyra solvers.");
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            RCP<const Thyra::LinearOpBase<SC> > thyraOp = ThyraUtils<SC,LO,GO,NO>::toThyra(crsOp.getCrsMatrix());

            Stratimikos::LinearSolverBuilder<SC> linearSolverBuilder;
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
            linearSolverBuilder.setPreconditioningStrategyFactory(abstractFactoryStd<Thyra::PreconditionerFactoryBase<SC>,Thyra::Ifpack2PreconditionerFactory<TCrsMatrix> >(),"Ifpack2");
#endif

            ParameterListPtr parameterList = sublist(ParameterList_,"Thyra");
            linearSolverBuilder.setParameterList(parameterList);

            RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lOWSFactory = linearSolverBuilder.createLinearSolveStrategy("");

            // AH 04/02/2020: Is this necessary?
            // RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
            // lOWSFactory->setOStream(out);
            // lOWSFactory->setVerbLevel(VERB_HIGH);

            LOWS_ = Thyra::linearOpWithSolve(*lOWSFactory,thyraOp);

#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Thyra");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelBlockPreconditioner")) {

            Teuchos::RCP< const Teuchos::Comm< int > > TC = K_->getMap()->getComm();
            Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > RepeatedMaps(1);
            UNVecPtr dofsPerNodeVector;
            ConstXMultiVectorPtrVecPtr nullSpaceBasisVec(1);
            Teuchos::ArrayRCP<DofOrdering> dofOrderings;
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > > dofsMapsVec = Teuchos::null;
            //Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > MainCoarseMapVector = Teuchos::null;

            FROSCH_ASSERT(ParameterList_->isParameter("Repeated Map Vector"),"Currently TwoLevelBlockPreconditioner cannot be constructed without Repeated Maps Vector ");
            FROSCH_ASSERT(ParameterList_->isParameter("DofsPerNode Vector"),"Currently, TwoLevelBlockPreconditioner cannot be constructed without DofsPerNode Vector.");
            FROSCH_ASSERT(ParameterList_->isParameter("DofOrdering Vector"),"Currently, TwoLevelBlockPreconditioner cannot be constructed without DofOrdering Vector.");

            if (ParameterList_->isParameter("Repeated Map Vector")) {
                RepeatedMaps = ExtractVectorFromParameterList<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > >(*ParameterList_,"Repeated Map Vector");
            }
            if (ParameterList_->isParameter("DofsPerNode Vector")) {
                dofsPerNodeVector = ExtractVectorFromParameterList<UN>(*ParameterList_,"DofsPerNode Vector");
            }
            if (ParameterList_->isParameter("DofOrdering Vector")) {
                dofOrderings = ExtractVectorFromParameterList<DofOrdering>(*ParameterList_,"DofOrdering Vector");
            }
            if (ParameterList_->isParameter("Dofs Maps Vector")) {
                dofsMapsVec = ExtractVectorFromParameterList<Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > >>(*ParameterList_,"Dofs Maps Vector");
            }

            FROSCH_ASSERT(RepeatedMaps.size()==dofsPerNodeVector.size(),"RepeatedMaps.size()!=dofsPerNodeVector.size()");
            FROSCH_ASSERT(RepeatedMaps.size()==dofOrderings.size(),"RepeatedMaps.size()!=dofOrderings.size()");
            TLBP = Teuchos::rcp(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(K_,ParameterList_));
            TLBP->initialize(ParameterList_->get("Dimension",3),
            dofsPerNodeVector,dofOrderings,
            ParameterList_->get("Overlap",1),
            RepeatedMaps,
            nullSpaceBasisVec);

            } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelPreconditioner")) {

            Teuchos::RCP< const Teuchos::Comm< int > > TC = K_->getMap()->getComm();
            Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > RepeatedMaps(1);
            Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > NodesMaps(1);

            UNVecPtr dofsPerNodeVector(1);
            Teuchos::ArrayRCP<DofOrdering> dofOrderings(1);
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > > dofsMapsVec = Teuchos::null;
            Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > MainCoarseMapVector = Teuchos::null;

            //FROSCH_ASSERT(ParameterList_->isParameter("Repeated Map Vector"),"Currently TwoLevelBlockPreconditioner cannot be constructed without Repeated Maps Vector ");
            if (ParameterList_->isParameter("Repeated Map Vector")) {
                RepeatedMaps = ExtractVectorFromParameterList<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > >(*ParameterList_,"Repeated Map Vector");
            }

            if (ParameterList_->isParameter("Nodes Map Vector")) {
                NodesMaps = ExtractVectorFromParameterList<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > >(*ParameterList_,"Nodes Map Vector");
            }

            if (ParameterList_->isParameter("DofsPerNode Vector")) {
                dofsPerNodeVector = ExtractVectorFromParameterList<UN>(*ParameterList_,"DofsPerNode Vector");
            }

            if (ParameterList_->isParameter("DofOrdering Vector")) {
                dofOrderings = ExtractVectorFromParameterList<DofOrdering>(*ParameterList_,"DofOrdering Vector");
            }

            if (ParameterList_->isParameter("Main Map Vector")) {
                MainCoarseMapVector = ExtractVectorFromParameterList<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > >(*ParameterList_,"Main Map Vector");
            }
            if (ParameterList_->isParameter("Dofs Maps Vector")) {
                dofsMapsVec = ExtractVectorFromParameterList<Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > >>(*ParameterList_,"Dofs Maps Vector");
            }

            ConstXMultiVectorPtr nodeList = null;
            GOVecPtr dirichletBoundaryDofs = null;
            ConstXMultiVectorPtr nullSpaceBasisVec = null;

            TLP = Teuchos::rcp(new TwoLevelPreconditioner<SC,LO,GO,NO>(K_,ParameterList_));
            TLP->initialize(ParameterList_->get("Dimension",3),
                            dofsPerNodeVector[0],
                            ParameterList_->get("Overlap",1),
                            nullSpaceBasisVec,
                            nodeList,
                            dofOrderings[0],
                            RepeatedMaps[0],
                            dofsMapsVec[0],
                            dirichletBoundaryDofs);
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
    }

    template<class SC,class LO,class GO,class NO>
    SubdomainSolver<SC,LO,GO,NO>::~SubdomainSolver()
    {
        FROSCH_TIMER_START_SUBDOMAINSOLVER(subdomainSolverTime,"SubdomainSolver::~SubdomainSolver");
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

#ifdef HAVE_SHYLU_DDFROSCH_THYRA
        LOWS_.reset();
#endif
    }

    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SUBDOMAINSOLVER(initializeTime,"SubdomainSolver::initialize");
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
            AMESOS_CHK_ERR(AmesosSolver_->SymbolicFactorization());
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos2")) {
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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("MueLu")) {
            IsInitialized_ = true;
            IsComputed_ = false;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Ifpack2")) {
            IsInitialized_ = true;
            IsComputed_ = false;
            Ifpack2Preconditioner_->initialize();
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {
            IsInitialized_ = true;
            IsComputed_ = false;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Thyra")) {
            // TODO: In the current implementation initialize() and compute() are part of apply() for Thyra.
            IsInitialized_ = true;
            IsComputed_ = false;
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelBlockPreconditioner")) {
            IsInitialized_ = true;
            IsComputed_ = false;
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelPreconditioner")) {
            IsInitialized_ = true;
            IsComputed_ = false;
        }  else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SUBDOMAINSOLVER(computeTime,"SubdomainSolver::compute");
        FROSCH_ASSERT(this->IsInitialized_,"ERROR: SubdomainSolver has to be initialized before calling compute()");
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos")) {
            IsComputed_ = true;
            AMESOS_CHK_ERR(AmesosSolver_->NumericFactorization());
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos2")) {
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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("MueLu")) {
            MueLuFactory_->SetupHierarchy(*MueLuHierarchy_);
            MueLuHierarchy_->IsPreconditioner(false);
            IsComputed_ = true;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Ifpack2")) {
            Ifpack2Preconditioner_->compute();
            IsComputed_ = true;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {
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
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Thyra")) {
            IsComputed_ = true;
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelBlockPreconditioner")) {
            TLBP->compute();
            IsComputed_ = true;
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelPreconditioner")) {
            TLP->compute();
            IsComputed_ = true;
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
        FROSCH_TIMER_START_SUBDOMAINSOLVER(applyTime,"SubdomainSolver::apply");
        FROSCH_ASSERT(IsComputed_,"!IsComputed_.");

#if defined(HAVE_SHYLU_DDFROSCH_AMESOS) && defined(HAVE_SHYLU_DDFROSCH_EPETRA)
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos")) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
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
#endif
        } else
#endif
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos2")) {
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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("MueLu")) {
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
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Ifpack2")) {
            const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
            TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
            *YTmp_ = y;
            const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(YTmp_.get());
            TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

            Ifpack2Preconditioner_->apply(*tpetraMultiVectorX,*tpetraMultiVectorY,mode,alpha,beta);
#endif
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {

            ConstXMultiVectorPtr xPtr = rcpFromRef(x);
            if (YTmp_.is_null()) YTmp_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
            BelosLinearProblem_->setProblem(YTmp_,xPtr);
            BelosSolverManager_->solve();
            y = *YTmp_;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Thyra")) {
            ConstXMultiVectorPtr xPtr = rcpFromRef(x);
            RCP<const Thyra::MultiVectorBase<SC> > thyraX = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xPtr);

            // AH 04/02/2020: Is there an easier way to do this? This seems to be close to the miminum in order to have two connected Thyra and Xpetra Multivectors
            if (ThyraYTmp_.is_null()) {
                XMultiVectorPtr yTmp = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
                ThyraYTmp_ = rcp_const_cast<Thyra::MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(yTmp));
            }
            const RCP<Tpetra::MultiVector<SC,LO,GO,NO> > yTmpTpMultVec = Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO>::getTpetraMultiVector(ThyraYTmp_);
            TEUCHOS_TEST_FOR_EXCEPT(is_null(yTmpTpMultVec));
            YTmp_ = rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(yTmpTpMultVec));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(YTmp_));

            Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*LOWS_, Thyra::NOTRANS, *thyraX, ThyraYTmp_.ptr());
            y = *YTmp_;

            /*
            const RCP<Tpetra::MultiVector<SC,LO,GO,NO> > yTpMultVec = Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO>::getTpetraMultiVector(rcpFromPtr(Y_inout));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(yTpMultVec));
            xY = rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(yTpMultVec));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(xY));
            */
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelBlockPreconditioner")) {
            if (YTmp_.is_null()) YTmp_ = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
            TLBP->apply(x,*YTmp_,Teuchos::NO_TRANS);
            y = *YTmp_;
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("TwoLevelPreconditioner")) {
            if (YTmp_.is_null()) YTmp_ = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),x.getNumVectors());
            TLP->apply(x,*YTmp_,Teuchos::NO_TRANS);
            y = *YTmp_;
        }  else {
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
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template<class SC,class LO,class GO,class NO>
    string SubdomainSolver<SC,LO,GO,NO>::description() const
    {
        return "Subdomain Solver"; // Add this->Description_;
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

    template<class SC,class LO,class GO,class NO>
    int SubdomainSolver<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr k,
                                                   bool reuseInitialize)
    {
        FROSCH_TIMER_START_SUBDOMAINSOLVER(resetMatrixTime,"SubdomainSolver::resetMatrix");
        K_ = k;
        FROSCH_ASSERT(!K_.is_null(),"FROSch::SubdomainSolver: K_ is null.");
        if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos")) {
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
            FROSCH_ASSERT(K_->getRowMap()->lib()==UseEpetra,"UnderlyingLib!=UseEpetra");
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
            // AH 10/18/2017: Dies könnten wir nach initialize() verschieben, oder?
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
            ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
            TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

            FROSCH_ASSERT(false,"FROSch::SubdomainSolver: resetMatrix() is not implemented for Amesos yet.");
#endif
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Amesos");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Amesos2")) {
            if (K_->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
                ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
                TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

                if (reuseInitialize) {
                    Amesos2SolverEpetra_->setA(epetraMat,Amesos2::SYMBFACT);
                } else {
                    Amesos2SolverEpetra_->setA(epetraMat,Amesos2::CLEAN);
                }
                /////////////////////////////////////////////
#else
                ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Epetra");
#endif
            } else if (K_->getRowMap()->lib()==UseTpetra) {
                const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
                const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
                TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

                if (reuseInitialize) {
                    Amesos2SolverTpetra_->setA(tpetraMat,Amesos2::SYMBFACT);
                } else {
                    Amesos2SolverTpetra_->setA(tpetraMat,Amesos2::CLEAN);
                }
               /////////////////////////////////////////////
            } else {
                FROSCH_ASSERT(false, "This can't happen. Either use Epetra or Tetra linear algebra stack.");
            }
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("MueLu")) {
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
            FROSCH_ASSERT(false,"FROSch::SubdomainSolver: resetMatrix() is not implemented for MueLu yet.");
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "MueLu");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Ifpack2")) {
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
            FROSCH_ASSERT(K_->getRowMap()->lib()==UseTpetra,"FROSch::SubdomainSolver: Ifpack2 is not compatible with Epetra.")

            // Convert matrix to Tpetra
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*K_);
            const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
            TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

            FROSCH_ASSERT(false,"FROSch::SubdomainSolver: resetMatrix() is not implemented for Ifpack2 yet.");
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Ifpack2");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
            FROSCH_ASSERT(false,"FROSch::SubdomainSolver: resetMatrix() is not implemented for Belos yet.");
#else
            ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Belos");
#endif
        } else if (!ParameterList_->get("SolverType","Amesos2").compare("Belos")) {
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
    FROSCH_ASSERT(false,"FROSch::SubdomainSolver: resetMatrix() is not implemented for Thyra yet.");
#else
    ThrowErrorMissingPackage("FROSch::SubdomainSolver", "Thyra");
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void SubdomainSolver<SC,LO,GO,NO>::residual(const XMultiVector & X,
                                                const XMultiVector & B,
                                                XMultiVector& R) const {
    SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one;
    apply(X,R);
    R.update(one,B,negone);
  }

}

#endif
