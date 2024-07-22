// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _THYRA_FROSCH_FACTORY_DEF_HPP
#define _THYRA_FROSCH_FACTORY_DEF_HPP

#include "Thyra_FROSchFactory_decl.hpp"

#include <FROSch_AlgebraicOverlappingPreconditioner_def.hpp>
#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>
#include <FROSch_OneLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelBlockPreconditioner_def.hpp>

#include <FROSch_TpetraPreconditioner_def.hpp>

namespace Thyra {

    using namespace FROSch;
    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    //Constructor
    template <class SC, class LO, class GO, class NO>
    FROSchFactory<SC,LO,GO,NO>::FROSchFactory()
    {

    }

    //-----------------------------------------------------------
    //Check Type -> so far redundant
    template <class SC, class LO, class GO, class NO>
    bool FROSchFactory<SC,LO,GO,NO>::isCompatible(const LinearOpSourceBase<SC>& fwdOpSrc) const
    {
        const ConstLinearOpBasePtr fwdOp = fwdOpSrc.getOp();
        if (ThyraUtils<SC,LO,GO,NO>::isEpetra(fwdOp)) {
            return true;
        } else if (ThyraUtils<SC,LO,GO,NO>::isTpetra(fwdOp)) {
            return true;
        } else {
            return false;
        }
    }

    //--------------------------------------------------------------
    //Create Default Prec -> Not used here (maybe somewhere else?)
    template<class SC, class LO, class GO , class NO>
    typename FROSchFactory<SC,LO,GO,NO>::PreconditionerBasePtr FROSchFactory<SC,LO,GO,NO>::createPrec() const
    {
        return rcp(new DefaultPreconditioner<SC>);
    }

    //-------------------------------------------------------------
    //Main Function to use FROSch as Prec
    template<class SC, class LO , class GO, class NO>
    void FROSchFactory<SC,LO,GO,NO>::initializePrec(const ConstLinearOpSourceBasePtr& fwdOpSrc,
                                                    PreconditionerBase<SC>* prec,
                                                    const ESupportSolveUse supportSolveUse) const
    {
        //PreCheck
        TEUCHOS_ASSERT(nonnull(fwdOpSrc));
        //TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
        TEUCHOS_ASSERT(prec);

        // Retrieve wrapped concrete Xpetra matrix from FwdOp
        const ConstLinearOpBasePtr fwdOp = fwdOpSrc->getOp();
        TEUCHOS_TEST_FOR_EXCEPT(is_null(fwdOp));

        // Check whether it is Epetra/Tpetra
        bool bIsEpetra  = ThyraUtils<SC,LO,GO,NO>::isEpetra(fwdOp);
        bool bIsTpetra  = ThyraUtils<SC,LO,GO,NO>::isTpetra(fwdOp);
        bool bIsBlocked = ThyraUtils<SC,LO,GO,NO>::isBlockedOperator(fwdOp);
        TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == true  && bIsTpetra == true));
        TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == bIsTpetra) && bIsBlocked == false);
        TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra != bIsTpetra) && bIsBlocked == true);

        // Check whether to use HalfPrecision
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
        const bool useHalfPrecision = paramList_->get<bool>("Use Half Precision", false);
#endif

        // Retrieve Matrix
        ConstXMatrixPtr A = ThyraUtils<SC,LO,GO,NO>::toXpetra(fwdOp);
        TEUCHOS_TEST_FOR_EXCEPT(is_null(A));

        CommPtr comm = A->getMap()->getComm();
        UnderlyingLib underlyingLib = A->getMap()->lib();

        // Retrieve concrete preconditioner object
        const Ptr<DefaultPreconditioner<SC> > defaultPrec = ptr(dynamic_cast<DefaultPreconditioner<SC> *>(prec));
        TEUCHOS_TEST_FOR_EXCEPT(is_null(defaultPrec));

        // extract preconditioner operator
        auto precOp = defaultPrec->getNonconstUnspecifiedPrecOp();

        // Abstract SchwarzPreconditioner
        RCP<SchwarzPreconditioner<SC,LO,GO,NO> > schwarzPreconditioner = null;
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
        RCP<SchwarzPreconditioner<HalfSC,LO,GO,NO> > halfSchwarzPreconditioner = null;
#endif
        const bool startingOver = (precOp.is_null() || !paramList_->isParameter("Recycling") || !paramList_->get("Recycling",true));

        if (startingOver) {
            FROSCH_ASSERT(paramList_->isParameter("FROSch Preconditioner Type"),"FROSch Preconditioner Type is not defined!");

            if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("AlgebraicOverlappingPreconditioner")) {
                // Extract the repeated map
                ConstXMapPtr repeatedMap = extractRepeatedMap(comm,underlyingLib);

                RCP<AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> > AOP(new AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>(A,paramList_));

                AOP->initialize(paramList_->get("Overlap",1),
                                repeatedMap);

                schwarzPreconditioner = AOP;
            } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("GDSWPreconditioner")) {
                // Extract the repeated map
                ConstXMapPtr repeatedMap = extractRepeatedMap(comm,underlyingLib);

                // Extract the coordinate list
                ConstXMultiVectorPtr coordinatesList = extractCoordinatesList(comm,underlyingLib);

                // Extract the dof ordering
                DofOrdering dofOrdering = NodeWise;
                if (!paramList_->get("DofOrdering","NodeWise").compare("NodeWise")) {
                    dofOrdering = NodeWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("DimensionWise")) {
                    dofOrdering = DimensionWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("Custom")) {
                    dofOrdering = Custom;
                } else {
                    FROSCH_ASSERT(false,"ERROR: Specify a valid DofOrdering.");
                }

                RCP<GDSWPreconditioner<SC,LO,GO,NO> > GP(new GDSWPreconditioner<SC,LO,GO,NO>(A,paramList_));

                GP->initialize(paramList_->get("Dimension",3),
                               paramList_->get("DofsPerNode",1),
                               dofOrdering,
                               paramList_->get("Overlap",1),
                               repeatedMap,
                               coordinatesList);

                schwarzPreconditioner = GP;
            } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("RGDSWPreconditioner")) {
                // Extract the repeated map
                ConstXMapPtr repeatedMap = extractRepeatedMap(comm,underlyingLib);

                // Extract the coordinate list
                ConstXMultiVectorPtr coordinatesList = extractCoordinatesList(comm,underlyingLib);

                // Extract the dof ordering
                DofOrdering dofOrdering = NodeWise;
                if (!paramList_->get("DofOrdering","NodeWise").compare("NodeWise")) {
                    dofOrdering = NodeWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("DimensionWise")) {
                    dofOrdering = DimensionWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("Custom")) {
                    dofOrdering = Custom;
                } else {
                    FROSCH_ASSERT(false,"ERROR: Specify a valid DofOrdering.");
                }
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
                if (useHalfPrecision)
                {
                    auto ANonConst = rcp_const_cast<XMatrix>(A);
                    auto halfA = Xpetra::convertToHalfPrecision(ANonConst);

                    RCP<RGDSWPreconditioner<HalfSC,LO,GO,NO> > halfRGP(new RGDSWPreconditioner<HalfSC,LO,GO,NO>(halfA,paramList_));

                    XMultiVectorPtr coordinatesListNonConst = (coordinatesList == null ? null : rcp_const_cast<XMultiVector>(coordinatesList));
                    auto halfCoordinatesList = rcp_dynamic_cast<const HalfPrecMultiVector> (Xpetra::convertToHalfPrecision(coordinatesListNonConst));
                    halfRGP->initialize(paramList_->get("Dimension",3),
                                        paramList_->get("DofsPerNode",1),
                                        dofOrdering,
                                        paramList_->get("Overlap",1),
                                        repeatedMap,
                                        halfCoordinatesList);

                    halfSchwarzPreconditioner = halfRGP;
                } else
#endif
                {
                    RCP<RGDSWPreconditioner<SC,LO,GO,NO> > RGP(new RGDSWPreconditioner<SC,LO,GO,NO>(A,paramList_));

                    RGP->initialize(paramList_->get("Dimension",3),
                                    paramList_->get("DofsPerNode",1),
                                    dofOrdering,
                                    paramList_->get("Overlap",1),
                                    repeatedMap,
                                    coordinatesList);

                    schwarzPreconditioner = RGP;
                }
            } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("OneLevelPreconditioner")) {
                // Extract the repeated map
                ConstXMapPtr repeatedMap = extractRepeatedMap(comm,underlyingLib);

                RCP<OneLevelPreconditioner<SC,LO,GO,NO> > OLP(new OneLevelPreconditioner<SC,LO,GO,NO>(A,paramList_));

                OLP->initialize(paramList_->get("Overlap",1),
                                repeatedMap);

                schwarzPreconditioner = OLP;
            } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelPreconditioner")) {
                // Extract the repeated map
                ConstXMapPtr repeatedMap = extractRepeatedMap(comm,underlyingLib);

                // Extract the null space
                ConstXMultiVectorPtr nullSpaceBasis = extractNullSpace(comm,underlyingLib);

                // Extract the coordinate list
                ConstXMultiVectorPtr coordinatesList = extractCoordinatesList(comm,underlyingLib);

                // Extract the dof ordering
                DofOrdering dofOrdering = NodeWise;
                if (!paramList_->get("DofOrdering","NodeWise").compare("NodeWise")) {
                    dofOrdering = NodeWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("DimensionWise")) {
                    dofOrdering = DimensionWise;
                } else if (!paramList_->get("DofOrdering","NodeWise").compare("Custom")) {
                    dofOrdering = Custom;
                } else {
                    FROSCH_ASSERT(false,"ERROR: Specify a valid DofOrdering.");
                }
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
                if (useHalfPrecision)
                {
                    auto ANonConst = rcp_const_cast<XMatrix>(A);
                    auto halfA = Xpetra::convertToHalfPrecision(ANonConst);

                    RCP<TwoLevelPreconditioner<HalfSC,LO,GO,NO> > halfTLP(new TwoLevelPreconditioner<HalfSC,LO,GO,NO>(halfA,paramList_));

                    XMultiVectorPtr coordinatesListNonConst = rcp_const_cast<XMultiVector>(coordinatesList);
                    XMultiVectorPtr nullSpaceBasisNonConst = rcp_const_cast<XMultiVector>(nullSpaceBasis);
                    auto halfCoordinatesList = (coordinatesListNonConst == null ? null : rcp_dynamic_cast<const HalfPrecMultiVector> (Xpetra::convertToHalfPrecision(coordinatesListNonConst)));
                    auto halfNullSpaceBasis = (nullSpaceBasisNonConst == null ? null : rcp_dynamic_cast<const HalfPrecMultiVector> (Xpetra::convertToHalfPrecision(nullSpaceBasisNonConst)));
                    halfTLP->initialize(paramList_->get("Dimension",3),
                                        paramList_->get("DofsPerNode",1),
                                        paramList_->get("Overlap",1),
                                        halfNullSpaceBasis,
                                        halfCoordinatesList,
                                        dofOrdering,
                                        repeatedMap);

                    halfSchwarzPreconditioner = halfTLP;
                } else
#endif
                {
                    RCP<TwoLevelPreconditioner<SC,LO,GO,NO> > TLP(new TwoLevelPreconditioner<SC,LO,GO,NO>(A,paramList_));

                    TLP->initialize(paramList_->get("Dimension",3),
                                    paramList_->get("DofsPerNode",1),
                                    paramList_->get("Overlap",1),
                                    nullSpaceBasis,
                                    coordinatesList,
                                    dofOrdering,
                                    repeatedMap);

                    schwarzPreconditioner = TLP;
                }
            } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelBlockPreconditioner")) {
                ConstXMapPtrVecPtr repeatedMaps = null;
                ConstXMultiVectorPtrVecPtr coordinatesList = null;
                UNVecPtr dofsPerNodeVector;
                DofOrderingVecPtr dofOrderings;

                FROSCH_ASSERT(paramList_->isParameter("DofsPerNode Vector"),"Currently, TwoLevelBlockPreconditioner cannot be constructed without DofsPerNode Vector.");
                FROSCH_ASSERT(paramList_->isParameter("DofOrdering Vector"),"Currently, TwoLevelBlockPreconditioner cannot be constructed without DofOrdering Vector.");
                // Extract the repeated map vector
                if (paramList_->isParameter("Repeated Map Vector")) {
                    XMapPtrVecPtr repeatedMapsTmp = ExtractVectorFromParameterList<XMapPtr>(*paramList_,"Repeated Map Vector");
                    XMultiVectorPtrVecPtr nodeListVecTmp = ExtractVectorFromParameterList<XMultiVectorPtr>(*paramList_,"Coordinates List Vector");
                    if (!repeatedMapsTmp.is_null()) {
                        repeatedMaps.resize(repeatedMapsTmp.size());
                        for (unsigned i=0; i<repeatedMaps.size(); i++) {
                            repeatedMaps[i] = repeatedMapsTmp[i].getConst();
                        }
                    }
                    // Extract the nodeList map vector
                    if (!nodeListVecTmp.is_null()) {
                      coordinatesList.resize(nodeListVecTmp.size());
                      for (unsigned i = 0; i<coordinatesList.size();i++) {
                        coordinatesList[i] = nodeListVecTmp[i].getConst();
                      }
                    }

                    FROSCH_ASSERT(!repeatedMaps.is_null(),"FROSch::FROSchFactory: repeatedMaps.is_null()");
                    // Extract the DofsPerNode  vector
                    dofsPerNodeVector = ExtractVectorFromParameterList<UN>(*paramList_,"DofsPerNode Vector");
                    // Extract the DofOrdering vector
                    dofOrderings = ExtractVectorFromParameterList<DofOrdering>(*paramList_,"DofOrdering Vector");
                } else {
                    FROSCH_ASSERT(false,"Currently, TwoLevelBlockPreconditioner cannot be constructed without Repeated Maps.");
                }

                FROSCH_ASSERT(repeatedMaps.size()==dofsPerNodeVector.size(),"RepeatedMaps.size()!=dofsPerNodeVector.size()");
                FROSCH_ASSERT(repeatedMaps.size()==dofOrderings.size(),"RepeatedMaps.size()!=dofOrderings.size()");

                RCP<TwoLevelBlockPreconditioner<SC,LO,GO,NO> > TLBP(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(A,paramList_));

                TLBP->initialize(paramList_->get("Dimension",3),
                                 dofsPerNodeVector,
                                 dofOrderings,
                                 paramList_->get("Overlap",1),
                                 coordinatesList,
                                 repeatedMaps);

                schwarzPreconditioner = TLBP;
            } else {
                FROSCH_ASSERT(false,"Thyra::FROSchFactory: Preconditioner Type is unknown.");
            }

            //-----------------------------------------------

            LinearOpBasePtr thyraPrecOp = null;

#if defined(HAVE_XPETRA_TPETRA)
            if (bIsTpetra) {
                RCP<Tpetra::Operator<SC,LO,GO,NO> > tpOp;
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
                if (halfSchwarzPreconditioner != null)
                {
                    halfSchwarzPreconditioner->compute();
                    auto halfSchwarzOperator = rcp (new HalfPrecOp (halfSchwarzPreconditioner));
                    RCP<TpetraPreconditioner<SC,LO,GO,NO> > TPreconditioner (new TpetraPreconditioner<SC,LO,GO,NO>(halfSchwarzOperator));
                    tpOp = rcp_dynamic_cast<Tpetra::Operator<SC,LO,GO,NO> >(TPreconditioner);
                } else
#endif
                {
                    schwarzPreconditioner->compute();
                    auto schwarzOperator = rcp_dynamic_cast<Xpetra::Operator<SC,LO,GO,NO> >(schwarzPreconditioner);
                    RCP<TpetraPreconditioner<SC,LO,GO,NO> > TPreconditioner (new TpetraPreconditioner<SC,LO,GO,NO>(schwarzOperator));
                    tpOp = rcp_dynamic_cast<Tpetra::Operator<SC,LO,GO,NO> >(TPreconditioner);
                }
                thyraPrecOp = Thyra::createLinearOp(tpOp);
            } else
#endif
            {
                RCP<Operator<SC,LO,GO,NO> > xpOp;
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
                if (halfSchwarzPreconditioner != null) {
                    halfSchwarzPreconditioner->compute();
                    xpOp = rcp_dynamic_cast<Operator<SC,LO,GO,NO> >(rcp (new HalfPrecOp (halfSchwarzPreconditioner)));
                } else
#endif
                {
                    schwarzPreconditioner->compute();
                    xpOp = rcp_dynamic_cast<Operator<SC,LO,GO,NO> >(schwarzPreconditioner);
                }

                ConstVectorSpaceBasePtr thyraRangeSpace  = ThyraUtils<SC,LO,GO,NO>::toThyra(schwarzPreconditioner->getRangeMap());
                ConstVectorSpaceBasePtr thyraDomainSpace = ThyraUtils<SC,LO,GO,NO>::toThyra(schwarzPreconditioner->getDomainMap());
                thyraPrecOp = fROSchLinearOp<SC,LO,GO,NO>(thyraRangeSpace,thyraDomainSpace,xpOp,bIsEpetra,bIsTpetra);
            }
            TEUCHOS_TEST_FOR_EXCEPT(is_null(thyraPrecOp));

            //Set SchwarzPreconditioner
            defaultPrec->initializeUnspecified(thyraPrecOp);
        } else {
            // cast to SchwarzPreconditioner
            RCP<Operator<SC,LO,GO,NO> > xpetraOp;
#if defined(HAVE_XPETRA_TPETRA)
            if (bIsTpetra) {
                RCP<Thyra::TpetraLinearOp<SC,LO,GO,NO> > thyra_precOp = rcp_dynamic_cast<Thyra::TpetraLinearOp<SC,LO,GO,NO>>(precOp, true);
                RCP<Tpetra::Operator<SC,LO,GO,NO> > tpOp = thyra_precOp->getTpetraOperator();
                RCP<TpetraPreconditioner<SC,LO,GO,NO> > tpetra_precOp = rcp_dynamic_cast<TpetraPreconditioner<SC,LO,GO,NO>>(tpOp, true);
                xpetraOp = rcp_dynamic_cast<Operator<SC,LO,GO,NO> >(tpetra_precOp->getSchwarzPreconditioner());
            } else
#endif
            {
                RCP<FROSchLinearOp<SC,LO,GO,NO> > fROSch_LinearOp = rcp_dynamic_cast<FROSchLinearOp<SC,LO,GO,NO> >(precOp,true);
                xpetraOp = fROSch_LinearOp->getXpetraOperator();
            }

#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
            if (useHalfPrecision)
            {
                if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("RGDSWPreconditioner") ||
                    !paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelPreconditioner")) {
                    auto ANonConst = rcp_const_cast<XMatrix>(A);
                    auto halfA = rcp_const_cast<const HalfPrecMatrix>(Xpetra::convertToHalfPrecision(ANonConst));
                    auto halfPrecOp = rcp_dynamic_cast<HalfPrecOp>(xpetraOp);
                    auto halfXpetraOp = rcp_dynamic_cast<Operator<HalfSC,LO,GO,NO>>(halfPrecOp->GetHalfPrecisionOperator());
                    RCP<OneLevelPreconditioner<HalfSC,LO,GO,NO> > TLP = rcp_dynamic_cast<OneLevelPreconditioner<HalfSC,LO,GO,NO> >(halfXpetraOp, true);
                    TLP->resetMatrix(halfA);

                    halfSchwarzPreconditioner = rcp_dynamic_cast<SchwarzPreconditioner<HalfSC,LO,GO,NO> >(halfXpetraOp, true);
                    halfSchwarzPreconditioner->compute();
                }
            } else
#endif
            {
                if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("AlgebraicOverlappingPreconditioner")) {
                    RCP<AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> > AOP = rcp_dynamic_cast<AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    AOP->resetMatrix(A);
                    schwarzPreconditioner = AOP;
                } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("GDSWPreconditioner")) {
                    RCP<AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> > GP = rcp_dynamic_cast<AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    GP->resetMatrix(A);
                    schwarzPreconditioner = GP;
                } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("RGDSWPreconditioner")) {
                    RCP<OneLevelPreconditioner<SC,LO,GO,NO> > RGP = rcp_dynamic_cast<OneLevelPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    RGP->resetMatrix(A);
                    schwarzPreconditioner = RGP;
                } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("OneLevelPreconditioner")) {
                    RCP<OneLevelPreconditioner<SC,LO,GO,NO> > OLP = rcp_dynamic_cast<OneLevelPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    OLP->resetMatrix(A);
                    schwarzPreconditioner = OLP;
                } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelPreconditioner")) {
                    RCP<OneLevelPreconditioner<SC,LO,GO,NO> > TLP = rcp_dynamic_cast<OneLevelPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    TLP->resetMatrix(A);
                    schwarzPreconditioner = TLP;
                } else if (!paramList_->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelBlockPreconditioner")) {
                    RCP<OneLevelPreconditioner<SC,LO,GO,NO> > TLBP = rcp_dynamic_cast<OneLevelPreconditioner<SC,LO,GO,NO> >(xpetraOp, true);
                    TLBP->resetMatrix(A);
                    schwarzPreconditioner = TLBP;
                } else {
                    FROSCH_ASSERT(false,"Thyra::FROSchFactory: Preconditioner Type is unknown.");
                }
                // recompute SchwarzPreconditioner
                schwarzPreconditioner->compute();
            }
        }
    }

    //-------------------------------------------------------------
    //uninitialize
    template <class SC, class LO, class GO, class NO>
    void FROSchFactory<SC,LO,GO,NO>::uninitializePrec(PreconditionerBase<SC>* prec,
                                                      ConstLinearOpSourceBasePtr* fwdOp,
                                                      ESupportSolveUse* supportSolveUse) const
    {
        TEUCHOS_ASSERT(prec);

        // Retrieve concrete preconditioner object
        const Ptr<DefaultPreconditioner<SC> > defaultPrec = ptr(dynamic_cast<DefaultPreconditioner<SC> *>(prec));
        TEUCHOS_TEST_FOR_EXCEPT(is_null(defaultPrec));

        if (fwdOp) {
            // TODO: Implement properly instead of returning default value
            *fwdOp = null;
        }

        if (supportSolveUse) {
            // TODO: Implement properly instead of returning default value
            *supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED;
        }

        defaultPrec->uninitialize();
    }
    //-----------------------------------------------------------------
    //Following Functione maybe needed later
    template <class SC, class LO, class GO, class NO>
    void FROSchFactory<SC,LO,GO,NO>::setParameterList(ParameterListPtr const & paramList)
    {
        TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
        paramList_ = paramList;
    }

    template<class SC, class LO,class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ParameterListPtr FROSchFactory<SC,LO,GO,NO>::unsetParameterList()
    {
        ParameterListPtr savedParamList = paramList_;
        paramList_ = null;
        return savedParamList;
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ParameterListPtr FROSchFactory<SC,LO,GO,NO>::getNonconstParameterList()
    {
        return paramList_;
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ConstParameterListPtr FROSchFactory<SC,LO,GO,NO>::getParameterList() const
    {
        return paramList_;
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ConstParameterListPtr FROSchFactory<SC,LO,GO,NO>::getValidParameters() const
    {
        static ConstParameterListPtr validPL;

        if (is_null(validPL))
        validPL = rcp(new ParameterList());

        return validPL;
    }

    template <class SC, class LO, class GO, class NO>
    string FROSchFactory<SC,LO,GO,NO>::description() const
    {
        return "FROSchFactory";
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ConstXMapPtr FROSchFactory<SC,LO,GO,NO>::extractRepeatedMap(CommPtr comm,
                                                                                                     UnderlyingLib lib) const
    {
        ConstXMapPtr repeatedMap = null;
        if (paramList_->isParameter("Repeated Map")) {
            repeatedMap = ExtractPtrFromParameterList<XMap>(*paramList_,"Repeated Map").getConst();
            if (repeatedMap.is_null()) {
                if (lib==UseTpetra) { // If coordinatesList.is_null(), we look for Tpetra/Epetra RCPs
                    RCP<const Tpetra::Map<LO,GO,NO> > repeatedMapTmp = ExtractPtrFromParameterList<const Tpetra::Map<LO,GO,NO> >(*paramList_,"Repeated Map");

                    RCP<const TpetraMap<LO,GO,NO> > xTpetraRepeatedMap(new const TpetraMap<LO,GO,NO>(repeatedMapTmp));
                    repeatedMap = rcp_dynamic_cast<ConstXMap>(xTpetraRepeatedMap);
                } else {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                    FROSCH_WARNING("FROSch::FROSchFactory",comm->getRank()==0,"Cannot retrieve Epetra objects from ParameterList. Use Xpetra instead.");
#endif
                }
            }
            FROSCH_ASSERT(!repeatedMap.is_null(),"FROSch::FROSchFactory: repeatedMap.is_null()");
        }
        return repeatedMap;
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ConstXMultiVectorPtr FROSchFactory<SC,LO,GO,NO>::extractCoordinatesList(CommPtr comm,
                                                                                                                 UnderlyingLib lib) const
    {
        ConstXMultiVectorPtr coordinatesList = null;
        if (paramList_->isParameter("Coordinates List")) {
            coordinatesList = ExtractPtrFromParameterList<XMultiVector>(*paramList_,"Coordinates List").getConst();
            if (coordinatesList.is_null()) {
                if (lib==UseTpetra) { // If coordinatesList.is_null(), we look for Tpetra/Epetra RCPs
                    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > coordinatesListTmp = ExtractPtrFromParameterList<Tpetra::MultiVector<SC,LO,GO,NO> >(*paramList_,"Coordinates List");

                    RCP<const Xpetra::TpetraMultiVector<SC,LO,GO,NO> > xTpetraCoordinatesList(new const Xpetra::TpetraMultiVector<SC,LO,GO,NO>(coordinatesListTmp));
                    coordinatesList = rcp_dynamic_cast<ConstXMultiVector>(xTpetraCoordinatesList);
                } else {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                    FROSCH_WARNING("FROSch::FROSchFactory",comm->getRank()==0,"Cannot retrieve Epetra objects from ParameterList. Use Xpetra instead.");
#endif
                }
            }
            FROSCH_ASSERT(!coordinatesList.is_null(),"FROSch::FROSchFactory: coordinatesList.is_null()");
        }
        return coordinatesList;
    }

    template <class SC, class LO, class GO, class NO>
    typename FROSchFactory<SC,LO,GO,NO>::ConstXMultiVectorPtr FROSchFactory<SC,LO,GO,NO>::extractNullSpace(CommPtr comm,
                                                                                                           UnderlyingLib lib) const
    {
        ConstXMultiVectorPtr nullSpaceBasis = null;
        if (paramList_->isParameter("Null Space")) {
            nullSpaceBasis = ExtractPtrFromParameterList<XMultiVector>(*paramList_,"Null Space").getConst();
            if (nullSpaceBasis.is_null()) {
                if (lib==UseTpetra) { // If nullSpaceBasis.is_null(), we look for Tpetra/Epetra RCPs
                    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > nullSpaceBasisTmp = ExtractPtrFromParameterList<Tpetra::MultiVector<SC,LO,GO,NO> >(*paramList_,"Null Space");

                    RCP<const Xpetra::TpetraMultiVector<SC,LO,GO,NO> > xTpetraNullSpaceBasis(new const Xpetra::TpetraMultiVector<SC,LO,GO,NO>(nullSpaceBasisTmp));
                    nullSpaceBasis = rcp_dynamic_cast<ConstXMultiVector>(xTpetraNullSpaceBasis);
                } else {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                    FROSCH_WARNING("FROSch::FROSchFactory",comm->getRank()==0,"Cannot retrieve Epetra objects from ParameterList. Use Xpetra instead.");
#endif
                }
            }
            FROSCH_ASSERT(!nullSpaceBasis.is_null(),"FROSch::FROSchFactory: nullSpaceBasis.is_null()");
        }
        return nullSpaceBasis;
    }

}
#endif
