// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_FROSCHPRECONDITIONER_DEF_HPP
#define _FROSCH_FROSCHPRECONDITIONER_DEF_HPP

#include <FROSch_FROSchPreconditioner_decl.hpp>
#include <FROSch_TwoLevelBlockPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int FROSchPreconditioner<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"FROSchPreconditioner::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int FROSchPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"FROSchPreconditioner::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::FROSchPreconditioner: !this->IsInitialized_");
        this->IsComputed_ = true;
        FROSchPreconditioner_->compute();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void FROSchPreconditioner<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                  XMultiVector &y,
                                                  ETransp mode,
                                                  SC alpha,
                                                  SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"FROSchPreconditioner::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::FROSchPreconditioner: !this->IsComputed_.");

        FROSchPreconditioner_->apply(x,y,mode,alpha,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int FROSchPreconditioner<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                        bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::FROSchPreconditioner: updateMatrix() is not implemented for the FROSchPreconditioner yet.");
    }

    template<class SC,class LO,class GO,class NO>
    FROSchPreconditioner<SC,LO,GO,NO>::FROSchPreconditioner(ConstXMatrixPtr k,
                                                            ParameterListPtr parameterList,
                                                            string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(FROSchPreconditionerTime,"FROSchPreconditioner::FROSchPreconditioner");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::FROSchPreconditioner: K_ is null.");

        if (!this->ParameterList_->get("Solver","TwoLevelPreconditioner").compare("TwoLevelPreconditioner")) {
            ParameterListPtr froschParameterList = sublist(this->ParameterList_,"FROSchPreconditioner");
            if (froschParameterList->isSublist(this->ParameterList_->get("Solver","TwoLevelPreconditioner"))) froschParameterList = sublist(froschParameterList,this->ParameterList_->get("Solver","TwoLevelPreconditioner"));
            // if (this->K_->getMap()->getComm()->getRank() == 0) froschParameterList->print(cout);

            ArrayRCP<RCP<const Map<LO,GO,NO> > > RepeatedMaps(1);
            ArrayRCP<RCP<const Map<LO,GO,NO> > > NodesMaps(1);

            UNVecPtr dofsPerNodeVector(1);
            ArrayRCP<DofOrdering> dofOrderings(1);
            ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > > dofsMapsVec(1);
            ArrayRCP<RCP<Map<LO,GO,NO> > > MainCoarseMapVector(1);

            if (froschParameterList->isParameter("Repeated Map Vector")) {
                RepeatedMaps = ExtractVectorFromParameterList<RCP<const Map<LO,GO,NO> > >(*froschParameterList,"Repeated Map Vector");
            }

            if (froschParameterList->isParameter("Nodes Map Vector")) {
                NodesMaps = ExtractVectorFromParameterList<RCP<const Map<LO,GO,NO> > >(*froschParameterList,"Nodes Map Vector");
            }

            if (froschParameterList->isParameter("DofsPerNode Vector")) {
                dofsPerNodeVector = ExtractVectorFromParameterList<UN>(*froschParameterList,"DofsPerNode Vector");
            }

            if (froschParameterList->isParameter("DofOrdering Vector")) {
                dofOrderings = ExtractVectorFromParameterList<DofOrdering>(*froschParameterList,"DofOrdering Vector");
            }

            if (froschParameterList->isParameter("Main Map Vector")) {
                MainCoarseMapVector = ExtractVectorFromParameterList<RCP<Map<LO,GO,NO> > >(*froschParameterList,"Main Map Vector");
            }
            if (froschParameterList->isParameter("Dofs Maps Vector")) {
                dofsMapsVec = ExtractVectorFromParameterList<ArrayRCP<RCP<const Map<LO,GO,NO> > >>(*froschParameterList,"Dofs Maps Vector");
            }

            ConstXMultiVectorPtr nodeList = null;
            GOVecPtr dirichletBoundaryDofs = null;
            ConstXMultiVectorPtr nullSpaceBasisVec = null;

            TwoLevelPreconditionerPtr twoLevelPreconditioner = rcp(new TwoLevelPreconditioner<SC,LO,GO,NO>(this->K_,froschParameterList));
            twoLevelPreconditioner->initialize(froschParameterList->get("Dimension",3),
                                               dofsPerNodeVector[0],
                                               froschParameterList->get("Overlap",1),
                                               nullSpaceBasisVec,
                                               nodeList,
                                               dofOrderings[0],
                                               RepeatedMaps[0],
                                               dofsMapsVec[0],
                                               dirichletBoundaryDofs);

            FROSchPreconditioner_ = twoLevelPreconditioner;
        } else if (!this->ParameterList_->get("Solver","TwoLevelPreconditioner").compare("TwoLevelBlockPreconditioner")) {
            ParameterListPtr froschParameterList = sublist(this->ParameterList_,"FROSchPreconditioner");
            if (froschParameterList->isSublist(this->ParameterList_->get("Solver","TwoLevelPreconditioner"))) froschParameterList = sublist(froschParameterList,this->ParameterList_->get("Solver","TwoLevelPreconditioner"));
            // if (this->K_->getMap()->getComm()->getRank() == 0) froschParameterList->print(cout);

            ArrayRCP<RCP<const Map<LO,GO,NO> > > RepeatedMaps(1);
            UNVecPtr dofsPerNodeVector;
            ConstXMultiVectorPtrVecPtr nullSpaceBasisVec = null;
            ArrayRCP<DofOrdering> dofOrderings = null;
            ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > > dofsMapsVec = null;

            FROSCH_ASSERT(froschParameterList->isParameter("Repeated Map Vector"),"FROSch::FROSchPreconditioner: Currently TwoLevelBlockPreconditioner cannot be constructed without Repeated Maps Vector ");
            FROSCH_ASSERT(froschParameterList->isParameter("DofsPerNode Vector"),"FROSch::FROSchPreconditioner: Currently, TwoLevelBlockPreconditioner cannot be constructed without DofsPerNode Vector.");
            FROSCH_ASSERT(froschParameterList->isParameter("DofOrdering Vector"),"FROSch::FROSchPreconditioner: Currently, TwoLevelBlockPreconditioner cannot be constructed without DofOrdering Vector.");

            if (froschParameterList->isParameter("Repeated Map Vector")) {
                RepeatedMaps = ExtractVectorFromParameterList<RCP<const Map<LO,GO,NO> > >(*froschParameterList,"Repeated Map Vector");
            }
            if (froschParameterList->isParameter("DofsPerNode Vector")) {
                dofsPerNodeVector = ExtractVectorFromParameterList<UN>(*froschParameterList,"DofsPerNode Vector");
            }
            if (froschParameterList->isParameter("DofOrdering Vector")) {
                dofOrderings = ExtractVectorFromParameterList<DofOrdering>(*froschParameterList,"DofOrdering Vector");
            }
            if (froschParameterList->isParameter("Dofs Maps Vector")) {
                dofsMapsVec = ExtractVectorFromParameterList<ArrayRCP<RCP<const Map<LO,GO,NO> > >>(*froschParameterList,"Dofs Maps Vector");
            }

            FROSCH_ASSERT(RepeatedMaps.size()==dofsPerNodeVector.size(),"FROSch::FROSchPreconditioner: RepeatedMaps.size()!=dofsPerNodeVector.size()");
            FROSCH_ASSERT(RepeatedMaps.size()==dofOrderings.size(),"FROSch::FROSchPreconditioner: RepeatedMaps.size()!=dofOrderings.size()");
            TwoLevelBlockPreconditionerPtr twoLevelBlockPreconditioner = rcp(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(this->K_,froschParameterList));
            twoLevelBlockPreconditioner->initialize(froschParameterList->get("Dimension",3),
                                                    dofsPerNodeVector,
                                                    dofOrderings,
                                                    froschParameterList->get("Overlap",1),
                                                    RepeatedMaps,
                                                    nullSpaceBasisVec);

            FROSchPreconditioner_ = twoLevelBlockPreconditioner;
        } else {
            FROSCH_ASSERT(false,"FROSch::FROSchPreconditioner: Solver unkown.");
        }
    }

}

#endif
