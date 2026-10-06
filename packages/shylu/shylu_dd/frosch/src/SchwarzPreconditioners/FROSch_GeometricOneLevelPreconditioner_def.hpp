// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GEOMETRICONELEVELPRECONDITIONER_DEF_HPP
#define _FROSCH_GEOMETRICONELEVELPRECONDITIONER_DEF_HPP

#include <FROSch_GeometricOneLevelPreconditioner_decl.hpp>
#include <FROSch_GeometricOverlappingOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;



    template <class SC,class LO,class GO,class NO>
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::GeometricOneLevelPreconditioner(
        ConstXMatrixPtr  k,
        GraphPtr         dualGraph,
        ParameterListPtr parameterList) 
    : SchwarzPreconditioner<SC,LO,GO,NO> (parameterList,k->getRangeMap()->getComm())
    , K_ (k)
    , SumOperator_ (new SumOperator<SC,LO,GO,NO>(k->getRangeMap()->getComm()))
    , MultiplicativeOperator_ (new MultiplicativeOperator<SC,LO,GO,NO>(k,parameterList))
    , OverlappingOperator_ ()
    {
        FROSCH_DETAILTIMER_START_LEVELID(oneLevelOptimizedPreconditionerTime, "GeometricOneLevelPreconditioner::GeometricOneLevelPreconditioner");

        if (!this->ParameterList_->get("OverlappingOperator Type", "GeometricOverlappingOperator").compare("GeometricOverlappingOperator")) {
            // Set the LevelID in the sublist
            int overlap = this->ParameterList_->get("Overlap",1);
            parameterList->sublist("GeometricOverlappingOperator").set("Level ID",this->LevelID_);
            OverlappingOperator_ = 
              GeometricOverlappingOperatorPtr(new GeometricOverlappingOperator<SC,LO,GO,NO>(k, overlap, dualGraph, sublist(parameterList,"GeometricOverlappingOperator")));
        } else {
            FROSCH_ASSERT(false,"Optimized operator type unkown.");
        }
        if (!this->ParameterList_->get("Level Combination", "Additive").compare("Multiplicative")) {
            UseMultiplicative_ = true;
        }
        if (UseMultiplicative_) {
            MultiplicativeOperator_->addOperator(OverlappingOperator_);
        } else {
            SumOperator_->addOperator(OverlappingOperator_);
        }

    }



    template <class SC,class LO,class GO,class NO>
    int
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::initialize(bool /*useDefaults*/)
    {
        FROSCH_TIMER_START_LEVELID(
            initializeTime, 
            "GeometricOneLevelPreconditioner::initialize");

        FROSCH_ASSERT(false, "GeometricOneLevelPreconditioner cannot be built without input parameters.");
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::initialize(XMapPtr overlappingMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime, "GeometricOneLevelPreconditioner::initialize");

        int ret = 0;

        if (!this->ParameterList_->get("OverlappingOperator Type", "GeometricOverlappingOperator").compare("GeometricOverlappingOperator")) {
            GeometricOverlappingOperatorPtr geometricOverlappingOperator = rcp_static_cast<GeometricOverlappingOperator<SC,LO,GO,NO> >(OverlappingOperator_);
            ret = geometricOverlappingOperator->initialize(overlappingMap);
        } else {
            FROSCH_ASSERT(false,"Optimized operator type unkown.");
        }
        return ret;
    }
 


    template <class SC,class LO,class GO,class NO>
    int
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime, "GeometricOneLevelPreconditioner::compute");

        FROSCH_ASSERT(false, "GeometricOneLevelPreconditioner cannot be computed without input parameters.");
    }


    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::compute(
        ConstXMatrixPtr neumannMatrix, 
        ConstXMatrixPtr robinMatrix)
    {
        FROSCH_TIMER_START_LEVELID(computeTime, "GeometricOneLevelPreconditioner::compute");

        int ret = 0;
        if (!this->ParameterList_->get("OverlappingOperator Type", "GeometricOverlappingOperator").compare("GeometricOverlappingOperator")) {
            GeometricOverlappingOperatorPtr geometricOverlappingOperator = rcp_static_cast<GeometricOverlappingOperator<SC,LO,GO,NO> >(OverlappingOperator_);
            ret = geometricOverlappingOperator->compute(neumannMatrix, robinMatrix);
        } else {
            FROSCH_ASSERT(false,"Optimized operator type unkown.");
        }
        return ret;
    }



    template <class SC,class LO,class GO,class NO>
    void GeometricOneLevelPreconditioner<SC,LO,GO,NO>::apply(
        const XMultiVector &x,
        XMultiVector       &y,
        ETransp             mode,
        SC                  alpha,
        SC                  beta) const
    {
        FROSCH_TIMER_START_LEVELID(applyTime, "GeometricOneLevelPreconditioner::apply");

        if (UseMultiplicative_) {
            return MultiplicativeOperator_->apply(x,y,true,mode,alpha,beta);
        }
        else{
            return SumOperator_->apply(x,y,true,mode,alpha,beta);
        }
    }



    template <class SC,class LO,class GO,class NO>
    const typename GeometricOneLevelPreconditioner<SC,LO,GO,NO>::ConstXMapPtr 
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }



    template <class SC,class LO,class GO,class NO>
    const typename GeometricOneLevelPreconditioner<SC,LO,GO,NO>::ConstXMapPtr 
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }



    template <class SC,class LO,class GO,class NO>
    void GeometricOneLevelPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                       const EVerbosityLevel verbLevel) const
    {
        if (UseMultiplicative_) {
            MultiplicativeOperator_->describe(out,verbLevel);
        }
        else{
            SumOperator_->describe(out,verbLevel);
        }

    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOneLevelPreconditioner<SC,LO,GO,NO>::communicateOverlappingTriangulation(
      XMultiVectorPtr                     nodeList,
      XMultiVectorTemplatePtr<long long>  elementList,
      XMultiVectorTemplatePtr<long long>  auxillaryList,
      XMultiVectorPtr                    &nodeListOverlapping,
      XMultiVectorTemplatePtr<long long> &elementListOverlapping,
      XMultiVectorTemplatePtr<long long> &auxillaryListOverlapping)

    {
        FROSCH_TIMER_START_LEVELID(computeTime, "GeometricOneLevelPreconditioner::compute");

        int ret = 0;
        if (!this->ParameterList_->get("OverlappingOperator Type", "GeometricOverlappingOperator").compare("GeometricOverlappingOperator")) {
            GeometricOverlappingOperatorPtr geometricOverlappingOperator = rcp_static_cast<GeometricOverlappingOperator<SC,LO,GO,NO> >(OverlappingOperator_);
            ret = geometricOverlappingOperator->communicateOverlappingTriangulation(
                      nodeList, elementList, auxillaryList, nodeListOverlapping, elementListOverlapping, auxillaryListOverlapping);
        } else {
            FROSCH_ASSERT(false,"Optimized operator type unkown.");
        }
        return ret;
    }



    template <class SC,class LO,class GO,class NO>
    string GeometricOneLevelPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "One-Level Optimized Schwarz Preconditioner";
    }



    template <class SC,class LO,class GO,class NO>
    int GeometricOneLevelPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_TIMER_START_LEVELID(resetMatrixTime,"OneLevelPreconditioner::resetMatrix");
        K_ = k;
        OverlappingOperator_->resetMatrix(K_);
        return 0;
    }

}

#endif // _FROSCH_GEOMETRICONELEVELPRECONDITIONER_DEF_HPP
