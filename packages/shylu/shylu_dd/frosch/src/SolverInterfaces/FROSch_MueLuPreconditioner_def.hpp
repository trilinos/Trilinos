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

#ifndef _FROSCH_MUELUPRECONDITIONER_DEF_HPP
#define _FROSCH_MUELUPRECONDITIONER_DEF_HPP

#include <FROSch_MueLuPreconditioner_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int MueLuPreconditioner<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"MueLuPreconditioner::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int MueLuPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"MueLuPreconditioner::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::MueLuPreconditioner: !this->IsInitialized_");
        MueLuFactory_->SetupHierarchy(*MueLuHierarchy_);
        MueLuHierarchy_->IsPreconditioner(false);
        this->IsComputed_ = true;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void MueLuPreconditioner<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                 XMultiVector &y,
                                                 ETransp mode,
                                                 SC alpha,
                                                 SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"MueLuPreconditioner::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::MueLuPreconditioner: !this->IsComputed_.");

        if (Y_.is_null()) Y_ = XMultiVectorFactory::Build(y.getMap(),y.getNumVectors());

        int mgridSweeps = this->ParameterList_->sublist("MueLu").get("mgridSweeps",-1);
        if (mgridSweeps>0) {
            MueLuHierarchy_->Iterate(x,*Y_,mgridSweeps);
        }
        else{
            typename ScalarTraits<SC>::magnitudeType tol = this->ParameterList_->sublist("MueLu").get("tol",1.e-6);
            MueLuHierarchy_->Iterate(x,*Y_,tol);
        }
        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int MueLuPreconditioner<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                       bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::MueLuPreconditioner: updateMatrix() is not implemented for the MueLuPreconditioner yet.");
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    MueLuPreconditioner<SC,LO,GO,NO>::MueLuPreconditioner(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList,
                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(MueLuPreconditionerTime,"MueLuPreconditioner::MueLuPreconditioner");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::MueLuPreconditioner: K_ is null.");

        MueLuFactory_ = rcp(new MueLu::ParameterListInterpreter<SC,LO,GO,NO>(this->ParameterList_->sublist("MueLu").sublist("MueLu Parameter")));
        RCP<XMultiVector> nullspace;

        if (!this->ParameterList_->sublist("MueLu").get("NullSpace","Laplace").compare("Laplace")) {
            nullspace = XMultiVectorFactory::Build(this->K_->getRowMap(), 1);
            nullspace->putScalar(1.);
        } else {
            FROSCH_ASSERT(false,"FROSch::MueLuPreconditioner: Only Laplacian null space is supported so far.");
        }
        MueLuHierarchy_ = MueLuFactory_->CreateHierarchy(); // Das vor den if block
        MueLuHierarchy_->GetLevel(0)->Set("A", Teuchos::rcp_const_cast<XMatrix>(this->K_)); // Das in den if block
        MueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace);
    }

}

#endif
