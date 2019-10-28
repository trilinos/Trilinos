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

#ifndef _FROSCH_SCHWARZOPERATOR_DEF_HPP
#define _FROSCH_SCHWARZOPERATOR_DEF_HPP

#include <FROSch_SchwarzOperator_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::SchwarzOperator(CommPtr comm) :
    MpiComm_ (comm),
    SerialComm_ (),
    K_ (),
    ParameterList_ (),
    Verbose_ (comm->getRank()==0),
    IsInitialized_ (false),
    IsComputed_ (false),
    LevelID_ (1)
    {
        SerialComm_ = createSerialComm<int>();
    }

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::SchwarzOperator(ConstXMatrixPtr k,
                                                  ParameterListPtr parameterList) :
    MpiComm_ (k->getRangeMap()->getComm()),
    SerialComm_ (),
    K_ (k),
    ParameterList_ (parameterList),
    Verbose_ (MpiComm_->getRank()==0),
    IsInitialized_ (false),
    IsComputed_ (false),
    LevelID_ (ParameterList_->get("Level ID",UN(1)))
    {
        FROSCH_ASSERT(getDomainMap()->isSameAs(*getRangeMap()),"SchwarzOperator assumes DomainMap==RangeMap");
        SerialComm_ = createSerialComm<int>();
    }

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::~SchwarzOperator()
    {

    }

    template<class SC,class LO,class GO,class NO>
    void SchwarzOperator<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                             XMultiVector &y,
                                             ETransp mode,
                                             SC alpha,
                                             SC beta) const
    {
        return apply(x,y,false,mode,alpha,beta);
    }

    template<class SC,class LO,class GO,class NO>
    typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr SchwarzOperator<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template<class SC,class LO,class GO,class NO>
    typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr SchwarzOperator<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template<class SC,class LO,class GO,class NO>
    bool SchwarzOperator<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_;
    }

    template<class SC,class LO,class GO,class NO>
    bool SchwarzOperator<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_;
    }

    template<class SC,class LO,class GO,class NO>
    int SchwarzOperator<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k) {
    // Maybe set IsComputed_ = false ? -> Go through code to be saver/cleaner
    // This function must be actively called by the user, and is only for recycling purposes.
    // The preconditioner is still computed and this point or the preconditioner was never computed and can now be computed with the matrix k now.
        K_ = k;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void SchwarzOperator<SC,LO,GO,NO>::residual(const XMultiVector & X,
                                                const XMultiVector & B,
                                                XMultiVector& R) const {
      SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one;
      apply(X,R);
      R.update(one,B,negone);
    }
}

#endif
