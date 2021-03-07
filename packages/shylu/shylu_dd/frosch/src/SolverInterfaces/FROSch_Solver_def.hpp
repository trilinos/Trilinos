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

#ifndef _FROSCH_SOLVER_DEF_HPP
#define _FROSCH_SOLVER_DEF_HPP

#include <FROSch_Solver_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    typename Solver<SC,LO,GO,NO>::ConstXMapPtr Solver<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template<class SC,class LO,class GO,class NO>
    typename Solver<SC,LO,GO,NO>::ConstXMapPtr Solver<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template<class SC,class LO,class GO,class NO>
    bool Solver<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_;
    }

    template<class SC,class LO,class GO,class NO>
    bool Solver<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_;
    }

    template<class SC,class LO,class GO,class NO>
    void Solver<SC,LO,GO,NO>::residual(const XMultiVector & x,
                                       const XMultiVector & b,
                                       XMultiVector& r) const
    {
        apply(x,r);
        r.update(ScalarTraits<SC>::one(),b,-ScalarTraits<SC>::one());
    }

    template<class SC,class LO,class GO,class NO>
    Solver<SC,LO,GO,NO>::Solver(ConstXMatrixPtr k,
                                ParameterListPtr parameterList,
                                string description) :
    K_ (k),
    ParameterList_ (parameterList),
    Description_ (description)
    {}

}

#endif
