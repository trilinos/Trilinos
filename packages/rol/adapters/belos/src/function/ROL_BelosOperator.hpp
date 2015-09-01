// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/** 
    \class Belos::OperatorTraits<Scalar,ROL::MultiVector<Scalar>,ROL::LinearOperator<Scalar>> 
    \brief Provides interface for using ROL::LinearOperator with Belos solvers
           (excluding block solvers).    
    \author Created by Greg von Winckel
*/


#ifndef ROL_BELOS_OPERATOR_HPP
#define ROL_BELOS_OPERATOR_HPP

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"

#include "ROL_MultiVector.hpp"
#include "ROL_LinearOperator.hpp"

namespace Belos {

    template<class Scalar>
    class OperatorTraits<Scalar,ROL::MultiVector<Scalar>,ROL::LinearOperator<Scalar> > {
        public:
            static void
            /// Note that ROL and Belos use reverse orderings
            /// for input and output vectors
            Apply(const ROL::LinearOperator<Scalar>& Op,
                  const ROL::MultiVector<Scalar> &X,
                        ROL::MultiVector<Scalar> &Y,
                  const ETrans trans = NOTRANS) {
                 Scalar tol=0;
                 int n = X.getNumberOfVectors();
                 for(int i=0;i<n;++i) {
                     Op.apply(*(Y.getVector(i)),*(X.getVector(i)),tol);    
                 }
            }

            static bool
            HasApplyTranspose(const ROL::LinearOperator<Scalar>& Op) {
                return false;
            }
    };

}

#endif

