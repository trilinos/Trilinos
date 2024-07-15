// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

