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

#include <iostream>

#include "Sacado.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_CompositeStepSQP.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace ROL;


/** \brief Objective function:  
    \f[f(x) = exp(x_1 x_2 x_3 x_4 x_5) + \frac{1}{2}*(x_1^3+x_2^3+1)^2 \f]
*/
template<class Real>
class Example_Objective {
    public:
        template<class ScalarT>
        ScalarT value(const Vector<ScalarT> &x, Real &tol);

};

template <class Real>
template <class ScalarT>
ScalarT Example_Objective<Real>::value(const Vector<ScalarT>& x, Real &tol) {
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    ScalarT x1 = (*xp)[0];
    ScalarT x2 = (*xp)[1];
    ScalarT x3 = (*xp)[2];
    ScalarT x4 = (*xp)[3];
    ScalarT x5 = (*xp)[4];

    ScalarT J = exp(x1*x2*x3*x4*x5) - 0.5 * pow( (pow(x1,3)+pow(x2,3)+1.0), 2);
    return J;  
}



/** \brief Equality constraint function with c_i(x) = 0, where :
    \f[ c_1(x) = x_1^2+x_2^2+x_3^2+x_4^2+x_5^2 - 10 \f]
    \f[ c_2(x) = x_2*x_3-5*x_4*x_5 \]
    \f[ c_3(x) = x_1^3 + x_2^3 + 1 \]
*/
template<class Real>
class Example_Constraint {
    public:
        int dim;
        
        template<class ScalarT>
        void value(Vector<ScalarT> &c, const Vector<ScalarT> &x, Real &tol);
};


template<class Real>
template<class ScalarT>
void Example_Constraint<Real>::value(Vector<ScalarT> &c, const Vector<ScalarT> &x, Real &tol) {
    Teuchos::RCP<std::vector<ScalarT> > cp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(c)).getVector());
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    ScalarT x1 = (*xp)[0];
    ScalarT x2 = (*xp)[1];
    ScalarT x3 = (*xp)[2];
    ScalarT x4 = (*xp)[3];
    ScalarT x5 = (*xp)[4];

    (*cp)[0] = x1*x1+x2*x2+x3*x3+x4*x4+x5*x5 - 10.0;
    (*cp)[1] = x2*x3 - 5.0*x4*x5;
    (*cp)[2] = x1*x1*x1 + x2*x2*x2 + 1.0;

}

