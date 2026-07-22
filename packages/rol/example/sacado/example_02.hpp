// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_StdVector.hpp"

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

    ROL::Ptr<const std::vector<ScalarT> > xp =
        (dynamic_cast<const StdVector<ScalarT>&>(x)).getVector();

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

    typedef std::vector<ScalarT> vector;
    typedef StdVector<ScalarT>   SV;
      
    

    ROL::Ptr<vector> cp =  dynamic_cast<SV&>(c).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    ScalarT x1 = (*xp)[0];
    ScalarT x2 = (*xp)[1];
    ScalarT x3 = (*xp)[2];
    ScalarT x4 = (*xp)[3];
    ScalarT x5 = (*xp)[4];

    (*cp)[0] = x1*x1+x2*x2+x3*x3+x4*x4+x5*x5 - 10.0;
    (*cp)[1] = x2*x3 - 5.0*x4*x5;
    (*cp)[2] = x1*x1*x1 + x2*x2*x2 + 1.0;

}

