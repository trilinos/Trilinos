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

//! \brief Verify the numerics using Sacado are working correctly 

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "FiniteElement.hpp"

#include <iostream>
#include <iomanip>


typedef double RealT;
typedef Sacado::Fad::SFad<RealT,1> FadType;


/* Examples of nonlinearities */

template<class ScalarT>
class Linear : public Nonlinearity<ScalarT> {
public:
    ScalarT operator() (const ScalarT &y, const ScalarT &y_x, const ScalarT &u, const ScalarT &x) {
         return y+2.0*u;
    } 
};

template<class ScalarT>
class Reaction : public Nonlinearity<ScalarT> {
public:
    ScalarT operator() (const ScalarT &y, const ScalarT &y_x, const ScalarT &u, const ScalarT &x) {
         return u*y;
    } 
};

template<class ScalarT>
class Advection : public Nonlinearity<ScalarT> {
public:
    ScalarT operator() (const ScalarT &y, const ScalarT &y_x, const ScalarT &u, const ScalarT &x) {
         return u*y_x;
    } 
};

template<class ScalarT>
class Burgers : public Nonlinearity<ScalarT> {
public:
    ScalarT operator() (const ScalarT &y, const ScalarT &y_x, const ScalarT &u, const ScalarT &x) {
         return y*y/2.0;
    } 
};







int main(int argc, char *argv[]) {

    // Set up MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprint > 0)
        outStream = Teuchos::rcp(&std::cout, false);
    else
        outStream = Teuchos::rcp(&bhs, false);

    int errorFlag = 0;

    int ni = atoi(argv[1]);
    int nq = atoi(argv[2]);
    int blk = atoi(argv[3]);
    
    bool diffTestFun = false;

    Teuchos::RCP<Nonlinearity<FadType> > func_fad = Teuchos::rcp(new Linear<FadType>);

    FiniteElement<FadType,RealT> fe_fad(ni,nq);

    Teuchos::RCP<std::vector<RealT> > x_rcp  = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp  = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > u_rcp  = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > f_rcp  = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > v_rcp  = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > jv_rcp = Teuchos::rcp( new std::vector<RealT>(ni, 0) );
    Teuchos::RCP<std::vector<RealT> > jac_rcp = Teuchos::rcp( new std::vector<RealT>(ni*ni, 0) );

    (*v_rcp)[0] = 1.0;
   
    StdVector<RealT> x(x_rcp);
    StdVector<RealT> y(y_rcp);
    StdVector<RealT> f(f_rcp);
    StdVector<RealT> u(u_rcp);
    StdVector<RealT> v(v_rcp);
    StdVector<RealT> jv(jv_rcp);
    StdVector<RealT> jac(jac_rcp);

    fe_fad.getInterpolationPoints(x);

    for(int i=0;i<ni;++i){ 
        (*u_rcp)[i] = 1.0;
        (*y_rcp)[i] = 1.0;
    } 

    Teuchos::RCP<std::vector<FadType> > y_fad_rcp   = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > v_fad_rcp   = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > jv_fad_rcp  = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > u_fad_rcp   = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > f_fad_rcp   = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > jac_fad_rcp = Teuchos::rcp( new std::vector<FadType> ); 

    y_fad_rcp->reserve(ni); 
    u_fad_rcp->reserve(ni); 
    v_fad_rcp->reserve(ni); 
    jv_fad_rcp->reserve(ni); 
    f_fad_rcp->reserve(ni); 
    jac_fad_rcp->reserve(ni*ni); 

    for(int i=0;i<ni;++i) {
        y_fad_rcp->push_back((*y_rcp)[i]);
        u_fad_rcp->push_back((*u_rcp)[i]);
        v_fad_rcp->push_back((*v_rcp)[i]);
        jv_fad_rcp->push_back(0);
        f_fad_rcp->push_back(0);
        for(int j=0;j<ni;++j){
            jac_fad_rcp->push_back(0);
        }
    }

    StdVector<FadType> y_fad(y_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> v_fad(v_fad_rcp);
    StdVector<FadType> jv_fad(jv_fad_rcp);
    StdVector<FadType> f_fad(f_fad_rcp);
    StdVector<FadType> jac_fad(jac_fad_rcp);
 
    fe_fad.vectorFunction(y_fad,u_fad,diffTestFun,func_fad,f_fad);
    fe_fad.applyJacobianBlock(y_fad,u_fad,diffTestFun,func_fad,v_fad,blk,jv_fad);
    fe_fad.buildJacobianBlock(y_fad,u_fad,diffTestFun,func_fad,blk,jac_fad);  

    for(int i=0;i<ni;++i) {
        (*f_rcp)[i] = (*f_fad_rcp)[i].val();
        (*jv_rcp)[i] = (*jv_fad_rcp)[i].val();
    }    

    printvec(*f_rcp,"f");
    printvec(*jv_rcp,"jv");

    std::cout << " Jacobian block " << blk << std::endl;
    for(int i=0;i<ni;++i) {
        for(int j=0;j<ni;++j) {
            (*jac_rcp)[i+ni*j] = (*jac_fad_rcp)[i+ni*j].val(); 
            std::cout << std::setw(14) <<(*jac_rcp)[i+ni*j];
        }  
       std::cout << std::endl;
    }   

   
 
    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

}
