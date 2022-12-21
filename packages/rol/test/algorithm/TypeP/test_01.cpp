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

/*! \file  test_23.cpp
    \brief Validate projected gradient algorithm.
*/

#include "ROL_TypeP_ProxGradientAlgorithm.hpp"
#include "ROL_l1Objective.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_StdObjective.hpp"

template<typename Real> 
class quad : public ROL::StdObjective<Real> { 
	private:
		int dim_; 
		std::vector<Real> a_, b_; 
	public: 
		quad(int dim) : dim_(dim) {
			a_.resize(dim); 
			b_.resize(dim); 

			for (int i = 0; i<dim; ++i){
				a_[i] = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				b_[i] = static_cast<Real>(2)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX) - static_cast<Real>(1);
			}
		}

		Real value(const std::vector<Real> &x, Real &tol){
			Real val(0); 
			for (int i = 0; i<dim_; ++i){
				val += static_cast<Real>(0.5)*a_[i]*x[i]*x[i] + b_[i]*x[i]; 
			}
			return val; 
		}

	 void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol){
		 
		 for (int i = 0; i<dim_; ++i){
			g[i] = a_[i]*x[i] + b_[i]; 
		 }
	 }
    void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol){
		 
		 for (int i = 0; i<dim_; ++i){
			hv[i] = a_[i]*v[i];  
		 }
	 }

		void getSolution(std::vector<Real> &x, const std::vector<Real> &wts, const std::vector<Real> &y) const {
			Real tmp(0); 

			for (int i = 0; i<dim_; ++i){
				tmp = -(b_[i] + wts[i])/a_[i]; 

				if (tmp > y[i]){
					x[i] = tmp;
			  }
				else{
					tmp = (wts[i] - b_[i])/a_[i]; 
					if (tmp < y[i]){
						x[i] = tmp; 
					}
					else{
						x[i] = y[i]; 
					}
				}
			}
		}
};

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    RealT tol = 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>());
    
    ROL::ParameterList list;
    list.sublist("Status Test").set("Gradient Tolerance",1e-7);
    list.sublist("Status Test").set("Constraint Tolerance",1e-8);
    list.sublist("Status Test").set("Step Tolerance",1e-12);
    list.sublist("Status Test").set("Iteration Limit", 250);
    list.sublist("Step").set("Type","Line Search");
    list.sublist("General").set("Output Level",iprint);
    
		int dim = 5; 
    ROL::Ptr<ROL::Vector<RealT>>     sol, wts, y;
    ROL::Ptr<ROL::Objective<RealT>>  nobj;
		ROL::Ptr<quad<RealT>> sobj; 
    ROL::Ptr<ROL::TypeP::ProxGradientAlgorithm<RealT>> algo;
    std::vector<RealT> data;
    RealT e1, e2, e3, e4, e5, err;

    *outStream << std::endl << "Diagonal LASSO" << std::endl << std::endl;
		ROL::Ptr<std::vector<RealT>> wtsP = ROL::makePtr<std::vector<RealT>>(dim); 
		ROL::Ptr<std::vector<RealT>> yP   = ROL::makePtr<std::vector<RealT>>(dim); 
		wts = ROL::makePtr<ROL::StdVector<RealT>>(wtsP);
    y   = ROL::makePtr<ROL::StdVector<RealT>>(yP);
    sol = ROL::makePtr<ROL::StdVector<RealT>>(dim);
		wts->randomize(0.0,1.0);
		y->randomize(-5.0, 5.0); 
		sol->zero(); 


		nobj = ROL::makePtr<ROL::l1Objective<RealT>>(wts,y); 
    sobj = ROL::makePtr<quad<RealT>>(dim); 

    algo = ROL::makePtr<ROL::TypeP::ProxGradientAlgorithm<RealT>>(list);
    algo->run(*sol,*sobj,*nobj,*outStream);

		std::vector<RealT> xstar(dim); 
		sobj->getSolution(xstar, *wtsP, *yP); 
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3] << "  x5 = " << data[4] << std::endl;
    e1 = (data[0]-xstar[0]);
    e2 = (data[1]-xstar[1]);
    e3 = (data[2]-xstar[2]);
    e4 = (data[3]-xstar[3]);
		e5 = (data[4]-xstar[4]); 
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    }
  
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}

