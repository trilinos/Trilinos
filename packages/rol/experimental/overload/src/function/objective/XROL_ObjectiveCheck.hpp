
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

#pragma once

namespace XROL {

/** \file  XROL_CheckObjective.hpp
    \brief Check gradient and Hessian using finite differences
           and check Hessian symmetry */

template<class X>
std::vector<std::vector<magnitude_t<X>>>
checkGradient( Objective<X>& obj, const X& x, const dual_t<X>& g,
             const X& d, std::ostream& os, Teuchos::ParameterList& parlist ) {

  using ROL::Finite_Difference_Arrays::shifts;
  using ROL::Finite_Difference_Arrays::weights;

  magnitude_t<X> eps      = ROL::ROL_EPSILON<magnitude_t<X>>();
  magnitude_t<X> sqrt_eps = std::sqrt(eps);

  auto objlist = parlist.sublist("Objective");
  auto gclist  = objlist.sublist("Gradient Check"); 
  
  bool printToStream = gclist.get("Print to Stream", true);

  index_t<X> numVals = 4;

  index_t<X> order    = static_cast<index_t<X>>(gclist.get("Finite Difference Order",1));
  index_t<X> numSteps = static_cast<index_t<X>>(gclist.get("Number of Steps",13));
  std::vector<magnitude_t<X>> steps(numSteps);
  
  if( gclist.isParameter("Finite Difference Steps") ) { 
    steps = Teuchos::getArrayFromStringParameter<magnitude_t<X>>(gclist,
              "Finite Difference Steps",numSteps).toVector();
  } else {
    double initialStepSize = gclist.get("Initial Step Size",1.0);
    double stepFactor = gclist.get("Step Reduction Factor",10.0);
    steps[0] = initialStepSize;
    for( index_t<X> i=1; i<numSteps; ++i) {
      steps[i] = steps[i-1]/stepFactor; 
    }
  }

  std::vector<magnitude_t<X>> tmp(numVals);
  std::vector<std::vector<magnitude_t<X>>> gCheck(numSteps, tmp);

  // Evaluate objective value at x. 
  obj.update(x);
  magnitude_t<X> val = obj.value(x,sqrt_eps);

  // Compute gradient at x.
  auto xnew = *clone(x);
  auto gtmp = *clone(g);

  //    obj.update(x); // Why do we need to update again ?
  obj.gradient(gtmp, x, sqrt_eps);

  // Evaluate the dual
  dual(xnew,gtmp);

  magnitude_t<X> dtg = dot(d,xnew);
  
  // Temporary vectors.

  for (index_t<X> i=0; i<numSteps; i++) {
    
    magnitude_t<X> eta = steps[i];

    set(xnew,x);

    // Compute gradient, finite-difference gradient, and absolute error.
    gCheck[i][0] = eta;
    gCheck[i][1] = dtg;

    gCheck[i][2] = weights[order-1][0] * val;

    for(index_t<X> j=0; j<order; ++j) {
      // Evaluate at x <- x+eta*c_i*d.
      axpy(xnew,eta*shifts[order-1][j], d);

      // Only evaluate at shifts where the weight is nonzero  
      if( weights[order-1][j+1] != 0 ) {
        obj.update(xnew);
        gCheck[i][2] += weights[order-1][j+1] * obj.value(xnew,sqrt_eps);
      }
    }

    gCheck[i][2] /= eta;

    gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);

    if (printToStream) {
      if (i==0) {
        os << std::right
           << std::setw(20) << "Step size"
           << std::setw(20) << "grad'*dir"
           << std::setw(20) << "FD approx"
           << std::setw(20) << "abs error"
           << "\n"
           << std::setw(20) << "---------"
           << std::setw(20) << "---------"
           << std::setw(20) << "---------"
           << std::setw(20) << "---------"
           << "\n";
      }
      os << std::scientific << std::setprecision(11) << std::right
         << std::setw(20) << gCheck[i][0]
         << std::setw(20) << gCheck[i][1]
         << std::setw(20) << gCheck[i][2]
         << std::setw(20) << gCheck[i][3]
         << "\n";
    }
  }
  return gCheck;   
}

template<class X>
std::vector<std::vector<magnitude_t<X>>>
checkHessVec( Objective<X>& obj, const X& x, const dual_t<X>& hv, const X& v, 
                   std::ostream &os, Teuchos::ParameterList& parlist ) {

  using ROL::Finite_Difference_Arrays::shifts;
  using ROL::Finite_Difference_Arrays::weights;

  magnitude_t<X> eps      = ROL::ROL_EPSILON<magnitude_t<X>>();
  magnitude_t<X> sqrt_eps = std::sqrt(eps);

  auto objlist = parlist.sublist("Objective");
  auto hvlist  = objlist.sublist("Hessian Check"); 
  
  bool printToStream = hvlist.get("Print to Stream", true);

  index_t<X> numVals = 4;

  index_t<X> order    = static_cast<index_t<X>>(hvlist.get("Finite Difference Order",1));
  index_t<X> numSteps = static_cast<index_t<X>>(hvlist.get("Number of Steps",10));

  std::vector<magnitude_t<X>> steps(numSteps);
  
  if( hvlist.isParameter("Finite Difference Steps") ) { 
    steps = Teuchos::getArrayFromStringParameter<magnitude_t<X>>(hvlist,
              "Finite Difference Steps",numSteps).toVector();
  } else {
    double initialStepSize = hvlist.get("Initial Step Size",1.0);
    double stepFactor = hvlist.get("Step Reduction Factor",10.0);
    steps[0] = initialStepSize;
    for( index_t<X> i=1; i<numSteps; ++i) {
      steps[i] = steps[i-1]/stepFactor; 
    }
  }

  std::vector<magnitude_t<X>> tmp(numVals);
  std::vector<std::vector<magnitude_t<X>>> hvCheck(numSteps, tmp);

  auto g = *clone(hv);
  obj.update(x);
  obj.gradient(g,x,sqrt_eps);
  
  auto Hv = *clone(hv);
  obj.hessVec(Hv,v,x,sqrt_eps);
  magnitude_t<X> normHv = norm(Hv);

  auto gdif = *clone(hv);
  auto gnew = *clone(hv);
  auto xnew = *clone(x);

  for (index_t<X> i=0; i<numSteps; i++) {

    magnitude_t<X> eta = steps[i]; 

    // Evaluate objective value at x+eta*d.
    set(xnew,x);

    set(gdif,g);
    scale(gdif,weights[order-1][0]);

    for(index_t<X> j=0; j<order; ++j) {
        
      // Evaluate at x <- x+eta*c_i*d.
      axpy(xnew,eta*shifts[order-1][j], v);

      // Only evaluate at shifts where the weight is nonzero  
      if( weights[order-1][j+1] != 0 ) {
          obj.update(xnew);
          obj.gradient(gnew, xnew, sqrt_eps); 
          axpy(gdif,weights[order-1][j+1],gnew);
      }
    }

    scale(gdif,1.0/eta);    

    // Compute norms of hessvec, finite-difference hessvec, and error.
    hvCheck[i][0] = eta;
    hvCheck[i][1] = normHv;
    hvCheck[i][2] = norm(gdif);
    axpy(gdif,-1.0, Hv);
    hvCheck[i][3] = norm(gdif);

    if (printToStream) {
      if (i==0) {
      os << std::right
         << std::setw(20) << "Step size"
         << std::setw(20) << "norm(Hess*vec)"
         << std::setw(20) << "norm(FD approx)"
         << std::setw(20) << "norm(abs error)"
         << "\n"
         << std::setw(20) << "---------"
         << std::setw(20) << "--------------"
         << std::setw(20) << "---------------"
         << std::setw(20) << "---------------"
         << "\n";
      }
      os << std::scientific << std::setprecision(11) << std::right
         << std::setw(20) << hvCheck[i][0]
         << std::setw(20) << hvCheck[i][1]
         << std::setw(20) << hvCheck[i][2]
         << std::setw(20) << hvCheck[i][3]
         << "\n";

    } // end if(printToStream)
  } // end for(index_t<X> i=0; i<numSteps; i++) 

  return hvCheck;
} // checkHessVec


template<class X>
std::vector<magnitude_t<X>>
checkHessSym( Objective<X>& obj, const X& x, const dual_t<X>& hv, const X& v, 
              const X &w, std::ostream &os, Teuchos::ParameterList& parlist ) {

  magnitude_t<X> eps      = ROL::ROL_EPSILON<magnitude_t<X>>();
  magnitude_t<X> sqrt_eps = std::sqrt(eps);

  auto h     = *clone(hv);
  auto hdual = *clone(x);

  obj.hessVec(h,v,x,sqrt_eps);

  dual(hdual,h);

  magnitude_t<X> wHv = dot(w,hdual);

  obj.hessVec(h,w,x,sqrt_eps);

  dual(hdual,h);
  magnitude_t<X> vHw = dot(v,hdual);

  std::vector<magnitude_t<X>> hsymCheck(3, 0);

  hsymCheck[0] = wHv;
  hsymCheck[1] = vHw;
  hsymCheck[2] = std::abs(vHw-wHv);
  
  auto objlist = parlist.sublist("Objective");
  auto hvlist  = objlist.sublist("Hessian Check"); 
  
  bool printToStream = hvlist.get("Print to Stream", true);


  if (printToStream) {
    os << std::right
       << std::setw(20) << "<w, H(x)v>"
       << std::setw(20) << "<v, H(x)w>"
       << std::setw(20) << "abs error"
       << "\n";
    os << std::scientific << std::setprecision(11) << std::right
       << std::setw(20) << hsymCheck[0]
       << std::setw(20) << hsymCheck[1]
       << std::setw(20) << hsymCheck[2]
       << "\n";
  }
  return hsymCheck;

} // checkHessSym

} // namespace XROL

