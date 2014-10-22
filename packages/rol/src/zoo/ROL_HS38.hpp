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

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 38th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS38_HPP
#define ROL_HS38_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"

namespace ROL {

  /** \brief W. Hock and K. Schittkowski 38th test function.
   */
  template<class Real>
  class Objective_HS38 : public Objective<Real> {
  public:
    Objective_HS38(void) {}

    Real value( const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      return 100.0 * std::pow((*ex)[1] - std::pow((*ex)[0],2.0),2.0) + std::pow(1.0-(*ex)[0],2.0) + 
              90.0 * std::pow((*ex)[3] - std::pow((*ex)[2],2.0),2.0) + std::pow(1.0-(*ex)[2],2.0) +
              10.1 * (std::pow((*ex)[1] - 1.0,2.0) + std::pow((*ex)[3]-1.0,2.0)) + 
              19.8 * ((*ex)[1] - 1.0) * ((*ex)[3] - 1.0);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > eg =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
      (*eg)[0] = -4.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) * (*ex)[0] - 2.0 * (1.0-(*ex)[0]);
      (*eg)[1] = 2.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) + 
                 2.0 * 10.1 * ((*ex)[1] - 1.0) + 19.8*((*ex)[3] - 1.0); 
      (*eg)[2] = -4.0 * 90.0 * ((*ex)[3] - std::pow((*ex)[2],2.0)) * (*ex)[2] - 2.0 * (1.0-(*ex)[2]);
      (*eg)[3] = 2.0 * 90.0 * ((*ex)[3] - std::pow((*ex)[2],2.0)) + 
                 2.0 * 10.1 * ((*ex)[3] - 1.0) + 19.8*((*ex)[1] - 1.0); 
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > ev =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > ehv =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
      Real h12 = -4.0 * 100.0 * (*ex)[0];
      Real h13 = 0.0;
      Real h14 = 0.0;
      Real h21 = -4.0 * 100.0 * (*ex)[0];
      Real h22 = 2.0 * 100.0 + 2.0 * 10.1;
      Real h23 = 0.0;
      Real h24 = 19.8;
      Real h31 = 0.0;
      Real h32 = 0.0;
      Real h33 = -4.0 * 90.0 * (*ex)[3] + 12.0 * 90.0 * std::pow((*ex)[2],2.0) + 2.0; 
      Real h34 = -4.0 * 90.0 * (*ex)[2];
      Real h41 = 0.0;
      Real h42 = 19.8;
      Real h43 = -4.0 * 90.0 * (*ex)[2];
      Real h44 = 2.0 * 90.0 + 2.0 * 10.1;

      (*ehv)[0] = h11 * (*ev)[0] + h12 * (*ev)[1] + h13 * (*ev)[2] + h14 * (*ev)[3];
      (*ehv)[1] = h21 * (*ev)[0] + h22 * (*ev)[1] + h23 * (*ev)[2] + h24 * (*ev)[3];
      (*ehv)[2] = h31 * (*ev)[0] + h32 * (*ev)[1] + h33 * (*ev)[2] + h34 * (*ev)[3];
      (*ehv)[3] = h41 * (*ev)[0] + h42 * (*ev)[1] + h43 * (*ev)[2] + h44 * (*ev)[3];
    } 
#endif
  };

  template<class Real>
  void getHS38( Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<BoundConstraint<Real> > &con, 
                Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
    // Resize Vectors
    n = 4;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_HS38<Real> );
    // Instantiate BoundConstraint
    std::vector<Real> l(n,-10.0);
    std::vector<Real> u(n,10.0);
    con = Teuchos::rcp( new StdBoundConstraint<Real>(l,u) );
    // Get Initial Guess
    (*x0p)[0] = -3.0;
    (*x0p)[1] = -1.0;
    (*x0p)[2] = -3.0;
    (*x0p)[3] = -1.0;
    // Get Solution
    (*xp)[0] = 1.0;
    (*xp)[1] = 1.0;
    (*xp)[2] = 1.0;
    (*xp)[3] = 1.0;
  }


}// End ROL Namespace

#endif
