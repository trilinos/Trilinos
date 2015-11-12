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
//         
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_09.cpp
    \brief Test of Primal Dual Interior Point KKT system
*/

#include "ROL_HS32.hpp" 
#include "ROL_InteriorPointPrimalDualResidual.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_GMRES.hpp"

//template<class Real>



typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>            vector;
  typedef ROL::Vector<RealT>            V;
  typedef ROL::StdVector<RealT>         SV;
//  typedef ROL::PartitionedVector<RealT> PV;

  typedef typename vector::size_type    uint;


  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = rcp(&std::cout, false);
  else
    outStream = rcp(&bhs, false);

  int errorFlag = 0;

  try {

    uint xo_dim = 3; // Dimension of optimization vectors
    uint ce_dim = 1; // Dimension of equality constraint
    uint ci_dim = 4; // Dimension of inequality constraint

    RealT left  = -1.0;
    RealT right =  1.0;

    // ----[ Full primal-dual vector ]----------------

    RCP<vector> xo_rcp = rcp( new vector(xo_dim,0.0) ); // opt
    RCP<vector> xs_rcp = rcp( new vector(ci_dim,0.0) ); // slack
    RCP<vector> xe_rcp = rcp( new vector(ce_dim,0.0) ); // equality multipliers
    RCP<vector> xi_rcp = rcp( new vector(ci_dim,0.0) ); // inequality multipliers

    RCP<V> xo = rcp( new SV(xo_rcp) ); 
    RCP<V> xs = rcp( new SV(xs_rcp) );
    RCP<V> xe = rcp( new SV(xe_rcp) );
    RCP<V> xi = rcp( new SV(xi_rcp) );     

    ROL::RandomizeVector(*xo,left,right);
    ROL::RandomizeVector(*xs,left,right);
    ROL::RandomizeVector(*xe,left,right);
    ROL::RandomizeVector(*xi,left,right);

    RCP<V> x = ROL::CreatePartitionedVector( xo, xs, xe, xi );
    

    // ----[ Full primal-dual direction vector ]------

    RCP<vector> vo_rcp = rcp( new vector(xo_dim,0.0) ); // opt
    RCP<vector> vs_rcp = rcp( new vector(ci_dim,0.0) ); // slack
    RCP<vector> ve_rcp = rcp( new vector(ce_dim,0.0) ); // equality multipliers
    RCP<vector> vi_rcp = rcp( new vector(ci_dim,0.0) ); // inequality multipliers
 
    RCP<V> vo = rcp( new SV(vo_rcp) );
    RCP<V> vs = rcp( new SV(vs_rcp) );
    RCP<V> ve = rcp( new SV(ve_rcp) );
    RCP<V> vi = rcp( new SV(vi_rcp) );     

    ROL::RandomizeVector(*vo,left,right);
    ROL::RandomizeVector(*vs,left,right);
    ROL::RandomizeVector(*ve,left,right);
    ROL::RandomizeVector(*vi,left,right);

    RCP<V> v = ROL::CreatePartitionedVector( vo, vs, ve, vi );


    // ----[ Full primal-dual residual vector ]------

    RCP<vector> ro_rcp = rcp( new vector(xo_dim,0.0) ); // opt
    RCP<vector> rs_rcp = rcp( new vector(ci_dim,0.0) ); // slack
    RCP<vector> re_rcp = rcp( new vector(ce_dim,0.0) ); // equality multipliers
    RCP<vector> ri_rcp = rcp( new vector(ci_dim,0.0) ); // inequality multipliers
 
    RCP<V> ro = rcp( new SV(vo_rcp) );
    RCP<V> rs = rcp( new SV(vs_rcp) );
    RCP<V> re = rcp( new SV(ve_rcp) );
    RCP<V> ri = rcp( new SV(vi_rcp) );     

    ROL::RandomizeVector(*ro,left,right);
    ROL::RandomizeVector(*rs,left,right);
    ROL::RandomizeVector(*re,left,right);
    ROL::RandomizeVector(*ri,left,right);

    RCP<V> r = ROL::CreatePartitionedVector( ro, rs, re, ri );

    // ----[ Primal-dual constraint ]-------

    RCP<ROL::Objective<RealT> > obj_hs32 = 
      rcp( new ROL::ZOO::Objective_HS32<RealT> );

    RCP<ROL::EqualityConstraint<RealT> > eqcon_hs32 = 
      rcp( new ROL::ZOO::EqualityConstraint_HS32<RealT> );
    
    RCP<ROL::EqualityConstraint<RealT> > incon_hs32 = 
      rcp( new ROL::ZOO::InequalityConstraint_HS32<RealT> );      


    *outStream << "Performing finite difference check on Primal-Dual KKT system" 
               << std::endl;

    using ROL::InteriorPoint::PrimalDualResidual;    

    PrimalDualResidual<RealT> con(obj_hs32,eqcon_hs32,incon_hs32, *x);
   
    con.checkApplyJacobian(*x,*v,*r,true,*outStream);  
    
  }

  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  

  return 0;
}


