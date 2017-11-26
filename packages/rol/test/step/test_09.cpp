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


  
  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makeSharedFromRef(std::cout);
  else
    outStream = ROL::makeSharedFromRef(bhs);

  int errorFlag = 0;

  try {

    uint xo_dim = 3; // Dimension of optimization vectors
    uint ce_dim = 1; // Dimension of equality constraint
    uint ci_dim = 4; // Dimension of inequality constraint

    RealT left  = -1.0;
    RealT right =  1.0;

    // ----[ Full primal-dual vector ]----------------

    ROL::SharedPointer<vector> xo_ptr = ROL::makeShared<vector>(xo_dim,0.0); // opt
    ROL::SharedPointer<vector> xs_ptr = ROL::makeShared<vector>(ci_dim,0.0); // slack
    ROL::SharedPointer<vector> xe_ptr = ROL::makeShared<vector>(ce_dim,0.0); // equality multipliers
    ROL::SharedPointer<vector> xi_ptr = ROL::makeShared<vector>(ci_dim,0.0); // inequality multipliers

    ROL::SharedPointer<V> xo = ROL::makeShared<SV>(xo_ptr); 
    ROL::SharedPointer<V> xs = ROL::makeShared<SV>(xs_ptr);
    ROL::SharedPointer<V> xe = ROL::makeShared<SV>(xe_ptr);
    ROL::SharedPointer<V> xi = ROL::makeShared<SV>(xi_ptr);     

    ROL::RandomizeVector(*xo,left,right);
    ROL::RandomizeVector(*xs,left,right);
    ROL::RandomizeVector(*xe,left,right);
    ROL::RandomizeVector(*xi,left,right);

    ROL::SharedPointer<V> x = ROL::CreatePartitionedVector( xo, xs, xe, xi );
    

    // ----[ Full primal-dual direction vector ]------

    ROL::SharedPointer<vector> vo_ptr = ROL::makeShared<vector>(xo_dim,0.0); // opt
    ROL::SharedPointer<vector> vs_ptr = ROL::makeShared<vector>(ci_dim,0.0); // slack
    ROL::SharedPointer<vector> ve_ptr = ROL::makeShared<vector>(ce_dim,0.0); // equality multipliers
    ROL::SharedPointer<vector> vi_ptr = ROL::makeShared<vector>(ci_dim,0.0); // inequality multipliers
 
    ROL::SharedPointer<V> vo = ROL::makeShared<SV>(vo_ptr);
    ROL::SharedPointer<V> vs = ROL::makeShared<SV>(vs_ptr);
    ROL::SharedPointer<V> ve = ROL::makeShared<SV>(ve_ptr);
    ROL::SharedPointer<V> vi = ROL::makeShared<SV>(vi_ptr);     

    ROL::RandomizeVector(*vo,left,right);
    ROL::RandomizeVector(*vs,left,right);
    ROL::RandomizeVector(*ve,left,right);
    ROL::RandomizeVector(*vi,left,right);

    ROL::SharedPointer<V> v = ROL::CreatePartitionedVector( vo, vs, ve, vi );


    // ----[ Full primal-dual residual vector ]------

    ROL::SharedPointer<vector> ro_ptr = ROL::makeShared<vector>(xo_dim,0.0); // opt
    ROL::SharedPointer<vector> rs_ptr = ROL::makeShared<vector>(ci_dim,0.0); // slack
    ROL::SharedPointer<vector> re_ptr = ROL::makeShared<vector>(ce_dim,0.0); // equality multipliers
    ROL::SharedPointer<vector> ri_ptr = ROL::makeShared<vector>(ci_dim,0.0); // inequality multipliers
 
    ROL::SharedPointer<V> ro = ROL::makeShared<SV>(vo_ptr);
    ROL::SharedPointer<V> rs = ROL::makeShared<SV>(vs_ptr);
    ROL::SharedPointer<V> re = ROL::makeShared<SV>(ve_ptr);
    ROL::SharedPointer<V> ri = ROL::makeShared<SV>(vi_ptr);     

    ROL::RandomizeVector(*ro,left,right);
    ROL::RandomizeVector(*rs,left,right);
    ROL::RandomizeVector(*re,left,right);
    ROL::RandomizeVector(*ri,left,right);

    ROL::SharedPointer<V> r = ROL::CreatePartitionedVector( ro, rs, re, ri );

    // ----[ Primal-dual constraint ]-------

    ROL::SharedPointer<ROL::Objective<RealT> > obj_hs32 = 
      ROL::makeShared<ROL::ZOO::Objective_HS32<RealT>>();

    ROL::SharedPointer<ROL::EqualityConstraint<RealT> > eqcon_hs32 = 
      ROL::makeShared<ROL::ZOO::EqualityConstraint_HS32<RealT>>();
    
    ROL::SharedPointer<ROL::EqualityConstraint<RealT> > incon_hs32 = 
      ROL::makeShared<ROL::ZOO::InequalityConstraint_HS32<RealT>>();      


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


