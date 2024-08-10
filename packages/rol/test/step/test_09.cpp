// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_09.cpp
    \brief Test of Primal Dual Interior Point KKT system
*/

#include "Teuchos_GlobalMPISession.hpp"

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
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    uint xo_dim = 3; // Dimension of optimization vectors
    uint ce_dim = 1; // Dimension of equality constraint
    uint ci_dim = 4; // Dimension of inequality constraint

    RealT left  = -1.0;
    RealT right =  1.0;

    // ----[ Full primal-dual vector ]----------------

    ROL::Ptr<vector> xo_ptr = ROL::makePtr<vector>(xo_dim,0.0); // opt
    ROL::Ptr<vector> xs_ptr = ROL::makePtr<vector>(ci_dim,0.0); // slack
    ROL::Ptr<vector> xe_ptr = ROL::makePtr<vector>(ce_dim,0.0); // equality multipliers
    ROL::Ptr<vector> xi_ptr = ROL::makePtr<vector>(ci_dim,0.0); // inequality multipliers

    ROL::Ptr<V> xo = ROL::makePtr<SV>(xo_ptr); 
    ROL::Ptr<V> xs = ROL::makePtr<SV>(xs_ptr);
    ROL::Ptr<V> xe = ROL::makePtr<SV>(xe_ptr);
    ROL::Ptr<V> xi = ROL::makePtr<SV>(xi_ptr);     

    ROL::RandomizeVector(*xo,left,right);
    ROL::RandomizeVector(*xs,left,right);
    ROL::RandomizeVector(*xe,left,right);
    ROL::RandomizeVector(*xi,left,right);

    ROL::Ptr<V> x = ROL::CreatePartitionedVector( xo, xs, xe, xi );
    

    // ----[ Full primal-dual direction vector ]------

    ROL::Ptr<vector> vo_ptr = ROL::makePtr<vector>(xo_dim,0.0); // opt
    ROL::Ptr<vector> vs_ptr = ROL::makePtr<vector>(ci_dim,0.0); // slack
    ROL::Ptr<vector> ve_ptr = ROL::makePtr<vector>(ce_dim,0.0); // equality multipliers
    ROL::Ptr<vector> vi_ptr = ROL::makePtr<vector>(ci_dim,0.0); // inequality multipliers
 
    ROL::Ptr<V> vo = ROL::makePtr<SV>(vo_ptr);
    ROL::Ptr<V> vs = ROL::makePtr<SV>(vs_ptr);
    ROL::Ptr<V> ve = ROL::makePtr<SV>(ve_ptr);
    ROL::Ptr<V> vi = ROL::makePtr<SV>(vi_ptr);     

    ROL::RandomizeVector(*vo,left,right);
    ROL::RandomizeVector(*vs,left,right);
    ROL::RandomizeVector(*ve,left,right);
    ROL::RandomizeVector(*vi,left,right);

    ROL::Ptr<V> v = ROL::CreatePartitionedVector( vo, vs, ve, vi );


    // ----[ Full primal-dual residual vector ]------

    ROL::Ptr<vector> ro_ptr = ROL::makePtr<vector>(xo_dim,0.0); // opt
    ROL::Ptr<vector> rs_ptr = ROL::makePtr<vector>(ci_dim,0.0); // slack
    ROL::Ptr<vector> re_ptr = ROL::makePtr<vector>(ce_dim,0.0); // equality multipliers
    ROL::Ptr<vector> ri_ptr = ROL::makePtr<vector>(ci_dim,0.0); // inequality multipliers
 
    ROL::Ptr<V> ro = ROL::makePtr<SV>(vo_ptr);
    ROL::Ptr<V> rs = ROL::makePtr<SV>(vs_ptr);
    ROL::Ptr<V> re = ROL::makePtr<SV>(ve_ptr);
    ROL::Ptr<V> ri = ROL::makePtr<SV>(vi_ptr);     

    ROL::RandomizeVector(*ro,left,right);
    ROL::RandomizeVector(*rs,left,right);
    ROL::RandomizeVector(*re,left,right);
    ROL::RandomizeVector(*ri,left,right);

    ROL::Ptr<V> r = ROL::CreatePartitionedVector( ro, rs, re, ri );

    // ----[ Primal-dual constraint ]-------

    ROL::Ptr<ROL::Objective<RealT> > obj_hs32 = 
      ROL::makePtr<ROL::ZOO::Objective_HS32<RealT>>();

    ROL::Ptr<ROL::EqualityConstraint<RealT> > eqcon_hs32 = 
      ROL::makePtr<ROL::ZOO::EqualityConstraint_HS32<RealT>>();
    
    ROL::Ptr<ROL::EqualityConstraint<RealT> > incon_hs32 = 
      ROL::makePtr<ROL::ZOO::InequalityConstraint_HS32<RealT>>();      


    *outStream << "Performing finite difference check on Primal-Dual KKT system" 
               << std::endl;

    using ROL::InteriorPoint::PrimalDualResidual;    

    PrimalDualResidual<RealT> con(obj_hs32,eqcon_hs32,incon_hs32, *x);
   
    con.checkApplyJacobian(*x,*v,*r,true,*outStream);  
    
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


