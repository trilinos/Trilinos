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

#include "Rocket.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

int main( int argc, char* argv[] ) {

   
  using vector = std::vector<double>;

  auto parlist = ROL::getParametersFromXmlFile( "Rocket.xml" ); 
  

  int     N  = parlist->get("Time Steps" ,      100   );  
  double  T  = parlist->get("Final Time",       20.0  );    
  double  g  = parlist->get("Gravity Constant", 9.8   );
  double  mr = parlist->get("Rocket Mass",      20.0  );
  double  mf = parlist->get("Fuel Mass",        100.0 );
  double  mu = parlist->get("Mass Penalty",     0.1   );
  double  ve = parlist->get("Exhaust Velocity", 1000.0);
  double  mt = mf+mr;    // Total mass
  double  dt = T/N;      // Time ste

  auto u_ptr = ROL::makePtr<vector>(N);
  auto u = ROL::makePtr<ROL::StdVector<double>>(u_ptr);
  auto l = u->dual().clone();

  auto z_ptr = ROL::makePtr<vector>(N,mf/T);
  auto z = ROL::makePtr<ROL::StdVector<double>>(z_ptr);


  // Trapezpoidal weights
  auto w_ptr = ROL::makePtr<vector>(N,dt);
  (*w_ptr)[0] *= 0.5; (*w_ptr)[N-1] *= 0.5;
  auto w = ROL::makePtr<ROL::StdVector<double>>(w_ptr);

  // Piecewise constant weights
  auto e_ptr = ROL::makePtr<vector>(N,dt);
  auto e = ROL::makePtr<ROL::StdVector<double>>(e_ptr);

  auto con = ROL::makePtr<Rocket::Constraint>( N, T, mt, mf, ve, g );

  double tol = 1.e-7; // Needed for solve

  // Compute solution for constant burn rate
  con->update_2(*z);
  con->solve(*l, *u, *z, tol);

  auto ucb_ptr = ROL::makePtr<vector>(N);
  auto ucb = ROL::makePtr<ROL::StdVector<double>>(ucb_ptr);
  ucb->set(*u);

  //  u->print(std::cout); 
  double htarg = w->dot(*u);  // Final height 

  auto obj  = ROL::makePtr<Rocket::Objective>( N, T, mt, htarg, mu, w );
  auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<double>>( obj, con, u, z, l );

  // Full space problem
  //  auto x = ROL::makePtr<ROL::Vector_SimOpt<double>>(u,z);
  //  ROL::OptimizationProblem<double> opt( obj, x, con, l );
  ROL::OptimizationProblem<double> opt( robj, z ); 
  ROL::OptimizationSolver<double> solver( opt, *parlist );
  solver.solve( std::cout );
  con->update_2(*z);
  con->solve(*l, *u, *z, tol);

  std::ofstream os("Rocket.txt");
  if( os.is_open() ) {
  
    os << std::setw(12) << "Time" 
       << std::setw(18) << "Burn_Rate_(const)" 
       << std::setw(16) << "Height_(const)" 
       << std::setw(16) << "Burn_Rate_(opt)" 
       << std::setw(16) << "Height_(opt)" << std::endl;

    vector hcb(N,0);
    vector hopt(N,0); 
    
    os << std::setw(12) << 0
       << std::setw(18) << 0
       << std::setw(16) << 0
       << std::setw(16) << 0
       << std::setw(16) << 0 << std::endl;

    os << std::setw(12) << 0
       << std::setw(18) << mf/T
       << std::setw(16) << 0
       << std::setw(16) << (*z_ptr)[0]
       << std::setw(16) << 0 << std::endl;

    hcb[0]  = 0.5*dt*(*ucb_ptr)[0];
    hopt[0] = 0.5*dt*(*u_ptr)[0];

    for( int i=1; i<N; ++i ) {
      hcb[i]  = hcb[i-1]  + 0.5*dt*( (*ucb_ptr)[i] + (*ucb_ptr)[i-1] );
      hopt[i] = hopt[i-1] + 0.5*dt*( (*u_ptr)[i]   + (*u_ptr)[i-1] );

      os << std::setw(12) << i*dt
         << std::setw(18) << mf/T
         << std::setw(16) << hcb[i]
         << std::setw(16) << (*z_ptr)[i-1]
         << std::setw(16) << hopt[i] << std::endl;

      os << std::setw(12) << i*dt
         << std::setw(18) << mf/T
         << std::setw(16) << hcb[i]
         << std::setw(16) << (*z_ptr)[i]
         << std::setw(16) << hopt[i] << std::endl;
    }
  }
  return 0;
}
