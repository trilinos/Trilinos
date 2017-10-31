#include "Rocket.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include <iostream>

int main( int argc, char* argv[] ) {

  using Teuchos::rcp; 
  using vector = std::vector<double>;

  auto parlist = rcp( new Teuchos::ParameterList() );
  Teuchos::updateParametersFromXmlFile("Rocket.xml",parlist.ptr()); 

  int     N  = parlist->get("Time Steps" ,      100   );  
  double  T  = parlist->get("Final Time",       20.0  );    
  double  g  = parlist->get("Gravity Constant", 9.8   );
  double  mr = parlist->get("Rocket Mass",      20.0  );
  double  mf = parlist->get("Fuel Mass",        100.0 );
  double  mu = parlist->get("Mass Penalty",     0.1   );
  double  ve = parlist->get("Exhaust Velocity", 1000.0);
  double  mt = mf+mr;    // Total mass
  double  dt = T/N;      // Time ste

  auto u_rcp = rcp( new vector(N) );
  auto u = rcp( new ROL::StdVector<double>(u_rcp) );
  auto l = u->dual().clone();

  auto z_rcp = rcp( new vector(N,mf/T) );
  auto z = rcp( new ROL::StdVector<double>(z_rcp) );


  // Trapezpoidal weights
  auto w_rcp = rcp( new vector(N,dt) );
  (*w_rcp)[0] *= 0.5; (*w_rcp)[N-1] *= 0.5;
  auto w = rcp( new ROL::StdVector<double>(w_rcp) );

  // Piecewise constant weights
  auto e_rcp = rcp( new vector(N,dt) );
  auto e = rcp( new ROL::StdVector<double>(e_rcp) );

  auto con = rcp( new Rocket::Constraint( N, T, mt, mf, ve, g ) );

  double tol = 1.e-7; // Needed for solve

  // Compute solution for constant burn rate
  con->update_2(*z);
  con->solve(*l, *u, *z, tol);
  //  u->print(std::cout); 
  double htarg = w->dot(*u);  // Final height 

  auto obj = rcp( new Rocket::Objective( N, T, mt, htarg, mu, w ) );
  auto robj = rcp( new ROL::Reduced_Objective_SimOpt<double>( obj, con, u, z, l ) );

  // Full space problem
  //  auto x = Teuchos::rcp( new ROL::Vector_SimOpt<double>(u,z) );
  //  ROL::OptimizationProblem<double> opt( obj, x, con, l );
  ROL::OptimizationProblem<double> opt( robj, z ); 
  ROL::OptimizationSolver<double> solver( opt, *parlist );
  solver.solve( std::cout );
  
  con->update_2(*z);
  con->solve(*l, *u, *z, tol);
  //  u->print(std::cout); 
  std::cout << "target height = " << htarg << 
               ", actual = " << w->dot(*u) << std::endl;
  std::cout << "Initial fuel mass   = " << mf << std::endl;
  std::cout << "Remaining fuel mass = " << mf-e->dot(*z) << std::endl;

  return 0;
}
