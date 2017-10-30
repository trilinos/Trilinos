#include "Rocket_Simple.hpp"
#include "ROL_RandomVector.hpp"

#include <iostream>

int main( int argc, char* argv[] ) {

  using Teuchos::rcp;

  int      N = 100;
  double   T = 1.0;
  double   g = 9.8;
  double  mr = 5;
  double  mf = 20;
  double  mt = mf+mr;
  double  ve = 10;
 

  auto u_rcp = rcp( new std::vector<double>(N) );
  auto u = rcp( new ROL::StdVector<double>(u_rcp) );

  auto z_rcp = rcp( new std::vector<double>(N) );
  auto z = rcp( new ROL::StdVector<double>(z_rcp) );
  
  auto vu = u->clone();
  auto vz = z->clone();

  auto du = u->clone();
  auto dz = z->clone();

  auto c = u->clone();

  ROL::RandomizeVector( *c,  .01,  1.0 );
  ROL::RandomizeVector( *u,  .01,  1.0 );
  ROL::RandomizeVector( *z,  0.01, 1.0 );
  ROL::RandomizeVector( *vu, .01,  1.0 );
  ROL::RandomizeVector( *vz, 0.01, 1.0 );
  ROL::RandomizeVector( *du, .01,  1.0 );
  ROL::RandomizeVector( *dz, .01,  1.0 );

//  double tol = 1.e-7;

  auto con = rcp( new Rocket::Constraint( N, T, mt, mf, ve, g ) );

  con->update_2(*z);

  con->checkSolve(*u,*z,*c,true,std::cout);

  con->checkApplyJacobian_1(*u, *z, *du, *vu, true, std::cout); 

  con->checkInverseJacobian_1(*vu,*du,*u,*z,true,std::cout);

  con->checkAdjointConsistencyJacobian_1(*c, *u, *vu, *z, true, std::cout);

  con->checkInverseAdjointJacobian_1(*u, *z, *du, *vu, true, std::cout); 
  
  con->checkApplyJacobian_2(*u,*z,*vz,*vu,true, std::cout);

  con->checkAdjointConsistencyJacobian_2(*c, *z, *vu, *dz, true, std::cout);
  


//  u->print(std::cout);



  return 0;
}
