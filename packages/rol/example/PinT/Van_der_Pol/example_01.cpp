// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iomanip>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"
#include "ROL_DynamicTrackingObjective.hpp"
#include "ROL_SerialObjective.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "VdP_DynamicConstraint.hpp"

using RealT = double;

int main( int argc, char* argv[] ) {

  using namespace ROL;

  using vector    = std::vector<RealT>;
  using SV        = StdVector<RealT>;
  using PV        = PartitionedVector<RealT>;
  using size_type = typename PV::size_type;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  auto outStream = makeStreamPtr( std::cout, argc > 1 );    
  int  errorFlag = 0;

  auto VdP_params = getParametersFromXmlFile( "VdP_Parameters.xml" );

  RealT T = VdP_params->get("Total Time", 2.0); 
  size_type Nt = static_cast<size_type>( VdP_params->get("Number of Time Steps",100) ); 

  auto u_initial = makePtr<vector>(2);
  (*u_initial)[0] = VdP_params->get("Initial Position",1.0);
  (*u_initial)[1] = VdP_params->get("Initial Velocity",1.0);
  
  auto initialCondition = makePtr<SV>(u_initial);

  auto uo = makePtr<SV>( makePtr<vector>(2) );
  auto un = makePtr<SV>( makePtr<vector>(2) );
  auto z  = makePtr<SV>( makePtr<vector>(1) );

  uo->randomize();
  un->randomize();
  z->randomize();
  
  // Tracking term is zero
  auto tracking = PV::create( *uo, Nt );
  tracking->zero();

  auto dyn_con = makePtr<VdP::DynamicConstraint<RealT>>();

  ValidateFunction<RealT> validator( 1, 13, 20, 11, true, *outStream);

  *outStream << std::string(80,'-') << std::endl;
  *outStream << "\n\nChecking DynamicConstraint:\n\n";

  auto timeStamps = TimeStamp<RealT>::make_uniform(0,T,{0.0,1.0},Nt);  

  DynamicConstraintCheck<RealT>::check( *dyn_con, validator, *uo, *un, *z, timeStamps->at(1),
    { 
      "applyJacobian_uo",
      "applyJacobian_un",
      "applyJacobian_z",
      "applyAdjointJacobian_uo",
      "applyAdjointJacobian_un",
      "applyAdjointJacobian_z",
      "applyInverseJacobian_un",
      "applyInverseAdjointJacobian_un",
      "applyAdjointHessian_uo_uo",
      "applyAdjointHessian_uo_z",
      "applyAdjointHessian_un_un",
      "applyAdjointHessian_un_z",
      "applyAdjointHessian_z_uo",
      "applyAdjointHessian_z_un" 
  } );


  auto U  = PV::create( *uo, Nt );
  auto Z  = PV::create( *z,  Nt );
  auto C  = U->clone();
  auto W  = U->clone();
  auto VU = U->clone();
  auto VZ = Z->clone();
  
  U->randomize();
  Z->randomize();
  C->randomize();
  W->randomize();
  VU->randomize();
  VZ->randomize();

  auto serial_con = make_SerialConstraint( dyn_con, *initialCondition, timeStamps );

  *outStream << std::string(80,'-') << std::endl;
  *outStream << "\nChecking SerialConstraint:\n";

  *outStream << "\n\ncheckApplyJacobian_1\n\n";
  serial_con->checkApplyJacobian_1( *U, *Z, *VU, *C, true, *outStream );

  *outStream << "\n\ncheckApplyJacobian_2\n\n";
  serial_con->checkApplyJacobian_2( *U, *Z, *VZ, *C, true, *outStream );

  *outStream << "\n\ncheckApplyAdjointHessian_11\n\n";
  serial_con->checkApplyAdjointHessian_11( *U, *Z, *W, *VU, *C, true, *outStream );

  serial_con->checkAdjointConsistencyJacobian_1( *W, *VU, *U, *Z, true, *outStream );
  serial_con->checkAdjointConsistencyJacobian_2( *W, *VZ, *U, *Z, true, *outStream );

  serial_con->checkInverseJacobian_1( *C, *VU, *U, *Z, true, *outStream );
  serial_con->checkInverseAdjointJacobian_1( *C, *VU, *U, *Z, true, *outStream );
 
  serial_con->checkSolve( *U, *Z, *C, true, *outStream );

  // Target is zero state
  auto U0 = partition(U->clone()); 
  U0->zero();

  // Control regularization parameter.
  RealT alpha = VdP_params->get("Control Penalty",1.0);;
 
  auto dyn_obj    = make_DynamicTrackingObjective(U0,alpha);

  *outStream << std::string(80,'-') << std::endl;
  *outStream << "\n\nChecking DynamicObjective:\n\n";

  DynamicObjectiveCheck<RealT>::check( *dyn_obj, validator, *uo, *un, *z, timeStamps->at(1),
    {
      "gradient_un",
      "gradient_z",
      "hessVec_un_un",
      "hessVec_uo_uo",
      "hessVec_z_z",
    } );

  auto serial_obj = make_SerialObjective( dyn_obj, *initialCondition, timeStamps );

  *outStream << std::string(80,'-') << std::endl;
  *outStream << "\nChecking SerialObjective:\n";
  
  *outStream << "\n\ncheckGradient_1\n\n";
  serial_obj->checkGradient_1( *U, *Z, *VU, true, *outStream ); 

  *outStream << "\n\ncheckGradient_2\n\n";
  serial_obj->checkGradient_2( *U, *Z, *VZ, true, *outStream ); 

  *outStream << "\n\ncheckHessVec_11\n\n";
  serial_obj->checkHessVec_11( *U, *Z, *VU, true, *outStream ); 

  *outStream << "\n\ncheckHessVec_22\n\n";
  serial_obj->checkHessVec_22( *U, *Z, *VZ, true, *outStream ); 

  // Initial guess of unit state and unit amplitude control
  U->setScalar(1.0);
  Z->setScalar(1.0);
  
  auto x       = make_Vector_SimOpt( U, Z );
  auto problem = make_OptimizationProblem( serial_obj, x, serial_con, W ); 
  auto solver  = make_OptimizationSolver( problem, VdP_params );

  solver->solve(*outStream);
  
  *outStream << "Optimal Solution\n\n";
  *outStream << std::setw(13) << "Time Step (k)"
             << std::setw(12)  << "t_k"
             << std::setw(16) << "Position"
             << std::setw(16) << "Velocity" 
             << std::setw(20) << "Damping Coefficient" << std::endl;
  *outStream << std::string(100,'-') << std::endl;
  for( size_type k=0; k<U->numVectors(); ++k ) {
    auto& u_k = *(static_cast<SV&>((*U)[k]).getVector());
    auto& z_k = *(static_cast<SV&>((*Z)[k]).getVector());
    *outStream << std::setw(13) << k 
               << std::setw(12) << std::setprecision(3) << timeStamps->at(k).t[0] 
               << std::setw(16) << std::setprecision(8) << u_k[0] 
               << std::setw(16) << std::setprecision(8) << u_k[1] 
               << std::setw(20) << std::setprecision(8) << z_k[0] << std::endl;
  }



  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;
}
