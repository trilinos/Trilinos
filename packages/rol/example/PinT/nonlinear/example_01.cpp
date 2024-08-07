// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve a control problem governed by a semilinear
           parabolic equation,
           \f[
              \min_{u,z} \;\sum_{n=1}^N frac{\delta_t}{2} (u_{n-1}-1)^2 (u_n-1)^2 (z_n-1)^2
           \f]
           subject to \f$u_0=0\f$ and
           \f[
               e^{u_n - u_{n-1} - \delta_t z_n} = 1.
           \f]
*/

#include "dynamicConstraint.hpp"
#include "dynamicObjective.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include <iostream>
//#include <fenv.h>

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;
  using uint  = std::vector<RealT>::size_type;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  // *** Example body.
  int errorFlag = 0;
  try {
    // Parse input parameter list
    ROL::Ptr<ROL::ParameterList> pl = ROL::getParametersFromXmlFile("input_ex01.xml");
    bool derivCheck = pl->get("Derivative Check",         true); // Check derivatives.
    uint nt         = pl->get("Temporal Discretization",   100); // Set temporal discretization.
    RealT T         = pl->get("End Time",                  1.0); // Set end time.
    RealT dt        = T/(static_cast<RealT>(nt)-1.0);

    // Initialize objective function.
    ROL::Ptr<ROL::DynamicObjective<RealT>> dyn_obj
      = ROL::makePtr<Objective_Nonlinear<RealT>>(*pl);
    ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con
      = ROL::makePtr<Constraint_Nonlinear<RealT>>(*pl);

    // Create control vectors.
    ROL::Ptr<ROL::StdVector<RealT>>         zk = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    ROL::Ptr<ROL::PartitionedVector<RealT>>  z = ROL::PartitionedVector<RealT>::create(*zk, nt);

    // Create initial state vector.
    ROL::Ptr<ROL::StdVector<RealT>> u0 = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);

    // Create constraint vector.
    ROL::Ptr<ROL::StdVector<RealT>> ck = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);

    // Construct reduced dynamic objective
    std::vector<ROL::TimeStamp<RealT>> timeStamp(nt);
    for( uint k=0; k<nt; ++k ) {
      timeStamp.at(k).t.resize(2);
      timeStamp.at(k).t.at(0) = k*dt;
      timeStamp.at(k).t.at(1) = (k+1)*dt;
    }
    ROL::ParameterList &rpl = pl->sublist("Reduced Dynamic Objective");
    ROL::Ptr<ROL::ReducedDynamicObjective<RealT>> obj
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl);

    // Check derivatives for dynamic interface and reduced dynamic objective
    if (derivCheck) {
      ROL::Ptr<ROL::PartitionedVector<RealT>> dz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::PartitionedVector<RealT>> hz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::StdVector<RealT>> uo = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
      ROL::Ptr<ROL::StdVector<RealT>> un = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
      zk->randomize();
      z->randomize();
      dz->randomize();
      hz->randomize();
      uo->randomize();
      un->randomize();
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicObjectiveCheck<RealT>::check(*dyn_obj,validate,*uo,*un,*zk);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*zk);
      obj->checkGradient(*z,*dz,true,*outStream);
      obj->checkHessVec(*z,*dz,true,*outStream);
      obj->checkHessSym(*z,*dz,*hz,true,*outStream);
    }

    // Set up optimization problem
    ROL::Ptr<ROL::ParameterList> pl_rol = ROL::getParametersFromXmlFile("input_rol.xml");
    z->zero();
    //ROL::OptimizationProblem<RealT> problem(obj,z,z_bnd);
    ROL::OptimizationProblem<RealT> problem(obj,z);
    ROL::OptimizationSolver<RealT> solver(problem,*pl_rol);
    solver.solve(*outStream);
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

