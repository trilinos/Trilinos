// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_parabolic.cpp
    \brief Shows how to solve a control problem governed by a semilinear
           parabolic equation,
           \f[
              \min_{u,z} \;\frac{1}{2} \int_0^T\int_0^1 (u(t,x)-\bar{u}(x))^2\,\mathrm{d}x\mathrm{d}t
                         \frac{\alpha}{2} \int_0^T z(t)^2\,\mathrm{d}t
           \f]
           subject to
           \f[
               u_t(t,x) - u_{xx}(t,x) + f(u(t,x),x) = z(t,x) \quad t\in (0,T], \; x\in (0,1)
           \f]
           with boundary conditions
           \f[
               u_x(t,0) = 0, \; u_x(t,1) = 0
           \f]
           and initial condition
           \f[
               u(0,x) = 0.
           \f]
*/

#include "example_parabolic_thyravec.hpp"

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
    ROL::Ptr<ROL::ParameterList> pl = ROL::getParametersFromXmlFile("example_parabolic.xml");
    ROL::ParameterList pl_prb = pl->sublist("Parabolic");
    bool derivCheck = pl_prb.get("Derivative Check",        true); // Check derivatives.
    uint nx         = pl_prb.get("Spatial Discretization",    64); // Set spatial discretization.
    uint nt         = pl_prb.get("Temporal Discretization",  100); // Set temporal discretization.
    RealT T         = pl_prb.get("End Time",                 1.0); // Set end time.
    RealT dt        = T/(static_cast<RealT>(nt)-1.0);

    // Initialize objective function.
    ROL::Ptr<ROL::DynamicObjective<RealT>> dyn_obj
      = ROL::makePtr<Objective_ParabolicControl<RealT>>(pl_prb);
    ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con
      = ROL::makePtr<Constraint_ParabolicControl<RealT>>(pl_prb);

    // Define vector spaces.
    ROL::Ptr<Thyra::VectorSpaceBase<RealT> > control_space = Thyra::defaultSpmdVectorSpace<RealT>(nx+2);
    ROL::Ptr<Thyra::VectorSpaceBase<RealT> > state_space = Thyra::defaultSpmdVectorSpace<RealT>(nx);

    // Create control vectors.
    ROL::Ptr<ROL::ThyraVector<RealT>>        zk = ROL::makePtr<ROL::ThyraVector<RealT>>(Thyra::createMember<RealT>(control_space));
    Thyra::put_scalar(0.0, (zk->getVector().ptr()));
    ROL::Ptr<ROL::PartitionedVector<RealT>>  z = ROL::PartitionedVector<RealT>::create(*zk, nt);

    // Create initial state vector.
    ROL::Ptr<ROL::ThyraVector<RealT>>        u0 = ROL::makePtr<ROL::ThyraVector<RealT>>(Thyra::createMember<RealT>(state_space));
    Thyra::put_scalar(0.0, (u0->getVector().ptr()));

    // Create constraint vector.
    ROL::Ptr<ROL::ThyraVector<RealT>>        ck = ROL::makePtr<ROL::ThyraVector<RealT>>(Thyra::createMember<RealT>(state_space));
    Thyra::put_scalar(0.0, (ck->getVector().ptr()));

    // Construct reduced dynamic objective
    std::vector<ROL::TimeStamp<RealT>> timeStamp(nt);
    for( uint k=0; k<nt; ++k ) {
      timeStamp.at(k).t.resize(2);
      timeStamp.at(k).t.at(0) = k*dt;
      timeStamp.at(k).t.at(1) = (k+1)*dt;
    }
    ROL::ParameterList &pl_reddyn = pl->sublist("Reduced Dynamic Objective");
    ROL::Ptr<ROL::ReducedDynamicObjective<RealT>> obj
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, pl_reddyn);

    // Check derivatives for dynamic interface and reduced dynamic objective
    if (derivCheck) {
      ROL::Ptr<ROL::PartitionedVector<RealT>> dz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::PartitionedVector<RealT>> hz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::ThyraVector<RealT>> uo = ROL::makePtr<ROL::ThyraVector<RealT>>(Thyra::createMember<RealT>(state_space));
      Thyra::put_scalar(0.0, (uo->getVector().ptr()));
      ROL::Ptr<ROL::ThyraVector<RealT>> un = ROL::makePtr<ROL::ThyraVector<RealT>>(Thyra::createMember<RealT>(state_space));
      Thyra::put_scalar(0.0, (un->getVector().ptr()));
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
    ROL::ParameterList pl_rol = pl->sublist("ROL");
    z->zero();
    ROL::OptimizationProblem<RealT> problem(obj,z);
    ROL::OptimizationSolver<RealT> solver(problem,pl_rol);
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

