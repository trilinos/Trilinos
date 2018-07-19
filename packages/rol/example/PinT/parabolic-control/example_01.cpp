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

/*! \file  example_01.cpp
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
#include <fenv.h>

int main(int argc, char *argv[]) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
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
    uint nx         = pl->get("Spatial Discretization",     64); // Set spatial discretization.
    uint nt         = pl->get("Temporal Discretization",   100); // Set temporal discretization.
    RealT T         = pl->get("End Time",                  1.0); // Set end time.
    RealT dt        = T/(static_cast<RealT>(nt)-1.0);

    // Initialize objective function.
    ROL::Ptr<ROL::DynamicObjective<RealT>> dyn_obj
      = ROL::makePtr<Objective_ParabolicControl<RealT>>(*pl);
    ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con
      = ROL::makePtr<Constraint_ParabolicControl<RealT>>(*pl);

    // Create control vectors.
    ROL::Ptr<ROL::StdVector<RealT>>         zk = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx+2));
    ROL::Ptr<ROL::PartitionedVector<RealT>>  z = ROL::PartitionedVector<RealT>::create(*zk, nt);

    // Create initial state vector.
    ROL::Ptr<ROL::StdVector<RealT>> u0 = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx,0.0));

    // Create constraint vector.
    ROL::Ptr<ROL::StdVector<RealT>> ck = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx,0.0));

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
      ROL::Ptr<ROL::StdVector<RealT>> uo = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx));
      ROL::Ptr<ROL::StdVector<RealT>> un = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx));
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

