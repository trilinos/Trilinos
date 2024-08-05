// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Bounds.hpp"
#include "Tanks_DynamicConstraint.hpp"
#include "Tanks_DynamicObjective.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include <iostream>

int main( int argc, char* argv[] ) {
  using RealT     = double;
  using size_type = std::vector<RealT>::size_type;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );
  int errorFlag  = 0;

  try {     // *** Example body.
    // Parse input parameter list
    ROL::Ptr<ROL::ParameterList> pl_ptr  = ROL::getParametersFromXmlFile("parameters_ex04.xml");
//    RealT height    = pl_ptr->get("Height of Tank",       10.0);  // Unused
//    RealT Qin00     = pl_ptr->get("Corner Inflow",       100.0);  // Unused
    RealT h_init    = pl_ptr->get("Initial Fluid Level",   2.0);
    RealT T         = pl_ptr->get("Total Time",           20.0);
    size_type Nt    = static_cast<size_type>(pl_ptr->get("Number of Time Stamps",100));
    size_type nrows = static_cast<size_type>(pl_ptr->get("Number of Rows"   ,3));
    size_type ncols = static_cast<size_type>(pl_ptr->get("Number of Columns",3));
    RealT dt = T/Nt;

    // Construct dynamic constraint and objective function
    ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con = Tanks::DynamicConstraint<RealT>::create(*pl_ptr);
    ROL::Ptr<ROL::DynamicObjective<RealT>>  dyn_obj = ROL::makePtr<Tanks::DynamicObjective<RealT>>(*pl_ptr);

    // Create control vectors
    ROL::Ptr<Tanks::ControlVector<RealT>>  zk  = Tanks::ControlVector<RealT>::create(*pl_ptr, "Control");
    ROL::Ptr<ROL::PartitionedVector<RealT>> z  = ROL::PartitionedVector<RealT>::create(*zk, Nt);
    ROL::Ptr<ROL::PartitionedVector<RealT>> dz = ROL::PartitionedVector<RealT>::create(*zk, Nt);
    ROL::Ptr<ROL::PartitionedVector<RealT>> hz = ROL::PartitionedVector<RealT>::create(*zk, Nt);
    for( size_type i=0; i<Nt; ++i ) {
      z->get(i)->randomize();
      dz->get(i)->randomize();
      hz->get(i)->randomize();
    }

    // Create control bounds
    //ROL::Ptr<Tanks::ControlVector<RealT>>   zk_lo = zk->clone( "Control Lower Bound (zk_lo)" );
    //ROL::Ptr<ROL::PartitionedVector<RealT>> z_lo  = ROL::PartitionedVector<RealT>::create(*zk_lo, Nt);
    //z_lo->zero();
    //ROL::Ptr<ROL::Bounds<RealT>>            z_bnd = ROL::makePtr<ROL::Bounds<RealT>>( *z_lo );

    // Create initial state vector
    ROL::Ptr<Tanks::StateVector<RealT>> u0 = Tanks::StateVector<RealT>::create(*pl_ptr, "Initial State");
    ROL::Ptr<Tanks::StateVector<RealT>> uo = Tanks::StateVector<RealT>::create(*pl_ptr, "Old State");
    ROL::Ptr<Tanks::StateVector<RealT>> un = Tanks::StateVector<RealT>::create(*pl_ptr, "New State");
    for( size_type i=0; i<nrows; ++i ) {
      for( size_type j=0; j<ncols; ++j ) {
        u0->h(i,j) = h_init;
      }
    }

    // Check dynamic objective derivatives
    uo->randomize();
    un->randomize();
    zk->randomize();
    ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
    ROL::DynamicObjectiveCheck<RealT>::check(*dyn_obj,validate,*uo,*un,*zk);
    ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*zk);

    // Create constraint vector
    ROL::Ptr<Tanks::StateVector<RealT>> ck = Tanks::StateVector<RealT>::create(*pl_ptr, "Constraint");

    // Construct reduced dynamic objective
    std::vector<ROL::TimeStamp<RealT>> timeStamp(Nt);
    for( size_type k=0; k<Nt; ++k ) {
      timeStamp.at(k).t.resize(2);
      timeStamp.at(k).t.at(0) = k*dt;
      timeStamp.at(k).t.at(1) = (k+1)*dt;
    }
    ROL::ParameterList &rpl = pl_ptr->sublist("Reduced Dynamic Objective");
    ROL::Ptr<ROL::ReducedDynamicObjective<RealT>> obj
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl);

    // Check gradient of reduced dynamic objective
    obj->checkGradient(*z,*dz,true,*outStream);
    obj->checkHessVec(*z,*dz,true,*outStream);
    obj->checkHessSym(*z,*dz,*hz,true,*outStream);

    // Set up optimization problem
    ROL::Ptr<ROL::ParameterList> rol_ptr  = ROL::getParametersFromXmlFile("rol-parameters.xml");
    z->zero();
    //ROL::OptimizationProblem<RealT> problem(obj,z,z_bnd);
    ROL::OptimizationProblem<RealT> problem(obj,z);
    ROL::OptimizationSolver<RealT> solver(problem,*rol_ptr);
    solver.solve(*outStream);
  }
  catch (std::logic_error &err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
