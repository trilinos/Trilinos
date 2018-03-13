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

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ODEConstraint_TimeSimOpt.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;
  typedef ROL::Ptr<const ROL::Vector<RealT>> CPtrVector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    double dt = 0.1;
    double tol = 1e-15;

    // allocate state vector
    std::vector<RealT> uo_data(3), un_data(3);
    uo_data[0] = 0.0; uo_data[1] = 2.0*M_PI; uo_data[2] = 1.0;
    un_data[0] = 0.0; un_data[1] = 0.0;      un_data[2] = 0.1;
    PtrVector uo_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(uo_data));
    PtrVector un_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(un_data));
    PtrVector u = ROL::makePtr<ROL::PartitionedVector<RealT>>(std::vector<PtrVector>({uo_vec,un_vec}));
    CPtrVector cu = u;
    PtrVector v_u = u->clone();
    PtrVector hv_u = u->clone();

    u->print(std::cout);

    // allocate control vector
    std::vector<RealT> z_data(1);
    PtrVector z = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(z_data));
    CPtrVector cz = z;
    PtrVector v_z = z->clone();
    PtrVector hv_z = z->clone();

    // allocate constraint vector
    std::vector<RealT> c_data(3);
    PtrVector c = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(c_data));
    CPtrVector cc  = c;
    PtrVector jv   = c->clone();
    PtrVector w    = c->clone();
    PtrVector w_c    = c->clone();
    PtrVector v_c  = c->clone();

    ROL::Ptr<ROL::Constraint_TimeSimOpt<RealT>> constraint = ROL::makePtr<ODE_Constraint<RealT>>(dt,2.0*M_PI);

    // ROL::RandomizeVector<RealT>(*u);
    ROL::RandomizeVector<RealT>(*z);
    ROL::RandomizeVector<RealT>(*v_u);
    ROL::RandomizeVector<RealT>(*v_z);
    ROL::RandomizeVector<RealT>(*v_c);
    ROL::RandomizeVector<RealT>(*jv);
    ROL::RandomizeVector<RealT>(*hv_u);
    ROL::RandomizeVector<RealT>(*hv_z);
    ROL::RandomizeVector<RealT>(*w);
    ROL::RandomizeVector<RealT>(*w_c);

    // check the solve
    *outStream << "Checking solve" << std::endl;
    constraint->checkSolve(*u,*z,*c,true,*outStream);

    // check the Jacobian_1
    *outStream << "Checking apply Jacobian 1" << std::endl;
    { 
      auto errors = constraint->checkApplyJacobian_1(*u,*z,*v_u,*jv,true,*outStream);
      if(errors[6][3] >= 1e-6)
        throw std::logic_error("Constraint apply jacobian 1 is incorrect");
    }

    // check the Jacobian_2
    *outStream << "Checking apply Jacobian 2" << std::endl;
    {
      auto errors = constraint->checkApplyJacobian_2(*u,*z,*v_z,*jv,true,*outStream);
      if(errors[6][3] >= 1e-6)
        throw std::logic_error("Constraint apply jacobian 2 is incorrect");
    }

    // check inverses Jacobian_1_new
    *outStream << "Checking apply Jacobian 2 new" << std::endl;
    {
      ROL::Vector<RealT> & v_un = *dynamic_cast<ROL::PartitionedVector<RealT>&>(*v_u).get(1);
      dynamic_cast<ROL::PartitionedVector<RealT>&>(*v_u).get(0)->scale(0.0);


      constraint->applyJacobian_1(*jv,*v_u,*u,*z,tol);

      PtrVector ijv = v_un.clone();
      constraint->applyInverseJacobian_1_new(*ijv,*jv,*uo_vec,*un_vec,
                                                      *z,tol);

      ijv->axpy(-1.0,v_un);

      *outStream << "Inverse Jacobian_1_new error norm = " << ijv->norm() << std::endl;
      if(ijv->norm() >= tol)
        throw std::logic_error("Inverse Jacobian_1_new not checking out");
    }

    // check the Adjoint Jacobian_1
    *outStream << "Checking apply Adjoint Jacobian 1" << std::endl;
    {
      auto error = constraint->checkAdjointConsistencyJacobian_1(*w,*v_u,*u,*z,true,*outStream);
      if(error >= 1e-8)
        throw std::logic_error("Constraint apply adjoint jacobian 1 is incorrect");
    }

    // check the Adjoint Jacobian_1
    *outStream << "Checking apply Adjoint Jacobian 2" << std::endl;
    {
      auto error = constraint->checkAdjointConsistencyJacobian_2(*w,*v_z,*u,*z,true,*outStream);
      if(error >= 1e-8)
        throw std::logic_error("Constraint apply adjoint jacobian 2 is incorrect");
    }

    // check the Adjoint Hessian_11
    *outStream << "Checking apply Adjoint Hessian_11" << std::endl;
    {
      auto errors = constraint->checkApplyAdjointHessian_11(*u,*z,*w_c,*v_u,*hv_u,true,*outStream);
      if(errors[6][3] >= 1e-5)
        throw std::logic_error("Constraint apply Adjoint Hessian 11 is incorrect");
    }

    // check the Adjoint Hessian_12
    *outStream << "Checking apply Adjoint Hessian_12" << std::endl;
    {
      auto errors = constraint->checkApplyAdjointHessian_12(*u,*z,*w_c,*v_u,*hv_z,true,*outStream);
      if(errors[6][3] >= 1e-5)
        throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
    }

    // check the Adjoint Hessian_21
    *outStream << "Checking apply Adjoint Hessian_21" << std::endl;
    {
      auto errors = constraint->checkApplyAdjointHessian_21(*u,*z,*w_c,*v_z,*hv_u,true,*outStream);
      if(errors[6][3] >= 1e-5)
        throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
    }

    // check the Adjoint Hessian_22
    *outStream << "Checking apply Adjoint Hessian_22" << std::endl;
    {
      auto errors = constraint->checkApplyAdjointHessian_22(*u,*z,*w_c,*v_z,*hv_z,true,*outStream);
      if(errors[6][3] >= 1e-5)
        throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
    }

    // check inverses adjoint Jacobian_1_new
    *outStream << "Checking apply Inverse Adjoint Jacobian 2 new" << std::endl;
    {
      dynamic_cast<ROL::PartitionedVector<RealT>&>(*v_u).get(0)->scale(0.0);

      PtrVector jv_u = un_vec->clone();

      constraint->applyAdjointJacobian_1_new(*jv,*v_c,
                                         *uo_vec,*un_vec,
                                         *z,tol);

      PtrVector iajv = c->clone();
      constraint->applyInverseAdjointJacobian_1_new(*iajv,*jv,
                                                    *uo_vec,*un_vec,
                                                    *z,tol);

      iajv->axpy(-1.0,*v_c);

      *outStream << "Inverse Adjoint Jacobian_1_new error norm = " << iajv->norm() << std::endl;
      if(iajv->norm() >= tol)
        throw std::logic_error("Inverse Adjoint Jacobian_1_new not checking out");
    }

    {
      ROL::RandomizeVector<RealT>(*u);
      ROL::RandomizeVector<RealT>(*z);

      *outStream << "Pre solve: solution =\n";
      u->print(*outStream);

      // using randomized vectors, print value
      constraint->value(*c,*cu,*cz,tol);
      *outStream << "Pre solve: constraint =\n";
      c->print(*outStream);

      // using randomized vectors, print value
      constraint->solve(*c,*u,*z,tol);
      *outStream << "Post solve: constraint =\n";
      c->print(*outStream);
      *outStream << "Post solve: solution =\n";
      u->print(*outStream);

      *outStream << "Post solve constraint norm = " << c->norm() << std::endl;
      if(c->norm() >= tol)
        throw std::logic_error("Constraint required accuracy not reached");
    }

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
