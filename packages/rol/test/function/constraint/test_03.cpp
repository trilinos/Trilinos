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

#include "ROL_Constraint_TimeSimOpt.hpp"
#include "ROL_Vector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"

typedef double RealT;

template <typename Real>
class ODE_Constraint : public ROL::Constraint_TimeSimOpt<Real> {
private:
  typedef ROL::StdVector<Real> VectorType;

  const std::vector<Real> & getVector( const ROL::Vector<Real> & x ) {
    return *dynamic_cast<const VectorType&>(x).getVector();
  }

  std::vector<Real> & getVector( ROL::Vector<Real>  & x ) {
    return *dynamic_cast<VectorType&>(x).getVector();
  }

  Real timestep_; 
  Real omega_; 

public:
   
  ODE_Constraint(double dt,double omega) : timestep_(dt), omega_(omega) {}

  virtual void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &u_old,
             const ROL::Vector<Real> &u_new,
             const ROL::Vector<Real> &z,
             Real &tol) override
  {
    auto & c_data = getVector(c);
    auto & uo_data = getVector(u_old);
    auto & un_data = getVector(u_new);
    auto & z_data = getVector(z);

    // solving (u,v are states, z is control)
    // 
    //    u' = v, v'=-omega^2 * u + z
    //    u(0) = 0
    //
    //    u(t) = sin(omega*t)
    //
    // using backward euler
    
    c_data[0] = un_data[0]-uo_data[0] - timestep_*(un_data[1] + z_data[0]);
    c_data[1] = un_data[1]-uo_data[1] + timestep_*omega_*omega_*un_data[0];
  }

  virtual void solve(ROL::Vector<Real> &c,
                     const ROL::Vector<Real> &u_old, ROL::Vector<Real> &u_new,
                     const ROL::Vector<Real> &z,
                     Real &tol) override
  {
    auto & uo_data = getVector(u_old);
    auto & un_data = getVector(u_new);
    auto & z_data = getVector(z);

    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    Real rhs_0 = uo_data[0]+timestep_*z_data[0];
    Real rhs_1 = uo_data[1];

    un_data[0] = gamma*(1.0 * rhs_0 +   b * rhs_1);
    un_data[1] = gamma*( -a * rhs_0 + 1.0 * rhs_1);

    value(c,u_old,u_new,z,tol);
  }

  virtual void applyJacobian_1_old(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_old,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vo_data = getVector(v_old);

    jv_data[0] = -vo_data[0];
    jv_data[1] = -vo_data[1];
  }

  virtual void applyJacobian_1_new(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_new,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vn_data = getVector(v_new);

    jv_data[0] = vn_data[0] - timestep_*vn_data[1];
    jv_data[1] = vn_data[1] + timestep_*omega_*omega_*vn_data[0];

             // [      1,   -dt ]
             // [ dt*w*w,     1 ]
  }

  virtual void applyInverseJacobian_1_new(ROL::Vector<Real> &ijv,
                                          const ROL::Vector<Real> &v_new,
                                          const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                          const ROL::Vector<Real> &z,
                                          Real &tol) override {
    auto & ijv_data = getVector(ijv);
    auto & v_data = getVector(v_new); 
      
    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    ijv_data[0] = gamma*(1.0 * v_data[0] +   b * v_data[1]);
    ijv_data[1] = gamma*( -a * v_data[0] + 1.0 * v_data[1]);
  }

  virtual void applyJacobian_2(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_new,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vn_data = getVector(v_new);

    jv_data[0] = -timestep_*vn_data[0];
    jv_data[1] = 0.0;
  }

  virtual void applyAdjointJacobian_1_old(ROL::Vector<Real> &ajv_old,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {
    auto & ajv_data = getVector(ajv_old);
    auto & v_data = getVector(dualv);

    ajv_data[0] = -v_data[0];
    ajv_data[1] = -v_data[1];
  }

  virtual void applyAdjointJacobian_1_new(ROL::Vector<Real> &ajv_new,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {

    auto & ajv_data = getVector(ajv_new);
    auto & v_data = getVector(dualv);

    ajv_data[0] = v_data[0] + timestep_*omega_*omega_*v_data[1];
    ajv_data[1] = v_data[1] - timestep_*v_data[0];
  }

  virtual void applyAdjointJacobian_2_time(ROL::Vector<Real> &ajv,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {
    auto & ajv_data = getVector(ajv);
    auto & v_data = getVector(dualv);

    ajv_data[0] = -timestep_*v_data[0];
  }

  virtual void applyInverseAdjointJacobian_1_new(ROL::Vector<Real> &iajv,
                                                 const ROL::Vector<Real> &v_new,
                                                 const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                                 const ROL::Vector<Real> &z,
                                                 Real &tol) override {
    auto & iajv_data = getVector(iajv);
    auto & v_data = getVector(v_new); 
      
    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    iajv_data[0] = gamma*(1.0 * v_data[0] -   a * v_data[1]);
    iajv_data[1] = gamma*(  b * v_data[0] + 1.0 * v_data[1]);
  }
};

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
    std::vector<RealT> uo_data(2), un_data(2);
    PtrVector uo_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(uo_data));
    PtrVector un_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(un_data));
    PtrVector u = ROL::makePtr<ROL::PartitionedVector<RealT>>(std::vector<PtrVector>({un_vec,uo_vec}));
    CPtrVector cu = u;
    PtrVector v_u = u->clone();


    // allocate control vector
    std::vector<RealT> z_data(1);
    PtrVector z = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(z_data));
    CPtrVector cz = z;
    PtrVector v_z = z->clone();

    // allocate constraint vector
    std::vector<RealT> c_data(2);
    PtrVector c = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(c_data));
    CPtrVector cc = c;
    PtrVector jv = c->clone();
    PtrVector w = c->clone();
    PtrVector v_c = c->clone();

    ROL::Ptr<ROL::Constraint_TimeSimOpt<RealT>> constraint = ROL::makePtr<ODE_Constraint<RealT>>(dt,2.0*M_PI);

    ROL::RandomizeVector<RealT>(*u);
    ROL::RandomizeVector<RealT>(*z);
    ROL::RandomizeVector<RealT>(*v_u);
    ROL::RandomizeVector<RealT>(*v_z);
    ROL::RandomizeVector<RealT>(*v_c);
    ROL::RandomizeVector<RealT>(*jv);
    ROL::RandomizeVector<RealT>(*w);

    // check the solve
    *outStream << "Checking solve" << std::endl;
    constraint->checkSolve(*u,*z,*c,true,*outStream);

    // check the Jacobian_1
    *outStream << "Checking apply Jacobian 1" << std::endl;
    { 
      auto errors = constraint->checkApplyJacobian_1(*u,*z,*v_u,*jv,true,*outStream);
      if(errors[6][3] >= 1e-8)
        throw std::logic_error("Constraint apply jacobian 1 is incorrect");
    }

    // check the Jacobian_2
    *outStream << "Checking apply Jacobian 2" << std::endl;
    {
      auto errors = constraint->checkApplyJacobian_2(*u,*z,*v_z,*jv,true,*outStream);
      if(errors[6][3] >= 1e-8)
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

    // check inverses adjoint Jacobian_1_new
    *outStream << "Checking apply Jacobian 2 new" << std::endl;
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
