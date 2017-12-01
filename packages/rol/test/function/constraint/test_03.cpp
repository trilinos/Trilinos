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

  void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &u_old,
             const ROL::Vector<Real> &u_new,
             const ROL::Vector<Real> &z_old,
             const ROL::Vector<Real> &z_new,
             Real &tol) override
  {
    auto c_data = getVector(c);
    auto uo_data = getVector(u_old);
    auto un_data = getVector(u_new);

    // solving:
    //    u' = v, v'=-omega^2 * u
    //    u(0) = 0
    //
    //    u(t) = sin(omega*t)
    //
    //    using backward euler
    
    c_data[0] = un_data[0]-uo_data[0] - timestep_*un_data[0];
    c_data[1] = un_data[1]-uo_data[1] + timestep_*omega_*omega_*un_data[0];

    std::cout << "Value = " 
              << c_data[0] << " " << un_data[0]<< " " << uo_data[0] << " " << un_data[0] << std::endl;
  }

  virtual void solve(ROL::Vector<Real> &c,
                     const ROL::Vector<Real> &u_old,
                     ROL::Vector<Real> &u_new,
                     const ROL::Vector<Real> &z_old,
                     const ROL::Vector<Real> &z_new,
                     Real &tol) override
  {
    auto uo_data = getVector(u_old);
    auto un_data = getVector(u_new);

    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    un_data[0] = gamma*(1.0 * uo_data[0] +   b * uo_data[1]);
    un_data[1] = gamma*( -a * uo_data[0] + 1.0 * uo_data[1]);

    value(c,u_old,u_new,z_old,z_new,tol);
  }
};

int main(int argc, char* argv[]) {

  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;

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
    double tol = 1e-16;

    // allocate state vector
    std::vector<RealT> uo_data(2), un_data(2);
    PtrVector uo_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(uo_data));
    PtrVector un_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(un_data));
    PtrVector u = ROL::makePtr<ROL::PartitionedVector<RealT>>(std::vector<PtrVector>({un_vec,uo_vec}));

    // allocate control vector
    std::vector<RealT> zo_data(2), zn_data(2);
    PtrVector zo_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(zo_data));
    PtrVector zn_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(zn_data));
    PtrVector z = ROL::makePtr<ROL::PartitionedVector<RealT>>(std::vector<PtrVector>({zn_vec,zo_vec}));

    // allocate constraint vector
    std::vector<RealT> c_data(2);
    PtrVector c = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(c_data));

    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> constraint = ROL::makePtr<ODE_Constraint<RealT>>(dt,2.0*M_PI);

    ROL::RandomizeVector<RealT>(*u);

    // check the solve
    // constraint->checkSolve(*u,*z,*c,*outStream);

    ROL::RandomizeVector<RealT>(*u);

    *outStream << "Pre solve: solution =\n";
    u->print(*outStream);

    // using randomized vectors, print value
    constraint->value(*c,*u,*z,tol);
    *outStream << "Pre solve: constrint =\n";
    c->print(*outStream);

    // using randomized vectors, print value
    constraint->solve(*c,*u,*z,tol);
    *outStream << "Post solve: constrint =\n";
    c->print(*outStream);
    *outStream << "Post solve: solution =\n";
    u->print(*outStream);

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
