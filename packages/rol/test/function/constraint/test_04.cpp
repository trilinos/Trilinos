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

#include "ROL_Constraint_PinTSimOpt.hpp"

typedef double RealT;

int main(int argc, char* argv[]) 
{

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

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  try {
    double dt = 0.1;
    double tol = 1e-15;

    if(numRanks!=3) {
      throw std::logic_error("Three processors are required to run this test!");
    }

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

    ROL::Ptr<ROL::Constraint_TimeSimOpt<RealT>> step_constraint = ROL::makePtr<ODE_Constraint<RealT>>(dt,2.0*M_PI);

    std::vector<std::vector<int>> sizes = { {6, 0},    // int(number of steps/3), remainder 
                                            {6, 1},
                                            {6, 2} };
    for(auto size : sizes) {
      ROL::Ptr<ROL::Constraint_PinTSimOpt<RealT>> pint_constraint 
        = ROL::makePtr<ROL::Constraint_PinTSimOpt<RealT>>(MPI_COMM_WORLD, 
                                                          step_constraint,
                                                          size[0]*3+size[1], // this is a "hard" number!
                                                          uo_vec, std::vector<int>({-1,0}),  
                                                          z,      std::vector<int>({0}));

      std::stringstream ss ;
      ss << "\n  Testing size = { " << size[0] << ", " << size[1] << " }" << std::endl;
      if(myRank==0) {
        if(pint_constraint->stepStart()!=0) 
          throw std::logic_error("Rank 0 step start is incorrect " + std::to_string(pint_constraint->stepStart())+ss.str());
        if(pint_constraint->stepEnd()!=size[0]+(size[1]>0 ? 1 : 0)) 
          throw std::logic_error("Rank 0 step end is incorrect " + std::to_string(pint_constraint->stepEnd())+ss.str());
      }
      else if(myRank==1) {
        if(pint_constraint->stepStart()!=size[0]+(size[1]>0 ? 1 : 0)) 
          throw std::logic_error("Rank 1 step start is incorrect " + std::to_string(pint_constraint->stepStart())+ss.str());
        if(pint_constraint->stepEnd()!=2*size[0]+ size[1])
          throw std::logic_error("Rank 1 step end is incorrect " + std::to_string(pint_constraint->stepEnd())+ss.str());
      }
      else if(myRank==2) {
        if(pint_constraint->stepStart()!=2*size[0]+size[1])
          throw std::logic_error("Rank 2 step start is incorrect " + std::to_string(pint_constraint->stepStart())+ss.str());
        if(pint_constraint->stepEnd()!=3*size[0]+size[1]) 
          throw std::logic_error("Rank 2 step end is incorrect " + std::to_string(pint_constraint->stepEnd())+ss.str());
      }
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  int errors = std::abs(errorFlag);
  MPI_Allreduce(&errors,&errorFlag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
