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

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"
#include "ROL_DynamicTrackingObjective.hpp"

#include "ROL_SerialConstraint.hpp"

#include "VdP_DynamicConstraint.hpp"


using RealT = double;

int main( int argc, char* argv[] ) {

  using ROL::Ptr;
  using ROL::makePtr;

  using vector    = std::vector<RealT>;
  using SV        = ROL::StdVector<RealT>;
  using PV        = ROL::PartitionedVector<RealT>;
  using size_type = typename PV::size_type;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );    
  int  errorFlag = 0;

  RealT T = 2.0;       // Total time
  size_type Nt = 100;  // Number of Time Steps


  auto uo = makePtr<SV>( makePtr<vector>(2) );
  auto un = makePtr<SV>( makePtr<vector>(2) );
  auto z  = makePtr<SV>( makePtr<vector>(1) );

  uo->randomize();
  un->randomize();
  z->randomize();
  
  // Tracking term is zero
  auto tracking = PV::create( *uo, Nt );
  tracking->zero();

  // Control regularization parameter is unity
  RealT alpha = 1.0;

  auto dyn_con = ROL::makePtr<VdP::DynamicConstraint<RealT>>();

  ROL::DynamicTrackingObjective<RealT> dyn_obj( tracking, alpha );

  ROL::ValidateFunction<RealT> validator( 1, 13, 20, 11, true, *outStream);

  // Check individual time step

  ROL::DynamicObjectiveCheck<RealT>::check( dyn_obj, validator, *uo, *un, *z, {
    "gradient_uo", "gradient_un", "gradient_z"  });

  ROL::DynamicConstraintCheck<RealT>::check( *dyn_con, validator, *uo, *un, *z, 
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
 
  auto timeStamps = ROL::TimeStamp<RealT>::make_uniform(0,T,{0.0,1.0},Nt);  

  *outStream << "\nChecking SerialConstraint:\n";

  ROL::SerialConstraint<RealT> serial_con( dyn_con, *uo, timeStamps );

  serial_con.checkApplyJacobian_1( *U, *Z, *VU, *C, true, *outStream );
  serial_con.checkApplyJacobian_2( *U, *Z, *VZ, *C, true, *outStream );

  serial_con.checkAdjointConsistencyJacobian_1( *W, *VU, *U, *Z, true, *outStream );
  serial_con.checkAdjointConsistencyJacobian_2( *W, *VU, *U, *Z, true, *outStream );



  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;
}
