/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <math.h>

#include <Kokkos_Value.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>

#include <Kokkos_Host.hpp>

#include <Kokkos_Host_macros.hpp>
#include <Element.hpp>
#include <CRSMatrixGatherFill.hpp>
#include <Dirichlet.hpp>
#include <CG_Solve.hpp>
#include <driver.hpp>
#include <Kokkos_Clear_macros.hpp>

namespace Test{

void test_Host(int beg, int end, int runs,int threads)
{
  if ( 0 < threads ) {
    Kokkos::Host::initialize( Kokkos::Host::SetThreadCount( threads ) );
    std::cout << std::endl << "\"Host with manually set threads = \" , " << threads << std::endl ;
    MiniImplTherm< float , Kokkos::Host >::driver( "Host-float" , beg , end , runs );
    MiniImplTherm< double, Kokkos::Host >::driver( "Host-double" , beg , end , runs );
    Kokkos::Host::finalize();
  }
  else {
    Kokkos::Host::initialize( Kokkos::Host::DetectAndUseCores() );
    threads = Kokkos::Host::detect_core_count();
    std::cout << std::endl << "\"Host with detected sequential threads = \" , " << threads << std::endl ;
    MiniImplTherm< float , Kokkos::Host >::driver( "Host-float" , beg , end , runs );
    MiniImplTherm< double, Kokkos::Host >::driver( "Host-double" , beg , end , runs );
    Kokkos::Host::finalize();
  }
} //test_host

}// namespace


