/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */


#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceHost_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//------------------------------------------------------------------------

namespace Test {

void run_test_host_hexgrad( int beg , int end )
{ Test::run_test_hexgrad< Kokkos::DeviceHost>( beg , end ); }

void run_test_host_gramschmidt( int beg , int end )
{ Test::run_test_gramschmidt< Kokkos::DeviceHost>( beg , end ); }

void run_test_tpi_hexgrad(int,int);
void run_test_tpi_gramschmidt(int,int);

void run_test_cuda_hexgrad(int,int);
void run_test_cuda_gramschmidt(int,int);

void run_test_tbb_hexgrad(int,int);
void run_test_tbb_gramschmidt(int,int);

void run_test_ferry_hexgrad(int,int);
void run_test_ferry_gramschmidt(int,int);

}

int main( int argc , char ** argv )
{
	Test::run_test_host_hexgrad( 10 , 20 );
	Test::run_test_tpi_hexgrad(  10 , 20 );
 	Test::run_test_cuda_hexgrad( 10 , 20 );
 	Test::run_test_tbb_hexgrad(  10 , 20 );
 	Test::run_test_ferry_hexgrad( 10 , 20 );
 
  	Test::run_test_host_gramschmidt( 10 , 20 );
  	Test::run_test_tpi_gramschmidt(  10 , 20 );
  	Test::run_test_cuda_gramschmidt( 10 , 20 );
 	Test::run_test_tbb_gramschmidt( 10 , 20 );
 	Test::run_test_ferry_gramschmidt(10 , 15 );

  return 0 ;
}

