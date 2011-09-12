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
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceTBB.hpp>
#include <Kokkos_DeviceTBB_ValueView.hpp>
#include <Kokkos_DeviceTBB_MultiVectorView.hpp>
#include <Kokkos_DeviceTBB_MDArrayView.hpp>
#include <Kokkos_DeviceTBB_ParallelFor.hpp>
#include <Kokkos_DeviceTBB_ParallelReduce.hpp>

#include <Kokkos_DeviceTBB_macros.hpp>
#include <assemble.hpp>
#include <CRSMatrixGatherFill.hpp>
#include <Dirichlet.hpp>
#include <CG_Solve.hpp>
#include <driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test{

  void test_TBB(int beg, int end , int runs,  int threads){

    std::cout << "Kokkos TBB - " << threads << " threads: " << std::endl;
    std::cout<<"Size , Mesh , Assemble , Gather , Dirichlet , Iteration , CG Time , Total"<<std::endl;
    Kokkos::DeviceTBB::initialize( threads );


    for(int i = beg ; i < end; i+=2)
     {
      double times[10], mins[10];
      for(int j = 0; j < runs; j++){

        run_kernel<Kokkos::DeviceTBB>(i,i,i,times);
      
        if(j == 0)
        {
          mins[0] = times[0];
          mins[1] = times[1];
          mins[2] = times[2];
          mins[3] = times[3];
          mins[4] = times[4];
          mins[5] = times[5];
          mins[6] = times[6];
        }
        for(int k = 0 ; k < 7 ; k++)
        {
          if(times[k] < mins[k])
            mins[k] = times[k];
        }
      }
      std::cout  <<i*i*i;
      for(int m = 0 ; m < 7; m++)
      {
        std::cout<<", "<<(times[m] * 1000)<<"ms";
      }
      std::cout<<std::endl;
    }
    Kokkos::DeviceTBB::finalize();

  }//test_TBB

}// namespace


