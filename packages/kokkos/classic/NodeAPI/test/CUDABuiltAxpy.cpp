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

#include "CUDABuiltAxpy.hpp"
#include "Kokkos_ThrustGPUNode.hpp"

#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
double mytimer(void)
{
  struct timeval tp;
  static long start=0, startu;
  if (!start)
    {
      gettimeofday(&tp, NULL);
      start = tp.tv_sec;
      startu = tp.tv_usec;
      return(0.0);
    }
  gettimeofday(&tp, NULL);
  return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
}
double t0;
#define TICK()  t0 = mytimer() // Use TICK and TOCK to time a code section
#define TOCK(t) t += mytimer() - t0

void initializeVectors(float *u, float *v);
void outputStats(double time4);

////////////////////////////////////////////////////////////////////////////////
void axpyGPU(float *y, const float*x, float alpha, int N)
{
  Teuchos::ParameterList pl;
  pl.set("Device Number",0);
  pl.set("Verbose",1);

  typedef Kokkos::ThrustGPUNode GPUNode;

  Teuchos::RCP<GPUNode> node = Teuchos::rcp(new GPUNode(pl));
  Teuchos::ArrayRCP<float> xbuf_dev = node->allocBuffer<float>(N);
  Teuchos::ArrayRCP<float> ybuf_dev = node->allocBuffer<float>(N);
  {
    Teuchos::ArrayRCP<float> xbuf_host = node->viewBufferNonConst<float>(Kokkos::WriteOnly,N,xbuf_dev);
    std::copy(x, x+N, xbuf_host.getRawPtr());
    Teuchos::ArrayRCP<float> ybuf_host = node->viewBufferNonConst<float>(Kokkos::WriteOnly,N,ybuf_dev);
    std::copy(y, y+N, ybuf_host.getRawPtr());
  }

  AddOp op;
  op.x = xbuf_dev.getRawPtr();
  op.y = ybuf_dev.getRawPtr();
  op.a = alpha;

  node->parallel_for<AddOp>(0,N,op);
}
////////////////////////////////////////////////////////////////////////////////

#define N 100000

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  float *u = new float[N];
  float *v = new float[N];
  float alpha = 2.3;
  double time4=0;

  initializeVectors(u,v);

  TICK();
  axpyGPU(u,v,alpha,N);
  TOCK(time4);

  outputStats(time4);  

  delete [] u; 
  delete [] v;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void initializeVectors(float *u, float *v)
{
  for(int i=0;i<N;i++)
    {
      u[i] = 1;
      v[i] = 2*i;
    }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void outputStats(double time4)
{
  std::cout << "Time (in s) for axpy computation:" << std::endl << std::endl;

  std::cout << "GPU time: " << time4 << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
