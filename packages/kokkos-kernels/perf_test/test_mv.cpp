/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>

#include <Kokkos_Core.hpp>
#include "Kokkos_MV.hpp"
#ifndef DEVICE
#define DEVICE 1
#endif
#if DEVICE==1
typedef Kokkos::Threads execution_space;
#else
typedef Kokkos::Cuda execution_space;
#endif

typedef double FLOAT;

#define EPSILON 1e-10

typedef MultiVectorDynamic<FLOAT,execution_space>::type mv_type;
typedef mv_type::HostMirror h_mv_type;
typedef Kokkos::View<FLOAT* ,Kokkos::LayoutLeft,execution_space >  vector_type ;
typedef Kokkos::View<FLOAT* ,Kokkos::LayoutLeft,Kokkos::Threads >  h2_vector_type ;
typedef vector_type::HostMirror h_vector_type;
typedef mv_type::size_type            size_type;

void test_mv_dot(int size, int numVecs, int loop)
{
  mv_type x("X",size,numVecs);
  mv_type y("Y",size,numVecs);
  mv_type r("R",size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = Kokkos::create_mirror_view(x);
  h_mv_type h_y = Kokkos::create_mirror_view(y);
  h_mv_type h_rh = Kokkos::create_mirror_view(r);
  h_mv_type h_rd = Kokkos::create_mirror_view(r);
  h_vector_type h_a = Kokkos::create_mirror_view(a);
  h2_vector_type h_a2("h2",numVecs);

  srand(17231);
  for(int k=0;k<numVecs;k++){
    h_a2(k) = 0;//(1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	h_a(k) = 0;
	for(int i=0; i<size;i++) {
	  h_x(i,k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	  h_y(i,k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	  h_a2(k)+= h_y(i,k)*h_x(i,k);
	}
  }

  Kokkos::deep_copy(x,h_x);
  Kokkos::deep_copy(y,h_y);
  Kokkos::deep_copy(a,h_a);
  h_vector_type h_b = h_a;

  MV_Dot(a,x,y);
  execution_space::fence();

  Kokkos::deep_copy(h_a,a);
  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k)-h_a2(k))/h_a2(k));
	if(fabs((h_a(k)-h_a2(k))/h_a2(k))>EPSILON) errors++;
  }

  timespec starttime,endtime;
  clock_gettime(CLOCK_REALTIME,&starttime);
  for(int i=0;i<loop;i++)
    MV_Dot(a,x,y);
  execution_space::fence();
  clock_gettime(CLOCK_REALTIME,&endtime);
  double time = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;


  printf("MV_Dot:       %6.2lf GB/s %8i Elements %3i Vectors %s\n",
		  2*size*numVecs*sizeof(FLOAT)*loop/time*1e-9, size, numVecs, errors==0?"PASSED":"FAILED");
}

void test_mv_add(int size, int numVecs, int loop)
{
  mv_type x("X",size,numVecs);
  mv_type y("Y",size,numVecs);
  mv_type r("R",size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = Kokkos::create_mirror_view(x);
  h_mv_type h_y = Kokkos::create_mirror_view(y);
  h_mv_type h_rh = Kokkos::create_mirror_view(r);
  h_mv_type h_rd = Kokkos::create_mirror_view(r);
  h_vector_type h_a = Kokkos::create_mirror_view(a);
  h_vector_type h_b("h_b",numVecs);

  srand(17231);
  for(int k=0;k<numVecs;k++){
    h_a(k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	for(int i=0; i<size;i++) {
	  h_x(i,k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	  h_y(i,k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	  h_rh(i,k) = h_a(k)*h_y(i,k) + h_a(k)*h_x(i,k);
	}
  }

  Kokkos::deep_copy(x,h_x);
  Kokkos::deep_copy(y,h_y);
  Kokkos::deep_copy(a,h_a);
  MV_Add(r,a,x,a,y);
  execution_space::fence();

  Kokkos::deep_copy(h_rd,r);
  for(int k=0;k<numVecs;k++){
	h_a(k) = 0;
	h_b(k) = 0;
    for(int i=0; i<size;i++) {
	  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
	  h_b(k)+= h_rh(i,k)*h_rh(i,k);
	}
  }

  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k))/h_b(k));
	if(fabs((h_a(k))/h_b(k))>EPSILON) errors++;
  }

  timespec starttime,endtime;
  clock_gettime(CLOCK_REALTIME,&starttime);
  for(int i=0;i<loop;i++)
    MV_Add(r,a,x,a,y);
  execution_space::fence();
  clock_gettime(CLOCK_REALTIME,&endtime);
  double time = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

  printf("MV_Add:       %6.2lf GB/s %8i Elements %3i Vectors %s\n",
		  3*size*numVecs*sizeof(FLOAT)*loop/time*1e-9, size, numVecs, errors==0?"PASSED":"FAILED");
}

void test_mv_mulscalar(int size, int numVecs, int loop)
{
  mv_type x("X",size,numVecs);
  mv_type r("R",size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = Kokkos::create_mirror_view(x);
  h_mv_type h_rh = Kokkos::create_mirror_view(r);
  h_mv_type h_rd = Kokkos::create_mirror_view(r);
  h_vector_type h_a = Kokkos::create_mirror_view(a);
  h_vector_type h_b("h_b",numVecs);

  srand(17231);
  for(int k=0;k<numVecs;k++){
    h_a(k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	for(int i=0; i<size;i++) {
	  h_x(i,k) = (1.0*(1.0*rand()/std::numeric_limits<unsigned int>::max())-0.5)*1;
	  h_rh(i,k) = h_a(k)*h_x(i,k);
	}
  }

  Kokkos::deep_copy(x,h_x);
  Kokkos::deep_copy(a,h_a);
  MV_MulScalar(r,a,x);
  execution_space::fence();

  Kokkos::deep_copy(h_rd,r);
  for(int k=0;k<numVecs;k++){
	h_a(k) = 0;
	h_b(k) = 0;
    for(int i=0; i<size;i++) {
	  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
	  h_b(k)+= h_rh(i,k)*h_rh(i,k);
	}
  }

  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k))/h_b(k));
	if(fabs((h_a(k))/h_b(k))>EPSILON) errors++;
  }

  timespec starttime,endtime;
  clock_gettime(CLOCK_REALTIME,&starttime);
  for(int i=0;i<loop;i++)
	  MV_MulScalar(r,a,x);
  execution_space::fence();
  clock_gettime(CLOCK_REALTIME,&endtime);
  double time = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

  printf("MV_MulScalar: %6.2lf GB/s %8i Elements %3i Vectors %s\n",
		  2*size*numVecs*sizeof(FLOAT)*loop/time*1e-9, size, numVecs, errors==0?"PASSED":"FAILED");
}


int main(int argc, char **argv)
{
 int size = 200000;
 int numVecs = 17;
 int loop = 100;
 int threads_per_numa=1;
 int device = 0;
 int numa=1;

 for(int i=0;i<argc;i++)
 {
  if((strcmp(argv[i],"-n")==0)) {loop=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-d")==0)) {device=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-s")==0)) {size=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-v")==0)) {numVecs=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--threads")==0)) {threads_per_numa=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--numa")==0)) {numa=atoi(argv[++i]); continue;}
 }

 Kokkos::InitArguments args_init;

 args_init.device_id = device;

 if(numa>1 || threads>1)
 {
   args_init.num_threads = numa*threads_per_numa;
   args_init.num_numa = numa;
 }

 Kokkos::initialize( args_init );
 {
   test_mv_dot(size,numVecs,loop);
   test_mv_add(size,numVecs,loop);
   test_mv_mulscalar(size,numVecs,loop);
 }
 Kokkos::finalize();

}

