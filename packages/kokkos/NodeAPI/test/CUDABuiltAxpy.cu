#include <Kokkos_ThrustGPUNode.hpp>
#include "Kokkos_ThrustGPUNode.cuh"

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

struct AddOp 
{
  float a;
  const float *x;
  float *y;

  inline KERNEL_PREFIX void execute(int i) const
  {
    y[i] += a*x[i];
  }
};

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
