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

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <iostream>
#include <functional>
#include <algorithm>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_Version.hpp"

#include "Kokkos_ThrustGPUNode.hpp"
#include "Kokkos_TPINode.hpp"

#include "PDP10_TestOps.hpp"

  struct CompStats {
    double seconds;
    double flops;
  };

  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;
  using std::scientific;
  using std::setprecision;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::null;
  using Teuchos::parameterList;

  int numPthreads;

  template <class NODE>
  RCP<NODE> getNode() {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  using Kokkos::SerialNode;
  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    static RCP<SerialNode> serialnode;
    if (serialnode == null) {
      Teuchos::ParameterList pl; 
      serialnode = rcp(new SerialNode(pl));
    }
    return serialnode;
  }

  using Kokkos::TPINode;
  template <>
  RCP<TPINode> getNode<TPINode>() {
    static RCP<TPINode> tpinode;
    if (tpinode == null) {
      Teuchos::ParameterList pl; 
      pl.set("Num Threads",numPthreads);
      pl.set("Verbose",1);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }

  using Kokkos::ThrustGPUNode;
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    static RCP<ThrustGPUNode> thrustnode;
    if (thrustnode == null) {
      Teuchos::ParameterList pl; 
      pl.set("Device Number",0);
      pl.set("Verbose",1);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }


  /////////////////////////////////////////////////////////
  template <class SCALAR, class NODE>
  void sumTest(RCP<Teuchos::Time> time, int N)
  {
    Teuchos::ArrayRCP<SCALAR> x;
    RCP<NODE> node = getNode<NODE>();
    Kokkos::ReadyBufferHelper<NODE> rbh(node);
    SCALAR result;
    {
      x = node->template allocBuffer<SCALAR>(N);
    }
    // set x[i] = 1, i=0:N-1
    {
      InitOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addNonConstBuffer<SCALAR>(x);
      rbh.end();
      node->parallel_for(0,N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      Teuchos::TimeMonitor localTimer(*time);
      SumOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addConstBuffer<SCALAR>(x);
      rbh.end();
      for (int i=0; i<3; ++i) {
        result = node->parallel_reduce(0,N,wdp);
      }
    }
    SCALAR expectedResult = (SCALAR)(N);
    TEUCHOS_TEST_FOR_EXCEPT(result != expectedResult);
    // compute sum x[i], i=1:N-2
    {
      Teuchos::TimeMonitor localTimer(*time);
      SumOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addConstBuffer<SCALAR>(x);
      rbh.end();
      for (int i=0; i<3; ++i) {
        result = node->parallel_reduce(1,N-1,wdp);
      }
    }
    expectedResult = (SCALAR)(N-2);
    TEUCHOS_TEST_FOR_EXCEPT(result != expectedResult);
    {
      x = Teuchos::null;
    }
    cout << "SumOp<" << Teuchos::TypeNameTraits<SCALAR>::name() << "> time: " << time->totalElapsedTime()/3.0 << std::endl;
  }


  /////////////////////////////////////////////////////////
  template <class Node>
  RCP<typename Kokkos::DefaultKernels<float,int,Node>::SparseOps::template bind_scalar<float>::other_type>
  gen_prob(RCP<Node> node, int N, size_t &totalNNZ)
  {
    typedef typename Kokkos::DefaultKernels<float,int,Node>::SparseOps   DSM;
    typedef typename DSM::template bind_scalar<float>::other_type       fDSM;
    typedef typename fDSM::template graph<int,Node>::graph_type                  GRPH;
    typedef typename fDSM::template matrix<float,int,Node>::matrix_type           MAT;
    // generate symmetric tridiagonal matrix
    RCP<GRPH> G = rcp(new GRPH(N,N,node,null));
    RCP<MAT>  A= rcp(new MAT(G,null));
    // allocate buffers for offsets, indices and values
    totalNNZ = 3*N - 2;
    ArrayRCP<size_t> offsets(N+1);
    ArrayRCP<int>    inds(totalNNZ);
    ArrayRCP<float>  vals(totalNNZ);
    {
      size_t NNZsofar = 0;
      offsets[0] = NNZsofar;
      inds[NNZsofar] = 0; inds[NNZsofar+1] =  1;
      vals[NNZsofar] = 2; vals[NNZsofar+1] = -1;
      NNZsofar += 2;
      for (int i=1; i != N-1; ++i) {
        offsets[i] = NNZsofar;
        inds[NNZsofar] = i-1; inds[NNZsofar+1] = i; inds[NNZsofar+2] = i+1;
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2; vals[NNZsofar+2] =  -1;
        NNZsofar += 3;
      }
      offsets[N-1] = NNZsofar;
      inds[NNZsofar] = N-2; inds[NNZsofar+1] = N-1;
      vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2;
      NNZsofar += 2;
      offsets[N]   = NNZsofar;
      TEUCHOS_TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
    }
    G->setStructure(offsets, inds);
    offsets = Teuchos::null;
    inds    = Teuchos::null;
    A->setValues(vals);
    vals    = Teuchos::null;
    fDSM::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,parameterList());
    RCP<fDSM> dsm = rcp(new fDSM(node));
    dsm->setGraphAndMatrix(G,A);
    return dsm;
  }


  /////////////////////////////////////////////////////////
  template <class Node>
  CompStats power_method(RCP<Teuchos::Time> time, int N, size_t niters, float tolerance, bool verbose)  
  {
    typedef Kokkos::MultiVector<float,Node>                        MV;
    typedef Kokkos::DefaultArithmetic< Kokkos::MultiVector<float,Node> > DMVA;
    RCP<Node> node = getNode<Node>();
    // create matrix
    size_t NNZ = 0;
    RCP<typename Kokkos::DefaultKernels<float,int,Node>::SparseOps::template bind_scalar<float>::other_type> A;
    A = gen_prob<Node>(node,N,NNZ);
    // Variables needed for iteration
    const float ONE  = 1.0f;
    const float ZERO = 0.0f;
  
    // create vectors
    Kokkos::MultiVector<float,Node> z(node), q(node), r(node);
    {
      ArrayRCP<float> zvals = node->template allocBuffer<float>(N),
                      qvals = node->template allocBuffer<float>(N),
                      rvals = node->template allocBuffer<float>(N);
      // allocate vectors, set pointers
      z.initializeValues(N,1,zvals,N);
      q.initializeValues(N,1,qvals,N);
      r.initializeValues(N,1,rvals,N);
    }
    Teuchos::ScalarTraits<float>::seedrandom(1234);
    DMVA::Random(z);
    DMVA::Init(q,ZERO);
    DMVA::Init(r,ZERO);
    float lambda = 0.0f;
    float normz, residual = 0.0f;
    // power iteration
    {
      Teuchos::TimeMonitor lcltimer(*time);
      for (size_t iter = 0; iter < niters; ++iter) {
        // cout << "z: "; printVec(z);
        normz = Teuchos::ScalarTraits<float>::squareroot( DMVA::Norm2Squared(z) );  // Compute |z|
        // cout << "|z|: " << normz << endl;
        DMVA::Scale(q, ONE/normz, (const MV&)z);                                    // Set q = z / |z|
        // cout << "q: "; printVec(q);
        A->template multiply<float,float>(Teuchos::NO_TRANS, ONE, (const MV&)q, z); // Compute z = A*q
      }
    }
    lambda = DMVA::Dot(q,z);                                                       // Approximate maximum eigenvalue: lamba = dot(q,z)
    DMVA::GESUM(r, ONE, (const MV&)z, -lambda, (const MV&)q, ZERO);
    residual = Teuchos::ScalarTraits<float>::squareroot(DMVA::Norm2Squared(r)) / lambda;
    cout << "PM final residual = " << residual << ", lambda = " << lambda << endl;
    CompStats cs;
    cs.seconds = time->totalElapsedTime();
    cs.flops   = niters * ( 2.0*N    // |z|
                          + 1.0*N    // q = z/|z|
                          + 2.0*NNZ  // A*q
                          );
    return cs;
  }


  /////////////////////////////////////////////////////////
  template <class Node>
  CompStats conjugate_gradient(RCP<Teuchos::Time> time, int N, size_t niters, float tolerance, bool verbose) 
  {
    typedef Kokkos::MultiVector<float,Node> MV;
    typedef Kokkos::DefaultArithmetic<MV> DMVA;
    RCP<Node> node = getNode<Node>();
    // create matrix
    size_t NNZ = 0;
    RCP<typename Kokkos::DefaultKernels<float,int,Node>::SparseOps::template bind_scalar<float>::other_type> A;
    A = gen_prob<Node>(node,N,NNZ);
    // Variables needed for iteration
    const float ONE  = 1.0f;
    const float ZERO = 0.0f;
    float residual = ZERO;
    float alpha, beta, r2, pAp, r2_old, norm0;
    // create vectors
    ArrayRCP<float> xvals = node->template allocBuffer<float>(N),
                    rvals = node->template allocBuffer<float>(N),
                    pvals = node->template allocBuffer<float>(N),
                    zvals = node->template allocBuffer<float>(N);
    // allocate vectors, set pointers
    MV x(node), r(node), p(node), z(node);
    x.initializeValues(N,1,xvals,N);
    r.initializeValues(N,1,rvals,N);
    p.initializeValues(N,1,pvals,N);
    z.initializeValues(N,1,zvals,N);
    // initial x = ones()
    // right hand size = random()
    // r = A*x - b = 0 - b = random()
    DMVA::Init(x,ZERO);
    Teuchos::ScalarTraits<float>::seedrandom(1234);
    DMVA::Random(r);
    r2 = DMVA::Norm2Squared(r);
    norm0 = Teuchos::ScalarTraits<float>::squareroot(r2);
    DMVA::Assign(p,(const MV&)r);
    // conjugate gradient iteration
    {
      Teuchos::TimeMonitor lcltimer(*time);
      for (size_t iter = 0; iter < niters; ++iter) {
        A->template multiply<float,float>(Teuchos::NO_TRANS, ONE, (const MV&)p, z);
        pAp = DMVA::Dot(p,z);                     // pAp = <p,z> = <p,A*p>
        alpha = r2 / pAp;                         // alpha = <r,r> / <p,A*p>
        DMVA::GESUM(x, alpha,(const MV&)p,ONE);   // x = x + alpha*p
        DMVA::GESUM(r,-alpha,(const MV&)p,ONE);   // r = r - alpha*p
        r2_old = r2;
        r2 = DMVA::Norm2Squared(r);               // <new r, new r>
        beta = r2 / r2_old;                       // beta = <old r, old r> / <new r, new r>
        DMVA::GESUM(p,ONE,(const MV&)r,beta);     // p = beta*p + r
      }
    }
    residual = Teuchos::ScalarTraits<float>::squareroot(r2) / norm0;
    cout << "CG final residual = " << residual << endl;
    CompStats cs;
    cs.seconds = time->totalElapsedTime();
    cs.flops   = niters * ( 2.0*NNZ // A*p
                          + 2.0*N   // <p,z>
                          + 3.0*N   // x update
                          + 3.0*N   // r update
                          + 2.0*N   // <r,r>
                          + 3.0*N   // p update
                          );
    return cs;
  }


  /////////////////////////////////////////////////////////
  // create a timer, for timing time.
  /////////////////////////////////////////////////////////
  inline RCP<Teuchos::Time> getNewTimer(const std::string &lbl) 
  {
    return Teuchos::TimeMonitor::getNewTimer(lbl);
  }


  /////////////////////////////////////////////////////////
  // for debugging
  /////////////////////////////////////////////////////////
  template <class Node>
  void printVec(const Kokkos::MultiVector<float,Node> &vec) {
    const int n = vec.getNumRows();
    ArrayRCP<const float> vals = vec.getValues(0);
    ArrayRCP<const float> vals_h = vec.getNode()->viewBuffer(n,vals);
    for (int i=0; i<n; ++i) {
      cout << "   " << vals_h[i];
    }
    cout << endl;
  }


  /////////////////////////////////////////////////////////
  // Do it!
  /////////////////////////////////////////////////////////
  int main(int argc, char **argv) {
    Teuchos::CommandLineProcessor cmdp(false,true);
    int numIters = 25;
    int N = 100000;
    numPthreads = 8;
    cmdp.setOption("N",  &N,  "Size of the grid.");
    cmdp.setOption("numPthreads",&numPthreads,"Number of Pthreads for use by TPI.");
    cmdp.setOption("numIters",&numIters,"Number of CG/PM iterations.h");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    cout << endl << Kokkos::Kokkos_Version() << endl;

    // 
    cout << "\nTesting SerialNode" << endl;
    sumTest<int,  SerialNode>(getNewTimer("sum int srl"  ),N);
    sumTest<float,SerialNode>(getNewTimer("sum float srl"),N);
    CompStats pm_serial = power_method      <SerialNode   >(getNewTimer("pm srl"),N,numIters,1e-5f,1);
    CompStats cg_serial = conjugate_gradient<SerialNode   >(getNewTimer("cg srl"),N,numIters,1e-5f,1);
    // 
    cout << "\nTesting TPINode" << endl;
    sumTest<int,  TPINode>(getNewTimer("sum int tpi"),N);
    sumTest<float,TPINode>(getNewTimer("sum float tpi"),N);
    CompStats pm_tpi = power_method      <TPINode      >(getNewTimer("pm tpi"),N,numIters,1e-5f,1);
    CompStats cg_tpi = conjugate_gradient<TPINode      >(getNewTimer("cg tpi"),N,numIters,1e-5f,1);
    // 
    cout << "\nTesting ThrustGPUNode" << endl;
    sumTest<int,  ThrustGPUNode>(getNewTimer("sum int gpu"),N);
    sumTest<float,ThrustGPUNode>(getNewTimer("sum float gpu"),N);
    CompStats pm_gpu = power_method      <ThrustGPUNode>(getNewTimer("pm gpu"),N,numIters,1e-5f,1);
    CompStats cg_gpu = conjugate_gradient<ThrustGPUNode>(getNewTimer("cg gpu"),N,numIters,1e-5f,1);

    // 
    // Timings, or it didn't happen...
    //
    cout.precision(3);
    cout << endl 
      << setw(6) << "test"    << "  " << setw(10) <<               "seconds"         << "  " << setw(10)<<               "flops"         << "  " << setw(13) << "mflops/second"                                                     << endl
      << setw(6) << "------"  << "  " << setw(10) <<               "----------"      << "  " << setw(10)<<               "----------"    << "  " << setw(13) << "-------------"                                                     << endl
      << setw(6) << "pm srl"  << "  " << setw(10) << scientific << pm_serial.seconds << "  " << setw(10)<< scientific << pm_serial.flops << "  " << setw(13) << fixed << pm_serial.flops/pm_serial.seconds/1.0e6 << endl
      << setw(6) << "pm tpi"  << "  " << setw(10) << scientific << pm_tpi.seconds    << "  " << setw(10)<< scientific << pm_tpi.flops    << "  " << setw(13) << fixed << pm_tpi.flops   /pm_tpi.seconds   /1.0e6 << endl
      << setw(6) << "pm gpu"  << "  " << setw(10) << scientific << pm_gpu.seconds    << "  " << setw(10)<< scientific << pm_gpu.flops    << "  " << setw(13) << fixed << pm_gpu.flops   /pm_gpu.seconds   /1.0e6 << endl
      << setw(6) << "------"  << "  " << setw(10) <<               "----------"      << "  " << setw(10)<<               "----------"    << "  " << setw(13) << "-------------"                                                     << endl
      << setw(6) << "cg srl"  << "  " << setw(10) << scientific << cg_serial.seconds << "  " << setw(10)<< scientific << cg_serial.flops << "  " << setw(13) << fixed << cg_serial.flops/cg_serial.seconds/1.0e6 << endl
      << setw(6) << "cg tpi"  << "  " << setw(10) << scientific << cg_tpi.seconds    << "  " << setw(10)<< scientific << cg_tpi.flops    << "  " << setw(13) << fixed << cg_tpi.flops   /cg_tpi.seconds   /1.0e6 << endl
      << setw(6) << "cg gpu"  << "  " << setw(10) << scientific << cg_gpu.seconds    << "  " << setw(10)<< scientific << cg_gpu.flops    << "  " << setw(13) << fixed << cg_gpu.flops   /cg_gpu.seconds   /1.0e6 << endl
      ;
    // cout << endl;
    // Teuchos::TimeMonitor::summarize();
    cout << "\nEnd Result: TEST PASSED" << endl;
    return 0;
  }
