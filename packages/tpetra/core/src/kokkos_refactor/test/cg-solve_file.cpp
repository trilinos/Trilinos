/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/
#pragma message "Start Compiling"
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Array.hpp>

// I/O for Harwell-Boeing files
#include <Tpetra_MatrixIO.hpp>

#include "Tpetra_Power_Method.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_KokkosRefactor_MultiVector.hpp"
#include "Tpetra_KokkosRefactor_Vector.hpp"
#include "Tpetra_KokkosRefactor_CrsMatrix.hpp"

#include <MatrixMarket_Tpetra.hpp>
#include <algorithm>
#include <functional>
#ifndef USE_KOKKOS_CLASSIC
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#else
#include <Kokkos_TPINode.hpp>
#endif
#include <impl/Kokkos_Timer.hpp>

template<class CrsMatrix, class Vector>
void cg_solve(Teuchos::RCP<CrsMatrix> A, Teuchos::RCP<Vector> b, Teuchos::RCP<Vector> x, int myproc) {
  typedef double ScalarType;
  typedef double magnitude_type;
  typedef typename CrsMatrix::local_ordinal_type LocalOrdinalType;
  Teuchos::RCP<Vector> r,p,Ap;
  int max_iter=200;
  double tolerance = 1e-8;
  r = Tpetra::createVector<ScalarType>(A->getRangeMap());
  p = Tpetra::createVector<ScalarType>(A->getRangeMap());
  Ap = Tpetra::createVector<ScalarType>(A->getRangeMap());

  int length = r->getLocalLength();
  for(int i = 0;i<length;i++) {
    x->replaceLocalValue(i,0);
    r->replaceLocalValue(i,1);
    Ap->replaceLocalValue(i,1);
  }

  magnitude_type normr = 0;
  magnitude_type rtrans = 0;
  magnitude_type oldrtrans = 0;

  LocalOrdinalType print_freq = max_iter/10;
  if (print_freq>50) print_freq = 50;
  if (print_freq<1)  print_freq = 1;

  ScalarType one = 1.0;
  ScalarType zero = 0.0;

  double dottime = 0;
  double addtime = 0;
  double matvectime = 0;

  Kokkos::Impl::Timer timer;
  p->update(1.0,*x,0.0,*x,0.0);
  addtime += timer.seconds(); timer.reset();


  A->apply(*p, *Ap);
  matvectime += timer.seconds(); timer.reset();

  r->update(1.0,*b,-1.0,*Ap,0.0);
  addtime += timer.seconds(); timer.reset();

  rtrans = r->dot(*r);
  dottime += timer.seconds(); timer.reset();

  normr = std::sqrt(rtrans);

  if (myproc == 0) {
    std::cout << "Initial Residual = "<< normr << std::endl;
  }

  magnitude_type brkdown_tol = std::numeric_limits<magnitude_type>::epsilon();

  int num_iters=0;
  for(LocalOrdinalType k=1; k <= max_iter && normr > tolerance; ++k) {
    if (k == 1) {
      p->update(1.0,*r,0.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    else {
      oldrtrans = rtrans;
      rtrans = r->dot(*r);
      dottime += timer.seconds(); timer.reset();
      magnitude_type beta = rtrans/oldrtrans;
      p->update(beta,*p,1.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    normr = std::sqrt(rtrans);
    if (myproc == 0 && (k%print_freq==0 || k==max_iter)) {
      std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
    }

    magnitude_type alpha = 0;
    magnitude_type p_ap_dot = 0;
    A->apply(*p, *Ap);
    matvectime += timer.seconds(); timer.reset();
    p_ap_dot = Ap->dot(*p);
    dottime += timer.seconds(); timer.reset();

   if (p_ap_dot < brkdown_tol) {
      if (p_ap_dot < 0 ) {
        std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
        return;
      }
      else brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans/p_ap_dot;


    x->update(1.0,*x,alpha,*p,0.0);
    r->update(1.0,*r,-alpha,*Ap,0.0);
    addtime += timer.seconds(); timer.reset();

    num_iters = k;
  }
  rtrans = r->dot(*r);

  normr = std::sqrt(rtrans);

  //if (myproc == 0) {
    std::cout << "Final Residual = "<< normr << " after " << num_iters << " iterations" << std::endl;
    std::cout << "WAXPBY time: " << addtime << std::endl;
    std::cout << "DOT    time: " << dottime << std::endl;
    std::cout << "MATVEC time: " << matvectime << std::endl;
    std::cout << "CG     time: " << matvectime+addtime+dottime << std::endl;
 // }

}


int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  //
  // Specify types used in this example
  //
  typedef double                                                  Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType            Magnitude;
  typedef int                                                     Ordinal;

#pragma message "HUHU?"
#ifndef USE_KOKKOS_CLASSIC
#ifdef COMPILE_CUDA
  typedef Kokkos::Compat::KokkosCudaWrapperNode                Node;
#pragma message "Compile CUDA"
#else
  #ifdef KOKKOS_HAVE_OPENMP
  #pragma message "Compile OpenMP"
    typedef Kokkos::Compat::KokkosOpenMPWrapperNode                Node;
  #else
  #ifdef KOKKOS_HAVE_PTHREAD
#pragma message "Compile Threads"
    typedef Kokkos::Compat::KokkosThreadsWrapperNode                Node;
    #endif
  #endif
#endif
#else
  typedef KokkosClassic::TPINode                Node;
#endif
  typedef Tpetra::MpiPlatform<Node>                            Platform;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal,Node>          CrsMatrix;
  using Teuchos::RCP;
  using Teuchos::tuple;


  //
  // Get example parameters from command-line processor
  //
  bool printMatrix = false;
  bool verbose = false;
  int niters = 100;
  int numthreads = 1;
  int numteams = 1;
  Magnitude tolerance = 1.0e-2;
  std::string filename("dc1.mtx");
  std::string filename_vector("dc1_b.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("numthreads",&numthreads,"Number of threads per thread team.");
  cmdp.setOption("numteams",&numteams,"Number of thread teams.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("filename_vector",&filename_vector,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("tolerance",&tolerance,"Relative residual tolerance used for solver.");
  cmdp.setOption("iterations",&niters,"Maximum number of iterations.");
  cmdp.setOption("printMatrix","noPrintMatrix",&printMatrix,"Print the full matrix after reading it.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  int verboseint = verbose?1:0;
  Teuchos::ParameterList params;
  params.set("Num Threads",numthreads,"Number of Threads per Threadteam");
  params.set("Num Teams",numteams,"Number of Threadteams");
  params.set("Verbose",verboseint,"Verbose output");

  //
  // Get the communicator and node
  //
  Node anode(params);
  RCP<Node>  node(&anode,false);

  Platform platform(node);
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();

  /*Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();*/
  const int myRank = comm->getRank();
  verbose = verbose || (myRank==0);
  //
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl << std::endl;
  }
  std::cout << "Comm info: " << *comm;


  // Read Tpetra::CrsMatrix from file
  //
  RCP<CrsMatrix> A = Tpetra::MatrixMarket::Reader<CrsMatrix>::readSparseFile(filename,comm,node);

  if (printMatrix) {
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> Vector;

  if ( A->getRangeMap() != A->getDomainMap() ) {
    throw std::runtime_error("TpetraExamples::powerMethod(): operator must have domain and range maps that are equivalent.");
  }
  // create three vectors, fill z with random numbers
  Teuchos::RCP<Vector> b, x;
  std::cout << "Create Vector\n";
  RCP<const CrsMatrix::map_type> map = A->getRangeMap();
  b = Tpetra::MatrixMarket::Reader<CrsMatrix>::readVectorFile(filename_vector,comm,node,map);
  x = Tpetra::createVector<Scalar>(A->getRangeMap());

  cg_solve(A,b,x,myRank);

  return 0;
}
