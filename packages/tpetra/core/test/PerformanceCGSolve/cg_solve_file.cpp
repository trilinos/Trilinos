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

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Version.hpp"

#include "TpetraUtils_MatrixGenerator.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLPerfTestArchive.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

#include <algorithm>
#include <functional>

template<class CrsMatrix, class Vector>
bool cg_solve (Teuchos::RCP<CrsMatrix> A, Teuchos::RCP<Vector> b, Teuchos::RCP<Vector> x, int myproc, double tolerance, int max_iter)
{
  using Teuchos::TimeMonitor;
  std::string addTimerName = "CG: axpby";
  std::string matvecTimerName = "CG: spmv";
  std::string dotTimerName = "CG: dot";
  static_assert (std::is_same<typename CrsMatrix::scalar_type, typename Vector::scalar_type>::value,
                 "The CrsMatrix and Vector template parameters must have the same scalar_type.");

  typedef typename Vector::scalar_type ScalarType;
  typedef typename Vector::mag_type magnitude_type;
  typedef typename Vector::local_ordinal_type LO;
  Teuchos::RCP<Vector> r,p,Ap;
  //int max_iter=200;
  //double tolerance = 1e-8;
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

  LO print_freq = max_iter/10;
  print_freq = std::min(print_freq, 50);
  print_freq = std::max(print_freq, 1);
  {
    TimeMonitor t(*TimeMonitor::getNewTimer(addTimerName));
    p->update(1.0,*x,0.0,*x,0.0);
  }
  {
    TimeMonitor t(*TimeMonitor::getNewTimer(matvecTimerName));
    A->apply(*p, *Ap);
  }
  {
    TimeMonitor t(*TimeMonitor::getNewTimer(addTimerName));
    r->update(1.0,*b,-1.0,*Ap,0.0);
  }
  {
    TimeMonitor t(*TimeMonitor::getNewTimer(dotTimerName));
    rtrans = r->dot(*r);
  }

  normr = std::sqrt(rtrans);

  if (myproc == 0) {
    std::cout << "Initial Residual = "<< normr << std::endl;
  }

  magnitude_type brkdown_tol = std::numeric_limits<magnitude_type>::epsilon();
  LO k;
  for(k=1; k <= max_iter && normr > tolerance; ++k) {
    if (k == 1) {
        TimeMonitor t(*TimeMonitor::getNewTimer(addTimerName));
        p->update(1.0,*r,0.0);
    }
    else {
      oldrtrans = rtrans;
      {
        TimeMonitor t(*TimeMonitor::getNewTimer(dotTimerName));
        rtrans = r->dot(*r);
      }
      {
        TimeMonitor t(*TimeMonitor::getNewTimer(addTimerName));
        magnitude_type beta = rtrans/oldrtrans;
        p->update(1.0,*r,beta);
      }
    }
    normr = std::sqrt(rtrans);
    if (myproc == 0 && (k%print_freq==0 || k==max_iter)) {
      std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
    }

    magnitude_type alpha = 0;
    magnitude_type p_ap_dot = 0;
    {
      TimeMonitor t(*TimeMonitor::getNewTimer(matvecTimerName));
      A->apply(*p, *Ap);
    }
    {
      TimeMonitor t(*TimeMonitor::getNewTimer(dotTimerName));
      p_ap_dot = Ap->dot(*p);
    }
    {
      TimeMonitor t(*TimeMonitor::getNewTimer(addTimerName));
      if (p_ap_dot < brkdown_tol) {
        if (p_ap_dot < 0 ) {
          std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
          return false;
        }
        else brkdown_tol = 0.1 * p_ap_dot;
      }
      alpha = rtrans / p_ap_dot;
      x->update(alpha,*p,1.0);
      r->update(-alpha,*Ap,1.0);
    }
  }
  {
    TimeMonitor t(*TimeMonitor::getNewTimer(dotTimerName));
    rtrans = r->dot(*r);
  }

  normr = std::sqrt(rtrans);
  return true;
}

namespace CGParams
{
  int nsize = 20;
  bool printMatrix = false;
  bool verbose = false;
  int niters = 100;
  double tolerance = 1.0e-2;
  std::string filename;
  std::string filename_vector;
  std::string testarchive("Tpetra_PerformanceTests.xml");
  std::string hostname;

  double tol_small = 0.05;
  double tol_large = 0.10;
}

template<class Node>
int run()
{
  using namespace CGParams;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::StackedTimer;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  typedef Tpetra::Vector<>::scalar_type                 Scalar;
  typedef typename Tpetra::Map<>::local_ordinal_type    LO;
  typedef typename Tpetra::Map<>::global_ordinal_type   GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>          crs_matrix_type;
  typedef Tpetra::Vector<Scalar,LO,GO,Node>             vec_type;
  typedef Tpetra::Map<LO,GO,Node>                       map_type;

  typedef typename vec_type::mag_type                   mag_type;

  //
  // Get the communicator and node
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  int myRank = comm->getRank();

  //
  // Say hello, print some communicator info
  //
  if (verbose) {
    if (myRank == 0) {
      cout << "Comm info: ";
    }
    cout << *comm;
  }


  // Read Tpetra::CrsMatrix from file
  //
  RCP<crs_matrix_type> A;
  if (! filename.empty ()) {
    A = Tpetra::MatrixMarket::Reader<crs_matrix_type>::readSparseFile (filename, comm);
  }
  else {
    A = Tpetra::Utils::MatrixGenerator<crs_matrix_type>::generate_miniFE_matrix (nsize, comm);
  }

  if (printMatrix) {
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream (Teuchos::rcpFromRef (cout));
    A->describe (*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    cout << endl << A->description() << endl << endl;
  }

  // This is a collective over the Maps' communicator.
  if (! A->getRangeMap ()->isSameAs (* (A->getDomainMap ()))) {
    throw std::runtime_error ("The matrix must have domain and range maps that are the same.");
  }

  // Either read the right-hand side b of the linear system from a
  // file, or generate it.
  RCP<const map_type> map = A->getRangeMap ();
  RCP<vec_type> b;
  if (nsize < 0) {
    typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
    b = reader_type::readVectorFile (filename_vector, map->getComm (),
                                     map);
  } else {
    typedef Tpetra::Utils::MatrixGenerator<crs_matrix_type> gen_type;
    b = gen_type::generate_miniFE_vector (nsize, map->getComm ()
                                         );
  }

  // The vector x on input is the initial guess for the CG solve.
  // On output, it is the approximate solution.
  RCP<vec_type> x (new vec_type (A->getDomainMap ()));

  // Solve the linear system Ax=b using CG.
  RCP<StackedTimer> timer = rcp(new StackedTimer("CG: global"));
  TimeMonitor::setStackedTimer(timer);

  bool success = cg_solve(A, b, x, myRank, tolerance, niters);
  timer->stopBaseTimer();
  StackedTimer::OutputOptions options;
  options.print_warnings = false;
  timer->report(std::cout, comm, options);

  std::string testBaseName = "Tpetra CGSolve ";
  if (Tpetra::Details::Behavior::cudaLaunchBlocking()) testBaseName += "CUDA_LAUNCH_BLOCKING ";

  auto xmlOut = timer->reportWatchrXML(testBaseName + std::to_string(comm->getSize()) + " ranks", comm);
  if(myRank == 0)
  {
    if(xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
    if(success)
      std::cout << "End Result: TEST PASSED\n";
    else
      std::cout << "End Result: TEST FAILED\n";
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
  using namespace Kokkos::Compat;
  using namespace CGParams;
  using default_exec = Tpetra::Details::DefaultTypes::execution_space;
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  //
  // Get example parameters from command-line processor
  //
  int numthreads = 1;
  int numgpus = 1;
  //Default for #gpus is $KOKKOS_NUM_DEVICES
  const char* rawKokkosNumDevices = std::getenv("KOKKOS_NUM_DEVICES");
  if(rawKokkosNumDevices)
    numgpus = std::atoi(rawKokkosNumDevices);
  int skipgpu = 999;

  bool useCuda = false;
  bool useOpenMP = false;
  bool useThreads = false;
  bool useSerial = false;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("numthreads",&numthreads,"Number of threads per thread team.");
  cmdp.setOption("numgpus",&numgpus,"Number of GPUs.");
  cmdp.setOption("skipgpu",&skipgpu,"Do not use this GPU.");
  cmdp.setOption("hostname",&hostname,"Override of hostname for PerfTest entry.");
  cmdp.setOption("testarchive",&testarchive,"Set filename for Performance Test archive.");
  cmdp.setOption("filename",&filename,"Filename for test matrix.");
  cmdp.setOption("filename_vector",&filename_vector,"Filename for test matrix vector.");
  cmdp.setOption("tolerance",&tolerance,"Relative residual tolerance used for solver.");
  cmdp.setOption("iterations",&niters,"Maximum number of iterations.");
  cmdp.setOption("printMatrix","noPrintMatrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("size",&nsize,"Generate miniFE matrix with X^3 elements.");
  cmdp.setOption("tol_small",&tol_small,"Tolerance for total CG-Time and final residual.");
  cmdp.setOption("tol_large",&tol_large,"Tolerance for individual times.");
  //Only provide the option to use Node types that are actually enabled in this build
#ifdef HAVE_TPETRA_INST_CUDA
  cmdp.setOption("cuda","no-cuda",&useCuda,"Use Cuda node");
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
  cmdp.setOption("openmp","no-openmp",&useOpenMP,"Use OpenMP node");
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
  cmdp.setOption("threads","no-threads",&useThreads,"Use Pthreads node");
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
  cmdp.setOption("serial","no-serial",&useSerial,"Use Serial node");
#endif
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  //If no node type was explicitly requested, use Tpetra's default node
  if(!useCuda && !useOpenMP && !useThreads && !useSerial)
  {
#ifdef HAVE_TPETRA_INST_CUDA
    if(std::is_same<default_exec, Kokkos::Cuda>::value)
    {
      std::cout << "No node specified in command-line args, so using default (Cuda)\n";
      useCuda = true;
    }
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
    if(std::is_same<default_exec, Kokkos::OpenMP>::value)
    {
      std::cout << "No node specified in command-line args, so using default (OpenMP)\n";
      useOpenMP = true;
    }
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
    if(std::is_same<default_exec, Kokkos::Threads>::value)
    {
      std::cout << "No node specified in command-line args, so using default (Pthreads)\n";
      useThreads = true;
    }
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
    if(std::is_same<default_exec, Kokkos::Serial>::value)
    {
      std::cout << "No node specified in command-line args, so using default (Serial)\n";
      useSerial = true;
    }
#endif
  }

  int myRank = 0;
#ifdef HAVE_MPI
  (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
#endif // HAVE_MPI

  Kokkos::InitArguments kokkosArgs;
  kokkosArgs.num_threads = numthreads;
  kokkosArgs.device_id = myRank % numgpus;
  kokkosArgs.skip_device = skipgpu;
  kokkosArgs.disable_warnings = !verbose;

  Kokkos::initialize (kokkosArgs);
  {
#ifdef HAVE_TPETRA_INST_CUDA
    if(useCuda)
      return run<KokkosCudaWrapperNode>();
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
    if(useOpenMP)
      return run<KokkosOpenMPWrapperNode>();
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
    if(useThreads)
      return run<KokkosThreadsWrapperNode>();
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
    if(useSerial)
      return run<KokkosSerialWrapperNode>();
#endif
    if(myRank == 0)
      std::cerr << "Error: no node type was enabled. CG was not run.\n";
  }
  Kokkos::finalize ();
  return 0;
}
