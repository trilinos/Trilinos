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

#ifndef CG_SOLVE_FILE_HPP
#define CG_SOLVE_FILE_HPP

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

#include <impl/Kokkos_Timer.hpp>

#include <algorithm>
#include <functional>


struct result_struct {
  double addtime,dottime,matvectime,final_residual;
  int niters;
  result_struct(double add, double dot, double matvec,int niter,double res):
    addtime(add),dottime(dot),matvectime(matvec),
    final_residual(res),niters(niter) {};
};


template<class Node>
Teuchos::XMLTestNode machine_configuration(Node node);

template<class Node>
Teuchos::XMLTestNode machine_configuration(Node node) {
  Teuchos::XMLTestNode config = Teuchos::PerfTest_MachineConfig();
  config.addString("KokkosNodeType",node->name());
  return config;
}


template<class CrsMatrix>
Teuchos::XMLTestNode test_entry(
    const std::string& filename_matrix,
    const std::string& filename_vector,
    int nsize, int mpi_ranks, int teams, int threads,
    Teuchos::RCP<CrsMatrix> A,
    result_struct results,int max_iters,int res_tolerance,
    double tol_small, double tol_large) {
  Teuchos::XMLTestNode configuration("TestConfiguration");
  configuration.addInt("MPI_Ranks",mpi_ranks);
  configuration.addInt("Teams",teams);
  configuration.addInt("Threads",threads);
  if(filename_matrix.empty()) {
    std::ostringstream strs;
    strs << "MiniFE-Generated " << nsize;
    configuration.addString("Matrix_File",strs.str());
  } else
    configuration.addString("Matrix_File",filename_matrix);

  if (filename_vector.empty ()) {
    std::ostringstream strs;
    strs << "MiniFE-Generated " << nsize;
    configuration.addString ("Vector_File", strs.str ());
  } else {
    configuration.addString ("Vector_File", filename_vector);
  }

  configuration.addInt("Matrix_Rows",A->getGlobalNumRows());
  configuration.addInt("Matrix_Cols",A->getGlobalNumCols());
  configuration.addInt("Matrix_NNZ",A->getGlobalNumEntries());
  configuration.addInt("Max_Iterations",max_iters);
  configuration.addInt("Tolerance",res_tolerance);


  Teuchos::XMLTestNode times("TestResults");
  times.addValueTolerance("Time_Add", Teuchos::ValueTolerance(results.addtime,tol_large));
  times.addValueTolerance("Time_Dot", Teuchos::ValueTolerance(results.dottime,tol_large));
  times.addValueTolerance("Time_MatVec", Teuchos::ValueTolerance(results.matvectime,tol_large));
  times.addValueTolerance("Time_CGSolve", Teuchos::ValueTolerance(results.matvectime+results.addtime+results.dottime,tol_small));
  times.addValueTolerance("Result_Iterations",Teuchos::ValueTolerance(results.niters,
                                                                      results.niters>0?results.niters-1:0,
                                                                      results.niters+1));
  times.addValueTolerance("Result_Final_Residual",Teuchos::ValueTolerance(results.final_residual,tol_small));

  //times.addString("Result_Residual", ValueTolerance(atof(argc[8]),4,6).as_string());

  Teuchos::XMLTestNode test("TPetra::CG-Solve");
  Teuchos::XMLTestNode entry("TestEntry");
  entry.addChild(configuration);
  entry.addChild(times);
  test.addChild(entry);
  return test;
}


template<class CrsMatrix, class Vector>
result_struct
cg_solve (Teuchos::RCP<CrsMatrix> A, Teuchos::RCP<Vector> b, Teuchos::RCP<Vector> x, int myproc)
{
  static_assert (std::is_same<typename CrsMatrix::scalar_type, typename Vector::scalar_type>::value,
                 "The CrsMatrix and Vector template parameters must have the same scalar_type.");

  typedef typename Vector::scalar_type ScalarType;
  typedef typename Vector::mag_type magnitude_type;
  typedef typename Vector::local_ordinal_type LO;
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

  LO print_freq = max_iter/10;
  if (print_freq>50) print_freq = 50;
  if (print_freq<1)  print_freq = 1;

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
  LO k;
  for(k=1; k <= max_iter && normr > tolerance; ++k) {
    if (k == 1) {
      p->update(1.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    else {
      oldrtrans = rtrans;
      rtrans = r->dot(*r);
      dottime += timer.seconds(); timer.reset();
      magnitude_type beta = rtrans/oldrtrans;
      p->update(1.0,*r,beta);
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
        return result_struct(0,0,0,0,0);
      }
      else brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans/p_ap_dot;


    x->update(alpha,*p,1.0);
    r->update(-alpha,*Ap,1.0);
    addtime += timer.seconds(); timer.reset();

  }
  rtrans = r->dot(*r);

  normr = std::sqrt(rtrans);


  return result_struct(addtime,dottime,matvectime,k-1,normr);
}

template<class Node>
int
run (int argc, char *argv[])
{
  using Teuchos::RCP;
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
  // Get example parameters from command-line processor
  //
  bool printMatrix = false;
  bool verbose = false;
  int niters = 100;
  int numthreads = 1;
  int numteams = 1;
  int numgpus = 1;
  int skipgpu = 999;
  int nsize = 20;
  mag_type tolerance = 1.0e-2;
  std::string filename;
  std::string filename_vector;
  std::string testarchive("Tpetra_PerformanceTests.xml");
  std::string hostname;

  double tol_small = 0.05;
  double tol_large = 0.10;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("numthreads",&numthreads,"Number of threads per thread team.");
  cmdp.setOption("numteams",&numteams,"Number of thread teams.");
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
  cmdp.setOption("tol_large",&tol_small,"Tolerance for individual times.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  int myRank = 0;
#ifdef HAVE_MPI
  (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
#endif // HAVE_MPI

  Kokkos::InitArguments kokkosArgs;
  kokkosArgs.num_threads = numthreads;
  kokkosArgs.num_numa = numteams; // ???
  kokkosArgs.device_id = myRank & numgpus;
  kokkosArgs.skip_device = skipgpu;
  kokkosArgs.disable_warnings = ! verbose;

  Kokkos::initialize (kokkosArgs);

  //
  // Get the communicator and node
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

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
    A = Tpetra::Utils::MatrixGenerator<crs_matrix_type>::generate_miniFE_matrix (nsize, comm, Teuchos::null);
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
                                     map->getNode (), map);
  } else {
    typedef Tpetra::Utils::MatrixGenerator<crs_matrix_type> gen_type;
    b = gen_type::generate_miniFE_vector (nsize, map->getComm (),
                                          map->getNode ());
  }

  // The vector x on input is the initial guess for the CG solve.
  // On output, it is the approximate solution.
  RCP<vec_type> x (new vec_type (A->getDomainMap ()));

  // Solve the linear system Ax=b using CG.
  result_struct results = cg_solve (A, b, x, myRank);

  // Print results.
  if (myRank == 0) {
    Teuchos::XMLTestNode machine_config =
      machine_configuration (map->getNode ());
    Teuchos::XMLTestNode test =
      test_entry (filename, filename_vector, nsize, comm->getSize (), numteams,
                  numthreads, A, results, niters, tolerance, tol_small,
                  tol_large);
    Teuchos::PerfTestResult comparison_result =
      Teuchos::PerfTest_CheckOrAdd_Test (machine_config, test, testarchive, hostname);
    switch (comparison_result) {
      case Teuchos::PerfTestPassed:
        cout << "PASSED" << endl;
        break;
      case Teuchos::PerfTestFailed:
        cout << "FAILED" << endl;
        break;
      case Teuchos::PerfTestNewMachine:
        cout << "PASSED. Adding new machine entry." << endl;
        break;
      case Teuchos::PerfTestNewConfiguration:
        cout << "PASSED. Adding new machine configuration." << endl;
        break;
      case Teuchos::PerfTestNewTest:
        cout << "PASSED. Adding new test entry." << endl;
        break;
      case Teuchos::PerfTestNewTestConfiguration:
        cout << "PASSED. Adding new test entry configuration." << endl;
        break;
      case Teuchos::PerfTestUpdatedTest:
        cout << "PASSED. Updating test entry." << endl;
        break;
    default:
      cout << "FAILED: Invalid comparison result." << endl;
    }
    if (verbose) {
      cout << test << endl;
    }
  }

  return EXIT_SUCCESS;
}

#endif
