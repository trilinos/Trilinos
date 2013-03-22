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

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Ifpack2_Parameters.hpp>

#include <Tpetra_HybridPlatform.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include "build_eigproblem.hpp"
#include "build_eigsolver.hpp"

using std::string;
using std::cout;
using std::endl;
using Teuchos::FancyOStream;
using Teuchos::ParameterList;
using Teuchos::ParameterXMLFileReader;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::fancyOStream;
using Teuchos::Comm;
using Teuchos::Time;
using Teuchos::Array;
using Teuchos::ScalarTraits;

namespace AnasaziTpetraHybridDriverTest {
  ParameterList plTestParams;
};

#ifdef HAVE_TPETRA_THREADED_MKL
#include <mkl_service.h>
template <class Node>
class ThreadedBlasKiller {
  public:
  static void kill() {
  }
};
template <>
class ThreadedBlasKiller<Kokkos::ThrustGPUNode> {
  public:
  static void kill() {
    mkl_set_dynamic(false);
    mkl_set_num_threads(1);
  }
};
#else
template <class Node>
class ThreadedBlasKiller {
  public:
  static void kill() {}
};
#endif

template <class Node>
class NodeDetails {
  public:
  static void printDetails(const RCP<FancyOStream> &fos, const RCP<const Node> &node) {}
  static void clear(const RCP<Node> &node) {}
};

template <>
class NodeDetails<Kokkos::ThrustGPUNode> {
  public:
  static void printDetails(const RCP<FancyOStream> &fos, const RCP<const Kokkos::ThrustGPUNode> &node) {
    node->printStatistics(fos);
  }
  static void clear(const RCP<Kokkos::ThrustGPUNode> &node) {
    node->clearStatistics();
  }
};

template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;

template <class Node>
class runTest {
  public:
  static void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) {
    using std::endl;
  
    ThreadedBlasKiller<Node>::kill();

    typedef double Scalar;
    typedef int LO; //LocalOrdinal
    typedef int GO; //GlobalOrdinal
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
    typedef Tpetra::Operator<Scalar,LO,GO,Node>    TOP;
    typedef Anasazi::Eigenproblem<Scalar,TMV,TOP>    AEP;
    typedef Anasazi::SolverManager<Scalar,TMV,TOP>   ASM;
    typedef Anasazi::MultiVecTraits<Scalar,TMV>      AMVT;

    AMVT::mvTimesMatAddMvTimer_ = Teuchos::TimeMonitor::getNewTimer("Anasazi/Tpetra::MvTimesMatAddMv()");
    AMVT::mvTransMvTimer_ = Teuchos::TimeMonitor::getNewTimer("Anasazi/Tpetra::MvTransMv()");

    const bool IAmRoot = (comm->getRank() == 0);

    RCP<FancyOStream> fos = fancyOStream(rcpFromRef(std::cout));
    fos->setShowProcRank(true);
    *fos << "Running test with Node == " << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << endl;

    int myWeight = myMachPL.get<int>("Node Weight",1);

    // The build_problem function is located in build_problem.hpp.
    // Note that build_problem calls build_precond and sets a preconditioner on the
    // linear-problem, if a preconditioner is specified.
    RCP<AEP> problem = build_eigproblem<Scalar,LO,GO,Node>(AnasaziTpetraHybridDriverTest::plTestParams, comm, node, myWeight);

    // The build_solver function is located in build_solver.hpp:
    RCP<ASM> solver = build_eigsolver<Scalar,TMV,TOP>(comm, AnasaziTpetraHybridDriverTest::plTestParams, problem);
    if (IAmRoot) *fos << Teuchos::typeName(*solver) << endl;

    NodeDetails<Node>::clear(node);

    Anasazi::ReturnType ret = solver->solve();

    if (IAmRoot) {
      *fos << "Converged in " << solver->getNumIters() << " iterations." << endl;
    }

    // RCP<TMV> R = rcp(new TMV(*problem->getRHS()));
    // problem->computeCurrResVec(&*R, &*problem->getLHS(), &*problem->getRHS());
    // Array<ScalarTraits<Scalar>::magnitudeType> norms(R->getNumVectors());
    // R->norm2(norms);
    // if (norms.size() < 1) {
    //   throw std::runtime_error("ERROR: norms.size()==0 indicates R->getNumVectors()==0.");
    // }
    // if (IAmRoot) {
    //   *fos << "2-Norm of 0th residual vec: " << norms[0] << endl;
    // }

    if (IAmRoot) {
      *fos << endl;
    }
    Teuchos::TimeMonitor::summarize(*fos,true,false,false);
    if (IAmRoot) {
      *fos << endl;
    }

    // print node details
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        NodeDetails<Node>::printDetails(fos,node.getConst());
      }
      comm->barrier();
    }
  }
};

int main(int argc, char*argv[])
{
  Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));

  const int numImages = comm->getSize(),
            myImageID = comm->getRank();
  const bool IAmRoot = (myImageID == 0);

  // Time timer("Total time for AnasaziTpetraHybridDriver test");
  // timer.start();

  // Get macine and test parameter files from command line
  string fnTestParams("");
  string fnMachine("");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("param-file", &fnTestParams, "XML Parameters file");
  cmdp.setOption("machine-file",&fnMachine,"Filename for XML machine description file.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Read XML machine file into a ParameterList
  if (IAmRoot) {
    cout << "Every proc machine parameters from: " << fnMachine << endl;
  }
  ParameterList plMachine = ParameterXMLFileReader(fnMachine).getParameters();

  // Read the contents of the xml file into a ParameterList.
  if (IAmRoot) {
    cout << "Every proc reading parameters from: " << fnTestParams << endl;
  }
  AnasaziTpetraHybridDriverTest::plTestParams = ParameterXMLFileReader(fnTestParams).getParameters();

  // Everything is node based from here on out; have the HybridPlatform launch us appropriately
  Tpetra::HybridPlatform platform(comm,plMachine);
  platform.runUserCode<runTest>();

  // timer.stop();
  // if (IAmRoot) {
  //   cout << "proc 0 total program time: " << timer.totalElapsedTime() << endl;
  // }
  
  return 0;
}
