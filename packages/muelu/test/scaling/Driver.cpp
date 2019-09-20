// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MatrixLoad.hpp>
#include <DriverCore.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#ifdef HAVE_MUELU_TPETRA
#include <BelosTpetraAdapter.hpp>    // => This header defines Belos::TpetraOp
#endif
#ifdef HAVE_MUELU_EPETRA
#include <BelosEpetraAdapter.hpp>    // => This header defines Belos::EpetraPrecOp
#endif
#endif


#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

#ifdef HAVE_MUELU_AMGX
#include <MueLu_AMGXOperator.hpp>
#include <MueLu_AMGX_Setup.hpp>
#endif
#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include "Xpetra_TpetraMultiVector.hpp"
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>
#endif

#ifdef HAVE_MUELU_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#endif


/*********************************************************************/

#ifdef HAVE_MUELU_TPETRA
#include "KokkosBlas1_abs_impl.hpp"
template<class RV, class XV, class SizeType>
void Temporary_Replacement_For_Kokkos_abs(const RV& R, const XV& X) {
  typedef typename XV::execution_space execution_space;
  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  KokkosBlas::Impl::V_Abs_Functor<RV, XV, SizeType> op (R, X);
  Kokkos::parallel_for (policy, op);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void equilibrateMatrix(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &Axpetra, std::string equilibrate) {
#include <MueLu_UseShortNames.hpp>
  using Tpetra::computeRowAndColumnOneNorms;
  using Tpetra::leftAndOrRightScaleCrsMatrix;
  bool equilibrate_1norm = (equilibrate == "1-norm");
  bool equilibrate_diag  = (equilibrate == "diag");
  bool equilibrate_no    = (equilibrate == "no");
  bool assumeSymmetric = false;
  typedef typename Tpetra::Details::EquilibrationInfo<typename Kokkos::ArithTraits<Scalar>::val_type,typename Node::device_type> equil_type;

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A = Utilities::Op2NonConstTpetraCrs(Axpetra);

  if(Axpetra->getRowMap()->lib() == Xpetra::UseTpetra) {
     equil_type equibResult_ = computeRowAndColumnOneNorms (*A, assumeSymmetric);
     if (equilibrate_1norm) {
        using device_type = typename Node::device_type;
        using mag_type = typename Kokkos::ArithTraits<Scalar>::mag_type;
        using mag_view_type = Kokkos::View<mag_type*, device_type>;
        using scalar_view_type = Kokkos::View<typename equil_type::val_type*, device_type>;

        mag_view_type rowDiagAbsVals ("rowDiagAbsVals",equibResult_.rowDiagonalEntries.extent (0));
        //        KokkosBlas::abs (rowDiagAbsVals, equibResult_.rowDiagonalEntries);
        Temporary_Replacement_For_Kokkos_abs<mag_view_type,scalar_view_type,LocalOrdinal>(rowDiagAbsVals, equibResult_.rowDiagonalEntries);

        mag_view_type colDiagAbsVals ("colDiagAbsVals",equibResult_.colDiagonalEntries.extent (0));

        //        KokkosBlas::abs (colDiagAbsVals, equibResult_.colDiagonalEntries);
        Temporary_Replacement_For_Kokkos_abs<mag_view_type,scalar_view_type,LocalOrdinal>(colDiagAbsVals, equibResult_.colDiagonalEntries);

        leftAndOrRightScaleCrsMatrix (*A, rowDiagAbsVals, colDiagAbsVals,
                                      true, true, equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
     }
     else if (equilibrate_diag) {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        leftAndOrRightScaleCrsMatrix (*A, equibResult_.rowNorms,
                                      colScalingFactors, true, true,
                                      equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
     else if (equilibrate_no) {
       // no-op
     }
     else
       throw std::runtime_error("Invalid 'equilibrate' option '"+equilibrate+"'");
  }
}
#endif


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::TimeMonitor;
  using Teuchos::ParameterList;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName       = "";                 clp.setOption("xml",                   &xmlFileName,       "read parameters from an xml file");
  std::string yamlFileName      = "";                 clp.setOption("yaml",                  &yamlFileName,      "read parameters from a yaml file");
  bool        printTimings      = true;               clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
  std::string timingsFormat     = "table-fixed";      clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
  int         writeMatricesOPT  = -2;                 clp.setOption("write",                 &writeMatricesOPT,  "write matrices to file (-1 means all; i>=0 means level i)");
  std::string dsolveType        = "belos", solveType; clp.setOption("solver",                &dsolveType,        "solve type: (none | cg | gmres | standalone | matvec)");
  std::string belosType         = "cg";               clp.setOption("belosType",             &belosType,         "belos solver type: (Pseudoblock CG | Block CG | Pseudoblock GMRES | Block GMRES | ...) see BelosSolverFactory.hpp for exhaustive list of solvers");
  double      dtol              = 1e-12, tol;         clp.setOption("tol",                   &dtol,              "solver convergence tolerance");
  bool        binaryFormat      = false;              clp.setOption("binary", "ascii",       &binaryFormat,      "print timings to screen");

  std::string rowMapFile;                             clp.setOption("rowmap",                &rowMapFile,        "map data file");
  std::string colMapFile;                             clp.setOption("colmap",                &colMapFile,        "colmap data file");
  std::string domainMapFile;                          clp.setOption("domainmap",             &domainMapFile,     "domainmap data file");
  std::string rangeMapFile;                           clp.setOption("rangemap",              &rangeMapFile,      "rangemap data file");
  std::string matrixFile;                             clp.setOption("matrix",                &matrixFile,        "matrix data file");
  std::string rhsFile;                                clp.setOption("rhs",                   &rhsFile,           "rhs data file");
  std::string coordFile;                              clp.setOption("coords",                &coordFile,         "coordinates data file");
  std::string nullFile;                               clp.setOption("nullspace",             &nullFile,          "nullspace data file");
  int         numRebuilds       = 0;                  clp.setOption("rebuild",               &numRebuilds,       "#times to rebuild hierarchy");
  int         numResolves       = 0;                  clp.setOption("resolve",               &numResolves,       "#times to redo solve");
  int         maxIts            = 200;                clp.setOption("its",                   &maxIts,            "maximum number of solver iterations");
  int         numVectors        = 1;                  clp.setOption("multivector",           &numVectors,        "number of rhs to solve simultaneously");
  bool        scaleResidualHist = true;               clp.setOption("scale", "noscale",      &scaleResidualHist, "scaled Krylov residual history");
  bool        solvePreconditioned = true;             clp.setOption("solve-preconditioned","no-solve-preconditioned", &solvePreconditioned, "use MueLu preconditioner in solve");  
  bool        useStackedTimer   = false;              clp.setOption("stacked-timer","no-stacked-timer", &useStackedTimer, "use stacked timer");

#ifdef HAVE_MUELU_TPETRA
  std::string equilibrate = "no" ;                    clp.setOption("equilibrate",           &equilibrate,       "equilibrate the system (no | diag | 1-norm)");
#endif
#ifdef HAVE_MUELU_CUDA
  bool profileSetup = false;                          clp.setOption("cuda-profile-setup", "no-cuda-profile-setup", &profileSetup, "enable CUDA profiling for setup");
  bool profileSolve = false;                          clp.setOption("cuda-profile-solve", "no-cuda-profile-solve", &profileSolve, "enable CUDA profiling for solve");
#else
  bool profileSetup = false;
  bool profileSolve = false;
#endif
  int  cacheSize = 0;                                 clp.setOption("cachesize",               &cacheSize,       "cache size (in KB)");
#ifdef HAVE_MPI
  int provideNodeComm = 0;                            clp.setOption("nodecomm",          &provideNodeComm,  "make the nodal communicator available w/ reduction factor X");
#endif

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error,
                             "Cannot provide both xml and yaml input files");

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  ParameterList paramList;
  auto inst = xpetraParameters.GetInstantiation();

  if (yamlFileName != "") {
    Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  } else {
    if (inst == Xpetra::COMPLEX_INT_INT)
      xmlFileName = (xmlFileName != "" ? xmlFileName : "scaling-complex.xml");
    else
      xmlFileName = (xmlFileName != "" ? xmlFileName : "scaling.xml");
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  }

  if (inst == Xpetra::COMPLEX_INT_INT && dsolveType == "belos") {
    belosType = "gmres";
    out << "WARNING: CG will not work with COMPLEX scalars, switching to GMRES"<<std::endl;
  }

#ifdef HAVE_MUELU_AMGX
//Initialize AMGX
MueLu::MueLu_AMGX_initialize();
MueLu::MueLu_AMGX_initialize_plugins();
#endif

  bool isDriver = paramList.isSublist("Run1");
  if (isDriver) {
    // update galeriParameters with the values from the XML file
    ParameterList& realParams = galeriParameters.GetParameterList();

    for (ParameterList::ConstIterator it = realParams.begin(); it != realParams.end(); it++) {
      const std::string& name = realParams.name(it);
      if (paramList.isParameter(name))
        realParams.setEntry(name, paramList.getEntry(name));
    }
  }

#ifdef HAVE_MPI
  // Generate the node-level communicator, if we want one
  Teuchos::RCP<const Teuchos::Comm<int> > nodeComm;
  int NodeId = comm->getRank();
  if(provideNodeComm) {
    nodeComm = MueLu::GenerateNodeComm(comm,NodeId,provideNodeComm);
    //    printf("DEBUG: Base rank %d => New, node %d, rank %d\n",comm->getRank(),NodeId,nodeComm->getRank());
  }
#endif

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;
#ifdef HAVE_MUELU_OPENMP
  std::string node_name = Node::name();
  if(!comm->getRank() && !node_name.compare("OpenMP/Wrapper"))
    galeriStream<<"OpenMP Max Threads = "<<omp_get_max_threads()<<std::endl;
#endif


  comm->barrier();
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  if(useStackedTimer)
    stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));


  RCP<Matrix>      A;
  RCP<const Map>   map;
  RCP<RealValuedMultiVector> coordinates;
  RCP<Xpetra::MultiVector<SC,LO,GO,NO> > nullspace;
  RCP<MultiVector> X, B;

  // Load the matrix off disk (or generate it via Galeri)
  MatrixLoad<SC,LO,GO,NO>(comm,lib,binaryFormat,matrixFile,rhsFile,rowMapFile,colMapFile,domainMapFile,rangeMapFile,coordFile,nullFile,map,A,coordinates,nullspace,X,B,numVectors,galeriParameters,xpetraParameters,galeriStream);
  comm->barrier();
  tm = Teuchos::null;

  // Do equilibration if requested
#ifdef HAVE_MUELU_TPETRA
  if(lib == Xpetra::UseTpetra) {
    equilibrateMatrix(A,equilibrate);
  }
#endif

  int numReruns = 1;
  if (paramList.isParameter("number of reruns"))
    numReruns = paramList.get<int>("number of reruns");

  for (int rerunCount = 1; rerunCount <= numReruns; rerunCount++) {
    bool stop = false;
    ParameterList mueluList, runList;
    const bool mustAlreadyExist = true;
    if (isDriver) {
      runList   = paramList.sublist("Run1",  mustAlreadyExist);
      mueluList = runList  .sublist("MueLu", mustAlreadyExist);
    } else {
      mueluList = paramList;
      stop = true;
    }

    if (nullspace.is_null()) {
      int blkSize = 1;
      if (mueluList.isSublist("Matrix")) {
        // Factory style parameter list
        const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
        if (operatorList.isParameter("PDE equations"))
          blkSize = operatorList.get<int>("PDE equations");

      } else if (paramList.isParameter("number of equations")) {
        // Easy style parameter list
        blkSize = paramList.get<int>("number of equations");
      }

      nullspace = MultiVectorFactory::Build(map, blkSize);
      for (int i = 0; i < blkSize; i++) {
        RCP<const Map> domainMap = A->getDomainMap();
        GO             indexBase = domainMap->getIndexBase();

        ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
        for (int j = 0; j < nsData.size(); j++) {
          GO GID = domainMap->getGlobalElement(j) - indexBase;

          if ((GID-i) % blkSize == 0)
            nsData[j] = Teuchos::ScalarTraits<SC>::one();
        }
      }
    }


    int runCount = 1;
    int   savedOut  = -1;
    FILE* openedOut = NULL;
    do {
      solveType = dsolveType;
      tol       = dtol;


      if (isDriver) {
        if (runList.isParameter("filename")) {
          // Redirect all output into a filename We have to redirect all output,
          // including printf's, therefore we cannot simply replace C++ cout
          // buffers, and have to use heavy machinary (dup2)
          std::string filename = runList.get<std::string>("filename");
          if (numReruns > 1)
            filename += "_run" + MueLu::toString(rerunCount);
          filename += (lib == Xpetra::UseEpetra ? ".epetra" : ".tpetra");

          savedOut  = dup(STDOUT_FILENO);
          openedOut = fopen(filename.c_str(), "w");
          dup2(fileno(openedOut), STDOUT_FILENO);
        }
        if (runList.isParameter("solver")) solveType = runList.get<std::string>("solver");
        if (runList.isParameter("tol"))    tol       = runList.get<double>     ("tol");
      }

      RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream& out2 = *fancy2;
      out2.setOutputToRootOnly(0);
      out2 << galeriStream.str();

      // =========================================================================
      // Preconditioner construction
      // =========================================================================
      bool useAMGX = mueluList.isParameter("use external multigrid package") && (mueluList.get<std::string>("use external multigrid package") == "amgx");
      bool useML   = mueluList.isParameter("use external multigrid package") && (mueluList.get<std::string>("use external multigrid package") == "ml");
#ifdef HAVE_MPI
      if(provideNodeComm && !useAMGX && !useML) {
        Teuchos::ParameterList& userParamList = mueluList.sublist("user data");
        userParamList.set("Node Comm",nodeComm);
      }
#endif
      out2<<"*********** MueLu ParameterList ***********"<<std::endl;
      out2<<mueluList;
      out2<<"*******************************************"<<std::endl;

      RCP<Hierarchy> H;
      RCP<Operator> Prec;
      bool preconditionerOK = true;
      try {
        comm->barrier();
        // Build the preconditioner numRebuilds+1 times
        MUELU_SWITCH_TIME_MONITOR(tm,"Driver: 2 - MueLu Setup");
        PreconditionerSetup(A,coordinates,nullspace,mueluList,profileSetup,useAMGX,useML,numRebuilds,H,Prec);

        comm->barrier();
        tm = Teuchos::null;
      }
      catch(const std::exception& e) {
        out2<<"MueLu_Driver: preconditioner setup crashed w/ message:"<<e.what()<<std::endl;
        H=Teuchos::null; Prec=Teuchos::null;
        preconditionerOK = false;
      }

      // =========================================================================
      // System solution (Ax = b)
      // =========================================================================
      if(preconditionerOK) {
        try {
          comm->barrier();
          if (writeMatricesOPT > -2) {
            tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.5 - Matrix output")));
            H->Write(writeMatricesOPT, writeMatricesOPT);
            tm = Teuchos::null;
          }

          // Solve the system numResolves+1 times
          SystemSolve(A,X,B,H,Prec,out2,solveType,belosType,profileSolve,useAMGX,useML,cacheSize,numResolves,scaleResidualHist,solvePreconditioned,maxIts,tol);

          comm->barrier();
        }
        catch(const std::exception& e) {
          out2<<"MueLu_Driver: solver crashed w/ message:"<<e.what()<<std::endl;
        }
      }
      else {
        out2<<"MueLu_Driver: Not solving system due to crash in preconditioner setup"<<std::endl;
      }


      tm = Teuchos::null;
      globalTimeMonitor = Teuchos::null;

      if (printTimings) {
        RCP<ParameterList> reportParams = rcp(new ParameterList);
        if (timingsFormat == "yaml") {
          reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
          reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
        }
        reportParams->set("How to merge timer sets",   "Union");
        reportParams->set("alwaysWriteLocal",          false);
        reportParams->set("writeGlobalStats",          true);
        reportParams->set("writeZeroTimers",           false);
        // FIXME: no "ignoreZeroTimers"

        const std::string filter = "";

        if (useStackedTimer) {
          Teuchos::StackedTimer::OutputOptions options;
          options.output_fraction = options.output_histogram = options.output_minmax = true;
          stacked_timer->report(out2, comm, options);
        }
        else {
          std::ios_base::fmtflags ff(out2.flags());
          if (timingsFormat == "table-fixed") out2 << std::fixed;
          else                                out2 << std::scientific;
          TimeMonitor::report(comm.ptr(), out, filter, reportParams);
          out2 << std::setiosflags(ff);
        }
      }

      TimeMonitor::clearCounters();
      out2 << std::endl;

      if (isDriver) {
        try {
          runList   = paramList.sublist("Run" + MueLu::toString(++runCount), mustAlreadyExist);
          mueluList = runList  .sublist("MueLu", mustAlreadyExist);
        } catch (Teuchos::Exceptions::InvalidParameterName& e) {
          stop = true;
        }
      }
      fflush(NULL);
      comm->barrier();
    } while (!stop);

    // Cleanup Output
    if (openedOut != NULL) {
      TEUCHOS_ASSERT(savedOut >= 0);
      dup2(savedOut, STDOUT_FILENO);
      fclose(openedOut);
      openedOut = NULL;
    }


  }//end reruns

#ifdef HAVE_MUELU_AMGX
// Finalize AMGX
MueLu::MueLu_AMGX_finalize_plugins();
MueLu::MueLu_AMGX_finalize();
#endif

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
