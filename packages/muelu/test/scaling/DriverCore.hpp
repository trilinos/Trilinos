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
#ifndef DRIVERCORE_HPP
#include <MueLu.hpp>

// Teuchos
#include <Teuchos_ScalarTraits.hpp>

// Xpetra
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Matrix.hpp>


// Belos
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
#include <Xpetra_TpetraOperator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
namespace BelosTpetra {
namespace Impl {
  extern void register_Cg (const bool verbose);
  extern void register_CgPipeline (const bool verbose);
  extern void register_CgSingleReduce (const bool verbose);
  extern void register_Gmres (const bool verbose);
  extern void register_GmresS (const bool verbose);
  extern void register_GmresSstep (const bool verbose);
  extern void register_GmresSingleReduce (const bool verbose);
} // namespace Impl
} // namespace BelosTpetra
#endif

#ifdef HAVE_MUELU_EPETRA
#include <BelosEpetraAdapter.hpp>    // => This header defines Belos::EpetraPrecOp
#endif
#endif

// Cuda
#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// AMGX
#ifdef HAVE_MUELU_AMGX
#include <MueLu_AMGXOperator.hpp>
#endif


#include <MueLu_CreateXpetraPreconditioner.hpp>



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
// A handy macro to switch time monitors in a StackedTimer-compatible way
#define MUELU_SWITCH_TIME_MONITOR(tm,timername)                                  \
  {tm=Teuchos::null; tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(timername)));}


//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
// Support for ML interface
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
#include <Xpetra_EpetraOperator.hpp>
#include <ml_MultiLevelPreconditioner.h>

// Helper functions for compilation purposes
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct ML_Wrapper{
  static void Generate_ML_MultiLevelPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & A, Teuchos::ParameterList & mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & mlopX) {
    throw std::runtime_error("Template parameter mismatch");
  }
};


template<class GlobalOrdinal>
struct ML_Wrapper<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> {
  static void Generate_ML_MultiLevelPreconditioner(Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& A,Teuchos::ParameterList & mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& mlopX) {
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    Teuchos::RCP<const Epetra_CrsMatrix> Aep   = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(A);
    Teuchos::RCP<Epetra_Operator> mlop  = Teuchos::rcp<Epetra_Operator>(new ML_Epetra::MultiLevelPreconditioner(*Aep,mueluList));
#if defined(HAVE_MUELU_BELOS)
    // NOTE: Belos needs the Apply() and AppleInverse() routines of ML swapped.  So...
    mlop = Teuchos::rcp<Belos::EpetraPrecOp>(new Belos::EpetraPrecOp(mlop));
#endif

    mlopX = Teuchos::rcp(new Xpetra::EpetraOperator<GO,NO>(mlop));
  }
};
#endif


//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
// This is a standard setup routine
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PreconditionerSetup(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & A,
                         Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node> > & coordinates,
                         Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & nullspace,
                         Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & material,
                         Teuchos::ParameterList & mueluList,
                         bool profileSetup,
                         bool useAMGX,
                         bool useML,
                         bool setNullSpace,
                         int numRebuilds,
                         Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & H,
                         Teuchos::RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & Prec) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();
  typedef typename Teuchos::ScalarTraits<SC>::coordinateType coordinate_type;
  typedef Xpetra::MultiVector<coordinate_type,LO,GO,NO> CoordinateMultiVector;
  // =========================================================================
  // Preconditioner construction
  // =========================================================================
#ifdef HAVE_MUELU_CUDA
  if(profileSetup) cudaProfilerStart();
#endif

   if(useML && lib != Xpetra::UseEpetra) throw std::runtime_error("Error: Cannot use ML on non-epetra matrices");

   for (int i = 0; i <= numRebuilds; i++) {
     A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
     if(useAMGX) {
#if defined(HAVE_MUELU_AMGX) && defined(HAVE_MUELU_TPETRA)
       RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > Ac      = Utilities::Op2NonConstTpetraCrs(A);
       RCP<Tpetra::Operator<SC,LO,GO,NO> > At       = Teuchos::rcp_dynamic_cast<Tpetra::Operator<SC,LO,GO,NO> >(Ac);
       RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > Top = MueLu::CreateTpetraPreconditioner(At, mueluList);
       Prec = Teuchos::rcp(new Xpetra::TpetraOperator<SC,LO,GO,NO>(Top));
#endif
     } else if(useML) {
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
       mueluList.remove("use external multigrid package");
       if(!coordinates.is_null()) {
         RCP<const Epetra_MultiVector> epetraCoord =  MueLu::Utilities<coordinate_type,LO,GO,NO>::MV2EpetraMV(coordinates);
         if(epetraCoord->NumVectors() > 0)  mueluList.set("x-coordinates",(*epetraCoord)[0]);
         if(epetraCoord->NumVectors() > 1)  mueluList.set("y-coordinates",(*epetraCoord)[1]);
         if(epetraCoord->NumVectors() > 2)  mueluList.set("z-coordinates",(*epetraCoord)[2]);
       }
       if(!material.is_null()) {
         RCP<const Epetra_MultiVector> epetraMat =  MueLu::Utilities<SC,LO,GO,NO>::MV2EpetraMV(material);
         mueluList.set("material coordinates",(*epetraMat)[0]);
       }
       ML_Wrapper<SC, LO, GO, NO>::Generate_ML_MultiLevelPreconditioner(A,mueluList,Prec);
#endif
     }
     else {
       Teuchos::Array<LO> lNodesPerDim(3, 10);
       Teuchos::ParameterList& userParamList = mueluList.sublist("user data");
       if(!coordinates.is_null())
         userParamList.set<RCP<CoordinateMultiVector> >("Coordinates", coordinates);
       if (!nullspace.is_null() && setNullSpace )
         userParamList.set<RCP<Xpetra::MultiVector<SC,LO,GO,NO>> >("Nullspace", nullspace);
       userParamList.set<Teuchos::Array<LO> >("Array<LO> lNodesPerDim", lNodesPerDim);
       H = MueLu::CreateXpetraPreconditioner(A, mueluList);
     }
   }
#ifdef HAVE_MUELU_CUDA
   if(profileSetup) cudaProfilerStop();
#endif
}


#if defined(HAVE_MUELU_EPETRA)

// Helper functions for compilation purposes
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct Matvec_Wrapper{
  static void UnwrapEpetra(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                           Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& X,
                           Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& B,
                           Teuchos::RCP<const Epetra_CrsMatrix>& Aepetra,
                           Teuchos::RCP<Epetra_MultiVector>& Xepetra,
                           Teuchos::RCP<Epetra_MultiVector>& Bepetra) {
    throw std::runtime_error("Template parameter mismatch");
  }
};


template<class GlobalOrdinal>
struct Matvec_Wrapper<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> {
  static void UnwrapEpetra(Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& A,
                           Teuchos::RCP<Xpetra::MultiVector<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& X,
                           Teuchos::RCP<Xpetra::MultiVector<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& B,
                           Teuchos::RCP<const Epetra_CrsMatrix>& Aepetra,
                           Teuchos::RCP<Epetra_MultiVector>& Xepetra,
                           Teuchos::RCP<Epetra_MultiVector>& Bepetra) {
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    Aepetra = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(A);
    Xepetra = Teuchos::rcp(& Xpetra::toEpetra(*X),false);
    Bepetra = Teuchos::rcp(& Xpetra::toEpetra(*B),false);
  }
};
#endif

//*************************************************************************************
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SystemSolve(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & A,
                 Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & X,
                 Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & B,
                 Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node > > & H,
                 Teuchos::RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node > > & Prec,
                 Teuchos::FancyOStream & out,
                 std::string solveType,
                 std::string belosType,
                 bool profileSolve,
                 bool useAMGX,
                 bool useML,
                 int cacheSize,
                 int numResolves,
                 bool scaleResidualHist,
                 bool solvePreconditioned,
                 int maxIts,
                 double tol) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero();

  // Cache clearing
  std::vector<int> tempVector;
  int min = 0, max = 10;
  int numInts = 0;
  if (cacheSize > 0) {
    cacheSize *= 1024; //convert to bytes
    numInts = cacheSize/sizeof(int) + 1;
    tempVector.resize(numInts);
  }

  // Get the raw matrices for matvec testing
#if defined(HAVE_MUELU_TPETRA)
      Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > Atpetra;
      Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO> > Xtpetra,Btpetra;
      if(lib==Xpetra::UseTpetra) {
        Atpetra = Utilities::Op2NonConstTpetraCrs(A);
        Xtpetra = rcp(& Xpetra::toTpetra(*X),false);
        Btpetra = rcp(& Xpetra::toTpetra(*B),false);
      }
#endif
#if defined(HAVE_MUELU_EPETRA)
      Teuchos::RCP<const Epetra_CrsMatrix> Aepetra;
      Teuchos::RCP<Epetra_MultiVector> Xepetra,Bepetra;
      if(lib==Xpetra::UseEpetra) {
        Matvec_Wrapper<SC,LO,GO,NO>::UnwrapEpetra(A,X,B,Aepetra,Xepetra,Bepetra);
      }
#endif

  for(int solveno = 0; solveno<=numResolves; solveno++) {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - LHS and RHS initialization")));
    X->putScalar(zero);
    tm = Teuchos::null;

    if (solveType == "none") {
      // Do not perform a solve
    } else if (solveType == "matvec") {
      // Just do matvecs
      tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 6 - Matvec")));
#if defined(HAVE_MUELU_TPETRA)
      if(lib==Xpetra::UseTpetra) Atpetra->apply(*Btpetra,*Xtpetra);
#endif
#if defined(HAVE_MUELU_EPETRA) && !defined(HAVE_MUELU_INST_COMPLEX_INT_INT) && !defined(HAVE_MUELU_INST_FLOAT_INT_INT)
      if(lib==Xpetra::UseEpetra) Aepetra->Apply(*Bepetra,*Xepetra);
#endif
      //clear the cache (and don't time it)
      tm = Teuchos::null;
      int ttt = rand();
      for (int i=0; i<numInts; ++i)
        tempVector[i] += (min + (ttt % static_cast<int>(max - min + 1)));
    } else if (solveType == "standalone") {
      tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Fixed Point Solve")));
#ifdef HAVE_MUELU_CUDA
      if(profileSolve) cudaProfilerStart();
#endif
      if(useAMGX) {
#if defined(HAVE_MUELU_AMGX) && defined(HAVE_MUELU_TPETRA)
        // Do a fixed-point iteraiton without convergence checks
        RCP<MultiVector> R = MultiVectorFactory::Build(X->getMap(), X->getNumVectors());
        for(int i=0; i<maxIts; i++) {
          Utilities::Residual(*A, *X, *B, *R);
          Prec->apply(*R,*X);
        }
#endif
      }
      else {
        H->IsPreconditioner(false);
        std::pair<LocalOrdinal,typename STS::magnitudeType> maxItsTol(maxIts, tol);
        H->Iterate(*B, *X, maxItsTol);
      }
#ifdef HAVE_MUELU_CUDA
      if(profileSolve) cudaProfilerStop();
#endif
    } else if (solveType == "belos") {
#ifdef HAVE_MUELU_BELOS
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Belos Solve")));
#ifdef HAVE_MUELU_CUDA
      if(profileSolve) cudaProfilerStart();
#endif
      // Operator and Multivector type that will be used with Belos
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;

      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
      Teuchos::RCP<OP> belosPrec; // Turns a MueLu::Hierarchy object into a Belos operator
      if(useAMGX) {
#if defined(HAVE_MUELU_AMGX) && defined(HAVE_MUELU_TPETRA)
        belosPrec = rcp(new Belos::XpetraOp <SC, LO, GO, NO>(Prec)); // Turns an Xpetra::Operator object into a Belos operator
#endif
      }
      else if(useML) {
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
        belosPrec = rcp(new Belos::XpetraOp <SC, LO, GO, NO>(Prec)); // Turns an Xpetra::Operator object into a Belos operator
#endif
      }
      else {
        H->IsPreconditioner(true);
        belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H)); // Turns a MueLu::Hierarchy object into a Belos operator
      }

      // Belos parameter list
      RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
      belosList->set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
      belosList->set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
      belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList->set("Output Frequency",      1);
      belosList->set("Output Style",          Belos::Brief);
      if (!scaleResidualHist)
        belosList->set("Implicit Residual Scaling", "None");

      int numIts;
      Belos::ReturnType ret = Belos::Unconverged;

#if defined(HAVE_MUELU_TPETRA)
      constexpr bool verbose = false;
      BelosTpetra::Impl::register_Cg (verbose);
      BelosTpetra::Impl::register_CgPipeline (verbose);
      BelosTpetra::Impl::register_CgSingleReduce (verbose);
      BelosTpetra::Impl::register_Gmres (verbose);
      BelosTpetra::Impl::register_GmresS (verbose);
      BelosTpetra::Impl::register_GmresSstep (verbose);
      BelosTpetra::Impl::register_GmresSingleReduce (verbose);

      try {

        TEUCHOS_TEST_FOR_EXCEPTION(lib!=Xpetra::UseTpetra, std::invalid_argument, "Need to use Tpetra backend in order to call a Belos Tpetra solver.");

        using tMV = Tpetra::MultiVector<SC, LO, GO, NO>;
        using tOP = Tpetra::Operator<SC, LO, GO, NO>;

        Teuchos::RCP<tOP> belosPrecTpetra;
        if(useAMGX) {
#if defined(HAVE_MUELU_AMGX) && defined(HAVE_MUELU_TPETRA)
          RCP<Xpetra::TpetraOperator<SC,LO,GO,NO> > xto = Teuchos::rcp_dynamic_cast<Xpetra::TpetraOperator<SC,LO,GO,NO> >(Prec);
          belosPrecTpetra = xto->getOperator();
#endif
        }
        else {
          belosPrecTpetra = rcp(new MueLu::TpetraOperator<SC,LO,GO,NO>(H));
        }

        // Construct a Belos LinearProblem object
        RCP<Belos::LinearProblem<SC, tMV, tOP> > belosProblem = rcp(new Belos::LinearProblem<SC, tMV, tOP>(Atpetra, Xtpetra, Btpetra));
        if(solvePreconditioned) belosProblem->setRightPrec(belosPrecTpetra);

        bool set = belosProblem->setProblem();
        if (set == false) {
          throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
        }

        // Create an iterative solver manager
        Belos::SolverFactory<SC, tMV, tOP> solverFactory;
        RCP< Belos::SolverManager<SC, tMV, tOP> > solver = solverFactory.create(belosType, belosList);
        solver->setProblem(belosProblem);

        // Perform solve
        ret = solver->solve();
        numIts = solver->getNumIters();

      } catch (std::invalid_argument&)
#endif
      {

        // Construct a Belos LinearProblem object
        RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
        if(solvePreconditioned) belosProblem->setRightPrec(belosPrec);

        bool set = belosProblem->setProblem();
        if (set == false) {
          throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
        }

        // Create an iterative solver manager
        Belos::SolverFactory<SC, MV, OP> solverFactory;
        RCP< Belos::SolverManager<SC, MV, OP> > solver = solverFactory.create(belosType, belosList);
        solver->setProblem(belosProblem);

        // Perform solve
        ret = solver->solve();
        numIts = solver->getNumIters();

      }

      // Get the number of iterations for this solve.
      out << "Number of iterations performed for this solve: " << numIts << std::endl;

      // Check convergence
      if (ret != Belos::Converged)
        out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
      else
        out << std::endl << "SUCCESS:  Belos converged!" << std::endl;
#ifdef HAVE_MUELU_CUDA
      if(profileSolve) cudaProfilerStop();
#endif
#endif //ifdef HAVE_MUELU_BELOS
    } else {
      throw MueLu::Exceptions::RuntimeError("Unknown solver type: \"" + solveType + "\"");
    }
  }// end resolves

}

#endif
