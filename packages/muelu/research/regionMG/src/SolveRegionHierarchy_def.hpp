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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SOLVEREGIONHIERARCHY_DEF_HPP
#define MUELU_SOLVEREGIONHIERARCHY_DEF_HPP

#include "SetupRegionHierarchy_def.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;

//! Recursive V-cycle in region fashion
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void vCycle(const int l, ///< ID of current level
            const int numLevels, ///< Total number of levels
            const std::string cycleType,
            RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & regHierarchy,
            RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& fineRegX, ///< solution
            RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > fineRegB, ///< right hand side
            Array<RCP<Teuchos::ParameterList> > smootherParams, ///< region smoother parameter list
            bool& zeroInitGuess,
            RCP<ParameterList> coarseSolverData = Teuchos::null,
            RCP<ParameterList> hierarchyData = Teuchos::null)
{
#include "MueLu_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  RCP<Level> level = regHierarchy->GetLevel(l);
  RCP<Matrix>    regMatrix      = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());
  RCP<const Map> regRowMap      = regMatrix->getRowMap();
  RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > regRowImporter = level->Get<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("rowImport");
  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regInterfaceScalings = level->Get<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("regInterfaceScalings");

  int cycleCount = 1;
  if(cycleType == "W" && l > 0) // W cycle and not on finest level
    cycleCount=2;
  if (l < numLevels - 1) { // fine or intermediate levels
    for(int cycle=0; cycle < cycleCount; cycle++){

//    std::cout << "level: " << l << std::endl;

      // extract data from hierarchy parameterlist
      std::string levelName("level" + std::to_string(l));
      ParameterList levelList;
      bool useCachedVectors = false;
      // if(Teuchos::nonnull(hierarchyData) &&  hierarchyData->isSublist(levelName)) {
      //   levelList = hierarchyData->sublist(levelName);
      //   if(levelList.isParameter("residual") && levelList.isParameter("solution")) {
      //     useCachedVectors = true;
      //   }
      // }

      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 1 - pre-smoother")));

      // pre-smoothing
      smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrix,
                    regRowMap, regRowImporter, zeroInitGuess);

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 2 - compute residual")));

      RCP<Vector> regRes;
      if(useCachedVectors) {
        regRes = levelList.get<RCP<Vector> >("residual");
      } else {
        regRes = VectorFactory::Build(regRowMap, true);
      }
      computeResidual(regRes, fineRegX, fineRegB, regMatrix, *smootherParams[l]);

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 3 - scale interface")));

      scaleInterfaceDOFs(regRes, regInterfaceScalings, true);

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 4 - create coarse vectors")));

      // Transfer to coarse level
      RCP<Vector> coarseRegX;
      RCP<Vector> coarseRegB;

      {
        RCP<Level> levelCoarse = regHierarchy->GetLevel(l+1);
        RCP<Matrix> regProlongCoarse = levelCoarse->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > regRowMapCoarse = regProlongCoarse->getColMap();
        // Get pre-communicated communication patterns for the fast MatVec
        const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = smootherParams[l+1]->get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
        const RCP<Import>        regionInterfaceImporter = smootherParams[l+1]->get<RCP<Import>>("Fast MatVec: interface importer");

        coarseRegX = VectorFactory::Build(regRowMapCoarse, true);
        coarseRegB = VectorFactory::Build(regRowMapCoarse, true);

        regProlongCoarse->apply(*regRes, *coarseRegB, Teuchos::TRANS, SC_ONE, SC_ZERO, true, regionInterfaceImporter, regionInterfaceLIDs);
          // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
          // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegB[j]->getMap()));
      }

      tm = Teuchos::null;
      bool coarseZeroInitGuess = true;

      // Call V-cycle recursively
      vCycle(l+1, numLevels, cycleType, regHierarchy,
             coarseRegX, coarseRegB,
             smootherParams, coarseZeroInitGuess, coarseSolverData, hierarchyData);

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 6 - transfer coarse to fine")));

      // Transfer coarse level correction to fine level
      RCP<Vector> regCorrection;
      {
        RCP<Level> levelCoarse = regHierarchy->GetLevel(l+1);
        RCP<Matrix>                                             regProlongCoarse     = levelCoarse->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
        // Get pre-communicated communication patterns for the fast MatVec
        const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = smootherParams[l]->get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
        const RCP<Import>        regionInterfaceImporter = smootherParams[l]->get<RCP<Import>>("Fast MatVec: interface importer");

        regCorrection = VectorFactory::Build(regRowMap, true);
        regProlongCoarse->apply(*coarseRegX, *regCorrection, Teuchos::NO_TRANS, SC_ONE, SC_ZERO, false, regionInterfaceImporter, regionInterfaceLIDs);
          // TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegX[j]->getMap()));
          // TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regCorrection[j]->getMap()));
      }

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 7 - add coarse grid correction")));

      // apply coarse grid correction
      fineRegX->update(SC_ONE, *regCorrection, SC_ONE);
      if (coarseZeroInitGuess) zeroInitGuess = true;

//    std::cout << "level: " << l << std::endl;

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: 8 - post-smoother")));

      // post-smoothing
      smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrix,
                    regRowMap, regRowImporter, zeroInitGuess);

      tm = Teuchos::null;

    }
  } else {

    // Coarsest grid solve

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(0);

    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("vCycle: * - coarsest grid solve")));

    // std::cout << "Applying coarse colver" << std::endl;

    const std::string coarseSolverType = coarseSolverData->get<std::string>("coarse solver type");
    if (coarseSolverType == "smoother") {
      smootherApply(smootherParams[l], fineRegX, fineRegB, regMatrix,
                  regRowMap, regRowImporter, zeroInitGuess);
    } else {
      zeroInitGuess = false;
      // First get the Xpetra vectors from region to composite format
      RCP<const Map> coarseRowMap = coarseSolverData->get<RCP<const Map> >("compCoarseRowMap");
      RCP<Vector> compX = VectorFactory::Build(coarseRowMap, true);
      RCP<Vector> compRhs = VectorFactory::Build(coarseRowMap, true);
      {
        RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(regInterfaceScalings->getMap());
        inverseInterfaceScaling->reciprocal(*regInterfaceScalings);
        fineRegB->elementWiseMultiply(SC_ONE, *fineRegB, *inverseInterfaceScaling, SC_ZERO);

        regionalToComposite(fineRegB, compRhs, regRowImporter);
      }

      if (coarseSolverType == "direct")
      {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)

        using DirectCoarseSolver = Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >;
        RCP<DirectCoarseSolver> coarseSolver = coarseSolverData->get<RCP<DirectCoarseSolver> >("direct solver object");

        TEUCHOS_TEST_FOR_EXCEPT_MSG(coarseRowMap->lib()!=Xpetra::UseTpetra,
            "Coarse solver requires Tpetra/Amesos2 stack.");
        TEUCHOS_ASSERT(!coarseSolver.is_null());

        // using Utilities = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

        // From here on we switch to Tpetra for simplicity
        // we could also implement a similar Epetra branch
        using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

      //    *fos << "Attempting to use Amesos2 to solve the coarse grid problem" << std::endl;
        RCP<Tpetra_MultiVector> tX = Utilities::MV2NonConstTpetraMV2(*compX);
        RCP<const Tpetra_MultiVector> tB = Utilities::MV2TpetraMV(compRhs);

        /* Solve!
         *
         * Calling solve() on the coarseSolver should just do a triangular solve, since symbolic
         * and numeric factorization are supposed to have happened during hierarchy setup.
         * Here, we just check if they're done and print message if not.
         *
         * We don't have to change the map of tX and tB since we have configured the Amesos2 solver
         * during its construction to work with non-continuous maps.
         */
        if (not coarseSolver->getStatus().symbolicFactorizationDone())
          *fos << "Symbolic factorization should have been done during hierarchy setup, "
              "but actually is missing. Anyway ... just do it right now." << std::endl;
        if (not coarseSolver->getStatus().numericFactorizationDone())
          *fos << "Numeric factorization should have been done during hierarchy setup, "
              "but actually is missing. Anyway ... just do it right now." << std::endl;
        coarseSolver->solve(tX.ptr(), tB.ptr());
#else
        *fos << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++\n"
             << "+ Coarse level direct solver requires Tpetra and Amesos2.   +\n"
             << "+ Skipping the coarse level solve.                          +\n"
             << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++"
             << std::endl;
#endif
      }
      else if (coarseSolverType == "amg") // use AMG as coarse level solver
      {
        const bool coarseSolverRebalance = coarseSolverData->get<bool>("coarse solver rebalance");

        // Extract the hierarchy from the coarseSolverData
        RCP<Hierarchy> amgHierarchy = coarseSolverData->get<RCP<Hierarchy>>("amg hierarchy object");

        // Run a single V-cycle
        if(coarseSolverRebalance==false){
          amgHierarchy->Iterate(*compRhs, *compX, 1, true);

        } else {
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
          RCP<const Import> rebalanceImporter = coarseSolverData->get<RCP<const Import> >("rebalanceImporter");

          // TODO: These vectors could be cached to improve performance
          RCP<Vector> rebalancedRhs = VectorFactory::Build(rebalanceImporter->getTargetMap());
          RCP<Vector> rebalancedX = VectorFactory::Build(rebalanceImporter->getTargetMap(), true);
          rebalancedRhs->doImport(*compRhs, *rebalanceImporter, Xpetra::INSERT);

          rebalancedRhs->replaceMap(rebalancedRhs->getMap()->removeEmptyProcesses());
          rebalancedX->replaceMap(rebalancedX->getMap()->removeEmptyProcesses());

          if(!amgHierarchy.is_null()){
            amgHierarchy->Iterate(*rebalancedRhs, *rebalancedX, 1, true);
          }

          rebalancedX->replaceMap(rebalanceImporter->getTargetMap());
          compX->doExport(*rebalancedX, *rebalanceImporter, Xpetra::INSERT);
#else
          amgHierarchy->Iterate(*compRhs, *compX, 1, true);
#endif
        }
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Unknown coarse solver type.");
      }

      // Transform back to region format
      RCP<Vector> quasiRegX;
      compositeToRegional(compX, quasiRegX, fineRegX,
                          regRowMap,
                          regRowImporter);

      tm = Teuchos::null;
    }
  }

  return;
} // vCycle

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void solveRegionProblem(const double tol, const bool scaleResidualHist, const int maxIts,
                        const std::string cycleType, const std::string convergenceLog,
                        RCP<Teuchos::ParameterList>& coarseSolverData,
                        Array<RCP<Teuchos::ParameterList> >& smootherParams,
                        RCP<Teuchos::ParameterList> hierarchyData,
                        RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & regHierarchy,
                        RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X,
                        RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& B) {

  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;
  using SC = Scalar;

  using Map    = Xpetra::Map<LO,GO,NO>;
  using Import = Xpetra::Import<LO,GO,NO>;
  using Matrix = Xpetra::Matrix<SC,LO,GO,NO>;
  using Vector = Xpetra::Vector<SC,LO,GO,NO>;
  using VectorFactory = Xpetra::VectorFactory<SC,LO,GO,NO>;

  using Level = MueLu::Level;

  using STS = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename STS::magnitudeType;
  // const Scalar zero = STS::zero();
  const Scalar one  = STS::one();

  // we start by extracting some basic data from the hierarchy
  const int numLevels = regHierarchy->GetNumLevels();
  RCP<Level> level0 = regHierarchy->GetLevel(0);
  RCP<Matrix> regMat  = level0->Get<RCP<Matrix> >("A");
  RCP<const Map> revisedRowMap  = regMat->getRowMap();
  RCP<Import> rowImport = level0->Get<RCP<Import> >("rowImport");
  RCP<const Map> dofMap = X->getMap();
  const int myRank = dofMap->getComm()->getRank();

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

  // We first use the non-level container variables to setup the fine grid problem.
  // This is ok since the initial setup just mimics the application and the outer
  // Krylov method.
  //
  // We switch to using the level container variables as soon as we enter the
  // recursive part of the algorithm.
  //

  // Composite residual vector
  RCP<Vector> compRes = VectorFactory::Build(dofMap, true);

  // transform composite vectors to regional layout
  Teuchos::RCP<Vector> quasiRegX;
  Teuchos::RCP<Vector> regX;
  compositeToRegional(X, quasiRegX, regX,
                      revisedRowMap, rowImport);

  RCP<Vector> quasiRegB;
  RCP<Vector> regB;
  compositeToRegional(B, quasiRegB, regB,
                      revisedRowMap, rowImport);

  RCP<Vector> regRes;
  regRes = VectorFactory::Build(revisedRowMap, true);

  /////////////////////////////////////////////////////////////////////////
  // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
  /////////////////////////////////////////////////////////////////////////

  // Prepare output of residual norm to file
  RCP<std::ofstream> log;
  if (myRank == 0)
    {
      log = rcp(new std::ofstream(convergenceLog.c_str()));
      (*log) << "# num procs = " << dofMap->getComm()->getSize() << "\n"
             << "# iteration | res-norm (scaled=" << scaleResidualHist << ")\n"
             << "#\n";
      *log << std::setprecision(16) << std::scientific;
    }

  // Print type of residual norm to the screen
  if (scaleResidualHist)
    out << "Using scaled residual norm." << std::endl;
  else
    out << "Using unscaled residual norm." << std::endl;


  // Richardson iterations
  magnitude_type normResIni = Teuchos::ScalarTraits<magnitude_type>::zero();
  const int old_precision = std::cout.precision();
  std::cout << std::setprecision(8) << std::scientific;
  int cycle = 0;

  Teuchos::RCP<Vector> regCorrect;
  regCorrect = VectorFactory::Build(revisedRowMap, true);
  for (cycle = 0; cycle < maxIts; ++cycle)
    {
      const Scalar SC_ZERO = Teuchos::ScalarTraits<SC>::zero();
      regCorrect->putScalar(SC_ZERO);
      // Get Stuff out of Hierarchy
      RCP<Level> level = regHierarchy->GetLevel(0);
      RCP<Vector> regInterfaceScalings = level->Get<RCP<Vector> >("regInterfaceScalings");
      // check for convergence
      {
        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////
        computeResidual(regRes, regX, regB, regMat, *smootherParams[0]);
        scaleInterfaceDOFs(regRes, regInterfaceScalings, true);

        compRes = VectorFactory::Build(dofMap, true);
        regionalToComposite(regRes, compRes, rowImport);

        typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();
        if(cycle == 0) { normResIni = normRes; }

        if (scaleResidualHist)
          normRes /= normResIni;

        // Output current residual norm to screen (on proc 0 only)
        out << cycle << "\t" << normRes << std::endl;
        if (myRank == 0)
          (*log) << cycle << "\t" << normRes << "\n";

        if (normRes < tol)
          break;
      }

      /////////////////////////////////////////////////////////////////////////
      // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
      /////////////////////////////////////////////////////////////////////////

      bool zeroInitGuess = true;
      scaleInterfaceDOFs(regRes, regInterfaceScalings, false);
      vCycle(0, numLevels, cycleType, regHierarchy,
             regCorrect, regRes,
             smootherParams, zeroInitGuess, coarseSolverData, hierarchyData);

      regX->update(one, *regCorrect, one);
    }
  out << "Number of iterations performed for this solve: " << cycle << std::endl;

  std::cout << std::setprecision(old_precision);
  std::cout.unsetf(std::ios::fixed | std::ios::scientific);
}

#endif // MUELU_SOLVEREGIONHIERARCHY_DEF_HPP
