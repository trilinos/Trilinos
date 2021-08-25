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
