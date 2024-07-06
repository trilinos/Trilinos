// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using SC = Tpetra::MultiVector<>::scalar_type;
  using crs_matrix_type = Tpetra::CrsMatrix<SC>;
  using map_type = Tpetra::Map<>;
  using OP = Tpetra::Operator<SC>;
  using MV = Tpetra::MultiVector<SC>;
  using LO = map_type::local_ordinal_type; // int LO;
  using GO = map_type::global_ordinal_type; // int GO;
  // typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto TeuchosComm = Tpetra::getDefaultComm ();

    int NumMyElements = 0;
    if (TeuchosComm->getRank()==0) {
      NumMyElements = 10;
    }
    Array<GO> uniqueMapArray(NumMyElements);
    for (LO i=0; i<uniqueMapArray.size(); i++) {
      uniqueMapArray[i] = i;
    }

    const Tpetra::global_size_t INVALID =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const map_type> UniqueMap = rcp (new map_type (INVALID, uniqueMapArray, 0, TeuchosComm));

    RCP<crs_matrix_type> K = rcp (new crs_matrix_type (UniqueMap, 10));

    for (LO i = 0; i < static_cast<LO>(UniqueMap->getLocalNumElements()); ++i) {
      LO numEntries = 10-i;
      Array<GO> indicesArray(numEntries);
      Array<SC> valuesArray(numEntries);
      for (LO k=0; k<numEntries; k++) {
        indicesArray[k] = k+i;
        valuesArray[k] = 1;
      }
      K->insertGlobalValues(UniqueMap->getGlobalElement(i), indicesArray(),
                            valuesArray());
    }
    K->fillComplete();

    //Solve with GMRES
    RCP<MV> Solution = rcp(new MV(UniqueMap,1));
    RCP<MV> RightHandSide = rcp(new MV(UniqueMap,1));

    Solution->putScalar(0.0);
    RightHandSide->putScalar(1.0);

    using Belos::LinearProblem;
    RCP<LinearProblem<SC,MV,OP> > belosLinearProblem =
      rcp (new LinearProblem<SC,MV,OP> (K, Solution, RightHandSide));
    belosLinearProblem->setProblem();

    RCP<ParameterList> solverParameterList = rcp(new ParameterList());
    solverParameterList->set("Convergence Tolerance",1.0e-12);
    solverParameterList->set("Verbosity",47);
    solverParameterList->set("Output Frequency",1);
    solverParameterList->set("Output Style",1);

    Belos::BlockGmresSolMgr<SC,MV,OP> solver(belosLinearProblem,
                                             solverParameterList);
    solver.solve();
  }

  return(EXIT_SUCCESS);
}
