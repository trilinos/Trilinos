// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_Parameters.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_MultiVector.hpp>

namespace Ifpack2 {

void getValidParameters(Teuchos::ParameterList &params) {
  using STS = Teuchos::ScalarTraits<double>;

  // params.clear();
  Teuchos::ParameterList empty;
  params = empty;

  // ============================================================ //
  // Parameters are reported from each used file in IFPACK2. Files //
  // are listed in alphabetical order, first all *.cpp, then *.hpp. //
  // Some options not very tested or documented anywhere          //
  // are not reported here.                                       //
  // ============================================================ //

  // Ifpack2_IlukGraph.hpp
  params.set("fact: iluk level-of-fill", 1);
  params.set("fact: iluk level-of-overlap", 0);

  // Ifpack2_Chebyshev
  params.set("chebyshev: max eigenvalue", STS::nan());
  params.set("chebyshev: ratio eigenvalue", STS::nan());
  params.set("chebyshev: min eigenvalue", 30.0);
  params.set("chebyshev: degree", 1);
  params.set("chebyshev: eigenvalue max iterations", 10);
  params.set("chebyshev: eigenvalue relative tolerance", 0.0);
  params.set("chebyshev: eigenvalue keep vector", false);
  params.set("chebyshev: assume matrix does not change", false);
  // params.set("chebyshev: operator inv diagonal",Teuchos::null);
  params.set("chebyshev: min diagonal value", STS::eps());
  params.set("chebyshev: zero starting solution", true);
  params.set("chebyshev: use native spmv", false);
  params.set("chebyshev: algorithm", "first");

  // Ifpack2_Amesos.cpp
  params.set("amesos: solver type", "Amesos_Klu");

  // Ifpack2_IC.cpp
  params.set("fact: level-of-fill", 1);
  params.set("fact: absolute threshold", 0.0);
  params.set("fact: relative threshold", 0.0);
  params.set("fact: drop tolerance", 0.0);

  // Ifpack2_ICT.cpp
  params.set("fact: ict level-of-fill", 1.0);
  params.set("fact: absolute threshold", 0.0);
  params.set("fact: relative threshold", 1.0);
  params.set("fact: relax value", 0.0);
  params.set("fact: drop tolerance", 0.0);

  // Ifpack2_ILU.cpp
  params.set("fact: level-of-fill", 0);
  params.set("fact: absolute threshold", 0.0);
  params.set("fact: relative threshold", 1.0);
  params.set("fact: relax value", 0.0);

  // Ifpack2_ILUT.cpp
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: absolute threshold", 0.0);
  params.set("fact: relative threshold", 1.0);
  params.set("fact: relax value", 0.0);
  params.set("fact: type", "serial");
  params.sublist("parallel ILUT options"); // FIXME this should be validated

  // Ifpack2_LocalSparseTriangularSolver.cpp
  params.set("trisolver: type", "Internal");
  params.set("trisolver: block size", 1);
  params.set("trisolver: reverse U", false);

  // Overlapping partitioner
  params.set("partitioner: local parts", 1);
  params.set("partitioner: overlap", 0);
  params.set("partitioner: print level", 0);

  // Ifpack2_Relaxation.cpp
  params.set("relaxation: container", "TriDi");
  params.set("relaxation: type", "Jacobi");
  params.set("relaxation: sweeps", 1);
  params.set("relaxation: direction", "forward");
  params.set("relaxation: damping factor", 1.0);
  params.set("relaxation: min diagonal value", 1.0);
  params.set("relaxation: zero starting solution", true);
  params.set("relaxation: backward mode", false);
  params.set("relaxation: use l1", false);
  params.set("relaxation: l1 eta", 1.5);
  params.set("relaxation: banded container superdiagonals", -1);
  params.set("relaxation: banded container subdiagonals", -1);
  params.set("relaxation: mtgs cluster size", 1);
  params.set("relaxation: mtgs coloring algorithm", "Default");
  params.set("relaxation: long row threshold", 0);

  // Ifpack2_SPARSKIT.cpp
  // ap 25 May 2016: all SPARSKIT for backwards compatibility ONLY
  params.set("fact: sparskit: lfil", 0);
  params.set("fact: sparskit: tol", 0.0);
  params.set("fact: sparskit: droptol", 0.0);
  params.set("fact: sparskit: permtol", 0.1);
  params.set("fact: sparskit: alph", 0.0);
  params.set("fact: sparskit: mbloc", -1);
  params.set("fact: sparskit: type", "ILUT");

  // Additive Schwarz preconditioner
  params.set("schwarz: compute condest",
             false); // mfh 24 Mar 2015: for backwards compatibility ONLY
  params.set("schwarz: combine mode", "ZERO"); // use string mode for this
  params.set("schwarz: use reordering", true);
  params.set("schwarz: filter singletons", false);
  params.set("schwarz: overlap level", 0);
  params.set("schwarz: num iterations", 1);

  params.set("subdomain solver name", "");
  params.set("inner solver name", "");
  params.set("schwarz: subdomain solver name", "");
  params.set("schwarz: inner solver name", "");

  Teuchos::ParameterList dummyListSubdomain;

  const std::vector<std::string> subdomainSolverParameterNames = {
      "inner preconditioner parameters", "subdomain solver parameters",
      "schwarz: inner preconditioner parameters",
      "schwarz: subdomain solver parameters"};
  for (auto &subdomainSolverParameterName : subdomainSolverParameterNames) {
    params.set(subdomainSolverParameterName, dummyListSubdomain);
    params.sublist(subdomainSolverParameterName).disableRecursiveValidation();
  }

  Teuchos::ParameterList dummyListReordering;
  params.set("schwarz: reordering list", dummyListReordering);
  // Ifpack2 doesn't attempt to validate options for Zoltan2
  params.sublist("schwarz: reordering list").disableRecursiveValidation();

  // Ifpack2_BlockRelaxation.hpp
  // params.set("relaxation: type", "Jacobi"); // already set
  // params.set("relaxation: sweeps", 1); // already set
  // params.get("relaxation: damping factor", 1.0); // already set
  // params.get("relaxation: zero starting solution", true); // already set
  params.set("partitioner: type", "greedy");
  params.set("zoltan2: algorithm", "phg");
  params.set("partitioner: local parts", 1);
  params.set("partitioner: overlap", 0);
  Teuchos::Array<Teuchos::ArrayRCP<int>> tmp0;
  params.set("partitioner: parts", tmp0);
  params.set("partitioner: maintain sparsity", false);
  params.set("block relaxation: decouple dofs", false);

  // Ifpack2_METISPartitioner.hpp
  // ap 25 May 2016: all METIS for backwards compatibility ONLY
  params.set("partitioner: use symmetric graph", true);

  // Ifpack2_Details_Amesos2Wrapper
  Teuchos::ParameterList dummyList;
  params.set("Amesos2", dummyList);
  params.sublist("Amesos2").disableRecursiveValidation();
  params.set("Amesos2 solver name", "KLU2");

  // Ifpack2_Details_UserPartitioner.hpp
  Teuchos::ArrayRCP<int> tmp;
  params.set("partitioner: map", tmp);

  // Ifpack2_LinePartitioner.hpp (FIXME)
  params.set("partitioner: line detection threshold", 0.0);
  params.set("partitioner: PDE equations", 1);
  Teuchos::RCP<Tpetra::MultiVector<>> dummy;
  params.set("partitioner: coordinates", dummy);

  // Ifpack2_Hypre.hpp
  params.set("hypre: Solver", "PCG");
  params.set("hypre: Preconditioner", "Euclid");
  params.set("hypre: SolveOrPrecondition", "Solver");
  params.sublist("hypre: Solver functions").disableRecursiveValidation();

  params.sublist("hypre: Preconditioner functions")
      .disableRecursiveValidation();
  params.sublist("Operators").disableRecursiveValidation();
  params.sublist("Coordinates").disableRecursiveValidation();
  params.set("hypre: Dump", false);
  params.set("hypre: SetPreconditioner", false);
  params.set("hypre: NumFunctions", 0);
}

} // namespace Ifpack2
