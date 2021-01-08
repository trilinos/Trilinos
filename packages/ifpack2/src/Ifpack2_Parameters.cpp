/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack2_Parameters.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace Ifpack2 {

void getValidParameters(Teuchos::ParameterList& params)
{
  using STS = Teuchos::ScalarTraits<double>;

  //params.clear();
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
  params.set("chebyshev: assume matrix does not change", false);
  // params.set("chebyshev: operator inv diagonal",Teuchos::null);
  params.set("chebyshev: min diagonal value", STS::eps());
  params.set("chebyshev: zero starting solution", true);

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

  // Ifpack2_SPARSKIT.cpp
  // ap 25 May 2016: all SPARSKIT for backwards compatibility ONLY
  params.set("fact: sparskit: lfil",    0);
  params.set("fact: sparskit: tol",     0.0);
  params.set("fact: sparskit: droptol", 0.0);
  params.set("fact: sparskit: permtol", 0.1);
  params.set("fact: sparskit: alph",    0.0);
  params.set("fact: sparskit: mbloc",   -1);
  params.set("fact: sparskit: type",    "ILUT");

  // Additive Schwarz preconditioner
  params.set("schwarz: compute condest", false); // mfh 24 Mar 2015: for backwards compatibility ONLY
  params.set("schwarz: combine mode", "ZERO"); // use string mode for this
  params.set("schwarz: use reordering", true);
  params.set("schwarz: filter singletons", false);
  params.set("schwarz: overlap level", 0);

  // Ifpack2_BlockRelaxation.hpp
  // params.set("relaxation: type", "Jacobi"); // already set
  // params.set("relaxation: sweeps", 1); // already set
  // params.get("relaxation: damping factor", 1.0); // already set
  // params.get("relaxation: zero starting solution", true); // already set
  params.set("partitioner: type", "greedy");
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
  params.set("Amesos2",dummyList);
  params.sublist("Amesos2").disableRecursiveValidation();
  params.set("Amesos2 solver name", "KLU2");

  // Ifpack2_Details_UserPartitioner.hpp
  Teuchos::ArrayRCP<int> tmp;
  params.set("partitioner: map", tmp);

  // Ifpack2_LinePartitioner.hpp (FIXME)
  params.set("partitioner: line detection threshold", 0.0);
  params.set("partitioner: PDE equations", 1);
  Teuchos::RCP<Tpetra::MultiVector<> > dummy;
  params.set("partitioner: coordinates",dummy);

  // Ifpack2_Hypre.hpp
  params.set("hypre: Solver", "PCG");
  params.set("hypre: Preconditioner", "Euclid");
  params.set("hypre: SolveOrPrecondition", "Solver");
  params.sublist("hypre: Solver functions").disableRecursiveValidation();

  params.sublist("hypre: Preconditioner functions").disableRecursiveValidation();
  params.set("hypre: SetPreconditioner", false);
  params.set("hypre: NumFunctions", 0);
}

}//namespace Ifpack2

