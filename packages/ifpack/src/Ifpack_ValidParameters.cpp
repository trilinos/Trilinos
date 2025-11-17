/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#include "Ifpack_ValidParameters.h"

Teuchos::ParameterList Ifpack_GetValidParameters()
{
  Teuchos::ParameterList List; // empty list

  // ============================================================ //
  // Parameters are reported from each used file in IFPACK. Files //
  // are listed in alphabetical order, first all *.cpp, then *.h. //
  // Some options not very tested or documented anywhere          //
  // are not reported here.                                       //
  // ============================================================ //
  
  // Ifpack_Amesos.cpp
  List.set("amesos: solver type", "Amesos_Klu");

  // Ifpack_IC.cpp
  List.set("fact: level-of-fill", (int)1);
  List.set("fact: absolute threshold", (double)0.0);
  List.set("fact: relative threshold", (double)0.0);
  List.set("fact: drop tolerance", (double)0.0);

  // Ifpack_ICT.cpp
  List.set("fact: ict level-of-fill", (double)1.0);
  List.set("fact: absolute threshold", (double)0.0);
  List.set("fact: relative threshold", (double)1.0);
  List.set("fact: relax value", (double)0.0);
  List.set("fact: drop tolerance", (double)0.0);

  // Ifpack_ILU.cpp
  List.set("fact: level-of-fill", (int)0);
  List.set("fact: absolute threshold", (double)0.0);
  List.set("fact: relative threshold", (double)1.0);
  List.set("fact: relax value", (double)0.0);

  // Ifpack_ILUT.cpp
  List.set("fact: ilut level-of-fill", (double)1.0);
  List.set("fact: absolute threshold", (double)0.0);
  List.set("fact: relative threshold", (double)1.0);
  List.set("fact: relax value", (double)0.0);

#ifdef HAVE_IFPACK_SUPERLU
  // Ifpack_SILU.cpp
  List.set("fact: drop tolerance",1e-4);
  List.set("fact: zero pivot threshold",1e-2);
  List.set("fact: maximum fill factor",10.0);
  List.set("fact: silu drop rule",9);
#endif

  // Ifpack_METISPartitioner.cpp
  List.set("partitioner: local parts", (int)1);
  List.set("partitioner: overlap", (int)0);
  List.set("partitioner: print level", (int)0);

  // Ifpack_PointRelaxation.cpp
  List.set("relaxation: type", "Jacobi");
  List.set("relaxation: sweeps", (int)1);
  List.set("relaxation: damping factor", (double)1.0);
  List.set("relaxation: min diagonal value", (double)1.0);
  List.set("relaxation: zero starting solution", true);
  List.set("relaxation: backward mode",false);
  List.set("relaxation: use l1",false);
  List.set("relaxation: l1 eta",(double)1.5);

  // Ifpack_SPARSKIT.cpp
  List.set("fact: sparskit: lfil", (int)0);
  List.set("fact: sparskit: tol", (double)0.0);
  List.set("fact: sparskit: droptol", (double)0.0);
  List.set("fact: sparskit: permtol", (double)0.1);
  List.set("fact: sparskit: alph", (double)0.0);
  List.set("fact: sparskit: mbloc", (int)(-1));
  List.set("fact: sparskit: type", ("ILUT"));

  // Ifpack_Chebyshev.cpp
  List.set("chebyshev: ratio eigenvalue", (double)30.0);
  List.set("chebyshev: min eigenvalue", (double)0.0);
  List.set("chebyshev: max eigenvalue", (double)-1.0);
  List.set("chebyshev: degree", (int)1);
  List.set("chebyshev: min diagonal value", (double)0.0);
  List.set("chebyshev: zero starting solution", (bool)true);

  // Additive Schwarz preconditioner
  List.set("schwarz: compute condest", true);
  List.set("schwarz: combine mode", "Zero"); // use std::string mode for this
  List.set("schwarz: reordering type", "none");
  List.set("schwarz: filter singletons", false);

  // Ifpack_BlockRelaxation.h
  // List.set("relaxation: type", "Jacobi"); // already set
  // List.set("relaxation: sweeps", 1); // already set
  // List.get("relaxation: damping factor", 1.0); // already set
  // List.get("relaxation: zero starting solution", true); // already set
  List.set("partitioner: type", "greedy");
  List.set("partitioner: local parts", (int)1);
  List.set("partitioner: overlap", (int)0);
  List.set("partitioner: keep singletons",false);

  // Ifpack_METISPartitioner.h
  List.set("partitioner: use symmetric graph", true);

  // Krylov smoother
  List.set("krylov: iterations",(int)5);
  List.set("krylov: tolerance",(double)0.001);
  List.set("krylov: solver",(int)1);
  List.set("krylov: preconditioner",(int)0);
  List.set("krylov: number of sweeps",(int)1);
  List.set("krylov: block size",(int)1);
  List.set("krylov: damping parameter",(double)1.0);
  List.set("krylov: zero starting solution",true);

  // Ifpack_Hypre.cpp
  List.set("hypre: Solver", "PCG");
  List.set("hypre: Preconditioner", "Euclid");
  List.set("hypre: SolveOrPrecondition", "Solver");
  List.sublist("hypre: Solver functions").disableRecursiveValidation();

  List.sublist("hypre: Preconditioner functions").disableRecursiveValidation();
  List.sublist("Operators").disableRecursiveValidation();
  List.sublist("Coordinates").disableRecursiveValidation();
  List.set("hypre: Dump", false);
  List.set("hypre: SetPreconditioner", false);
  List.set("hypre: NumFunctions", 0);


  return(List);
}

