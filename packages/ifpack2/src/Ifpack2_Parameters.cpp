/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

namespace Ifpack2 {

void getValidParameters(Teuchos::ParameterList& params)
{
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
  params.set("fact: iluk level-of-fill", (int)1);
  params.set("fact: iluk level-of-overlap", (int)0);

  // Ifpack2_Amesos.cpp
  params.set("amesos: solver type", "Amesos_Klu");

  // Ifpack2_IC.cpp
  params.set("fact: level-of-fill", (int)1);
  params.set("fact: absolute threshold", (double)0.0);
  params.set("fact: relative threshold", (double)0.0);
  params.set("fact: drop tolerance", (double)0.0);

  // Ifpack2_ICT.cpp
  params.set("fact: ict level-of-fill", (double)1.0);
  params.set("fact: absolute threshold", (double)0.0);
  params.set("fact: relative threshold", (double)1.0);
  params.set("fact: relax value", (double)0.0);
  params.set("fact: drop tolerance", (double)0.0);

  // Ifpack2_ILU.cpp
  params.set("fact: level-of-fill", (int)0);
  params.set("fact: absolute threshold", (double)0.0);
  params.set("fact: relative threshold", (double)1.0);
  params.set("fact: relax value", (double)0.0);

  // Ifpack2_ILUT.cpp
  params.set("fact: ilut level-of-fill", (double)1.0);
  params.set("fact: absolute threshold", (double)0.0);
  params.set("fact: relative threshold", (double)1.0);
  params.set("fact: relax value", (double)0.0);

  // Ifpack2_METISPartitioner.cpp
  params.set("partitioner: local parts", (int)1);
  params.set("partitioner: overlap", (int)0);
  params.set("partitioner: print level", (int)0);

  // Ifpack2_Relaxation.cpp
  params.set("relaxation: type", "Jacobi");
  params.set("relaxation: sweeps", (int)1);
  params.set("relaxation: damping factor", (double)1.0);
  params.set("relaxation: min diagonal value", (double)1.0);
  params.set("relaxation: zero starting solution", true);
  params.set("relaxation: backward mode",false);
  params.set("relaxation: use l1",false);
  params.set("relaxation: l1 eta",(double)1.5);

  // Ifpack2_SPARSKIT.cpp
  params.set("fact: sparskit: lfil", (int)0);
  params.set("fact: sparskit: tol", (double)0.0);
  params.set("fact: sparskit: droptol", (double)0.0);
  params.set("fact: sparskit: permtol", (double)0.1);
  params.set("fact: sparskit: alph", (double)0.0);
  params.set("fact: sparskit: mbloc", (int)(-1));
  params.set("fact: sparskit: type", ("ILUT"));

  // Additive Schwarz preconditioner
  params.set("schwarz: compute condest", false); // mfh 24 Mar 2015: for backwards compatibility ONLY
  params.set("schwarz: combine mode", "ZERO"); // use string mode for this
  params.set("schwarz: use reordering", true);
  params.set("schwarz: filter singletons", false);
  params.set("schwarz: overlap level", (int)0);

  // Ifpack2_BlockRelaxation.hpp
  // params.set("relaxation: type", "Jacobi"); // already set
  // params.set("relaxation: sweeps", 1); // already set
  // params.get("relaxation: damping factor", 1.0); // already set
  // params.get("relaxation: zero starting solution", true); // already set
  params.set("partitioner: type", "greedy");
  params.set("partitioner: local parts", (int)1);
  params.set("partitioner: overlap", (int)0);

  // Krylov smoother
  params.set("krylov: iteration type",(int)1);
  params.set("krylov: number of iterations",(int)5);
  params.set("krylov: residual tolerance",(double)0.001);
  params.set("krylov: block size",(int)1);
  params.set("krylov: zero starting solution", true);
  params.set("krylov: preconditioner type",(int)1);

  // Ifpack2_METISPartitioner.hpp
  params.set("partitioner: use symmetric graph", true);

  // Ifpack2_Details_UserPartitioner.hpp
  Teuchos::ArrayRCP<int> tmp;
  params.set("partitioner: map", tmp);
}

}//namespace Ifpack2

