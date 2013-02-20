/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack2_Parameters.hpp"

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
  params.set("schwarz: compute condest", true);
  params.set("schwarz: combine mode", "Zero"); // use string mode for this
  params.set("schwarz: reordering type", "none");
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
  params.set("krylov: number of iterations",(int)5);
  params.set("krylov: residual tolerance",(double)0.001);
  params.set("krylov: block size",(int)1);
  params.set("krylov: zero starting solution", true);
  params.set("krylov: preconditioner type",(int)1);

  // Ifpack2_METISPartitioner.hpp
  params.set("partitioner: use symmetric graph", true);
}

}//namespace Ifpack2

