#include "Ifpack_ValidParameters.h"

#ifdef HAVE_IFPACK_TEUCHOS

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

  // Ifpack_SPARSKIT.cpp
  List.set("fact: sparskit: lfil", (int)0);
  List.set("fact: sparskit: tol", (double)0.0);
  List.set("fact: sparskit: droptol", (double)0.0);
  List.set("fact: sparskit: permtol", (double)0.1);
  List.set("fact: sparskit: alph", (double)0.0);
  List.set("fact: sparskit: mbloc", (int)(-1));
  List.set("fact: sparskit: type", ("ILUT"));

  // Additive Schwarz preconditioner
  List.set("schwarz: compute condest", true);
  List.set("schwarz: combine mode", "Zero"); // use string mode for this
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

  // Ifpack_METISPartitioner.h
  List.set("partitioner: use symmetric graph", true);

  return(List);
}

#endif
